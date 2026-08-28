#!/usr/bin/env python3
"""Generate and submit the production Tet2025 CellBouncer SLURM workflow.

The working scientific execution graph is retained for demultiplexing,
contamination estimation, ploidy/refinement, post-hoc audit,
identity reconciliation, and mitochondrial analysis.  Outputs now live with
their source libraries under ``mapping_output/Tet_2025_Multiome-RNA_N``;
cross-library products live under ``mapping_output/aggregate_library_analysis``.

``AMBIENT_PLOTS`` is a plot-only stage for one to eight selected contamination
conditions.  It creates descriptive cross-library plots and the installed
CellBouncer ``plot/contam.R`` summaries without performing condition ranking.
CONTAM and AMBIENT_PLOTS can select original demux identities, reconciled
identities, or both for a paired identity-sensitivity analysis.

``AMBIENT_VALIDATE`` performs the production identity/ambient control.  For
each library with an applied reconciliation, it fits one deterministic
empty-droplet profile from the applied production roster, freezes that profile,
and rescales the original and reconciled assignments against it.  This isolates
the assignment effect without the discarded presence-screen experiment or
profile-refitting confounding.

``AMBIENT_SWAP_TEST`` discovers supported unexpected identity populations from
the validated reconciliation event table.  It selects exact-identity events
that already cleared the event-mass, biological-admissibility, occupancy, and
evidence policy.  It then independently requires QC-passing NN tetraploid
support for every diploid-to-heterotypic proposal and reapplies the event-mass
floor after cell-level eligibility filtering.  Eligible assignments are
evaluated with both a shared frozen empty-drop profile (G/H) and independently
fitted profiles (J/K).  ``identity_reconciliation_figures.py`` owns all
analysis and figures, including combined and per-identity views.  Candidate-
level views do not multiply estimator jobs by event count, and no donor or
library candidate list is hard-coded or supplied manually.

``GEX_AMBIENT`` reuses completed per-cell contamination estimates to infer an
ambient gene-expression profile with CellBouncer's ``ambient_rna_gex`` path.
It can consume an explicit cluster file, reuse a cluster column in the H5AD
configured for ``PLOIDY_NN``, or create RNA-only numbered Leiden nuisance
clusters directly from every barcode in each library's raw filtered MEX.
Automatic clustering uses a permissive anchor-and-transfer design rather than
the production QC exclusion mask.  The C++ fit receives an exact cohort-subset
MEX containing only cells with expression, clusters, and valid contamination
rates.  Versioned gene-distribution summaries test whether ambient RNA is
concentrated in a few genes or broadly distributed.

``CLEANUP_RESULTS`` is a standalone maintenance stage for the entire pipeline.
It keeps the newest completed result for each versioned aggregate workflow,
keeps the newest valid per-library generated profile/result generation, removes
stale figures from the retained plot tree, and retains only the newest log and
generated-script attempt for each logical job.  Canonical in-place products and
shared inputs (CONDF, condition indexes, reconciliation decisions, posthoc,
hybrid, ploidy, tetra-refine, mitochondrial outputs, and migration journals)
are never cleanup targets.

The default condition is the current focused production candidate,
``IND_CK_RF_SX0_GATED_RFREE_PFIT``.  Every historical condition remains
explicitly selectable through ``--conditions`` or the named condition sets.
Dry-run generation remains the default; pass ``--submit`` to launch jobs.

When CONDF or DEMUX has pending work, the default managed-daemon lifecycle
creates run-scoped main/HET/species shared-memory segments on the configured
daemon nodes, gates consumers on an atomic readiness marker, and releases only
those run-owned segments after every consumer settles.  Use
``--no-manage-vcf-daemons`` only for the historical externally managed mode.
"""

# Legacy revision prose intentionally removed.  Source history remains in git;
# this file documents the current production contract only.

import argparse
import atexit
import csv
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
from datetime import datetime
from pathlib import Path


# AMBIENT_PLOTS uses only standard-library imports during orchestration.  The
# numerical plotting stack is loaded inside its compute-node worker.
AMBIENT_PLOT_DEFAULT_CONDITION = "IND_CK_RF_SX0_GATED_RFREE_PFIT"
ORCHESTRATOR_RELEASE = "2026-08-27-identity-reconciliation-round1"
CANDIDATE_AXIS_STAGE = "IDENTITY_CANDIDATE_AXIS"
CANDIDATE_AXIS_MIN_EVIDENCE = 10
CANDIDATE_AXIS_MIN_EVIDENCE_SOURCE = (
    "AUDITED_TETRA_SCORE_CALLS_V3_DEFAULT_10"
)
AMBIENT_PROFILE_REQUIRED_VERSION = "2.6-adaptive-simplex"
AMBIENT_CONDITION_SHORT_NAMES = {
    "IND_CK_RF_SX0_GATED_RFREE_PFIT": "CK_SX0_GATED",
}
AMBIENT_PLOT_DIRNAME = "ambient_rna"
AMBIENT_CONTAMINATION_SUBDIR = "contamination"
AMBIENT_R_MODULE = "R/4.5.1"
AMBIENT_PYTHON_MODULES = ("miniforge/3", "genomics-base/latest")
AMBIENT_DEPLOYED_CONTAM_R = "/nvme/software/packages/cellbouncer/dev/plot/contam.R"
CONTAM_ASSIGNMENT_SOURCES = ("demux", "reconciled")
IDENTITY_AMBIENT_SELECTOR = "reconciliation-four-arm"
IDENTITY_AMBIENT_DIRNAME = "reconciliation_four_arm"
IDENTITY_AMBIENT_CANDIDATE_SET = "applied"
IDENTITY_AMBIENT_PLAN_VERSION = "identity_ambient_comparison_plan_V4"
_IDENTITY_AMBIENT_PLAN_CACHE = {}
IDENTITY_AMBIENT_ARMS = {
    "demux_original": {
        "arm": "A",
        "label": "A: demux assignments + original roster",
        "short_label": "A demux/original",
        "assignment_basis": "demux",
        "roster_basis": "original",
    },
    "demux_augmented": {
        "arm": "B",
        "label": "B: demux assignments + augmented roster",
        "short_label": "B demux/augmented",
        "assignment_basis": "demux",
        "roster_basis": "augmented",
    },
    "reconciled_augmented": {
        "arm": "C",
        "label": "C: reconciled assignments + augmented roster",
        "short_label": "C reconciled/augmented",
        "assignment_basis": "reconciled",
        "roster_basis": "augmented",
    },
    "reconciled_replacement": {
        "arm": "D",
        "label": "D: reconciled assignments + replacement-only roster",
        "short_label": "D reconciled/replacement",
        "assignment_basis": "reconciled",
        "roster_basis": "replacement",
        "optional": True,
    },
}
IDENTITY_AMBIENT_ARM_ORDER = tuple(IDENTITY_AMBIENT_ARMS)
AMBIENT_VALIDATION_DIRNAME = "ambient_validation"
AMBIENT_VALIDATION_PLAN_VERSION = "fixed_profile_identity_validation_V1"
AMBIENT_VALIDATION_PROFILE_MAX_EMPTY = 50000
AMBIENT_VALIDATION_PROFILE_SEED = 42
AMBIENT_VALIDATION_PROFILE_STARTS = 2
AMBIENT_VALIDATION_FIXED_ARMS = {
    "demux_fixed_empty": {
        "arm": "E",
        "label": "E: demux assignments + frozen empty-drop profile",
        "short_label": "E demux/fixed-empty",
        "source_arm": "demux_augmented",
    },
    "reconciled_fixed_empty": {
        "arm": "F",
        "label": "F: reconciled assignments + frozen empty-drop profile",
        "short_label": "F reconciled/fixed-empty",
        "source_arm": "reconciled_augmented",
    },
}
AMBIENT_SWAP_TEST_DIRNAME = "ambient_swap_test"
AMBIENT_SWAP_TEST_PLAN_VERSION = "auto_discovered_fixed_profile_swap_test_V3"
AMBIENT_SWAP_HETEROTYPIC_NN_MIN_PROB = 0.90
AMBIENT_SWAP_TEST_ARMS = {
    "original_fixed": {
        "arm": "G",
        "label": "G: original assignments + frozen candidate-roster profile",
        "short_label": "G original/fixed-profile",
        "assignment_key": "original_assignments",
        "profile_mode": "fixed",
    },
    "candidate_fixed": {
        "arm": "H",
        "label": "H: candidate-overlay assignments + same frozen profile",
        "short_label": "H proposed/fixed-profile",
        "assignment_key": "candidate_assignments",
        "profile_mode": "fixed",
    },
    "original_joint": {
        "arm": "J",
        "label": "J: original assignments + independently fitted profile",
        "short_label": "J original/joint-profile",
        "assignment_key": "original_assignments",
        "profile_mode": "joint",
    },
    "candidate_joint": {
        "arm": "K",
        "label": "K: candidate-overlay assignments + independently fitted profile",
        "short_label": "K proposed/joint-profile",
        "assignment_key": "candidate_assignments",
        "profile_mode": "joint",
    },
}

# =============================================================================
# Paths (all absolute, hardcoded)
# =============================================================================

PROJECT_ROOT = "/mnt/beegfs/tetmultiome_rna_mapped"
BEEGFS_ROOT = os.path.join(PROJECT_ROOT, "mapping_output")
AGGREGATE_ROOT = os.path.join(BEEGFS_ROOT, "aggregate_library_analysis")
CONDITION_INDEX_ROOT = os.path.join(AGGREGATE_ROOT, "condition_index")
CONDF_DIR = os.path.join(AGGREGATE_ROOT, "condf")
SOFTWARE_BIN = "/nvme/software/packages/cellbouncer/dev/bin"
SOFTWARE_PLOT = "/nvme/software/packages/cellbouncer/dev/plot"
IDENTITY_RECONCILIATION_FIGURES = os.path.join(
    SOFTWARE_PLOT, "identity_reconciliation_figures.py")
CONTAM_R = os.path.join(SOFTWARE_PLOT, "contam.R")
PANEL_METADATA = "/mnt/beegfs/tetmultiome_rna_mapped/Misc_Metadata/panel_metadata.tsv"
EXPECTED_LINES_DIR = "/mnt/beegfs/tetmultiome_rna_mapped/Misc_Metadata"
LIBRARY_CONVERSIONS_XLSX = os.path.join(EXPECTED_LINES_DIR, "Library_conversions.xlsx")
# Identity reconciliation uses the authoritative project-wide workbook colocated
# with the mapping outputs; 2025_LineMeta defines the global biological universe.
IDENTITY_LIBRARY_CONVERSIONS_XLSX = os.path.join(BEEGFS_ROOT, "Library_conversions.xlsx")

# Cross-library roots.  Per-library writers use libN directory symlinks from
# these roots into each library's demux_nomito analysis folders; aggregate
# files stay physically under aggregate_library_analysis.
EXPECTED_POOL_METADATA = os.path.join(EXPECTED_LINES_DIR, "pool_combinations.tsv")
AUDIT_ROOT = os.path.join(AGGREGATE_ROOT, "posthoc")
HYBRID_ROOT = os.path.join(AGGREGATE_ROOT, "hybrid")
MT_FUSION_ROOT = os.path.join(AGGREGATE_ROOT, "mitochondrial")
AMBIENT_PLOT_ROOT = os.path.join(AGGREGATE_ROOT, "ambient_rna")
GEX_AMBIENT_ROOT = os.path.join(AGGREGATE_ROOT, "gex_ambient")
FIGURE_ROOT = os.path.join(AGGREGATE_ROOT, "figures")
DEFAULT_AUDIT_ROOT = AUDIT_ROOT
DEFAULT_HYBRID_ROOT = HYBRID_ROOT
DEFAULT_MT_FUSION_ROOT = MT_FUSION_ROOT

# Production NoMito panel family. CONDF/DEMUX consume the nuclear panels through
# vcf_loader_daemon shared-memory segments, while MT_FUSION reads its compact
# mitochondrial BCF and site manifest directly from BeeGFS. Keep the direct
# source paths here so the daemon-backed segments and mt panel are tied to one
# explicit production panel family.
NOMITO_PANEL_DIR = "/mnt/beegfs/home/b/vcfdownsample/Downsample_ATAC_Species_poolInformer/NoMito"
VCF_SOURCE_PATHS = {
    "interindiv_20M": os.path.join(NOMITO_PANEL_DIR, "tet.vars.downsampled_20M.bcf"),
    "interindiv_het_10M": os.path.join(NOMITO_PANEL_DIR, "tet.vars.het_10M.bcf"),
    "species_20M": os.path.join(NOMITO_PANEL_DIR, "tet.vars.species_20M.bcf"),
}
MT_PANEL_FILES = {
    "vcf": os.path.join(NOMITO_PANEL_DIR, "tet.vars.mt_fusion_ratio.bcf"),
    "vcf_index": os.path.join(NOMITO_PANEL_DIR, "tet.vars.mt_fusion_ratio.bcf.csi"),
    "site_manifest": os.path.join(NOMITO_PANEL_DIR, "tet.vars.mt_fusion_ratio.site_manifest.tsv"),
    "haplotype_groups": os.path.join(NOMITO_PANEL_DIR, "tet.vars.mt_fusion_ratio.haplotype_groups.tsv"),
    "haplotype_pairwise": os.path.join(NOMITO_PANEL_DIR, "tet.vars.mt_fusion_ratio.haplotype_pairwise.tsv"),
    "pair_audit": os.path.join(NOMITO_PANEL_DIR, "tet.vars.mt_fusion_ratio.pair_audit.tsv"),
    "sample_audit": os.path.join(NOMITO_PANEL_DIR, "tet.vars.mt_fusion_ratio.sample_audit.tsv"),
    "sites_bed": os.path.join(NOMITO_PANEL_DIR, "tet.vars.mt_fusion_ratio.sites.bed"),
}

# The trained model remains at its versioned source location.  New inference,
# refinement, reconciliation, and audit outputs no longer write there.
PLOIDY_MODEL_ROOT = os.path.join(
    PROJECT_ROOT, "ploidy_classifier", "retrain_nomito_20260814")
IDENTITY_RECONCILIATION_ROOT = os.path.join(
    AGGREGATE_ROOT, "identity_reconciliation")
TETRA_REFINE_ROOT = os.path.join(AGGREGATE_ROOT, "tetra_refine")
PLOIDY_CALLS_ROOT = os.path.join(AGGREGATE_ROOT, "ploidy")
DEFAULT_TETRA_REFINE_ROOT = TETRA_REFINE_ROOT
PLOIDY_NN_H5AD = os.path.join(BEEGFS_ROOT, "h5a5_outs", "unfiltered_normed_tetmultiome_rna.h5ad")
PLOIDY_NN_WEIGHTS = os.path.join(PLOIDY_MODEL_ROOT, "model", "ploidy_nn_weights.pt")
PLOIDY_NN_MODULE = "ploidy-inference/latest"
IDENTITY_AUDIT_ROOT = AUDIT_ROOT
ATAC_MAPPING_ROOT = "/mnt/beegfs/tetmultiome_atac/mapping_output"
REFINE_CONTAM_CONDITION = "IND_CK_RF_SX0_GATED_RFREE_PFIT"
REFINE_EXTERNAL_PLOIDY_MIN_PROB = 0.90
PROCESS_SCRIPTS_DIR = os.path.dirname(os.path.abspath(__file__))
DEMUX_OUTPUT_ROOT = None
ALLOWED_IDENTITIES = None

# The focused set is the production default.  The complete historical registry
# remains available explicitly as ``--condition-set all``.

LIB_PREFIX = "Tet_2025_Multiome-RNA_"
DEMUX_SUBDIR = "demux_nomito"

# Condf file paths (library-independent, generated by Stage 1 CONDF)
CONDF_PATHS = {
    "interindiv_20M": os.path.join(CONDF_DIR, "demux_input_20M.condf"),
    "interindiv_het_10M": os.path.join(CONDF_DIR, "demux_input_HET_10M.condf"),
    "species_20M": os.path.join(CONDF_DIR, "demux_input_species_20M.condf"),
}

# Nuclear BCF/VCF source files live in NOMITO_PANEL_DIR and are loaded into
# shared memory by vcf_loader_daemon. Pipeline stages access them via
# --shared_vcf/--shared_het_vcf/--species_shared_vcf rather than opening these
# paths directly. VCF_SOURCE_PATHS records the exact production source files.

# Shared memory segment names (from vcf_loader_daemon running on pika and char).
SHARED_VCF = {
    "interindiv_20M": "/demux_vcf_20M",
    "interindiv_het_10M": "/demux_vcf_20M_het",
    "species_20M": "/demux_vcf_20M_species",
}

# SLURM defaults
SLURM_PARTITION = "compute"
SLURM_TIME = "7-00:00:00"
# Node-local /dev/shm panels make CONDF and DEMUX daemon-bound. Only
# those two generated stages receive this nodelist; every daemon-free stage is
# still unpinned. Override with --daemon-nodes, or pass an empty value to disable
# pinning intentionally after provisioning the requested segments elsewhere.
DAEMON_NODELIST = "pika,char"
MANAGE_VCF_DAEMONS = True
VCF_DAEMON_MEM_GB = 64
VCF_DAEMON_CREATE_WAIT_SECONDS = 1800
VCF_DAEMON_CLEANUP_WAIT_SECONDS = 300
VCF_DAEMON_STATE_ROOT = os.path.join(AGGREGATE_ROOT, "vcf_daemon_state")
MANAGED_VCF_RUN_ID = None
MANAGED_VCF_READY_FILE = None

MODULES = [
    "miniforge/3",
    "htslib/1.20",
    "cellbouncer/dev",
]

ALL_LIBRARIES = list(range(1, 41))
PROFILE_RESTARTS = 8
CK_CANDIDATE_SPLIT = 0.05
CK_SX0_ABBREV = "IND_CK_RF_SX0_RFREE_PFIT"
CK_SX0_GATED_ABBREV = "IND_CK_RF_SX0_GATED_RFREE_PFIT"
CK_SX025_ABBREV = "IND_CK_RF_SX025_RFREE_PFIT"
CK_SX025_GATED_ABBREV = "IND_CK_RF_SX025_GATED_RFREE_PFIT"
CK_SX050_ABBREV = "IND_CK_RF_SX050_RFREE_PFIT"
CK_SX075_ABBREV = "IND_CK_RF_SX075_RFREE_PFIT"
CK_SX1_ABBREV = "IND_CK_RF_SX1_RFREE_PFIT"
# Backward-compatible name retained for the historical gated condition.
CK_GEOMETRY_GATED_ABBREV = CK_SX025_GATED_ABBREV
CK_GEOMETRY_GATE_VERSION = "CK_GEOMETRY_GATE_V1"
CK_GEOMETRY_BASE_STRENGTH = 0.25
CK_GEOMETRY_FALLBACK_STRENGTH = 1.0
CK_GEOMETRY_PARENT_AXIS_ALPHA_MIN = 0.80
CK_GEOMETRY_AMBIENT_ORTHOGONAL_NORM_MAX = 0.10
CK_GEOMETRY_PARENT_MASS_MIN = 0.90
LEGACY2C_MODE = "legacy2c"
LEGACY2C_MIN_SIGNAL_GAP = 0.10


# =============================================================================
# Library metadata (from Tet2025_Libraries.tsv)
# =============================================================================

LIB_INFO = {
    1:  {"batch": 1, "lib_type": "Standard",  "tets": 12, "dips": 0},
    2:  {"batch": 1, "lib_type": "Standard",  "tets": 12, "dips": 0},
    3:  {"batch": 1, "lib_type": "Standard",  "tets": 10, "dips": 0},
    4:  {"batch": 1, "lib_type": "Standard",  "tets": 10, "dips": 0},
    5:  {"batch": 1, "lib_type": "Pool",      "tets": 84, "dips": 0},
    6:  {"batch": 1, "lib_type": "Pool",      "tets": 84, "dips": 0},
    7:  {"batch": 1, "lib_type": "Diploid",   "tets": 0,  "dips": 19},
    8:  {"batch": 1, "lib_type": "Diploid",   "tets": 0,  "dips": 21},
    9:  {"batch": 2, "lib_type": "Standard",  "tets": 12, "dips": 0},
    10: {"batch": 2, "lib_type": "Standard",  "tets": 12, "dips": 0},
    11: {"batch": 2, "lib_type": "Standard",  "tets": 10, "dips": 2},
    12: {"batch": 2, "lib_type": "Diploid",   "tets": 1,  "dips": 19},
    13: {"batch": 2, "lib_type": "Standard",  "tets": 12, "dips": 0},
    14: {"batch": 2, "lib_type": "Standard",  "tets": 12, "dips": 0},
    15: {"batch": 2, "lib_type": "Pool",      "tets": 77, "dips": 15},
    16: {"batch": 2, "lib_type": "Pool",      "tets": 77, "dips": 15},
    17: {"batch": 2, "lib_type": "Standard",  "tets": 12, "dips": 0},
    18: {"batch": 2, "lib_type": "Standard",  "tets": 12, "dips": 0},
    19: {"batch": 2, "lib_type": "Standard",  "tets": 10, "dips": 2},
    20: {"batch": 2, "lib_type": "Diploid",   "tets": 1,  "dips": 19},
    21: {"batch": 2, "lib_type": "Standard",  "tets": 12, "dips": 0},
    22: {"batch": 2, "lib_type": "Standard",  "tets": 12, "dips": 0},
    23: {"batch": 2, "lib_type": "Pool",      "tets": 77, "dips": 15},
    24: {"batch": 2, "lib_type": "Pool",      "tets": 77, "dips": 15},
    25: {"batch": 3, "lib_type": "Standard",  "tets": 10, "dips": 6},
    26: {"batch": 3, "lib_type": "Standard",  "tets": 12, "dips": 0},
    27: {"batch": 3, "lib_type": "Standard",  "tets": 11, "dips": 0},
    28: {"batch": 3, "lib_type": "Standard",  "tets": 10, "dips": 1},
    29: {"batch": 3, "lib_type": "Standard",  "tets": 12, "dips": 6},
    30: {"batch": 3, "lib_type": "Diploid",   "tets": 8,  "dips": 15},
    31: {"batch": 3, "lib_type": "Pool",      "tets": 46, "dips": 0},
    32: {"batch": 3, "lib_type": "Pool",      "tets": 46, "dips": 0},
    33: {"batch": 3, "lib_type": "Standard",  "tets": 11, "dips": 6},
    34: {"batch": 3, "lib_type": "Standard",  "tets": 12, "dips": 0},
    35: {"batch": 3, "lib_type": "Standard",  "tets": 13, "dips": 0},
    36: {"batch": 3, "lib_type": "Standard",  "tets": 9,  "dips": 1},
    37: {"batch": 3, "lib_type": "Standard",  "tets": 11, "dips": 5},
    38: {"batch": 3, "lib_type": "Diploid",   "tets": 7,  "dips": 15},
    39: {"batch": 3, "lib_type": "Pool",      "tets": 46, "dips": 0},
    40: {"batch": 3, "lib_type": "Pool",      "tets": 46, "dips": 0},
}


# =============================================================================
# Condition registry: full selectable estimator roster.
# =============================================================================
#
# Each condition is a dict with:
#   abbrev:               Short name used in directory names and job names
#   mode:                 1, "1+SR", "legacy2c", 2, 3, 4, or 5
#   flags:                Extra CLI flags for tet_contam_estimate (beyond the
#                         base flags that mode implies)
#   needs_empty_drops:    False, "individual", or "species"
#   needs_species_counts: bool (only for mode 1+SR conditions)
#
# Base flags per mode (always included, not in "flags" column):
#   Mode 1:    --interindividual -X {el}
#   Mode 1+SR: --interindividual -X {el} --species_regularize -P {pm}
#              (explicit joint-evidence mode; disabled from automatic runs unless requested)
#   Mode 2:    --interspecies -X {el}
#   Mode 3:    --interindividual -X {el} --warm_start {ed_ind}
#   Mode 4:    --interspecies -X {el} --warm_start {ed_sp}
#   Mode 5:    --interspecies -X {el} --fixed_ambient {ed_sp}

CONDITIONS = [
    # --- Mode 1: Interindividual, no warm start (8 conditions) ---
    {"abbrev": "IND_BASE",    "mode": 1, "flags": "",
     "needs_empty_drops": False, "needs_species_counts": False},
    {"abbrev": "IND_LOO",     "mode": 1, "flags": "--loo",
     "needs_empty_drops": False, "needs_species_counts": False},
    {"abbrev": "IND_AD",      "mode": 1, "flags": "--adaptive_prior",
     "needs_empty_drops": False, "needs_species_counts": False},
    {"abbrev": "IND_AD_LOO",  "mode": 1, "flags": "--adaptive_prior --loo",
     "needs_empty_drops": False, "needs_species_counts": False},
    {"abbrev": "IND_RF",      "mode": 1, "flags": "--r_feedback",
     "needs_empty_drops": False, "needs_species_counts": False},
    {"abbrev": "IND_RF_LOO",  "mode": 1, "flags": "--r_feedback --loo",
     "needs_empty_drops": False, "needs_species_counts": False},
    {"abbrev": "IND_RA",      "mode": 1, "flags": "--r_feedback --adaptive_prior",
     "needs_empty_drops": False, "needs_species_counts": False},
    {"abbrev": "IND_RA_LOO",  "mode": 1, "flags": "--r_feedback --adaptive_prior --loo",
     "needs_empty_drops": False, "needs_species_counts": False},

    # --- Mode 1 + species regularization (2 conditions) ---
    {"abbrev": "IND_AD_LOO_SR",  "mode": "1+SR",
     "flags": "--adaptive_prior --loo --species_counts {sc}",
     "needs_empty_drops": False, "needs_species_counts": True},
    {"abbrev": "IND_RA_LOO_SR",  "mode": "1+SR",
     "flags": "--r_feedback --adaptive_prior --loo --species_counts {sc}",
     "needs_empty_drops": False, "needs_species_counts": True},

    # --- Legacy compiled two-component compatibility comparators (2 conditions) ---
    {"abbrev": "LEGACY2C_CLASSIC", "mode": LEGACY2C_MODE,
     "legacy2c_variant": "classic", "flags": "",
     "needs_empty_drops": False, "needs_species_counts": False},
    {"abbrev": "LEGACY2C_TET_AWARE", "mode": LEGACY2C_MODE,
     "legacy2c_variant": "tet_aware", "flags": "",
     "needs_empty_drops": False, "needs_species_counts": False},

    # --- Mode 2: Interspecies, no warm start (6 conditions) ---
    {"abbrev": "SP_BASE",     "mode": 2, "flags": "",
     "needs_empty_drops": False, "needs_species_counts": False},
    {"abbrev": "SP_LOO",      "mode": 2, "flags": "--loo",
     "needs_empty_drops": False, "needs_species_counts": False},
    {"abbrev": "SP_AD",       "mode": 2, "flags": "--adaptive_prior",
     "needs_empty_drops": False, "needs_species_counts": False},
    {"abbrev": "SP_AD_LOO",   "mode": 2, "flags": "--adaptive_prior --loo",
     "needs_empty_drops": False, "needs_species_counts": False},
    {"abbrev": "SP_RF",       "mode": 2, "flags": "--r_feedback",
     "needs_empty_drops": False, "needs_species_counts": False},
    {"abbrev": "SP_RF_LOO",   "mode": 2, "flags": "--r_feedback --loo",
     "needs_empty_drops": False, "needs_species_counts": False},
    {"abbrev": "SP_FIXED_EMPTY", "mode": 5, "flags": "",
     "needs_empty_drops": "species_fixed", "needs_species_counts": False},

    # --- Mode 3: Interindividual + warm start (4 conditions) ---
    {"abbrev": "IND_WS",          "mode": 3, "flags": "",
     "needs_empty_drops": "individual", "needs_species_counts": False},
    {"abbrev": "IND_WS_LOO",      "mode": 3, "flags": "--loo",
     "needs_empty_drops": "individual", "needs_species_counts": False},
    {"abbrev": "IND_WS_AD_LOO",   "mode": 3, "flags": "--adaptive_prior --loo",
     "needs_empty_drops": "individual", "needs_species_counts": False},
    {"abbrev": "IND_WS_RA_LOO",   "mode": 3, "flags": "--r_feedback --adaptive_prior --loo",
     "needs_empty_drops": "individual", "needs_species_counts": False},

    # --- Mode 4: Interspecies + warm start (3 conditions) ---
    {"abbrev": "SP_WS",          "mode": 4, "flags": "",
     "needs_empty_drops": "species", "needs_species_counts": False},
    {"abbrev": "SP_WS_LOO",      "mode": 4, "flags": "--loo",
     "needs_empty_drops": "species", "needs_species_counts": False},
    {"abbrev": "SP_WS_AD_LOO",   "mode": 4, "flags": "--adaptive_prior --loo",
     "needs_empty_drops": "species", "needs_species_counts": False},
]

# Complete 16-cell candidate-keyed production factorial ported from the
# frozen CK_FACTORIAL_DESIGN.  The four ordinary
# estimator mechanisms vary independently while candidate-keyed rows stay on.
CK_FACTORIAL_DESIGN = {
    "IND_CK":               {"warm_start": False, "adaptive_prior": False, "r_feedback": False, "loo": False},
    "IND_LOO_CK":           {"warm_start": False, "adaptive_prior": False, "r_feedback": False, "loo": True},
    "IND_AD_CK":            {"warm_start": False, "adaptive_prior": True,  "r_feedback": False, "loo": False},
    "IND_AD_LOO_CK":        {"warm_start": False, "adaptive_prior": True,  "r_feedback": False, "loo": True},
    "IND_RF_CK":            {"warm_start": False, "adaptive_prior": False, "r_feedback": True,  "loo": False},
    "IND_RF_LOO_CK":        {"warm_start": False, "adaptive_prior": False, "r_feedback": True,  "loo": True},
    "IND_RA_CK":            {"warm_start": False, "adaptive_prior": True,  "r_feedback": True,  "loo": False},
    "IND_RA_LOO_CK":        {"warm_start": False, "adaptive_prior": True,  "r_feedback": True,  "loo": True},
    "IND_WS_CK":            {"warm_start": True,  "adaptive_prior": False, "r_feedback": False, "loo": False},
    "IND_WS_LOO_CK":        {"warm_start": True,  "adaptive_prior": False, "r_feedback": False, "loo": True},
    "IND_WS_AD_CK":         {"warm_start": True,  "adaptive_prior": True,  "r_feedback": False, "loo": False},
    "IND_WS_AD_LOO_CK":     {"warm_start": True,  "adaptive_prior": True,  "r_feedback": False, "loo": True},
    "IND_WS_RF_CK":         {"warm_start": True,  "adaptive_prior": False, "r_feedback": True,  "loo": False},
    "IND_WS_RF_LOO_CK":     {"warm_start": True,  "adaptive_prior": False, "r_feedback": True,  "loo": True},
    "IND_WS_RA_CK":         {"warm_start": True,  "adaptive_prior": True,  "r_feedback": True,  "loo": False},
    "IND_WS_RA_LOO_CK":     {"warm_start": True,  "adaptive_prior": True,  "r_feedback": True,  "loo": True},
}


def _build_ck_condition(abbrev, spec):
    flags = []
    if spec["r_feedback"]:
        flags.append("--r_feedback")
    if spec["adaptive_prior"]:
        flags.append("--adaptive_prior")
    if spec["loo"]:
        flags.append("--loo")
    flags.extend([
        "--candidate_keyed_rows",
        "--candidate_keyed_split",
        f"{CK_CANDIDATE_SPLIT:.2f}",
    ])
    return {
        "abbrev": abbrev,
        "mode": 3 if spec["warm_start"] else 1,
        "flags": " ".join(flags),
        "needs_empty_drops": "individual" if spec["warm_start"] else False,
        "needs_species_counts": False,
        "candidate_keyed": True,
    }


def _build_final_ck_condition(abbrev, source_exclusion_strength, role):
    """Build one frozen fitted-profile/free-r CK finalist."""
    return {
        "abbrev": abbrev,
        "mode": 1,
        "flags": (
            "--r_feedback --candidate_keyed_rows "
            f"--candidate_keyed_split {CK_CANDIDATE_SPLIT:.2f} "
            f"--source_exclusion_strength {float(source_exclusion_strength):.2f}"
        ),
        "needs_empty_drops": False,
        "needs_species_counts": False,
        "candidate_keyed": True,
        "condition_role": role,
        "source_exclusion_strength": float(source_exclusion_strength),
        "final_ck_condition": True,
    }


def _build_final_ck_gated_condition(abbrev, base_strength, role):
    """Build one frozen CK_GEOMETRY_GATE_V1 finalist."""
    return {
        "abbrev": abbrev,
        "mode": 1,
        "runner": "geometry_gated",
        "flags": (
            "--r_feedback --candidate_keyed_rows "
            f"--candidate_keyed_split {CK_CANDIDATE_SPLIT:.2f}"
        ),
        "needs_empty_drops": False,
        "needs_species_counts": False,
        "candidate_keyed": True,
        "condition_role": role,
        "base_source_exclusion_strength": float(base_strength),
        "fallback_source_exclusion_strength": CK_GEOMETRY_FALLBACK_STRENGTH,
        "geometry_gate_version": CK_GEOMETRY_GATE_VERSION,
        "geometry_gate_parent_axis_alpha_min": CK_GEOMETRY_PARENT_AXIS_ALPHA_MIN,
        "geometry_gate_ambient_orthogonal_norm_max": CK_GEOMETRY_AMBIENT_ORTHOGONAL_NORM_MAX,
        "geometry_gate_parent_mass_min": CK_GEOMETRY_PARENT_MASS_MIN,
        "final_ck_condition": True,
    }


_CK_CONDITIONS_BY_ABBREV = {
    abbrev: _build_ck_condition(abbrev, spec)
    for abbrev, spec in CK_FACTORIAL_DESIGN.items()
}
_CK_NO_WS_ABBREVS = [
    abbrev for abbrev, spec in CK_FACTORIAL_DESIGN.items()
    if not spec["warm_start"]
]
_CK_WS_ABBREVS = [
    abbrev for abbrev, spec in CK_FACTORIAL_DESIGN.items()
    if spec["warm_start"]
]

# Keep job-generation order aligned with the ordinary estimator families:
# no-warm-start CK follows Mode 1 and warm-start CK follows Mode 3.
_BASE_CONDITIONS = CONDITIONS
CONDITIONS = []
for condition in _BASE_CONDITIONS:
    CONDITIONS.append(condition)
    if condition["abbrev"] == "IND_RA_LOO":
        CONDITIONS.extend(
            _CK_CONDITIONS_BY_ABBREV[name] for name in _CK_NO_WS_ABBREVS)
        CONDITIONS.append({
            "abbrev": CK_GEOMETRY_GATED_ABBREV,
            "mode": 1,
            "runner": "geometry_gated",
            "flags": (
                "--r_feedback --candidate_keyed_rows "
                f"--candidate_keyed_split {CK_CANDIDATE_SPLIT:.2f}"
            ),
            "needs_empty_drops": False,
            "needs_species_counts": False,
            "candidate_keyed": True,
            "condition_role": "validation_candidate",
            "base_source_exclusion_strength": CK_GEOMETRY_BASE_STRENGTH,
            "fallback_source_exclusion_strength": CK_GEOMETRY_FALLBACK_STRENGTH,
            "geometry_gate_version": CK_GEOMETRY_GATE_VERSION,
            "geometry_gate_parent_axis_alpha_min": CK_GEOMETRY_PARENT_AXIS_ALPHA_MIN,
            "geometry_gate_ambient_orthogonal_norm_max": CK_GEOMETRY_AMBIENT_ORTHOGONAL_NORM_MAX,
            "geometry_gate_parent_mass_min": CK_GEOMETRY_PARENT_MASS_MIN,
        })
    elif condition["abbrev"] == "IND_WS_RA_LOO":
        CONDITIONS.extend(
            _CK_CONDITIONS_BY_ABBREV[name] for name in _CK_WS_ABBREVS)

# Preserve the complete historical sequence above. Append only the six final CK
# conditions that were not already present; SX025_GATED remains the original
# historical entry and is not duplicated or redefined.
FINAL_CK_NEW_CONDITIONS = [
    _build_final_ck_condition(CK_SX0_ABBREV, 0.00, "all_source_anchor"),
    _build_final_ck_gated_condition(
        CK_SX0_GATED_ABBREV, 0.00, "lead_single_rate_candidate"),
    _build_final_ck_condition(
        CK_SX025_ABBREV, 0.25, "ungated_gate_mechanism_comparator"),
    _build_final_ck_condition(
        CK_SX050_ABBREV, 0.50, "individual_profile_comparator"),
    _build_final_ck_condition(
        CK_SX075_ABBREV, 0.75, "balanced_profile_comparator"),
    _build_final_ck_condition(
        CK_SX1_ABBREV, 1.00, "species_profile_and_gate_endpoint"),
]
CONDITIONS.extend(FINAL_CK_NEW_CONDITIONS)

# The historical SX025_GATED registry entry remains byte-for-byte unchanged.
# Its final CK role is declared only in the evaluator/reporting metadata below,
# not by rewriting the historical condition definition.

# Build lookup by abbreviation
COND_BY_ABBREV = {c["abbrev"]: c for c in CONDITIONS}
ALL_CONDITION_ABBREVS = [c["abbrev"] for c in CONDITIONS]

FINAL_CK_ABBREVS = [
    CK_SX0_ABBREV,
    CK_SX0_GATED_ABBREV,
    CK_SX025_ABBREV,
    CK_SX025_GATED_ABBREV,
    CK_SX050_ABBREV,
    CK_SX075_ABBREV,
    CK_SX1_ABBREV,
]
CK_RATE_ABBREVS = [
    CK_SX0_ABBREV,
    CK_SX0_GATED_ABBREV,
    CK_SX025_ABBREV,
    CK_SX025_GATED_ABBREV,
    CK_SX1_ABBREV,
]
CK_PROFILE_ABBREVS = [
    CK_SX0_ABBREV,
    CK_SX050_ABBREV,
    CK_SX075_ABBREV,
    CK_SX1_ABBREV,
]
CK_MINIMAL_ABBREVS = [
    CK_SX0_ABBREV,
    CK_SX0_GATED_ABBREV,
    CK_SX025_GATED_ABBREV,
]

# Efficient real-library comparison set. The complete historical registry
# remains the default. This explicit subset retains the final CK roster plus
# the historical anchors needed to interpret CK: IND_BASE; matched RF/source-
# exclusion controls; the prior best warm-start CK/non-CK pair; both compiled
# legacy two-component comparators; and fitted, RF, hard-exclusion, and fixed-
# empty native-species controls.
REAL_LIBRARY_COMPARISON_ABBREVS = [
    "IND_BASE",
    "IND_RF",
    "IND_RF_LOO",
    "IND_WS_RA_LOO",
    "IND_WS_RA_LOO_CK",
    "LEGACY2C_CLASSIC",
    "LEGACY2C_TET_AWARE",
    "SP_BASE",
    "SP_RF",
    "SP_RF_LOO",
    "SP_FIXED_EMPTY",
    *FINAL_CK_ABBREVS,
]

# Named sets are additive conveniences only. They never redefine or delete the
# cumulative full registry.
CONDITION_SETS = {
    "focused": [CK_SX0_GATED_ABBREV],
    "all": ALL_CONDITION_ABBREVS,
    "real-library-comparison": REAL_LIBRARY_COMPARISON_ABBREVS,
    "ck-finalists": FINAL_CK_ABBREVS,
    "ck-rate": CK_RATE_ABBREVS,
    "ck-profile": CK_PROFILE_ABBREVS,
    "ck-minimal": CK_MINIMAL_ABBREVS,
}
DEFAULT_CONDITION_SET = "focused"


def resolve_condition_set(name):
    """Resolve one named condition roster."""
    if name not in CONDITION_SETS:
        raise ValueError(
            f"unknown condition set {name!r}; available: "
            + ", ".join(sorted(CONDITION_SETS))
        )
    return [COND_BY_ABBREV[abbrev] for abbrev in CONDITION_SETS[name]]


# =============================================================================
# Path helpers
# =============================================================================

def get_lib_dir(lib_num):
    """Return the base directory for a library under BEEGFS_ROOT."""
    return os.path.join(BEEGFS_ROOT, f"{LIB_PREFIX}{lib_num}")


def get_demux_mapping_root():
    """Return the mapping root that owns the selected DEMUX outputs."""
    return DEMUX_OUTPUT_ROOT or BEEGFS_ROOT


def get_demux_dir(lib_num):
    """Return the demux output directory for a library.

    With --demux-output-root, inputs still come from BEEGFS_ROOT while outputs
    mirror the normal mapping layout under the selected root:
      <root>/Tet_2025_Multiome-RNA_<N>/demux_nomito/
    """
    if DEMUX_OUTPUT_ROOT:
        return os.path.join(
            DEMUX_OUTPUT_ROOT, f"{LIB_PREFIX}{lib_num}", DEMUX_SUBDIR)
    return os.path.join(get_lib_dir(lib_num), DEMUX_SUBDIR)


def get_demux_prefix(lib_num):
    """Return the demux output prefix (without extension) for a library."""
    return os.path.join(get_demux_dir(lib_num), f"lib{lib_num}_demuxed")


def get_local_condf_path(lib_num, key="interindiv_20M"):
    """Return the self-contained CONDF copy stored in a library demux directory."""
    filenames = {
        "interindiv_20M": "demux_input_20M.condf",
        "interindiv_het_10M": "demux_input_HET_10M.condf",
        "species_20M": "demux_input_species_20M.condf",
    }
    if key not in filenames:
        raise KeyError(f"unknown CONDF key: {key}")
    return os.path.join(get_demux_dir(lib_num), filenames[key])


def get_filtered_barcodes(lib_num):
    """Return path to the CellRanger filtered barcodes file for a library."""
    return os.path.join(
        get_lib_dir(lib_num), "filtered", "barcodes.tsv.gz")


def get_filtered_features(lib_num):
    """Return the CellRanger filtered-matrix feature file."""
    return os.path.join(get_lib_dir(lib_num), "filtered", "features.tsv.gz")


def get_filtered_matrix(lib_num):
    """Return the CellRanger filtered sparse count matrix."""
    return os.path.join(get_lib_dir(lib_num), "filtered", "matrix.mtx.gz")


def _gex_safe_name(value):
    """Return a stable, path-safe GEX analysis label."""
    label = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value)).strip("._")
    if not label:
        raise ValueError("GEX ambient analysis name has no safe characters")
    return label


def get_gex_cluster_path(lib_num, args):
    """Return the manual or automatically generated per-library clusters."""
    if args.gex_cluster_source in {"auto", "h5ad"}:
        marker_path = (os.path.abspath(args.gex_marker_genes)
                       if args.gex_marker_genes else None)
        marker_digest = None
        if marker_path and os.path.isfile(marker_path):
            digest = hashlib.sha1()
            with open(marker_path, "rb") as handle:
                for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                    digest.update(chunk)
            marker_digest = digest.hexdigest()
        signature_payload = {
            "source": args.gex_cluster_source,
            "h5ad": (os.path.abspath(args.ploidy_nn_h5ad)
                     if args.gex_cluster_source == "h5ad" else None),
            "column": args.gex_h5ad_cluster_column,
            "marker_genes": marker_path,
            "marker_genes_sha1": marker_digest,
            "feature_type": args.gex_feature_type,
            "data_features": args.gex_auto_data_features,
            "min_genes": args.gex_auto_min_genes,
            "neighbors": args.gex_auto_neighbors,
            "pca_components": args.gex_auto_pca_components,
            "resolutions": args.gex_auto_resolution_grid,
            "stability_repeats": args.gex_auto_stability_repeats,
            "exclude_regex": args.gex_auto_exclude_genes_regex,
            "seed": args.gex_auto_random_seed,
            "minimum_cluster_cells": args.gex_min_cluster_cells,
        }
        signature = hashlib.sha1(
            json.dumps(signature_payload, sort_keys=True).encode("utf-8")
        ).hexdigest()[:12]
        token = f"{args.gex_cluster_source}.{signature}"
        return os.path.join(
            get_gex_ambient_summary_dir(args), "clusters",
            f"lib{lib_num}.{token}.tsv")
    if args.gex_cluster_source != "manual":
        raise ValueError(f"unknown GEX cluster source: {args.gex_cluster_source}")
    if not args.gex_clusters_template:
        return None
    try:
        return os.path.abspath(args.gex_clusters_template.format(
            lib=lib_num, lib_num=lib_num, library=f"lib{lib_num}"))
    except (KeyError, ValueError) as exc:
        raise ValueError(
            "--gex-clusters-template supports only {lib}, {lib_num}, and "
            f"{{library}}: {exc}") from exc


def get_gex_auto_cluster_qc_path(lib_num, args):
    return get_gex_cluster_path(lib_num, args) + ".qc.tsv"


def get_gex_auto_cluster_contract_path(lib_num, args):
    return get_gex_cluster_path(lib_num, args) + ".contract.json"


def gex_auto_cluster_outputs_complete(lib_num, args):
    return all(check_file_exists(path) for path in (
        get_gex_cluster_path(lib_num, args),
        get_gex_auto_cluster_qc_path(lib_num, args),
        get_gex_auto_cluster_contract_path(lib_num, args),
    ))


def get_gex_ambient_run_dir(lib_num, cond_abbrev, assignment_source, args):
    """Return the isolated, source-adjacent directory for one GEX fit."""
    return os.path.join(
        get_contam_dir(lib_num, cond_abbrev, assignment_source),
        "gex_ambient", _gex_safe_name(args.gex_analysis_name))


def get_gex_ambient_prefix(lib_num, cond_abbrev, assignment_source, args):
    run_dir = get_gex_ambient_run_dir(
        lib_num, cond_abbrev, assignment_source, args)
    return os.path.join(run_dir, f"lib{lib_num}_demuxed")


def get_gex_ambient_summary_dir(args):
    return os.path.join(
        os.path.abspath(args.gex_ambient_root),
        _gex_safe_name(args.gex_analysis_name))


def gex_ambient_outputs_complete(lib_num, cond_abbrev, assignment_source, args):
    prefix = get_gex_ambient_prefix(
        lib_num, cond_abbrev, assignment_source, args)
    return all(check_file_exists(path) for path in (
        prefix + ".gex_profile",
        prefix + "_mtx/barcodes.tsv.gz",
        prefix + "_mtx/features.tsv.gz",
        prefix + "_mtx/matrix.mtx.gz",
        prefix + ".gex_run_contract.json",
        prefix + ".eligible_clusters.tsv",
        prefix + ".eligible_mtx/barcodes.tsv.gz",
        prefix + ".eligible_mtx/matrix.mtx.gz",
        prefix + ".gex_input_qc.tsv",
    ))


_EXPECTED_LINES_WORKBOOK_CACHE = None


def get_expected_lines_path(lib_num):
    """Return the canonical orchestrator-owned individual expected-lines path."""
    return os.path.join(get_lib_dir(lib_num), f"lib{lib_num}_expected_lines.txt")


def _load_expected_lines_workbook():
    """Load and cache the authoritative expected-lines workbook conversion data."""
    global _EXPECTED_LINES_WORKBOOK_CACHE
    if _EXPECTED_LINES_WORKBOOK_CACHE is not None:
        return _EXPECTED_LINES_WORKBOOK_CACHE

    if not os.path.isfile(LIBRARY_CONVERSIONS_XLSX):
        raise FileNotFoundError(
            f"Library conversions workbook missing: {LIBRARY_CONVERSIONS_XLSX}")

    try:
        import pandas as pd
    except ImportError as e:
        raise RuntimeError(
            "pandas is required to read Library_conversions.xlsx; "
            "run the orchestrator in the standard miniforge environment") from e

    try:
        convert_df = pd.read_excel(
            LIBRARY_CONVERSIONS_XLSX, sheet_name="convert")
        libs_df = pd.read_excel(
            LIBRARY_CONVERSIONS_XLSX, sheet_name="libs")
    except Exception as e:
        raise RuntimeError(
            f"failed to read expected-lines workbook {LIBRARY_CONVERSIONS_XLSX}: {e}") from e

    missing_convert = [
        column for column in ("VCF_ID", "Library")
        if column not in convert_df.columns
    ]
    missing_libs = [
        column for column in ("lib", "line")
        if column not in libs_df.columns
    ]
    if missing_convert or missing_libs:
        details = []
        if missing_convert:
            details.append(
                "convert sheet missing column(s): " + ", ".join(missing_convert))
        if missing_libs:
            details.append(
                "libs sheet missing column(s): " + ", ".join(missing_libs))
        raise ValueError(
            f"invalid expected-lines workbook {LIBRARY_CONVERSIONS_XLSX}: "
            + "; ".join(details))

    library_to_vcf = {}
    for _, row in convert_df.iterrows():
        vcf_id = row.get("VCF_ID")
        library = row.get("Library")
        if pd.isna(vcf_id) or pd.isna(library):
            continue
        library_to_vcf[str(library)] = str(vcf_id)

    _EXPECTED_LINES_WORKBOOK_CACHE = (pd, libs_df, library_to_vcf)
    return _EXPECTED_LINES_WORKBOOK_CACHE


def generate_expected_lines(lib_num):
    """Generate one library's individual expected-lines file from the workbook.

    This is the same conversion used by setup_and_submit_all_parallel_demux_het.sh:
    the ``convert`` sheet maps ``Library`` names to ``VCF_ID`` values, while the
    ``libs`` sheet supplies each library's ``line`` entries. Fusion components
    separated by ``x`` are converted independently, sorted, and joined with ``+``.
    The unique converted identities are then written in sorted order.
    """
    pd, libs_df, library_to_vcf = _load_expected_lines_workbook()
    output_file = get_expected_lines_path(lib_num)

    lib_rows = libs_df[libs_df["lib"] == lib_num]
    if lib_rows.empty:
        raise ValueError(
            f"no entries in workbook libs sheet for library {lib_num}: "
            f"{LIBRARY_CONVERSIONS_XLSX}")

    def convert_line_to_vcf(line_value):
        if pd.isna(line_value):
            return None
        line_text = str(line_value).strip()
        parts = line_text.split("x") if "x" in line_text else [line_text]
        vcf_parts = []
        for part in parts:
            part = part.strip()
            if part in library_to_vcf:
                vcf_parts.append(library_to_vcf[part])
            else:
                return None
        return "+".join(sorted(vcf_parts))

    vcf_lines = []
    for _, row in lib_rows.iterrows():
        line_value = row.get("line", None)
        if pd.isna(line_value):
            continue
        vcf_line = convert_line_to_vcf(line_value)
        if vcf_line:
            vcf_lines.append(vcf_line)

    vcf_lines = sorted(set(vcf_lines))
    if not vcf_lines:
        raise ValueError(
            f"no expected identities could be generated for library {lib_num} "
            f"from {LIBRARY_CONVERSIONS_XLSX}")

    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    payload = "".join(f"{line}\n" for line in vcf_lines)
    current = None
    try:
        with open(output_file, "r") as existing:
            current = existing.read()
    except OSError:
        pass

    if current != payload:
        staged = output_file + ".new"
        with open(staged, "w") as out:
            out.write(payload)
        os.replace(staged, output_file)

    return output_file


def get_expected_lines(lib_num, create=True):
    """Return the canonical expected-lines path, generating it by default."""
    path = get_expected_lines_path(lib_num)
    if create:
        return generate_expected_lines(lib_num)
    if not os.path.isfile(path) or os.path.getsize(path) == 0:
        return None
    return path


def get_individual_ambient_candidates_path(lib_num):
    """Return the explicit individual ambient-candidate roster path."""
    return os.path.join(
        get_demux_dir(lib_num), f"lib{lib_num}_ambient_candidates.txt")


def derive_individual_ambient_candidates(individual_expected, candidates_path):
    """Derive the per-library individual donor roster from expected lines.

    Expected-lines rows may be singlets or fusion identities such as A+B.  CK
    profile rows are keyed to individual ambient donors, so the roster contains
    the unique constituent tokens only, one per line.  The file is written next
    to the library's demux products; no temporary-directory path is used.
    """
    if not individual_expected or not os.path.isfile(individual_expected):
        raise FileNotFoundError(
            f"individual expected-lines file missing: {individual_expected}")

    candidates = set()
    with open(individual_expected, "r") as src:
        for raw in src:
            line = raw.split("#", 1)[0].strip()
            if not line:
                continue
            identity = line.split("\t", 1)[0].strip()
            for token in identity.split("+"):
                token = token.strip()
                if token:
                    candidates.add(token)

    if not candidates:
        raise ValueError(
            f"no individual ambient candidates could be derived from {individual_expected}")

    os.makedirs(os.path.dirname(candidates_path), exist_ok=True)
    payload = "".join(f"{candidate}\n" for candidate in sorted(candidates))
    current = None
    try:
        with open(candidates_path, "r") as existing:
            current = existing.read()
    except OSError:
        pass
    if current != payload:
        staged = candidates_path + ".new"
        with open(staged, "w") as out:
            out.write(payload)
        os.replace(staged, candidates_path)
    return candidates_path


def get_individual_ambient_candidates(lib_num, create=True):
    """Return/create the explicit individual ambient-candidate roster."""
    path = get_individual_ambient_candidates_path(lib_num)
    if create:
        return derive_individual_ambient_candidates(get_expected_lines(lib_num), path)
    if not os.path.isfile(path) or os.path.getsize(path) == 0:
        return None
    return path


def _identity_components(identity):
    """Return the nonempty singlet components of one receiver identity."""
    return tuple(
        token.strip() for token in str(identity or "").strip().split("+")
        if token.strip()
    )


def _canonical_identity(identity):
    """Return the CellBouncer donor identity used by controlled comparisons."""
    components = _identity_components(identity)
    if not components:
        return ""
    # Identity reconciliation uses A+A to retain post-hoc biological ploidy.
    # CellBouncer assignments encode donor genotype, not biological ploidy:
    # homotypic A+A is therefore A, while a real heterotypic A+B remains A+B.
    # Never send a repeated-donor biological-state label to the estimator.
    if len(set(components)) == 1:
        return components[0]
    return "+".join(sorted(components))


def _read_identity_lines(path):
    identities = []
    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            line = raw.split("#", 1)[0].strip()
            if not line:
                continue
            identity = line.split("\t", 1)[0].split(None, 1)[0].strip()
            if identity:
                identities.append(identity)
    return list(dict.fromkeys(identities))


def _open_text_auto(path):
    return (
        gzip.open(path, "rt", encoding="utf-8", errors="replace")
        if str(path).endswith(".gz") else
        open(path, "r", encoding="utf-8", errors="replace")
    )


def _read_assignment_identity_map(path):
    assignments = {}
    with _open_text_auto(path) as handle:
        for line_number, raw in enumerate(handle, start=1):
            parts = raw.strip().split()
            if not parts:
                continue
            if parts[0].lower() == "barcode":
                continue
            if len(parts) < 2:
                raise ValueError(
                    f"assignment line {line_number} has fewer than two fields: {path}")
            assignments[parts[0]] = parts[1]
    if not assignments:
        raise ValueError(f"no assignments found in {path}")
    return assignments


def _read_assignment_rows(path):
    """Read a four-column CellBouncer assignment file keyed by barcode."""
    rows = {}
    with _open_text_auto(path) as handle:
        for line_number, raw in enumerate(handle, start=1):
            parts = raw.strip().split()
            if not parts or parts[0].lower() == "barcode":
                continue
            if len(parts) != 4:
                raise ValueError(
                    f"assignment line {line_number} does not have four fields: {path}")
            if parts[0] in rows:
                raise ValueError(
                    f"duplicate assignment barcode {parts[0]} in {path}")
            rows[parts[0]] = tuple(parts)
    if not rows:
        raise ValueError(f"no assignments found in {path}")
    return rows


def _truthy(value):
    return str(value or "").strip().lower() in {"1", "true", "t", "yes", "y"}


def _write_if_changed(path, payload):
    """Atomically update a generated comparison input without mtime churn."""
    path = os.path.abspath(path)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    try:
        with open(path, "r", encoding="utf-8") as handle:
            if handle.read() == payload:
                return path
    except OSError:
        pass
    staged = f"{path}.new.{os.getpid()}"
    with open(staged, "w", encoding="utf-8") as handle:
        handle.write(payload)
    os.replace(staged, path)
    return path


def _identity_ambient_candidate_set(candidate_set=None):
    value = candidate_set or IDENTITY_AMBIENT_CANDIDATE_SET
    if value not in {"applied", "exploratory"}:
        raise ValueError(
            "reconciliation candidate set must be applied or exploratory: "
            f"{value}")
    return value


def get_identity_ambient_plan_dir(lib_num, candidate_set=None):
    candidate_set = _identity_ambient_candidate_set(candidate_set)
    return os.path.join(
        get_demux_dir(lib_num), IDENTITY_AMBIENT_DIRNAME,
        candidate_set, "plan")


def get_identity_ambient_context_path(lib_num, candidate_set=None):
    return os.path.join(
        get_identity_ambient_plan_dir(lib_num, candidate_set),
        f"lib{lib_num}.comparison_context.json")


def get_identity_ambient_scrutiny_path(lib_num, candidate_set=None):
    return os.path.join(
        get_identity_ambient_plan_dir(lib_num, candidate_set),
        f"lib{lib_num}.scrutiny_cells.tsv")


def get_identity_ambient_candidate_provenance_path(
        lib_num, candidate_set=None):
    return os.path.join(
        get_identity_ambient_plan_dir(lib_num, candidate_set),
        f"lib{lib_num}.ambient_candidate_provenance.tsv")


def get_identity_ambient_reconciled_assignments_path(lib_num, candidate_set=None):
    """Return the fixed-universe reconciliation overlay used by arms C/D."""
    return os.path.join(
        get_identity_ambient_plan_dir(lib_num, candidate_set),
        f"lib{lib_num}.comparison_reconciled.assignments")


def _identity_ambient_roster_path(
        lib_num, roster_basis, kind, candidate_set=None):
    if roster_basis not in {"augmented", "replacement"}:
        raise ValueError(f"generated roster basis is not supported: {roster_basis}")
    if kind not in {"receiver_lines", "ambient_candidates"}:
        raise ValueError(f"unknown identity/ambient roster kind: {kind}")
    return os.path.join(
        get_identity_ambient_plan_dir(lib_num, candidate_set),
        f"lib{lib_num}.{roster_basis}_{kind}.txt")


def _identity_ambient_file_record(path):
    path = os.path.abspath(path)
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    status = os.stat(path)
    return {
        "path": path,
        "size_bytes": int(status.st_size),
        "sha256": digest.hexdigest(),
    }


def _selected_reconciliation_candidate(row, candidate_set):
    if _truthy(row.get("reassignment_applied")):
        return True
    if candidate_set == "applied":
        return False
    if candidate_set != "exploratory":
        raise ValueError(
            f"unknown reconciliation candidate set: {candidate_set}")
    proposed = str(row.get("proposed_donor_genotype", "")).strip()
    original = str(row.get("original_demux_assignment", "")).strip()
    current = str(row.get("current_refined_assignment", "")).strip() or original
    if not proposed or _canonical_identity(proposed) in {
            _canonical_identity(original), _canonical_identity(current)}:
        return False
    if str(row.get("proposed_droplet_state", "")).strip() != "SINGLE_CELL":
        return False
    admissibility = str(
        row.get("proposed_biological_admissibility", "")).strip()
    if admissibility not in {
            "BIOLOGICAL_SINGLE_CELL_ALLOWED", "SINGLET_IDENTITY_CANDIDATE"}:
        return False
    confidence = str(row.get("decision_confidence", "")).strip()
    action = str(row.get("final_action", "")).strip()
    return (
        confidence in {"DECISIVE", "STRONG_NOT_AUTOAPPLIED"}
        and action not in {"", "KEEP"}
    )


def prepare_identity_ambient_comparison(lib_num, candidate_set=None):
    """Create the fixed rosters and barcode strata for the four-arm analysis.

    The plan is derived from validated reconciliation decisions plus the actual
    reconciled assignment export. It never changes the demux count/CONDF bundle.
    Arm D is marked not applicable unless an affected original identity is fully
    displaced and removing it changes the ambient donor roster.
    """
    candidate_set = _identity_ambient_candidate_set(candidate_set)
    cache_key = (int(lib_num), candidate_set)
    if cache_key in _IDENTITY_AMBIENT_PLAN_CACHE:
        return dict(_IDENTITY_AMBIENT_PLAN_CACHE[cache_key])

    expected_path = get_expected_lines(lib_num)
    original_candidates_path = get_individual_ambient_candidates(lib_num)
    demux_assignments_path = get_demux_prefix(lib_num) + ".assignments"
    demux_samples_path = get_demux_prefix(lib_num) + ".samples"
    reconciled_assignments_path = get_reconciled_assignments_path(lib_num)
    decisions_path = get_reconciled_cells_path(lib_num)
    required = (
        expected_path, original_candidates_path, demux_assignments_path,
        reconciled_assignments_path, decisions_path, demux_samples_path,
    )
    missing = [path for path in required if not path or not check_file_exists(path)]
    if missing:
        raise FileNotFoundError(
            "four-arm comparison requires completed demux and validated identity "
            "reconciliation inputs; missing: " + ", ".join(str(x) for x in missing))

    input_records = {
        "original_receiver_lines": _identity_ambient_file_record(expected_path),
        "original_ambient_candidates": _identity_ambient_file_record(
            original_candidates_path),
        "demux_assignments": _identity_ambient_file_record(
            demux_assignments_path),
        "validated_reconciled_assignments": _identity_ambient_file_record(
            reconciled_assignments_path),
        "reconciled_cells": _identity_ambient_file_record(decisions_path),
        "demux_samples": _identity_ambient_file_record(demux_samples_path),
    }
    plan_fingerprint_payload = {
        "plan_version": IDENTITY_AMBIENT_PLAN_VERSION,
        "candidate_set": candidate_set,
        "inputs": input_records,
    }
    plan_fingerprint = hashlib.sha256(json.dumps(
        plan_fingerprint_payload, sort_keys=True,
        separators=(",", ":")).encode("utf-8")).hexdigest()

    original_receivers = _read_identity_lines(expected_path)
    original_donors = set(_read_identity_lines(original_candidates_path))
    demux_rows = _read_assignment_rows(demux_assignments_path)
    reconciled_rows = _read_assignment_rows(reconciled_assignments_path)
    demux_assignments = {
        barcode: row[1] for barcode, row in demux_rows.items()}
    validated_reconciled_assignments = {
        barcode: row[1] for barcode, row in reconciled_rows.items()}
    unexpected_reconciled = sorted(set(reconciled_rows) - set(demux_rows))
    if unexpected_reconciled:
        raise ValueError(
            "validated reconciliation contains barcodes absent from the demux "
            "assignment universe (first 10): "
            + ", ".join(unexpected_reconciled[:10]))

    # Keep the cell universe identical in A/B/C/D. The validated compatibility
    # export intentionally excludes technical multiplets, so using it directly
    # would confound assignment changes with barcode removal. Arms C/D instead
    # overlay every validated single-cell identity onto the original demux rows
    # and preserve all other rows unchanged.
    reconciled_assignments = dict(demux_assignments)
    reconciled_assignments.update(validated_reconciled_assignments)
    comparison_assignments_path = (
        get_identity_ambient_reconciled_assignments_path(
            lib_num, candidate_set))
    comparison_lines = []
    for barcode, demux_row in demux_rows.items():
        reconciled_row = reconciled_rows.get(barcode)
        # Reconciliation may describe a homotetraploid as A+A. That notation
        # remains in the reconciliation records, but the CellBouncer comparison
        # ledger must carry its donor identity A and type S.
        biological_identity = (
            reconciled_row[1] if reconciled_row is not None else demux_row[1])
        identity = _canonical_identity(biological_identity)
        if not identity:
            raise ValueError(
                f"empty reconciled genotype identity for barcode {barcode}")
        assignment_type = "D" if "+" in identity else "S"
        # A/B/C/D are a controlled identity/roster experiment. CellBouncer
        # uses assignment LLRs as statistical weights, so carrying a changed
        # reconciliation score into C/D would confound the identity contrast.
        # Preserve the original demux score for every barcode in every arm.
        score = demux_row[3]
        comparison_lines.append(
            "\t".join((barcode, identity, assignment_type, score)))
    _write_if_changed(
        comparison_assignments_path, "\n".join(comparison_lines) + "\n")

    with _open_text_auto(decisions_path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required_fields = {
            "barcode", "original_demux_assignment", "proposed_donor_genotype",
            "reconciled_donor_genotype", "proposed_droplet_state",
            "proposed_biological_admissibility", "final_action",
            "decision_confidence", "reassignment_applied",
            "explicit_multiplet_evidence", "occupancy_resolution_status",
            "library_exchange_evidence_eligible", "event_id",
        }
        missing_fields = required_fields - set(reader.fieldnames or [])
        if missing_fields:
            raise ValueError(
                f"reconciliation decision schema missing {sorted(missing_fields)}: "
                f"{decisions_path}")
        decision_rows = {}
        for row in reader:
            barcode = str(row.get("barcode", "")).strip()
            if not barcode:
                raise ValueError(
                    f"reconciliation decision row has an empty barcode: {decisions_path}")
            if barcode in decision_rows:
                raise ValueError(
                    f"duplicate reconciliation decision barcode {barcode}: "
                    f"{decisions_path}")
            decision_rows[barcode] = row

    decision_only_barcodes = sorted(set(decision_rows) - set(demux_assignments))
    if decision_only_barcodes:
        raise ValueError(
            "reconciliation decisions contain barcodes absent from the demux "
            "assignment universe (first 10): "
            + ", ".join(decision_only_barcodes[:10]))
    for barcode, row in decision_rows.items():
        recorded_original = str(
            row.get("original_demux_assignment", "")).strip()
        if (not recorded_original or
                _canonical_identity(recorded_original) !=
                _canonical_identity(demux_assignments[barcode])):
            raise ValueError(
                f"reconciliation decision/demux original identity mismatch for "
                f"{barcode}: {recorded_original!r} != "
                f"{demux_assignments[barcode]!r}")
    missing_decisions = sorted(set(validated_reconciled_assignments) - set(decision_rows))
    if missing_decisions:
        raise ValueError(
            "validated reconciliation assignments lack decision rows (first 10): "
            + ", ".join(missing_decisions[:10]))
    for barcode, identity in validated_reconciled_assignments.items():
        recorded_final = str(
            decision_rows[barcode].get("reconciled_donor_genotype", "")).strip()
        if (_canonical_identity(recorded_final) !=
                _canonical_identity(identity)):
            raise ValueError(
                f"reconciliation decision/export final identity mismatch for "
                f"{barcode}: {recorded_final!r} != {identity!r}")
    for barcode, row in decision_rows.items():
        if not _truthy(row.get("reassignment_applied")):
            continue
        if barcode not in validated_reconciled_assignments:
            raise ValueError(
                f"applied reconciliation decision is absent from the validated "
                f"single-cell assignment export: {barcode}")

    selected_rows = {}
    selected_source_identities = set()
    selected_candidates = set()
    for barcode, row in decision_rows.items():
        if not barcode or not _selected_reconciliation_candidate(row, candidate_set):
            continue
        original = (
            str(row.get("original_demux_assignment", "")).strip()
            or demux_assignments.get(barcode, "")
        )
        proposed = str(row.get("proposed_donor_genotype", "")).strip()
        reconciled = (
            str(row.get("reconciled_donor_genotype", "")).strip()
            or reconciled_assignments.get(barcode, "")
        )
        selected_rows[barcode] = row
        if original:
            selected_source_identities.add(_canonical_identity(original))
        candidate_identities = (
            (reconciled,) if candidate_set == "applied" else
            (proposed, reconciled))
        for identity in candidate_identities:
            if identity:
                selected_candidates.add(identity)

    # Every identity actually present after reconciliation must propagate even
    # if its decision row is not selected by the candidate filter.
    final_receivers = list(dict.fromkeys(reconciled_assignments.values()))
    final_receiver_keys = {_canonical_identity(x) for x in final_receivers}
    original_receiver_keys = {_canonical_identity(x) for x in original_receivers}
    actual_new_receivers = [
        identity for identity in final_receivers
        if _canonical_identity(identity) not in original_receiver_keys
    ]
    selected_candidates.update(actual_new_receivers)

    # Fixed groups: all cells from an affected source population, all cells
    # finally assigned to candidate identities, and every selected proposal row.
    candidate_keys = {_canonical_identity(x) for x in selected_candidates}
    scrutiny = set(selected_rows)
    scrutiny.update(
        barcode for barcode, identity in demux_assignments.items()
        if _canonical_identity(identity) in selected_source_identities)
    scrutiny.update(
        barcode for barcode, identity in reconciled_assignments.items()
        if _canonical_identity(identity) in candidate_keys)

    # Full displacement is deliberately strict: every original cell must be in
    # the reconciled export, none may retain the source identity, and all must
    # converge on one replacement identity.
    displaced = {}
    for source_key in sorted(selected_source_identities):
        source_barcodes = [
            barcode for barcode, identity in demux_assignments.items()
            if _canonical_identity(identity) == source_key
        ]
        if not source_barcodes or any(
                barcode not in reconciled_assignments for barcode in source_barcodes):
            continue
        replacements = {
            _canonical_identity(reconciled_assignments[barcode])
            for barcode in source_barcodes
        }
        if source_key not in replacements and len(replacements) == 1:
            displaced[source_key] = next(iter(replacements))

    augmented_receivers = list(original_receivers)
    for identity in final_receivers:
        if _canonical_identity(identity) not in {
                _canonical_identity(x) for x in augmented_receivers}:
            augmented_receivers.append(identity)
    augmented_donors = set(original_donors)
    for identity in (*final_receivers, *selected_candidates):
        augmented_donors.update(_identity_components(identity))

    # Production safety invariant: the applied roster may gain a donor only
    # through an applied reconciliation row or an identity actually present in
    # the validated final assignment export. Cell-level alternative candidates
    # are not library-level evidence that a donor was physically present.
    applied_added_donors = set()
    for row in selected_rows.values():
        if not _truthy(row.get("reassignment_applied")):
            continue
        applied_added_donors.update(_identity_components(
            str(row.get("reconciled_donor_genotype", "")).strip()))
    final_receiver_donors = set()
    for identity in final_receivers:
        final_receiver_donors.update(_identity_components(identity))
    if candidate_set == "applied":
        ungrounded = sorted(
            donor for donor in augmented_donors - original_donors
            if donor not in applied_added_donors and
            donor not in final_receiver_donors)
        if ungrounded:
            raise ValueError(
                f"lib{lib_num} applied ambient roster contains donors without "
                "an applied reconciliation or validated final receiver: "
                + ", ".join(ungrounded))

    replacement_receivers = [
        identity for identity in augmented_receivers
        if _canonical_identity(identity) not in displaced
    ]
    retained_receiver_donors = set()
    for identity in replacement_receivers:
        retained_receiver_donors.update(_identity_components(identity))
    displaced_source_components = set()
    for source_identity in displaced:
        displaced_source_components.update(_identity_components(source_identity))
    removable_displaced_components = (
        displaced_source_components - retained_receiver_donors)
    # Arm D removes only donor components made obsolete by a fully displaced
    # source. Keep selected-but-not-final comparison candidates so C->D is not
    # confounded by unrelated candidate removal.
    replacement_donors = (
        set(augmented_donors) - removable_displaced_components)
    replacement_eligible = bool(removable_displaced_components)

    sample_names = set(_read_identity_lines(get_demux_prefix(lib_num) + ".samples"))
    unknown = sorted(augmented_donors - sample_names)
    if unknown:
        raise ValueError(
            f"lib{lib_num} reconciliation ambient candidates absent from the "
            f"demux/VCF sample universe: {', '.join(unknown)}")

    plan_dir = get_identity_ambient_plan_dir(lib_num, candidate_set)
    augmented_receiver_path = _identity_ambient_roster_path(
        lib_num, "augmented", "receiver_lines", candidate_set)
    augmented_candidates_path = _identity_ambient_roster_path(
        lib_num, "augmented", "ambient_candidates", candidate_set)
    replacement_receiver_path = _identity_ambient_roster_path(
        lib_num, "replacement", "receiver_lines", candidate_set)
    replacement_candidates_path = _identity_ambient_roster_path(
        lib_num, "replacement", "ambient_candidates", candidate_set)
    _write_if_changed(
        augmented_receiver_path,
        "".join(f"{identity}\n" for identity in augmented_receivers))
    _write_if_changed(
        augmented_candidates_path,
        "".join(f"{identity}\n" for identity in sorted(augmented_donors)))
    if replacement_eligible:
        _write_if_changed(
            replacement_receiver_path,
            "".join(f"{identity}\n" for identity in replacement_receivers))
        _write_if_changed(
            replacement_candidates_path,
            "".join(f"{identity}\n" for identity in sorted(replacement_donors)))

    # Human-auditable donor provenance. This deliberately records roster-level
    # evidence rather than merely listing the candidates emitted to the solver.
    # In applied mode every added donor must show applied or validated-final
    # support; exploratory-only alternatives are labelled explicitly.
    provenance_evidence = {}
    for barcode, row in selected_rows.items():
        provenance_fields = (
            ("reconciled_donor_genotype",)
            if candidate_set == "applied" else
            ("proposed_donor_genotype", "reconciled_donor_genotype"))
        identities = {
            str(row.get(field, "")).strip()
            for field in provenance_fields
            if str(row.get(field, "")).strip()
        }
        components = set()
        for identity in identities:
            components.update(_identity_components(identity))
        for donor in components:
            evidence = provenance_evidence.setdefault(donor, {
                "barcodes": set(), "applied_barcodes": set(),
                "unapplied_barcodes": set(), "confidences": set(),
                "event_ids": set(), "source_identities": set(),
            })
            evidence["barcodes"].add(barcode)
            if _truthy(row.get("reassignment_applied")):
                evidence["applied_barcodes"].add(barcode)
            else:
                evidence["unapplied_barcodes"].add(barcode)
            confidence = str(row.get("decision_confidence", "")).strip()
            if confidence:
                evidence["confidences"].add(confidence)
            event_id = str(row.get("event_id", "")).strip()
            if event_id:
                evidence["event_ids"].add(event_id)
            evidence["source_identities"].update(identities)

    final_receiver_cell_counts = {donor: 0 for donor in augmented_donors}
    for identity in reconciled_assignments.values():
        for donor in _identity_components(identity):
            if donor in final_receiver_cell_counts:
                final_receiver_cell_counts[donor] += 1

    provenance_path = get_identity_ambient_candidate_provenance_path(
        lib_num, candidate_set)
    provenance_lines = [
        "donor\troster_status\tinclusion_reason\tselected_rows"
        "\tapplied_rows\tunapplied_exploratory_rows\tfinal_receiver_cells"
        "\tdecision_confidences\tevent_ids\tsource_identities"
    ]
    for donor in sorted(augmented_donors):
        evidence = provenance_evidence.get(donor, {})
        applied_count = len(evidence.get("applied_barcodes", ()))
        unapplied_count = len(evidence.get("unapplied_barcodes", ()))
        final_count = final_receiver_cell_counts.get(donor, 0)
        if donor in original_donors:
            roster_status = "original"
            inclusion_reason = "original_library_ambient_candidate"
        elif applied_count:
            roster_status = "added"
            inclusion_reason = "applied_reconciliation"
        elif final_count:
            roster_status = "added"
            inclusion_reason = "validated_final_receiver"
        else:
            roster_status = "added_exploratory"
            inclusion_reason = "unapplied_exploratory_candidate"
        provenance_lines.append("\t".join((
            donor,
            roster_status,
            inclusion_reason,
            str(len(evidence.get("barcodes", ()))),
            str(applied_count),
            str(unapplied_count),
            str(final_count),
            ";".join(sorted(evidence.get("confidences", ()))) or "NA",
            ";".join(sorted(evidence.get("event_ids", ()))) or "NA",
            ";".join(sorted(evidence.get("source_identities", ()))) or "NA",
        )))
    _write_if_changed(provenance_path, "\n".join(provenance_lines) + "\n")

    scrutiny_rows = []
    all_barcodes = sorted(set(demux_assignments) | set(reconciled_assignments))
    for barcode in all_barcodes:
        demux_identity = demux_assignments.get(barcode, "")
        reconciled_identity = reconciled_assignments.get(barcode, "")
        changed = bool(
            demux_identity and reconciled_identity and
            _canonical_identity(demux_identity) != _canonical_identity(reconciled_identity))
        selected = barcode in selected_rows
        scrutinized = barcode in scrutiny
        stratum = (
            "changed_target" if changed else
            "scrutinized_other" if scrutinized else
            "background"
        )
        row = selected_rows.get(barcode, {})
        scrutiny_rows.append({
            "barcode": barcode,
            "demux_identity": demux_identity or "NA",
            "reconciled_identity": reconciled_identity or "NA",
            "proposed_identity": str(
                row.get("proposed_donor_genotype", "")).strip() or "NA",
            "selected_candidate_row": int(selected),
            "reassignment_applied": int(
                _truthy(row.get("reassignment_applied"))),
            "changed": int(changed),
            "scrutinized": int(scrutinized),
            "stratum": stratum,
            "transition": f"{demux_identity or 'NA'} -> {reconciled_identity or 'NA'}",
            "event_id": str(row.get("event_id", "")).strip() or "NA",
            "decision_confidence": str(
                row.get("decision_confidence", "")).strip() or "NA",
        })
    scrutiny_path = get_identity_ambient_scrutiny_path(lib_num, candidate_set)
    output = [
        "barcode\tdemux_identity\treconciled_identity\tproposed_identity"
        "\tselected_candidate_row\treassignment_applied\tchanged\tscrutinized"
        "\tstratum\ttransition\tevent_id\tdecision_confidence"
    ]
    output.extend("\t".join(str(row[key]) for key in (
        "barcode", "demux_identity", "reconciled_identity", "proposed_identity",
        "selected_candidate_row", "reassignment_applied", "changed", "scrutinized",
        "stratum", "transition", "event_id", "decision_confidence"))
        for row in scrutiny_rows)
    _write_if_changed(scrutiny_path, "\n".join(output) + "\n")

    plan_artifact_paths = {
        "reconciled_assignments": comparison_assignments_path,
        "scrutiny_cells": scrutiny_path,
        "ambient_candidate_provenance": provenance_path,
        "augmented_receiver_lines": augmented_receiver_path,
        "augmented_ambient_candidates": augmented_candidates_path,
    }
    if replacement_eligible:
        plan_artifact_paths.update({
            "replacement_receiver_lines": replacement_receiver_path,
            "replacement_ambient_candidates": replacement_candidates_path,
        })
    plan_artifact_records = {
        name: _identity_ambient_file_record(path)
        for name, path in plan_artifact_paths.items()
    }

    context = {
        "schema_version": 2,
        "plan_version": IDENTITY_AMBIENT_PLAN_VERSION,
        "plan_fingerprint": plan_fingerprint,
        "input_records": input_records,
        "plan_artifact_records": plan_artifact_records,
        "library": int(lib_num),
        "candidate_set": candidate_set,
        "demux_assignments": demux_assignments_path,
        "validated_reconciled_assignments": reconciled_assignments_path,
        "reconciled_assignments": comparison_assignments_path,
        "reconciled_cells": decisions_path,
        "scrutiny_cells": scrutiny_path,
        "ambient_candidate_provenance": provenance_path,
        "original_receiver_lines": expected_path,
        "original_ambient_candidates": original_candidates_path,
        "augmented_receiver_lines": augmented_receiver_path,
        "augmented_ambient_candidates": augmented_candidates_path,
        "replacement_receiver_lines": (
            replacement_receiver_path if replacement_eligible else ""),
        "replacement_ambient_candidates": (
            replacement_candidates_path if replacement_eligible else ""),
        "n_original_receivers": len(original_receivers),
        "n_original_candidates": len(original_donors),
        "n_augmented_receivers": len(augmented_receivers),
        "n_augmented_candidates": len(augmented_donors),
        "n_replacement_receivers": (
            len(replacement_receivers) if replacement_eligible else 0),
        "n_replacement_candidates": (
            len(replacement_donors) if replacement_eligible else 0),
        "n_scrutinized": len(scrutiny),
        "n_changed": sum(int(row["changed"]) for row in scrutiny_rows),
        "n_comparison_barcodes": len(reconciled_assignments),
        "n_validated_single_cell_barcodes": len(
            validated_reconciled_assignments),
        "n_preserved_non_single_cell_barcodes": (
            len(reconciled_assignments) -
            len(validated_reconciled_assignments)),
        "assignment_score_basis": "original_demux_all_arms",
        "n_background": sum(
            row["stratum"] == "background" for row in scrutiny_rows),
        "selected_candidate_identities": sorted(
            selected_candidates, key=_canonical_identity),
        "added_candidate_components": sorted(augmented_donors - original_donors),
        "fully_displaced_identity_map": displaced,
        "removed_candidate_components": sorted(
            removable_displaced_components),
        "replacement_arm_eligible": replacement_eligible,
        "replacement_arm_skip_reason": (
            "" if replacement_eligible else
            "no fully displaced source changed the ambient donor roster"),
    }
    context_path = get_identity_ambient_context_path(lib_num, candidate_set)
    _write_if_changed(
        context_path, json.dumps(context, indent=2, sort_keys=True) + "\n")
    context["context_path"] = context_path
    context["plan_dir"] = plan_dir
    _IDENTITY_AMBIENT_PLAN_CACHE[cache_key] = dict(context)
    return dict(context)


def identity_ambient_arm_inputs(lib_num, assignment_source):
    """Return concrete assignment/receiver/donor inputs for one arm."""
    if assignment_source not in IDENTITY_AMBIENT_ARMS:
        raise ValueError(f"not a reconciliation four-arm source: {assignment_source}")
    context = prepare_identity_ambient_comparison(lib_num)
    spec = IDENTITY_AMBIENT_ARMS[assignment_source]
    roster_basis = spec["roster_basis"]
    if roster_basis == "original":
        receiver = context["original_receiver_lines"]
        candidates = context["original_ambient_candidates"]
    elif roster_basis == "augmented":
        # Arm B changes only ambient donors; Arm C also admits actual final
        # receiver identities. This preserves the planned factorial contrast.
        receiver = (
            context["original_receiver_lines"]
            if assignment_source == "demux_augmented" else
            context["augmented_receiver_lines"]
        )
        candidates = context["augmented_ambient_candidates"]
    else:
        if not context["replacement_arm_eligible"]:
            return None
        receiver = context["replacement_receiver_lines"]
        candidates = context["replacement_ambient_candidates"]
    return {
        **spec,
        "assignment_source": assignment_source,
        "assignments": (
            context["demux_assignments"]
            if spec["assignment_basis"] == "demux" else
            context["reconciled_assignments"]),
        "receiver_lines": receiver,
        "ambient_candidates": candidates,
        "scrutiny_cells": context["scrutiny_cells"],
        "context_path": context["context_path"],
        "candidate_set": context["candidate_set"],
        "plan_fingerprint": context["plan_fingerprint"],
        "replacement_arm_eligible": context["replacement_arm_eligible"],
        "replacement_arm_skip_reason": context["replacement_arm_skip_reason"],
    }


def load_individual_to_species_map():
    """Load panel metadata as individual -> species and species-name set."""
    if not os.path.isfile(PANEL_METADATA):
        raise FileNotFoundError(f"panel_metadata not found: {PANEL_METADATA}")
    mapping = {}
    species = set()
    with open(PANEL_METADATA, "r") as f:
        header = f.readline().rstrip("\n").split("\t")
        try:
            indiv_col = header.index("indiv_id")
            species_col = header.index("species")
        except ValueError:
            # Fallback for two-column files without the exact header names.
            indiv_col, species_col = 0, 1
            if len(header) >= 2 and header[0] != "indiv_id":
                mapping[header[indiv_col]] = header[species_col]
                species.add(header[species_col])
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) <= max(indiv_col, species_col):
                continue
            indiv = parts[indiv_col].strip()
            sp = parts[species_col].strip()
            if indiv and sp:
                mapping[indiv] = sp
                species.add(sp)
    return mapping, species


def collapse_expected_lines_to_species(individual_expected, species_expected):
    """Create a species-native expected-lines file from individual expected lines.

    IMPORTANT: this is a compatibility generator for the current native species
    estimator, which accepts one- or two-component species identities.  It must
    never write the pseudo-species label Hy into the native species expected-lines
    file.  Hy / Chinobo identities are folded into native B/C evidence.

    Current compatibility behavior:
      Hy / Chinobo-mCherry        -> B+C
      Hy+O / Chinobo-mCherry+O    -> B+O and C+O
      Hy+B                        -> B and B+C
      Hy+C                        -> C and B+C
      Hy+H                        -> B+H and C+H

    The function also emits native component singlets as ambient-source entries
    so the cell-based species ambient profile is allowed to estimate every native
    species component present in expected identities.  This is a temporary bridge
    until the C++ species-native estimator supports weighted multi-species
    endogenous composition vectors, e.g. Hy+O = {B:0.25, C:0.25, O:0.5}.
    """
    indiv_to_species, species_names = load_individual_to_species_map()

    native_species = {"B", "C", "H", "O"}
    hybrid_species_labels = {"Hy", "HY", "Hybrid", "hybrid", "Chinobo", "chinobo"}

    def is_chinobo_label(x):
        xl = str(x).lower()
        return "chinobo" in xl or xl in {"hy", "hybrid"}

    def token_species_options(tok):
        """Return native species options for one diploid token.

        Normal diploid identities return one native species.  The Chinobo/Hy
        hybrid returns two native options, B and C, because it carries one bonobo
        and one chimp haploid genome.
        """
        tok = tok.strip()
        if not tok:
            return set()
        if tok in native_species:
            return {tok}
        if tok in hybrid_species_labels or is_chinobo_label(tok):
            return {"B", "C"}
        if tok in species_names:
            # If panel_metadata itself uses Hy as the species, fold it here.
            if tok in hybrid_species_labels or is_chinobo_label(tok):
                return {"B", "C"}
            if tok in native_species:
                return {tok}
            raise ValueError(
                f"Species token '{tok}' is not a native species label B/C/H/O; "
                f"native interspecies mode cannot use pseudo-species labels")
        if tok in indiv_to_species:
            sp = indiv_to_species[tok]
            if sp in native_species:
                return {sp}
            if sp in hybrid_species_labels or is_chinobo_label(sp) or is_chinobo_label(tok):
                return {"B", "C"}
            raise ValueError(
                f"Cannot map individual '{tok}' with species '{sp}' to native B/C/H/O "
                f"using {PANEL_METADATA}; source file: {individual_expected}")
        raise ValueError(
            f"Cannot map expected-lines token '{tok}' to a native species composition using "
            f"{PANEL_METADATA}; source file: {individual_expected}")

    def pair_label(a, b):
        if a == b:
            return a
        x, y = sorted((a, b))
        return f"{x}+{y}"

    out_entries = set()
    ambient_source_species = set()

    with open(individual_expected, "r") as f:
        for raw in f:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            toks = [x.strip() for x in line.split("+") if x.strip()]
            if not toks:
                continue
            option_sets = [token_species_options(tok) for tok in toks]
            for opts in option_sets:
                ambient_source_species.update(opts)

            if len(option_sets) == 1:
                opts = sorted(option_sets[0])
                if len(opts) == 1:
                    out_entries.add(opts[0])
                else:
                    # Hy / Chinobo-mCherry alone: represent as B+C for the
                    # current two-component species estimator.
                    out_entries.add(pair_label(opts[0], opts[1]))
            elif len(option_sets) == 2:
                # Current estimator compatibility: if either side is Hy/B+C,
                # project to the valid pairwise native combinations.  The final
                # biological model must replace this with weighted vectors.
                for a in sorted(option_sets[0]):
                    for b in sorted(option_sets[1]):
                        out_entries.add(pair_label(a, b))
            else:
                # Do not silently invent a B+C+O-style line here because the
                # current C++ species contamination estimator is not known to
                # support arbitrary weighted multi-species endogenous identities.
                # Project all component combinations pairwise as a conservative
                # compatibility fallback.
                flat = sorted(set().union(*option_sets))
                for i, a in enumerate(flat):
                    out_entries.add(a)
                    for b in flat[i + 1:]:
                        out_entries.add(pair_label(a, b))

    # Critical for ambient estimation: include every native species component as
    # an explicit source species.  Without this, B+C expected identities can be
    # legal while B is still absent from the ambient-source universe.
    out_entries.update(ambient_source_species)

    # Refuse to write pseudo-species labels to native expected-lines.
    bad = sorted(x for x in out_entries if "Hy" in x or "Chinobo" in x)
    if bad:
        raise ValueError(f"Refusing to write non-native species labels to {species_expected}: {bad}")

    os.makedirs(os.path.dirname(species_expected), exist_ok=True)
    with open(species_expected, "w") as out:
        out.write("# species-native expected lines generated by orchestrate_tetraploid.py\n")
        out.write(f"# source: {individual_expected}\n")
        out.write("# Hy / Chinobo identities folded to native B/C species components.\n")
        out.write("# Expected-lines remain pairwise/source-compatible for candidate restriction;\n")
        out.write("# tet_contam_estimate uses .assignments + -P panel metadata for weighted hybrid cell composition.\n")
        for entry in sorted(out_entries):
            out.write(entry + "\n")
    return species_expected


def get_species_expected_lines(lib_num, create=True):
    """Return species-native expected lines path, creating it from individual expected lines if needed."""
    indiv_expected = get_expected_lines(lib_num)
    if indiv_expected is None:
        return None
    species_expected = os.path.join(get_demux_dir(lib_num), f"lib{lib_num}_expected_species_lines.txt")
    if create:
        try:
            collapse_expected_lines_to_species(indiv_expected, species_expected)
        except Exception as e:
            raise RuntimeError(f"failed to create species expected-lines for lib{lib_num}: {e}") from e
    elif not os.path.isfile(species_expected):
        return None
    return species_expected


def get_bam_path(lib_num):
    """Return path to the BAM file for a library.

    Checks multiple known locations in order of preference.
    """
    lib_dir = get_lib_dir(lib_num)
    candidates = [
        os.path.join(lib_dir, "gex.bam"),
        os.path.join(lib_dir, "Aligned.sortedByCoord.out.bam"),
        os.path.join(lib_dir, "outs", "possorted_bam.bam"),
    ]
    for path in candidates:
        if os.path.isfile(path):
            return path
    return None


def _identity_ambient_contam_dir(
        lib_num, cond_abbrev, candidate_set, assignment_source):
    return os.path.join(
        get_demux_dir(lib_num), "contamination", cond_abbrev,
        IDENTITY_AMBIENT_DIRNAME, str(candidate_set), assignment_source)


def get_contam_dir(lib_num, cond_abbrev, assignment_source="demux"):
    """Return the output directory for a specific condition/library pair."""
    base = os.path.join(get_demux_dir(lib_num), "contamination", cond_abbrev)
    if assignment_source == "demux":
        return base
    if assignment_source == "reconciled":
        return os.path.join(base, "reconciled")
    if assignment_source in IDENTITY_AMBIENT_ARMS:
        return _identity_ambient_contam_dir(
            lib_num, cond_abbrev, _identity_ambient_candidate_set(),
            assignment_source)
    raise ValueError(f"unknown contamination assignment source: {assignment_source}")


def get_contam_prefix(lib_num, cond_abbrev, assignment_source="demux"):
    """Return the output prefix for a specific condition/library pair."""
    return os.path.join(get_contam_dir(lib_num, cond_abbrev, assignment_source),
                        f"lib{lib_num}_demuxed")


def _identity_ambient_required_output_suffixes(
        cond_abbrev, include_assignment_sidecar=True):
    required = [
        ".contam_rate", ".contam_prof", ".allele_ratio",
        ".contam_diagnostics.tsv", ".cell_source_profile.tsv",
        ".profile_fit_diagnostics.tsv", ".condf_coverage.tsv",
        ".run_contract.json", ".decontam.assignments",
    ]
    if include_assignment_sidecar:
        required.append(".assignments")
    cond = COND_BY_ABBREV.get(cond_abbrev)
    if cond and cond.get("runner") == "geometry_gated":
        required.append(".geometry_gate_audit.tsv")
    return tuple(required)


def _reusable_legacy_arm_a_prefix(lib_num, cond_abbrev):
    """Return a complete pre-applied Arm A prefix, if one exists.

    Arm A uses original demux assignments and the original ambient roster, so
    its numerical inputs are independent of candidate-set selection. Older
    releases nevertheless stored it below the candidate-set directory. This
    read-only resolver lets the applied plan reuse that completed calculation.
    """
    current_inputs = identity_ambient_arm_inputs(lib_num, "demux_original")
    required = _identity_ambient_required_output_suffixes(
        cond_abbrev, include_assignment_sidecar=False)
    for legacy_set in ("supported", "exploratory"):
        prefix = os.path.join(
            _identity_ambient_contam_dir(
                lib_num, cond_abbrev, legacy_set, "demux_original"),
            f"lib{lib_num}_demuxed")
        if not all(check_file_exists(prefix + suffix) for suffix in required):
            continue
        contract_path = prefix + ".identity_ambient_arm.tsv"
        if not check_file_exists(contract_path):
            continue
        try:
            with open(contract_path, "r", encoding="utf-8", newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
        except (OSError, csv.Error):
            continue
        if len(rows) != 1:
            continue
        row = rows[0]
        if not (
                row.get("library") == f"lib{lib_num}" and
                row.get("condition") == cond_abbrev and
                row.get("arm") == "A" and
                row.get("arm_key") == "demux_original" and
                row.get("assignment_basis") == "demux" and
                row.get("roster_basis") == "original" and
                os.path.abspath(row.get("assignment_path", "")) ==
                os.path.abspath(current_inputs["assignments"]) and
                os.path.abspath(row.get("receiver_lines", "")) ==
                os.path.abspath(current_inputs["receiver_lines"]) and
                os.path.abspath(row.get("ambient_candidates", "")) ==
                os.path.abspath(current_inputs["ambient_candidates"]) and
                row.get("assignment_update_mode") == "iterative_frozen" and
                row.get("assignment_score_basis") ==
                "original_demux_all_arms"):
            continue
        return prefix
    return None


def materialize_reusable_identity_ambient_arm_a(lib_num, cond_abbrev):
    """Link a validated legacy Arm A bundle into the current applied plan.

    Only immutable numerical artifacts are linked. A fresh lightweight arm
    contract records the current plan/context, so plots use the current applied
    scrutiny strata without pretending the old broad roster was part of Arm A.
    """
    if _identity_ambient_candidate_set() != "applied":
        return None
    target_prefix = get_contam_prefix(
        lib_num, cond_abbrev, "demux_original")
    if all(check_file_exists(target_prefix + suffix) for suffix in
           _identity_ambient_required_output_suffixes(cond_abbrev)) and \
            check_file_exists(target_prefix + ".identity_ambient_arm.tsv"):
        return None
    source_prefix = _reusable_legacy_arm_a_prefix(lib_num, cond_abbrev)
    if source_prefix is None:
        return None

    target_dir = os.path.dirname(target_prefix)
    target_name = os.path.basename(target_prefix)
    os.makedirs(target_dir, exist_ok=True)
    existing = list(Path(target_dir).glob(target_name + ".*"))
    if existing:
        raise RuntimeError(
            "cannot reuse candidate-independent Arm A because the applied "
            "target already contains partial artifacts: " + target_prefix)
    source_dir = Path(os.path.dirname(source_prefix))
    source_name = os.path.basename(source_prefix)
    for source in source_dir.glob(source_name + ".*"):
        if source.name.endswith(".identity_ambient_arm.tsv"):
            continue
        if not (source.is_file() or source.is_symlink()):
            continue
        target = Path(target_dir) / source.name.replace(
            source_name, target_name, 1)
        os.symlink(os.path.abspath(source), target)

    arm_inputs = identity_ambient_arm_inputs(lib_num, "demux_original")
    assignment_sidecar = Path(target_prefix + ".assignments")
    if not assignment_sidecar.exists() and not assignment_sidecar.is_symlink():
        os.symlink(
            os.path.abspath(arm_inputs["assignments"]), assignment_sidecar)
    contract_path = target_prefix + ".identity_ambient_arm.tsv"
    contract = (
        "library\tcondition\tarm\tarm_key\tassignment_basis\troster_basis"
        "\tcandidate_set\tplan_fingerprint\tassignment_path\treceiver_lines"
        "\tambient_candidates\tscrutiny_cells\tcontext_path"
        "\tassignment_update_mode\tassignment_score_basis\n"
        f"lib{lib_num}\t{cond_abbrev}\tA\tdemux_original\tdemux\toriginal"
        f"\t{arm_inputs['candidate_set']}\t{arm_inputs['plan_fingerprint']}"
        f"\t{arm_inputs['assignments']}\t{arm_inputs['receiver_lines']}"
        f"\t{arm_inputs['ambient_candidates']}\t{arm_inputs['scrutiny_cells']}"
        f"\t{arm_inputs['context_path']}\titerative_frozen"
        "\toriginal_demux_all_arms\n"
    )
    _write_if_changed(contract_path, contract)
    if not all(check_file_exists(target_prefix + suffix) for suffix in
               _identity_ambient_required_output_suffixes(cond_abbrev)):
        raise RuntimeError(
            "candidate-independent Arm A reuse produced an incomplete bundle: "
            + target_prefix)
    return source_prefix


def resolve_contam_assignment_sources(value):
    """Expand the user-facing contamination assignment-source selection."""
    if value == "both":
        return list(CONTAM_ASSIGNMENT_SOURCES)
    if value in CONTAM_ASSIGNMENT_SOURCES:
        return [value]
    if value == IDENTITY_AMBIENT_SELECTOR:
        return list(IDENTITY_AMBIENT_ARM_ORDER)
    raise ValueError(f"unknown contamination assignment source: {value}")


def get_contam_assignments_path(lib_num, assignment_source):
    """Return the assignment file used by one contamination-estimator run."""
    if assignment_source in IDENTITY_AMBIENT_ARMS:
        inputs = identity_ambient_arm_inputs(lib_num, assignment_source)
        if inputs is None:
            raise ValueError(
                f"{assignment_source} is not applicable to lib{lib_num}")
        return inputs["assignments"]
    if assignment_source == "demux":
        return get_demux_prefix(lib_num) + ".assignments"
    if assignment_source == "reconciled":
        return get_reconciled_assignments_path(lib_num)
    raise ValueError(f"unknown contamination assignment source: {assignment_source}")


def contam_source_applicable(lib_num, assignment_source):
    """Return (applicable, reason) for a selected contamination source/arm."""
    if assignment_source != "reconciled_replacement":
        return True, ""
    inputs = identity_ambient_arm_inputs(lib_num, assignment_source)
    if inputs is None:
        context = prepare_identity_ambient_comparison(lib_num)
        return False, context["replacement_arm_skip_reason"]
    return True, ""


def get_empty_drops_indiv(lib_num):
    """Return path to individual-level empty drops ambient profile.

    Output lives at the raw prefix since tet_ambient_profile runs on raw
    (unfiltered) barcode counts.
    """
    out_dir = get_demux_dir(lib_num)
    raw_prefix = os.path.join(out_dir, f"lib{lib_num}_raw")
    return raw_prefix + ".contam_prof_empty"


def get_empty_drops_species(lib_num):
    """Return path to species-level empty drops ambient profile.

    Output lives at the raw prefix since tet_ambient_profile runs on raw
    (unfiltered) barcode counts.
    """
    out_dir = get_demux_dir(lib_num)
    raw_prefix = os.path.join(out_dir, f"lib{lib_num}_raw")
    return raw_prefix + ".species_prof_empty"


def get_species_counts(lib_num):
    """Return path to species-diagnostic counts file."""
    return get_demux_prefix(lib_num) + ".species_counts"


def get_species_condf(lib_num):
    """Return path to species condf file alongside species_counts."""
    return get_demux_prefix(lib_num) + ".species_condf"


def get_library_analysis_dir(lib_num, name):
    """Return a shallow per-library analysis folder inside demux_nomito."""
    return os.path.join(get_demux_dir(lib_num), name)


def get_audit_root():
    """Return the root for post-hoc QC/swap-audit outputs."""
    return AUDIT_ROOT


def get_hybrid_root():
    """Return the root for post-hoc hybrid reconciliation outputs."""
    return HYBRID_ROOT


def get_mt_fusion_root():
    """Return the root for mitochondrial fusion-ratio outputs."""
    return MT_FUSION_ROOT


def get_mt_rna_ambient_variant(args):
    """Return the MT output namespace for one optional RNA ambient input.

    Historical MT paths stay byte-for-byte unchanged when the covariate is
    disabled.  Enabled runs are isolated by condition and assignment source;
    controlled four-arm runs additionally include the reconciliation candidate
    set because ``applied`` and ``exploratory`` are different scientific plans.
    """
    if args is None or not args.mt_rna_ambient_condition:
        return None
    condition = str(args.mt_rna_ambient_condition)
    source = str(args.mt_rna_ambient_assignment_source)
    leaf = source
    if source in IDENTITY_AMBIENT_ARMS:
        leaf = f"{source}__{_identity_ambient_candidate_set()}"
    leaf += f"__mt-{args.mt_assignment_source}"
    return os.path.join("rna_ambient", condition, leaf)


def get_mt_lib_dir(lib_num, args=None):
    base = os.path.join(get_mt_fusion_root(), f"lib{lib_num}")
    variant = get_mt_rna_ambient_variant(args)
    return os.path.join(base, variant) if variant else base


def get_mt_prefix(lib_num, args=None):
    return os.path.join(get_mt_lib_dir(lib_num, args), f"lib{lib_num}")


def get_mt_ratio_path(lib_num, args=None):
    return get_mt_prefix(lib_num, args) + ".mt_ratio.tsv"


def get_mt_profile_path(lib_num, args=None):
    return get_mt_prefix(lib_num, args) + ".mt_profile.tsv.gz"


def get_mt_qc_path(lib_num, args=None):
    return get_mt_prefix(lib_num, args) + ".mt_qc.tsv"


def get_refined_simple_assignments_path(lib_num):
    return os.path.join(TETRA_REFINE_ROOT, f"lib{lib_num}", f"lib{lib_num}.assignments_refined")


def _expand_mt_template(template, lib_num):
    if not template:
        return None
    return template.format(lib=lib_num, lib_num=lib_num, library=f"lib{lib_num}")


def get_mt_bam_path(lib_num, args):
    if args.mt_bam_template:
        return _expand_mt_template(args.mt_bam_template, lib_num)
    return get_bam_path(lib_num)


def get_mt_library_id(lib_num, args):
    return _expand_mt_template(args.mt_library_id_template, lib_num)


def get_mt_rna_ambient_path(lib_num, args):
    """Return the selected RNA contamination-rate covariate for MT analysis."""
    if not args.mt_rna_ambient_condition:
        return None
    return get_contam_prefix(
        lib_num, args.mt_rna_ambient_condition,
        args.mt_rna_ambient_assignment_source) + ".contam_rate"


def get_mt_rna_ambient_provenance_path(lib_num, args):
    return get_mt_prefix(lib_num, args) + ".rna_ambient_source.json"


def get_mt_rna_ambient_job_id(lib_num, args, contam_job_ids):
    """Return the exact freshly submitted CONTAM job consumed by MT, if any."""
    if not args.mt_rna_ambient_condition:
        return None
    return contam_job_ids.get((
        args.mt_rna_ambient_condition,
        args.mt_rna_ambient_assignment_source,
        lib_num,
    ))


def get_mt_rna_ambient_provenance(lib_num, args):
    """Return the exact current RNA-ambient input contract for one MT run.

    Four-arm inputs are accepted only when their completion contract matches
    the current reconciliation plan. File digests make a same-plan estimator
    rerun visible to MT rather than silently reusing ratios made from old rates.
    """
    rate_path = get_mt_rna_ambient_path(lib_num, args)
    if not rate_path:
        return None
    if not check_file_exists(rate_path):
        raise FileNotFoundError(
            f"RNA ambient contam_rate missing or empty: {rate_path}")
    source = args.mt_rna_ambient_assignment_source
    provenance = {
        "schema_version": 1,
        "library": int(lib_num),
        "condition": args.mt_rna_ambient_condition,
        "assignment_source": source,
        "mt_assignment_source": args.mt_assignment_source,
        "mt_assignments": _identity_ambient_file_record(
            get_mt_assignments_path(lib_num, args)),
        "rna_ambient_qc_max": args.mt_rna_ambient_max,
        "candidate_set": "",
        "plan_fingerprint": "",
        "contam_rate": _identity_ambient_file_record(rate_path),
        "arm_contract": None,
    }
    if source in IDENTITY_AMBIENT_ARMS:
        if not check_output_exists(
                lib_num, args.mt_rna_ambient_condition, source):
            raise ValueError(
                f"RNA ambient arm {source} is incomplete or belongs to a stale "
                f"comparison plan for lib{lib_num}")
        arm_inputs = identity_ambient_arm_inputs(lib_num, source)
        if arm_inputs is None:
            raise ValueError(
                f"RNA ambient arm {source} is not applicable to lib{lib_num}")
        arm_contract = get_contam_prefix(
            lib_num, args.mt_rna_ambient_condition,
            source) + ".identity_ambient_arm.tsv"
        provenance.update({
            "candidate_set": arm_inputs["candidate_set"],
            "plan_fingerprint": arm_inputs["plan_fingerprint"],
            "arm_contract": _identity_ambient_file_record(arm_contract),
        })
    return provenance


def mt_outputs_complete(lib_num, args):
    """Return whether MT outputs match the currently selected RNA covariate."""
    outputs = (
        get_mt_ratio_path(lib_num, args),
        get_mt_profile_path(lib_num, args),
        get_mt_qc_path(lib_num, args),
    )
    if not all(check_file_exists(path) for path in outputs):
        return False
    if not args.mt_rna_ambient_condition:
        # Preserve the historical completion contract exactly when RNA
        # contamination is not requested.
        return True
    provenance_path = get_mt_rna_ambient_provenance_path(lib_num, args)
    if not check_file_exists(provenance_path):
        return False
    try:
        with open(provenance_path, "r", encoding="utf-8") as handle:
            observed = json.load(handle)
        return observed == get_mt_rna_ambient_provenance(lib_num, args)
    except (OSError, ValueError, TypeError, json.JSONDecodeError):
        return False


def _normalize_mt_manifest_library_id(value):
    """Return a numeric library ID for manifest keys like 7, 07, or lib7."""
    text = str(value).strip()
    if text.lower().startswith("lib"):
        text = text[3:]
    if not text.isdigit():
        return None
    return int(text)


def get_mt_manifest_library_numbers(site_manifest):
    """Return library numbers that have at least one row in the MT ratio manifest."""
    libraries = set()
    with open(site_manifest, "r", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or "library_id" not in reader.fieldnames:
            raise ValueError(
                f"MT site manifest is missing required library_id column: {site_manifest}")
        for row in reader:
            lib_num = _normalize_mt_manifest_library_id(row.get("library_id", ""))
            if lib_num is not None:
                libraries.add(lib_num)
    return libraries


def get_identity_reconciliation_root():
    return IDENTITY_RECONCILIATION_ROOT


def get_reconciled_cells_path(lib_num):
    return os.path.join(
        get_identity_reconciliation_root(), "decisions",
        f"lib{lib_num}.reconciled_cells.tsv.gz")


def get_reconciled_assignments_path(lib_num):
    return os.path.join(
        get_identity_reconciliation_root(), "decisions",
        f"lib{lib_num}.reconciled_single_cell.assignments")


def get_mt_assignments_path(lib_num, args, prefer_future_refined=False):
    # Final MT analysis should consume the adjudicated single-cell identity product.
    # Explicit demux/refined modes remain available only for deliberate comparison runs.
    source = args.mt_assignment_source
    if source in ("auto", "reconciled"):
        return get_reconciled_assignments_path(lib_num)
    if source == "demux":
        return get_demux_prefix(lib_num) + ".assignments"
    if source == "refined":
        return get_refined_simple_assignments_path(lib_num)
    raise ValueError(f"unknown MT assignment source: {source}")


def get_mt_population_prefix(args):
    if args.mt_population_prefix:
        return os.path.abspath(args.mt_population_prefix)
    base = os.path.join(get_mt_fusion_root(), "population")
    variant = get_mt_rna_ambient_variant(args)
    if variant:
        base = os.path.join(base, variant)
    return os.path.join(base, "mt_population")


def get_identity_subdir(args, name):
    return os.path.join(os.path.abspath(args.identity_reconciliation_root), name)


def get_identity_candidate_path(lib_num, args):
    return os.path.join(
        get_identity_subdir(args, "candidates"),
        f"lib{lib_num}.identity_candidates.tsv.gz")


def get_identity_score_pair_path(lib_num, args):
    return os.path.join(
        get_identity_subdir(args, "score_pairs"),
        f"lib{lib_num}.reconciliation_score_pairs.tsv.gz")


def get_identity_score_pair_summary_path(lib_num, args):
    return os.path.join(
        get_identity_subdir(args, "score_pairs"),
        f"lib{lib_num}.reconciliation_score_pair_summary.tsv")


def get_identity_nuclear_score_path(lib_num, args):
    return os.path.join(
        get_identity_subdir(args, "nuclear"),
        f"lib{lib_num}.identity_hypothesis_scores.tsv.gz")


def get_identity_targeted_hypothesis_score_path(lib_num, args):
    return os.path.join(
        get_identity_subdir(args, "nuclear"),
        f"lib{lib_num}.reconciliation_pair_hypothesis_scores.tsv.gz")


def get_identity_targeted_fold_score_path(lib_num, args):
    return os.path.join(
        get_identity_subdir(args, "nuclear"),
        f"lib{lib_num}.reconciliation_pair_site_fold_scores.tsv.gz")


def get_identity_probability_path(lib_num, args):
    return os.path.join(
        get_identity_subdir(args, "nuclear"),
        f"lib{lib_num}.identity_pair_probabilities.tsv.gz")


def get_identity_probability_provenance_path(lib_num, args):
    return os.path.join(
        get_identity_subdir(args, "nuclear"),
        f"lib{lib_num}.identity_pair_probability_provenance.tsv")


def get_identity_score_output_root(args):
    if args.identity_score_output_root:
        return os.path.abspath(args.identity_score_output_root)
    requested = {
        token.strip().upper()
        for token in str(args.stage or "").split(",") if token.strip()
    }
    dirname = (
        "reconciliation_swap_score_summary_v6_1_frozen_reuse"
        if "IDENTITY_SCORE_AGGREGATE_ONLY" in requested
        else "reconciliation_swap_score_summary_v6_1"
    )
    return get_identity_subdir(args, dirname)


def get_mt_identity_score_path(lib_num, args):
    return os.path.join(
        get_identity_subdir(args, "mt"),
        f"lib{lib_num}.mt_identity_scores.tsv.gz")


def get_atac_identity_score_path(lib_num, args):
    return os.path.join(
        get_identity_subdir(args, "atac"),
        f"lib{lib_num}.atac_identity_scores.tsv.gz")


def get_audit_lib_dir(lib_num):
    return os.path.join(get_audit_root(), f"lib{lib_num}")


def get_hybrid_lib_dir(lib_num):
    return os.path.join(get_hybrid_root(), f"lib{lib_num}")


def get_call_qc_path(lib_num):
    return os.path.join(get_audit_lib_dir(lib_num), f"lib{lib_num}.call_qc.tsv.gz")


def get_species_qc_path(lib_num):
    return os.path.join(get_audit_lib_dir(lib_num), f"lib{lib_num}.species_qc.tsv")


def get_swap_report_path(lib_num):
    return os.path.join(get_audit_lib_dir(lib_num), f"lib{lib_num}.swap_report.tsv")


def get_hybrid_summary_path(lib_num):
    return os.path.join(get_hybrid_lib_dir(lib_num), f"lib{lib_num}.hybrid_posthoc_summary.tsv")


def get_ploidy_calls_root():
    return PLOIDY_CALLS_ROOT


def get_ploidy_calls_path(lib_num):
    return os.path.join(get_ploidy_calls_root(), f"lib{lib_num}.ploidy_calls_nn.tsv")


def get_unexpected_component_nn_root():
    """Analysis-only NN cross-check outputs under the current audit root."""
    return os.path.join(get_audit_root(), "unexpected_component_nn")


def has_unexpected_component_signal(lib_num):
    """Return True only for libraries currently flagged with a foreign component."""
    path = get_swap_report_path(lib_num)
    if not os.path.isfile(path):
        return False
    try:
        with open(path, "r", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            row = next(reader, None)
        if not row:
            return False
        return (
            str(row.get("audit_verdict", "")).strip() == "UNEXPECTED_COMPONENT_SIGNAL"
            and str(row.get("top_unexpected_component_supported", "")).strip() not in {"", "NA", "."}
        )
    except Exception:
        return False


def get_refined_assignments_path(lib_num):
    return os.path.join(TETRA_REFINE_ROOT, f"lib{lib_num}", f"lib{lib_num}.refined_assignments")


def get_refine_summary_path(lib_num):
    return os.path.join(TETRA_REFINE_ROOT, f"lib{lib_num}", f"lib{lib_num}.refine_summary")


def get_default_vcf_samples_for_audit(lib_num):
    """Use the demux sample list as the audit sample universe for that library."""
    return get_demux_prefix(lib_num) + ".samples"


def get_expected_pool_metadata():
    return EXPECTED_POOL_METADATA


def get_allowed_identities():
    return ALLOWED_IDENTITIES


def get_log_dir():
    """Return the aggregate orchestration log directory."""
    return os.path.join(AGGREGATE_ROOT, "logs")


def get_script_dir():
    """Return the aggregate generated-sbatch directory."""
    return os.path.join(AGGREGATE_ROOT, "slurm_scripts")


def publish_figure_shortcut(analysis, target):
    """Atomically publish one stable aggregate figure-directory shortcut."""
    names = {
        "ambient_plots": "ambient_plots",
        "ambient_validation": "ambient_validation",
        "ambient_swap_test": "ambient_swap_test",
    }
    analysis = str(analysis).strip().lower()
    if analysis not in names:
        raise ValueError(f"unknown figure analysis: {analysis}")
    aggregate_root = Path(AGGREGATE_ROOT).resolve()
    figure_root = Path(FIGURE_ROOT)
    target = Path(target)
    if not target.is_dir() or target.is_symlink():
        raise RuntimeError(f"figure target is not a real directory: {target}")
    target = target.resolve()
    try:
        target.relative_to(aggregate_root)
    except ValueError as exc:
        raise RuntimeError(
            f"figure target escaped aggregate_library_analysis: {target}") from exc
    figure_root.mkdir(parents=True, exist_ok=True)
    if figure_root.is_symlink() or not figure_root.is_dir():
        raise RuntimeError(f"unsafe central figure root: {figure_root}")
    link = figure_root / names[analysis]
    if link.exists() and not link.is_symlink():
        raise RuntimeError(
            f"refusing to replace non-symlink figure index entry: {link}")
    temporary = figure_root / (
        f".{link.name}.tmp.{os.environ.get('SLURM_JOB_ID', os.getpid())}")
    if temporary.exists() or temporary.is_symlink():
        temporary.unlink()
    relative_target = os.path.relpath(target, figure_root.resolve())
    temporary.symlink_to(relative_target, target_is_directory=True)
    os.replace(temporary, link)
    print(f"Published figure index: {link} -> {relative_target}")


def get_ambient_plot_run_dir(args, cond_list):
    """Return a stable, collision-resistant directory for one plot selection."""
    conditions = [condition["abbrev"] for condition in cond_list]
    assignment_sources = resolve_contam_assignment_sources(
        args.contam_assignment_source)
    if args.ambient_plot_label:
        label = _ambient_condition_slug(args.ambient_plot_label)
        if args.contam_assignment_source == IDENTITY_AMBIENT_SELECTOR:
            label = _ambient_condition_slug(
                f"{label}__{args.reconciliation_candidate_set}")
    elif (len(conditions) == 1 and
          args.contam_assignment_source == IDENTITY_AMBIENT_SELECTOR):
        label = _ambient_condition_slug(
            f"{conditions[0]}__reconciliation-four-arm-"
            f"{args.reconciliation_candidate_set}")
    elif len(conditions) == 1 and assignment_sources == ["demux"]:
        label = _ambient_condition_slug(conditions[0])
    elif len(conditions) == 1 and assignment_sources == ["reconciled"]:
        label = _ambient_condition_slug(f"{conditions[0]}__reconciled")
    elif len(conditions) == 1 and assignment_sources == ["demux", "reconciled"]:
        label = _ambient_condition_slug(
            f"{conditions[0]}__demux-vs-reconciled")
    else:
        selection = [
            f"{condition}\t{source}"
            for condition in conditions for source in assignment_sources
        ]
        if args.contam_assignment_source == IDENTITY_AMBIENT_SELECTOR:
            selection.append(
                "candidate_set\t" + args.reconciliation_candidate_set)
        digest = hashlib.sha1("\n".join(selection).encode("utf-8")).hexdigest()[:10]
        candidate_suffix = (
            f"_{args.reconciliation_candidate_set}"
            if args.contam_assignment_source == IDENTITY_AMBIENT_SELECTOR else "")
        label = f"selected_{len(conditions) * len(assignment_sources)}{candidate_suffix}_{digest}"
    return os.path.join(os.path.abspath(args.ambient_plot_root), label)


def _ensure_directory_index(link_path, target_path):
    """Create one safe directory symlink used by an unchanged helper.

    The helper-facing trees contain links only; the files themselves are
    written beneath the corresponding library.  Existing real directories are
    never replaced automatically because they may contain an unmigrated run.
    """
    link_path = os.path.abspath(link_path)
    target_path = os.path.abspath(target_path)
    os.makedirs(target_path, exist_ok=True)
    os.makedirs(os.path.dirname(link_path), exist_ok=True)

    if os.path.lexists(link_path):
        if os.path.islink(link_path):
            if os.path.realpath(link_path) == os.path.realpath(target_path):
                return
            raise RuntimeError(
                f"compatibility link points somewhere else: {link_path} -> "
                f"{os.readlink(link_path)} (expected {target_path})")
        raise RuntimeError(
            f"refusing to replace existing non-symlink path: {link_path}; "
            "run the migration guide first")

    relative_target = os.path.relpath(target_path, os.path.dirname(link_path))
    os.symlink(relative_target, link_path)


def prepare_output_layout(lib_nums, cond_list, stages, args):
    """Create the shallow production folders and helper compatibility indexes."""
    os.makedirs(AGGREGATE_ROOT, exist_ok=True)
    os.makedirs(get_log_dir(), exist_ok=True)
    os.makedirs(get_script_dir(), exist_ok=True)

    condition_names = {c["abbrev"] for c in cond_list}
    if "TETRA_REFINE" in stages and args.refine_contam_condition:
        condition_names.add(args.refine_contam_condition)
    if "HYBRID" in stages:
        for name in (
                args.hybrid_individual_condition,
                args.hybrid_species_condition,
                args.hybrid_fixed_species_condition):
            if name:
                condition_names.add(name)

    for lib_num in lib_nums:
        if {"CONTAM", "AMBIENT_PLOTS", "GEX_AMBIENT", "TETRA_REFINE", "HYBRID"} & set(stages):
            for condition in sorted(condition_names):
                _ensure_directory_index(
                    os.path.join(CONDITION_INDEX_ROOT, condition, f"lib{lib_num}"),
                    get_contam_dir(lib_num, condition))

        if (AUDIT_ROOT == DEFAULT_AUDIT_ROOT and
                ({"POSTHOC", "POSTHOC_SUMMARY", "UNEXPECTED_COMPONENT_NN",
                  "IDENTITY_RECONCILIATION", "IDENTITY_RECONCILE_ONLY",
                  "IDENTITY_SCORE"} & set(stages))):
            _ensure_directory_index(
                os.path.join(AUDIT_ROOT, f"lib{lib_num}"),
                get_library_analysis_dir(lib_num, "posthoc"))

        if HYBRID_ROOT == DEFAULT_HYBRID_ROOT and "HYBRID" in stages:
            _ensure_directory_index(
                os.path.join(HYBRID_ROOT, f"lib{lib_num}"),
                get_library_analysis_dir(lib_num, "hybrid"))

        if TETRA_REFINE_ROOT == DEFAULT_TETRA_REFINE_ROOT and (
                {"TETRA_REFINE", "POSTHOC", "IDENTITY_RECONCILIATION",
                 "IDENTITY_RECONCILE_ONLY", "IDENTITY_SCORE"} & set(stages)):
            _ensure_directory_index(
                os.path.join(TETRA_REFINE_ROOT, f"lib{lib_num}"),
                get_library_analysis_dir(lib_num, "tetra_refine"))

        if MT_FUSION_ROOT == DEFAULT_MT_FUSION_ROOT and (
                {"MT_FUSION", "MT_POPULATION"} & set(stages)):
            _ensure_directory_index(
                os.path.join(MT_FUSION_ROOT, f"lib{lib_num}"),
                get_library_analysis_dir(lib_num, "mt_fusion"))

    for root in (AUDIT_ROOT, HYBRID_ROOT, MT_FUSION_ROOT,
                 TETRA_REFINE_ROOT, PLOIDY_CALLS_ROOT,
                 IDENTITY_RECONCILIATION_ROOT, CONDITION_INDEX_ROOT):
        os.makedirs(root, exist_ok=True)
    if {"AMBIENT_PLOTS", "AMBIENT_SWAP_TEST"} & set(stages):
        os.makedirs(os.path.abspath(args.ambient_plot_root), exist_ok=True)
    if "GEX_AMBIENT" in stages:
        os.makedirs(os.path.abspath(args.gex_ambient_root), exist_ok=True)


def module_block(modules=None):
    """Return a module load block for sbatch scripts."""
    if modules is None:
        modules = MODULES
    lines = ["module purge"]
    lines.extend(f"module load {m}" for m in modules)
    lines.extend([
        'echo "Runtime host: $(hostname)"',
        'echo "SLURM job: ${SLURM_JOB_ID:-not-set}"',
        'free -h || true',
    ])
    return "\n".join(lines)


def _vcf_daemon_safe_token(value):
    """Return a short token safe in POSIX-shm names and SLURM job names."""
    token = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value)).strip("_.-")
    return token[:72] or "run"


def get_vcf_daemon_nodes():
    """Return the explicit nodes that can run daemon-backed CONDF/DEMUX."""
    nodes = []
    for raw in DAEMON_NODELIST.split(","):
        node = raw.strip()
        if node and node not in nodes:
            if not re.fullmatch(r"[A-Za-z0-9_.-]+", node):
                raise ValueError(f"invalid daemon node name: {node!r}")
            nodes.append(node)
    return nodes


def configure_managed_vcf_segments():
    """Select collision-free shared-memory names for this orchestration run."""
    global MANAGED_VCF_RUN_ID, MANAGED_VCF_READY_FILE
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    MANAGED_VCF_RUN_ID = _vcf_daemon_safe_token(f"{stamp}_{os.getpid()}")
    base = f"/tetvcf_{MANAGED_VCF_RUN_ID}"
    SHARED_VCF["interindiv_20M"] = base
    SHARED_VCF["interindiv_het_10M"] = base + "_het"
    SHARED_VCF["species_20M"] = base + "_species"
    MANAGED_VCF_READY_FILE = f"/dev/shm/{base[1:]}.ready.tsv"
    return MANAGED_VCF_RUN_ID


def get_vcf_daemon_state_dir(run_id=None):
    token = run_id or MANAGED_VCF_RUN_ID
    if not token:
        raise RuntimeError("managed VCF daemon run ID has not been configured")
    return os.path.join(VCF_DAEMON_STATE_ROOT, _vcf_daemon_safe_token(token))


def resolve_vcf_daemon_reference_bam(lib_nums):
    """Choose the first selected BAM as the shared chromosome/TID authority."""
    for lib_num in lib_nums:
        bam = get_bam_path(lib_num)
        if bam and check_file_exists(bam):
            return os.path.abspath(bam)
    return None


def vcf_segment_wait_block(required_segments):
    """Return the compute-node readiness gate for daemon-backed consumers."""
    required = [str(segment) for segment in required_segments]
    tests = " && ".join(
        f'[[ -s "/dev/shm/{segment.lstrip("/")}" ]]' for segment in required)
    if MANAGED_VCF_READY_FILE:
        ready_test = f'[[ -s "{MANAGED_VCF_READY_FILE}" ]] && {tests}'
        detail = (
            f'echo "  Ready marker: {MANAGED_VCF_READY_FILE}" >&2\n'
            f'echo "  Required segments: {", ".join(required)}" >&2')
        return f'''echo "Waiting for orchestrator-managed VCF daemon readiness on $(hostname)"
VCF_WAIT_DEADLINE=$(( $(date +%s) + {VCF_DAEMON_CREATE_WAIT_SECONDS} ))
while ! ( {ready_test} ); do
    if [[ $(date +%s) -ge "$VCF_WAIT_DEADLINE" ]]; then
        echo "ERROR: timed out waiting for the managed VCF daemon on $(hostname)" >&2
        {detail}
        exit 1
    fi
    sleep 5
done
echo "Managed VCF daemon ready on $(hostname): {', '.join(required)}"'''
    return f'''if ! ( {tests} ); then
    echo "ERROR: required externally managed VCF segment(s) missing on $(hostname): {', '.join(required)}" >&2
    exit 1
fi'''


def generate_vcf_daemon_holder_script(node, reference_bam, run_id):
    """Generate one node-pinned foreground VCF daemon holder.

    The holder owns run-scoped shm names, retains the daemon as its child, and
    stops it from an EXIT trap after the cleanup watcher plants TEARDOWN.
    """
    if not MANAGED_VCF_READY_FILE:
        raise RuntimeError("managed VCF ready-file path is not configured")
    state_dir = get_vcf_daemon_state_dir(run_id)
    log_dir = os.path.join(get_log_dir(), "vcf_daemon", run_id)
    script_dir = get_script_dir()
    os.makedirs(state_dir, exist_ok=True)
    os.makedirs(log_dir, exist_ok=True)
    os.makedirs(script_dir, exist_ok=True)
    safe_node = _vcf_daemon_safe_token(node)
    safe_run = _vcf_daemon_safe_token(run_id)
    script_path = os.path.join(
        script_dir, f"vcf_daemon_holder_{safe_run}_{safe_node}.sbatch")
    daemon = os.path.join(SOFTWARE_BIN, "vcf_loader_daemon")
    base = SHARED_VCF["interindiv_20M"]
    main_seg = f"/dev/shm/{SHARED_VCF['interindiv_20M'].lstrip('/')}"
    het_seg = f"/dev/shm/{SHARED_VCF['interindiv_het_10M'].lstrip('/')}"
    species_seg = f"/dev/shm/{SHARED_VCF['species_20M'].lstrip('/')}"
    sentinel = os.path.join(state_dir, "TEARDOWN")
    node_state = os.path.join(state_dir, f"{safe_node}.holder.tsv")

    script = f'''#!/bin/bash
#SBATCH --job-name=tetvcf_{safe_run[-24:]}_{safe_node}
#SBATCH --output={log_dir}/holder_{safe_node}_%j.out
#SBATCH --error={log_dir}/holder_{safe_node}_%j.err
#SBATCH --time={SLURM_TIME}
#SBATCH --cpus-per-task=1
#SBATCH --mem={VCF_DAEMON_MEM_GB}G
#SBATCH --partition={SLURM_PARTITION}
#SBATCH --nodes=1
#SBATCH --nodelist={node}
#SBATCH --chdir={AGGREGATE_ROOT}

set -euo pipefail
{module_block()}

DAEMON="{daemon}"
BASE="{base}"
READY="{MANAGED_VCF_READY_FILE}"
SENTINEL="{sentinel}"
NODE_STATE="{node_state}"
CHILD_PID=""
OWNS_SEGMENTS=0

stop_child() {{
    if [[ -n "${{CHILD_PID:-}}" ]] && kill -0 "$CHILD_PID" 2>/dev/null; then
        echo "Stopping owned VCF daemon child $CHILD_PID on $(hostname)" >&2
        kill "$CHILD_PID" 2>/dev/null || true
        for _ in $(seq 1 120); do
            kill -0 "$CHILD_PID" 2>/dev/null || break
            sleep 1
        done
        if kill -0 "$CHILD_PID" 2>/dev/null; then
            echo "Force-killing VCF daemon child $CHILD_PID after cleanup wait" >&2
            kill -9 "$CHILD_PID" 2>/dev/null || true
        fi
        wait "$CHILD_PID" 2>/dev/null || true
    fi
}}

cleanup() {{
    set +e
    if [[ "${{OWNS_SEGMENTS:-0}}" == "1" ]]; then
        stop_child
        if command -v timeout >/dev/null 2>&1; then
            timeout 30 "$DAEMON" --destroy --name "$BASE" >/dev/null 2>&1 || true
        else
            "$DAEMON" --destroy --name "$BASE" >/dev/null 2>&1 || true
        fi
        rm -f "$READY" "$NODE_STATE"
    fi
}}
trap cleanup EXIT
trap 'exit 143' TERM INT

echo "================================================================"
echo "MANAGED VCF DAEMON HOLDER"
echo "  Run: {run_id}"
echo "  Node: $(hostname)"
echo "  Base segment: $BASE"
echo "  Ready marker: $READY"
echo "  Main BCF: {VCF_SOURCE_PATHS['interindiv_20M']}"
echo "  HET BCF: {VCF_SOURCE_PATHS['interindiv_het_10M']}"
echo "  Species BCF: {VCF_SOURCE_PATHS['species_20M']}"
echo "  Reference BAM: {reference_bam}"
echo "================================================================"

for path in "$DAEMON" "{VCF_SOURCE_PATHS['interindiv_20M']}" "{VCF_SOURCE_PATHS['interindiv_het_10M']}" "{VCF_SOURCE_PATHS['species_20M']}" "{reference_bam}"; do
    [[ -s "$path" ]] || {{ echo "ERROR: required daemon input missing or empty: $path" >&2; exit 1; }}
done
if [[ -f "$SENTINEL" ]]; then
    echo "Teardown was requested before this holder started; no VCF segments will be created"
    exit 0
fi
for path in "{main_seg}" "{het_seg}" "{species_seg}" "$READY"; do
    if [[ -e "$path" ]]; then
        echo "ERROR: run-scoped VCF artifact already exists before holder start: $path" >&2
        exit 1
    fi
done

OWNS_SEGMENTS=1
"$DAEMON" \
    --vcf "{VCF_SOURCE_PATHS['interindiv_20M']}" \
    --het_vcf "{VCF_SOURCE_PATHS['interindiv_het_10M']}" \
    --species_vcf "{VCF_SOURCE_PATHS['species_20M']}" \
    --bam "{reference_bam}" \
    --name "$BASE" \
    --qual 50 \
    --foreground \
    --ready-file "$READY" &
CHILD_PID=$!

DEADLINE=$(( $(date +%s) + {VCF_DAEMON_CREATE_WAIT_SECONDS} ))
while [[ ! -s "$READY" ]]; do
    if ! kill -0 "$CHILD_PID" 2>/dev/null; then
        echo "ERROR: vcf_loader_daemon child exited before readiness" >&2
        wait "$CHILD_PID" || true
        exit 1
    fi
    if [[ $(date +%s) -ge "$DEADLINE" ]]; then
        echo "ERROR: VCF daemon did not become ready within {VCF_DAEMON_CREATE_WAIT_SECONDS}s" >&2
        exit 1
    fi
    sleep 2
done

for path in "{main_seg}" "{het_seg}" "{species_seg}"; do
    [[ -s "$path" ]] || {{ echo "ERROR: readiness marker exists but segment is missing/empty: $path" >&2; exit 1; }}
done
printf 'node\t%s\nholder_jobid\t%s\ndaemon_pid\t%s\nready_file\t%s\n' \
    "$(hostname)" "$SLURM_JOB_ID" "$CHILD_PID" "$READY" > "${{NODE_STATE}}.tmp"
mv -f "${{NODE_STATE}}.tmp" "$NODE_STATE"
echo "VCF daemon is ready; holding until teardown sentinel: $SENTINEL"

while [[ ! -f "$SENTINEL" ]]; do
    if ! kill -0 "$CHILD_PID" 2>/dev/null; then
        echo "ERROR: owned VCF daemon child exited before teardown" >&2
        wait "$CHILD_PID" || true
        exit 1
    fi
    [[ -s "$READY" ]] || {{ echo "ERROR: VCF ready marker disappeared before teardown" >&2; exit 1; }}
    sleep 5
done

echo "Teardown sentinel observed; holder exiting and cleaning its run-scoped segments"
'''
    with open(script_path, "w") as handle:
        handle.write(script)
    os.chmod(script_path, 0o755)
    return script_path


def generate_vcf_daemon_cleanup_script(run_id, holder_job_ids, consumer_job_ids):
    """Generate a watcher that tears holders down after all consumers settle."""
    state_dir = get_vcf_daemon_state_dir(run_id)
    log_dir = os.path.join(get_log_dir(), "vcf_daemon", run_id)
    script_dir = get_script_dir()
    os.makedirs(state_dir, exist_ok=True)
    os.makedirs(log_dir, exist_ok=True)
    os.makedirs(script_dir, exist_ok=True)
    script_path = os.path.join(
        script_dir, f"vcf_daemon_cleanup_{_vcf_daemon_safe_token(run_id)}.sbatch")
    holders = " ".join(str(x) for x in holder_job_ids if str(x).isdigit())
    consumers = " ".join(str(x) for x in consumer_job_ids if str(x).isdigit())
    sentinel = os.path.join(state_dir, "TEARDOWN")
    script = f'''#!/bin/bash
#SBATCH --job-name=tetvcf_cleanup
#SBATCH --output={log_dir}/cleanup_%j.out
#SBATCH --error={log_dir}/cleanup_%j.err
#SBATCH --time={SLURM_TIME}
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --partition={SLURM_PARTITION}
#SBATCH --nodes=1
#SBATCH --chdir={AGGREGATE_ROOT}

set -euo pipefail
{module_block()}

HOLDERS=({holders})
CONSUMERS=({consumers})
SENTINEL="{sentinel}"
HOST_FAILED=0

job_live() {{
    local state
    state="$(squeue -h -j "$1" -o '%T' 2>/dev/null || true)"
    [[ "$state" =~ (PENDING|RUNNING|CONFIGURING|COMPLETING|RESIZING) ]]
}}

job_invalid_dependency() {{
    local reason
    reason="$(squeue -h -j "$1" -o '%r' 2>/dev/null || true)"
    [[ "$reason" == *DependencyNeverSatisfied* || "$reason" == *InvalidDependency* ]]
}}

echo "Watching VCF consumers: ${{CONSUMERS[*]}}"
echo "Watching VCF holders: ${{HOLDERS[*]}}"
mkdir -p "{state_dir}"

while true; do
    ACTIVE=0
    for holder in "${{HOLDERS[@]}}"; do
        if ! job_live "$holder"; then
            echo "ERROR: VCF holder job $holder ended before consumer completion" >&2
            HOST_FAILED=1
        fi
    done
    if [[ "$HOST_FAILED" -ne 0 ]]; then
        for consumer in "${{CONSUMERS[@]}}"; do
            job_live "$consumer" && scancel "$consumer" 2>/dev/null || true
        done
        break
    fi
    for consumer in "${{CONSUMERS[@]}}"; do
        if job_invalid_dependency "$consumer"; then
            echo "Canceling consumer $consumer with unsatisfiable dependency" >&2
            scancel "$consumer" 2>/dev/null || true
            continue
        fi
        job_live "$consumer" && ACTIVE=1
    done
    [[ "$ACTIVE" -eq 0 ]] && break
    sleep 10
done

touch "$SENTINEL"
echo "All VCF consumers settled; teardown sentinel planted: $SENTINEL"

DEADLINE=$(( $(date +%s) + {VCF_DAEMON_CLEANUP_WAIT_SECONDS} ))
while true; do
    LIVE=0
    for holder in "${{HOLDERS[@]}}"; do
        job_live "$holder" && LIVE=1
    done
    [[ "$LIVE" -eq 0 ]] && break
    if [[ $(date +%s) -ge "$DEADLINE" ]]; then
        echo "WARNING: holder cleanup exceeded {VCF_DAEMON_CLEANUP_WAIT_SECONDS}s; using scancel backstop" >&2
        for holder in "${{HOLDERS[@]}}"; do
            job_live "$holder" && scancel "$holder" 2>/dev/null || true
        done
        SCANCEL_DEADLINE=$(( $(date +%s) + 120 ))
        while true; do
            STILL_LIVE=0
            for holder in "${{HOLDERS[@]}}"; do
                job_live "$holder" && STILL_LIVE=1
            done
            [[ "$STILL_LIVE" -eq 0 ]] && break
            [[ $(date +%s) -ge "$SCANCEL_DEADLINE" ]] && break
            sleep 2
        done
        break
    fi
    sleep 5
done

sleep 5
for holder in "${{HOLDERS[@]}}"; do
    if job_live "$holder"; then
        echo "ERROR: VCF holder $holder remained live after teardown" >&2
        HOST_FAILED=1
    fi
done

ALL_JOBS="$(printf '%s\n' "${{HOLDERS[@]}}" "${{CONSUMERS[@]}}" | awk 'NF' | paste -sd, -)"
if [[ -n "$ALL_JOBS" ]]; then
    sacct -j "$ALL_JOBS" --format=JobID%-14,JobName%-32,NodeList%-12,State%-20,ExitCode%-10,Elapsed%-12,MaxRSS%-12 2>/dev/null || true
fi

if [[ "$HOST_FAILED" -ne 0 ]]; then
    echo "ERROR: managed VCF daemon lifecycle ended with a holder failure" >&2
    exit 1
fi
echo "Managed VCF daemon lifecycle completed and run-scoped segments were released"
'''
    with open(script_path, "w") as handle:
        handle.write(script)
    os.chmod(script_path, 0o755)
    return script_path


def resolve_process_script(script_name):
    """Resolve POSTHOC/HYBRID helper scripts from the installed CellBouncer bin.

    SLURM jobs may run on nodes that cannot see the login-node source tree.
    Therefore generated sbatch files must use the deployed shared install path,
    not /home/b/... or the directory containing this orchestrator.
    """
    installed = os.path.join(SOFTWARE_BIN, script_name)
    if os.path.isfile(installed):
        return installed
    return os.path.join("/__MISSING_CELLBOUNCER_HELPER_SCRIPT__", script_name)


def required_posthoc_scripts():
    return [
        resolve_process_script("swap_audit_prepare.py"),
        resolve_process_script("swap_audit_summarize.py"),
        resolve_process_script("swap_audit_aggregate.py"),
    ]


def required_hybrid_scripts():
    return [resolve_process_script("hybrid_posthoc_reconcile.py")]


def required_mt_population_scripts():
    return [resolve_process_script("mt_population_structure.py")]


def resolve_geometry_gate_script():
    """Resolve the gated helper only from compute-visible shared paths."""
    candidates = [
        os.path.join(SOFTWARE_BIN, "geometry_gated_contam_estimate.py"),
        os.path.join(PROCESS_SCRIPTS_DIR, "geometry_gated_contam_estimate.py"),
    ]
    seen = set()
    for candidate in candidates:
        normalized = os.path.abspath(candidate)
        if normalized in seen:
            continue
        seen.add(normalized)
        if not normalized.startswith(("/mnt/beegfs/", "/nvme/")):
            continue
        if os.path.isfile(normalized):
            return normalized
    return os.path.join(
        "/__MISSING_CELLBOUNCER_HELPER_SCRIPT__",
        "geometry_gated_contam_estimate.py",
    )


def required_geometry_gate_scripts():
    return [resolve_geometry_gate_script()]


def resolve_contam_r_script():
    """Resolve contam.R from the deployed tree, then a source checkout."""
    candidates = (
        CONTAM_R,
        os.path.abspath(os.path.join(PROCESS_SCRIPTS_DIR, "..", "plot", "contam.R")),
        os.path.abspath(os.path.join(PROCESS_SCRIPTS_DIR, "..", "Plot", "contam.R")),
    )
    for candidate in candidates:
        if os.path.isfile(candidate):
            return candidate
    return CONTAM_R


def required_ploidy_nn_scripts():
    return [resolve_process_script("run_ploidy_nn_inference.py")]


def required_tetra_refine_scripts():
    return [resolve_process_script("run_tetra_refine_for_library.py")]


def species_profile_alias_block(out_prefix):
    """Shell block that standardizes species-native contam output naming.

    tet_contam_estimate writes the native mixture profile to .contam_prof.
    Downstream orchestration wants the unambiguous species alias .species_prof.
    Create the alias after successful species-native runs without requiring a
    C++ rewrite.
    """
    return (
        f'if [[ -f "{out_prefix}.contam_prof" && ! -e "{out_prefix}.species_prof" ]]; then\n'
        f'    ln -sf "{out_prefix}.contam_prof" "{out_prefix}.species_prof"\n'
        f'    echo "Created species profile alias: {out_prefix}.species_prof -> {out_prefix}.contam_prof"\n'
        f'fi'
    )


# =============================================================================
# File existence checks
# =============================================================================

def check_file_exists(path):
    """Return True if file exists and is non-empty."""
    return bool(path) and os.path.isfile(path) and os.path.getsize(path) > 0


def identity_validation_summary_passes(path):
    """Require a populated validation ledger whose every check passed."""
    if not check_file_exists(path):
        return False, "missing"
    try:
        with open(path, newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        if not rows:
            return False, "empty"
        failed = [
            row for row in rows
            if str(row.get("status", "")).strip().upper() != "PASS"
            or int(float(row.get("n_failures", -1))) != 0
        ]
        if failed:
            return False, ",".join(
                str(row.get("check", "UNKNOWN")) for row in failed[:5]
            )
        return True, "PASS"
    except Exception as exc:
        return False, f"unreadable: {exc}"


def demux_outputs_complete(lib_num, individual_only=False):
    """Return whether the complete selected DEMUX bundle is nonempty."""
    prefix = get_demux_prefix(lib_num)
    required = [prefix + suffix for suffix in (
        ".counts", ".assignments", ".summary", ".diagnostics.gz",
        ".runner_ups.gz", ".samples")]
    if not individual_only:
        raw_prefix = os.path.join(get_demux_dir(lib_num), f"lib{lib_num}_raw")
        required.extend(prefix + suffix for suffix in (
            ".condf", ".species_counts", ".species_condf",
            ".species_samples", ".species_assignments"))
        required.extend(raw_prefix + suffix for suffix in (
            ".counts", ".samples", ".condf", ".species_counts",
            ".species_condf", ".species_samples"))
    return all(check_file_exists(path) for path in required)


def tetra_refine_outputs_complete(lib_num):
    """Return whether every downstream-facing refinement product is ready."""
    return all(check_file_exists(path) for path in (
        get_refined_assignments_path(lib_num),
        get_refined_simple_assignments_path(lib_num),
        get_refine_summary_path(lib_num),
    ))


def check_upstream_ready(lib_num, cond, create_derived=True):
    """Check all prerequisites for running a condition on a library.

    Returns (ready, missing_list) where missing_list has human-readable
    descriptions of what is missing.
    """
    missing = []

    # Native modes intentionally validate only the native bundle they consume.
    prefix = get_demux_prefix(lib_num)
    if cond["mode"] in (2, 4, 5):
        for ext in [".species_counts", ".species_condf", ".species_samples", ".species_assignments"]:
            if not check_file_exists(prefix + ext):
                missing.append(f"species-native demux output: {prefix}{ext}")
        # Native species estimator now consumes original individual assignments
        # to derive weighted species-composition vectors for hybrid fusions.
        for ext in [".assignments", ".samples"]:
            if not check_file_exists(prefix + ext):
                missing.append(f"individual-native assignment metadata for weighted species composition: {prefix}{ext}")
        if not check_file_exists(PANEL_METADATA):
            missing.append(f"panel_metadata for weighted species composition: {PANEL_METADATA}")
    else:
        for ext in [".counts", ".condf", ".assignments", ".samples"]:
            if not check_file_exists(prefix + ext):
                missing.append(f"individual-native demux output: {prefix}{ext}")

    # Native modes need expected_lines at the matching resolution.
    if cond["mode"] in (2, 4, 5):
        try:
            el = get_species_expected_lines(lib_num)
        except Exception as e:
            el = None
            missing.append(f"species expected_lines file ({e})")
        if el is None:
            missing.append("species expected_lines file")
    else:
        el = get_expected_lines(lib_num)
        if el is None:
            missing.append("expected_lines file")

    # Every individual-native estimator receives an explicit library-derived
    # donor roster. CK requires it, and applying it to matched non-CK conditions
    # keeps the candidate universe identical across comparisons.
    if cond["mode"] in (1, 3, "1+SR", LEGACY2C_MODE):
        try:
            candidates = get_individual_ambient_candidates(
                lib_num, create=create_derived)
        except Exception as e:
            missing.append(f"individual ambient candidates ({e})")
        else:
            if candidates is None or not check_file_exists(candidates):
                missing.append("individual ambient candidates file")

    # Only explicit joint/regularization mode needs panel_metadata at solve time.
    if cond["mode"] == "1+SR":
        if not check_file_exists(PANEL_METADATA):
            missing.append(f"panel_metadata: {PANEL_METADATA}")
        if cond["needs_species_counts"]:
            for ext in [".species_counts", ".species_condf", ".species_samples"]:
                if not check_file_exists(prefix + ext):
                    missing.append(f"joint species artifact: {prefix}{ext}")

    # Warm start conditions need empty drops output
    if cond["needs_empty_drops"] == "individual":
        if not check_file_exists(get_empty_drops_indiv(lib_num)):
            missing.append(f"empty_drops_indiv: {get_empty_drops_indiv(lib_num)}")
    elif cond["needs_empty_drops"] in ("species", "species_fixed"):
        if not check_file_exists(get_empty_drops_species(lib_num)):
            missing.append(f"empty_drops_species: {get_empty_drops_species(lib_num)}")

    return (len(missing) == 0, missing)


def check_output_exists(lib_num, cond_abbrev, assignment_source="demux"):
    """Check if outputs already exist for a condition/library pair.

    Species-level conditions are not complete unless the native .species_prof
    exists too. This prevents old SP/SP_WS runs with only .contam_rate from
    being treated as finished after the tet_contam_estimate write-gate fix.
    """
    prefix = get_contam_prefix(lib_num, cond_abbrev, assignment_source)
    if assignment_source in IDENTITY_AMBIENT_ARMS:
        required = (
            ".contam_rate", ".contam_prof", ".allele_ratio",
            ".contam_diagnostics.tsv", ".cell_source_profile.tsv",
            ".profile_fit_diagnostics.tsv", ".condf_coverage.tsv",
            ".run_contract.json", ".decontam.assignments",
            ".identity_ambient_arm.tsv",
        )
        cond = COND_BY_ABBREV.get(cond_abbrev)
        if cond and cond.get("runner") == "geometry_gated":
            required = required + (".geometry_gate_audit.tsv",)
        if not all(check_file_exists(prefix + suffix) for suffix in required):
            return False
        try:
            arm_inputs = identity_ambient_arm_inputs(
                lib_num, assignment_source)
            if arm_inputs is None:
                return False
            with open(
                    prefix + ".identity_ambient_arm.tsv", "r",
                    encoding="utf-8", newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            return (
                len(rows) == 1 and
                rows[0].get("library") == f"lib{lib_num}" and
                rows[0].get("condition") == cond_abbrev and
                rows[0].get("arm") == arm_inputs["arm"] and
                rows[0].get("arm_key") == assignment_source and
                rows[0].get("assignment_basis") ==
                arm_inputs["assignment_basis"] and
                rows[0].get("roster_basis") == arm_inputs["roster_basis"] and
                rows[0].get("candidate_set") == arm_inputs["candidate_set"] and
                rows[0].get("plan_fingerprint") ==
                arm_inputs["plan_fingerprint"] and
                os.path.abspath(rows[0].get("assignment_path", "")) ==
                os.path.abspath(arm_inputs["assignments"]) and
                os.path.abspath(rows[0].get("receiver_lines", "")) ==
                os.path.abspath(arm_inputs["receiver_lines"]) and
                os.path.abspath(rows[0].get("ambient_candidates", "")) ==
                os.path.abspath(arm_inputs["ambient_candidates"]) and
                os.path.abspath(rows[0].get("scrutiny_cells", "")) ==
                os.path.abspath(arm_inputs["scrutiny_cells"]) and
                os.path.abspath(rows[0].get("context_path", "")) ==
                os.path.abspath(arm_inputs["context_path"]) and
                rows[0].get("assignment_update_mode") ==
                "iterative_frozen" and
                rows[0].get("assignment_score_basis") ==
                "original_demux_all_arms")
        except (OSError, ValueError, csv.Error):
            return False
    if not check_file_exists(prefix + ".contam_rate"):
        return False

    if not check_file_exists(prefix + ".contam_prof"):
        return False

    cond = COND_BY_ABBREV.get(cond_abbrev)
    if cond and cond.get("runner") == "geometry_gated":
        return all(check_file_exists(prefix + suffix) for suffix in (
            ".contam_rate",
            ".contam_prof",
            ".allele_ratio",
            ".contam_diagnostics.tsv",
            ".profile_fit_diagnostics.tsv",
            ".condf_coverage.tsv",
            ".run_contract.json",
            ".geometry_gate_audit.tsv",
        ))
    if cond and cond.get("final_ck_condition") and cond.get("runner") != "geometry_gated":
        return all(check_file_exists(prefix + suffix) for suffix in (
            ".contam_rate",
            ".contam_prof",
            ".allele_ratio",
            ".contam_diagnostics.tsv",
            ".profile_fit_diagnostics.tsv",
            ".condf_coverage.tsv",
            ".run_contract.json",
        ))
    if cond and cond["mode"] == LEGACY2C_MODE:
        return (check_file_exists(prefix + ".decontam.assignments") and
                check_file_exists(prefix + ".legacy2c_candidates.txt") and
                check_file_exists(prefix + ".condf_coverage.tsv") and
                check_file_exists(prefix + ".profile_fit_diagnostics.tsv") and
                check_file_exists(prefix + ".legacy2c_diagnostics.tsv") and
                check_file_exists(prefix + ".run_contract.json"))
    if cond and cond["mode"] in (2, 4, 5, "1+SR"):
        # Native species runs write .contam_prof; the orchestrator also creates
        # .species_prof as an explicit alias for downstream clarity. Accept the
        # primary profile when diagnosing older completed runs whose alias was
        # not created yet.
        return (check_file_exists(prefix + ".species_prof") or
                check_file_exists(prefix + ".contam_prof"))

    return True


# =============================================================================
# Stage 1 CONDF: .condf generation
# =============================================================================

def generate_condf_script(vcf_key, force=False):
    """Generate sbatch script for .condf generation for one VCF set.

    Returns (script_path, condf_output_path).
    """
    condf_path = CONDF_PATHS[vcf_key]
    job_name = f"condf_{vcf_key}"
    log_dir = get_log_dir()
    script_dir = get_script_dir()

    binary = os.path.join(SOFTWARE_BIN, "demux_parallel")

    # Build validation block: wait for this run's node-local daemon readiness.
    shm_segment = SHARED_VCF[vcf_key]
    validation_block = vcf_segment_wait_block([shm_segment])

    # Build skip/cleanup behavior. A forced CONDF refresh removes the old
    # artifact before demux_parallel runs so the result is genuinely rebuilt,
    # rather than relying on the binary to overwrite an existing file.
    skip_block = ""
    force_cleanup_block = ""
    if force:
        force_cleanup_block = (
            f'echo "Removing stale CONDF before regeneration: {condf_path}"\n'
            f'rm -f "{condf_path}"\n'
        )
    else:
        skip_block = (
            f'if [[ -s "{condf_path}" ]]; then\n'
            f'    echo "✅ .condf already exists: {condf_path}"\n'
            f'    echo "Skipping (use --regenerate-condf or --force to regenerate)"\n'
            f'    exit 0\n'
            f'fi\n'
        )

    shm_name = SHARED_VCF[vcf_key]
    # demux_parallel --dump_conditional appends ".condf" to the -o prefix,
    # so strip the .condf extension from condf_path to get the correct output.
    condf_prefix = condf_path
    if condf_prefix.endswith(".condf"):
        condf_prefix = condf_prefix[:-6]
    command = f'{binary} --dump_conditional -o "{condf_prefix}" --shared_vcf {shm_name} -t 160'

    # Output validation
    out_val = (
        f'if [[ ! -s "{condf_path}" ]]; then\n'
        f'    echo "❌ ERROR: Expected output missing: {condf_path}"\n'
        f'    exit 1\n'
        f'fi\n'
        f'echo "✅ .condf generated: {condf_path}"'
    )

    daemon_node_line = f"#SBATCH --nodelist={DAEMON_NODELIST}" if DAEMON_NODELIST else ""
    condf_cpus = 95 if MANAGED_VCF_READY_FILE else 96

    script = f"""#!/bin/bash
#SBATCH --job-name={job_name}
#SBATCH --output={log_dir}/{job_name}_%j.out
#SBATCH --error={log_dir}/{job_name}_%j.err
#SBATCH --time={SLURM_TIME}
#SBATCH --cpus-per-task={condf_cpus}
#SBATCH --mem=512G
#SBATCH --partition={SLURM_PARTITION}
#SBATCH --nodes=1
{daemon_node_line}

set -euo pipefail
{module_block()}

echo "================================================================"
echo "STAGE 1 CONDF: .condf generation"
echo "  VCF set: {vcf_key}"
echo "  Shared VCF: {SHARED_VCF[vcf_key]}"
echo "  Output: {condf_path}"
echo "  Started: $(date)"
echo "  Node: $(hostname)"
echo "================================================================"
echo ""

# Ensure output directory exists
mkdir -p "{CONDF_DIR}"

# Validate shared memory segment
{validation_block}

# Skip if output already exists, or remove it for an explicit refresh
{skip_block}{force_cleanup_block}
echo "Running: {command}"
echo ""

{command}

echo ""
# Validate output
{out_val}
echo "Finished: $(date)"
"""

    os.makedirs(log_dir, exist_ok=True)
    os.makedirs(script_dir, exist_ok=True)
    script_path = os.path.join(script_dir, f"{job_name}.sbatch")
    with open(script_path, "w") as f:
        f.write(script)
    os.chmod(script_path, stat.S_IRWXU | stat.S_IRGRP | stat.S_IXGRP |
             stat.S_IROTH | stat.S_IXOTH)
    return script_path, condf_path


# =============================================================================
# Stage 2 DEMUX: Demux (per library)
# =============================================================================

def generate_demux_script(lib_num, condf_job_ids=None, force=False, use_het_vcf=True, individual_only=False):
    """Generate sbatch script for demux of one library.

    Two passes within a single job:
      Pass 1: filtered barcodes, 80 threads (produces .counts, .assignments, etc.)
      Pass 2: raw barcodes, 32 threads, --skip_assignment (produces raw .counts for empties)

    Returns script_path.
    """
    job_name = f"demux_lib{lib_num}"
    log_dir = get_log_dir()
    script_dir = get_script_dir()

    binary = os.path.join(SOFTWARE_BIN, "demux_parallel")
    daemon = os.path.join(SOFTWARE_BIN, "vcf_loader_daemon")
    bam = get_bam_path(lib_num)
    out_dir = get_demux_dir(lib_num)
    prefix_filtered = get_demux_prefix(lib_num)
    prefix_raw = os.path.join(out_dir, f"lib{lib_num}_raw")
    daemon_ready_block = vcf_segment_wait_block(
        [SHARED_VCF["interindiv_20M"]]
        if individual_only else [
            SHARED_VCF["interindiv_20M"],
            *([SHARED_VCF["interindiv_het_10M"]] if use_het_vcf else []),
            SHARED_VCF["species_20M"],
        ])
    demux_exclusive_line = (
        "" if MANAGED_VCF_READY_FILE else "#SBATCH --exclusive")

    filtered_bc = get_filtered_barcodes(lib_num)
    el = get_expected_lines(lib_num)

    if individual_only:
        # Focused panel-causality mode: preserve the standard filtered-barcode
        # individual assignment command and omit all optional/secondary work.
        dep_line = ""
        if condf_job_ids:
            dep_ids = ":".join(str(j) for j in condf_job_ids)
            dep_line = f"#SBATCH --dependency=afterok:{dep_ids}"

        val_lines = []
        for vf in [bam, filtered_bc, el]:
            val_lines.append(
                f'if [[ ! -s "{vf}" ]]; then\n'
                f'    echo "ERROR: Required file missing or empty: {vf}"\n'
                f'    exit 1\n'
                f'fi')
        validation_block = "\n".join(val_lines)

        skip_block = ""
        if not force:
            skip_block = (
                f'if [[ -s "{prefix_filtered}.counts" ]] && '
                f'[[ -s "{prefix_filtered}.assignments" ]] && '
                f'[[ -s "{prefix_filtered}.summary" ]] && '
                f'[[ -s "{prefix_filtered}.diagnostics.gz" ]] && '
                f'[[ -s "{prefix_filtered}.runner_ups.gz" ]] && '
                f'[[ -s "{prefix_filtered}.samples" ]]; then\n'
                f'    echo "Complete individual-only demux outputs already exist for lib{lib_num}"\n'
                f'    echo "Skipping (use --force to rerun)"\n'
                f'    exit 0\n'
                f'fi\n'
            )

        count_policy = "--force_recount" if force else "--reuse_counts"
        shm_interindiv = SHARED_VCF["interindiv_20M"]
        pass1_cmd = (
            f'{binary} -b "{bam}" -o "{prefix_filtered}" '
            f'--shared_vcf {shm_interindiv} -f '
            f'--barcodes "{filtered_bc}" '
            f'-I "{el}" '
            f'--dump_selection_audit '
            f'{count_policy} '
            f'-t 80'
        )

        expected = [
            f"{prefix_filtered}.counts",
            f"{prefix_filtered}.assignments",
            f"{prefix_filtered}.summary",
            f"{prefix_filtered}.diagnostics.gz",
            f"{prefix_filtered}.runner_ups.gz",
            f"{prefix_filtered}.samples",
        ]
        out_check_lines = []
        for ef in expected:
            out_check_lines.append(
                f'if [[ ! -s "{ef}" ]]; then\n'
                f'    echo "ERROR: Expected output missing or empty: {ef}"\n'
                f'    MISSING=$((MISSING + 1))\n'
                f'fi')
        out_check_block = "\n".join(out_check_lines)

        info = LIB_INFO.get(lib_num, {})
        daemon_node_line = (
            f"#SBATCH --nodelist={DAEMON_NODELIST}" if DAEMON_NODELIST else ""
        )

        script = f"""#!/bin/bash
#SBATCH --job-name={job_name}
#SBATCH --output={log_dir}/{job_name}_%j.out
#SBATCH --error={log_dir}/{job_name}_%j.err
#SBATCH --time={SLURM_TIME}
#SBATCH --cpus-per-task=90
#SBATCH --mem=1000G
#SBATCH --partition={SLURM_PARTITION}
#SBATCH --nodes=1
{demux_exclusive_line}
{daemon_node_line}
{dep_line}

set -euo pipefail
{module_block()}

echo "================================================================"
echo "STAGE 2 DEMUX: Individual-only demux lib{lib_num}"
echo "  Batch: {info.get('batch', '?')}  Type: {info.get('lib_type', '?')}"
echo "  Tets: {info.get('tets', '?')}  Dips: {info.get('dips', '?')}"
echo "  Started: $(date)"
echo "  Node: $(hostname)"
echo "  Main individual panel: enabled"
echo "  HET panel: disabled"
echo "  Species panel: disabled"
echo "  Raw-barcode pass: disabled"
echo "  Pileup sidecars: disabled"
echo "  Selection audit: enabled"
echo "================================================================"
echo ""

mkdir -p "{out_dir}"

{daemon_ready_block}

{validation_block}

{skip_block}echo "=== FILTERED INDIVIDUAL DEMUX (80 threads) ==="
echo "Running: {pass1_cmd}"
echo ""

{pass1_cmd}

echo ""
echo "=== OUTPUT VALIDATION ==="
MISSING=0
{out_check_block}
if [[ "${{MISSING}}" -gt 0 ]]; then
    echo "ERROR: lib{lib_num} individual-only demux failed: ${{MISSING}} output files missing"
    exit 1
fi

echo "lib{lib_num} individual-only demux completed successfully"
echo "Finished: $(date)"
"""

        os.makedirs(log_dir, exist_ok=True)
        os.makedirs(script_dir, exist_ok=True)
        script_path = os.path.join(script_dir, f"{job_name}.sbatch")
        with open(script_path, "w") as f:
            f.write(script)
        os.chmod(script_path, stat.S_IRWXU | stat.S_IRGRP | stat.S_IXGRP |
                 stat.S_IROTH | stat.S_IXOTH)
        return script_path

    # SLURM dependency line
    dep_line = ""
    if condf_job_ids:
        dep_ids = ":".join(str(j) for j in condf_job_ids)
        dep_line = f"#SBATCH --dependency=afterok:{dep_ids}"

    # Validation block. expected_lines is REQUIRED for demux assignment so the
    # demuxer constrains its candidate identity set to the legal pool. Without
    # -i/-I, demux_parallel picks from all 28*28 combo identities and frequently
    # assigns cells to out-of-pool combos, which downstream causes the contam
    # estimator to attribute residual signal from misassigned cells to JOS3C1
    # (or other singletons) as phantom ambient mass.
    local_condf_interindiv = get_local_condf_path(lib_num, "interindiv_20M")
    condf_interindiv_source = (
        local_condf_interindiv
        if os.path.isfile(local_condf_interindiv) and os.path.getsize(local_condf_interindiv) > 0
        else CONDF_PATHS["interindiv_20M"]
    )
    val_files = [bam, filtered_bc, el, PANEL_METADATA, condf_interindiv_source]
    val_lines = []
    for vf in val_files:
        val_lines.append(
            f'if [[ ! -s "{vf}" ]]; then\n'
            f'    echo "❌ ERROR: Required file missing or empty: {vf}"\n'
            f'    exit 1\n'
            f'fi')
    validation_block = "\n".join(val_lines)

    # Skip check. A partial bundle is never treated as complete.
    skip_block = ""
    if not force:
        skip_block = (
            f'if [[ -s "{prefix_filtered}.counts" ]] && '
            f'[[ -s "{prefix_filtered}.assignments" ]] && '
            f'[[ -s "{prefix_filtered}.summary" ]] && '
            f'[[ -s "{prefix_filtered}.diagnostics.gz" ]] && '
            f'[[ -s "{prefix_filtered}.runner_ups.gz" ]] && '
            f'[[ -s "{prefix_filtered}.samples" ]] && '
            f'[[ -s "{prefix_filtered}.condf" ]] && '
            f'[[ -s "{prefix_filtered}.species_counts" ]] && '
            f'[[ -s "{prefix_filtered}.species_condf" ]] && '
            f'[[ -s "{prefix_filtered}.species_samples" ]] && '
            f'[[ -s "{prefix_filtered}.species_assignments" ]] && '
            f'[[ -s "{prefix_raw}.counts" ]] && '
            f'[[ -s "{prefix_raw}.samples" ]] && '
            f'[[ -s "{prefix_raw}.condf" ]] && '
            f'[[ -s "{prefix_raw}.species_counts" ]] && '
            f'[[ -s "{prefix_raw}.species_condf" ]] && '
            f'[[ -s "{prefix_raw}.species_samples" ]]; then\n'
            f'    echo "✅ Complete demux outputs already exist for lib{lib_num}"\n'
            f'    echo "Skipping (use --force to rerun)"\n'
            f'    exit 0\n'
            f'fi\n'
        )

    shm_interindiv = SHARED_VCF["interindiv_20M"]
    shm_het = SHARED_VCF["interindiv_het_10M"]
    shm_species = SHARED_VCF["species_20M"]
    het_vcf_arg = f"--shared_het_vcf {shm_het} " if use_het_vcf else ""
    count_policy = "--force_recount" if force else "--reuse_counts"

    # Pass 1 always requests the five-column selection audit used downstream.
    # Nuclear identity scoring consumes the pileup sidecars. They can only be
    # produced during a real BAM counting pass, so request them on forced full
    # recounts;
    # demux_parallel v2.13 refuses --dump_pileup on the reuse path fail-closed.
    pileup_arg = f'--dump_pileup "{prefix_filtered}" ' if force else ""

    # Pass 1: filtered barcodes with species output and heterozygous-site
    # diagnostics (uses shared memory VCF daemon)
    # -f (--disable_conditional): skip .condf computation since Stage 1 CONDF
    # already generated them; .condf files are used downstream by tet_contam_estimate
    # and tet_ambient_profile, not by demux_parallel during counting.
    # -I expected_lines: restrict identity candidates to the library's expected
    # pool. parse_idfile in demux_parallel parses singlet rows (bare identity
    # names) and combo rows (A+B) from the same file, so one -I covers both.
    # -I (idfile_doublet_given) also triggers a post-round filter_identities
    # step that drops first-round assignments to singlet identities that were
    # not explicitly listed in the file. Without this restriction,
    # demux_parallel considers all 28*28=784 combo identities, frequently
    # produces out-of-pool assignments, and corrupts downstream ambient
    # profile estimation.
    pass1_cmd = (
        f'{binary} -b "{bam}" -o "{prefix_filtered}" '
        f'--shared_vcf {shm_interindiv} '
        f'{het_vcf_arg}-f '
        f'--barcodes "{filtered_bc}" '
        f'-I "{el}" '
        f'--species_shared_vcf {shm_species} '
        f'--species_counts_output '
        f'--species_assignment_output '
        f'--panel_metadata "{PANEL_METADATA}" '
        f'--species_panel_mode count_only '
        f'--dump_selection_audit '
        f'{pileup_arg}'
        f'{count_policy} '
        f'-t 80'
    )

    # Pass 2: raw barcodes, skip assignment (uses shared memory VCF daemon)
    # Species counts with count_only mode are needed for empty-drops estimation
    # downstream.  Even though --skip_assignment suppresses assignment writing,
    # keep -I on this demux_parallel invocation so every demux/counting command
    # receives the same legal-pool restriction and future changes cannot silently
    # re-enable an unrestricted raw assignment path.
    # Thread count kept low (8) because dual-panel counting (interindividual +
    # species) allocates per-thread count structures for both panels; 32 threads
    # exceeds 1TB on high-read libraries. 8 threads keeps peak memory under 500GB.
    pass2_cmd = (
        f'{binary} -b "{bam}" -o "{prefix_raw}" '
        f'--shared_vcf {shm_interindiv} -f '
        f'-I "{el}" '
        f'--species_shared_vcf {shm_species} '
        f'--species_counts_output '
        f'--panel_metadata "{PANEL_METADATA}" '
        f'--species_panel_mode count_only '
        f'--skip_assignment '
        f'{count_policy} '
        f'-t 8'
    )

    # Keep the production demux directory self-contained. Stage 1 may still
    # generate library-independent CONDF files in CONDF_DIR, but after DEMUX
    # each selected panel-family CONDF is copied into the library's demux_nomito
    # directory. The filtered/raw .condf links are relative local links, so
    # moving or retiring the shared generation directory cannot break them.
    condf_interindiv = condf_interindiv_source
    local_condf_interindiv = get_local_condf_path(lib_num, "interindiv_20M")
    local_condf_het = get_local_condf_path(lib_num, "interindiv_het_10M")
    local_condf_species = get_local_condf_path(lib_num, "species_20M")
    shared_condf_het = CONDF_PATHS["interindiv_het_10M"]
    shared_condf_species = CONDF_PATHS["species_20M"]
    symlink_block = (
        f'echo "=== POST-DEMUX LOCAL CONDF BUNDLE ==="\n'
        f'if [[ ! -s "{local_condf_interindiv}" ]]; then cp -p "{condf_interindiv}" "{local_condf_interindiv}"; fi\n'
        f'if [[ -s "{shared_condf_het}" && ! -s "{local_condf_het}" ]]; then cp -p "{shared_condf_het}" "{local_condf_het}"; fi\n'
        f'if [[ -s "{shared_condf_species}" && ! -s "{local_condf_species}" ]]; then cp -p "{shared_condf_species}" "{local_condf_species}"; fi\n'
        f'ln -sfn "demux_input_20M.condf" "{prefix_filtered}.condf"\n'
        f'ln -sfn "demux_input_20M.condf" "{prefix_raw}.condf"\n'
        f'ln -sfn "$(basename \"{prefix_filtered}.samples\")" "{prefix_raw}.samples"\n'
        f'if [[ -s "{prefix_filtered}.species_samples" && ! -e "{prefix_raw}.species_samples" ]]; then ln -sfn "$(basename \"{prefix_filtered}.species_samples\")" "{prefix_raw}.species_samples"; fi\n'
        f'echo "Local CONDF copies and relative prefix links are ready"'
    )

    # Output validation
    expected = [
        f"{prefix_filtered}.counts",
        f"{prefix_filtered}.assignments",
        f"{prefix_filtered}.summary",
        f"{prefix_filtered}.diagnostics.gz",
        f"{prefix_filtered}.runner_ups.gz",
        f"{prefix_filtered}.samples",
        f"{prefix_filtered}.condf",
        f"{prefix_filtered}.species_counts",
        f"{prefix_filtered}.species_condf",
        f"{prefix_filtered}.species_samples",
        f"{prefix_filtered}.species_assignments",
        f"{prefix_raw}.counts",
        f"{prefix_raw}.samples",
        f"{prefix_raw}.condf",
        f"{prefix_raw}.species_counts",
        f"{prefix_raw}.species_condf",
        f"{prefix_raw}.species_samples",
    ]
    # Forced full recounts must publish the pileup sidecars required by nuclear
    # identity scoring; their absence is a job failure, not a warning.
    if force:
        expected.append(f"{prefix_filtered}.pileup_sites.tsv.gz")
        expected.append(f"{prefix_filtered}.pileup_obs.tsv.gz")
        expected.append(f"{prefix_filtered}.pileup_molecules.tsv.gz")
    out_check_lines = []
    for ef in expected:
        out_check_lines.append(
            f'if [[ ! -s "{ef}" ]]; then\n'
            f'    echo "❌ ERROR: Expected output missing or empty: {ef}"\n'
            f'    MISSING=$((MISSING + 1))\n'
            f'fi')
    out_check_block = "\n".join(out_check_lines)

    info = LIB_INFO.get(lib_num, {})

    daemon_node_line = f"#SBATCH --nodelist={DAEMON_NODELIST}" if DAEMON_NODELIST else ""

    script = f"""#!/bin/bash
#SBATCH --job-name={job_name}
#SBATCH --output={log_dir}/{job_name}_%j.out
#SBATCH --error={log_dir}/{job_name}_%j.err
#SBATCH --time={SLURM_TIME}
#SBATCH --cpus-per-task=90
#SBATCH --mem=1000G
#SBATCH --partition={SLURM_PARTITION}
#SBATCH --nodes=1
{demux_exclusive_line}
{daemon_node_line}
{dep_line}

set -euo pipefail
{module_block()}

echo "================================================================"
echo "STAGE 2 DEMUX: Demux lib{lib_num}"
echo "  Batch: {info.get('batch', '?')}  Type: {info.get('lib_type', '?')}"
echo "  Tets: {info.get('tets', '?')}  Dips: {info.get('dips', '?')}"
echo "  Started: $(date)"
echo "  Node: $(hostname)"
echo "  HET VCF panel: {'enabled' if use_het_vcf else 'disabled'}"
echo "  Selection audit: enabled"
echo "  Pileup sidecars: {'enabled (forced recount)' if force else 'disabled (reuse path)'}"
echo "================================================================"
echo ""

mkdir -p "{out_dir}"

{daemon_ready_block}

# Validate inputs
{validation_block}

# Skip if outputs exist
{skip_block}
echo "=== PASS 1: Filtered barcodes (80 threads) ==="
echo "Running: {pass1_cmd}"
echo ""

{pass1_cmd}

echo ""
echo "=== PASS 2: Raw barcodes (8 threads, skip_assignment) ==="
echo "Running: {pass2_cmd}"
echo ""

{pass2_cmd}

echo ""
{symlink_block}

echo ""
echo "=== OUTPUT VALIDATION ==="
MISSING=0
{out_check_block}
if [[ "${{MISSING}}" -gt 0 ]]; then
    echo "❌ ERROR: lib{lib_num} demux failed: ${{MISSING}} output files missing"
    exit 1
fi

echo "✅ lib{lib_num} demux completed successfully"
echo "Finished: $(date)"
"""

    os.makedirs(log_dir, exist_ok=True)
    os.makedirs(script_dir, exist_ok=True)
    script_path = os.path.join(script_dir, f"{job_name}.sbatch")
    with open(script_path, "w") as f:
        f.write(script)
    os.chmod(script_path, stat.S_IRWXU | stat.S_IRGRP | stat.S_IXGRP |
             stat.S_IROTH | stat.S_IXOTH)
    return script_path


# =============================================================================
# Stage 3 EMPTY_DROPS: Empty drops ambient profile (per library)
# =============================================================================

def generate_empty_drops_script(lib_num, demux_job_id=None, force=False):
    """Generate sbatch script for empty drops estimation on one library.

    Runs tet_ambient_profile twice: --interindividual then --interspecies.
    Uses the RAW barcode counts (from DEMUX Pass 2) so that empty droplets
    are present in the .counts file. The tool partitions barcodes into
    "cell" (in filtered list) and "empty" (not in filtered list).

    Returns script_path.
    """
    job_name = f"empty_lib{lib_num}"
    log_dir = get_log_dir()
    script_dir = get_script_dir()

    binary = os.path.join(SOFTWARE_BIN, "tet_ambient_profile")
    out_dir = get_demux_dir(lib_num)
    raw_prefix = os.path.join(out_dir, f"lib{lib_num}_raw")
    prefix_filtered = get_demux_prefix(lib_num)
    condf_interindiv = raw_prefix + ".condf"
    # Native species empty-drop estimation must use the per-library
    # species-native condf generated by demux_parallel. The centralized
    # CONDF_PATHS["species_20M"] artifact is the old individual-shaped condf
    # from --dump_conditional and will fail the dimensional guard.
    condf_species = raw_prefix + ".species_condf"
    filtered_bc = get_filtered_barcodes(lib_num)
    el = get_expected_lines(lib_num)
    if el is None:
        raise FileNotFoundError(f"lib{lib_num}: expected_lines file required for EMPTY_DROPS")
    species_el = get_species_expected_lines(lib_num)
    if species_el is None:
        raise FileNotFoundError(f"lib{lib_num}: species expected_lines file required for EMPTY_DROPS")
    decompressed_bc = os.path.join(out_dir, "filtered_barcodes.tsv")

    out_indiv = raw_prefix + ".contam_prof_empty"
    out_species = raw_prefix + ".species_prof_empty"

    # SLURM dependency
    dep_line = ""
    if demux_job_id is not None:
        dep_line = f"#SBATCH --dependency=afterok:{demux_job_id}"

    # Self-healing repair uses only the self-contained demux directory.
    # The migrated production layout stores demux_input_20M.condf locally and
    # both filtered/raw prefix links resolve to that file.
    local_condf_interindiv = get_local_condf_path(lib_num, "interindiv_20M")
    symlink_repair_block = (
        f'FILTERED_PREFIX="{prefix_filtered}"\n'
        f'RAW_PREFIX="{raw_prefix}"\n'
        f'LOCAL_CONDF="{local_condf_interindiv}"\n'
        f'if [[ -s "${{LOCAL_CONDF}}" && ! -e "${{FILTERED_PREFIX}}.condf" ]]; then\n'
        f'    ln -sfn "demux_input_20M.condf" "${{FILTERED_PREFIX}}.condf"\n'
        f'    echo "Created missing local link: ${{FILTERED_PREFIX}}.condf -> demux_input_20M.condf"\n'
        f'fi\n'
        f'if [[ -s "${{LOCAL_CONDF}}" && ! -e "${{RAW_PREFIX}}.condf" ]]; then\n'
        f'    ln -sfn "demux_input_20M.condf" "${{RAW_PREFIX}}.condf"\n'
        f'    echo "Created missing local link: ${{RAW_PREFIX}}.condf -> demux_input_20M.condf"\n'
        f'fi\n'
        f'if [[ -f "${{FILTERED_PREFIX}}.species_condf" && ! -e "${{RAW_PREFIX}}.species_condf" ]]; then\n'
        f'    ln -sfn "${{FILTERED_PREFIX##*/}}.species_condf" "${{RAW_PREFIX}}.species_condf"\n'
        f'    echo "Created missing symlink: ${{RAW_PREFIX}}.species_condf -> ${{FILTERED_PREFIX}}.species_condf"\n'
        f'fi\n'
        f'if [[ -f "${{FILTERED_PREFIX}}.species_samples" && ! -e "${{RAW_PREFIX}}.species_samples" ]]; then\n'
        f'    ln -sfn "${{FILTERED_PREFIX##*/}}.species_samples" "${{RAW_PREFIX}}.species_samples"\n'
        f'    echo "Created missing symlink: ${{RAW_PREFIX}}.species_samples -> ${{FILTERED_PREFIX}}.species_samples"\n'
        f'fi\n'
        f'if [[ -f "${{FILTERED_PREFIX}}.samples" && ! -e "${{RAW_PREFIX}}.samples" ]]; then\n'
        f'    ln -sfn "${{FILTERED_PREFIX##*/}}.samples" "${{RAW_PREFIX}}.samples"\n'
        f'    echo "Created missing symlink: ${{RAW_PREFIX}}.samples -> ${{FILTERED_PREFIX}}.samples"\n'
        f'fi'
    )

    # Validation block: check raw prefix files (produced by DEMUX Pass 2 + symlinks)
    val_files = [
        raw_prefix + ".counts",
        raw_prefix + ".condf",
        raw_prefix + ".samples",
        raw_prefix + ".species_counts",
        raw_prefix + ".species_condf",
        raw_prefix + ".species_samples",
        filtered_bc,
        el,
        species_el,
        PANEL_METADATA,
    ]
    val_lines = []
    for vf in val_files:
        val_lines.append(
            f'if [[ ! -f "{vf}" ]]; then\n'
            f'    echo "❌ ERROR: Required file missing: {vf}"\n'
            f'    exit 1\n'
            f'fi')
    validation_block = "\n".join(val_lines)

    # Skip check
    skip_block = ""
    if not force:
        skip_block = (
            f'if [[ -f "{out_indiv}" ]] && [[ -f "{out_species}" ]]; then\n'
            f'    echo "✅ Empty drops outputs already exist for lib{lib_num}"\n'
            f'    echo "Skipping (use --force to rerun)"\n'
            f'    exit 0\n'
            f'fi\n'
        )

    # Decompress filtered barcodes for stable/reproducible input path; tet_ambient_profile can now also read gzipped barcode files directly
    decompress_block = (
        f'echo "Decompressing filtered barcodes..."\n'
        f'zcat "{filtered_bc}" > "{decompressed_bc}"\n'
        f'echo "  Wrote {decompressed_bc}"'
    )

    # tet_ambient_profile uses --ids/-i for the same expected-lines restriction
    # that demux_parallel receives via -I.  Do not pass uppercase -I here unless
    # the C++ CLI is extended to accept it; current production builds expose
    # --ids/-i for this tool.
    cmd_indiv = (
        f'{binary} -o "{raw_prefix}" '
        f'--filtered_barcodes "{decompressed_bc}" '
        f'--interindividual '
        f'--condf "{condf_interindiv}" '
        f'--ids "{el}" '
        f'-T "${{SLURM_CPUS_PER_TASK}}"'
    )

    cmd_species = (
        f'{binary} -o "{raw_prefix}" '
        f'--filtered_barcodes "{decompressed_bc}" '
        f'--interspecies '
        f'--condf "{condf_species}" '
        f'--ids "{species_el}" '
        f'--panel_metadata "{PANEL_METADATA}" '
        f'-T "${{SLURM_CPUS_PER_TASK}}"'
    )

    script = f"""#!/bin/bash
#SBATCH --job-name={job_name}
#SBATCH --output={log_dir}/{job_name}_%j.out
#SBATCH --error={log_dir}/{job_name}_%j.err
#SBATCH --time={SLURM_TIME}
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --partition={SLURM_PARTITION}
#SBATCH --nodes=1
{dep_line}

set -euo pipefail
{module_block()}

echo "================================================================"
echo "STAGE 3 EMPTY_DROPS: Empty drops estimation lib{lib_num}"
echo "  Raw prefix: {raw_prefix}"
echo "  Individual expected lines: {el}"
echo "  Species expected lines: {species_el}"
echo "  Started: $(date)"
echo "  Node: $(hostname)"
echo "================================================================"
echo ""

# Repair missing symlinks if DEMUX crashed before its symlink step
{symlink_repair_block}

echo ""
# Validate inputs
{validation_block}

# Skip if outputs exist
{skip_block}
# Decompress filtered barcodes for a stable local input path
{decompress_block}

echo ""
echo "=== Interindividual ambient profile ==="
echo "Running: {cmd_indiv}"
echo ""

{cmd_indiv}

echo ""
echo "=== Interspecies ambient profile ==="
echo "Running: {cmd_species}"
echo ""

{cmd_species}

echo ""
# Validate outputs
MISSING=0
if [[ ! -f "{out_indiv}" ]]; then
    echo "❌ ERROR: Expected output missing: {out_indiv}"
    MISSING=$((MISSING + 1))
fi
if [[ ! -f "{out_species}" ]]; then
    echo "❌ ERROR: Expected output missing: {out_species}"
    MISSING=$((MISSING + 1))
fi
if [[ "${{MISSING}}" -gt 0 ]]; then
    echo "❌ ERROR: lib{lib_num} empty drops failed: ${{MISSING}} output files missing"
    exit 1
fi

echo "✅ lib{lib_num} empty drops completed successfully"
echo "Finished: $(date)"
"""

    os.makedirs(log_dir, exist_ok=True)
    os.makedirs(script_dir, exist_ok=True)
    script_path = os.path.join(script_dir, f"{job_name}.sbatch")
    with open(script_path, "w") as f:
        f.write(script)
    os.chmod(script_path, stat.S_IRWXU | stat.S_IRGRP | stat.S_IXGRP |
             stat.S_IROTH | stat.S_IXOTH)
    return script_path


# =============================================================================
# Stage 4 CONTAM: Contamination estimation (per condition x library)
# =============================================================================

def build_contam_command(lib_num, cond, assignment_source="demux",
                         out_prefix_override=None, fixed_ambient=None,
                         input_overrides=None):
    """Build the selected contamination-estimator command for a condition/library.

    Returns the full command string.
    """
    mode = cond["mode"]
    binary_name = "legacy2c_contam_estimate" if mode == LEGACY2C_MODE else "tet_contam_estimate"
    binary = os.path.join(SOFTWARE_BIN, binary_name)
    prefix = get_demux_prefix(lib_num)
    arm_inputs = dict(input_overrides) if input_overrides else (
        identity_ambient_arm_inputs(lib_num, assignment_source)
        if assignment_source in IDENTITY_AMBIENT_ARMS else None
    )
    if assignment_source in IDENTITY_AMBIENT_ARMS and arm_inputs is None:
        raise ValueError(
            f"{assignment_source} is not applicable to lib{lib_num}")
    el = (
        arm_inputs["receiver_lines"] if arm_inputs else
        get_expected_lines(lib_num)
    )
    individual_candidates = (
        arm_inputs["ambient_candidates"] if arm_inputs else
        get_individual_ambient_candidates(lib_num)
        if cond["mode"] in (1, 3, "1+SR", LEGACY2C_MODE) else None
    )
    species_el = get_species_expected_lines(lib_num) if cond["mode"] in (2, 4, 5) else None
    out_prefix = (
        os.path.abspath(out_prefix_override) if out_prefix_override else
        get_contam_prefix(lib_num, cond["abbrev"], assignment_source)
    )
    if (not input_overrides and assignment_source == "reconciled" and
            cond["mode"] in (1, 3, "1+SR", LEGACY2C_MODE)):
        # Reconciliation may identify a cross-library swap whose true receiver
        # identity is absent from the original library roster. The job builds
        # this receiver-only union while leaving ambient_candidates restricted
        # to the original library donor roster.
        el = out_prefix + ".reconciled_expected_lines.txt"

    extra = cond["flags"]
    if arm_inputs:
        # Controlled A-D comparisons must hold supplied receiver identities
        # fixed while retaining the estimator's normal outer convergence.
        extra = (extra + " --freeze_assignments").strip()
    if fixed_ambient:
        extra = (
            extra + f' --fixed_ambient "{os.path.abspath(fixed_ambient)}"'
            ' --fixed_ambient_basis library'
        ).strip()

    # Resolve template placeholders in extra flags
    sc = get_species_counts(lib_num)
    extra = extra.replace("{sc}", f'"{sc}"')

    # Build base command by mode
    condf_interindiv = prefix + ".condf"
    species_condf = get_species_condf(lib_num)
    if cond.get("runner") == "geometry_gated":
        helper = resolve_geometry_gate_script()
        estimator_args = (
            f'-o "{out_prefix}" --interindividual -X "{el}" '
            f'--ambient_candidates "{individual_candidates}" '
            f'--condf "{condf_interindiv}" '
            f'--strict_condf --run_class production '
            f'--assignments_basis library --expected_lines_basis library '
            f'--ambient_candidates_basis library '
            f'--condition_key "{cond["abbrev"]}" '
            f'--run_contract "{out_prefix}.run_contract.json" '
            f'--profile_restarts {PROFILE_RESTARTS}'
        )
        if extra.strip():
            estimator_args += f" {extra}"
        estimator_args += ' -T "${SLURM_CPUS_PER_TASK}"'
        return (
            f'python3 "{helper}" '
            f'--estimator-binary "{binary}" '
            f'--condition-key "{cond["abbrev"]}" '
            f'--gate-version "{cond["geometry_gate_version"]}" '
            f'--base-strength {cond["base_source_exclusion_strength"]:.17g} '
            f'--fallback-strength {cond["fallback_source_exclusion_strength"]:.17g} '
            f'--parent-axis-alpha-threshold {cond["geometry_gate_parent_axis_alpha_min"]:.17g} '
            f'--ambient-orthogonal-norm-threshold {cond["geometry_gate_ambient_orthogonal_norm_max"]:.17g} '
            f'--parent-mass-threshold {cond["geometry_gate_parent_mass_min"]:.17g} '
            f'-- {estimator_args}'
        )
    if mode == 1:
        cmd = (f'{binary} -o "{out_prefix}" --interindividual -X "{el}" '
               f'--ambient_candidates "{individual_candidates}" '
               f'--condf "{condf_interindiv}" '
               f'--profile_restarts {PROFILE_RESTARTS}')
    elif mode == "1+SR":
        cmd = (f'{binary} -o "{out_prefix}" --interindividual -X "{el}" '
               f'--ambient_candidates "{individual_candidates}" '
               f'--condf "{condf_interindiv}" '
               f'--species_regularize --allow_legacy_species_regularize -P "{PANEL_METADATA}" '
               f'--profile_restarts {PROFILE_RESTARTS}')
    elif mode == LEGACY2C_MODE:
        variant = cond["legacy2c_variant"]
        cmd = (f'{binary} --variant {variant} '
               f'--output_prefix "{out_prefix}" '
               f'--receiver_lines "{el}" '
               f'--ambient_candidates "{individual_candidates}" '
               f'--strict_condf '
               f'--run_class production '
               f'--assignments_basis library '
               f'--expected_lines_basis library '
               f'--ambient_candidates_basis library '
               f'--run_contract "{out_prefix}.run_contract.json"')
        if variant == "tet_aware":
            cmd += f' --min_signal_gap {LEGACY2C_MIN_SIGNAL_GAP:.2f}'
    elif mode == 2:
        cmd = (f'{binary} -o "{out_prefix}" --interspecies '
               f'-X "{species_el}" '
               f'--species_condf "{species_condf}" '
               f'-P "{PANEL_METADATA}" '
               f'--profile_restarts {PROFILE_RESTARTS}')
    elif mode == 3:
        if fixed_ambient:
            # A frozen profile replaces, rather than supplements, the Mode-3
            # warm start; the estimator intentionally rejects both together.
            cmd = (f'{binary} -o "{out_prefix}" --interindividual -X "{el}" '
                   f'--ambient_candidates "{individual_candidates}" '
                   f'--condf "{condf_interindiv}" '
                   f'--profile_restarts {PROFILE_RESTARTS}')
        else:
            ed_ind = get_empty_drops_indiv(lib_num)
            cmd = (f'{binary} -o "{out_prefix}" --interindividual -X "{el}" '
                   f'--ambient_candidates "{individual_candidates}" '
                   f'--condf "{condf_interindiv}" '
                   f'--warm_start "{ed_ind}" '
                   f'--profile_restarts {PROFILE_RESTARTS}')
    elif mode == 4:
        ed_sp = get_empty_drops_species(lib_num)
        cmd = (f'{binary} -o "{out_prefix}" --interspecies '
               f'-X "{species_el}" '
               f'--species_condf "{species_condf}" '
               f'-P "{PANEL_METADATA}" '
               f'--warm_start "{ed_sp}" '
               f'--profile_restarts {PROFILE_RESTARTS}')
    elif mode == 5:
        ed_sp = get_empty_drops_species(lib_num)
        cmd = (f'{binary} -o "{out_prefix}" --interspecies '
               f'-X "{species_el}" '
               f'--species_condf "{species_condf}" '
               f'-P "{PANEL_METADATA}" '
               f'--fixed_ambient "{ed_sp}" '
               f'--profile_restarts {PROFILE_RESTARTS}')
    else:
        raise ValueError(f"Unknown mode: {mode}")

    if cond.get("final_ck_condition") or arm_inputs:
        cmd += (
            f' --strict_condf --run_class production'
            f' --assignments_basis library'
            f' --expected_lines_basis library'
            f' --ambient_candidates_basis library'
            f' --run_contract "{out_prefix}.run_contract.json"'
        )
        if cond["mode"] == 3:
            cmd += ' --warm_start_basis library'

    if extra.strip():
        cmd += f" {extra}"

    cmd += ' -T "${SLURM_CPUS_PER_TASK}"'
    return cmd


def generate_contam_script(lib_num, cond, dep_job_id=None, force=False,
                           assignment_source="demux"):
    """Generate sbatch script for contamination estimation.

    dep_job_id should be the DEMUX stage job (for non-warm-start conditions) or
    EMPTY_DROPS stage job (for warm-start conditions).

    Returns script_path.
    """
    abbrev = cond["abbrev"]
    arm_inputs = (
        identity_ambient_arm_inputs(lib_num, assignment_source)
        if assignment_source in IDENTITY_AMBIENT_ARMS else None
    )
    if assignment_source in IDENTITY_AMBIENT_ARMS and arm_inputs is None:
        raise ValueError(
            f"cannot generate non-applicable {assignment_source} for lib{lib_num}")
    assignment_basis = (
        arm_inputs["assignment_basis"] if arm_inputs else assignment_source)
    source_suffix = "" if assignment_source == "demux" else f"_{assignment_source}"
    job_name = f"qc_{abbrev}{source_suffix}_lib{lib_num}"
    if arm_inputs:
        # Keep SLURM's display name compact while preserving the descriptive
        # output/log directory names.
        job_name = f"qc_{abbrev}_arm{arm_inputs['arm']}_lib{lib_num}"
        log_dir = os.path.join(
            get_log_dir(), abbrev, IDENTITY_AMBIENT_DIRNAME,
            _identity_ambient_candidate_set(), assignment_source)
    else:
        log_dir = os.path.join(get_log_dir(), abbrev, assignment_source)
    script_dir = get_script_dir()

    out_dir = get_contam_dir(lib_num, abbrev, assignment_source)
    out_prefix = get_contam_prefix(lib_num, abbrev, assignment_source)
    arm_contract_path = out_prefix + ".identity_ambient_arm.tsv"
    prefix = get_demux_prefix(lib_num)
    selected_assignments = get_contam_assignments_path(
        lib_num, assignment_source)
    demux_assignments = prefix + ".assignments"

    # SLURM dependency
    dep_line = ""
    if dep_job_id is not None:
        dep_line = f"#SBATCH --dependency=afterok:{dep_job_id}"

    # Build validation block for upstream files
    val_files = []
    if cond["mode"] in (2, 4, 5):
        for ext in [".species_counts", ".species_condf", ".species_samples", ".species_assignments", ".samples"]:
            val_files.append(prefix + ext)
        val_files.append(selected_assignments)
        val_files.append(PANEL_METADATA)
    else:
        for ext in [".counts", ".condf", ".samples"]:
            val_files.append(prefix + ext)
        val_files.append(selected_assignments)
    if assignment_basis == "reconciled":
        # Reconciled exports use NA for unchanged-cell scores. The
        # contamination model requires a positive numeric confidence for every
        # row, so the job recovers the original demux LLR before launching the
        # estimator. This is also needed by native species comparison runs.
        val_files.append(demux_assignments)
    if cond["mode"] in (2, 4, 5):
        el = get_species_expected_lines(lib_num)
    else:
        el = (
            arm_inputs["receiver_lines"] if arm_inputs else
            get_expected_lines(lib_num))
    if el:
        val_files.append(el)
    if cond["mode"] in (1, 3, "1+SR", LEGACY2C_MODE):
        val_files.append(
            arm_inputs["ambient_candidates"] if arm_inputs else
            get_individual_ambient_candidates(lib_num))
    if arm_inputs:
        val_files.extend([
            arm_inputs["scrutiny_cells"], arm_inputs["context_path"],
        ])

    if cond["mode"] == "1+SR":
        val_files.append(PANEL_METADATA)
        if cond["needs_species_counts"]:
            val_files.append(get_species_counts(lib_num))
            val_files.append(get_species_condf(lib_num))
            val_files.append(prefix + ".species_samples")
    if cond["needs_empty_drops"] == "individual":
        val_files.append(get_empty_drops_indiv(lib_num))
    elif cond["needs_empty_drops"] in ("species", "species_fixed"):
        val_files.append(get_empty_drops_species(lib_num))

    val_lines = []
    for vf in val_files:
        val_lines.append(
            f'if [[ ! -f "{vf}" ]]; then\n'
            f'    echo "❌ ERROR: Required upstream file missing: {vf}"\n'
            f'    exit 1\n'
            f'fi')
    validation_block = "\n".join(val_lines)

    # Skip check
    skip_block = ""
    early_arm_completion_block = ""
    if not force:
        if cond.get("runner") == "geometry_gated":
            required = (
                ".contam_rate", ".contam_prof", ".allele_ratio",
                ".contam_diagnostics.tsv", ".profile_fit_diagnostics.tsv",
                ".condf_coverage.tsv", ".run_contract.json",
                ".geometry_gate_audit.tsv",
            )
            tests = " && ".join(
                f'[[ -s "{out_prefix}{suffix}" ]]' for suffix in required)
            skip_block = (
                f'if {tests}; then\n'
                f'    echo "✅ Complete geometry-gated outputs already exist: {out_prefix}"\n'
                f'    echo "Skipping (use --force to rerun)"\n'
                f'    exit 0\n'
                f'fi\n'
            )
        elif cond.get("final_ck_condition"):
            required = (
                ".contam_rate", ".contam_prof", ".allele_ratio",
                ".contam_diagnostics.tsv", ".profile_fit_diagnostics.tsv",
                ".condf_coverage.tsv", ".run_contract.json",
            )
            tests = " && ".join(
                f'[[ -s "{out_prefix}{suffix}" ]]' for suffix in required)
            skip_block = (
                f'if {tests}; then\n'
                f'    echo "✅ Complete final CK outputs already exist: {out_prefix}"\n'
                f'    echo "Skipping (use --force to rerun)"\n'
                f'    exit 0\n'
                f'fi\n'
            )
        elif cond["mode"] == LEGACY2C_MODE:
            skip_block = (
                f'if [[ -f "{out_prefix}.contam_rate" ]] && '
                f'[[ -f "{out_prefix}.contam_prof" ]] && '
                f'[[ -f "{out_prefix}.decontam.assignments" ]] && '
                f'[[ -f "{out_prefix}.legacy2c_candidates.txt" ]] && '
                f'[[ -f "{out_prefix}.condf_coverage.tsv" ]] && '
                f'[[ -f "{out_prefix}.profile_fit_diagnostics.tsv" ]] && '
                f'[[ -f "{out_prefix}.legacy2c_diagnostics.tsv" ]] && '
                f'[[ -f "{out_prefix}.run_contract.json" ]]; then\n'
                f'    echo "✅ Complete legacy2c outputs already exist: {out_prefix}"\n'
                f'    echo "Skipping (use --force to rerun)"\n'
                f'    exit 0\n'
                f'fi\n'
            )
        elif cond["mode"] in (2, 4, 5, "1+SR"):
            skip_block = (
                f'if [[ -f "{out_prefix}.contam_rate" ]] && '
                f'([[ -f "{out_prefix}.species_prof" ]] || [[ -f "{out_prefix}.contam_prof" ]]); then\n'
                f'    echo "✅ Outputs already exist: {out_prefix}.contam_rate and native species profile"\n'
                f'    echo "Skipping (use --force to rerun)"\n'
                f'    exit 0\n'
                f'fi\n'
            )
        else:
            skip_block = (
                f'if [[ -f "{out_prefix}.contam_rate" ]] && '
                f'[[ -f "{out_prefix}.contam_prof" ]]; then\n'
                f'    echo "✅ Outputs already exist: {out_prefix}.contam_rate and .contam_prof"\n'
                f'    echo "Skipping (use --force to rerun)"\n'
                f'    exit 0\n'
                f'fi\n'
            )
    if arm_inputs and not force:
        required = [
            ".contam_rate", ".contam_prof", ".allele_ratio",
            ".contam_diagnostics.tsv", ".cell_source_profile.tsv",
            ".profile_fit_diagnostics.tsv", ".condf_coverage.tsv",
            ".run_contract.json", ".decontam.assignments",
            ".identity_ambient_arm.tsv",
        ]
        if cond.get("runner") == "geometry_gated":
            required.append(".geometry_gate_audit.tsv")
        tests = " && ".join(
            f'[[ -s "{out_prefix}{suffix}" ]]' for suffix in required)
        any_artifact = " || ".join(
            f'[[ -e "{out_prefix}{suffix}" ]]' for suffix in required)
        early_arm_completion_block = f'''
if {tests}; then
    if awk -F '\t' -v candidate_set='{arm_inputs['candidate_set']}' -v fingerprint='{arm_inputs['plan_fingerprint']}' '
        NR == 1 {{ for (i=1; i<=NF; i++) column[$i]=i; next }}
        NR == 2 {{
            rows++
            ok = (column["candidate_set"] && column["plan_fingerprint"] &&
                  column["assignment_update_mode"] &&
                  column["assignment_score_basis"] &&
                  $(column["candidate_set"]) == candidate_set &&
                  $(column["plan_fingerprint"]) == fingerprint &&
                  $(column["assignment_update_mode"]) == "iterative_frozen" &&
                  $(column["assignment_score_basis"]) == "original_demux_all_arms")
        }}
        NR > 2 {{ rows++ }}
        END {{ exit !(ok && rows == 1) }}
    ' "{arm_contract_path}"; then
        echo "✅ Complete four-arm output already exists for plan {arm_inputs['plan_fingerprint']}: {out_prefix}"
        echo "Skipping (use --force to rerun)"
        exit 0
    fi
    echo "❌ ERROR: existing four-arm output belongs to a different comparison plan" >&2
    echo "Current plan: {arm_inputs['plan_fingerprint']} ({arm_inputs['candidate_set']})" >&2
    echo "Rerun with --force to replace only this isolated arm output." >&2
    exit 1
fi
if {any_artifact}; then
    echo "❌ ERROR: incomplete/stale four-arm artifacts exist without a matching completion contract: {out_prefix}" >&2
    echo "Rerun with --force to replace only this isolated arm output." >&2
    exit 1
fi
'''
        # Four-arm completion is checked before output-side assignment links
        # are created, so a dry reuse cannot mutate a completed result bundle.
        skip_block = ""

    # Expose only the native input bundle needed by the selected mode. A
    # reconciliation export cannot be linked directly because unchanged rows
    # may carry "NA" in the numeric LLR column. Four-arm comparison overlays
    # already carry the original demux LLR for every cell; legacy reconciled
    # runs retain their historical score-recovery behavior.
    if assignment_basis == "reconciled":
        assignment_input_block = f'''prepared_assignments="{out_prefix}.assignments"
prepared_tmp="${{prepared_assignments}}.tmp.${{SLURM_JOB_ID:-$$}}"
rm -f "${{prepared_tmp}}"
if ! awk '
BEGIN {{ FS="[[:space:]]+"; OFS="\\t" }}
function numeric(x) {{
    return x ~ /^[-+]?([0-9]+([.][0-9]*)?|[.][0-9]+)([eE][-+]?[0-9]+)?$/
}}
NR == FNR {{
    if (NF >= 4 && numeric($4) && ($4 + 0) > 0) original_llr[$1] = $4
    next
}}
NF != 4 {{
    printf "ERROR: reconciled assignment line %d has %d columns; expected 4\\n", FNR, NF > "/dev/stderr"
    bad = 1
    next
}}
{{
    score = $4
    if (!numeric(score) || (score + 0) <= 0) {{
        if (!($1 in original_llr)) {{
            printf "ERROR: no positive original demux LLR for reconciled barcode %s (line %d)\\n", $1, FNR > "/dev/stderr"
            bad = 1
            next
        }}
        score = original_llr[$1]
        inherited++
    }} else {{
        reconciled_score++
    }}
    # The compatibility export can use A+A to retain a post-hoc biological
    # homotetraploid classification. CellBouncer assignments encode donor
    # genotype, not biological ploidy, so pass A/S rather than A+A/D. Preserve
    # true heterotypic A+B identities as D.
    model_identity = $2
    n_components = split($2, identity_components, "+")
    if (n_components == 2 && identity_components[1] == identity_components[2]) {{
        model_identity = identity_components[1]
    }}
    assignment_type = (index(model_identity, "+") > 0) ? "D" : "S"
    if ($2 != model_identity || $3 != assignment_type) normalized_type++
    print $1, model_identity, assignment_type, score
    written++
}}
END {{
    if (written == 0) {{
        print "ERROR: no usable reconciled assignment rows were produced" > "/dev/stderr"
        bad = 1
    }}
    if (bad) exit 42
    printf "Prepared %d reconciled assignments: %d reconciliation scores, %d inherited demux LLRs, %d genotype-type normalizations\\n", written, reconciled_score, inherited, normalized_type > "/dev/stderr"
}}
' "{demux_assignments}" "{selected_assignments}" > "${{prepared_tmp}}"; then
    rm -f "${{prepared_tmp}}"
    exit 1
fi
mv -f "${{prepared_tmp}}" "${{prepared_assignments}}"'''
        if (assignment_source == "reconciled" and
                cond["mode"] in (1, 3, "1+SR", LEGACY2C_MODE)):
            reconciled_expected = out_prefix + ".reconciled_expected_lines.txt"
            assignment_input_block += f'''
prepared_expected="{reconciled_expected}"
expected_tmp="${{prepared_expected}}.tmp.${{SLURM_JOB_ID:-$$}}"
rm -f "${{expected_tmp}}"
{{
    awk 'NF && $1 !~ /^#/ {{ print $1 }}' "{el}"
    awk 'NF >= 2 {{ print $2 }}' "${{prepared_assignments}}"
}} | awk 'NF && !seen[$0]++' > "${{expected_tmp}}"
if [[ ! -s "${{expected_tmp}}" ]]; then
    echo "ERROR: reconciled receiver expected-lines roster is empty" >&2
    rm -f "${{expected_tmp}}"
    exit 1
fi
mv -f "${{expected_tmp}}" "${{prepared_expected}}"
echo "Prepared reconciled receiver roster: $(wc -l < "${{prepared_expected}}") identities" >&2'''
    else:
        assignment_input_block = (
            f'ln -sf "{selected_assignments}" "{out_prefix}.assignments"')

    symlink_cmds = []
    if cond["mode"] in (2, 4, 5):
        for ext in [".species_assignments", ".species_counts", ".species_condf", ".species_samples"]:
            src = prefix + ext
            dst = out_prefix + ext
            symlink_cmds.append(f'ln -sf "{src}" "{dst}"')
        # Also expose the original individual assignments/sample list to
        # tet_contam_estimate.  Native species mode uses these only to derive
        # weighted endogenous species composition vectors, e.g.
        # Chinobo-mCherry+JOS3C1 -> B:0.25,C:0.25,O:0.5.
        symlink_cmds.append(assignment_input_block)
        src = prefix + ".samples"
        dst = out_prefix + ".samples"
        symlink_cmds.append(f'ln -sf "{src}" "{dst}"')
    else:
        symlink_cmds.append(assignment_input_block)
        for ext in [".counts", ".condf", ".samples"]:
            src = prefix + ext
            dst = out_prefix + ext
            symlink_cmds.append(f'ln -sf "{src}" "{dst}"')
        if cond["mode"] == "1+SR":
            for ext in [".species_counts", ".species_condf", ".species_samples"]:
                src = prefix + ext
                dst = out_prefix + ext
                symlink_cmds.append(f'ln -sf "{src}" "{dst}"')
    symlink_block = "\n".join(symlink_cmds)

    command = build_contam_command(
        lib_num, cond, assignment_source=assignment_source)

    # --force must remove stale estimator outputs before launching.
    # Both estimator binaries refuse to overwrite existing .contam_prof/.contam_rate
    # when it is running from already-counted demux files, so orchestrator-level
    # force has to clean only output artifacts while preserving input symlinks
    # such as .counts/.condf/.samples/.species_counts/.species_condf.
    force_cleanup_block = ""
    if force:
        cleanup_targets = [
            out_prefix + ".contam_prof",
            out_prefix + ".contam_rate",
            out_prefix + ".species_prof",
            out_prefix + ".allele_ratios",
            out_prefix + ".allele_ratios.tsv",
            out_prefix + ".allele_ratios.tsv.gz",
            out_prefix + ".cell_contam.tsv",
            out_prefix + ".cell_contam.tsv.gz",
            out_prefix + ".diagnostics.tsv",
            out_prefix + ".diagnostics.tsv.gz",
            out_prefix + ".pass1.contam_prof",
            out_prefix + ".pass1.contam_rate",
            out_prefix + ".pass1.decontam.assignments",
            out_prefix + ".pass1.species_prof",
            out_prefix + ".model_fit.tsv",
            out_prefix + ".profile_fit_diagnostics.tsv",
            out_prefix + ".condf_coverage.tsv",
            out_prefix + ".contam_diagnostics.tsv",
            out_prefix + ".contam_diagnostics.tsv.gz",
            out_prefix + ".class_residuals.tsv",
            out_prefix + ".class_residuals.tsv.gz",
            out_prefix + ".cell_source_profile.tsv",
            out_prefix + ".cell_source_profile.tsv.gz",
            out_prefix + ".run_contract.json",
            out_prefix + ".decontam.assignments",
            out_prefix + ".legacy2c_candidates.txt",
            out_prefix + ".legacy2c_diagnostics.tsv",
            out_prefix + ".geometry_gate_audit.tsv",
            out_prefix + ".allele_ratio",
            out_prefix + ".identity_ambient_arm.tsv",
        ]
        quoted = " ".join(f'"{x}"' for x in cleanup_targets)
        endpoint_cleanup = ""
        if cond.get("runner") == "geometry_gated":
            endpoint_cleanup = (
                f'rm -f "{out_prefix}.geometry_base_endpoint".*\n'
                f'rm -f "{out_prefix}.geometry_fallback_endpoint".*\n'
            )
        force_cleanup_block = (
            'echo "=== FORCE CLEANUP: stale contamination-estimator outputs ==="\n'
            f'rm -f {quoted}\n'
            f'{endpoint_cleanup}'
        )

    # Output validation
    expected_outputs = [out_prefix + ".contam_rate", out_prefix + ".contam_prof"]
    if cond.get("runner") == "geometry_gated":
        expected_outputs.extend([
            out_prefix + ".allele_ratio",
            out_prefix + ".contam_diagnostics.tsv",
            out_prefix + ".profile_fit_diagnostics.tsv",
            out_prefix + ".condf_coverage.tsv",
            out_prefix + ".run_contract.json",
            out_prefix + ".geometry_gate_audit.tsv",
        ])
    if cond.get("final_ck_condition") and cond.get("runner") != "geometry_gated":
        expected_outputs.extend([
            out_prefix + ".allele_ratio",
            out_prefix + ".contam_diagnostics.tsv",
            out_prefix + ".profile_fit_diagnostics.tsv",
            out_prefix + ".condf_coverage.tsv",
            out_prefix + ".run_contract.json",
        ])
    if cond["mode"] == LEGACY2C_MODE:
        expected_outputs.extend([
            out_prefix + ".decontam.assignments",
            out_prefix + ".legacy2c_candidates.txt",
            out_prefix + ".condf_coverage.tsv",
            out_prefix + ".profile_fit_diagnostics.tsv",
            out_prefix + ".legacy2c_diagnostics.tsv",
            out_prefix + ".run_contract.json",
        ])
    if cond["mode"] in (2, 4, 5, "1+SR"):
        expected_outputs.append(out_prefix + ".species_prof")
    if arm_inputs:
        expected_outputs.extend([
            out_prefix + ".allele_ratio",
            out_prefix + ".contam_diagnostics.tsv",
            out_prefix + ".cell_source_profile.tsv",
            out_prefix + ".profile_fit_diagnostics.tsv",
            out_prefix + ".condf_coverage.tsv",
            out_prefix + ".run_contract.json",
            out_prefix + ".decontam.assignments",
            arm_contract_path,
        ])
        if cond.get("runner") == "geometry_gated":
            expected_outputs.append(out_prefix + ".geometry_gate_audit.tsv")
        expected_outputs = list(dict.fromkeys(expected_outputs))
    out_check_lines = []
    output_test = "-s" if arm_inputs else "-f"
    for ef in expected_outputs:
        out_check_lines.append(
            f'if [[ ! {output_test} "{ef}" ]]; then\n'
            f'    echo "❌ ERROR: Expected output missing or empty: {ef}"\n'
            f'    MISSING=$((MISSING + 1))\n'
            f'fi')
    out_check_block = "\n".join(out_check_lines)

    info = LIB_INFO.get(lib_num, {})
    mode_desc = {1: "Interindividual", "1+SR": "Interindividual+SR",
                 LEGACY2C_MODE: "Legacy two-component compatibility",
                 2: "Interspecies", 3: "Interindividual+WS",
                 4: "Interspecies+WS", 5: "Interspecies+FixedEmpty"}
    if cond.get("runner") == "geometry_gated":
        mode_desc[1] = "Interindividual CK geometry-gated"

    # For native species outputs, create an explicit .species_prof alias after
    # tet_contam_estimate writes the native profile to .contam_prof.
    species_alias_block = ""
    if cond["mode"] in (2, 4, 5, "1+SR"):
        species_alias_block = species_profile_alias_block(out_prefix)

    arm_plan_guard_function = ""
    arm_plan_guard_call = ""
    arm_validation_block = ""
    arm_echo_block = ""
    if arm_inputs:
        assignment_context_key = (
            "demux_assignments"
            if arm_inputs["assignment_basis"] == "demux" else
            "reconciled_assignments")
        if assignment_source in {"demux_original", "demux_augmented"}:
            receiver_context_key = "original_receiver_lines"
        elif assignment_source == "reconciled_augmented":
            receiver_context_key = "augmented_receiver_lines"
        else:
            receiver_context_key = "replacement_receiver_lines"
        candidate_context_key = (
            "original_ambient_candidates"
            if assignment_source == "demux_original" else
            "augmented_ambient_candidates"
            if assignment_source in {
                "demux_augmented", "reconciled_augmented"} else
            "replacement_ambient_candidates")
        arm_plan_guard_function = f'''
validate_identity_ambient_plan() {{
    python3 - "{arm_inputs['context_path']}" \
        "{arm_inputs['candidate_set']}" "{arm_inputs['plan_fingerprint']}" \
        "{assignment_context_key}" "{selected_assignments}" \
        "{receiver_context_key}" "{arm_inputs['receiver_lines']}" \
        "{candidate_context_key}" "{arm_inputs['ambient_candidates']}" \
        "{arm_inputs['scrutiny_cells']}" <<'PY'
import hashlib
import json
import os
import sys

(context_path, candidate_set, fingerprint, assignment_key, assignments,
 receiver_key, receivers, candidate_key, candidates, scrutiny) = sys.argv[1:]
with open(context_path, "r", encoding="utf-8") as handle:
    context = json.load(handle)
if context.get("candidate_set") != candidate_set:
    raise SystemExit("ERROR: four-arm context candidate set changed after submission")
if context.get("plan_fingerprint") != fingerprint:
    raise SystemExit("ERROR: four-arm context fingerprint changed after submission")
for key, expected in (
    (assignment_key, assignments),
    (receiver_key, receivers),
    (candidate_key, candidates),
    ("scrutiny_cells", scrutiny),
):
    observed = str(context.get(key, ""))
    if os.path.abspath(observed) != os.path.abspath(expected):
        raise SystemExit(
            f"ERROR: four-arm context {{key}} path changed: "
            f"{{observed!r}} != {{expected!r}}")

def digest(path):
    value = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()

for group in ("input_records", "plan_artifact_records"):
    records = context.get(group)
    if not isinstance(records, dict) or not records:
        raise SystemExit(f"ERROR: four-arm context lacks {{group}}")
    for name, record in records.items():
        if not isinstance(record, dict):
            raise SystemExit(f"ERROR: malformed {{group}} record {{name}}")
        path = str(record.get("path", ""))
        if (not os.path.isfile(path) or
                os.path.getsize(path) != int(record.get("size_bytes", -1)) or
                digest(path) != record.get("sha256")):
            raise SystemExit(
                f"ERROR: four-arm plan input/artifact changed: {{name}}: {{path}}")
PY
}}
'''
        arm_plan_guard_call = "validate_identity_ambient_plan"
        arm_echo_block = (
            f'echo "  Comparison arm: {arm_inputs["arm"]} '
            f'({arm_inputs["label"]})"\n'
            f'echo "  Receiver roster: {arm_inputs["receiver_lines"]}"\n'
            f'echo "  Ambient donor roster: {arm_inputs["ambient_candidates"]}"'
        )
        arm_validation_block = f'''
validate_identity_ambient_plan
echo "=== FOUR-ARM CONTRACT VALIDATION ==="
comparison_tmp=$(mktemp -d "${{SLURM_TMPDIR:-/tmp}}/identity-ambient-${{SLURM_JOB_ID:-$$}}.XXXXXX")
trap 'rm -rf "$comparison_tmp"' EXIT
awk 'NF >= 2 {{print $1 "\\t" $2}}' "{out_prefix}.assignments" | LC_ALL=C sort > "$comparison_tmp/input.tsv"
awk 'NF >= 2 {{print $1 "\\t" $2}}' "{out_prefix}.decontam.assignments" | LC_ALL=C sort > "$comparison_tmp/final.tsv"
if ! cmp -s "$comparison_tmp/input.tsv" "$comparison_tmp/final.tsv"; then
    echo "❌ ERROR: --freeze_assignments invariant failed; output identities differ from the supplied arm assignments" >&2
    diff -u "$comparison_tmp/input.tsv" "$comparison_tmp/final.tsv" | head -80 >&2 || true
    exit 1
fi
if ! awk '
NR == FNR {{ if (NF && $1 !~ /^#/) allowed[$1] = 1; next }}
NF >= 2 && !($1 in allowed) {{
    printf "ERROR: contamination profile donor %s is absent from arm roster\\n", $1 > "/dev/stderr"
    bad = 1
}}
NF >= 2 {{ observed[$1] = 1 }}
END {{
    for (source in allowed) if (!(source in observed)) {{
        printf "ERROR: arm donor %s is absent from contamination profile\\n", source > "/dev/stderr"
        bad = 1
    }}
    exit bad
}}
' "{arm_inputs['ambient_candidates']}" "{out_prefix}.contam_prof"; then
    exit 1
fi
python3 - "{out_prefix}.run_contract.json" "{out_prefix}.assignments" "{arm_inputs['receiver_lines']}" "{arm_inputs['ambient_candidates']}" "{demux_assignments}" "{out_prefix}.contam_rate" <<'PY'
import json
import math
import os
import sys

contract_path, assignments, receivers, candidates, demux_assignments, rate_path = sys.argv[1:]
with open(contract_path, "r", encoding="utf-8") as handle:
    contract = json.load(handle)
endpoint = contract
if endpoint.get("assignment_update_mode") != "iterative_frozen":
    raise SystemExit(
        "ERROR: run contract does not prove iterative frozen assignments: "
        + repr(endpoint.get("assignment_update_mode")))
if endpoint.get("freeze_assignments") is not True:
    raise SystemExit("ERROR: run contract lacks freeze_assignments=true")
if endpoint.get("production_contract_pass") is not True:
    raise SystemExit("ERROR: run contract lacks production_contract_pass=true")
for basis_field in (
        "assignments_basis", "expected_lines_basis",
        "ambient_candidates_basis"):
    if endpoint.get(basis_field) != "library":
        raise SystemExit(
            f"ERROR: run contract {{basis_field}} is not library: "
            + repr(endpoint.get(basis_field)))
def recorded_path(name):
    value = endpoint.get(name)
    return value.get("path", "") if isinstance(value, dict) else str(value or "")
for name, expected in (
    ("assignments", assignments),
    ("expected_lines", receivers),
    ("ambient_candidates", candidates),
):
    observed = recorded_path(name)
    if os.path.abspath(observed) != os.path.abspath(expected):
        raise SystemExit(
            f"ERROR: run contract {{name}} path mismatch: {{observed!r}} != {{expected!r}}")

def assignment_scores(path):
    result = {{}}
    with open(path, "r", encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, 1):
            parts = raw.split()
            if not parts:
                continue
            if len(parts) != 4:
                raise SystemExit(
                    f"ERROR: malformed assignment row {{line_number}} in {{path}}")
            try:
                score = float(parts[3])
            except ValueError as exc:
                raise SystemExit(
                    f"ERROR: nonnumeric assignment score for {{parts[0]}} in {{path}}") from exc
            if not math.isfinite(score):
                raise SystemExit(
                    f"ERROR: nonfinite assignment score for {{parts[0]}} in {{path}}")
            result[parts[0]] = score
    return result

arm_scores = assignment_scores(assignments)
demux_scores = assignment_scores(demux_assignments)
if set(arm_scores) != set(demux_scores):
    raise SystemExit("ERROR: four-arm assignment barcode universe differs from demux")
for barcode, expected_score in demux_scores.items():
    if not math.isclose(
            arm_scores[barcode], expected_score, rel_tol=1e-12, abs_tol=1e-12):
        raise SystemExit(
            f"ERROR: four-arm LLR weight differs from demux for {{barcode}}: "
            f"{{arm_scores[barcode]}} != {{expected_score}}")

rate_barcodes = set()
with open(rate_path, "r", encoding="utf-8") as handle:
    for line_number, raw in enumerate(handle, 1):
        parts = raw.split()
        if not parts:
            continue
        if len(parts) < 2:
            raise SystemExit(
                f"ERROR: malformed contamination-rate row {{line_number}} in {{rate_path}}")
        if parts[0] in rate_barcodes:
            raise SystemExit(
                f"ERROR: duplicate contamination-rate barcode {{parts[0]}}")
        try:
            rate = float(parts[1])
        except ValueError as exc:
            raise SystemExit(
                f"ERROR: nonnumeric contamination rate for {{parts[0]}}") from exc
        if not math.isfinite(rate):
            raise SystemExit(
                f"ERROR: nonfinite contamination rate for {{parts[0]}}")
        rate_barcodes.add(parts[0])
if rate_barcodes != set(arm_scores):
    missing = sorted(set(arm_scores) - rate_barcodes)[:10]
    extra = sorted(rate_barcodes - set(arm_scores))[:10]
    raise SystemExit(
        f"ERROR: contamination-rate barcode universe mismatch; "
        f"missing={{missing}}, extra={{extra}}")
PY
arm_tmp="{arm_contract_path}.tmp.${{SLURM_JOB_ID:-$$}}"
printf 'library\tcondition\tarm\tarm_key\tassignment_basis\troster_basis\tcandidate_set\tplan_fingerprint\tassignment_path\treceiver_lines\tambient_candidates\tscrutiny_cells\tcontext_path\tassignment_update_mode\tassignment_score_basis\n' > "$arm_tmp"
printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \\
    'lib{lib_num}' '{abbrev}' '{arm_inputs['arm']}' '{assignment_source}' \\
    '{arm_inputs['assignment_basis']}' '{arm_inputs['roster_basis']}' \\
    '{arm_inputs['candidate_set']}' '{arm_inputs['plan_fingerprint']}' \\
    '{selected_assignments}' '{arm_inputs['receiver_lines']}' \\
    '{arm_inputs['ambient_candidates']}' '{arm_inputs['scrutiny_cells']}' \\
    '{arm_inputs['context_path']}' 'iterative_frozen' \\
    'original_demux_all_arms' >> "$arm_tmp"
mv -f "$arm_tmp" "{arm_contract_path}"
echo "✅ Four-arm identity/roster contract validated"
'''

    script = f"""#!/bin/bash
#SBATCH --job-name={job_name}
#SBATCH --output={log_dir}/{job_name}_%j.out
#SBATCH --error={log_dir}/{job_name}_%j.err
#SBATCH --time={SLURM_TIME}
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --partition={SLURM_PARTITION}
#SBATCH --nodes=1
#SBATCH --chdir={AGGREGATE_ROOT}
{dep_line}

set -euo pipefail
{module_block()}

{arm_plan_guard_function}

echo "================================================================"
echo "STAGE 4 CONTAM: Contamination estimation"
echo "  Condition: {abbrev}"
echo "  Assignment source: {assignment_source}"
echo "  Assignment basis: {assignment_basis}"
{arm_echo_block}
echo "  Assignments: {selected_assignments}"
echo "  Mode: {mode_desc.get(cond['mode'], cond['mode'])}"
echo "  Library: lib{lib_num} (Batch {info.get('batch', '?')}, {info.get('lib_type', '?')})"
echo "  Started: $(date)"
echo "  Node: $(hostname)"
echo "  Output prefix: {out_prefix}"
echo "================================================================"
echo ""

# Create output directory and symlinks
mkdir -p "{out_dir}"
{arm_plan_guard_call}
{early_arm_completion_block}
{symlink_block}

# Validate upstream files
{validation_block}

# Skip if output exists
{skip_block}
{force_cleanup_block}
echo "Running: {command}"
echo ""

{command}

{species_alias_block}

{arm_validation_block}

echo ""
echo "=== OUTPUT VALIDATION ==="
MISSING=0
{out_check_block}
if [[ "${{MISSING}}" -gt 0 ]]; then
    echo "❌ ERROR: {abbrev}/lib{lib_num} failed: ${{MISSING}} output files missing"
    exit 1
fi

echo "✅ {abbrev}/lib{lib_num} completed successfully"
echo "Finished: $(date)"
"""

    os.makedirs(log_dir, exist_ok=True)
    os.makedirs(script_dir, exist_ok=True)
    script_path = os.path.join(script_dir, f"{job_name}.sbatch")
    with open(script_path, "w") as f:
        f.write(script)
    os.chmod(script_path, stat.S_IRWXU | stat.S_IRGRP | stat.S_IXGRP |
             stat.S_IROTH | stat.S_IXOTH)
    return script_path


# =============================================================================
# Stage 5 PLOIDY_NN: optional NN ploidy inference
# =============================================================================

def _lib_range_for_ploidy_nn(lib_nums):
    """Return the comma-separated syntax expected by run_ploidy_nn_inference.py."""
    return ",".join(str(x) for x in sorted(set(lib_nums)))


def generate_ploidy_nn_script(lib_nums, force=False, args=None):
    """Generate one sbatch script for NN ploidy inference over requested libraries."""
    libs = sorted(set(lib_nums))
    label = f"lib{libs[0]}" if len(libs) == 1 else f"lib{libs[0]}_{libs[-1]}"
    job_name = f"ploidy_nn_{label}"
    log_dir = os.path.join(get_log_dir(), "PLOIDY_NN")
    script_dir = get_script_dir()

    helper = resolve_process_script("run_ploidy_nn_inference.py")
    h5ad = args.ploidy_nn_h5ad if args else PLOIDY_NN_H5AD
    weights = args.ploidy_nn_weights if args else PLOIDY_NN_WEIGHTS
    scaler = weights[:-3] + "_scaler.npz" if weights.endswith(".pt") else weights + "_scaler.npz"
    out_dir = args.ploidy_calls_root if args else get_ploidy_calls_root()
    ploidy_module = args.ploidy_nn_module if args else PLOIDY_NN_MODULE
    qc_only = bool(args.ploidy_nn_qc_only) if args else False
    lib_range = _lib_range_for_ploidy_nn(libs)

    cmd = (
        f'python3 "{helper}" '
        f'--h5ad "{h5ad}" '
        f'--weights "{weights}" '
        f'--output_dir "{out_dir}" '
        f'--lib_range "{lib_range}"'
    )
    if qc_only:
        cmd += " --qc_only"
    if force:
        cmd += " --force"

    output_checks = "\n".join(
        f'if [[ ! -f "{os.path.join(out_dir, f"lib{lib}.ploidy_calls_nn.tsv")}" ]]; then\n'
        f'    echo "❌ ERROR: Expected ploidy output missing: {os.path.join(out_dir, f"lib{lib}.ploidy_calls_nn.tsv")}"\n'
        f'    MISSING=$((MISSING + 1))\n'
        f'fi'
        for lib in libs
    )

    ploidy_module_block = ""
    if ploidy_module:
        ploidy_module_block = f"""
# Load the module prerequisite and then the cluster-provided ploidy inference
# module. Do not activate conda here; Lmod requires miniforge/3 as a parent
# module for ploidy-inference/latest on this cluster.
module load miniforge/3
module load "{ploidy_module}"
"""

    script = f"""#!/bin/bash
#SBATCH --job-name={job_name}
#SBATCH --output={log_dir}/{job_name}_%j.out
#SBATCH --error={log_dir}/{job_name}_%j.err
#SBATCH --time={SLURM_TIME}
#SBATCH --cpus-per-task=16
#SBATCH --mem=256G
#SBATCH --partition={SLURM_PARTITION}
#SBATCH --nodes=1

set -euo pipefail
module purge
{ploidy_module_block}

mkdir -p "{log_dir}" "{out_dir}"

echo "================================================================"
echo "STAGE 5 PLOIDY_NN: optional homotypic/tetraploid NN ploidy inference"
echo "  Libraries: {lib_range}"
echo "  h5ad: {h5ad}"
echo "  weights: {weights}"
echo "  output_dir: {out_dir}"
echo "  helper: {helper}"
echo "  Started: $(date)"
echo "  Node: $(hostname)"
echo "================================================================"

for f in "{helper}" "{h5ad}" "{weights}" "{scaler}"; do
    if [[ ! -f "$f" ]]; then
        echo "❌ ERROR: Required PLOIDY_NN input missing: $f"
        exit 1
    fi
done


echo "Running: {cmd}"
{cmd}

MISSING=0
{output_checks}
if [[ "$MISSING" -gt 0 ]]; then
    exit 1
fi

echo "✅ PLOIDY_NN completed for libraries: {lib_range}"
echo "Finished: $(date)"
"""
    os.makedirs(log_dir, exist_ok=True)
    os.makedirs(script_dir, exist_ok=True)
    script_path = os.path.join(script_dir, f"{job_name}.sbatch")
    with open(script_path, "w") as f:
        f.write(script)
    os.chmod(script_path, stat.S_IRWXU | stat.S_IRGRP | stat.S_IXGRP |
             stat.S_IROTH | stat.S_IXOTH)
    return script_path


# =============================================================================
# Stage 6 TETRA_REFINE: optional biological assignment refinement
# =============================================================================

def generate_tetra_refine_script(lib_num, dep_job_ids=None, force=False, args=None):
    """Generate sbatch script for tetra_refine for one library."""
    job_name = f"tetra_refine_lib{lib_num}"
    log_dir = os.path.join(get_log_dir(), "TETRA_REFINE")
    script_dir = get_script_dir()
    helper = resolve_process_script("run_tetra_refine_for_library.py")
    output_root = TETRA_REFINE_ROOT
    ploidy_root = args.ploidy_calls_root if args else get_ploidy_calls_root()
    contam_condition = args.refine_contam_condition if args else REFINE_CONTAM_CONDITION
    external_min_prob = args.refine_external_ploidy_min_prob if args else REFINE_EXTERNAL_PLOIDY_MIN_PROB
    demux_prefix = get_demux_prefix(lib_num)
    expected_lines = get_expected_lines(lib_num)

    deps = [str(x) for x in (dep_job_ids or []) if x]
    dep_line = f"#SBATCH --dependency=afterok:{':'.join(deps)}" if deps else ""

    cmd = (
        f'python3 "{helper}" '
        f'--lib {lib_num} '
        f'--demux_prefix "{demux_prefix}" '
        f'--expected "{expected_lines}" '
        f'--output_root "{output_root}" '
        f'--ploidy_calls_root "{ploidy_root}" '
        f'--quant_root "{CONDITION_INDEX_ROOT}" '
        f'--contam_condition "{contam_condition or ""}" '
        f'--external_ploidy_min_prob {external_min_prob}'
    )
    if contam_condition:
        cmd += (
            f' --contam_rate "{get_contam_prefix(lib_num, contam_condition)}'
            '.contam_rate"')
    if args and args.require_refine_external_ploidy:
        cmd += " --require_external_ploidy"
    if args and args.require_refine_contam_rate:
        cmd += " --require_contam_rate"
    if args and args.refine_scoring_only:
        cmd += " --scoring_only"

    refined = get_refined_assignments_path(lib_num)
    refined_simple = get_refined_simple_assignments_path(lib_num)
    summary = get_refine_summary_path(lib_num)

    force_clean = ""
    if force:
        out_prefix = os.path.join(output_root, f"lib{lib_num}", f"lib{lib_num}")
        force_clean = f"""
# --force means remove stale tetra_refine outputs only, not demux inputs.
rm -f "{out_prefix}.refined_assignments" \
      "{out_prefix}.refined_changed" \
      "{out_prefix}.assignments_refined" \
      "{out_prefix}.refine_summary"
"""

    script = f"""#!/bin/bash
#SBATCH --job-name={job_name}
#SBATCH --output={log_dir}/{job_name}_%j.out
#SBATCH --error={log_dir}/{job_name}_%j.err
#SBATCH --time={SLURM_TIME}
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --partition={SLURM_PARTITION}
#SBATCH --nodes=1
{dep_line}

set -euo pipefail
{module_block()}
export PATH="{SOFTWARE_BIN}:$PATH"

mkdir -p "{log_dir}" "{os.path.join(output_root, f'lib{lib_num}')}"

echo "================================================================"
echo "STAGE 6 TETRA_REFINE: optional biological assignment refinement"
echo "  Library: lib{lib_num}"
echo "  demux_prefix: {demux_prefix}"
echo "  expected_lines: {expected_lines}"
echo "  output_root: {output_root}"
echo "  ploidy_calls_root: {ploidy_root}"
echo "  contam_condition: {contam_condition}"
echo "  external_ploidy_min_prob: {external_min_prob}"
echo "  helper: {helper}"
echo "  Started: $(date)"
echo "  Node: $(hostname)"
echo "================================================================"

if [[ ! -f "{helper}" ]]; then
    echo "❌ ERROR: Required TETRA_REFINE helper script missing: {helper}"
    exit 1
fi
if [[ ! -x "{os.path.join(SOFTWARE_BIN, 'tetra_refine')}" ]]; then
    echo "❌ ERROR: tetra_refine binary missing or not executable: {os.path.join(SOFTWARE_BIN, 'tetra_refine')}"
    exit 1
fi
{force_clean}

echo "Running: {cmd}"
{cmd}

if [[ ! -f "{refined}" ]]; then
    echo "❌ ERROR: Expected tetra_refine output missing: {refined}"
    exit 1
fi
if [[ ! -f "{summary}" ]]; then
    echo "❌ ERROR: Expected tetra_refine summary missing: {summary}"
    exit 1
fi
if [[ ! -f "{refined_simple}" ]]; then
    echo "❌ ERROR: Expected simple refined assignments missing: {refined_simple}"
    exit 1
fi

echo "✅ TETRA_REFINE output: {refined}"
echo "Finished: $(date)"
"""
    os.makedirs(log_dir, exist_ok=True)
    os.makedirs(script_dir, exist_ok=True)
    script_path = os.path.join(script_dir, f"{job_name}.sbatch")
    with open(script_path, "w") as f:
        f.write(script)
    os.chmod(script_path, stat.S_IRWXU | stat.S_IRGRP | stat.S_IXGRP |
             stat.S_IROTH | stat.S_IXOTH)
    return script_path


# =============================================================================
# Stage 7 POSTHOC: tetra_score_calls + sample-swap audit
# =============================================================================

def generate_posthoc_script(lib_num, dep_job_id=None, force=False):
    """Generate sbatch script for per-cell scoring and swap audit for one library."""
    job_name = f"posthoc_lib{lib_num}"
    log_dir = os.path.join(get_log_dir(), "POSTHOC")
    script_dir = get_script_dir()
    audit_root = get_audit_root()
    audit_lib_dir = get_audit_lib_dir(lib_num)
    prefix = get_demux_prefix(lib_num)
    samples = get_default_vcf_samples_for_audit(lib_num)
    expected_meta = get_expected_pool_metadata()
    allowed = get_allowed_identities()

    dep_line = f"#SBATCH --dependency=afterok:{dep_job_id}" if dep_job_id is not None else ""

    val_files = [
        prefix + ".counts", prefix + ".condf", prefix + ".samples",
        prefix + ".assignments", prefix + ".diagnostics.gz", prefix + ".runner_ups.gz",
        prefix + ".species_counts", prefix + ".species_condf", prefix + ".species_samples",
        samples, expected_meta, PANEL_METADATA,
    ]
    if allowed:
        val_files.append(allowed)
    val_lines = []
    for vf in val_files:
        val_lines.append(
            f'if [[ ! -f "{vf}" ]]; then\n'
            f'    echo "❌ ERROR: Required POSTHOC file missing: {vf}"\n'
            f'    exit 1\n'
            f'fi')
    validation_block = "\n".join(val_lines)

    skip_block = ""
    if not force:
        skip_block = (
            f'if [[ -f "{get_call_qc_path(lib_num)}" ]] && '
            f'[[ -f "{get_species_qc_path(lib_num)}" ]] && '
            f'[[ -f "{get_swap_report_path(lib_num)}" ]]; then\n'
            f'    echo "✅ POSTHOC outputs already exist for lib{lib_num}"\n'
            f'    echo "Skipping (use --force to rerun)"\n'
            f'    exit 0\n'
            f'fi\n'
        )

    allowed_args = ""
    if allowed:
        allowed_args = f' --allowed_identities "{allowed}" --allowed_identities_source file'

    swap_audit_prepare = resolve_process_script("swap_audit_prepare.py")
    swap_audit_summarize = resolve_process_script("swap_audit_summarize.py")
    swap_audit_aggregate = resolve_process_script("swap_audit_aggregate.py")

    prep_cmd = (
        f'python3 "{swap_audit_prepare}" '
        f'--vcf_samples "{samples}" '
        f'--expected_pool_metadata "{expected_meta}" '
        f'--demux_root "{get_demux_dir(lib_num)}" '
        f'--audit_root "{audit_root}" '
        f'--refined_assignments_root "{TETRA_REFINE_ROOT}" '
        f'--libraries lib{lib_num} '
        f'--panel_metadata "{PANEL_METADATA}" '
        f'--overwrite{allowed_args}'
    )
    summarize_cmd = (
        f'python3 "{swap_audit_summarize}" '
        f'--lib lib{lib_num} '
        f'--capabilities_manifest "{os.path.join(audit_lib_dir, f"lib{lib_num}.capabilities.json")}" '
        f'--audit_root "{audit_root}" '
        f'--expected_metadata "{expected_meta}" '
        f'--panel_metadata "{PANEL_METADATA}" '
        f'--out_prefix "{os.path.join(audit_lib_dir, f"lib{lib_num}")}" '
        f'--overwrite'
    )
    capabilities_manifest = os.path.join(audit_lib_dir, f"lib{lib_num}.capabilities.json")
    per_lib_score_commands = os.path.join(audit_lib_dir, f"lib{lib_num}.run_tetra_score_calls_commands.sh")
    per_lib_audit_commands = os.path.join(audit_lib_dir, f"lib{lib_num}.run_audit_demux_parallel_commands.sh")

    script = f"""#!/bin/bash
#SBATCH --job-name={job_name}
#SBATCH --output={log_dir}/{job_name}_%j.out
#SBATCH --error={log_dir}/{job_name}_%j.err
#SBATCH --time={SLURM_TIME}
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --partition={SLURM_PARTITION}
#SBATCH --nodes=1
{dep_line}

set -euo pipefail
{module_block()}
export PATH="{SOFTWARE_BIN}:$PATH"
export OMP_NUM_THREADS="${{SLURM_CPUS_PER_TASK}}"

mkdir -p "{log_dir}" "{audit_root}"

echo "================================================================"
echo "STAGE 7 POSTHOC: tetra_score_calls + sample-swap audit"
echo "  Library: lib{lib_num}"
echo "  Demux prefix: {prefix}"
echo "  Audit root: {audit_root}"
echo "  Started: $(date)"
echo "  Node: $(hostname)"
echo "================================================================"

# Validate helper scripts resolved by orchestrator at job-generation time
for helper in "{swap_audit_prepare}" "{swap_audit_summarize}" "{swap_audit_aggregate}"; do
    if [[ ! -f "$helper" ]]; then
        echo "❌ ERROR: Required helper script missing: $helper"
        exit 1
    fi
done

{validation_block}

{skip_block}

echo "=== Preparing audit manifest and commands ==="
echo "Running: {prep_cmd}"
{prep_cmd}

manifest="{capabilities_manifest}"
score_commands="{per_lib_score_commands}"
audit_commands="{per_lib_audit_commands}"

if [[ ! -s "$manifest" ]]; then
    echo "❌ ERROR: swap_audit_prepare did not create capabilities manifest: $manifest"
    echo "POSTHOC preflight failed; refusing to run tetra_score_calls, audit demux, summarize, or aggregate."
    exit 1
fi

if [[ ! -s "$score_commands" ]]; then
    echo "❌ ERROR: swap_audit_prepare did not create per-library tetra_score_calls commands: $score_commands"
    exit 1
fi

if [[ ! -s "$audit_commands" ]]; then
    echo "❌ ERROR: swap_audit_prepare did not create per-library audit demux commands: $audit_commands"
    exit 1
fi

echo "=== Running tetra_score_calls ==="
bash "$score_commands"

echo "=== Running demux reuse-count swap audit ==="
bash "$audit_commands"

echo "=== Summarizing swap audit ==="
echo "Running: {summarize_cmd}"
{summarize_cmd}

MISSING=0
for f in \
    "{get_call_qc_path(lib_num)}" \
    "{get_species_qc_path(lib_num)}" \
    "{get_swap_report_path(lib_num)}"; do
    if [[ ! -f "$f" ]]; then
        echo "❌ ERROR: Expected POSTHOC output missing: $f"
        MISSING=$((MISSING + 1))
    else
        echo "✅ $f"
    fi
done
if [[ "$MISSING" -gt 0 ]]; then
    exit 1
fi

echo "Finished POSTHOC lib{lib_num}: $(date)"
"""
    os.makedirs(log_dir, exist_ok=True)
    os.makedirs(script_dir, exist_ok=True)
    script_path = os.path.join(script_dir, f"{job_name}.sbatch")
    with open(script_path, "w") as f:
        f.write(script)
    os.chmod(script_path, stat.S_IRWXU | stat.S_IRGRP | stat.S_IXGRP |
             stat.S_IROTH | stat.S_IXOTH)
    return script_path


# =============================================================================
# POSTHOC_SUMMARY: re-summarize existing audit outputs only
# =============================================================================

def generate_posthoc_summary_script(lib_num, dep_job_id=None):
    """Generate a small SLURM job that reruns only swap_audit_summarize.py."""
    job_name = f"posthoc_summary_lib{lib_num}"
    log_dir = os.path.join(get_log_dir(), "POSTHOC_SUMMARY")
    script_dir = get_script_dir()
    audit_root = get_audit_root()
    audit_lib_dir = get_audit_lib_dir(lib_num)
    expected_meta = get_expected_pool_metadata()
    manifest = os.path.join(audit_lib_dir, f"lib{lib_num}.capabilities.json")
    out_prefix = os.path.join(audit_lib_dir, f"lib{lib_num}")
    swap_audit_summarize = resolve_process_script("swap_audit_summarize.py")

    dep_line = f"#SBATCH --dependency=afterok:{dep_job_id}" if dep_job_id is not None else ""
    summarize_cmd = (
        f'python3 "{swap_audit_summarize}" '
        f'--lib lib{lib_num} '
        f'--capabilities_manifest "{manifest}" '
        f'--audit_root "{audit_root}" '
        f'--expected_metadata "{expected_meta}" '
        f'--panel_metadata "{PANEL_METADATA}" '
        f'--out_prefix "{out_prefix}" '
        f'--overwrite'
    )
    swap_report = get_swap_report_path(lib_num)

    script = f"""#!/bin/bash
#SBATCH --job-name={job_name}
#SBATCH --output={log_dir}/{job_name}_%j.out
#SBATCH --error={log_dir}/{job_name}_%j.err
#SBATCH --time={SLURM_TIME}
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --partition={SLURM_PARTITION}
#SBATCH --nodes=1
{dep_line}

set -euo pipefail
{module_block()}
export PATH="{SOFTWARE_BIN}:$PATH"
mkdir -p "{log_dir}" "{audit_lib_dir}"

echo "================================================================"
echo "POSTHOC_SUMMARY: existing audit outputs -> refreshed verdicts"
echo "  Library: lib{lib_num}"
echo "  Audit root: {audit_root}"
echo "  Started: $(date)"
echo "  Node: $(hostname)"
echo "================================================================"

if [[ ! -f "{swap_audit_summarize}" ]]; then
    echo "ERROR: Required summary helper missing: {swap_audit_summarize}"
    exit 1
fi
if [[ ! -f "{manifest}" ]]; then
    echo "ERROR: Capabilities manifest missing: {manifest}"
    exit 1
fi

echo "Running summary-only refresh. No demux_parallel or tetra_score_calls will run."
echo "Command: {summarize_cmd}"
{summarize_cmd} &
summary_pid=$!
while kill -0 "$summary_pid" 2>/dev/null; do
    echo "  lib{lib_num} summary still running: $(date)"
    sleep 60
done
wait "$summary_pid"

if [[ ! -f "{swap_report}" ]]; then
    echo "ERROR: Expected refreshed swap report missing: {swap_report}"
    exit 1
fi

echo "POSTHOC_SUMMARY complete for lib{lib_num}: $(date)"
echo "  {swap_report}"
"""

    os.makedirs(log_dir, exist_ok=True)
    os.makedirs(script_dir, exist_ok=True)
    script_path = os.path.join(script_dir, f"{job_name}.sbatch")
    with open(script_path, "w") as f:
        f.write(script)
    os.chmod(script_path, stat.S_IRWXU | stat.S_IRGRP | stat.S_IXGRP |
             stat.S_IROTH | stat.S_IXOTH)
    return script_path


# =============================================================================
# UNEXPECTED_COMPONENT_NN: genotype-alert cells x existing NN ploidy calls
# =============================================================================

def generate_unexpected_component_nn_script(lib_num, dep_job_id=None):
    """Generate one analysis-only SLURM job for a library NN cross-check."""
    job_name = f"unexpected_component_nn_lib{lib_num}"
    log_dir = os.path.join(get_log_dir(), "UNEXPECTED_COMPONENT_NN")
    script_dir = get_script_dir()
    audit_root = get_audit_root()
    out_root = get_unexpected_component_nn_root()
    libname = f"lib{lib_num}"
    helper = resolve_process_script("analyze_unexpected_component_nn.py")
    dep_line = f"#SBATCH --dependency=afterok:{dep_job_id}" if dep_job_id is not None else ""

    swap_report = get_swap_report_path(lib_num)
    swap_scores = os.path.join(get_audit_lib_dir(lib_num), f"{libname}.swap_scores.tsv")
    call_qc = get_call_qc_path(lib_num)
    nn_calls = get_ploidy_calls_path(lib_num)
    summary_out = os.path.join(out_root, libname, f"{libname}.unexpected_component_nn_summary.tsv")

    cmd = (
        f'python3 "{helper}" '
        f'--lib {libname} '
        f'--audit-root "{audit_root}" '
        f'--nn-root "{get_ploidy_calls_root()}" '
        f'--out-root "{out_root}"'
    )

    script = f"""#!/bin/bash
#SBATCH --job-name={job_name}
#SBATCH --output={log_dir}/{job_name}_%j.out
#SBATCH --error={log_dir}/{job_name}_%j.err
#SBATCH --time={SLURM_TIME}
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --partition={SLURM_PARTITION}
#SBATCH --nodes=1
{dep_line}

set -euo pipefail
{module_block()}
export PATH="{SOFTWARE_BIN}:$PATH"
mkdir -p "{log_dir}" "{out_root}/{libname}"

echo "================================================================"
echo "UNEXPECTED_COMPONENT_NN: POSTHOC genotype signal x NN P(tet)"
echo "  Library: {libname}"
echo "  Started: $(date)"
echo "  Node: $(hostname)"
echo "================================================================"

for f in "{helper}" "{swap_report}" "{swap_scores}" "{call_qc}" "{nn_calls}"; do
    if [[ ! -f "$f" ]]; then
        echo "ERROR: Required input missing: $f"
        exit 1
    fi
done

echo "Analysis only: no demux, NN inference, tetra_refine, or POSTHOC recomputation."
echo "Running: {cmd}"
{cmd}

if [[ ! -f "{summary_out}" ]]; then
    echo "ERROR: Expected summary missing: {summary_out}"
    exit 1
fi

echo "UNEXPECTED_COMPONENT_NN complete for {libname}: $(date)"
"""

    os.makedirs(log_dir, exist_ok=True)
    os.makedirs(script_dir, exist_ok=True)
    script_path = os.path.join(script_dir, f"{job_name}.sbatch")
    with open(script_path, "w") as f:
        f.write(script)
    os.chmod(script_path, stat.S_IRWXU | stat.S_IRGRP | stat.S_IXGRP |
             stat.S_IROTH | stat.S_IXOTH)
    return script_path


def generate_unexpected_component_nn_aggregate_script(dep_job_ids=None, libraries=None):
    """Aggregate per-library NN cross-check outputs after all jobs finish."""
    job_name = "unexpected_component_nn_aggregate"
    log_dir = os.path.join(get_log_dir(), "UNEXPECTED_COMPONENT_NN")
    script_dir = get_script_dir()
    out_root = get_unexpected_component_nn_root()
    helper = resolve_process_script("analyze_unexpected_component_nn.py")
    deps = [str(x) for x in (dep_job_ids or []) if x]
    dep_line = f"#SBATCH --dependency=afterok:{':'.join(deps)}" if deps else ""
    lib_spec = " ".join(str(x) for x in (libraries or [])) or "1-40"
    cmd = (
        f'python3 "{helper}" --aggregate '
        f'--libraries "{lib_spec}" '
        f'--audit-root "{get_audit_root()}" '
        f'--nn-root "{get_ploidy_calls_root()}" '
        f'--out-root "{out_root}"'
    )
    aggregate_out = os.path.join(out_root, "all_libraries.unexpected_component_nn_summary.tsv")

    script = f"""#!/bin/bash
#SBATCH --job-name={job_name}
#SBATCH --output={log_dir}/{job_name}_%j.out
#SBATCH --error={log_dir}/{job_name}_%j.err
#SBATCH --time={SLURM_TIME}
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --partition={SLURM_PARTITION}
#SBATCH --nodes=1
{dep_line}

set -euo pipefail
{module_block()}
export PATH="{SOFTWARE_BIN}:$PATH"
mkdir -p "{log_dir}" "{out_root}"

if [[ ! -f "{helper}" ]]; then
    echo "ERROR: Required helper missing: {helper}"
    exit 1
fi

echo "=== Aggregate unexpected-component x NN cross-check ==="
echo "Running: {cmd}"
{cmd}

if [[ ! -f "{aggregate_out}" ]]; then
    echo "ERROR: Expected aggregate missing: {aggregate_out}"
    exit 1
fi

echo "Aggregate complete: {aggregate_out}"
"""

    os.makedirs(log_dir, exist_ok=True)
    os.makedirs(script_dir, exist_ok=True)
    script_path = os.path.join(script_dir, f"{job_name}.sbatch")
    with open(script_path, "w") as f:
        f.write(script)
    os.chmod(script_path, stat.S_IRWXU | stat.S_IRGRP | stat.S_IXGRP |
             stat.S_IROTH | stat.S_IXOTH)
    return script_path


# =============================================================================
# POSTHOC final aggregation: one shared job after per-library POSTHOC jobs
# =============================================================================

def generate_posthoc_aggregate_script(dep_job_ids=None):
    """Generate one final all-library swap-audit aggregation job."""
    job_name = "posthoc_aggregate"
    log_dir = os.path.join(get_log_dir(), "POSTHOC")
    script_dir = get_script_dir()
    audit_root = get_audit_root()
    swap_audit_aggregate = resolve_process_script("swap_audit_aggregate.py")

    deps = [str(x) for x in (dep_job_ids or []) if x]
    dep_line = f"#SBATCH --dependency=afterok:{':'.join(deps)}" if deps else ""

    aggregate_cmd = (
        f'python3 "{swap_audit_aggregate}" '
        f'--audit_root "{audit_root}" '
        f'--refined_assignments_root "{TETRA_REFINE_ROOT}" '
        f'--overwrite'
    )
    aggregate_output = os.path.join(audit_root, "all_libraries.swap_audit.tsv")

    script = f"""#!/bin/bash
#SBATCH --job-name={job_name}
#SBATCH --output={log_dir}/{job_name}_%j.out
#SBATCH --error={log_dir}/{job_name}_%j.err
#SBATCH --time={SLURM_TIME}
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --partition={SLURM_PARTITION}
#SBATCH --nodes=1
{dep_line}

set -euo pipefail
{module_block()}
export PATH="{SOFTWARE_BIN}:$PATH"
mkdir -p "{log_dir}" "{audit_root}"

if [[ ! -f "{swap_audit_aggregate}" ]]; then
    echo "ERROR: Required POSTHOC aggregation helper missing: {swap_audit_aggregate}"
    exit 1
fi

echo "=== Final all-library POSTHOC aggregation ==="
echo "Running: {aggregate_cmd}"
{aggregate_cmd}

if [[ ! -f "{aggregate_output}" ]]; then
    echo "ERROR: Expected POSTHOC aggregate missing: {aggregate_output}"
    exit 1
fi

echo "POSTHOC aggregate complete: {aggregate_output}"
"""

    os.makedirs(log_dir, exist_ok=True)
    os.makedirs(script_dir, exist_ok=True)
    script_path = os.path.join(script_dir, f"{job_name}.sbatch")
    with open(script_path, "w") as f:
        f.write(script)
    os.chmod(script_path, stat.S_IRWXU | stat.S_IRGRP | stat.S_IXGRP |
             stat.S_IROTH | stat.S_IXOTH)
    return script_path


# =============================================================================
# Stage 6 HYBRID: post-hoc independent individual/species reconciliation
# =============================================================================

def generate_hybrid_script(lib_num, dep_job_ids=None, force=False,
                           individual_condition="IND_RF",
                           species_condition="SP_RF",
                           fixed_species_condition="SP_FIXED_EMPTY"):
    """Generate sbatch script for explicit post-hoc hybrid reconciliation.

    This is not the legacy species_regularize path and not a C++ joint likelihood.
    It combines outputs from independently-run native individual and species paths.
    """
    job_name = f"hybrid_lib{lib_num}"
    log_dir = os.path.join(get_log_dir(), "HYBRID")
    script_dir = get_script_dir()
    out_dir = get_hybrid_lib_dir(lib_num)
    out_file = get_hybrid_summary_path(lib_num)

    dep_line = ""
    if dep_job_ids:
        dep_ids = ":".join(str(j) for j in dep_job_ids if j is not None)
        if dep_ids:
            dep_line = f"#SBATCH --dependency=afterok:{dep_ids}"

    val_files = [
        get_swap_report_path(lib_num),
        get_call_qc_path(lib_num),
        get_species_qc_path(lib_num),
        get_contam_prefix(lib_num, individual_condition) + ".contam_rate",
        get_contam_prefix(lib_num, species_condition) + ".contam_rate",
        get_contam_prefix(lib_num, species_condition) + ".species_prof",
    ]
    if fixed_species_condition:
        fixed_prefix = get_contam_prefix(lib_num, fixed_species_condition)
        val_files.extend([fixed_prefix + ".contam_rate", fixed_prefix + ".species_prof"])
    val_lines = []
    for vf in val_files:
        val_lines.append(
            f'if [[ ! -f "{vf}" ]]; then\n'
            f'    echo "❌ ERROR: Required HYBRID input missing: {vf}"\n'
            f'    exit 1\n'
            f'fi')
    validation_block = "\n".join(val_lines)

    skip_block = ""
    if not force:
        skip_block = (
            f'if [[ -f "{out_file}" ]]; then\n'
            f'    echo "✅ HYBRID output already exists for lib{lib_num}: {out_file}"\n'
            f'    echo "Skipping (use --force to rerun)"\n'
            f'    exit 0\n'
            f'fi\n'
        )

    hybrid_reconcile = resolve_process_script("hybrid_posthoc_reconcile.py")

    # Repair explicit .species_prof aliases if native species jobs completed under
    # a binary that wrote only .contam_prof.
    repair_aliases = []
    for cond_name in [species_condition, fixed_species_condition]:
        if cond_name:
            repair_aliases.append(species_profile_alias_block(get_contam_prefix(lib_num, cond_name)))
    repair_species_aliases_block = "\n".join(repair_aliases)

    cmd = (
        f'python3 "{hybrid_reconcile}" '
        f'--lib lib{lib_num} '
        f'--audit_root "{get_audit_root()}" '
        f'--quant_root "{CONDITION_INDEX_ROOT}" '
        f'--individual_condition "{individual_condition}" '
        f'--species_condition "{species_condition}" '
        f'--fixed_species_condition "{fixed_species_condition}" '
        f'--out "{out_file}"'
    )

    script = f"""#!/bin/bash
#SBATCH --job-name={job_name}
#SBATCH --output={log_dir}/{job_name}_%j.out
#SBATCH --error={log_dir}/{job_name}_%j.err
#SBATCH --time={SLURM_TIME}
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --partition={SLURM_PARTITION}
#SBATCH --nodes=1
{dep_line}

set -euo pipefail
{module_block()}
export PATH="{SOFTWARE_BIN}:$PATH"
mkdir -p "{log_dir}" "{out_dir}"

echo "================================================================"
echo "STAGE 8 HYBRID: post-hoc reconciliation"
echo "  Library: lib{lib_num}"
echo "  Individual condition: {individual_condition}"
echo "  Species condition: {species_condition}"
echo "  Fixed species comparator: {fixed_species_condition}"
echo "  Output: {out_file}"
echo "  Started: $(date)"
echo "  Node: $(hostname)"
echo "================================================================"

# Validate helper script resolved by orchestrator at job-generation time
if [[ ! -f "{hybrid_reconcile}" ]]; then
    echo "❌ ERROR: Required helper script missing: {hybrid_reconcile}"
    exit 1
fi

# Repair native species profile aliases before validation
{repair_species_aliases_block}

{validation_block}

{skip_block}

echo "Running: {cmd}"
{cmd}

if [[ ! -f "{out_file}" ]]; then
    echo "❌ ERROR: Expected HYBRID output missing: {out_file}"
    exit 1
fi

echo "✅ HYBRID output: {out_file}"
echo "Finished HYBRID lib{lib_num}: $(date)"
"""
    os.makedirs(log_dir, exist_ok=True)
    os.makedirs(script_dir, exist_ok=True)
    script_path = os.path.join(script_dir, f"{job_name}.sbatch")
    with open(script_path, "w") as f:
        f.write(script)
    os.chmod(script_path, stat.S_IRWXU | stat.S_IRGRP | stat.S_IXGRP |
             stat.S_IROTH | stat.S_IXOTH)
    return script_path


# =============================================================================
# Stage 9 MT_FUSION: per-cell mitochondrial ratios
# =============================================================================

def generate_mt_fusion_script(lib_num, dep_job_ids=None, force=False, args=None,
                              prefer_future_refined=False):
    """Generate one per-library mt_fusion_ratio SLURM job."""
    job_name = f"mt_fusion_lib{lib_num}"
    log_dir = os.path.join(get_log_dir(), "MT_FUSION")
    mt_variant = get_mt_rna_ambient_variant(args)
    if mt_variant:
        log_dir = os.path.join(log_dir, mt_variant)
    script_dir = get_script_dir()
    out_dir = get_mt_lib_dir(lib_num, args)
    out_prefix = get_mt_prefix(lib_num, args)
    binary = os.path.join(SOFTWARE_BIN, "mt_fusion_ratio")
    bam = get_mt_bam_path(lib_num, args)
    assignments = get_mt_assignments_path(
        lib_num, args, prefer_future_refined=prefer_future_refined)
    library_id = get_mt_library_id(lib_num, args)
    rna_ambient_path = get_mt_rna_ambient_path(lib_num, args)
    rna_ambient_provenance_path = (
        get_mt_rna_ambient_provenance_path(lib_num, args)
        if rna_ambient_path else None)
    rna_ambient_arm_contract = None
    rna_ambient_candidate_set = ""
    rna_ambient_plan_fingerprint = ""
    rna_ambient_arm = ""
    rna_ambient_assignment_basis = ""
    rna_ambient_roster_basis = ""
    if (rna_ambient_path and
            args.mt_rna_ambient_assignment_source in IDENTITY_AMBIENT_ARMS):
        arm_inputs = identity_ambient_arm_inputs(
            lib_num, args.mt_rna_ambient_assignment_source)
        if arm_inputs is None:
            raise ValueError(
                f"RNA ambient arm {args.mt_rna_ambient_assignment_source} "
                f"is not applicable to lib{lib_num}")
        rna_ambient_candidate_set = arm_inputs["candidate_set"]
        rna_ambient_plan_fingerprint = arm_inputs["plan_fingerprint"]
        rna_ambient_arm = arm_inputs["arm"]
        rna_ambient_assignment_basis = arm_inputs["assignment_basis"]
        rna_ambient_roster_basis = arm_inputs["roster_basis"]
        rna_ambient_arm_contract = get_contam_prefix(
            lib_num, args.mt_rna_ambient_condition,
            args.mt_rna_ambient_assignment_source,
        ) + ".identity_ambient_arm.tsv"

    deps = [str(x) for x in (dep_job_ids or []) if x]
    dep_line = f"#SBATCH --dependency=afterok:{':'.join(deps)}" if deps else ""

    self_learn_rho = (
        args.mt_rho_mode == "shrink" and
        args.mt_rho_reference is None and
        args.mt_pooled_rho is None
    )
    scratch_prefix = os.path.join(out_dir, f"lib{lib_num}.rho_learning_free")
    scratch_ratio = scratch_prefix + ".mt_ratio.tsv"
    learned_rho_reference = os.path.join(out_dir, f"lib{lib_num}.rho_reference.tsv")

    base_args = (
        f'"{binary}" '
        f'--bam "{bam}" '
        f'--vcf "{args.mt_vcf}" '
        f'--assignments "{assignments}" '
        f'--library_id "{library_id}" '
        f'--site_manifest "{args.mt_site_manifest}" '
        f'--assay_mode {args.mt_assay_mode} '
        f'--likelihood {args.mt_likelihood} '
        f'--single_parent_epsilon {args.mt_single_parent_epsilon:.17g} '
        f'--rho_low_information_molecules {args.mt_rho_low_information_molecules} '
        f'--site_influence_mode {args.mt_site_influence_mode} '
        f'--threads "${{SLURM_CPUS_PER_TASK}}"'
    )

    common_optional_args = ""
    if args.mt_site_calibration:
        common_optional_args += f' --site_calibration "{args.mt_site_calibration}"'
    if args.mt_site_calibration_stratum_template:
        common_optional_args += (
            ' --site_calibration_stratum '
            f'"{_expand_mt_template(args.mt_site_calibration_stratum_template, lib_num)}"'
        )
    if args.mt_mask_bed:
        common_optional_args += f' --mt_mask_bed "{args.mt_mask_bed}"'
    if args.mt_atac_include_singletons:
        common_optional_args += ' --atac_include_singletons'
    if rna_ambient_path:
        common_optional_args += (
            f' --rna_ambient_fraction_file "{rna_ambient_path}"'
            f' --rna_ambient_condition "{args.mt_rna_ambient_condition}"'
        )
        if args.mt_rna_ambient_max is not None:
            common_optional_args += (
                f' --rna_ambient_qc_max {args.mt_rna_ambient_max:.17g}')

    auto_empty_barcodes = (
        args.mt_ambient_qc_max is not None and
        not args.mt_ambient_none and
        not args.mt_ambient_profile_template and
        not args.mt_empty_barcodes_template
    )
    auto_empty_barcodes_path = os.path.join(
        out_dir, f"lib{lib_num}.mt_empty_barcodes.txt")
    filtered_barcodes = get_filtered_barcodes(lib_num)

    final_ambient_args = ""
    if args.mt_ambient_profile_template:
        final_ambient_args += (
            f' --ambient_profile "{_expand_mt_template(args.mt_ambient_profile_template, lib_num)}"'
        )
    elif args.mt_empty_barcodes_template:
        final_ambient_args += (
            f' --empty_barcodes "{_expand_mt_template(args.mt_empty_barcodes_template, lib_num)}"'
        )
    elif auto_empty_barcodes:
        final_ambient_args += f' --empty_barcodes "{auto_empty_barcodes_path}"'
    else:
        # No mt-specific ambient source is implied by the nuclear contamination
        # outputs. Default explicitly to the estimator's no-ambient mode when
        # ambient QC is not requested.
        final_ambient_args += ' --ambient_none'
    if args.mt_ambient_qc_max is not None:
        final_ambient_args += f' --ambient_qc_max {args.mt_ambient_qc_max:.17g}'

    final_base_cmd = (
        base_args +
        f' --output_prefix "{out_prefix}"' +
        ' --write_profile_grid' +
        f' --profile_grid_step {args.mt_profile_grid_step:.17g}' +
        common_optional_args +
        final_ambient_args
    )

    scratch_free_cmd = (
        base_args +
        f' --output_prefix "{scratch_prefix}"' +
        common_optional_args +
        ' --rho_mode free --ambient_none'
    )

    canonical_cmd = ""
    if not self_learn_rho:
        canonical_cmd = (
            final_base_cmd +
            f' --rho_mode {args.mt_rho_mode}' +
            f' --rho_prior_strength {args.mt_rho_prior_strength:.17g}'
        )
        if args.mt_rho_reference:
            canonical_cmd += f' --rho_reference "{args.mt_rho_reference}"'
        if args.mt_pooled_rho is not None:
            canonical_cmd += f' --pooled_rho {args.mt_pooled_rho:.17g}'

    val_files = [binary, bam, args.mt_vcf, assignments, args.mt_site_manifest]
    if args.mt_ambient_profile_template:
        val_files.append(_expand_mt_template(args.mt_ambient_profile_template, lib_num))
    if args.mt_empty_barcodes_template:
        val_files.append(_expand_mt_template(args.mt_empty_barcodes_template, lib_num))
    if auto_empty_barcodes:
        val_files.append(filtered_barcodes)
    if args.mt_site_calibration:
        val_files.append(args.mt_site_calibration)
    if args.mt_rho_reference:
        val_files.append(args.mt_rho_reference)
    if args.mt_mask_bed:
        val_files.append(args.mt_mask_bed)
    if rna_ambient_path:
        val_files.append(rna_ambient_path)
    if rna_ambient_arm_contract:
        val_files.append(rna_ambient_arm_contract)
    validation_block = "\n".join(
        f'if [[ ! -s "{path}" ]]; then echo "ERROR: Required MT_FUSION input missing or empty: {path}"; exit 1; fi'
        for path in val_files
    )

    rna_ambient_contract_function = ""
    rna_ambient_write_block = ""
    if rna_ambient_path:
        arm_contract_arg = rna_ambient_arm_contract or ""
        rna_ambient_qc_max_token = (
            "NA" if args.mt_rna_ambient_max is None else
            format(args.mt_rna_ambient_max, ".17g"))
        rna_ambient_contract_function = f'''
validate_rna_ambient_provenance() {{
    python3 - "$1" "{rna_ambient_path}" "{arm_contract_arg}" \
        "{rna_ambient_provenance_path}" "{args.mt_rna_ambient_condition}" \
        "{args.mt_rna_ambient_assignment_source}" \
        "{rna_ambient_candidate_set}" "{rna_ambient_plan_fingerprint}" \
        "lib{lib_num}" "{rna_ambient_arm}" \
        "{rna_ambient_assignment_basis}" "{rna_ambient_roster_basis}" \
        "{rna_ambient_qc_max_token}" "{args.mt_assignment_source}" \
        "{assignments}" <<'PY'
import csv
import hashlib
import json
import os
import sys

(mode, rate_path, arm_contract_path, provenance_path, condition, source,
 candidate_set, plan_fingerprint, library, arm, assignment_basis,
 roster_basis, qc_max_token, mt_assignment_source, mt_assignments) = sys.argv[1:]

def record(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    status = os.stat(path)
    return {{
        "path": os.path.abspath(path),
        "size_bytes": int(status.st_size),
        "sha256": digest.hexdigest(),
    }}

arm_record = None
if arm_contract_path:
    with open(arm_contract_path, "r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if len(rows) != 1:
        raise SystemExit("RNA ambient arm contract must contain exactly one row")
    row = rows[0]
    expected = {{
        "library": library,
        "condition": condition,
        "arm": arm,
        "arm_key": source,
        "assignment_basis": assignment_basis,
        "roster_basis": roster_basis,
        "candidate_set": candidate_set,
        "plan_fingerprint": plan_fingerprint,
        "assignment_update_mode": "iterative_frozen",
        "assignment_score_basis": "original_demux_all_arms",
    }}
    for field, value in expected.items():
        if row.get(field) != value:
            raise SystemExit(
                f"RNA ambient arm contract {{field}} mismatch: "
                f"{{row.get(field)!r}} != {{value!r}}")
    arm_record = record(arm_contract_path)

payload = {{
    "schema_version": 1,
    "library": int(library[3:]),
    "condition": condition,
    "assignment_source": source,
    "mt_assignment_source": mt_assignment_source,
    "mt_assignments": record(mt_assignments),
    "rna_ambient_qc_max": (
        None if qc_max_token == "NA" else float(qc_max_token)),
    "candidate_set": candidate_set,
    "plan_fingerprint": plan_fingerprint,
    "contam_rate": record(rate_path),
    "arm_contract": arm_record,
}}
if mode == "check":
    try:
        with open(provenance_path, "r", encoding="utf-8") as handle:
            observed = json.load(handle)
    except (OSError, json.JSONDecodeError) as exc:
        raise SystemExit(f"RNA ambient provenance is unavailable: {{exc}}")
    if observed != payload:
        raise SystemExit("RNA ambient provenance is stale")
elif mode == "write":
    tmp = f"{{provenance_path}}.tmp.{{os.getpid()}}"
    with open(tmp, "w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)
        handle.write("\\n")
    os.replace(tmp, provenance_path)
else:
    raise SystemExit(f"unknown provenance mode: {{mode}}")
PY
}}
'''
        rna_ambient_write_block = "validate_rna_ambient_provenance write"

    skip_block = ""
    if not force:
        if rna_ambient_path:
            outputs = (
                get_mt_ratio_path(lib_num, args),
                get_mt_profile_path(lib_num, args),
                get_mt_qc_path(lib_num, args),
            )
            complete = " && ".join(
                f'[[ -s "{path}" ]]' for path in outputs)
            any_artifact = " || ".join(
                f'[[ -e "{path}" ]]' for path in
                (*outputs, rna_ambient_provenance_path))
            skip_block = f'''if {complete}; then
    if validate_rna_ambient_provenance check; then
        echo "MT_FUSION outputs already match the selected RNA ambient source for lib{lib_num}; skipping"
        exit 0
    fi
    echo "ERROR: existing MT_FUSION outputs have stale RNA ambient provenance; rerun with --force" >&2
    exit 1
fi
if {any_artifact}; then
    echo "ERROR: incomplete MT_FUSION outputs/provenance exist; rerun with --force" >&2
    exit 1
fi
'''
        else:
            skip_block = (
                f'if [[ -s "{get_mt_ratio_path(lib_num, args)}" ]] && '
                f'[[ -s "{get_mt_profile_path(lib_num, args)}" ]] && '
                f'[[ -s "{get_mt_qc_path(lib_num, args)}" ]]; then\n'
                f'    echo "MT_FUSION outputs already exist for lib{lib_num}; skipping (use --force to rerun)"\n'
                f'    exit 0\n'
                f'fi\n'
            )

    force_cleanup = ""
    if force:
        force_cleanup = (
            f'rm -f "{out_prefix}.mt_ratio.tsv" "{out_prefix}.mt_profile.tsv.gz" '
            f'"{out_prefix}.mt_qc.tsv" "{out_prefix}.mt_site_counts.tsv" '
            f'"{out_prefix}.mt_ambient_profile.tsv"\n'
        )
        if rna_ambient_provenance_path:
            force_cleanup += f'rm -f "{rna_ambient_provenance_path}"\n'
        if self_learn_rho:
            force_cleanup += (
                f'rm -f "{learned_rho_reference}" "{scratch_prefix}".mt_*\n'
            )

    auto_empty_barcodes_block = ""
    if auto_empty_barcodes:
        auto_empty_barcodes_block = f"""
echo "=== AUTO MT NON-CELL BARCODE DERIVATION ==="
AUTO_EMPTY_BARCODES="{auto_empty_barcodes_path}"
MT_BAM="{bam}"
FILTERED_BARCODES="{filtered_barcodes}"
TMP_SUFFIX="${{SLURM_JOB_ID:-$$}}"
CHR_M_BARCODES_TMP="${{AUTO_EMPTY_BARCODES}}.chrM.${{TMP_SUFFIX}}.tmp"
FILTERED_BARCODES_TMP="${{AUTO_EMPTY_BARCODES}}.filtered.${{TMP_SUFFIX}}.tmp"
EMPTY_BARCODES_TMP="${{AUTO_EMPTY_BARCODES}}.new.${{TMP_SUFFIX}}.tmp"
trap 'rm -f "$CHR_M_BARCODES_TMP" "$FILTERED_BARCODES_TMP" "$EMPTY_BARCODES_TMP"' EXIT

samtools view "$MT_BAM" chrM \\
    | awk '{{for (i = 12; i <= NF; ++i) {{if ($i ~ /^CB:Z:/) {{cb = substr($i, 6); if (cb != "") print cb; break}}}}}}' \\
    | LC_ALL=C sort -u > "$CHR_M_BARCODES_TMP"

gzip -cd "$FILTERED_BARCODES" \\
    | awk '{{sub(/\\r$/, "", $1); if ($1 != "") print $1}}' \\
    | LC_ALL=C sort -u > "$FILTERED_BARCODES_TMP"

LC_ALL=C comm -23 "$CHR_M_BARCODES_TMP" "$FILTERED_BARCODES_TMP" > "$EMPTY_BARCODES_TMP"

CHR_M_BARCODE_COUNT=$(wc -l < "$CHR_M_BARCODES_TMP")
FILTERED_BARCODE_COUNT=$(wc -l < "$FILTERED_BARCODES_TMP")
EMPTY_BARCODE_COUNT=$(wc -l < "$EMPTY_BARCODES_TMP")
echo "  Unique chrM CB barcodes: $CHR_M_BARCODE_COUNT"
echo "  Filtered-cell barcodes: $FILTERED_BARCODE_COUNT"
echo "  Resulting non-cell barcodes: $EMPTY_BARCODE_COUNT"

if [[ ! -s "$EMPTY_BARCODES_TMP" ]]; then
    echo "ERROR: Auto-generated MT non-cell barcode list is empty for lib{lib_num}"
    exit 1
fi

mv -f "$EMPTY_BARCODES_TMP" "$AUTO_EMPTY_BARCODES"
rm -f "$CHR_M_BARCODES_TMP" "$FILTERED_BARCODES_TMP"
trap - EXIT
echo "  Wrote: $AUTO_EMPTY_BARCODES"
echo ""
"""

    rho_learning_block = ""
    scratch_cleanup_block = ""
    if self_learn_rho:
        rho_learning_block = f"""
echo "=== SAME-LIBRARY RHO REFERENCE LEARNING ==="
SCRATCH_PREFIX="{scratch_prefix}"
SCRATCH_RATIO="{scratch_ratio}"
RHO_REFERENCE="{learned_rho_reference}"
rm -f "$RHO_REFERENCE" "$SCRATCH_PREFIX".mt_*

echo "Running scratch free-rho fit: {scratch_free_cmd}"
{scratch_free_cmd}

if [[ ! -s "$SCRATCH_RATIO" ]]; then
    echo "ERROR: Rho-learning scratch ratio output missing or empty: $SCRATCH_RATIO"
    exit 1
fi

echo "Building same-library rho reference: $RHO_REFERENCE"
set +e
python3 - "$SCRATCH_RATIO" "$RHO_REFERENCE" "{args.mt_assay_mode}" "{library_id}" <<'PY'
import csv
import math
import statistics
import sys
from collections import defaultdict

ratio_path, out_path, assay_mode, library_id = sys.argv[1:5]

by_pair = defaultdict(list)
all_rho = []

with open(ratio_path, newline="") as handle:
    reader = csv.DictReader(handle, delimiter="\\t")
    fieldnames = set(reader.fieldnames or [])
    required = {{
        "canonical_parent1",
        "canonical_parent2",
        "converged",
    }}
    missing = required - fieldnames
    rho_column = (
        "rho" if "rho" in fieldnames else
        "overdispersion_rho" if "overdispersion_rho" in fieldnames else
        None
    )
    if missing or rho_column is None:
        if rho_column is None:
            missing = set(missing)
            missing.add("rho/overdispersion_rho")
        raise SystemExit(
            "rho-learning scratch output missing required column(s): "
            + ", ".join(sorted(missing))
        )

    for row in reader:
        converged = str(row.get("converged", "")).strip().lower()
        if converged not in {{"1", "true", "yes"}}:
            continue

        try:
            rho = float(row[rho_column])
        except (TypeError, ValueError):
            continue

        if not math.isfinite(rho):
            continue

        p1 = str(row["canonical_parent1"]).strip()
        p2 = str(row["canonical_parent2"]).strip()
        if not p1 or not p2:
            continue

        by_pair[(p1, p2)].append(rho)
        all_rho.append(rho)

if not all_rho:
    raise SystemExit(2)

with open(out_path, "w", newline="") as handle:
    writer = csv.writer(handle, delimiter="\\t", lineterminator="\\n")
    writer.writerow([
        "assay_mode",
        "library_id",
        "canonical_parent1",
        "canonical_parent2",
        "rho",
    ])

    for p1, p2 in sorted(by_pair):
        writer.writerow([
            assay_mode,
            library_id,
            p1,
            p2,
            format(statistics.median(by_pair[(p1, p2)]), ".17g"),
        ])

    writer.writerow([
        assay_mode,
        library_id,
        "*",
        "*",
        format(statistics.median(all_rho), ".17g"),
    ])
PY
RHO_BUILD_STATUS=$?
set -e

if [[ "$RHO_BUILD_STATUS" -eq 0 ]]; then
    if [[ ! -s "$RHO_REFERENCE" ]]; then
        echo "ERROR: Same-library rho reference builder returned success but output is missing/empty: $RHO_REFERENCE"
        exit 1
    fi
    echo "  Learned same-library rho reference: $RHO_REFERENCE"
    FINAL_RHO_ARGS=(--rho_mode shrink --rho_reference "$RHO_REFERENCE" --rho_prior_strength {args.mt_rho_prior_strength:.17g})
elif [[ "$RHO_BUILD_STATUS" -eq 2 ]]; then
    rm -f "$RHO_REFERENCE"
    echo "WARNING: No eligible converged finite-rho cells for lib{lib_num}; final canonical fit will use rho_mode=free."
    echo "WARNING: No cross-library rho reference will be borrowed."
    FINAL_RHO_ARGS=(--rho_mode free --rho_prior_strength {args.mt_rho_prior_strength:.17g})
else
    echo "ERROR: Same-library rho reference construction failed with exit code $RHO_BUILD_STATUS"
    exit "$RHO_BUILD_STATUS"
fi

echo "Running final canonical MT_FUSION with rho args: ${{FINAL_RHO_ARGS[*]}}"
{final_base_cmd} "${{FINAL_RHO_ARGS[@]}}"
"""
        scratch_cleanup_block = f"""
echo "Final canonical MT_FUSION outputs validated; removing noncanonical rho-learning scratch outputs."
rm -f "{scratch_prefix}".mt_*
"""
    else:
        rho_learning_block = f"""
echo "Running: {canonical_cmd}"
{canonical_cmd}
"""

    script = f"""#!/bin/bash
#SBATCH --job-name={job_name}
#SBATCH --output={log_dir}/{job_name}_%j.out
#SBATCH --error={log_dir}/{job_name}_%j.err
#SBATCH --time={SLURM_TIME}
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --partition={SLURM_PARTITION}
#SBATCH --nodes=1
{dep_line}

set -euo pipefail
{module_block()}
module load samtools/1.20

mkdir -p "{log_dir}" "{out_dir}"

{rna_ambient_contract_function}

echo "================================================================"
echo "STAGE 9 MT_FUSION: per-cell mitochondrial parental ratios"
echo "  Library: lib{lib_num}"
echo "  Manifest library_id: {library_id}"
echo "  Assignments: {assignments}"
echo "  Assay: {args.mt_assay_mode}"
echo "  RNA ambient condition: {args.mt_rna_ambient_condition or 'disabled'}"
echo "  RNA ambient rates: {rna_ambient_path or 'disabled'}"
echo "  Output prefix: {out_prefix}"
echo "  Started: $(date)"
echo "  Node: $(hostname)"
echo "================================================================"

{validation_block}
{skip_block}{force_cleanup}
{auto_empty_barcodes_block}
{rho_learning_block}
for f in "{get_mt_ratio_path(lib_num, args)}" "{get_mt_profile_path(lib_num, args)}" "{get_mt_qc_path(lib_num, args)}"; do
    if [[ ! -s "$f" ]]; then
        echo "ERROR: Expected MT_FUSION output missing or empty: $f"
        exit 1
    fi
done
{rna_ambient_write_block}
{scratch_cleanup_block}
echo "MT_FUSION complete for lib{lib_num}: $(date)"
"""
    os.makedirs(log_dir, exist_ok=True)
    os.makedirs(script_dir, exist_ok=True)
    script_path = os.path.join(script_dir, f"{job_name}.sbatch")
    with open(script_path, "w") as f:
        f.write(script)
    os.chmod(script_path, stat.S_IRWXU | stat.S_IRGRP | stat.S_IXGRP |
             stat.S_IROTH | stat.S_IXOTH)
    return script_path

def generate_mt_population_script(lib_nums, dep_job_ids=None, force=False, args=None):
    """Generate the final MT population job from reconciliation-owned metadata.

    The identity-reconciliation rich cell table is the cell-level authority.
    mt_population_structure derives the mt parent pair from the ratio output and
    uses a uniquely resolved reconciled_uid as the physical fusion replicate.
    No workbook is opened or rebuilt by this stage.
    """
    libs = sorted(set(lib_nums))
    job_name = "mt_population"
    log_dir = os.path.join(get_log_dir(), "MT_FUSION")
    mt_variant = get_mt_rna_ambient_variant(args)
    if mt_variant:
        log_dir = os.path.join(log_dir, mt_variant, "population")
    script_dir = get_script_dir()
    helper = resolve_process_script("mt_population_structure.py")
    out_prefix = get_mt_population_prefix(args)
    out_dir = os.path.dirname(out_prefix)
    deps = [str(x) for x in (dep_job_ids or []) if x]
    dep_line = f"#SBATCH --dependency=afterok:{':'.join(deps)}" if deps else ""

    ratio_paths = [get_mt_ratio_path(lib, args) for lib in libs]
    profile_paths = [get_mt_profile_path(lib, args) for lib in libs]
    reconciled_paths = [get_reconciled_cells_path(lib) for lib in libs]
    ratio_args = " ".join(f'"{x}"' for x in ratio_paths)
    profile_args = " ".join(f'"{x}"' for x in profile_paths)
    reconciled_args = " ".join(f'"{x}"' for x in reconciled_paths)
    cmd = (
        f'python3 "{helper}" '
        f'--ratio-tsv {ratio_args} '
        f'--profile-tsv {profile_args} '
        f'--reconciled-cells {reconciled_args} '
        f'--group-by "{args.mt_group_by}" '
        f'--output-prefix "{out_prefix}" '
        f'--meaningful-sd {args.mt_meaningful_sd:.17g} '
        f'--min-cells {args.mt_min_cells} '
        f'--delta-bic {args.mt_delta_bic:.17g} '
        f'--min-component-fraction {args.mt_min_component_fraction:.17g} '
        f'--min-component-cells {args.mt_min_component_cells} '
        f'--min-component-separation {args.mt_min_component_separation:.17g} '
        f'--membership-threshold {args.mt_membership_threshold:.17g} '
        f'--min-confident-fraction {args.mt_min_confident_fraction:.17g}'
    )
    if args.mt_rna_ambient_max is not None:
        cmd += f' --rna-ambient-max {args.mt_rna_ambient_max:.17g}'

    required = [helper] + ratio_paths + profile_paths + reconciled_paths
    validation_block = "\n".join(
        f'if [[ ! -s "{path}" ]]; then echo "ERROR: Required MT population input missing or empty: {path}"; exit 1; fi'
        for path in required
    )

    groups_out = out_prefix + ".mt_population_groups.tsv"
    cells_out = out_prefix + ".mt_population_cells.tsv"
    exclusions_out = out_prefix + ".mt_population_exclusions.tsv"
    params_out = out_prefix + ".mt_population_run_parameters.tsv"
    ambient_associations_out = (
        out_prefix + ".mt_population_ambient_associations.tsv")
    skip_block = ""
    if not force:
        skip_block = (
            f'if [[ -s "{groups_out}" ]] && [[ -s "{cells_out}" ]] && '
            f'[[ -s "{exclusions_out}" ]] && [[ -s "{params_out}" ]] && '
            f'[[ -s "{ambient_associations_out}" ]]; then\n'
            f'    echo "MT population outputs already exist; skipping (use --force to rerun)"\n'
            f'    exit 0\n'
            f'fi\n'
        )
    force_cleanup = ""
    if force:
        force_cleanup = f'rm -f "{out_prefix}.mt_population_"*.tsv\n'

    script = f"""#!/bin/bash
#SBATCH --job-name={job_name}
#SBATCH --output={log_dir}/{job_name}_%j.out
#SBATCH --error={log_dir}/{job_name}_%j.err
#SBATCH --time={SLURM_TIME}
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --partition={SLURM_PARTITION}
#SBATCH --nodes=1
{dep_line}

set -euo pipefail
module purge
module load miniforge/3
module load genomics-base/latest
mkdir -p "{log_dir}" "{out_dir}"

echo "================================================================"
echo "STAGE 10 MT_POPULATION: reconciliation-aware mitochondrial population structure"
echo "  Libraries: {' '.join('lib'+str(x) for x in libs)}"
echo "  Identity reconciliation root: {get_identity_reconciliation_root()}"
echo "  Group by: {args.mt_group_by}"
echo "  Output prefix: {out_prefix}"
echo "  Started: $(date)"
echo "  Node: $(hostname)"
echo "================================================================"

{validation_block}
{skip_block}{force_cleanup}
echo "Running: {cmd}"
{cmd}

for f in "{groups_out}" "{cells_out}" "{exclusions_out}" "{params_out}" "{ambient_associations_out}"; do
    if [[ ! -s "$f" ]]; then
        echo "ERROR: Expected MT population output missing or empty: $f"
        exit 1
    fi
done

echo "MT population structure complete: $(date)"
"""
    os.makedirs(log_dir, exist_ok=True)
    os.makedirs(script_dir, exist_ok=True)
    script_path = os.path.join(script_dir, f"{job_name}.sbatch")
    with open(script_path, "w") as f:
        f.write(script)
    os.chmod(script_path, stat.S_IRWXU | stat.S_IRGRP | stat.S_IXGRP |
             stat.S_IROTH | stat.S_IXOTH)
    return script_path


# =============================================================================
# IDENTITY_RECONCILIATION: post-hoc identity reconciliation and final evidence
# =============================================================================

def _identity_helper(name):
    return resolve_process_script(name)


def _candidate_axis_paths(args):
    root = os.path.abspath(
        args.identity_candidate_axis_output_root or os.path.join(
            args.identity_reconciliation_root, "candidate_axis"))
    aggregate = (
        root if os.path.isfile(os.path.join(
            root, "candidate_axis_cell_scores.tsv.gz"))
        else os.path.join(root, "aggregate"))
    return {
        "root": root,
        "pairs": os.path.join(root, "score_pairs"),
        "scorer": os.path.join(root, "scorer"),
        "aggregate": aggregate,
        "logs": os.path.join(root, "logs"),
    }


def _candidate_axis_open_tsv(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    return opener(path, "rt", newline="")


def _candidate_axis_finalized_event_keys(lib_num, args):
    decisions = os.path.join(
        get_identity_subdir(args, "decisions"),
        f"lib{lib_num}.reconciled_cells.tsv.gz",
    )
    keys = set()
    with _candidate_axis_open_tsv(decisions) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            "event_id", "event_class", "final_action",
            "proposed_donor_genotype", "original_demux_assignment",
        }
        missing = required - set(reader.fieldnames or ())
        if missing:
            raise ValueError(
                f"candidate-axis decision schema missing {sorted(missing)}: {decisions}"
            )
        for row in reader:
            event_class = str(row.get("event_class", "")).strip().upper()
            event_id = str(row.get("event_id", "")).strip()
            proposal = _canonical_identity(row.get("proposed_donor_genotype", ""))
            original = _canonical_identity(row.get("original_demux_assignment", ""))
            if (
                str(row.get("final_action", "")).strip().upper()
                == "REASSIGN_GENOTYPE"
                and event_class in {
                    "LIKELY_UNEXPECTED_INTACT_BIOLOGICAL_LINE",
                    "LIKELY_UNEXPECTED_SINGLET_POPULATION",
                }
                and event_id and proposal and original and proposal != original
                and "HOMOTET" not in event_class
                and "OCCUPANCY" not in event_class
                and "PLOIDY" not in event_class
            ):
                keys.add((f"lib{lib_num}", event_id, proposal))
    return keys


def _candidate_axis_has_finalized_event(lib_num, args):
    return bool(_candidate_axis_finalized_event_keys(lib_num, args))


def _candidate_axis_frozen_parameters(args, lib_num, has_finalized_event=None):
    if has_finalized_event is None:
        has_finalized_event = _candidate_axis_has_finalized_event(lib_num, args)
    probability = get_identity_probability_path(lib_num, args)
    provenance = get_identity_probability_provenance_path(lib_num, args)
    poor_fit = float(args.identity_poor_fit_residual)
    if not math.isfinite(poor_fit) or not 0 <= poor_fit <= 1:
        raise ValueError("frozen V3 poor_fit_residual is outside [0,1]")
    if not has_finalized_event:
        return {
            "error_ref": "0.001", "error_alt": "0.001",
            "min_evidence": str(CANDIDATE_AXIS_MIN_EVIDENCE),
            "min_evidence_source": "NOT_APPLICABLE_ZERO_EVENT",
            "poor_fit_residual": format(poor_fit, ".17g"),
            "applicability": "NOT_APPLICABLE_ZERO_EVENT",
        }
    error_ref = error_alt = 0.001
    if os.path.isfile(probability):
        errors_ref = set()
        errors_alt = set()
        with _candidate_axis_open_tsv(probability) as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if (not reader.fieldnames or
                    not {"error_ref", "error_alt"} <= set(reader.fieldnames)):
                raise ValueError(
                    f"frozen V3 table lacks error_ref/error_alt: {probability}")
            n_rows = 0
            for row in reader:
                n_rows += 1
                for field, destination in (
                        ("error_ref", errors_ref), ("error_alt", errors_alt)):
                    text = str(row.get(field, "")).strip()
                    if not text or text.upper() == "NA":
                        raise ValueError(f"frozen V3 row is missing {field}")
                    value = float(text)
                    if not math.isfinite(value):
                        raise ValueError(
                            f"frozen V3 {field} contains a nonfinite value")
                    destination.add(value)
            if n_rows == 0:
                raise ValueError(
                    f"frozen V3 probability table is empty: {probability}")
        if len(errors_ref) != 1 or len(errors_alt) != 1:
            raise ValueError(
                "frozen V3 error_ref/error_alt values are not unanimous")
        error_ref = next(iter(errors_ref))
        error_alt = next(iter(errors_alt))
        if not (0 <= error_ref <= 1 and 0 <= error_alt <= 1 and
                error_ref + error_alt < 1):
            raise ValueError(
                "frozen V3 error pair is outside the transform domain")
    if os.path.isfile(provenance):
        with _candidate_axis_open_tsv(provenance) as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        if len(rows) == 1:
            poor_fit_text = str(rows[0].get("poor_fit_residual", "")).strip()
            if poor_fit_text and poor_fit_text.upper() != "NA":
                poor_fit = float(poor_fit_text)
    if not math.isfinite(poor_fit) or not 0 <= poor_fit <= 1:
        raise ValueError("frozen V3 poor_fit_residual is outside [0,1]")
    return {
        "error_ref": format(error_ref, ".17g"),
        "error_alt": format(error_alt, ".17g"),
        "min_evidence": str(CANDIDATE_AXIS_MIN_EVIDENCE),
        "min_evidence_source": CANDIDATE_AXIS_MIN_EVIDENCE_SOURCE,
        "poor_fit_residual": format(poor_fit, ".17g"),
        "applicability": "APPLICABLE",
    }


def _identity_deps(dep_job_ids):
    deps = [str(x) for x in (dep_job_ids or []) if x]
    return f"#SBATCH --dependency=afterok:{':'.join(deps)}" if deps else ""


def _write_identity_sbatch(
        job_name, body, dep_job_ids=None, cpus=2, mem="8G",
        time=SLURM_TIME, modules=None, commands=()):
    log_dir = os.path.join(get_log_dir(), "IDENTITY_RECONCILIATION")
    script_dir = get_script_dir()
    dep_line = _identity_deps(dep_job_ids)
    if modules is None:
        modules = (
            "miniforge/3", "genomics-base/latest",
            "htslib/1.20", "cellbouncer/dev")
    module_lines = "\n".join(f"module load {module}" for module in modules)
    command_checks = "\n".join(
        f"command -v {shlex.quote(command)} >/dev/null"
        for command in commands)
    script = f"""#!/bin/bash
#SBATCH --job-name={job_name}
#SBATCH --output={log_dir}/{job_name}_%j.out
#SBATCH --error={log_dir}/{job_name}_%j.err
#SBATCH --time={SLURM_TIME}
#SBATCH --cpus-per-task={cpus}
#SBATCH --mem={mem}
#SBATCH --partition={SLURM_PARTITION}
#SBATCH --nodes=1
{dep_line}

set -euo pipefail
module purge
{module_lines}
module list 2>&1
{command_checks}
export PATH="{SOFTWARE_BIN}:$PATH"
echo "Runtime host: $(hostname)"
echo "SLURM job: ${{SLURM_JOB_ID:-not-set}}"
free -h || true
{body}
"""
    os.makedirs(log_dir, exist_ok=True)
    os.makedirs(script_dir, exist_ok=True)
    path = os.path.join(script_dir, f"{job_name}.sbatch")
    with open(path, "w") as handle:
        handle.write(script)
    os.chmod(path, stat.S_IRWXU | stat.S_IRGRP | stat.S_IXGRP | stat.S_IROTH | stat.S_IXOTH)
    return path


def generate_identity_metadata_script(lib_nums, args, dep_job_ids=None):
    root = os.path.abspath(args.identity_reconciliation_root)
    metadata = get_identity_subdir(args, "metadata")
    helper = _identity_helper("identity_reconciliation.py")
    distinguishability = os.path.join(metadata, "nuclear_panel_distinguishability.tsv")
    panel_utility = os.path.join(SOFTWARE_BIN, "nuclear_panel_distinguishability")
    body = f'''mkdir -p "{root}" "{metadata}"
echo "IDENTITY_METADATA: workbook -> deterministic genotype/UID contracts"
python3 "{helper}" metadata --xlsx "{IDENTITY_LIBRARY_CONVERSIONS_XLSX}" --panel-metadata "{PANEL_METADATA}" --output-root "{metadata}"
test -s "{os.path.join(metadata, 'library_expected_genotypes.tsv')}"
test -s "{os.path.join(metadata, 'global_biological_lines.tsv')}"
test -s "{os.path.join(metadata, 'global_donors.tsv')}"
echo "IDENTITY_METADATA: deployed NoMito nuclear-panel distinguishability context"
"{panel_utility}" --vcf "{VCF_SOURCE_PATHS['interindiv_20M']}" --output "{distinguishability}"
test -s "{distinguishability}"
'''
    return _write_identity_sbatch(
        "identity_metadata", body, dep_job_ids, cpus=1, mem="4G")


def generate_identity_candidates_script(lib_num, args, dep_job_ids=None):
    """Generate candidate hypotheses for one library only.

    Candidate construction is library-local once the shared metadata contract
    exists. Run the helper in a private per-library output directory because it
    also emits an all_libraries.identity_events_candidates.tsv summary; direct
    parallel writes to the shared candidate root would otherwise race.
    """
    out = get_identity_subdir(args, "candidates")
    job_out = os.path.join(out, ".per_library", f"lib{lib_num}")
    helper = _identity_helper("identity_reconciliation.py")
    refined_arg = f'--refined-root "{TETRA_REFINE_ROOT}"' if TETRA_REFINE_ROOT else ""
    body = f"""mkdir -p "{out}" "{job_out}"
echo "IDENTITY_CANDIDATES lib{lib_num}: sparse hypotheses only; no ambient-RNA inputs"
python3 "{helper}" candidates --libraries {lib_num} --demux-root "{get_demux_mapping_root()}" --audit-root "{args.identity_audit_root}" {refined_arg} --nn-root "{PLOIDY_CALLS_ROOT}" --metadata-root "{get_identity_subdir(args, 'metadata')}" --output-root "{job_out}" --top-k {args.identity_candidate_top_k} --max-candidates {args.identity_max_candidates}

# Publish only this library's products into the canonical shared candidate root.
# The helper's all-library event summary remains isolated for the aggregate job.
for f in "{job_out}"/lib{lib_num}.*; do
    [[ -e "$f" ]] || continue
    mv -f "$f" "{out}/"
done
test -s "{get_identity_candidate_path(lib_num, args)}"
"""
    return _write_identity_sbatch(
        f"identity_candidates_lib{lib_num}", body, dep_job_ids,
        cpus=1, mem="12G")


def generate_identity_candidates_aggregate_script(lib_nums, args, dep_job_ids=None):
    """Combine per-library candidate-event summaries after parallel generation."""
    out = get_identity_subdir(args, "candidates")
    per_root = os.path.join(out, ".per_library")
    final_path = os.path.join(out, "all_libraries.identity_events_candidates.tsv")
    libs = " ".join(str(x) for x in lib_nums)
    body = f"""mkdir -p "{out}"
echo "IDENTITY_CANDIDATES_AGGREGATE: combine per-library event summaries"
tmp="{final_path}.tmp.$$"
rm -f "$tmp"
first=1
for n in {libs}; do
    f="{per_root}/lib${{n}}/all_libraries.identity_events_candidates.tsv"
    if [[ ! -s "$f" ]]; then
        echo "ERROR: missing per-library candidate event summary: $f"
        exit 1
    fi
    if [[ "$first" -eq 1 ]]; then
        cat "$f" > "$tmp"
        first=0
    else
        tail -n +2 "$f" >> "$tmp"
    fi
done
mv -f "$tmp" "{final_path}"
test -s "{final_path}"
"""
    return _write_identity_sbatch(
        "identity_candidates_aggregate", body, dep_job_ids,
        cpus=1, mem="2G", time="1:00:00")

def generate_identity_doublet_context_script(lib_nums, args, dep_job_ids=None):
    out = get_identity_subdir(args, "doublet_context")
    helper = _identity_helper("identity_reconciliation.py")
    libs = " ".join(str(x) for x in lib_nums)
    dragon = os.path.join(SOFTWARE_BIN, "doublet_dragon")
    body = f'''mkdir -p "{out}"
echo "IDENTITY_DOUBLET_CONTEXT: Doublet Dragon on diploid-resolvable subset; real biological fusions excluded"
python3 "{helper}" doublet-context --libraries {libs} --metadata-root "{get_identity_subdir(args, 'metadata')}" --demux-root "{get_demux_mapping_root()}" --doublet-dragon "{dragon}" --output-root "{out}"
'''
    return _write_identity_sbatch(
        "identity_doublet_context", body, dep_job_ids, cpus=1, mem="4G")


def generate_identity_nuclear_script(lib_num, args, dep_job_ids=None):
    outdir = get_identity_subdir(args, "nuclear")
    prefix = get_demux_prefix(lib_num)
    cand = get_identity_candidate_path(lib_num, args)
    score = get_identity_nuclear_score_path(lib_num, args)
    folds = os.path.join(outdir, f"lib{lib_num}.identity_site_fold_scores.tsv.gz")
    cmd = (
        f'"{os.path.join(SOFTWARE_BIN, "tetra_score_calls")}" '
        f'--counts "{prefix}.counts" --samples "{prefix}.samples" --assignments "{prefix}.assignments" '
        f'--candidate_manifest "{cand}" --output "{score}" --site-fold-output "{folds}" '
        f'--pileup-sites "{prefix}.pileup_sites.tsv.gz" --pileup-observations "{prefix}.pileup_obs.tsv.gz" '
        f'--pileup-molecules "{prefix}.pileup_molecules.tsv.gz" '
        f'--site-folds {args.identity_site_folds}')
    body = f'''mkdir -p "{outdir}"
echo "NUCLEAR_IDENTITY_SCORE lib{lib_num}: reconciliation hypothesis likelihoods and site folds"
echo "Pair probabilities are deliberately deferred until reconciliation nominates an exact swap."
{cmd}
test -s "{score}"
'''
    return _write_identity_sbatch(
        f"identity_nuclear_lib{lib_num}", body, dep_job_ids,
        cpus=1, mem="32G")


def generate_identity_score_pair_manifest_script(
        lib_nums, args, dep_job_ids=None):
    helper = _identity_helper("identity_reconciliation.py")
    pairs = get_identity_subdir(args, "score_pairs")
    decisions = get_identity_subdir(args, "decisions")
    metadata = get_identity_subdir(args, "metadata")
    validation = os.path.join(
        get_identity_subdir(args, "validation"), "validation_summary.tsv")
    libraries = " ".join(str(value) for value in lib_nums)
    combined = os.path.join(
        pairs, "all_libraries.reconciliation_score_pair_summary.tsv")
    body = f'''mkdir -p "{pairs}"
echo "IDENTITY_SCORE_PAIRS: original allowed demux winner versus reconciliation-nominated swap only"
python3 "{helper}" score-pairs --libraries {libraries} --decisions-root "{decisions}" --metadata-root "{metadata}" --validation-summary "{validation}" --output-root "{pairs}"
test -s "{combined}"
'''
    return _write_identity_sbatch(
        "identity_score_pairs", body, dep_job_ids,
        cpus=1, mem="4G")


def generate_identity_probability_score_script(
        lib_num, args, dep_job_ids=None):
    outdir = get_identity_subdir(args, "nuclear")
    prefix = get_demux_prefix(lib_num)
    pair_manifest = get_identity_score_pair_path(lib_num, args)
    pair_summary = get_identity_score_pair_summary_path(lib_num, args)
    hypothesis_scores = get_identity_targeted_hypothesis_score_path(
        lib_num, args)
    fold_scores = get_identity_targeted_fold_score_path(lib_num, args)
    probability = get_identity_probability_path(lib_num, args)
    provenance = get_identity_probability_provenance_path(lib_num, args)
    scorer = os.path.join(SOFTWARE_BIN, "tetra_score_calls")
    command = (
        f'"{scorer}" '
        f'--counts "{prefix}.counts" --samples "{prefix}.samples" '
        f'--assignments "{prefix}.assignments" '
        f'--candidate_manifest "{pair_manifest}" '
        f'--output "{hypothesis_scores}" --site-fold-output "{fold_scores}" '
        f'--pileup-sites "{prefix}.pileup_sites.tsv.gz" '
        f'--pileup-observations "{prefix}.pileup_obs.tsv.gz" '
        f'--pileup-molecules "{prefix}.pileup_molecules.tsv.gz" '
        f'--probability-output "{probability}" '
        f'--probability-resamples {args.identity_probability_resamples} '
        f'--probability-seed {args.identity_probability_seed + lib_num} '
        f'--poor-fit-residual {args.identity_poor_fit_residual} '
        f'--site-folds {args.identity_site_folds}')
    body = f'''mkdir -p "{outdir}"
test -s "{pair_manifest}"
test -s "{pair_summary}"
pair_count=$(awk -F '\\t' '
  NR == 1 {{ for (i = 1; i <= NF; ++i) if ($i == "n_score_pairs") column = i; next }}
  NR == 2 {{ if (!column) exit 2; print $column; found = 1 }}
  END {{ if (!found) exit 3 }}
' "{pair_summary}")
if [[ "$pair_count" == "0" ]]; then
  echo "IDENTITY_PROBABILITY_SCORE lib{lib_num}: no reconciliation-nominated swaps; not scored"
  rm -f "{hypothesis_scores}" "{fold_scores}" "{probability}" "{provenance}"
  exit 0
fi
echo "IDENTITY_PROBABILITY_SCORE lib{lib_num}: $pair_count exact original-versus-nominated-swap pairs"
if [[ -s "{prefix}.pileup_molecules.tsv.gz" ]]; then
  echo "Scorer may orient QC to molecule evidence; v6.1 retains those rows as diagnostic and unevaluable"
else
  echo "Molecule sidecar absent; v6.1 frozen site-only QC alignment is available"
fi
{command}
test -s "{hypothesis_scores}"
test -s "{probability}"
hash_or_na() {{
  if [[ -s "$1" ]]; then sha256sum "$1" | awk '{{print $1}}'; else printf 'NA'; fi
}}
provenance_tmp="{provenance}.tmp.$$"
{{
  printf '%s\t' \
    'library' 'pair_manifest' 'pair_manifest_sha256' 'pair_summary' \
    'pair_summary_sha256' 'probability_output' 'probability_output_sha256' \
    'counts_path' 'counts_sha256' 'samples_path' 'samples_sha256' \
    'assignments_path' 'assignments_sha256' 'pileup_sites_path' \
    'pileup_sites_sha256' 'pileup_observations_path' \
    'pileup_observations_sha256' 'pileup_molecules_path' \
    'pileup_molecules_sha256' 'scorer_binary' 'scorer_binary_sha256' \
    'probability_resamples' 'probability_seed' 'poor_fit_residual' 'site_folds'
  printf '%s\n' 'schema_version'
  printf '%s\t' \
    'lib{lib_num}' '{pair_manifest}' "$(hash_or_na '{pair_manifest}')" \
    '{pair_summary}' "$(hash_or_na '{pair_summary}')" \
    '{probability}' "$(hash_or_na '{probability}')" \
    '{prefix}.counts' "$(hash_or_na '{prefix}.counts')" \
    '{prefix}.samples' "$(hash_or_na '{prefix}.samples')" \
    '{prefix}.assignments' "$(hash_or_na '{prefix}.assignments')" \
    '{prefix}.pileup_sites.tsv.gz' "$(hash_or_na '{prefix}.pileup_sites.tsv.gz')" \
    '{prefix}.pileup_obs.tsv.gz' "$(hash_or_na '{prefix}.pileup_obs.tsv.gz')" \
    '{prefix}.pileup_molecules.tsv.gz' "$(hash_or_na '{prefix}.pileup_molecules.tsv.gz')" \
    '{scorer}' "$(hash_or_na '{scorer}')" '{args.identity_probability_resamples}' \
    '{args.identity_probability_seed + lib_num}' '{args.identity_poor_fit_residual}' \
    '{args.identity_site_folds}'
  printf '%s\n' 'identity_pair_probability_provenance_v1'
}} > "$provenance_tmp"
mv -f "$provenance_tmp" "{provenance}"
test -s "{provenance}"
'''
    return _write_identity_sbatch(
        f"identity_probability_lib{lib_num}", body, dep_job_ids,
        cpus=4, mem="32G", time="1-00:00:00")


def generate_identity_probability_aggregate_script(
        lib_nums, args, dep_job_ids=None, include_ambient=True,
        reuse_frozen=False):
    helper = _identity_helper("identity_probability_aggregate.py")
    output_root = get_identity_score_output_root(args)
    libraries = " ".join(str(value) for value in lib_nums)
    ambient_arg = ""
    if include_ambient and args.identity_score_ambient_root:
        ambient_arg = (
            f'--ambient-root "{os.path.abspath(args.identity_score_ambient_root)}" ')
    if args.identity_score_ambient_condition:
        ambient_arg += (
            f'--ambient-condition "{args.identity_score_ambient_condition}" ')
    reuse_arg = "--reuse-frozen-probabilities " if reuse_frozen else ""
    scorer = os.path.join(SOFTWARE_BIN, "tetra_score_calls")
    body = f'''mkdir -p "{output_root}"
echo "IDENTITY_PROBABILITY_AGGREGATE v6.1: frozen site-aligned QC -> separated event scopes"
python3 "{helper}" --identity-root "{os.path.abspath(args.identity_reconciliation_root)}" --libraries {libraries} --output-root "{output_root}" {ambient_arg}{reuse_arg}--scorer-binary "{scorer}"
test -s "{os.path.join(output_root, 'cell_swap_identity_scores.tsv.gz')}"
test -s "{os.path.join(output_root, 'library_swap_event_summary.tsv')}"
test -s "{os.path.join(output_root, 'library_swap_event_scope_summary.tsv')}"
test -s "{os.path.join(output_root, 'cross_library_swap_summary.tsv')}"
test -s "{os.path.join(output_root, 'identity_score_provenance.tsv')}"
test -s "{os.path.join(output_root, 'plot_data', 'plotting_handoff.md')}"
'''
    return _write_identity_sbatch(
        "identity_probability_aggregate", body, dep_job_ids,
        cpus=8, mem="128G", time="4:00:00")


def _write_candidate_axis_sbatch(
        job_name, body, args, dep_job_ids=None, logical_dependencies=None,
        cpus=2, mem="8G", modules=("miniforge/3", "genomics-base/latest"),
        commands=("python3",)):
    paths = _candidate_axis_paths(args)
    os.makedirs(paths["logs"], exist_ok=True)
    dependency = _identity_deps(dep_job_ids)
    logical = ",".join(logical_dependencies or ()) or "NONE"
    module_lines = "\n".join(f"module load {name}" for name in modules)
    command_lines = "\n".join(
        f"command -v {name} >/dev/null" for name in commands)
    python_check = (
        "python3 -c 'import argparse,csv,gzip,math,pathlib'"
        if "python3" in commands else "")
    script = f'''#!/bin/bash
#SBATCH --job-name={job_name}
#SBATCH --output={paths["logs"]}/{job_name}_%j.out
#SBATCH --error={paths["logs"]}/{job_name}_%j.err
#SBATCH --time={SLURM_TIME}
#SBATCH --cpus-per-task={cpus}
#SBATCH --mem={mem}
#SBATCH --partition={SLURM_PARTITION}
#SBATCH --nodes=1
{dependency}
# Candidate-axis logical-dependency: {logical}

set -euo pipefail
export LC_ALL=C
export LANG=C
module purge
{module_lines}
module list
{command_lines}
{python_check}
echo "Runtime host: $(hostname)"
echo "SLURM job: ${{SLURM_JOB_ID:-not-set}}"
{body}
'''
    path = os.path.join(paths["logs"], f"{job_name}.sbatch")
    with open(path, "w", encoding="utf-8") as handle:
        handle.write(script)
    os.chmod(
        path,
        stat.S_IRWXU | stat.S_IRGRP | stat.S_IXGRP |
        stat.S_IROTH | stat.S_IXOTH,
    )
    syntax = subprocess.run(
        ["bash", "-n", path], capture_output=True, text=True, check=False)
    if syntax.returncode != 0:
        raise RuntimeError(
            f"candidate-axis sbatch syntax validation failed for {path}: "
            f"{syntax.stderr.strip()}")
    return path


def generate_identity_candidate_axis_pairs_script(
        lib_num, args, dep_job_ids=None):
    paths = _candidate_axis_paths(args)
    helper = _identity_helper("identity_reconciliation.py")
    prefix = get_demux_prefix(lib_num)
    parameters = _candidate_axis_frozen_parameters(args, lib_num)
    frozen_v3_lines = []
    for option, path in (
            ("frozen-v3-probability",
             get_identity_probability_path(lib_num, args)),
            ("frozen-v3-provenance",
             get_identity_probability_provenance_path(lib_num, args))):
        if os.path.isfile(path):
            frozen_v3_lines.append(
                f"  --{option} {shlex.quote(path)} \\")
    frozen_v3 = "\n".join(frozen_v3_lines)
    if frozen_v3:
        frozen_v3 += "\n"
    validation = os.path.join(
        get_identity_subdir(args, "validation"), "validation_summary.tsv")
    v6_source = _identity_helper("identity_probability_aggregate.py")
    frozen_v6 = ""
    if args.identity_candidate_axis_v6_3_root:
        v6_root = os.path.abspath(args.identity_candidate_axis_v6_3_root)
        v6_cell = os.path.join(v6_root, "cell_swap_identity_scores.tsv.gz")
        v6_review = os.path.join(v6_root, "review_only_identity_scores.tsv.gz")
        if os.path.isfile(v6_cell) and os.path.isfile(v6_review):
            frozen_v6 = (
                f' --frozen-v6-3-cell {shlex.quote(v6_cell)}'
                f' --frozen-v6-3-review-only {shlex.quote(v6_review)}'
            )
    targeted = ""
    if args.identity_candidate_axis_event_id:
        targeted = (
            f' --event-id {shlex.quote(args.identity_candidate_axis_event_id)}'
            f' --proposed-identity {shlex.quote(args.identity_candidate_axis_proposal)}'
        )
    body = f'''mkdir -p {shlex.quote(paths["pairs"])}
python3 -B {shlex.quote(helper)} candidate-axis-pairs \\
  --libraries {lib_num} \\
  --decisions-root {shlex.quote(get_identity_subdir(args, "decisions"))} \\
  --metadata-root {shlex.quote(get_identity_subdir(args, "metadata"))} \\
  --validation-summary {shlex.quote(validation)} \\
  --output-root {shlex.quote(paths["pairs"])} \\
  --samples {shlex.quote(prefix + ".samples")} \\
  --pileup-sites {shlex.quote(prefix + ".pileup_sites.tsv.gz")} \\
  --pileup-observations {shlex.quote(prefix + ".pileup_obs.tsv.gz")} \\
{frozen_v3}  --v6-3-source {shlex.quote(v6_source)} \\
  --min-evidence {parameters["min_evidence"]} \\
  --min-evidence-source {shlex.quote(parameters["min_evidence_source"])} \\
  --poor-fit-residual {parameters["poor_fit_residual"]}{targeted}{frozen_v6}
test -s {shlex.quote(os.path.join(paths["pairs"], f"lib{lib_num}.candidate_axis_pair_source_audit.tsv"))}
test -s {shlex.quote(os.path.join(paths["pairs"], f"lib{lib_num}.candidate_axis_pairs.tsv.gz"))}
test -s {shlex.quote(os.path.join(paths["pairs"], f"lib{lib_num}.candidate_axis_pair_exclusions.tsv.gz"))}
test -s {shlex.quote(os.path.join(paths["pairs"], f"lib{lib_num}.candidate_axis_pair_summary.tsv"))}
'''
    return _write_candidate_axis_sbatch(
        f"IDENTITY_CANDIDATE_AXIS_PAIRS_lib{lib_num}", body, args, dep_job_ids,
        cpus=1, mem="4G")


def generate_identity_candidate_axis_score_script(
        lib_num, args, dep_job_ids=None):
    paths = _candidate_axis_paths(args)
    aggregate_helper = _identity_helper("identity_candidate_axis_aggregate.py")
    scorer = os.path.join(SOFTWARE_BIN, "tetra_score_calls")
    prefix = get_demux_prefix(lib_num)
    parameters = _candidate_axis_frozen_parameters(args, lib_num)
    pair_manifest = os.path.join(
        paths["pairs"], f"lib{lib_num}.candidate_axis_pairs.tsv.gz")
    pair_summary = os.path.join(
        paths["pairs"], f"lib{lib_num}.candidate_axis_pair_summary.tsv")
    raw_output = os.path.join(
        paths["scorer"], f"lib{lib_num}.candidate_axis_scores.tsv.gz")
    provenance = os.path.join(
        paths["scorer"], f"lib{lib_num}.candidate_axis_score_provenance.tsv")
    body = f'''mkdir -p {shlex.quote(paths["scorer"])}
test -s {shlex.quote(pair_manifest)}
test -s {shlex.quote(pair_summary)}
temp_parent="${{SLURM_TMPDIR:-/dev/shm}}"
if [[ "$temp_parent" != /* || ! -d "$temp_parent" || ! -w "$temp_parent" ]]; then
  echo "ERROR: candidate-axis temp parent must be an absolute writable directory: $temp_parent" >&2
  exit 1
fi
temp_parent=$(realpath -e -- "$temp_parent")
job_tmp=$(mktemp -d -- "$temp_parent/identity_candidate_axis_${{SLURM_JOB_ID:-manual}}_XXXXXX")
if [[ "$(dirname -- "$job_tmp")" != "$temp_parent" || ! -d "$job_tmp" ]]; then
  echo "ERROR: unsafe candidate-axis temporary child: $job_tmp" >&2
  exit 1
fi
cleanup_candidate_axis_tmp() {{
  if [[ -n "${{job_tmp:-}}" && -d "$job_tmp" && "$(dirname -- "$job_tmp")" == "$temp_parent" && "$(basename -- "$job_tmp")" == identity_candidate_axis_* ]]; then
    rm -rf -- "$job_tmp"
  fi
}}
trap cleanup_candidate_axis_tmp EXIT
score_command=(
  {shlex.quote(scorer)}
  --samples {shlex.quote(prefix + ".samples")}
  --candidate_manifest {shlex.quote(pair_manifest)}
  --pileup-sites {shlex.quote(prefix + ".pileup_sites.tsv.gz")}
  --pileup-observations {shlex.quote(prefix + ".pileup_obs.tsv.gz")}
  --candidate-axis-output {shlex.quote(raw_output)}
  --candidate-axis-temp-dir "$job_tmp"
  --libname lib{lib_num}
  --error_ref {parameters["error_ref"]}
  --error_alt {parameters["error_alt"]}
  --min_evidence {parameters["min_evidence"]}
  --poor-fit-residual {parameters["poor_fit_residual"]}
)
printf -v score_command_text '%q ' "${{score_command[@]}}"
provenance_common=(
  --output {shlex.quote(provenance)}
  --library lib{lib_num}
  --command "$score_command_text"
  --input {shlex.quote("pair_manifest=" + pair_manifest)}
  --input {shlex.quote("pair_summary=" + pair_summary)}
  --input {shlex.quote("samples=" + prefix + ".samples")}
  --input {shlex.quote("pileup_sites=" + prefix + ".pileup_sites.tsv.gz")}
  --input {shlex.quote("pileup_observations=" + prefix + ".pileup_obs.tsv.gz")}
  --scorer-binary {shlex.quote(scorer)}
  --temp-root "$job_tmp"
  --error-ref {parameters["error_ref"]}
  --error-alt {parameters["error_alt"]}
  --min-evidence {parameters["min_evidence"]}
  --poor-fit-residual {parameters["poor_fit_residual"]}
)
"${{score_command[@]}}"
python3 -B {shlex.quote(aggregate_helper)} score-provenance \\
  "${{provenance_common[@]}}" --raw-output {shlex.quote(raw_output)}
test -s {shlex.quote(raw_output)}
test -s {shlex.quote(provenance)}
'''
    return _write_candidate_axis_sbatch(
        f"IDENTITY_CANDIDATE_AXIS_SCORE_lib{lib_num}", body, args, dep_job_ids,
        logical_dependencies=(f"IDENTITY_CANDIDATE_AXIS_PAIRS_lib{lib_num}",),
        cpus=1, mem="64G",
        modules=("miniforge/3", "genomics-base/latest", "htslib/1.20", "cellbouncer/dev"),
        commands=(
            "tetra_score_calls", "python3", "realpath", "mktemp",
            "dirname", "basename", "rm"))


def generate_identity_candidate_axis_aggregate_script(
        lib_nums, args, dep_job_ids=None):
    paths = _candidate_axis_paths(args)
    helper = _identity_helper("identity_candidate_axis_aggregate.py")
    v6_arg = ""
    if args.identity_candidate_axis_v6_3_root:
        v6_root = os.path.abspath(args.identity_candidate_axis_v6_3_root)
        if all(os.path.isfile(os.path.join(v6_root, name)) for name in (
                "cell_swap_identity_scores.tsv.gz",
                "review_only_identity_scores.tsv.gz")):
            v6_arg = " --v6-3-root " + shlex.quote(v6_root)
    targeted = ""
    if args.identity_candidate_axis_event_id:
        targeted = (
            f' --event-id {shlex.quote(args.identity_candidate_axis_event_id)}'
            f' --proposed-identity {shlex.quote(args.identity_candidate_axis_proposal)}'
        )
    libraries = " ".join(str(value) for value in lib_nums)
    body = f'''python3 -B {shlex.quote(helper)} \\
  --identity-root {shlex.quote(os.path.abspath(args.identity_reconciliation_root))} \\
  --libraries {libraries} \\
  --pair-root {shlex.quote(paths["pairs"])} \\
  --score-root {shlex.quote(paths["scorer"])} \\
  --output-root {shlex.quote(paths["aggregate"])}{targeted}{v6_arg}
for output in \\
  candidate_axis_cell_scores.tsv.gz \\
  candidate_axis_event_summary.tsv \\
  candidate_axis_event_strata_summary.tsv \\
  candidate_axis_event_scope_summary.tsv \\
  candidate_axis_applied_vs_retained_contrast.tsv \\
  candidate_axis_pair_exclusions.tsv.gz \\
  candidate_axis_provenance.tsv \\
  candidate_axis_run_summary.tsv \\
  candidate_axis_output_manifest.tsv; do
  test -s {shlex.quote(paths["aggregate"])}"/$output"
done
'''
    return _write_candidate_axis_sbatch(
        "IDENTITY_CANDIDATE_AXIS_AGGREGATE", body, args, dep_job_ids,
        logical_dependencies=tuple(
            [f"IDENTITY_CANDIDATE_AXIS_PAIRS_lib{value}" for value in lib_nums]
            + [f"IDENTITY_CANDIDATE_AXIS_SCORE_lib{value}" for value in lib_nums]
        ),
        cpus=1, mem="4G")


def generate_identity_mt_script(lib_num, args, dep_job_ids=None):
    outdir = get_identity_subdir(args, "mt")
    outprefix = os.path.join(outdir, f"lib{lib_num}")
    score = get_mt_identity_score_path(lib_num, args)
    cand = get_identity_candidate_path(lib_num, args)
    bam = get_bam_path(lib_num) or os.path.join(get_lib_dir(lib_num), "gex.bam")
    binary = os.path.join(SOFTWARE_BIN, "mt_identity_score")
    optional = _identity_helper("identity_reconciliation.py")
    required = "1" if args.identity_require_mt else "0"
    cmd = (
        f'"{binary}" --bam "{bam}" --vcf "{args.mt_panel}" --candidate_manifest "{cand}" '
        f'--library_id "lib{lib_num}" --site_manifest "{args.mt_site_manifest}" '
        f'--haplotype_groups "{args.mt_haplotype_groups}" --haplotype_pairwise "{args.mt_haplotype_pairwise}" '
        f'--output_prefix "{outprefix}" --threads "${{SLURM_CPUS_PER_TASK}}"')
    body = f'''mkdir -p "{outdir}"
required={required}
missing=0
for f in "{binary}" "{bam}" "{args.mt_panel}" "{args.mt_site_manifest}" "{args.mt_haplotype_pairwise}"; do
    [[ -s "$f" ]] || missing=1
done
if [[ "$missing" -ne 0 ]]; then
    if [[ "$required" -eq 1 ]]; then echo "ERROR: required mitochondrial identity evidence input missing"; exit 1; fi
    python3 "{optional}" optional-status --candidate-manifest "{cand}" --modality mt --status MT_UNAVAILABLE --output "{score}"
    exit 0
fi
echo "MT_IDENTITY_SCORE lib{lib_num}: one RNA chrM pass; ambient RNA is not an input"
{cmd}
test -s "{score}"
'''
    return _write_identity_sbatch(f"identity_mt_lib{lib_num}", body, dep_job_ids, cpus=4, mem="64G", time="1-00:00:00")


def _identity_expand_path(template, lib_num):
    if not template:
        return ""
    return template.format(lib=lib_num, lib_num=lib_num, library=f"lib{lib_num}")


def get_identity_atac_bam(lib_num, args):
    if args.identity_atac_bam_template:
        return _identity_expand_path(args.identity_atac_bam_template, lib_num)
    return os.path.join(args.identity_atac_root, f"Tet_2025_Multiome-ATAC_{lib_num}", "atac.bam")


def generate_identity_atac_script(lib_num, args, dep_job_ids=None):
    outdir = get_identity_subdir(args, "atac")
    cand = get_identity_candidate_path(lib_num, args)
    final_score = get_atac_identity_score_path(lib_num, args)
    optional = _identity_helper("identity_reconciliation.py")
    if args.identity_evidence_mode == "rna":
        body = f'''mkdir -p "{outdir}"
python3 "{optional}" optional-status --candidate-manifest "{cand}" --modality atac --status ATAC_NOT_REQUESTED --output "{final_score}"
'''
        return _write_identity_sbatch(f"identity_atac_not_requested_lib{lib_num}", body, dep_job_ids, cpus=1, mem="2G", time="1:00:00")

    atac_bam = get_identity_atac_bam(lib_num, args)
    atac_panel = args.atac_panel or ""
    demux_prefix = get_demux_prefix(lib_num)
    count_prefix = os.path.join(outdir, f"lib{lib_num}.atac_count")
    raw_score = os.path.join(outdir, f"lib{lib_num}.atac_identity_scores.raw.tsv.gz")
    qc = count_prefix + ".atac_qc.tsv"
    explicit_map = _identity_expand_path(args.atac_barcode_map, lib_num) if args.atac_barcode_map else ""
    auto_map = ""
    rna_w = ""
    atac_w = ""
    if not explicit_map and args.rna_barcode_whitelist and args.atac_barcode_whitelist:
        auto_map = os.path.join(outdir, f"lib{lib_num}.atac_to_rna_barcodes.tsv")
        rna_w = _identity_expand_path(args.rna_barcode_whitelist, lib_num)
        atac_w = _identity_expand_path(args.atac_barcode_whitelist, lib_num)

    main_panel = VCF_SOURCE_PATHS["interindiv_20M"]
    rna_bam = get_bam_path(lib_num) or os.path.join(get_lib_dir(lib_num), "gex.bam")
    filtered = get_filtered_barcodes(lib_num)
    demux = os.path.join(SOFTWARE_BIN, "demux_parallel")
    scorer = os.path.join(SOFTWARE_BIN, "tetra_score_calls")
    finalizer = _identity_helper("identity_reconciliation.py")
    het_arg = f'--atac_het_vcf "{args.atac_het_panel}"' if args.atac_het_panel else ""
    base_cmd = (
        f'"{demux}" -b "{rna_bam}" -o "{count_prefix}" --vcf "{main_panel}" -f '
        f'--barcodes "{filtered}" --atac_bam "{atac_bam}" --atac_vcf "{atac_panel}" '
        f'{het_arg} --reuse_counts --skip_assignment -t "${{SLURM_CPUS_PER_TASK}}"'
    )
    direct_or_mapped = ""
    final_map_arg = ""
    if explicit_map:
        final_map_arg = f'--barcode-map "{explicit_map}"'
        direct_or_mapped = f'''if [[ ! -s "{explicit_map}" ]]; then
    python3 "{optional}" optional-status --candidate-manifest "{cand}" --modality atac --status ATAC_UNAVAILABLE --output "{final_score}"
    exit 0
fi
if ! {base_cmd} --atac_barcode_map "{explicit_map}"; then
    python3 "{optional}" optional-status --candidate-manifest "{cand}" --modality atac --status ATAC_UNAVAILABLE --output "{final_score}"
    exit 0
fi'''
    elif auto_map:
        final_map_arg = f'--barcode-map "{auto_map}"'
        direct_or_mapped = f'''# Auto mode: try the expected direct RNA barcode namespace first. If that
# overlap is inadequate, deterministically construct the paired-whitelist map
# and retry in mapped mode. Either failure remains optional evidence.
if ! {base_cmd}; then
    rm -f "{count_prefix}.atac.counts" "{qc}"
    if [[ ! -s "{rna_w}" || ! -s "{atac_w}" ]] || \\
       ! python3 "{_identity_helper("identity_reconciliation.py")}" atac-barcode-map --rna-whitelist "{rna_w}" --atac-whitelist "{atac_w}" --output "{auto_map}"; then
        python3 "{optional}" optional-status --candidate-manifest "{cand}" --modality atac --status ATAC_UNAVAILABLE --output "{final_score}"
        exit 0
    fi
    if ! {base_cmd} --atac_barcode_map "{auto_map}"; then
        python3 "{optional}" optional-status --candidate-manifest "{cand}" --modality atac --status ATAC_UNAVAILABLE --output "{final_score}"
        exit 0
    fi
fi'''
    else:
        direct_or_mapped = f'''if ! {base_cmd}; then
    python3 "{optional}" optional-status --candidate-manifest "{cand}" --modality atac --status ATAC_UNAVAILABLE --output "{final_score}"
    exit 0
fi'''

    body = f'''mkdir -p "{outdir}"
if [[ -z "{atac_panel}" || ! -s "{atac_bam}" || ! -s "{main_panel}" || ! -s "{filtered}" ]]; then
    python3 "{optional}" optional-status --candidate-manifest "{cand}" --modality atac --status ATAC_UNAVAILABLE --output "{final_score}"
    exit 0
fi
ln -sfn "{demux_prefix}.counts" "{count_prefix}.counts"
ln -sfn "{demux_prefix}.samples" "{count_prefix}.samples"
{direct_or_mapped}
if ! "{scorer}" --counts "{count_prefix}.atac.counts" --samples "{demux_prefix}.samples" --assignments "{demux_prefix}.assignments" --candidate_manifest "{cand}" --score-prefix atac --output "{raw_score}"; then
    python3 "{optional}" optional-status --candidate-manifest "{cand}" --modality atac --status ATAC_INSUFFICIENT --output "{final_score}"
    exit 0
fi
python3 "{finalizer}" atac-finalize --scores "{raw_score}" --qc "{qc}" --output "{final_score}" {final_map_arg}
test -s "{final_score}"
'''
    return _write_identity_sbatch(f"identity_atac_lib{lib_num}", body, dep_job_ids, cpus=8, mem="64G", time="1-00:00:00")


def generate_identity_reconcile_script(lib_nums, args, dep_job_ids=None):
    helper = _identity_helper("identity_reconciliation.py")
    libs = " ".join(str(x) for x in lib_nums)
    decisions = get_identity_subdir(args, "decisions")
    reports = get_identity_subdir(args, "reports")
    auto = "--auto-apply" if args.identity_auto_apply else ""
    atac_arg = (
        f'--atac-root "{get_identity_subdir(args, "atac")}" '
        if args.identity_evidence_mode == "rna-atac" else ""
    )
    body = f'''mkdir -p "{decisions}" "{reports}"
echo "IDENTITY_RECONCILE: two-pass event context then one decision pass"
python3 "{helper}" reconcile --libraries {libs} --metadata-root "{get_identity_subdir(args, 'metadata')}" --candidate-root "{get_identity_subdir(args, 'candidates')}" --doublet-context-root "{get_identity_subdir(args, 'doublet_context')}" --nuclear-root "{get_identity_subdir(args, 'nuclear')}" --mt-root "{get_identity_subdir(args, 'mt')}" {atac_arg}--nn-root "{PLOIDY_CALLS_ROOT}" --refined-root "{TETRA_REFINE_ROOT}" --refined-optional --audit-root "{args.identity_audit_root}" --demux-root "{get_demux_mapping_root()}" --panel-metadata "{PANEL_METADATA}" --panel-distinguishability "{os.path.join(get_identity_subdir(args, 'metadata'), 'nuclear_panel_distinguishability.tsv')}" --output-root "{decisions}" --reports-root "{reports}" --evidence-mode {args.identity_evidence_mode} --write-reconciled-assignments {auto}
test -s "{os.path.join(decisions, 'all_libraries.reconciled_cells.tsv.gz')}"
'''
    return _write_identity_sbatch(
        "identity_reconcile", body, dep_job_ids, cpus=1, mem="16G")


def generate_identity_validate_script(lib_nums, args, dep_job_ids=None):
    helper = _identity_helper("validate_identity_reconciliation.py")
    libs = " ".join(str(x) for x in lib_nums)
    validation = get_identity_subdir(args, "validation")
    atac_root = (
        get_identity_subdir(args, "atac")
        if args.identity_evidence_mode == "rna-atac" else "")
    body = f'''mkdir -p "{validation}"
echo "IDENTITY_VALIDATE: structural invariants; unresolved scientific evidence is nonfatal"
python3 "{helper}" --libraries {libs} --demux-root "{get_demux_mapping_root()}" --metadata-root "{get_identity_subdir(args, 'metadata')}" --candidate-root "{get_identity_subdir(args, 'candidates')}" --doublet-context-root "{get_identity_subdir(args, 'doublet_context')}" --nuclear-root "{get_identity_subdir(args, 'nuclear')}" --mt-root "{get_identity_subdir(args, 'mt')}" --atac-root "{atac_root}" --decisions-root "{get_identity_subdir(args, 'decisions')}" --reports-root "{get_identity_subdir(args, 'reports')}" --evidence-mode {args.identity_evidence_mode} --output-root "{validation}"
'''
    return _write_identity_sbatch(
        "identity_validate", body, dep_job_ids, cpus=1, mem="12G")


def generate_identity_finalize_script(
        lib_nums, args, candidate_axis_root="", frozen_ambient_root="",
        four_arm_root="", dep_job_ids=None):
    helper = _identity_helper("identity_reconciliation.py")
    libraries = " ".join(str(value) for value in lib_nums)
    aggregate = os.path.join(
        os.path.abspath(args.identity_reconciliation_root), "aggregate")
    final_assignments = os.path.join(
        os.path.abspath(args.identity_reconciliation_root), "final_assignments")
    optional = ""
    for option, value in (
            ("candidate-axis-root", candidate_axis_root),
            ("frozen-ambient-root", frozen_ambient_root),
            ("four-arm-root", four_arm_root),
            ("review-input", args.identity_review_input)):
        if value:
            optional += f" --{option} {shlex.quote(os.path.abspath(value))}"
    assignment_tests = "\n".join(
        f'test -s {shlex.quote(os.path.join(final_assignments, f"lib{value}.reconciled.assignments"))}'
        for value in lib_nums)
    body = f'''mkdir -p {shlex.quote(aggregate)} {shlex.quote(final_assignments)}
echo "IDENTITY_FINALIZE: assignment-noninterfering evidence joins and canonical ledgers"
python3 -B {shlex.quote(helper)} finalize \\
  --libraries {libraries} \\
  --demux-root {shlex.quote(get_demux_mapping_root())} \\
  --candidate-root {shlex.quote(get_identity_subdir(args, "candidates"))} \\
  --decisions-root {shlex.quote(get_identity_subdir(args, "decisions"))} \\
  --evidence-mode {shlex.quote(args.identity_evidence_mode)} \\
  --run-id {shlex.quote(ORCHESTRATOR_RELEASE)} \\
  --output-root {shlex.quote(aggregate)}{optional}
for output in \\
  identity_reconciliation_final_cells.tsv.gz \\
  identity_reconciliation_candidate_audit.tsv.gz \\
  identity_reconciliation_final_events.tsv \\
  identity_reconciliation_review_queue.tsv.gz \\
  identity_reconciliation_run_summary.tsv; do
  test -s {shlex.quote(aggregate)}/"$output"
done
{assignment_tests}
'''
    return _write_identity_sbatch(
        "identity_finalize", body, dep_job_ids, cpus=1, mem="64G",
        modules=("miniforge/3",),
        commands=("python3",))


def _identity_candidate_axis_outputs_complete(args, eligible_libraries):
    if not eligible_libraries:
        return False
    paths = _candidate_axis_paths(args)
    required = (
        "candidate_axis_cell_scores.tsv.gz",
        "candidate_axis_event_summary.tsv",
        "candidate_axis_event_strata_summary.tsv",
        "candidate_axis_event_scope_summary.tsv",
        "candidate_axis_applied_vs_retained_contrast.tsv",
        "candidate_axis_pair_exclusions.tsv.gz",
        "candidate_axis_run_summary.tsv",
    )
    if not all(check_file_exists(os.path.join(paths["aggregate"], name))
               for name in required):
        return False
    passing_libraries = set()
    run_summary_path = os.path.join(
        paths["aggregate"], "candidate_axis_run_summary.tsv")
    with _candidate_axis_open_tsv(run_summary_path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not {"library", "accounting_status"} <= set(reader.fieldnames or ()):
            return False
        for row in reader:
            match = re.fullmatch(
                r"(?:lib)?(\d+)", str(row.get("library", "")).strip(),
                flags=re.I)
            if (match and
                    str(row.get("accounting_status", "")).strip().upper()
                    == "PASS"):
                passing_libraries.add(int(match.group(1)))
    if not set(eligible_libraries) <= passing_libraries:
        return False
    expected_keys = set()
    for lib_num in eligible_libraries:
        expected_keys.update(_candidate_axis_finalized_event_keys(lib_num, args))
    represented_keys = set()
    for name in (
            "candidate_axis_cell_scores.tsv.gz",
            "candidate_axis_pair_exclusions.tsv.gz"):
        path = os.path.join(paths["aggregate"], name)
        with _candidate_axis_open_tsv(path) as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            required_fields = {
                "library", "selected_supported_event_id",
                "selected_supported_event_proposal",
            }
            if not required_fields <= set(reader.fieldnames or ()):
                return False
            for row in reader:
                match = re.fullmatch(
                    r"(?:lib)?(\d+)", str(row.get("library", "")).strip(),
                    flags=re.I)
                event_id = str(
                    row.get("selected_supported_event_id", "")).strip()
                proposal = _canonical_identity(
                    row.get("selected_supported_event_proposal", ""))
                if match and event_id and proposal:
                    represented_keys.add(
                        (f"lib{int(match.group(1))}", event_id, proposal))
    return bool(expected_keys) and expected_keys <= represented_keys


def _identity_worker_configure(args):
    global AUDIT_ROOT, HYBRID_ROOT, MT_FUSION_ROOT
    global IDENTITY_RECONCILIATION_ROOT, TETRA_REFINE_ROOT, PLOIDY_CALLS_ROOT
    global EXPECTED_POOL_METADATA, ALLOWED_IDENTITIES, PANEL_METADATA
    global DEMUX_OUTPUT_ROOT, CONDF_DIR, CONDF_PATHS
    global IDENTITY_AMBIENT_CANDIDATE_SET
    AUDIT_ROOT = args.audit_root
    HYBRID_ROOT = args.hybrid_root
    MT_FUSION_ROOT = os.path.abspath(args.mt_output_root)
    IDENTITY_RECONCILIATION_ROOT = os.path.abspath(
        args.identity_reconciliation_root)
    TETRA_REFINE_ROOT = args.refined_assignments_root
    PLOIDY_CALLS_ROOT = args.ploidy_calls_root
    EXPECTED_POOL_METADATA = args.expected_pool_metadata
    ALLOWED_IDENTITIES = args.allowed_identities
    PANEL_METADATA = os.path.abspath(args.panel_metadata)
    DEMUX_OUTPUT_ROOT = (
        os.path.abspath(args.demux_output_root)
        if args.demux_output_root else None)
    CONDF_DIR = os.path.abspath(args.condf_dir)
    CONDF_PATHS = {
        "interindiv_20M": os.path.join(CONDF_DIR, "demux_input_20M.condf"),
        "interindiv_het_10M": os.path.join(
            CONDF_DIR, "demux_input_HET_10M.condf"),
        "species_20M": os.path.join(CONDF_DIR, "demux_input_species_20M.condf"),
    }
    IDENTITY_AMBIENT_CANDIDATE_SET = "applied"
    args.reconciliation_candidate_set = "applied"
    args.contam_assignment_source = IDENTITY_AMBIENT_SELECTOR
    args.identity_candidate_axis_event_id = None
    args.identity_candidate_axis_proposal = None
    if not args.identity_candidate_axis_output_root:
        args.identity_candidate_axis_output_root = os.path.join(
            IDENTITY_RECONCILIATION_ROOT, "candidate_axis")


def _identity_selected_event_libraries(lib_nums, args):
    """Return event-bearing libraries restricted to the explicit selection."""
    selected = set(int(value) for value in lib_nums)
    with _candidate_axis_open_tsv(get_identity_event_path(args)) as handle:
        return {
            int(match.group(1))
            for row in csv.DictReader(handle, delimiter="\t")
            for match in [re.fullmatch(
                r"(?:lib)?(\d+)", str(row.get("library", "")).strip(),
                flags=re.I)]
            if match and int(match.group(1)) in selected
        }


def _identity_four_arm_plot_outputs_missing(bundle):
    """Return selected libraries lacking a completed contam.R task product."""
    missing = set()
    try:
        with open(
                bundle["r_manifest"], "r", encoding="utf-8",
                newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
    except (KeyError, OSError, ValueError, csv.Error):
        return set()
    for row in rows:
        try:
            lib_num = int(row["library"])
            prefix = row["output_prefix"]
            status = os.path.join(
                bundle["plot_root"], "contam_r", "status",
                f"{row['condition_slug']}.lib{lib_num}.tsv")
            if not all(check_file_exists(prefix + suffix) for suffix in (
                    ".contam.pdf", ".contam.png", ".contam.stats")):
                missing.add(lib_num)
                continue
            with open(
                    status, "r", encoding="utf-8", newline="") as handle:
                status_rows = list(csv.DictReader(handle, delimiter="\t"))
            if (len(status_rows) != 1 or
                    status_rows[0].get("status", "").upper() != "PASS"):
                missing.add(lib_num)
        except (KeyError, OSError, ValueError, csv.Error):
            if row.get("library") and str(row["library"]).isdigit():
                missing.add(int(row["library"]))
    return missing


def identity_final_evidence_only_worker(args, lib_nums):
    """Submit subset aggregates and finalization from completed evidence only."""
    condition = COND_BY_ABBREV[AMBIENT_PLOT_DEFAULT_CONDITION]
    missing = {}

    eligible_axis = [
        lib_num for lib_num in lib_nums
        if _candidate_axis_has_finalized_event(lib_num, args)
    ]
    candidate_axis_root = _candidate_axis_paths(args)["root"]
    if eligible_axis and not _identity_candidate_axis_outputs_complete(
            args, eligible_axis):
        for lib_num in eligible_axis:
            missing.setdefault(lib_num, []).append("candidate-axis aggregate")
    if not eligible_axis:
        candidate_axis_root = ""

    event_libraries = _identity_selected_event_libraries(lib_nums, args)
    four_arm_root = ""
    four_arm_bundle = None
    missing_four_arm_r = set()
    if event_libraries:
        four_arm_ready = True
        for lib_num in sorted(event_libraries):
            try:
                context = prepare_identity_ambient_comparison(
                    lib_num, candidate_set="applied")
            except (FileNotFoundError, OSError, ValueError):
                missing.setdefault(lib_num, []).append("four-arm inputs")
                four_arm_ready = False
                continue
            for arm_key in IDENTITY_AMBIENT_ARM_ORDER:
                if (arm_key == "reconciled_replacement" and
                        not context["replacement_arm_eligible"]):
                    continue
                if not check_output_exists(
                        lib_num, condition["abbrev"], arm_key):
                    missing.setdefault(lib_num, []).append(
                        f"four-arm {arm_key}")
        if four_arm_ready:
            four_arm_root = get_ambient_plot_run_dir(args, [condition])
            four_arm_bundle = ambient_generate_plot_job_bundle(
                orchestrator_path=os.path.abspath(__file__),
                mapping_root=get_demux_mapping_root(),
                aggregate_root=AGGREGATE_ROOT,
                plot_root=four_arm_root,
                libraries=sorted(event_libraries),
                conditions=[condition["abbrev"]],
                assignment_sources=IDENTITY_AMBIENT_ARM_ORDER,
                script_dir=get_script_dir(),
                log_dir=os.path.join(
                    get_log_dir(), "IDENTITY_RECONCILIATION", "four_arm"),
                partition=SLURM_PARTITION,
                identity_ambient_candidate_set="applied",
            )
            missing_four_arm_r = _identity_four_arm_plot_outputs_missing(
                four_arm_bundle)

    frozen_ambient_root = ""
    frozen_ambient_bundle = None
    try:
        discovery = discover_ambient_swap_events(args)
    except ValueError as exc:
        if "no supported exact-identity" not in str(exc):
            raise
        discovery = None
    if discovery is not None:
        plans = {}
        for lib_num in lib_nums:
            try:
                plan = prepare_ambient_swap_test_plan(
                    lib_num, args, discovery=discovery)
            except (FileNotFoundError, OSError, ValueError):
                missing.setdefault(lib_num, []).append(
                    "controlled ambient inputs")
                continue
            plans[lib_num] = plan
            if not plan.get("applicable"):
                continue
            if not ambient_validation_profile_outputs_complete(plan):
                missing.setdefault(lib_num, []).append(
                    "controlled ambient profile")
            if not ambient_swap_arm_outputs_complete(
                    lib_num, condition, plan):
                missing.setdefault(lib_num, []).append(
                    "controlled ambient arms")
        if any(plan.get("applicable") for plan in plans.values()):
            frozen_ambient_bundle = generate_ambient_swap_aggregate_script(
                args, [condition], plans, discovery, arm_job_ids=None)
            frozen_ambient_root = frozen_ambient_bundle["output_root"]

    if missing:
        detail = "; ".join(
            f"lib{lib_num} ({', '.join(dict.fromkeys(reasons))})"
            for lib_num, reasons in sorted(missing.items()))
        raise ValueError(
            "IDENTITY_FINAL_EVIDENCE_ONLY selected-library evidence "
            "incomplete: " + detail)

    aggregate_jobs = []
    if four_arm_bundle is not None:
        if missing_four_arm_r:
            four_arm_r_job = _ambient_submit_sbatch(
                four_arm_bundle["r_sbatch"], dependency_job_ids=None)
            print(
                "IDENTITY_FINAL_EVIDENCE_ONLY: selected-library contam.R "
                f"array job {four_arm_r_job}; independent of aggregation")
        aggregate_jobs.append(submit_job(four_arm_bundle["aggregate_sbatch"]))
    if frozen_ambient_bundle is not None:
        aggregate_jobs.append(submit_job(frozen_ambient_bundle["script"]))

    finalizer = generate_identity_finalize_script(
        lib_nums, args, candidate_axis_root=candidate_axis_root,
        frozen_ambient_root=frozen_ambient_root,
        four_arm_root=four_arm_root, dep_job_ids=aggregate_jobs)
    finalizer_job = submit_job(finalizer)
    print(
        "IDENTITY_FINAL_EVIDENCE_ONLY: submitted fresh selected-library "
        f"aggregates and canonical finalizer job {finalizer_job}")
    return 0


def resolve_identity_finalize_only_evidence_roots(
        lib_nums, args, require_completed=False):
    """Resolve completed evidence roots without generating or submitting work."""
    condition = COND_BY_ABBREV[AMBIENT_PLOT_DEFAULT_CONDITION]
    path_args = argparse.Namespace(**vars(args))
    path_args.contam_assignment_source = IDENTITY_AMBIENT_SELECTOR
    path_args.reconciliation_candidate_set = "applied"
    path_args.identity_candidate_axis_output_root = os.path.join(
        os.path.abspath(path_args.identity_reconciliation_root),
        "candidate_axis")

    eligible_axis = [
        lib_num for lib_num in lib_nums
        if _candidate_axis_has_finalized_event(lib_num, path_args)
    ]
    candidate_axis_root = (
        _candidate_axis_paths(path_args)["root"] if eligible_axis else "")
    if (require_completed and eligible_axis and
            not _identity_candidate_axis_outputs_complete(
                path_args, eligible_axis)):
        raise ValueError(
            "completed candidate-axis aggregate does not cover every "
            "selected finalized event/proposal")

    event_libraries = _identity_selected_event_libraries(
        lib_nums, path_args)
    four_arm_root = (
        get_ambient_plot_run_dir(path_args, [condition])
        if event_libraries else "")
    if require_completed and four_arm_root:
        four_required = (
            "reconciliation_planned_contrast_cells.tsv",
            "reconciliation_planned_contrasts.tsv",
            "reconciliation_exact_donor_burden.tsv",
        )
        missing_four = [
            name for name in four_required
            if not check_file_exists(os.path.join(
                four_arm_root, "data", name))]
        if missing_four:
            raise ValueError(
                "completed four-arm aggregate is missing: "
                + ", ".join(missing_four))

    frozen_ambient_root = ""
    try:
        discovery = discover_ambient_swap_events(path_args)
    except ValueError as exc:
        if "no supported exact-identity" not in str(exc):
            raise
        discovery = None
    if (discovery is not None and
            set(int(value) for value in lib_nums) &
            set(discovery["by_library"])):
        discovery_id = ambient_swap_discovery_id(discovery)
        candidate_root = os.path.join(
            get_ambient_swap_test_output_root(
                path_args, [condition], discovery_id),
            "identity_reconciliation_figures")
        controlled_required = (
            "ambient_swap_cell_contrasts.tsv.gz",
            "ambient_swap_candidate_cells.tsv",
            "ambient_swap_library_applicability.tsv",
        )
        controlled_paths = [
            os.path.join(candidate_root, "data", name)
            for name in controlled_required]
        complete = all(check_file_exists(path) for path in controlled_paths)
        partial = any(os.path.lexists(path) for path in controlled_paths)
        if complete:
            frozen_ambient_root = candidate_root
        elif require_completed and partial:
            raise ValueError(
                "controlled ambient aggregate is incomplete: "
                + ", ".join(
                    name for name, path in zip(
                        controlled_required, controlled_paths)
                    if not check_file_exists(path)))

    return {
        "candidate_axis_root": candidate_axis_root,
        "four_arm_root": four_arm_root,
        "frozen_ambient_root": frozen_ambient_root,
    }


def identity_final_evidence_worker(payload):
    """Expand validated, data-dependent reconciliation evidence work.

    This is an internal scheduler bridge, not a user-facing scientific stage.
    It exists so zero-event libraries never receive artificial candidate-axis
    or reconciliation-specific ambient jobs.
    """
    data = json.loads(payload)
    args = argparse.Namespace(**data["args"])
    lib_nums = [int(value) for value in data["libraries"]]
    _identity_worker_configure(args)
    if data.get("reuse_only"):
        return identity_final_evidence_only_worker(args, lib_nums)
    condition = COND_BY_ABBREV[AMBIENT_PLOT_DEFAULT_CONDITION]

    eligible_axis = [
        lib_num for lib_num in lib_nums
        if _candidate_axis_has_finalized_event(lib_num, args)
    ]
    candidate_axis_root = _candidate_axis_paths(args)["root"]
    candidate_axis_aggregate_job = None
    completed_axis = (
        _identity_candidate_axis_outputs_complete(args, eligible_axis)
        if eligible_axis else False)
    if eligible_axis and not completed_axis:
        if os.path.lexists(candidate_axis_root):
            raise ValueError(
                "candidate-axis root exists but is incomplete or does not "
                "contain every current finalized event/proposal key; use a "
                "new empty --identity-candidate-axis-output-root instead of "
                f"overwriting frozen evidence: {candidate_axis_root}")
        score_jobs = []
        for lib_num in eligible_axis:
            pair_script = generate_identity_candidate_axis_pairs_script(
                lib_num, args)
            pair_job = submit_job(pair_script)
            score_script = generate_identity_candidate_axis_score_script(
                lib_num, args, [pair_job])
            score_jobs.append(submit_job(score_script))
        aggregate_script = generate_identity_candidate_axis_aggregate_script(
            eligible_axis, args, score_jobs)
        candidate_axis_aggregate_job = submit_job(aggregate_script)
        print(
            "IDENTITY_FINAL_EVIDENCE: candidate axis submitted for "
            f"{len(eligible_axis)} event-bearing eligible libraries; "
            f"aggregate job {candidate_axis_aggregate_job}")
    elif eligible_axis:
        print(
            "IDENTITY_FINAL_EVIDENCE: reusing compatible completed candidate "
            f"axis root {candidate_axis_root}")
    else:
        candidate_axis_root = ""
        print("IDENTITY_FINAL_EVIDENCE: no eligible candidate-axis events")

    event_libraries = _identity_selected_event_libraries(lib_nums, args)

    four_arm_root = ""
    four_arm_aggregate_job = None
    if event_libraries:
        arm_jobs = []
        for lib_num in sorted(event_libraries):
            context = prepare_identity_ambient_comparison(
                lib_num, candidate_set="applied")
            for arm_key in IDENTITY_AMBIENT_ARM_ORDER:
                if (arm_key == "reconciled_replacement" and
                        not context["replacement_arm_eligible"]):
                    continue
                if (not args.force and check_output_exists(
                        lib_num, condition["abbrev"], arm_key)):
                    print(
                        "IDENTITY_FINAL_EVIDENCE: reusing compatible "
                        f"four-arm output {arm_key}/lib{lib_num}")
                    continue
                out_prefix = get_contam_prefix(
                    lib_num, condition["abbrev"], arm_key)
                stale_arm_output = any(
                    os.path.lexists(out_prefix + suffix)
                    for suffix in (
                        ".contam_rate", ".contam_prof", ".allele_ratio",
                        ".contam_diagnostics.tsv",
                        ".cell_source_profile.tsv",
                        ".profile_fit_diagnostics.tsv", ".condf_coverage.tsv",
                        ".run_contract.json", ".decontam.assignments",
                        ".identity_ambient_arm.tsv",
                        ".geometry_gate_audit.tsv",
                    ))
                script = generate_contam_script(
                    lib_num, condition,
                    force=(args.force or stale_arm_output),
                    assignment_source=arm_key)
                arm_jobs.append(submit_job(script))
        four_arm_root = get_ambient_plot_run_dir(args, [condition])
        bundle = ambient_generate_plot_job_bundle(
            orchestrator_path=os.path.abspath(__file__),
            mapping_root=get_demux_mapping_root(),
            aggregate_root=AGGREGATE_ROOT,
            plot_root=four_arm_root,
            libraries=sorted(event_libraries),
            conditions=[condition["abbrev"]],
            assignment_sources=IDENTITY_AMBIENT_ARM_ORDER,
            script_dir=get_script_dir(),
            log_dir=os.path.join(
                get_log_dir(), "IDENTITY_RECONCILIATION", "four_arm"),
            partition=SLURM_PARTITION,
            identity_ambient_candidate_set="applied",
        )
        submitted = ambient_submit_plot_job_bundle(
            bundle, upstream_job_ids=arm_jobs, submit=True)
        four_arm_aggregate_job = submitted["aggregate_job_id"]
        print(
            "IDENTITY_FINAL_EVIDENCE: A/B/C and applicable D ambient work "
            f"submitted for {len(event_libraries)} event-bearing libraries; "
            f"aggregate job {four_arm_aggregate_job}")
    else:
        print("IDENTITY_FINAL_EVIDENCE: all selected libraries are zero-event")

    frozen_ambient_root = ""
    frozen_ambient_aggregate_job = None
    try:
        discovery = discover_ambient_swap_events(args)
    except ValueError as exc:
        if "no supported exact-identity" not in str(exc):
            raise
        discovery = None
    if discovery is not None:
        plans = {}
        profile_jobs = {}
        arm_jobs = []
        for lib_num in lib_nums:
            plan = prepare_ambient_swap_test_plan(
                lib_num, args, discovery=discovery)
            plans[lib_num] = plan
            if not plan.get("applicable"):
                continue
            if (not args.force and
                    ambient_validation_profile_outputs_complete(plan)):
                print(
                    "IDENTITY_FINAL_EVIDENCE: reusing compatible "
                    f"ambient-swap profile lib{lib_num}")
            else:
                profile_script = generate_ambient_swap_profile_script(
                    lib_num, plan, force=args.force)
                profile_jobs[lib_num] = submit_job(profile_script)
        for lib_num, plan in plans.items():
            if not plan.get("applicable"):
                continue
            if (not args.force and
                    ambient_validation_profile_outputs_complete(plan) and
                    ambient_swap_arm_outputs_complete(
                        lib_num, condition, plan)):
                print(
                    "IDENTITY_FINAL_EVIDENCE: reusing compatible "
                    f"ambient-swap arms lib{lib_num}")
                continue
            arm_script = generate_ambient_swap_arm_script(
                lib_num, condition, plan,
                dep_job_id=profile_jobs.get(lib_num), force=args.force)
            arm_jobs.append(submit_job(arm_script))
        if any(plan.get("applicable") for plan in plans.values()):
            bundle = generate_ambient_swap_aggregate_script(
                args, [condition], plans, discovery, arm_job_ids=arm_jobs)
            frozen_ambient_root = bundle["output_root"]
            frozen_ambient_aggregate_job = submit_job(bundle["script"])
            print(
                "IDENTITY_FINAL_EVIDENCE: controlled same-profile current/"
                "proposal ambient work submitted; aggregate job "
                f"{frozen_ambient_aggregate_job}")

    dependencies = [
        job for job in (
            candidate_axis_aggregate_job, four_arm_aggregate_job,
            frozen_ambient_aggregate_job)
        if job
    ]
    finalizer = generate_identity_finalize_script(
        lib_nums, args, candidate_axis_root=candidate_axis_root,
        frozen_ambient_root=frozen_ambient_root,
        four_arm_root=four_arm_root, dep_job_ids=dependencies)
    finalizer_job = submit_job(finalizer)
    print(f"IDENTITY_FINAL_EVIDENCE: canonical finalizer job {finalizer_job}")
    return 0


def generate_identity_final_evidence_planner_script(
        lib_nums, args, dep_job_ids=None, reuse_only=False):
    worker_args = dict(vars(args))
    if args.identity_evidence_mode == "rna":
        for field in (
                "identity_atac_root", "identity_atac_bam_template",
                "atac_panel", "atac_het_panel", "atac_barcode_map",
                "rna_barcode_whitelist", "atac_barcode_whitelist"):
            worker_args[field] = None
    payload = json.dumps({
        "libraries": list(lib_nums),
        "args": worker_args,
        "reuse_only": bool(reuse_only),
    }, sort_keys=True, separators=(",", ":"))
    stage_name = (
        "IDENTITY_FINAL_EVIDENCE_ONLY" if reuse_only
        else "IDENTITY_FINAL_EVIDENCE")
    body = f'''echo "{stage_name}: enumerate validated events and submit only applicable evidence jobs"
python3 -B {shlex.quote(os.path.abspath(__file__))} \\
  --_identity-final-evidence-worker {shlex.quote(payload)}
'''
    return _write_identity_sbatch(
        stage_name.lower(), body, dep_job_ids,
        cpus=1, mem="4G",
        modules=("miniforge/3", "genomics-base/latest"),
        commands=("python3", "sbatch"))


# =============================================================================
# AMBIENT_SWAP_TEST: automatically discovered fixed-profile candidate evaluation
# =============================================================================

_AMBIENT_SWAP_DISCOVERY_CACHE = {}
AMBIENT_SWAP_SUPPORTED_EVENT_CLASSES = {
    "LIKELY_UNEXPECTED_INTACT_BIOLOGICAL_LINE",
    "LIKELY_UNEXPECTED_SINGLET_POPULATION",
}
AMBIENT_SWAP_SUPPORTED_EVENT_CONFIDENCES = {"STRONG", "DECISIVE"}


def get_identity_event_path(args):
    return os.path.join(
        get_identity_subdir(args, "decisions"),
        "all_libraries.identity_events.tsv")


def get_identity_reconciliation_manifest_path(args):
    return os.path.join(
        get_identity_subdir(args, "decisions"),
        "reconciliation_manifest.json")


def discover_ambient_swap_events(args):
    """Discover supported unexpected populations from reconciliation output.

    The reconciliation event table is authoritative.  This stage does not
    reproduce nuclear thresholds and does not consume a hand-maintained list.
    It selects only exact-identity, library-unexpected biological populations
    that the upstream occupancy-aware policy classified STRONG or DECISIVE.
    """
    event_path = get_identity_event_path(args)
    manifest_path = get_identity_reconciliation_manifest_path(args)
    for path in (event_path, manifest_path):
        if not check_file_exists(path):
            raise FileNotFoundError(
                f"ambient swap discovery input missing: {path}")
    cache_key = (
        event_path, os.path.getsize(event_path), os.stat(event_path).st_mtime_ns,
        manifest_path, os.path.getsize(manifest_path),
        os.stat(manifest_path).st_mtime_ns,
    )
    cached = _AMBIENT_SWAP_DISCOVERY_CACHE.get(cache_key)
    if cached is not None:
        return {
            **cached,
            "by_library": {
                library: [dict(row) for row in rows]
                for library, rows in cached["by_library"].items()
            },
        }
    with open(manifest_path, "r", encoding="utf-8") as handle:
        manifest = json.load(handle)
    if not manifest.get("occupancy_aware_reconciliation"):
        raise ValueError(
            "AMBIENT_SWAP_TEST requires occupancy-aware reconciliation output")
    if not manifest.get("library_exchange_inference"):
        raise ValueError(
            "AMBIENT_SWAP_TEST requires reconciliation library-exchange inference")
    min_cells = int(manifest.get("event_min_cells", 0))
    if min_cells <= 0:
        raise ValueError(
            "reconciliation manifest lacks a valid event_min_cells policy")

    by_library = {}
    audit_rows = []
    with _open_text_auto(event_path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            "event_id", "library", "event_scope", "unexpected_component",
            "n_implicated_cells", "event_class", "event_confidence",
            "event_reason", "singlet_library_relationship",
        }
        missing = required - set(reader.fieldnames or ())
        if missing:
            raise ValueError(
                f"identity event schema missing {sorted(missing)}: {event_path}")
        for row in reader:
            library_text = str(row.get("library", "")).strip()
            match = re.fullmatch(r"(?:lib)?(\d+)", library_text,
                                 flags=re.IGNORECASE)
            if not match:
                raise ValueError(
                    f"identity event has invalid library: {library_text!r}")
            library = int(match.group(1))
            raw_identity = str(row.get("unexpected_component", "")).strip()
            raw_components = _identity_components(raw_identity)
            identity = _canonical_identity(raw_identity)
            try:
                n_cells = int(float(str(
                    row.get("n_implicated_cells", "0")).strip() or 0))
            except ValueError as exc:
                raise ValueError(
                    f"identity event {row.get('event_id')} has invalid "
                    "n_implicated_cells") from exc
            reasons = []
            if str(row.get("event_scope", "")).strip() != "EXACT_IDENTITY":
                reasons.append("not_exact_identity")
            if str(row.get("event_class", "")).strip() not in \
                    AMBIENT_SWAP_SUPPORTED_EVENT_CLASSES:
                reasons.append("unsupported_event_class")
            if str(row.get("event_confidence", "")).strip() not in \
                    AMBIENT_SWAP_SUPPORTED_EVENT_CONFIDENCES:
                reasons.append("below_supported_event_confidence")
            if n_cells < min_cells:
                reasons.append("below_event_mass_policy")
            if (len(raw_components) > 1 and
                    len(set(raw_components)) == 1):
                reasons.append("homotet_ploidy_event_not_sample_swap")
            if not identity:
                reasons.append("empty_candidate_identity")
            selected = not reasons
            audit = {
                **row,
                "library_num": library,
                "candidate_identity": identity,
                "n_implicated_cells_int": n_cells,
                "ambient_swap_selected": selected,
                "ambient_swap_exclusion_reason": ";".join(reasons) or "NA",
            }
            audit_rows.append(audit)
            if selected:
                by_library.setdefault(library, []).append(audit)
    selected_events = sum(len(rows) for rows in by_library.values())
    if selected_events == 0:
        raise ValueError(
            "reconciliation contains no supported exact-identity unexpected "
            "events for ambient swap testing")
    result = {
        "event_path": event_path,
        "event_record": _identity_ambient_file_record(event_path),
        "manifest_path": manifest_path,
        "manifest_record": _identity_ambient_file_record(manifest_path),
        "manifest": manifest,
        "by_library": by_library,
        "audit_rows": audit_rows,
        "selected_event_count": selected_events,
        "selected_library_count": len(by_library),
        "event_min_cells": min_cells,
    }
    _AMBIENT_SWAP_DISCOVERY_CACHE[cache_key] = result
    return {
        **result,
        "by_library": {
            library: [dict(row) for row in rows]
            for library, rows in by_library.items()
        },
    }


def ambient_swap_discovery_id(discovery):
    policy = re.sub(
        r"[^A-Za-z0-9_.-]+", "_",
        str(discovery["manifest"].get("policy_version", "policy")))
    generation_basis = "\n".join((
        discovery["event_record"]["sha256"],
        AMBIENT_SWAP_TEST_PLAN_VERSION,
        f"heterotypic_nn_min={AMBIENT_SWAP_HETEROTYPIC_NN_MIN_PROB:.6f}",
        "heterotypic_nn_qc_required=true",
        "post_gate_event_mass_required=true",
    ))
    generation_hash = hashlib.sha256(
        generation_basis.encode("utf-8")).hexdigest()[:10]
    return f"auto_{policy[:40]}_{generation_hash}"


def get_ambient_swap_test_root(lib_num, discovery_id):
    return os.path.join(
        get_demux_dir(lib_num), AMBIENT_SWAP_TEST_DIRNAME, discovery_id)


def get_ambient_swap_test_arm_prefix(
        lib_num, condition, arm_key, discovery_id):
    if arm_key not in AMBIENT_SWAP_TEST_ARMS:
        raise ValueError(f"unknown ambient swap-test arm: {arm_key}")
    return os.path.join(
        get_ambient_swap_test_root(lib_num, discovery_id), "conditions",
        str(condition), arm_key, f"lib{lib_num}_demuxed")


def get_ambient_swap_test_output_root(args, cond_list, discovery_id):
    conditions = [condition["abbrev"] for condition in cond_list]
    short_conditions = [
        AMBIENT_CONDITION_SHORT_NAMES.get(condition, condition)
        for condition in conditions]
    if len(conditions) == 1:
        label = _ambient_condition_slug(
            f"{short_conditions[0]}__ambient-swap-test")
    else:
        digest = hashlib.sha1(
            "\n".join(conditions).encode("utf-8")).hexdigest()[:10]
        label = f"selected_{len(conditions)}__ambient-swap-test_{digest}"
    return os.path.join(
        os.path.abspath(args.ambient_plot_root), label, discovery_id)


def _read_reconciliation_decision_rows(path):
    with _open_text_auto(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            "barcode", "original_demux_assignment",
            "current_refined_assignment", "current_ploidy",
            "proposed_donor_genotype", "proposed_droplet_state",
            "proposed_biological_ploidy",
            "proposed_biological_admissibility", "final_action",
            "decision_confidence", "reassignment_applied",
            "nn_prob_tetraploid", "nn_qc_pass",
            "explicit_multiplet_evidence", "occupancy_resolution_status",
            "library_exchange_evidence_eligible", "event_id",
        }
        missing = required - set(reader.fieldnames or ())
        if missing:
            raise ValueError(
                f"reconciliation decision schema missing {sorted(missing)}: "
                f"{path}")
        rows = {}
        for row in reader:
            barcode = str(row.get("barcode", "")).strip()
            if not barcode:
                raise ValueError(
                    f"reconciliation decision row has empty barcode: {path}")
            if barcode in rows:
                raise ValueError(
                    f"duplicate reconciliation barcode {barcode}: {path}")
            rows[barcode] = row
    return rows


def _ambient_swap_diploid_to_heterotypic(row):
    """Return whether one proposal crosses diploid -> heterotypic tetraploid."""
    proposed_components = _identity_components(
        row.get("proposed_donor_genotype", ""))
    if (len(proposed_components) != 2 or
            len(set(proposed_components)) != 2):
        return False

    current_ploidy = str(row.get("current_ploidy", "")).strip().upper()
    if current_ploidy in {
            "T", "TET", "TETRAPLOID", "HOMOTYPIC", "HETEROTYPIC"}:
        return False
    if current_ploidy in {"D", "DIP", "DIPLOID"}:
        return True

    current_identity = (
        row.get("current_refined_assignment", "") or
        row.get("original_demux_assignment", ""))
    current_components = _identity_components(current_identity)
    return len(set(current_components)) == 1


def _ambient_swap_row_eligibility(row):
    """Return the upstream-policy basis for testing one proposed assignment."""
    proposed = _canonical_identity(row.get("proposed_donor_genotype", ""))
    original = _canonical_identity(
        row.get("current_refined_assignment", "") or
        row.get("original_demux_assignment", ""))
    if not proposed or proposed == original:
        return False, "no_identity_change"
    if str(row.get("proposed_droplet_state", "")).strip() != "SINGLE_CELL":
        return False, "not_single_cell"
    if str(row.get("proposed_biological_admissibility", "")).strip() not in {
            "BIOLOGICAL_SINGLE_CELL_ALLOWED", "SINGLET_IDENTITY_CANDIDATE"}:
        return False, "biologically_inadmissible"
    if _truthy(row.get("explicit_multiplet_evidence")):
        return False, "explicit_multiplet_evidence"
    if str(row.get("occupancy_resolution_status", "")).strip() in {
            "DONOR_PAIR_CELLULAR_ORIGIN_UNRESOLVED",
            "HOMOTET_VS_SAME_DONOR_DOUBLET_UNRESOLVED"}:
        return False, "occupancy_unresolved"
    if _ambient_swap_diploid_to_heterotypic(row):
        if str(row.get("proposed_biological_ploidy", "")).strip().upper() \
                != "TETRAPLOID":
            return False, "heterotypic_proposal_lacks_tetraploid_metadata"
        if not _truthy(row.get("nn_qc_pass")):
            return False, "diploid_to_heterotypic_nn_qc_failed"
        try:
            nn_prob_tetraploid = float(str(
                row.get("nn_prob_tetraploid", "")).strip())
        except (TypeError, ValueError):
            return False, "diploid_to_heterotypic_nn_missing"
        if (not math.isfinite(nn_prob_tetraploid) or
                nn_prob_tetraploid < AMBIENT_SWAP_HETEROTYPIC_NN_MIN_PROB):
            return False, "diploid_to_heterotypic_nn_below_threshold"
    if _truthy(row.get("reassignment_applied")):
        return True, "upstream_reassignment_applied"
    if _truthy(row.get("library_exchange_evidence_eligible")):
        return True, "upstream_library_exchange_eligible"
    confidence = str(row.get("decision_confidence", "")).strip()
    if confidence == "DECISIVE":
        return True, "upstream_decisive_not_applied"
    if confidence == "STRONG_NOT_AUTOAPPLIED":
        return True, "upstream_strong_not_autoapplied"
    return False, "not_supported_by_upstream_policy"


def prepare_ambient_swap_test_plan(lib_num, args, discovery=None):
    """Build one active-current-versus-proposed overlay per library."""
    discovery = discovery or discover_ambient_swap_events(args)
    discovery_id = ambient_swap_discovery_id(discovery)
    requested_rows = list(discovery["by_library"].get(int(lib_num), ()))
    requested = {
        _canonical_identity(row["candidate_identity"]): row
        for row in requested_rows
    }
    if len(requested) != len(requested_rows):
        raise ValueError(
            f"lib{lib_num} reconciliation emitted duplicate supported exact-"
            "identity events")
    if not requested:
        return {
            "schema_version": 1,
            "plan_version": AMBIENT_SWAP_TEST_PLAN_VERSION,
            "library": int(lib_num),
            "discovery_id": discovery_id,
            "event_path": discovery["event_path"],
            "applicable": False,
            "skip_reason": (
                "reconciliation found no supported unexpected exact-identity "
                "event in this library"),
            "candidate_identities": [],
        }

    demux_prefix = get_demux_prefix(lib_num)
    original_assignments_path = demux_prefix + ".assignments"
    samples_path = demux_prefix + ".samples"
    decisions_path = get_reconciled_cells_path(lib_num)
    original_receivers_path = get_expected_lines(lib_num)
    original_donors_path = get_individual_ambient_candidates(lib_num)
    required = (
        original_assignments_path, samples_path, decisions_path,
        original_receivers_path, original_donors_path,
    )
    missing = [path for path in required if not check_file_exists(path)]
    if missing:
        raise FileNotFoundError(
            "ambient swap test requires completed demux and reconciliation; "
            "missing: " + ", ".join(missing))

    demux_rows = _read_assignment_rows(original_assignments_path)
    decision_rows = _read_reconciliation_decision_rows(decisions_path)
    unknown_barcodes = sorted(set(decision_rows) - set(demux_rows))
    if unknown_barcodes:
        raise ValueError(
            "reconciliation ledger contains barcodes absent from demux "
            "assignments (first 10): " + ", ".join(unknown_barcodes[:10]))

    selected = {}
    excluded_counts = {
        identity: {} for identity in requested
    }
    matched_rows = {identity: 0 for identity in requested}
    requested_event_ids = {
        identity: str(row.get("event_id", "")).strip()
        for identity, row in requested.items()
    }
    for barcode, row in decision_rows.items():
        proposed = _canonical_identity(
            row.get("proposed_donor_genotype", ""))
        if proposed not in requested:
            continue
        if (str(row.get("event_id", "")).strip() !=
                requested_event_ids[proposed]):
            continue
        matched_rows[proposed] += 1
        eligible, basis = _ambient_swap_row_eligibility(row)
        if not eligible:
            counts = excluded_counts[proposed]
            counts[basis] = counts.get(basis, 0) + 1
            continue
        recorded_original = _canonical_identity(
            row.get("original_demux_assignment", ""))
        observed_original = _canonical_identity(demux_rows[barcode][1])
        if not recorded_original or recorded_original != observed_original:
            raise ValueError(
                f"lib{lib_num} reconciliation/demux identity mismatch for "
                f"{barcode}: {recorded_original!r} != {observed_original!r}")
        comparison_current = _canonical_identity(
            row.get("current_refined_assignment", "")) or observed_original
        selected[barcode] = {
            "row": row,
            "original_identity": comparison_current,
            "candidate_identity": proposed,
            "inclusion_basis": basis,
        }

    # Event mass was established upstream before this stage's independent
    # cell-level NN/occupancy gates. Reapply the same floor to the surviving
    # cells so a handful of eligible outliers cannot keep an identity alive.
    post_gate_counts = {identity: 0 for identity in requested}
    for item in selected.values():
        post_gate_counts[item["candidate_identity"]] += 1
    below_mass = {
        identity for identity, count in post_gate_counts.items()
        if 0 < count < discovery["event_min_cells"]
    }
    if below_mass:
        for barcode in list(selected):
            identity = selected[barcode]["candidate_identity"]
            if identity not in below_mass:
                continue
            del selected[barcode]
            counts = excluded_counts[identity]
            reason = "below_event_mass_after_cell_eligibility"
            counts[reason] = counts.get(reason, 0) + 1

    if not selected:
        details = "; ".join(
            f"{identity}: matched={matched_rows[identity]}, "
            f"eligible={post_gate_counts[identity]}"
            for identity in sorted(requested))
        return {
            "schema_version": 1,
            "plan_version": AMBIENT_SWAP_TEST_PLAN_VERSION,
            "library": int(lib_num),
            "discovery_id": discovery_id,
            "event_path": discovery["event_path"],
            "applicable": False,
            "skip_reason": (
                "no discovered-event cell passed upstream reconciliation "
                f"eligibility ({details})"),
            "candidate_identities": sorted(requested),
        }

    selected_candidates = sorted({
        item["candidate_identity"] for item in selected.values()})
    original_receivers = _read_identity_lines(original_receivers_path)
    receiver_keys = {_canonical_identity(x) for x in original_receivers}
    receiver_lines = list(original_receivers)
    for identity in selected_candidates:
        if identity not in receiver_keys:
            receiver_lines.append(identity)
            receiver_keys.add(identity)
    production_donors = set(_read_identity_lines(original_donors_path))
    for identity in selected_candidates:
        production_donors.update(_identity_components(identity))
    samples = set(_read_identity_lines(samples_path))
    unknown_donors = sorted(production_donors - samples)
    if unknown_donors:
        raise ValueError(
            f"lib{lib_num} discovered candidate components absent from .samples: "
            + ", ".join(unknown_donors))

    root = get_ambient_swap_test_root(lib_num, discovery_id)
    plan_dir = os.path.join(root, "plan")
    profile_dir = os.path.join(root, "fixed_profile")
    paths = {
        "current_assignments": os.path.join(
            plan_dir, f"lib{lib_num}.comparison_current_overlay.assignments"),
        "candidate_assignments": os.path.join(
            plan_dir, f"lib{lib_num}.candidate_overlay.assignments"),
        "receiver_lines": os.path.join(
            plan_dir, f"lib{lib_num}.receiver_lines.txt"),
        "production_roster": os.path.join(
            plan_dir, f"lib{lib_num}.production_roster.txt"),
        "candidate_cells": os.path.join(
            plan_dir, f"lib{lib_num}.candidate_cells.tsv"),
        "cell_strata": os.path.join(
            plan_dir, f"lib{lib_num}.cell_strata.tsv"),
        "event_discovery": os.path.join(
            plan_dir, f"lib{lib_num}.event_discovery.tsv"),
    }

    current_lines = []
    overlay_lines = []
    selected_originals = {
        item["original_identity"] for item in selected.values()}
    for barcode, demux_row in demux_rows.items():
        decision = decision_rows.get(barcode, {})
        current_identity = (
            _canonical_identity(decision.get("current_refined_assignment", ""))
            or _canonical_identity(demux_row[1]))
        current_type = "D" if "+" in current_identity else "S"
        current_lines.append("\t".join((
            barcode, current_identity, current_type, demux_row[3])))
        identity = (
            selected[barcode]["candidate_identity"]
            if barcode in selected else current_identity)
        assignment_type = "D" if "+" in identity else "S"
        overlay_lines.append("\t".join((
            barcode, identity, assignment_type, demux_row[3])))
    _write_if_changed(
        paths["current_assignments"],
        "\n".join(current_lines) + "\n")
    _write_if_changed(
        paths["candidate_assignments"],
        "\n".join(overlay_lines) + "\n")
    _write_if_changed(
        paths["receiver_lines"],
        "".join(f"{identity}\n" for identity in receiver_lines))
    _write_if_changed(
        paths["production_roster"],
        "".join(f"{identity}\n" for identity in sorted(production_donors)))

    candidate_header = (
        "barcode\toriginal_identity\tcandidate_identity\tevent_id"
        "\tdecision_confidence\tfinal_action\treassignment_applied"
        "\tlibrary_exchange_evidence_eligible\tinclusion_basis"
        "\tcurrent_ploidy\tproposed_biological_ploidy"
        "\tnn_prob_tetraploid\tnn_qc_pass")
    candidate_lines = [candidate_header]
    basis_counts = {identity: {} for identity in requested}
    applied_counts = {identity: 0 for identity in requested}
    selected_counts = {identity: 0 for identity in requested}
    for barcode in sorted(selected):
        item = selected[barcode]
        row = item["row"]
        identity = item["candidate_identity"]
        selected_counts[identity] += 1
        basis = item["inclusion_basis"]
        basis_counts[identity][basis] = basis_counts[identity].get(basis, 0) + 1
        applied_counts[identity] += int(_truthy(row.get("reassignment_applied")))
        candidate_lines.append("\t".join((
            barcode,
            item["original_identity"],
            identity,
            str(row.get("event_id", "")).strip() or "NA",
            str(row.get("decision_confidence", "")).strip() or "NA",
            str(row.get("final_action", "")).strip() or "NA",
            "TRUE" if _truthy(row.get("reassignment_applied")) else "FALSE",
            "TRUE" if _truthy(row.get(
                "library_exchange_evidence_eligible")) else "FALSE",
            basis,
            str(row.get("current_ploidy", "")).strip() or "NA",
            str(row.get("proposed_biological_ploidy", "")).strip() or "NA",
            str(row.get("nn_prob_tetraploid", "")).strip() or "NA",
            "TRUE" if _truthy(row.get("nn_qc_pass")) else "FALSE",
        )))
    _write_if_changed(
        paths["candidate_cells"], "\n".join(candidate_lines) + "\n")

    strata_lines = [
        "barcode\toriginal_identity\tcandidate_identity\tstratum"]
    for barcode, demux_row in demux_rows.items():
        decision = decision_rows.get(barcode, {})
        original = (
            _canonical_identity(decision.get("current_refined_assignment", ""))
            or _canonical_identity(demux_row[1]))
        if barcode in selected:
            candidate = selected[barcode]["candidate_identity"]
            stratum = "candidate_cell"
        elif original in selected_originals:
            candidate = "NA"
            stratum = "source_control"
        else:
            candidate = "NA"
            stratum = "background"
        strata_lines.append("\t".join((
            barcode, original, candidate, stratum)))
    _write_if_changed(
        paths["cell_strata"], "\n".join(strata_lines) + "\n")

    match_lines = [
        "library\tcandidate_identity\tevent_id\tevent_class"
        "\tevent_confidence\tn_implicated_cells\tevent_min_cells_policy"
        "\tmatched_decision_rows\tselected_cells\tapplied_selected_cells"
        "\tinclusion_bases\texclusion_reasons\tevent_reason"]
    for identity in sorted(requested):
        event_row = requested[identity]
        match_lines.append("\t".join((
            f"lib{lib_num}",
            identity,
            str(event_row.get("event_id", "")).strip() or "NA",
            str(event_row.get("event_class", "")).strip() or "NA",
            str(event_row.get("event_confidence", "")).strip() or "NA",
            str(event_row.get("n_implicated_cells_int", 0)),
            str(discovery["event_min_cells"]),
            str(matched_rows[identity]),
            str(selected_counts[identity]),
            str(applied_counts[identity]),
            ";".join(
                f"{key}:{value}" for key, value in
                sorted(basis_counts[identity].items())) or "NA",
            ";".join(
                f"{key}:{value}" for key, value in
                sorted(excluded_counts[identity].items())) or "NA",
            str(event_row.get("event_reason", "")).strip() or "NA",
        )))
    _write_if_changed(
        paths["event_discovery"], "\n".join(match_lines) + "\n")

    raw_prefix = os.path.join(get_demux_dir(lib_num), f"lib{lib_num}_raw")
    inputs = {
        "identity_events": discovery["event_record"],
        "reconciliation_manifest": discovery["manifest_record"],
        "demux_original_assignments": _identity_ambient_file_record(
            original_assignments_path),
        "comparison_current_assignments": _identity_ambient_file_record(
            paths["current_assignments"]),
        "reconciled_cells": _identity_ambient_file_record(decisions_path),
        "raw_counts": _ambient_validation_large_file_record(
            raw_prefix + ".counts"),
        "raw_condf": _identity_ambient_file_record(raw_prefix + ".condf"),
        "samples": _identity_ambient_file_record(samples_path),
        "filtered_barcodes": _identity_ambient_file_record(
            get_filtered_barcodes(lib_num)),
    }
    rosters = {
        key: _identity_ambient_file_record(path)
        for key, path in paths.items()
    }
    payload = {
        "plan_version": AMBIENT_SWAP_TEST_PLAN_VERSION,
        "library": int(lib_num),
        "discovery_id": discovery_id,
        "applicable": True,
        "profile_method": {
            "max_empty": AMBIENT_VALIDATION_PROFILE_MAX_EMPTY,
            "seed": AMBIENT_VALIDATION_PROFILE_SEED,
            "bootstrap": 0,
            "profile_starts": AMBIENT_VALIDATION_PROFILE_STARTS,
            "roster": "original_donors_plus_discovered_candidate_components",
        },
        "cell_eligibility": {
            "diploid_to_heterotypic_nn_gate": True,
            "nn_qc_required": True,
            "nn_prob_tetraploid_min": (
                AMBIENT_SWAP_HETEROTYPIC_NN_MIN_PROB),
            "event_mass_reapplied_after_cell_gates": True,
            "event_min_cells": discovery["event_min_cells"],
        },
        "inputs": inputs,
        "rosters": rosters,
    }
    fingerprint = hashlib.sha256(json.dumps(
        payload, sort_keys=True, separators=(",", ":")).encode(
            "utf-8")).hexdigest()
    profile_tag = fingerprint[:12]
    plan = {
        **payload,
        **paths,
        "schema_version": 1,
        "plan_fingerprint": fingerprint,
        "root": root,
        "plan_dir": plan_dir,
        "profile_dir": profile_dir,
        "event_path": discovery["event_path"],
        "reconciliation_manifest": discovery["manifest_path"],
        "original_assignments": paths["current_assignments"],
        "raw_prefix": raw_prefix,
        "demux_prefix": demux_prefix,
        "fixed_profile": os.path.join(
            profile_dir,
            f"lib{lib_num}.discovered_candidates_50k_{profile_tag}.contam_prof"),
        "fixed_diagnostics": os.path.join(
            profile_dir,
            f"lib{lib_num}.discovered_candidates_50k_{profile_tag}.diagnostics.tsv"),
        "candidate_identities": selected_candidates,
        "requested_candidate_identities": sorted(requested),
        "n_selected_cells": len(selected),
        "n_source_controls": sum(
            1 for line in strata_lines[1:] if line.endswith("\tsource_control")),
        "production_roster_list": sorted(production_donors),
        "skip_reason": "",
    }
    plan_path = os.path.join(
        plan_dir, f"lib{lib_num}.ambient_swap_test_plan.json")
    plan["plan_path"] = plan_path
    _write_if_changed(
        plan_path, json.dumps(plan, indent=2, sort_keys=True) + "\n")
    return plan


# =============================================================================
# AMBIENT_VALIDATE: production fixed-profile assignment validation
# =============================================================================

def get_ambient_validation_root(lib_num, candidate_set="applied"):
    return os.path.join(
        get_demux_dir(lib_num), IDENTITY_AMBIENT_DIRNAME,
        _identity_ambient_candidate_set(candidate_set),
        AMBIENT_VALIDATION_DIRNAME)


def get_ambient_validation_fixed_prefix(
        lib_num, condition, arm_key, candidate_set="applied"):
    if arm_key not in AMBIENT_VALIDATION_FIXED_ARMS:
        raise ValueError(f"unknown ambient validation arm: {arm_key}")
    return os.path.join(
        get_ambient_validation_root(lib_num, candidate_set),
        "conditions", str(condition), arm_key, f"lib{lib_num}_demuxed")


def _ambient_validation_large_file_record(path):
    """Record a large immutable payload without an expensive content scan."""
    path = os.path.abspath(path)
    stat_result = os.stat(path)
    return {
        "path": path,
        "size_bytes": int(stat_result.st_size),
        "mtime_ns": int(stat_result.st_mtime_ns),
        "sha256": "",
        "fingerprint_basis": "size_mtime",
    }


def prepare_ambient_validation_plan(lib_num, candidate_set="applied"):
    """Freeze the applied targets and production donor roster for E/F.

    Libraries without an applied donor addition are valid but not applicable;
    callers skip them rather than failing an all-library production run.
    """
    candidate_set = _identity_ambient_candidate_set(candidate_set)
    if candidate_set != "applied":
        raise ValueError(
            "AMBIENT_VALIDATE always uses the applied reconciliation set")
    comparison = prepare_identity_ambient_comparison(
        lib_num, candidate_set=candidate_set)
    targets = sorted(set(comparison.get("added_candidate_components", ())))
    if not targets:
        return {
            "schema_version": 1,
            "plan_version": AMBIENT_VALIDATION_PLAN_VERSION,
            "library": int(lib_num),
            "candidate_set": candidate_set,
            "applicable": False,
            "skip_reason": "no applied reconciliation added a donor component",
            "targets_list": [],
            "identity_ambient_plan_fingerprint": comparison["plan_fingerprint"],
        }
    production = sorted(set(_read_identity_lines(
        comparison["augmented_ambient_candidates"])))
    samples_path = get_demux_prefix(lib_num) + ".samples"
    samples = _read_identity_lines(samples_path)
    unknown = sorted(set(production) - set(samples))
    if unknown:
        raise ValueError(
            f"lib{lib_num} production donors absent from .samples: "
            + ", ".join(unknown))
    root = get_ambient_validation_root(lib_num, candidate_set)
    plan_dir = os.path.join(root, "plan")
    profile_dir = os.path.join(root, "fixed_profile")
    paths = {
        "targets": os.path.join(plan_dir, f"lib{lib_num}.targets.txt"),
        "production_roster": os.path.join(
            plan_dir, f"lib{lib_num}.production_roster.txt"),
    }
    values = {
        "targets": targets,
        "production_roster": production,
    }
    for key, path in paths.items():
        _write_if_changed(path, "".join(f"{value}\n" for value in values[key]))

    payload = {
        "plan_version": AMBIENT_VALIDATION_PLAN_VERSION,
        "library": int(lib_num),
        "candidate_set": candidate_set,
        "applicable": True,
        "identity_ambient_plan_fingerprint": comparison["plan_fingerprint"],
        "profile_method": {
            "max_empty": AMBIENT_VALIDATION_PROFILE_MAX_EMPTY,
            "seed": AMBIENT_VALIDATION_PROFILE_SEED,
            "bootstrap": 0,
            "profile_starts": AMBIENT_VALIDATION_PROFILE_STARTS,
            "roster": "applied_production_augmented",
        },
        "inputs": {
            "raw_counts": _ambient_validation_large_file_record(os.path.join(
                get_demux_dir(lib_num), f"lib{lib_num}_raw.counts")),
            "raw_condf": _identity_ambient_file_record(os.path.join(
                get_demux_dir(lib_num), f"lib{lib_num}_raw.condf")),
            "samples": _identity_ambient_file_record(samples_path),
            "filtered_barcodes": _identity_ambient_file_record(
                get_filtered_barcodes(lib_num)),
            "scrutiny_cells": _identity_ambient_file_record(
                comparison["scrutiny_cells"]),
        },
        "rosters": {
            key: _identity_ambient_file_record(path)
            for key, path in paths.items()
        },
    }
    fingerprint = hashlib.sha256(json.dumps(
        payload, sort_keys=True, separators=(",", ":")).encode(
            "utf-8")).hexdigest()
    profile_tag = fingerprint[:12]
    plan = {
        **payload,
        **paths,
        "schema_version": 1,
        "plan_fingerprint": fingerprint,
        "root": root,
        "plan_dir": plan_dir,
        "profile_dir": profile_dir,
        "context_path": comparison["context_path"],
        "scrutiny_cells": comparison["scrutiny_cells"],
        "raw_prefix": os.path.join(get_demux_dir(lib_num), f"lib{lib_num}_raw"),
        "demux_prefix": get_demux_prefix(lib_num),
        "fixed_profile": os.path.join(
            profile_dir,
            f"lib{lib_num}.production_50k_{profile_tag}.contam_prof"),
        "fixed_diagnostics": os.path.join(
            profile_dir,
            f"lib{lib_num}.production_50k_{profile_tag}.diagnostics.tsv"),
        "targets_list": targets,
        "production_roster_list": production,
        "skip_reason": "",
    }
    plan_path = os.path.join(
        plan_dir, f"lib{lib_num}.fixed_profile_validation_plan.json")
    plan["plan_path"] = plan_path
    _write_if_changed(
        plan_path, json.dumps(plan, indent=2, sort_keys=True) + "\n")
    return plan


def ambient_validation_profile_outputs_complete(plan):
    if not all(check_file_exists(path) for path in (
            plan["fixed_profile"], plan["fixed_diagnostics"])):
        return False
    metrics = {}
    try:
        with open(plan["fixed_diagnostics"], "r", encoding="utf-8") as handle:
            for row in csv.reader(handle, delimiter="\t"):
                if len(row) >= 2 and row[0] != "metric":
                    metrics[row[0]] = row[1]
        successful_starts = int(float(
            metrics.get("profile_starts_successful", "0")))
    except (OSError, TypeError, ValueError):
        return False
    return (
        successful_starts >= AMBIENT_VALIDATION_PROFILE_STARTS and
        metrics.get("bulk_log_likelihood_valid", "").lower() == "true"
    )


def _ambient_validation_plan_guard_shell(plan):
    return f'''python3 - "{plan['plan_path']}" "{plan['plan_fingerprint']}" <<'PY'
import hashlib
import json
import os
import sys

path, expected_fingerprint = sys.argv[1:]
with open(path, "r", encoding="utf-8") as handle:
    plan = json.load(handle)
if plan.get("plan_fingerprint") != expected_fingerprint:
    raise SystemExit("ERROR: ambient validation plan fingerprint changed")
for section in ("inputs", "rosters"):
    records = plan.get(section, {{}})
    if not records:
        raise SystemExit(f"ERROR: ambient validation plan lacks {{section}} records")
    for name, record in records.items():
        record_path = str(record.get("path", ""))
        if not os.path.isfile(record_path):
            raise SystemExit(f"ERROR: ambient validation {{section}} file missing: {{name}}: {{record_path}}")
        if os.path.getsize(record_path) != int(record.get("size_bytes", -1)):
            raise SystemExit(f"ERROR: ambient validation {{section}} size changed: {{name}}")
        expected_mtime = record.get("mtime_ns")
        if (expected_mtime is not None and
                os.stat(record_path).st_mtime_ns != int(expected_mtime)):
            raise SystemExit(f"ERROR: ambient validation {{section}} mtime changed: {{name}}")
        expected_digest = str(record.get("sha256", ""))
        if not expected_digest:
            continue
        digest = hashlib.sha256()
        with open(record_path, "rb") as source:
            for block in iter(lambda: source.read(1024 * 1024), b""):
                digest.update(block)
        if digest.hexdigest() != expected_digest:
            raise SystemExit(f"ERROR: ambient validation {{section}} digest changed: {{name}}")
PY'''


def generate_ambient_validation_profile_script(lib_num, plan, force=False):
    """Generate the single production-roster empty-drop profile job."""
    binary = os.path.join(SOFTWARE_BIN, "tet_ambient_profile")
    job_name = f"ambV_profile_lib{lib_num}"
    log_dir = os.path.join(
        get_log_dir(), "AMBIENT_VALIDATE", f"lib{lib_num}")
    script_path = os.path.join(get_script_dir(), f"{job_name}.sbatch")
    os.makedirs(log_dir, exist_ok=True)
    os.makedirs(get_script_dir(), exist_ok=True)
    os.makedirs(plan["profile_dir"], exist_ok=True)
    skip_block = ""
    if not force and ambient_validation_profile_outputs_complete(plan):
        skip_block = (
            'echo "Matching fixed empty-drop profile already exists; skipping"\n'
            "exit 0")
    production_run = f'''echo "Fitting deterministic production empty-drop profile"
rm -f "{plan['fixed_profile']}" "{plan['fixed_diagnostics']}"
"{binary}" -o "{plan['raw_prefix']}" \\
  --filtered_barcodes "{get_filtered_barcodes(lib_num)}" \\
  --interindividual --condf "{plan['raw_prefix']}.condf" \\
  --ids "{plan['production_roster']}" \\
  --max_empty {AMBIENT_VALIDATION_PROFILE_MAX_EMPTY} \\
  --seed {AMBIENT_VALIDATION_PROFILE_SEED} \\
  --profile_starts {AMBIENT_VALIDATION_PROFILE_STARTS} \\
  --bootstrap 0 --output "{plan['fixed_profile']}" \\
  --diagnostics "{plan['fixed_diagnostics']}" \\
  -T "${{SLURM_CPUS_PER_TASK}}"'''
    script = f'''#!/bin/bash
#SBATCH --job-name={job_name}
#SBATCH --output={log_dir}/{job_name}_%j.out
#SBATCH --error={log_dir}/{job_name}_%j.err
#SBATCH --chdir={AGGREGATE_ROOT}
#SBATCH --time={SLURM_TIME}
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --partition={SLURM_PARTITION}
#SBATCH --nodes=1

set -euo pipefail
{module_block(AMBIENT_PYTHON_MODULES + ("htslib/1.20", "cellbouncer/dev"))}
mkdir -p "{plan['profile_dir']}"
{_ambient_validation_plan_guard_shell(plan)}
{skip_block}

"{binary}" --help 2>&1 | grep -F -- "--profile_starts" >/dev/null || {{
    echo "ERROR: deployed tet_ambient_profile lacks --profile_starts" >&2
    exit 2
}}
"{binary}" --version 2>&1 | grep -F "v{AMBIENT_PROFILE_REQUIRED_VERSION}" >/dev/null || {{
    echo "ERROR: deployed tet_ambient_profile must be v{AMBIENT_PROFILE_REQUIRED_VERSION}" >&2
    exit 2
}}

{production_run}
awk '
    NR == FNR {{if (NF && $1 !~ /^#/) allowed[$1]=1; next}}
    NF >= 2 {{
        if (!($1 in allowed)) {{
            printf "ERROR: non-production donor %s entered fixed profile\\n", $1 > "/dev/stderr"
            bad=1
        }}
        observed[$1]=1
        mass += $2
    }}
    END {{
        for (donor in allowed) if (!(donor in observed)) {{
            printf "ERROR: production donor %s is absent from fixed profile\\n", donor > "/dev/stderr"
            bad=1
        }}
        if (!(mass > 0.999 && mass < 1.001)) {{
            printf "ERROR: fixed profile mass sum is %.12g, expected 1\\n", mass > "/dev/stderr"
            bad=1
        }}
        exit bad
    }}
' "{plan['production_roster']}" "{plan['fixed_profile']}"

test -s "{plan['fixed_profile']}"
test -s "{plan['fixed_diagnostics']}"
awk -F '\t' '$1 == "profile_starts_successful" && $2 >= {AMBIENT_VALIDATION_PROFILE_STARTS} {{ok=1}} END {{exit !ok}}' \
  "{plan['fixed_diagnostics']}"
echo "AMBIENT_VALIDATE fixed profile complete for lib{lib_num}"
'''
    _ambient_write_text(script_path, script, executable=True)
    return script_path


def _ambient_validation_output_cleanup(prefix):
    suffixes = (
        ".contam_prof", ".contam_rate", ".allele_ratio",
        ".contam_diagnostics.tsv", ".cell_source_profile.tsv",
        ".profile_fit_diagnostics.tsv", ".condf_coverage.tsv",
        ".run_contract.json", ".decontam.assignments", ".model_fit.tsv",
        ".geometry_gate_audit.tsv", ".pass1.contam_prof",
        ".pass1.contam_rate", ".pass1.decontam.assignments",
        ".fixed_ambient_validation_arm.tsv",
    )
    return " ".join(shlex.quote(prefix + suffix) for suffix in suffixes)


def ambient_validation_fixed_outputs_complete(prefix, cond):
    required = (
        ".contam_prof", ".contam_rate", ".allele_ratio",
        ".contam_diagnostics.tsv", ".cell_source_profile.tsv",
        ".profile_fit_diagnostics.tsv", ".condf_coverage.tsv",
        ".run_contract.json", ".decontam.assignments", ".model_fit.tsv",
        ".fixed_ambient_validation_arm.tsv",
    )
    if cond.get("runner") == "geometry_gated":
        required += (".geometry_gate_audit.tsv",)
    return all(check_file_exists(prefix + suffix) for suffix in required)


def generate_ambient_validation_fixed_script(
        lib_num, cond, args, plan, dep_job_id=None, force=False):
    if cond["mode"] not in (1, 3):
        raise ValueError(
            f"AMBIENT_VALIDATE supports individual modes 1/3; "
            f"{cond['abbrev']} has mode {cond['mode']}")
    job_name = f"ambV_fixed_{cond['abbrev']}_lib{lib_num}"
    log_dir = os.path.join(
        get_log_dir(), "AMBIENT_VALIDATE", cond["abbrev"], f"lib{lib_num}")
    script_path = os.path.join(get_script_dir(), f"{job_name}.sbatch")
    os.makedirs(log_dir, exist_ok=True)
    os.makedirs(get_script_dir(), exist_ok=True)
    dep_line = (
        f"#SBATCH --dependency=afterok:{dep_job_id}"
        if dep_job_id else "")
    blocks = []
    for arm_key, fixed_spec in AMBIENT_VALIDATION_FIXED_ARMS.items():
        source_arm = fixed_spec["source_arm"]
        arm_inputs = identity_ambient_arm_inputs(lib_num, source_arm)
        prefix = get_ambient_validation_fixed_prefix(
            lib_num, cond["abbrev"], arm_key,
            args.reconciliation_candidate_set)
        command = build_contam_command(
            lib_num, cond, assignment_source=source_arm,
            out_prefix_override=prefix,
            fixed_ambient=plan["fixed_profile"])
        required_suffixes = [
            ".contam_prof", ".contam_rate", ".allele_ratio",
            ".contam_diagnostics.tsv", ".cell_source_profile.tsv",
            ".profile_fit_diagnostics.tsv", ".condf_coverage.tsv",
            ".run_contract.json", ".decontam.assignments", ".model_fit.tsv",
        ]
        if cond.get("runner") == "geometry_gated":
            required_suffixes.append(".geometry_gate_audit.tsv")
        completion_tests = " && ".join(
            f'[[ -s "{prefix}{suffix}" ]]'
            for suffix in required_suffixes)
        required_output_tests = "\n".join(
            f'test -s "{prefix}{suffix}"'
            for suffix in required_suffixes)
        cleanup = (
            f'rm -f {_ambient_validation_output_cleanup(prefix)}\n'
            f'rm -f "{prefix}.geometry_base_endpoint".* '
            f'"{prefix}.geometry_fallback_endpoint".*')
        contract = prefix + ".fixed_ambient_validation_arm.tsv"
        run_body = f'''{cleanup}
ln -sfn "{plan['demux_prefix']}.counts" "{prefix}.counts"
ln -sfn "{plan['demux_prefix']}.condf" "{prefix}.condf"
ln -sfn "{plan['demux_prefix']}.samples" "{prefix}.samples"
ln -sfn "{arm_inputs['assignments']}" "{prefix}.assignments"
{command}
grep -q $'^ambient_profile_fixed\\ttrue$' "{prefix}.model_fit.tsv" || {{
    echo "ERROR: fixed ambient profile was not retained for {fixed_spec['arm']}" >&2
    exit 3
}}
fixed_profile_sha256=$(sha256sum "{plan['fixed_profile']}" | awk '{{print $1}}')
printf '%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n' \
  'library' 'condition' 'arm' 'arm_key' 'fixed_profile' 'fixed_profile_sha256' 'validation_plan_fingerprint' \
  'lib{lib_num}' '{cond['abbrev']}' '{fixed_spec['arm']}' '{arm_key}' \
  '{plan['fixed_profile']}' "$fixed_profile_sha256" '{plan['plan_fingerprint']}' > "{contract}"
{required_output_tests}
'''
        if force:
            execution = run_body
        else:
            execution = f'''if {completion_tests} && \
   [[ -s "{prefix}.fixed_ambient_validation_arm.tsv" ]] && \
   grep -q $'^ambient_profile_fixed\\ttrue$' "{prefix}.model_fit.tsv"; then
    current_fixed_sha256=$(sha256sum "{plan['fixed_profile']}" | awk '{{print $1}}')
    if awk -F '\\t' -v sha="$current_fixed_sha256" '
        NR == 1 {{for (i=1; i<=NF; i++) col[$i]=i; next}}
        NR == 2 {{
            rows++
            ok=(col["fixed_profile_sha256"] &&
                $(col["fixed_profile_sha256"]) == sha)
        }}
        NR > 2 {{rows++}}
        END {{exit !(ok && rows == 1)}}
    ' "{contract}"; then
        echo "Complete {fixed_spec['short_label']} output uses the same frozen profile; skipping"
    else
        echo "Existing {fixed_spec['short_label']} output uses a stale profile; rebuilding"
{run_body}
    fi
else
{run_body}
fi
'''
        block = f'''echo "Running {fixed_spec['label']}"
mkdir -p "{os.path.dirname(prefix)}"
{execution}
'''
        blocks.append(block)
    script = f'''#!/bin/bash
#SBATCH --job-name={job_name}
#SBATCH --output={log_dir}/{job_name}_%j.out
#SBATCH --error={log_dir}/{job_name}_%j.err
#SBATCH --chdir={AGGREGATE_ROOT}
#SBATCH --time={SLURM_TIME}
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --partition={SLURM_PARTITION}
#SBATCH --nodes=1
{dep_line}

set -euo pipefail
{module_block()}
test -s "{plan['fixed_profile']}"
{_ambient_validation_plan_guard_shell(plan)}
{chr(10).join(blocks)}
echo "AMBIENT_VALIDATE frozen-profile rescoring complete for {cond['abbrev']} lib{lib_num}"
'''
    _ambient_write_text(script_path, script, executable=True)
    return script_path


def _ambient_validation_read_two_column(path):
    result = {}
    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            parts = raw.split()
            if len(parts) < 2 or parts[0].lower() in {"barcode", "cell"}:
                continue
            try:
                value = float(parts[1])
            except ValueError:
                continue
            if math.isfinite(value):
                result[parts[0]] = value
    return result


def ambient_validation_worker(spec_path):
    """Aggregate fixed-profile E/F assignment contrasts across libraries."""
    with open(spec_path, "r", encoding="utf-8") as handle:
        spec = json.load(handle)
    import numpy as np
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    output_root = Path(spec["output_root"])
    data_dir = output_root / "data"
    figure_dir = output_root / "figures"
    data_dir.mkdir(parents=True, exist_ok=True)
    figure_dir.mkdir(parents=True, exist_ok=True)
    contrast_rows = []
    scatter_payload = []

    for library_spec in spec["libraries"]:
        library = int(library_spec["library"])
        with open(library_spec["scrutiny_cells"], "r", encoding="utf-8") as handle:
            scrutiny = {
                row["barcode"]: row["stratum"]
                for row in csv.DictReader(handle, delimiter="\t")}
        for condition, paths in library_spec["conditions"].items():
            left = _ambient_validation_read_two_column(paths["E"])
            right = _ambient_validation_read_two_column(paths["F"])
            common = sorted(set(left) & set(right))
            if not common:
                raise RuntimeError(
                    f"lib{library} {condition}: E/F have no paired cells")
            for stratum in (
                    "changed_target", "scrutinized_other", "background", "all"):
                cells = [
                    barcode for barcode in common
                    if stratum == "all" or
                    scrutiny.get(barcode, "background") == stratum]
                if not cells:
                    continue
                lv = np.asarray([left[cell] for cell in cells], dtype=float)
                rv = np.asarray([right[cell] for cell in cells], dtype=float)
                delta = rv - lv
                contrast_rows.append({
                    "library": f"lib{library}", "condition": condition,
                    "comparison": "fixed_profile_E_to_F", "stratum": stratum,
                    "n_paired": len(cells),
                    "median_demux_c": float(np.median(lv)),
                    "median_reconciled_c": float(np.median(rv)),
                    "median_delta_c": float(np.median(delta)),
                    "q25_delta_c": float(np.quantile(delta, 0.25)),
                    "q75_delta_c": float(np.quantile(delta, 0.75)),
                    "fraction_reconciled_lower": float(np.mean(delta < 0)),
                    "fraction_equal": float(np.mean(delta == 0)),
                })
                if stratum != "all":
                    scatter_payload.extend((
                        f"lib{library}", condition, stratum,
                        float(left[cell]), float(right[cell])) for cell in cells)

    def write_tsv(path, rows, fields=None):
        fields = list(fields or (list(rows[0]) if rows else ()))
        if not fields:
            raise RuntimeError(f"no TSV schema available for {path}")
        with open(path, "w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t",
                                    lineterminator="\n")
            writer.writeheader(); writer.writerows(rows)

    write_tsv(
        data_dir / "fixed_profile_assignment_contrasts.tsv", contrast_rows,
        ("library", "condition", "comparison", "stratum", "n_paired",
         "median_demux_c", "median_reconciled_c", "median_delta_c",
         "q25_delta_c", "q75_delta_c", "fraction_reconciled_lower",
         "fraction_equal"))

    applicability_rows = list(spec.get("applicability", ()))
    if applicability_rows:
        write_tsv(
            data_dir / "library_applicability.tsv", applicability_rows,
            ("library", "applicable", "targets", "reason", "fixed_profile",
             "plan_fingerprint"))

    if scatter_payload:
        fig, ax = plt.subplots(figsize=(7, 6.5))
        palette = {"changed_target": "#CC6677", "scrutinized_other": "#4477AA", "background": "#BBBBBB"}
        for stratum in ("background", "scrutinized_other", "changed_target"):
            rows = [row for row in scatter_payload if row[2] == stratum]
            if not rows:
                continue
            ax.scatter([row[3] for row in rows], [row[4] for row in rows],
                       s=7 if stratum == "background" else 16,
                       alpha=0.20 if stratum == "background" else 0.60,
                       color=palette[stratum], label=stratum.replace("_", " "))
        limits = ax.get_xlim(); upper = max(limits[1], ax.get_ylim()[1])
        ax.plot([0, upper], [0, upper], color="#333333", linewidth=0.8,
                linestyle="--")
        ax.set_xlim(0, upper); ax.set_ylim(0, upper)
        ax.set_xlabel("E: demux assignments, frozen empty-drop profile (c)")
        ax.set_ylabel("F: reconciled assignments, same frozen profile (c)")
        ax.set_title("All applicable libraries: fixed-profile assignment contrast")
        ax.legend(frameon=False)
        fig.tight_layout()
        for suffix in ("pdf", "png"):
            fig.savefig(figure_dir / f"fixed_profile_assignment_scatter.{suffix}",
                        dpi=220, bbox_inches="tight")
        plt.close(fig)

        def kde_curve(values, grid):
            values = np.asarray(values, dtype=float)
            values = values[np.isfinite(values)]
            if values.size == 0:
                return np.zeros_like(grid)
            spread = (
                float(np.std(values, ddof=1)) if values.size > 1 else 0.0)
            bandwidth = max(
                0.005, 1.06 * spread * values.size ** (-0.2))
            density = np.zeros_like(grid)
            normalizer = bandwidth * math.sqrt(2.0 * math.pi) * values.size
            for start in range(0, values.size, 2000):
                chunk = values[start:start + 2000]
                scaled = (grid[:, None] - chunk[None, :]) / bandwidth
                density += np.exp(-0.5 * scaled * scaled).sum(axis=1)
            return density / normalizer

        library_conditions = sorted({
            (row[0], row[1]) for row in scatter_payload})
        for library_name, condition in library_conditions:
            strata = [
                stratum for stratum in (
                    "changed_target", "scrutinized_other", "background")
                if any(row[0] == library_name and row[1] == condition and
                       row[2] == stratum for row in scatter_payload)]
            if not strata:
                continue
            values_all = [
                value for row in scatter_payload
                if row[0] == library_name and row[1] == condition
                for value in (row[3], row[4])]
            upper = min(1.0, max(0.10, float(np.quantile(values_all, 0.995)) * 1.10))
            grid = np.linspace(0.0, upper, 320)
            fig, axes = plt.subplots(
                len(strata), 1, figsize=(8.5, 2.55 * len(strata) + 1.3),
                sharex=True, squeeze=False)
            for ax, stratum in zip(axes.ravel(), strata):
                rows = [
                    row for row in scatter_payload
                    if row[0] == library_name and row[1] == condition and
                    row[2] == stratum]
                left = [row[3] for row in rows]
                right = [row[4] for row in rows]
                ax.plot(grid, kde_curve(left, grid), color="#4477AA",
                        linewidth=1.8, label="E demux / fixed empty")
                ax.plot(grid, kde_curve(right, grid), color="#CC6677",
                        linewidth=1.8, label="F reconciled / same profile")
                ax.set_ylabel("Density")
                ax.set_title(
                    f"{stratum.replace('_', ' ')} (n={len(rows):,})",
                    fontsize=10, loc="left")
                ax.grid(alpha=0.18)
            axes[0, 0].legend(frameon=False, fontsize=8, ncol=2)
            axes[-1, 0].set_xlabel("Estimated contamination fraction (c)")
            fig.suptitle(
                f"{library_name}: assignment-only fixed-profile KDE\n{condition}",
                fontsize=12, y=0.995)
            fig.tight_layout(rect=(0, 0, 1, 0.95))
            safe_condition = _ambient_condition_slug(condition)
            for suffix in ("pdf", "png"):
                fig.savefig(
                    figure_dir /
                    f"{library_name}_{safe_condition}_fixed_profile_kde_by_stratum.{suffix}",
                    dpi=220, bbox_inches="tight")
            plt.close(fig)

    for condition in sorted({row["condition"] for row in contrast_rows}):
        condition_rows = [
            row for row in contrast_rows
            if row["condition"] == condition and row["stratum"] != "all"]
        libraries = sorted(
            {row["library"] for row in condition_rows},
            key=lambda value: int(value.removeprefix("lib")))
        if not libraries:
            continue
        y = np.arange(len(libraries))
        fig, ax = plt.subplots(
            figsize=(9.5, max(5.0, 0.34 * len(libraries) + 2.0)))
        palette = {
            "changed_target": "#CC6677",
            "scrutinized_other": "#4477AA",
            "background": "#999999",
        }
        labels = {
            "changed_target": "changed target cells",
            "scrutinized_other": "scrutinized unchanged cells",
            "background": "background cells",
        }
        for stratum in (
                "background", "scrutinized_other", "changed_target"):
            x = []
            yy = []
            for index, library_name in enumerate(libraries):
                row = next((
                    item for item in condition_rows
                    if item["library"] == library_name and
                    item["stratum"] == stratum), None)
                if row is not None:
                    x.append(row["median_delta_c"])
                    yy.append(index)
            if x:
                ax.scatter(x, yy, s=38, color=palette[stratum],
                           label=labels[stratum], zorder=3)
        ax.axvline(0, color="#333333", linewidth=0.9, linestyle="--")
        ax.set_yticks(y)
        ax.set_yticklabels(libraries)
        ax.invert_yaxis()
        ax.set_xlabel("Median paired Δc (reconciled − demux)")
        ax.set_title(f"Fixed-profile assignment effect across libraries\n{condition}")
        ax.grid(axis="x", alpha=0.18)
        ax.legend(frameon=False, fontsize=8)
        fig.tight_layout()
        safe_condition = _ambient_condition_slug(condition)
        for suffix in ("pdf", "png"):
            fig.savefig(
                figure_dir /
                f"all_libraries_{safe_condition}_fixed_profile_delta.{suffix}",
                dpi=220, bbox_inches="tight")
        plt.close(fig)

    summary = {
        "schema_version": 1,
        "plan_version": AMBIENT_VALIDATION_PLAN_VERSION,
        "applicable_libraries": sum(
            row.get("applicable") == "yes"
            for row in spec.get("applicability", ())),
        "not_applicable_libraries": sum(
            row.get("applicable") != "yes"
            for row in spec.get("applicability", ())),
        "contrast_rows": len(contrast_rows),
        "interpretation": (
            "E-to-F holds one production-roster empty-drop profile fixed. "
            "A lower F contamination fraction therefore supports the reconciled "
            "assignment without confounding from profile refitting. Changed target "
            "cells are the primary evidence; scrutinized-other and background cells "
            "are negative-control strata."),
    }
    _ambient_write_json(
        output_root / "fixed_profile_validation_summary.json", summary)
    print(f"AMBIENT_VALIDATE wrote results under {output_root}")


def generate_ambient_validation_aggregate_script(
        args, cond_list, plans, fixed_job_ids=None):
    output_root = os.path.join(
        get_ambient_plot_run_dir(args, cond_list), "fixed_profile_validation")
    libraries = []
    applicability = []
    for lib_num, plan in sorted(plans.items()):
        applicable = bool(plan.get("applicable"))
        applicability.append({
            "library": f"lib{lib_num}",
            "applicable": "yes" if applicable else "no",
            "targets": ";".join(plan.get("targets_list", ())) or "NA",
            "reason": plan.get("skip_reason", "") or "applied donor addition",
            "fixed_profile": plan.get("fixed_profile", "") or "NA",
            "plan_fingerprint": plan.get("plan_fingerprint", "") or "NA",
        })
        if not applicable:
            continue
        conditions = {}
        for cond in cond_list:
            condition = cond["abbrev"]
            conditions[condition] = {
                "E": get_ambient_validation_fixed_prefix(
                    lib_num, condition, "demux_fixed_empty") + ".contam_rate",
                "F": get_ambient_validation_fixed_prefix(
                    lib_num, condition, "reconciled_fixed_empty") + ".contam_rate",
            }
        libraries.append({
            "library": lib_num,
            "targets": plan["targets_list"],
            "scrutiny_cells": plan["scrutiny_cells"],
            "fixed_profile": plan["fixed_profile"],
            "plan_fingerprint": plan["plan_fingerprint"],
            "conditions": conditions,
        })
    spec = {
        "worker_kind": "ambient_validation",
        "output_root": output_root,
        "libraries": libraries,
        "applicability": applicability,
    }
    slug = hashlib.sha256(json.dumps(
        spec, sort_keys=True).encode("utf-8")).hexdigest()[:10]
    spec_path = os.path.join(
        get_script_dir(), f"ambient_validate_{slug}.worker.json")
    _ambient_write_json(spec_path, spec)
    job_name = f"ambV_aggregate_{slug}"
    log_dir = os.path.join(get_log_dir(), "AMBIENT_VALIDATE")
    script_path = os.path.join(get_script_dir(), f"{job_name}.sbatch")
    dep_line = _ambient_dependency_line(fixed_job_ids)
    script = f'''#!/bin/bash
#SBATCH --job-name={job_name}
#SBATCH --output={log_dir}/{job_name}_%j.out
#SBATCH --error={log_dir}/{job_name}_%j.err
#SBATCH --chdir={AGGREGATE_ROOT}
#SBATCH --time={SLURM_TIME}
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --partition={SLURM_PARTITION}
#SBATCH --nodes=1
{dep_line}

set -euo pipefail
module purge
{chr(10).join(f'module load {module}' for module in AMBIENT_PYTHON_MODULES)}
echo "Runtime host: $(hostname)"
echo "SLURM job: ${{SLURM_JOB_ID:-not-set}}"
free -h || true
export MPLCONFIGDIR="${{SLURM_TMPDIR:-/tmp}}/matplotlib-${{SLURM_JOB_ID:-$$}}"
mkdir -p "$MPLCONFIGDIR"
python3 "{os.path.abspath(__file__)}" --_ambient-plots-worker "{spec_path}"
test -s "{os.path.join(output_root, 'fixed_profile_validation_summary.json')}"
test -d "{os.path.join(output_root, 'figures')}"
python3 "{os.path.abspath(__file__)}" --_publish-figure-shortcut \
  ambient_validation "{os.path.join(output_root, 'figures')}"
'''
    os.makedirs(log_dir, exist_ok=True)
    _ambient_write_text(script_path, script, executable=True)
    return {
        "script": script_path, "spec": spec_path,
        "output_root": output_root,
    }


def generate_ambient_swap_profile_script(lib_num, plan, force=False):
    """Generate one candidate-roster empty-drop profile job."""
    binary = os.path.join(SOFTWARE_BIN, "tet_ambient_profile")
    job_name = f"ambS_profile_lib{lib_num}"
    log_dir = os.path.join(
        get_log_dir(), "AMBIENT_SWAP_TEST", f"lib{lib_num}")
    script_path = os.path.join(get_script_dir(), f"{job_name}.sbatch")
    os.makedirs(log_dir, exist_ok=True)
    os.makedirs(get_script_dir(), exist_ok=True)
    os.makedirs(plan["profile_dir"], exist_ok=True)
    skip_block = ""
    if not force and ambient_validation_profile_outputs_complete(plan):
        skip_block = (
            'echo "Matching candidate-roster profile already exists; skipping"\n'
            "exit 0")
    script = f'''#!/bin/bash
#SBATCH --job-name={job_name}
#SBATCH --output={log_dir}/{job_name}_%j.out
#SBATCH --error={log_dir}/{job_name}_%j.err
#SBATCH --chdir={AGGREGATE_ROOT}
#SBATCH --time={SLURM_TIME}
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --partition={SLURM_PARTITION}
#SBATCH --nodes=1

set -euo pipefail
{module_block(AMBIENT_PYTHON_MODULES + ("htslib/1.20", "cellbouncer/dev"))}
mkdir -p "{plan['profile_dir']}"
{_ambient_validation_plan_guard_shell(plan)}
{skip_block}

"{binary}" --help 2>&1 | grep -F -- "--profile_starts" >/dev/null || {{
    echo "ERROR: deployed tet_ambient_profile lacks --profile_starts" >&2
    exit 2
}}
"{binary}" --version 2>&1 | grep -F "v{AMBIENT_PROFILE_REQUIRED_VERSION}" >/dev/null || {{
    echo "ERROR: deployed tet_ambient_profile must be v{AMBIENT_PROFILE_REQUIRED_VERSION}" >&2
    exit 2
}}

echo "Fitting deterministic candidate-roster empty-drop profile"
rm -f "{plan['fixed_profile']}" "{plan['fixed_diagnostics']}"
"{binary}" -o "{plan['raw_prefix']}" \
  --filtered_barcodes "{get_filtered_barcodes(lib_num)}" \
  --interindividual --condf "{plan['raw_prefix']}.condf" \
  --ids "{plan['production_roster']}" \
  --max_empty {AMBIENT_VALIDATION_PROFILE_MAX_EMPTY} \
  --seed {AMBIENT_VALIDATION_PROFILE_SEED} \
  --profile_starts {AMBIENT_VALIDATION_PROFILE_STARTS} \
  --bootstrap 0 --output "{plan['fixed_profile']}" \
  --diagnostics "{plan['fixed_diagnostics']}" \
  -T "${{SLURM_CPUS_PER_TASK}}"

awk '
    NR == FNR {{if (NF && $1 !~ /^#/) allowed[$1]=1; next}}
    NF >= 2 {{
        if (!($1 in allowed)) {{
            printf "ERROR: non-roster donor %s entered fixed profile\\n", $1 > "/dev/stderr"
            bad=1
        }}
        observed[$1]=1
        mass += $2
    }}
    END {{
        for (donor in allowed) if (!(donor in observed)) {{
            printf "ERROR: roster donor %s is absent from fixed profile\\n", donor > "/dev/stderr"
            bad=1
        }}
        if (!(mass > 0.999 && mass < 1.001)) {{
            printf "ERROR: fixed profile mass sum is %.12g, expected 1\\n", mass > "/dev/stderr"
            bad=1
        }}
        exit bad
    }}
' "{plan['production_roster']}" "{plan['fixed_profile']}"

test -s "{plan['fixed_profile']}"
test -s "{plan['fixed_diagnostics']}"
awk -F '\t' '$1 == "profile_starts_successful" && $2 >= {AMBIENT_VALIDATION_PROFILE_STARTS} {{ok=1}} END {{exit !ok}}' \
  "{plan['fixed_diagnostics']}"
echo "AMBIENT_SWAP_TEST fixed profile complete for lib{lib_num}"
'''
    _ambient_write_text(script_path, script, executable=True)
    return script_path


def ambient_swap_arm_outputs_complete(lib_num, cond, plan):
    """Return whether all G/H/J/K outputs match the current frozen plan."""
    if not check_file_exists(plan["fixed_profile"]):
        return False
    digest = hashlib.sha256()
    try:
        with open(plan["fixed_profile"], "rb") as handle:
            for block in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(block)
        fixed_profile_sha256 = digest.hexdigest()
        for arm_key, arm_spec in AMBIENT_SWAP_TEST_ARMS.items():
            prefix = get_ambient_swap_test_arm_prefix(
                lib_num, cond["abbrev"], arm_key, plan["discovery_id"])
            required_suffixes = [
                ".contam_prof", ".contam_rate", ".allele_ratio",
                ".contam_diagnostics.tsv", ".cell_source_profile.tsv",
                ".profile_fit_diagnostics.tsv", ".condf_coverage.tsv",
                ".run_contract.json", ".decontam.assignments",
                ".model_fit.tsv", ".ambient_swap_test_arm.tsv",
            ]
            if cond.get("runner") == "geometry_gated":
                required_suffixes.append(".geometry_gate_audit.tsv")
            if not all(check_file_exists(prefix + suffix)
                       for suffix in required_suffixes):
                return False

            expected_fixed = (
                "true" if arm_spec["profile_mode"] == "fixed" else "false")
            with open(
                    prefix + ".model_fit.tsv", "r",
                    encoding="utf-8", errors="replace") as handle:
                if not any(
                        line.rstrip("\n") ==
                        f"ambient_profile_fixed\t{expected_fixed}"
                        for line in handle):
                    return False

            with open(
                    prefix + ".ambient_swap_test_arm.tsv", "r",
                    encoding="utf-8", newline="") as handle:
                reader = csv.DictReader(handle, delimiter="\t")
                fields = set(reader.fieldnames or ())
                rows = list(reader)
            if len(rows) != 1:
                return False
            row = rows[0]
            profile_sha_field = (
                "reference_fixed_profile_sha256"
                if "reference_fixed_profile_sha256" in fields
                else "fixed_profile_sha256")
            if (
                    row.get("library") != f"lib{lib_num}" or
                    row.get("condition") != cond["abbrev"] or
                    row.get("arm") != arm_spec["arm"] or
                    row.get("arm_key") != arm_key or
                    row.get("profile_mode") != arm_spec["profile_mode"] or
                    row.get(profile_sha_field) != fixed_profile_sha256 or
                    row.get("swap_test_plan_fingerprint") !=
                    plan["plan_fingerprint"]):
                return False
    except (OSError, ValueError, csv.Error):
        return False
    return True


def generate_ambient_swap_arm_script(
        lib_num, cond, plan, dep_job_id=None, force=False):
    if cond["mode"] not in (1, 3):
        raise ValueError(
            f"AMBIENT_SWAP_TEST supports individual modes 1/3; "
            f"{cond['abbrev']} has mode {cond['mode']}")
    job_name = f"ambS_arms_{cond['abbrev']}_lib{lib_num}"
    log_dir = os.path.join(
        get_log_dir(), "AMBIENT_SWAP_TEST", cond["abbrev"],
        f"lib{lib_num}")
    script_path = os.path.join(get_script_dir(), f"{job_name}.sbatch")
    os.makedirs(log_dir, exist_ok=True)
    dep_line = (
        f"#SBATCH --dependency=afterok:{dep_job_id}"
        if dep_job_id else "")
    blocks = []
    input_overrides = {
        "receiver_lines": plan["receiver_lines"],
        "ambient_candidates": plan["production_roster"],
    }
    for arm_key, arm_spec in AMBIENT_SWAP_TEST_ARMS.items():
        prefix = get_ambient_swap_test_arm_prefix(
            lib_num, cond["abbrev"], arm_key, plan["discovery_id"])
        assignments = plan[arm_spec["assignment_key"]]
        profile_mode = arm_spec["profile_mode"]
        fixed_ambient = (
            plan["fixed_profile"] if profile_mode == "fixed" else None)
        expected_fixed_flag = "true" if profile_mode == "fixed" else "false"
        command = build_contam_command(
            lib_num, cond, assignment_source="demux",
            out_prefix_override=prefix,
            fixed_ambient=fixed_ambient,
            input_overrides=input_overrides)
        required_suffixes = [
            ".contam_prof", ".contam_rate", ".allele_ratio",
            ".contam_diagnostics.tsv", ".cell_source_profile.tsv",
            ".profile_fit_diagnostics.tsv", ".condf_coverage.tsv",
            ".run_contract.json", ".decontam.assignments",
            ".model_fit.tsv",
        ]
        if cond.get("runner") == "geometry_gated":
            required_suffixes.append(".geometry_gate_audit.tsv")
        completion_tests = " && ".join(
            f'[[ -s "{prefix}{suffix}" ]]'
            for suffix in required_suffixes)
        required_output_tests = "\n".join(
            f'test -s "{prefix}{suffix}"'
            for suffix in required_suffixes)
        cleanup = (
            f'rm -f {_ambient_validation_output_cleanup(prefix)} '
            f'"{prefix}.ambient_swap_test_arm.tsv"\n'
            f'rm -f "{prefix}.geometry_base_endpoint".* '
            f'"{prefix}.geometry_fallback_endpoint".*')
        contract = prefix + ".ambient_swap_test_arm.tsv"
        run_body = f'''{cleanup}
ln -sfn "{plan['demux_prefix']}.counts" "{prefix}.counts"
ln -sfn "{plan['demux_prefix']}.condf" "{prefix}.condf"
ln -sfn "{plan['demux_prefix']}.samples" "{prefix}.samples"
ln -sfn "{assignments}" "{prefix}.assignments"
{command}
grep -q $'^ambient_profile_fixed\\t{expected_fixed_flag}$' "{prefix}.model_fit.tsv" || {{
    echo "ERROR: unexpected ambient-profile mode for {arm_spec['arm']} ({profile_mode})" >&2
    exit 3
}}
fixed_profile_sha256=$(sha256sum "{plan['fixed_profile']}" | awk '{{print $1}}')
printf '%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n' \
  'library' 'condition' 'arm' 'arm_key' 'profile_mode' 'reference_fixed_profile' 'reference_fixed_profile_sha256' 'swap_test_plan_fingerprint' \
  'lib{lib_num}' '{cond['abbrev']}' '{arm_spec['arm']}' '{arm_key}' '{profile_mode}' \
  '{plan['fixed_profile']}' "$fixed_profile_sha256" '{plan['plan_fingerprint']}' > "{contract}"
{required_output_tests}
'''
        if force:
            execution = run_body
        else:
            execution = f'''if {completion_tests} && \
   [[ -s "{contract}" ]] && \
   grep -q $'^ambient_profile_fixed\\t{expected_fixed_flag}$' "{prefix}.model_fit.tsv"; then
    current_fixed_sha256=$(sha256sum "{plan['fixed_profile']}" | awk '{{print $1}}')
    if awk -F '\t' -v sha="$current_fixed_sha256" \
        -v mode="{profile_mode}" -v fingerprint="{plan['plan_fingerprint']}" '
        NR == 1 {{for (i=1; i<=NF; i++) col[$i]=i; next}}
        NR == 2 {{
            rows++
            profile_col=(col["reference_fixed_profile_sha256"] ?
                         col["reference_fixed_profile_sha256"] :
                         col["fixed_profile_sha256"])
            mode_ok=(col["profile_mode"] ?
                     $(col["profile_mode"]) == mode : mode == "fixed")
            ok=(mode_ok && profile_col && $(profile_col) == sha &&
                col["swap_test_plan_fingerprint"] &&
                $(col["swap_test_plan_fingerprint"]) == fingerprint)
        }}
        NR > 2 {{rows++}}
        END {{exit !(ok && rows == 1)}}
    ' "{contract}"; then
        echo "Complete {arm_spec['short_label']} output uses the same plan; skipping"
    else
        echo "Existing {arm_spec['short_label']} output uses a stale plan; rebuilding"
{run_body}
    fi
else
{run_body}
fi
'''
        blocks.append(f'''echo "Running {arm_spec['label']}"
mkdir -p "{os.path.dirname(prefix)}"
{execution}
''')
    script = f'''#!/bin/bash
#SBATCH --job-name={job_name}
#SBATCH --output={log_dir}/{job_name}_%j.out
#SBATCH --error={log_dir}/{job_name}_%j.err
#SBATCH --chdir={AGGREGATE_ROOT}
#SBATCH --time={SLURM_TIME}
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --partition={SLURM_PARTITION}
#SBATCH --nodes=1
{dep_line}

set -euo pipefail
{module_block()}
test -s "{plan['fixed_profile']}"
{_ambient_validation_plan_guard_shell(plan)}
{chr(10).join(blocks)}
echo "AMBIENT_SWAP_TEST rescoring complete for {cond['abbrev']} lib{lib_num}"
'''
    _ambient_write_text(script_path, script, executable=True)
    return script_path



def generate_ambient_swap_aggregate_script(
        args, cond_list, plans, discovery, arm_job_ids=None):
    discovery_id = ambient_swap_discovery_id(discovery)
    run_root = get_ambient_swap_test_output_root(
        args, cond_list, discovery_id)
    output_root = os.path.join(run_root, "identity_reconciliation_figures")
    libraries = []
    applicability = []
    for lib_num, plan in sorted(plans.items()):
        applicable = bool(plan.get("applicable"))
        applicability.append({
            "library": f"lib{lib_num}",
            "applicable": "yes" if applicable else "no",
            "candidates": ";".join(
                plan.get("candidate_identities", ())) or "NA",
            "reason": plan.get("skip_reason", "") or "supported events discovered",
            "plan_fingerprint": plan.get("plan_fingerprint", "") or "NA",
        })
        if not applicable:
            continue
        conditions = {}
        for cond in cond_list:
            condition = cond["abbrev"]
            conditions[condition] = {
                arm_spec["arm"]: get_ambient_swap_test_arm_prefix(
                    lib_num, condition, arm_key, discovery_id)
                    + ".contam_rate"
                for arm_key, arm_spec in AMBIENT_SWAP_TEST_ARMS.items()
            }
        libraries.append({
            "library": lib_num,
            "candidate_cells": plan["candidate_cells"],
            "cell_strata": plan["cell_strata"],
            "event_discovery": plan["event_discovery"],
            "fixed_profile": plan["fixed_profile"],
            "plan_fingerprint": plan["plan_fingerprint"],
            "conditions": conditions,
        })
    spec = {
        "schema_version": "ambient_swap_figure_input_V2",
        "output_root": output_root,
        "run_root": run_root,
        "condition_short_names": {
            condition["abbrev"]: AMBIENT_CONDITION_SHORT_NAMES.get(
                condition["abbrev"], condition["abbrev"])
            for condition in cond_list
        },
        "identity_event_path": discovery["event_path"],
        "reconciliation_manifest": discovery["manifest_path"],
        "discovery_id": discovery_id,
        "libraries": libraries,
        "applicability": applicability,
    }
    slug = hashlib.sha256(json.dumps(
        spec, sort_keys=True).encode("utf-8")).hexdigest()[:10]
    spec_path = os.path.join(
        get_script_dir(), f"ambient_swap_test_{slug}.worker.json")
    _ambient_write_json(spec_path, spec)
    job_name = f"ambS_aggregate_{slug}"
    log_dir = os.path.join(get_log_dir(), "AMBIENT_SWAP_TEST")
    script_path = os.path.join(get_script_dir(), f"{job_name}.sbatch")
    dep_line = _ambient_dependency_line(arm_job_ids)
    script = f'''#!/bin/bash
#SBATCH --job-name={job_name}
#SBATCH --output={log_dir}/{job_name}_%j.out
#SBATCH --error={log_dir}/{job_name}_%j.err
#SBATCH --chdir={AGGREGATE_ROOT}
#SBATCH --time={SLURM_TIME}
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --partition={SLURM_PARTITION}
#SBATCH --nodes=1
{dep_line}

set -euo pipefail
module purge
{chr(10).join(f'module load {module}' for module in AMBIENT_PYTHON_MODULES)}
echo "Runtime host: $(hostname)"
echo "SLURM job: ${{SLURM_JOB_ID:-not-set}}"
free -h || true
export MPLCONFIGDIR="${{SLURM_TMPDIR:-/tmp}}/matplotlib-${{SLURM_JOB_ID:-$$}}"
mkdir -p "$MPLCONFIGDIR"
STAGING_ROOT="{output_root}.staging.${{SLURM_JOB_ID:-$$}}"
python3 - "{run_root}" "{output_root}" "$STAGING_ROOT" <<'PY'
import shutil
import sys
from pathlib import Path

run_root = Path(sys.argv[1]).resolve()
output_root = Path(sys.argv[2])
staging_root = Path(sys.argv[3])
for candidate in (output_root, staging_root):
    if candidate.parent.resolve() != run_root:
        raise SystemExit("ERROR: swap-test figure path escaped the current run root")
    if candidate.is_symlink():
        raise SystemExit("ERROR: refusing a symlinked swap-test figure path")
if staging_root.exists():
    shutil.rmtree(staging_root)
staging_root.mkdir(parents=True)
PY
python3 "{IDENTITY_RECONCILIATION_FIGURES}" \
  --ambient-swap-spec "{spec_path}" \
  --output-dir "$STAGING_ROOT"
test -s "$STAGING_ROOT/ambient_swap_figure_summary.json"
test -s "$STAGING_ROOT/data/ambient_swap_figure_manifest.tsv"
python3 - "{run_root}" "{output_root}" "$STAGING_ROOT" <<'PY'
import os
import shutil
import sys
from pathlib import Path

run_root = Path(sys.argv[1]).resolve()
output_root = Path(sys.argv[2])
staging_root = Path(sys.argv[3])
for candidate in (output_root, staging_root):
    if candidate.parent.resolve() != run_root or candidate.is_symlink():
        raise SystemExit("ERROR: unsafe swap-test figure promotion path")
backup_root = output_root.with_name(
    output_root.name + ".previous." + os.environ.get("SLURM_JOB_ID", str(os.getpid())))
if backup_root.exists():
    if backup_root.is_symlink():
        raise SystemExit("ERROR: refusing a symlinked swap-test figure backup")
    shutil.rmtree(backup_root)
had_previous = output_root.exists()
if had_previous:
    os.replace(output_root, backup_root)
try:
    os.replace(staging_root, output_root)
except Exception:
    if had_previous and backup_root.exists() and not output_root.exists():
        os.replace(backup_root, output_root)
    raise
if backup_root.exists():
    shutil.rmtree(backup_root)
PY
test -s "{output_root}/ambient_swap_figure_summary.json"
test -s "{output_root}/data/ambient_swap_figure_manifest.tsv"
python3 "{os.path.abspath(__file__)}" --_publish-figure-shortcut \
  ambient_swap_test "{output_root}"
'''
    os.makedirs(log_dir, exist_ok=True)
    _ambient_write_text(script_path, script, executable=True)
    return {
        "script": script_path,
        "spec": spec_path,
        "output_root": output_root,
    }


# =============================================================================
# CLEANUP_RESULTS: retain newest completed products across all stage families
# =============================================================================

def _cleanup_nonempty_file(path):
    try:
        return path.is_file() and path.stat().st_size > 0
    except OSError:
        return False


def _cleanup_profile_rank(generation):
    """Return (completion rank, newest marker mtime) for a swap generation."""
    diagnostics = sorted(
        (generation / "fixed_profile").glob("*.diagnostics.tsv"))
    valid_profile = False
    marker_mtime = generation.stat().st_mtime_ns
    for path in diagnostics:
        if not _cleanup_nonempty_file(path):
            continue
        metrics = {}
        try:
            with path.open("r", encoding="utf-8", errors="replace") as handle:
                for row in csv.reader(handle, delimiter="\t"):
                    if len(row) >= 2 and row[0] != "metric":
                        metrics[row[0]] = row[1]
            valid_profile = (
                int(float(metrics.get("profile_starts_successful", "0"))) >=
                AMBIENT_VALIDATION_PROFILE_STARTS and
                metrics.get("bulk_log_likelihood_valid", "").lower() ==
                "true")
        except (OSError, TypeError, ValueError):
            valid_profile = False
        if valid_profile:
            marker_mtime = max(marker_mtime, path.stat().st_mtime_ns)
            break
    if not valid_profile:
        return 0, marker_mtime

    arm_names = set()
    for contract in generation.glob(
            "conditions/*/*/*.ambient_swap_test_arm.tsv"):
        if not _cleanup_nonempty_file(contract):
            continue
        arm_names.add(contract.parent.name)
        marker_mtime = max(marker_mtime, contract.stat().st_mtime_ns)
    required = set(AMBIENT_SWAP_TEST_ARMS)
    return (2 if required.issubset(arm_names) else 1), marker_mtime


def cleanup_results_worker(spec_path):
    """Remove obsolete generated results while preserving canonical state.

    The cleanup target set is deliberately constructed from known result
    contracts.  It never recursively treats the aggregate root as disposable.
    Every directory deletion is first atomically renamed inside its existing
    parent, so an interrupted cleanup cannot leave a partially visible result.
    """
    import shutil

    with open(spec_path, "r", encoding="utf-8") as handle:
        spec = json.load(handle)

    aggregate_root = Path(spec["aggregate_root"]).resolve()
    mapping_root = Path(spec["mapping_root"]).resolve()
    ambient_root = Path(spec["ambient_root"]).resolve()
    figure_root = Path(spec["figure_root"])
    log_root = Path(spec["log_root"]).resolve()
    script_root = Path(spec["script_root"]).resolve()
    libraries = sorted({int(value) for value in spec["libraries"]})
    job_id = str(os.environ.get("SLURM_JOB_ID", os.getpid()))
    records = []
    kept = []

    if aggregate_root == mapping_root or aggregate_root.parent != mapping_root:
        raise RuntimeError(
            f"unsafe aggregate cleanup root: {aggregate_root}")
    expected_names = {
        "ambient_rna", "gex_ambient", "condition_index", "condf", "figures", "hybrid",
        "identity_reconciliation", "logs", "mitochondrial", "ploidy",
        "posthoc", "slurm_scripts", "tetra_refine", "vcf_daemon_state",
    }
    if aggregate_root.name != "aggregate_library_analysis":
        raise RuntimeError(
            f"cleanup root has unexpected basename: {aggregate_root}")

    # Refuse to race any other active job from this orchestrator.  Unrelated
    # mapping workflows are allowed to continue.
    relevant_job = re.compile(
        r"^(?:amb|gex|qc_|condf|demux|empty|identity|mt_|tetra|posthoc|hybrid|"
        r"ploidy|ucnn|vcfD)", re.IGNORECASE)
    try:
        result = subprocess.run(
            ["squeue", "-h", "-u", os.environ.get("USER", ""),
             "-o", "%A|%j|%T"],
            check=False, capture_output=True, text=True)
    except OSError as exc:
        raise RuntimeError(f"cannot verify active SLURM jobs: {exc}") from exc
    if result.returncode != 0:
        raise RuntimeError(
            "cannot verify active SLURM jobs: "
            + (result.stderr.strip() or f"squeue exit {result.returncode}"))
    conflicts = []
    for line in result.stdout.splitlines():
        fields = line.split("|", 2)
        if len(fields) != 3:
            continue
        active_id, active_name, active_state = fields
        if active_id == job_id:
            continue
        if relevant_job.match(active_name):
            conflicts.append(
                f"{active_id} {active_name} {active_state}")
    if conflicts:
        raise RuntimeError(
            "refusing result cleanup while relevant jobs are active: "
            + "; ".join(conflicts[:20]))

    def is_within(path, root):
        try:
            path.resolve().relative_to(root.resolve())
            return True
        except (OSError, ValueError):
            return False

    def remove_file(path, scope, reason):
        path = Path(path)
        if path.is_symlink() or not path.is_file():
            return
        # Per-library fixed-profile files live under mapping_root; callers only
        # pass explicit current-plan sibling candidates from fixed_profile/.
        allowed = (ambient_root, log_root, script_root, mapping_root)
        if not any(is_within(path, root) for root in allowed):
            raise RuntimeError(f"file cleanup escaped allowlist: {path}")
        path.unlink()
        records.append((scope, "removed_file", str(path), reason))

    def remove_tree(path, allowed_root, scope, reason):
        path = Path(path)
        allowed_root = Path(allowed_root).resolve()
        if path.is_symlink() or not path.is_dir():
            return
        if path.resolve() == allowed_root or not is_within(path, allowed_root):
            raise RuntimeError(f"directory cleanup escaped allowlist: {path}")
        tombstone = path.with_name(
            f".cleanup_delete_{job_id}_{path.name}")
        suffix = 0
        while tombstone.exists() or tombstone.is_symlink():
            suffix += 1
            tombstone = path.with_name(
                f".cleanup_delete_{job_id}_{suffix}_{path.name}")
        os.replace(path, tombstone)
        shutil.rmtree(tombstone)
        records.append((scope, "removed_tree", str(path), reason))

    def newest(items):
        return max(items, key=lambda item: (item[0], str(item[1]))) if items else None

    # ---- Aggregate figure/result generations ----
    plot_candidates = []
    validation_candidates = []
    swap_candidates = []
    top_children = []
    if ambient_root.is_dir() and not ambient_root.is_symlink():
        top_children = sorted(
            child for child in ambient_root.iterdir()
            if child.is_dir() and not child.is_symlink())
        for top in top_children:
            plot_marker = top / "data" / "plot_manifest.tsv"
            if _cleanup_nonempty_file(plot_marker):
                plot_candidates.append((plot_marker.stat().st_mtime_ns, top))
            validation_marker = (
                top / "fixed_profile_validation" /
                "fixed_profile_validation_summary.json")
            if _cleanup_nonempty_file(validation_marker):
                validation_candidates.append(
                    (validation_marker.stat().st_mtime_ns, top))
            for generation in sorted(top.glob("auto_*")):
                if generation.is_symlink() or not generation.is_dir():
                    continue
                summary = (generation / "identity_reconciliation_figures" /
                           "ambient_swap_figure_summary.json")
                manifest = (generation / "identity_reconciliation_figures" /
                            "data" / "ambient_swap_figure_manifest.tsv")
                if (_cleanup_nonempty_file(summary) and
                        _cleanup_nonempty_file(manifest)):
                    swap_candidates.append((
                        max(summary.stat().st_mtime_ns,
                            manifest.stat().st_mtime_ns),
                        generation, top))

    latest_plot = newest(plot_candidates)
    latest_validation = newest(validation_candidates)
    latest_swap = max(
        swap_candidates,
        key=lambda item: (item[0], str(item[1])),
        default=None)
    keep_tops = set()
    if latest_plot:
        keep_tops.add(latest_plot[1].resolve())
        kept.append(("AMBIENT_PLOTS", str(latest_plot[1])))
    if latest_validation:
        keep_tops.add(latest_validation[1].resolve())
        kept.append(("AMBIENT_VALIDATE", str(latest_validation[1])))
    keep_swap_generation = None
    if latest_swap:
        keep_swap_generation = latest_swap[1].resolve()
        keep_tops.add(latest_swap[2].resolve())
        kept.append(("AMBIENT_SWAP_TEST", str(latest_swap[1])))

    recognized_names = tuple(
        sorted(set(spec.get("condition_names", ())), key=len, reverse=True))

    def recognized_ambient_root(top):
        return (
            _cleanup_nonempty_file(
                top / "data" / "ambient_plot_worker_spec.json") or
            (top / "fixed_profile_validation").exists() or
            (top / "contam_r").exists() or
            (top / "kde").exists() or
            "ambient-swap-test" in top.name or
            top.name.startswith("selected_") or
            any(top.name.startswith(name) for name in recognized_names)
        )

    if top_children and not (latest_plot or latest_validation or latest_swap):
        raise RuntimeError(
            "ambient_rna contains directories but no completed result marker; "
            "refusing to delete any aggregate result")

    for top in top_children:
        if not recognized_ambient_root(top):
            kept.append(("UNRECOGNIZED", str(top)))
            continue
        if top.resolve() not in keep_tops:
            remove_tree(
                top, ambient_root, "aggregate_results",
                "older aggregate result generation")
            continue
        if keep_swap_generation is not None:
            for generation in sorted(top.glob("auto_*")):
                if (generation.is_dir() and not generation.is_symlink() and
                        generation.resolve() != keep_swap_generation):
                    remove_tree(
                        generation, ambient_root, "aggregate_results",
                        "older AMBIENT_SWAP_TEST generation")

    # Remove stale PDF/PNG files from the retained AMBIENT_PLOTS tree.  The
    # current Python and contam.R manifests define the complete allowed set.
    if latest_plot and latest_plot[1].is_dir():
        current_plot = latest_plot[1]
        allowed_figures = set()
        manifest = current_plot / "data" / "plot_manifest.tsv"
        try:
            with manifest.open("r", encoding="utf-8", newline="") as handle:
                for row in csv.DictReader(handle, delimiter="\t"):
                    value = str(row.get("path", "")).strip()
                    if not value:
                        continue
                    path = Path(value)
                    if not path.is_absolute():
                        path = current_plot / path
                    if is_within(path, current_plot):
                        allowed_figures.add(path.resolve())
        except OSError:
            pass
        task_manifest = current_plot / "data" / "contam_r_tasks.tsv"
        if _cleanup_nonempty_file(task_manifest):
            with task_manifest.open(
                    "r", encoding="utf-8", newline="") as handle:
                for row in csv.DictReader(handle, delimiter="\t"):
                    prefix = str(row.get("output_prefix", "")).strip()
                    if not prefix:
                        continue
                    for suffix in (".contam.pdf", ".contam.png"):
                        path = Path(prefix + suffix)
                        if is_within(path, current_plot):
                            allowed_figures.add(path.resolve())
        for root_text, dirnames, filenames in os.walk(current_plot):
            dirnames[:] = [
                name for name in dirnames
                if not Path(root_text, name).is_symlink()]
            for filename in filenames:
                if not filename.lower().endswith((".pdf", ".png")):
                    continue
                path = Path(root_text, filename)
                if path.resolve() not in allowed_figures:
                    remove_file(
                        path, "aggregate_figures",
                        "not present in newest plot manifests")

    # Abandoned staging/backup directories are generated only by figure jobs.
    if ambient_root.is_dir():
        stale_dirs = []
        for root_text, dirnames, _ in os.walk(ambient_root):
            for name in list(dirnames):
                if (name == ".staging" or ".staging." in name or
                        ".previous." in name or
                        name.startswith(".cleanup_delete_")):
                    path = Path(root_text, name)
                    if not path.is_symlink():
                        stale_dirs.append(path)
        for path in sorted(stale_dirs, key=lambda item: len(item.parts), reverse=True):
            if path.exists():
                remove_tree(
                    path, ambient_root, "aggregate_staging",
                    "abandoned staging or backup directory")

    # Refresh the short, stable central figure index after obsolete aggregate
    # generations have been removed.  These are relative symlinks, not copies.
    if figure_root.parent.resolve() != aggregate_root:
        raise RuntimeError(f"unsafe central figure root: {figure_root}")
    figure_root.mkdir(parents=True, exist_ok=True)
    if figure_root.is_symlink() or not figure_root.is_dir():
        raise RuntimeError(f"unsafe central figure root: {figure_root}")
    figure_targets = {
        "ambient_plots": latest_plot[1] if latest_plot else None,
        "ambient_validation": (
            latest_validation[1] / "fixed_profile_validation" / "figures"
            if latest_validation else None),
        "ambient_swap_test": (
            latest_swap[1] / "identity_reconciliation_figures"
            if latest_swap else None),
    }
    for name, target in figure_targets.items():
        link = figure_root / name
        if target is None or not target.is_dir():
            if link.is_symlink():
                link.unlink()
                records.append((
                    "figure_index", "removed_link", str(link),
                    "no completed current figure target"))
            continue
        target = target.resolve()
        if not is_within(target, aggregate_root):
            raise RuntimeError(f"figure target escaped aggregate root: {target}")
        if link.exists() and not link.is_symlink():
            raise RuntimeError(
                f"refusing to replace non-symlink figure index entry: {link}")
        temporary = figure_root / f".{name}.tmp.{job_id}"
        if temporary.exists() or temporary.is_symlink():
            temporary.unlink()
        temporary.symlink_to(
            os.path.relpath(target, figure_root.resolve()),
            target_is_directory=True)
        os.replace(temporary, link)
        kept.append((f"FIGURE_INDEX:{name}", str(target)))

    # ---- Per-library generated profiles and versioned swap results ----
    for lib_num in libraries:
        demux = (mapping_root / f"{spec['library_prefix']}{lib_num}" /
                 spec["demux_subdir"])
        swap_root = demux / AMBIENT_SWAP_TEST_DIRNAME
        if swap_root.is_dir() and not swap_root.is_symlink():
            generations = [
                path for path in swap_root.glob("auto_*")
                if path.is_dir() and not path.is_symlink()]
            ranked = [(_cleanup_profile_rank(path), path) for path in generations]
            keep_generation = max(
                ranked,
                key=lambda item: (
                    item[0][0], item[0][1], str(item[1])),
                default=None)
            if keep_generation:
                kept.append((
                    f"lib{lib_num}:AMBIENT_SWAP_TEST",
                    str(keep_generation[1])))
                for _, path in ranked:
                    if path.resolve() != keep_generation[1].resolve():
                        remove_tree(
                            path, swap_root,
                            f"lib{lib_num}_swap_results",
                            "older per-library swap-test generation")

        for candidate_set in ("applied", "exploratory", "supported"):
            validation_root = (
                demux / IDENTITY_AMBIENT_DIRNAME / candidate_set /
                AMBIENT_VALIDATION_DIRNAME)
            fixed_dir = validation_root / "fixed_profile"
            plan_dir = validation_root / "plan"
            if not fixed_dir.is_dir() or fixed_dir.is_symlink():
                continue
            protected = set()
            for plan_path in sorted(plan_dir.glob("*.json")):
                try:
                    with plan_path.open("r", encoding="utf-8") as handle:
                        plan = json.load(handle)
                except (OSError, json.JSONDecodeError):
                    continue
                for key in ("fixed_profile", "fixed_diagnostics"):
                    value = str(plan.get(key, "")).strip()
                    if value:
                        path = Path(value)
                        if is_within(path, fixed_dir):
                            protected.add(path.resolve())
            if not protected:
                continue
            for path in sorted(fixed_dir.iterdir()):
                if (path.is_file() and not path.is_symlink() and
                        path.suffix in {".contam_prof", ".tsv"} and
                        path.resolve() not in protected):
                    remove_file(
                        path, f"lib{lib_num}_fixed_profiles",
                        "not referenced by the current validation plan")

    # ---- All-stage logs: retain newest attempt per logical job/task ----
    if log_root.is_dir() and not log_root.is_symlink():
        attempts = {}
        one_job = re.compile(r"^(.*)_(\d+)$")
        array_job = re.compile(r"^(.*)_(\d+)_(\d+)$")
        for root_text, dirnames, filenames in os.walk(log_root):
            dirnames[:] = [
                name for name in dirnames
                if not Path(root_text, name).is_symlink()]
            relative_root = str(Path(root_text).resolve().relative_to(log_root))
            for filename in filenames:
                extension = Path(filename).suffix.lower()
                if extension not in {".out", ".err"}:
                    continue
                stem = filename[:-len(extension)]
                match = array_job.match(stem)
                if match:
                    logical = (relative_root, match.group(1), match.group(3))
                    attempt = match.group(2)
                else:
                    match = one_job.match(stem)
                    if not match:
                        continue
                    logical = (relative_root, match.group(1), "")
                    attempt = match.group(2)
                path = Path(root_text, filename)
                key = (logical, attempt)
                attempts.setdefault(key, []).append(path)
        by_logical = {}
        for (logical, attempt), paths in attempts.items():
            score = max(path.stat().st_mtime_ns for path in paths)
            by_logical.setdefault(logical, []).append(
                (score, int(attempt), paths))
        for logical, candidates in by_logical.items():
            winner = max(candidates, key=lambda item: (item[0], item[1]))
            for candidate in candidates:
                if candidate is winner:
                    continue
                for path in candidate[2]:
                    remove_file(
                        path, "logs",
                        "older attempt for logical job " + "/".join(logical))

    # ---- Generated orchestration artifacts: retain newest hashed variant ----
    if script_root.is_dir() and not script_root.is_symlink():
        grouped = {}
        hash_token = re.compile(r"(?<![0-9A-Fa-f])[0-9A-Fa-f]{8,64}(?![0-9A-Fa-f])")
        for path in sorted(script_root.iterdir()):
            if path.is_symlink() or not path.is_file():
                continue
            if ".tmp." in path.name:
                remove_file(
                    path, "slurm_scripts", "abandoned temporary artifact")
                continue
            if path.suffix.lower() not in {".sbatch", ".json", ".tsv"}:
                continue
            normalized = hash_token.sub("HASH", path.name)
            grouped.setdefault(normalized, []).append(path)
        for normalized, paths in grouped.items():
            if len(paths) < 2:
                continue
            winner = max(
                paths, key=lambda path: (path.stat().st_mtime_ns, path.name))
            for path in paths:
                if path != winner:
                    remove_file(
                        path, "slurm_scripts",
                        f"older generated artifact family {normalized}")

    protected_roots = sorted(expected_names & {
        child.name for child in aggregate_root.iterdir()
        if child.exists() or child.is_symlink()})
    summary_path = aggregate_root / "cleanup_results_summary.tsv"
    tmp_summary = summary_path.with_name(
        f".{summary_path.name}.tmp.{job_id}")
    with tmp_summary.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(("scope", "action", "path", "reason"))
        for row in records:
            writer.writerow(row)
        for scope, path in kept:
            writer.writerow((scope, "kept", path, "newest completed result"))
    os.replace(tmp_summary, summary_path)

    state = {
        "schema_version": 1,
        "orchestrator_release": ORCHESTRATOR_RELEASE,
        "completed_at": datetime.now().astimezone().isoformat(),
        "removed_count": len(records),
        "kept_results": [
            {"scope": scope, "path": path} for scope, path in kept],
        "protected_aggregate_roots": protected_roots,
        "summary": str(summary_path),
    }
    state_path = aggregate_root / "cleanup_results_state.json"
    tmp_state = state_path.with_name(f".{state_path.name}.tmp.{job_id}")
    with tmp_state.open("w", encoding="utf-8") as handle:
        json.dump(state, handle, indent=2, sort_keys=True)
        handle.write("\n")
    os.replace(tmp_state, state_path)
    print(f"CLEANUP_RESULTS removed {len(records)} obsolete item(s)")
    print(f"Cleanup ledger: {summary_path}")
    print("Canonical in-place stage outputs and shared inputs were preserved")


def generate_cleanup_results_script(lib_nums):
    """Generate one all-stage cleanup job for the production result tree."""
    payload = {
        "schema_version": 1,
        "aggregate_root": AGGREGATE_ROOT,
        "mapping_root": BEEGFS_ROOT,
        "ambient_root": AMBIENT_PLOT_ROOT,
        "figure_root": FIGURE_ROOT,
        "log_root": get_log_dir(),
        "script_root": get_script_dir(),
        "libraries": sorted(set(int(value) for value in lib_nums)),
        "library_prefix": LIB_PREFIX,
        "demux_subdir": DEMUX_SUBDIR,
        "condition_names": ALL_CONDITION_ABBREVS,
    }
    spec_path = os.path.join(get_script_dir(), "cleanup_results.worker.json")
    _ambient_write_json(spec_path, payload)
    log_dir = os.path.join(get_log_dir(), "CLEANUP_RESULTS")
    script_path = os.path.join(get_script_dir(), "cleanup_results.sbatch")
    script = f'''#!/bin/bash
#SBATCH --job-name=cleanup_results
#SBATCH --output={log_dir}/cleanup_results_%j.out
#SBATCH --error={log_dir}/cleanup_results_%j.err
#SBATCH --chdir={AGGREGATE_ROOT}
#SBATCH --time={SLURM_TIME}
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --partition={SLURM_PARTITION}
#SBATCH --nodes=1

set -euo pipefail
{module_block(AMBIENT_PYTHON_MODULES)}
python3 "{os.path.abspath(__file__)}" \
  --_cleanup-results-worker "{spec_path}"
test -s "{os.path.join(AGGREGATE_ROOT, 'cleanup_results_summary.tsv')}"
test -s "{os.path.join(AGGREGATE_ROOT, 'cleanup_results_state.json')}"
'''
    os.makedirs(log_dir, exist_ok=True)
    _ambient_write_text(script_path, script, executable=True)
    return {"script": script_path, "spec": spec_path}


# =============================================================================
# Job submission
# =============================================================================

def submit_job(script_path):
    """Submit one sbatch script and return its job ID.

    Use Slurm's machine-readable response and retry once for a transient
    scheduler/transport failure.  A second failure is real: the requested job
    was not submitted, so dependent work cannot safely continue.
    """
    last_detail = ""
    for attempt in (1, 2):
        try:
            result = subprocess.run(
                ["sbatch", "--parsable", f"--chdir={AGGREGATE_ROOT}",
                 script_path],
                capture_output=True, text=True, check=False)
        except OSError as exc:
            last_detail = str(exc)
            result = None

        if result is not None and result.returncode == 0:
            # --parsable returns JOBID or JOBID;cluster.
            job_id = result.stdout.strip().split(";", 1)[0]
            if job_id.isdigit():
                return job_id
            last_detail = f"unexpected sbatch response: {result.stdout.strip()!r}"
        elif result is not None:
            last_detail = (
                result.stderr.strip() or result.stdout.strip()
                or f"exit {result.returncode}")

        if attempt == 1:
            print(f"  WARNING: sbatch attempt 1 failed for {script_path}: {last_detail}")
            print("  Retrying once...")

    raise RuntimeError(f"sbatch did not submit {script_path}: {last_detail}")


# =============================================================================
# GEX_AMBIENT: gene-expression ambient profile and concentration diagnostics
# =============================================================================

def gex_add_internal_worker_argument(parser):
    parser.add_argument(
        "--_gex-ambient-summary-worker", metavar="SPEC_JSON", default=None,
        help=argparse.SUPPRESS)
    parser.add_argument(
        "--_gex-auto-cluster-worker", metavar="SPEC_JSON", default=None,
        help=argparse.SUPPRESS)
    parser.add_argument(
        "--_gex-cluster-intersection-worker", metavar="SPEC_JSON", default=None,
        help=argparse.SUPPRESS)


def gex_maybe_run_internal_worker(args):
    intersection_spec = getattr(
        args, "_gex_cluster_intersection_worker", None)
    if intersection_spec:
        gex_run_cluster_intersection_worker(intersection_spec)
        return True
    cluster_spec = getattr(args, "_gex_auto_cluster_worker", None)
    if cluster_spec:
        gex_run_auto_cluster_worker(cluster_spec)
        return True
    spec_path = getattr(args, "_gex_ambient_summary_worker", None)
    if not spec_path:
        return False
    gex_run_summary_worker(spec_path)
    return True


def _gex_extract_library_number(value):
    text = str(value).strip()
    if text.lower().startswith("lib"):
        text = text[3:]
    else:
        text = text.split("_")[-1]
    return int(text)


def _gex_canonical_barcode(value):
    """Match the barcode convention already used by PLOIDY_NN output."""
    return str(value).strip().split("-", 1)[0]


def gex_run_auto_cluster_worker(spec_path):
    """Create RNA-only numbered nuisance clusters or reuse an H5AD column."""
    import gzip
    import numpy as np

    try:
        import scipy.io
        import scipy.sparse
        from sklearn.decomposition import TruncatedSVD
        from sklearn.metrics import adjusted_rand_score, silhouette_score
        from sklearn.neighbors import KNeighborsClassifier, NearestNeighbors
    except ImportError as exc:
        raise RuntimeError(
            "GEX automatic clustering requires scipy and scikit-learn "
            "from the configured genomics-base environment") from exc

    with open(spec_path, "r", encoding="utf-8") as handle:
        spec = json.load(handle)
    h5ad_path = spec["h5ad"]
    lib_num = int(spec["library"])
    source_mode = str(spec["source_mode"])
    n_components = int(spec["pca_components"])
    seed = int(spec["random_seed"])
    output = Path(spec["output"])
    qc_output = Path(spec["qc_output"])
    contract_output = Path(spec["contract_output"])

    with gzip.open(spec["mex_barcodes"], "rt", encoding="utf-8",
                   errors="replace") as handle:
        mex_barcode_list = [
            _gex_canonical_barcode(line) for line in handle if line.strip()]
    if not mex_barcode_list:
        raise RuntimeError(f"filtered MEX barcode file is empty: {spec['mex_barcodes']}")
    if len(set(mex_barcode_list)) != len(mex_barcode_list):
        raise RuntimeError(
            f"lib{lib_num}: filtered MEX contains duplicate canonical barcodes")
    mex_index = {barcode: idx for idx, barcode in enumerate(mex_barcode_list)}
    mex_barcodes = set(mex_barcode_list)

    adata = None
    libraries = np.asarray([], dtype=int)
    if source_mode == "auto":
        # The per-library filtered MEX is authoritative. Do not inherit prior
        # notebook or NN-H5AD cell filtering in automatic mode.
        selected = np.arange(len(mex_barcode_list), dtype=int)
        selected_barcodes = np.asarray(mex_barcode_list, dtype=object)
        mex_columns = selected.copy()
        n_h5ad_library_cells = 0
    elif source_mode == "h5ad":
        try:
            import anndata as ad
        except ImportError as exc:
            raise RuntimeError(
                "--gex-cluster-source h5ad requires anndata from the configured "
                "genomics-base environment") from exc
        adata = ad.read_h5ad(h5ad_path, backed="r")
        if "library" not in adata.obs.columns:
            raise RuntimeError(f"NN H5AD has no obs['library'] column: {h5ad_path}")
        libraries = np.asarray([
            _gex_extract_library_number(value)
            for value in adata.obs["library"].astype(str).values], dtype=int)
        obs_names = np.asarray(adata.obs_names.astype(str))
        canonical = np.asarray([
            _gex_canonical_barcode(value) for value in obs_names], dtype=object)
        selected = np.flatnonzero(
            (libraries == lib_num) & np.isin(canonical, list(mex_barcodes)))
        if selected.size < max(2, int(spec["minimum_cells"])):
            raise RuntimeError(
                f"lib{lib_num}: only {selected.size} NN-H5AD cells overlap the "
                f"filtered MEX; minimum is {spec['minimum_cells']}")
        selected_barcodes = canonical[selected]
        if len(set(selected_barcodes.tolist())) != selected.size:
            raise RuntimeError(
                f"lib{lib_num}: NN H5AD contains duplicate canonical barcodes")
        mex_columns = np.asarray(
            [mex_index[barcode] for barcode in selected_barcodes], dtype=int)
        n_h5ad_library_cells = int((libraries == lib_num).sum())
    else:
        raise RuntimeError(f"unknown GEX cluster source mode: {source_mode}")

    def read_features():
        rows = []
        opener = gzip.open if str(spec["mex_features"]).endswith(".gz") else open
        with opener(spec["mex_features"], "rt", encoding="utf-8",
                    errors="replace") as handle:
            for line in handle:
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 2:
                    continue
                rows.append((fields[1], fields[2] if len(fields) > 2 else ""))
        return rows

    def load_expression():
        feature_rows = read_features()
        matrix_opener = (
            gzip.open if str(spec["mex_matrix"]).endswith(".gz") else open)
        with matrix_opener(spec["mex_matrix"], "rb") as handle:
            matrix = scipy.io.mmread(handle)
        matrix = scipy.sparse.csr_matrix(matrix)
        if matrix.shape[0] != len(feature_rows):
            raise RuntimeError(
                f"lib{lib_num}: MEX feature rows ({len(feature_rows)}) do not "
                f"match matrix rows ({matrix.shape[0]})")
        if matrix.shape[1] != len(mex_barcode_list):
            raise RuntimeError(
                f"lib{lib_num}: MEX barcode rows ({len(mex_barcode_list)}) do "
                f"not match matrix columns ({matrix.shape[1]})")
        gex_rows = [
            idx for idx, (_, feature_type) in enumerate(feature_rows)
            if not feature_type or feature_type == spec["feature_type"]]
        if not gex_rows:
            raise RuntimeError(
                f"lib{lib_num}: no {spec['feature_type']!r} features in MEX")
        genes = np.asarray([feature_rows[idx][0] for idx in gex_rows],
                           dtype=object)
        expression = matrix[gex_rows, :][:, mex_columns].T.tocsr()
        expression.sum_duplicates()
        expression.eliminate_zeros()
        return expression, genes

    def load_marker_genes():
        path = spec.get("marker_genes")
        if not path:
            return set()
        markers = set()
        with open(path, "r", encoding="utf-8", errors="replace") as handle:
            for line in handle:
                line = line.strip()
                if line and not line.startswith("#"):
                    markers.add(re.split(r"[\t, ]+", line, maxsplit=1)[0])
        return markers

    def make_embedding(expression, genes, anchor_mask):
        anchor_expression = expression[anchor_mask, :].tocsr()
        detected = np.diff(anchor_expression.tocsc().indptr)
        min_detected = max(3, int(math.ceil(anchor_expression.shape[0] * 0.001)))
        technical = np.asarray([
            bool(re.search(spec["clustering_exclude_regex"], str(gene)))
            if spec["clustering_exclude_regex"] else False
            for gene in genes], dtype=bool)
        eligible = (detected >= min_detected) & ~technical
        eligible_idx = np.flatnonzero(eligible)
        if eligible_idx.size < 2:
            raise RuntimeError(
                f"lib{lib_num}: fewer than two nontechnical expressed genes "
                "remain for automatic clustering")

        # Poisson deviance ranks departures from constant relative expression
        # directly on raw UMI counts.  Unlike dispersion-HVG selection it has no
        # log-normalization/pseudocount-dependent mean-variance fit.
        cell_totals = np.asarray(anchor_expression.sum(axis=1)).ravel()
        gene_totals = np.asarray(anchor_expression.sum(axis=0)).ravel()
        grand_total = float(gene_totals.sum())
        if grand_total <= 0 or np.any(cell_totals <= 0):
            raise RuntimeError(
                f"lib{lib_num}: zero-count anchor encountered in raw MEX")
        deviance = np.full(expression.shape[1], -np.inf, dtype=np.float64)
        anchor_csc = anchor_expression[:, eligible_idx].tocsc()
        for local_idx, gene_idx in enumerate(eligible_idx):
            start, end = anchor_csc.indptr[local_idx:local_idx + 2]
            rows = anchor_csc.indices[start:end]
            values = anchor_csc.data[start:end].astype(np.float64, copy=False)
            proportion = gene_totals[gene_idx] / grand_total
            expected = cell_totals[rows] * proportion
            good = (values > 0) & (expected > 0)
            if np.any(good):
                deviance[gene_idx] = 2.0 * np.sum(
                    values[good] * np.log(values[good] / expected[good]))

        data_feature_count = int(spec["data_features"])
        if data_feature_count <= 0 or data_feature_count >= eligible_idx.size:
            chosen = set(eligible_idx.tolist())
            feature_selection = "all_detected_nontechnical_genes"
        else:
            ranked = eligible_idx[np.argsort(deviance[eligible_idx])[::-1]]
            chosen = set(ranked[:data_feature_count].tolist())
            feature_selection = "poisson_deviance"

        markers = load_marker_genes()
        markers_present = set()
        if markers:
            for idx, gene in enumerate(genes):
                if str(gene) in markers and detected[idx] > 0 and not technical[idx]:
                    chosen.add(idx)
                    markers_present.add(str(gene))
        feature_idx = np.asarray(sorted(chosen), dtype=int)
        if feature_idx.size < 2:
            raise RuntimeError(f"lib{lib_num}: too few clustering features")

        cell_totals_all = np.asarray(expression.sum(axis=1)).ravel()
        normalized = expression[:, feature_idx].astype(np.float32, copy=True)
        normalized = scipy.sparse.diags(
            10000.0 / np.maximum(cell_totals_all, 1.0)) @ normalized
        normalized = normalized.tocsr()
        normalized.data = np.log1p(normalized.data)
        max_components = min(
            n_components, int(anchor_mask.sum()) - 1, feature_idx.size - 1)
        if max_components < 2:
            raise RuntimeError(f"lib{lib_num}: too few dimensions for clustering")
        reducer = TruncatedSVD(n_components=max_components, random_state=seed)
        reducer.fit(normalized[anchor_mask, :])
        result = np.asarray(reducer.transform(normalized), dtype=np.float32)
        return result, feature_selection, feature_idx, markers, markers_present

    expression = None
    genes = None
    embedding = None
    embedding_source = "not_required"
    feature_selection = "not_required"
    feature_idx = np.asarray([], dtype=int)
    markers_requested = set()
    markers_present = set()
    cell_gene_counts = None

    raw_labels = None
    if source_mode == "h5ad":
        column = str(spec["h5ad_cluster_column"])
        if column not in adata.obs.columns:
            raise RuntimeError(
                f"NN H5AD has no obs[{column!r}] column: {h5ad_path}")
        raw_labels = np.asarray(adata.obs[column].astype(object).values[selected],
                                dtype=object)
        label_present = np.asarray([
            value is not None and str(value).strip() not in {"", "nan", "None", "<NA>"}
            for value in raw_labels], dtype=bool)
        if not np.any(label_present):
            raise RuntimeError(
                f"lib{lib_num}: H5AD column {column!r} has no usable labels")
        label_names = sorted({str(value) for value in raw_labels[label_present]})
        label_to_int = {name: idx for idx, name in enumerate(label_names)}
        labels = np.full(selected.size, -1, dtype=int)
        labels[label_present] = [
            label_to_int[str(value)] for value in raw_labels[label_present]]
        anchor_mask = label_present
        candidate_metrics = []
        chosen_resolution = None
        chosen_stability = None
        chosen_silhouette = None
        chosen_modularity = None
        if not np.all(label_present):
            expression, genes = load_expression()
            cell_gene_counts = np.diff(expression.indptr)
            embedding, feature_selection, feature_idx, markers_requested, markers_present = (
                make_embedding(expression, genes, anchor_mask))
            embedding_source = "raw_mex_log1p_truncated_svd_for_label_transfer"
            classifier = KNeighborsClassifier(
                n_neighbors=min(15, int(anchor_mask.sum())), weights="distance")
            classifier.fit(embedding[anchor_mask], labels[anchor_mask])
            labels[~anchor_mask] = classifier.predict(embedding[~anchor_mask])
            confidence = np.ones(selected.size, dtype=float)
            confidence[~anchor_mask] = np.max(
                classifier.predict_proba(embedding[~anchor_mask]), axis=1)
        else:
            confidence = np.ones(selected.size, dtype=float)
    elif source_mode == "auto":
        expression, genes = load_expression()
        cell_gene_counts = np.diff(expression.indptr)
        anchor_mask = cell_gene_counts >= int(spec["minimum_genes"])
        if anchor_mask.sum() < max(2, int(spec["minimum_cells"]) * 2):
            raise RuntimeError(
                f"lib{lib_num}: only {int(anchor_mask.sum())} cells meet the "
                f"permissive {spec['minimum_genes']}-gene anchor floor")
        embedding, feature_selection, feature_idx, markers_requested, markers_present = (
            make_embedding(expression, genes, anchor_mask))
        embedding_source = "raw_mex_poisson_deviance_log1p_truncated_svd"
        if (embedding.ndim != 2 or embedding.shape[0] != selected.size or
                embedding.shape[1] < 2 or not np.isfinite(embedding).all()):
            raise RuntimeError(
                f"lib{lib_num}: invalid clustering embedding shape {embedding.shape}")

        try:
            import igraph as ig
            import leidenalg
        except ImportError as exc:
            raise RuntimeError(
                "automatic GEX clustering requires python-igraph and leidenalg "
                "from the configured genomics-base environment") from exc

        anchor_embedding = embedding[anchor_mask]
        n_neighbors = min(
            int(spec["neighbors"]), anchor_embedding.shape[0] - 1)
        if n_neighbors < 2:
            raise RuntimeError(f"lib{lib_num}: too few anchor cells for neighbors")
        neighbor_model = NearestNeighbors(
            n_neighbors=n_neighbors + 1, metric="euclidean", n_jobs=-1)
        neighbor_model.fit(anchor_embedding)
        distances, neighbors = neighbor_model.kneighbors(anchor_embedding)
        edge_weights = {}
        for row in range(anchor_embedding.shape[0]):
            for distance, col in zip(distances[row, 1:], neighbors[row, 1:]):
                left, right = sorted((int(row), int(col)))
                if left == right:
                    continue
                weight = 1.0 / (1.0 + float(distance))
                edge_weights[(left, right)] = max(
                    weight, edge_weights.get((left, right), 0.0))
        graph = ig.Graph(
            n=anchor_embedding.shape[0], edges=list(edge_weights), directed=False)
        graph.es["weight"] = list(edge_weights.values())

        resolutions = [
            float(value) for value in str(spec["resolution_grid"]).split(",")
            if value.strip()]
        repeats = int(spec["stability_repeats"])
        candidate_metrics = []
        candidate_labels = {}
        rng = np.random.default_rng(seed)
        silhouette_indices = np.arange(anchor_embedding.shape[0])
        if silhouette_indices.size > 5000:
            silhouette_indices = np.sort(rng.choice(
                silhouette_indices, size=5000, replace=False))
        for resolution in resolutions:
            repeat_labels = []
            for repeat in range(repeats):
                partition = leidenalg.find_partition(
                    graph, leidenalg.RBConfigurationVertexPartition,
                    weights="weight", resolution_parameter=resolution,
                    seed=seed + repeat, n_iterations=-1)
                repeat_labels.append(np.asarray(partition.membership, dtype=int))
            first = repeat_labels[0]
            cluster_counts = np.bincount(first)
            n_clusters = int(cluster_counts.size)
            stability_pairs = [
                adjusted_rand_score(repeat_labels[left], repeat_labels[right])
                for left in range(repeats)
                for right in range(left + 1, repeats)]
            stability = float(np.mean(stability_pairs)) if stability_pairs else 1.0
            silhouette = float("nan")
            if 1 < n_clusters < silhouette_indices.size:
                silhouette = float(silhouette_score(
                    anchor_embedding[silhouette_indices],
                    first[silhouette_indices], metric="euclidean"))
            modularity = float(graph.modularity(first.tolist(), weights="weight"))
            metric = {
                "resolution": resolution, "clusters": n_clusters,
                "minimum_cluster_cells": int(cluster_counts.min()),
                "stability_ari": stability, "silhouette": silhouette,
                "modularity": modularity,
            }
            candidate_metrics.append(metric)
            candidate_labels[resolution] = first

        valid = [
            item for item in candidate_metrics
            if item["clusters"] >= 2 and
            item["minimum_cluster_cells"] >= int(spec["minimum_cells"]) and
            math.isfinite(item["silhouette"])]
        if not valid:
            raise RuntimeError(
                f"lib{lib_num}: no Leiden resolution produced at least two "
                f"clusters with >= {spec['minimum_cells']} anchor cells each; "
                "inspect the cluster QC or expand --gex-auto-resolution-grid")
        best_stability = max(item["stability_ari"] for item in valid)
        stable = [
            item for item in valid
            if item["stability_ari"] >= best_stability - 0.05]
        chosen = max(
            stable,
            key=lambda item: (
                item["silhouette"], item["modularity"],
                -item["clusters"], -item["resolution"]))
        chosen_resolution = float(chosen["resolution"])
        chosen_stability = float(chosen["stability_ari"])
        chosen_silhouette = float(chosen["silhouette"])
        chosen_modularity = float(chosen["modularity"])
        anchor_labels = candidate_labels[chosen_resolution]
        labels = np.full(selected.size, -1, dtype=int)
        labels[anchor_mask] = anchor_labels
        confidence = np.ones(selected.size, dtype=float)
        if np.any(~anchor_mask):
            classifier = KNeighborsClassifier(
                n_neighbors=min(15, int(anchor_mask.sum())), weights="distance")
            classifier.fit(anchor_embedding, anchor_labels)
            labels[~anchor_mask] = classifier.predict(embedding[~anchor_mask])
            confidence[~anchor_mask] = np.max(
                classifier.predict_proba(embedding[~anchor_mask]), axis=1)
    else:
        raise RuntimeError(f"unknown GEX cluster source mode: {source_mode}")

    if np.any(labels < 0):
        raise RuntimeError(f"lib{lib_num}: not every selected cell received a cluster")
    n_clusters = len(set(labels.tolist()))
    if n_clusters < 2:
        raise RuntimeError(f"lib{lib_num}: clustering produced fewer than two groups")

    # Deterministic arbitrary numbering.  H5AD-column labels are ordered by
    # their source names; automatic Leiden labels are ordered by centroids.
    if source_mode == "auto":
        centroids = []
        for label in sorted(set(labels.tolist())):
            members = embedding[labels == label]
            centroids.append((tuple(np.mean(members, axis=0).tolist()), label))
        label_order = [old for _, old in sorted(centroids)]
    else:
        label_order = sorted(set(labels.tolist()))
    old_to_new = {old: new for new, old in enumerate(label_order)}
    labels = np.asarray([old_to_new[int(label)] for label in labels], dtype=int)
    n_clusters = len(old_to_new)

    output.parent.mkdir(parents=True, exist_ok=True)
    cluster_tmp = output.with_name(output.name + f".tmp.{os.getpid()}")
    with cluster_tmp.open("w", encoding="utf-8") as handle:
        for barcode, label in sorted(zip(selected_barcodes, labels)):
            handle.write(f"{barcode}\t{label}\n")
    os.replace(cluster_tmp, output)

    counts = {int(label): int((labels == label).sum())
              for label in sorted(set(labels.tolist()))}
    _gex_atomic_tsv(
        qc_output,
        ["library", "cluster", "n_cells", "n_anchor_cells",
         "n_transferred_cells", "median_assignment_confidence",
         "minimum_assignment_confidence", "n_h5ad_library_cells",
         "n_mex_barcodes", "n_overlap_cells", "source_mode",
         "embedding_source", "embedding_dimensions", "feature_selection",
         "selected_features", "marker_genes_requested",
         "marker_genes_present", "realized_clusters", "chosen_resolution",
         "stability_ari", "silhouette", "modularity", "random_seed"],
        ({
            "library": lib_num, "cluster": label, "n_cells": count,
            "n_anchor_cells": int(((labels == label) & anchor_mask).sum()),
            "n_transferred_cells": int(((labels == label) & ~anchor_mask).sum()),
            "median_assignment_confidence": float(np.median(confidence[labels == label])),
            "minimum_assignment_confidence": float(np.min(confidence[labels == label])),
            "n_h5ad_library_cells": n_h5ad_library_cells,
            "n_mex_barcodes": len(mex_barcode_list),
            "n_overlap_cells": int(selected.size),
            "source_mode": source_mode, "embedding_source": embedding_source,
            "embedding_dimensions": int(embedding.shape[1]) if embedding is not None else 0,
            "feature_selection": feature_selection,
            "selected_features": int(feature_idx.size),
            "marker_genes_requested": len(markers_requested),
            "marker_genes_present": len(markers_present),
            "realized_clusters": n_clusters,
            "chosen_resolution": chosen_resolution,
            "stability_ari": chosen_stability,
            "silhouette": chosen_silhouette,
            "modularity": chosen_modularity, "random_seed": seed,
        } for label, count in counts.items()))
    contract = dict(spec)
    contract.update({
        "schema_version": 2, "source_mode": source_mode,
        "embedding_source": embedding_source,
        "embedding_dimensions": int(embedding.shape[1]) if embedding is not None else 0,
        "feature_selection": feature_selection,
        "selected_features": int(feature_idx.size),
        "marker_genes_requested": sorted(markers_requested),
        "marker_genes_present": sorted(markers_present),
        "n_h5ad_library_cells": n_h5ad_library_cells,
        "n_mex_barcodes": len(mex_barcode_list),
        "n_overlap_cells": int(selected.size), "realized_clusters": n_clusters,
        "n_anchor_cells": int(anchor_mask.sum()),
        "n_transferred_cells": int((~anchor_mask).sum()),
        "candidate_metrics": candidate_metrics,
        "chosen_resolution": chosen_resolution,
        "cluster_labels_are_semantic": False,
    })
    _ambient_write_json(contract_output, contract)
    print(
        f"GEX_AMBIENT {source_mode} clusters: lib{lib_num}, {selected.size} "
        f"cells, {n_clusters} numbered groups -> {output}")


def generate_gex_auto_cluster_script(lib_num, args, force=False):
    """Generate one RNA/H5AD-derived numbered-cluster job per library."""
    output = get_gex_cluster_path(lib_num, args)
    qc_output = get_gex_auto_cluster_qc_path(lib_num, args)
    contract_output = get_gex_auto_cluster_contract_path(lib_num, args)
    summary_dir = get_gex_ambient_summary_dir(args)
    spec_dir = os.path.join(summary_dir, "cluster_specs")
    os.makedirs(spec_dir, exist_ok=True)
    spec_path = os.path.join(spec_dir, f"lib{lib_num}.auto_cluster_spec.json")
    _ambient_write_json(spec_path, {
        "schema_version": 2,
        "library": lib_num,
        "source_mode": args.gex_cluster_source,
        "h5ad": (os.path.abspath(args.ploidy_nn_h5ad)
                 if args.gex_cluster_source == "h5ad" else None),
        "mex_barcodes": get_filtered_barcodes(lib_num),
        "mex_features": get_filtered_features(lib_num),
        "mex_matrix": get_filtered_matrix(lib_num),
        "feature_type": args.gex_feature_type,
        "output": output,
        "qc_output": qc_output,
        "contract_output": contract_output,
        "h5ad_cluster_column": args.gex_h5ad_cluster_column,
        "marker_genes": os.path.abspath(args.gex_marker_genes)
        if args.gex_marker_genes else None,
        "data_features": args.gex_auto_data_features,
        "minimum_genes": args.gex_auto_min_genes,
        "neighbors": args.gex_auto_neighbors,
        "pca_components": args.gex_auto_pca_components,
        "resolution_grid": args.gex_auto_resolution_grid,
        "stability_repeats": args.gex_auto_stability_repeats,
        "clustering_exclude_regex": args.gex_auto_exclude_genes_regex,
        "random_seed": args.gex_auto_random_seed,
        "minimum_cells": args.gex_min_cluster_cells,
        "cluster_labels_are_semantic": False,
    })
    analysis = _gex_safe_name(args.gex_analysis_name)
    log_dir = os.path.join(get_log_dir(), "GEX_AMBIENT", analysis, "clusters")
    script_path = os.path.join(
        get_script_dir(), f"gex_auto_clusters_{analysis}_lib{lib_num}.sbatch")
    force_note = "regenerate" if force else "create"
    script = f'''#!/bin/bash
#SBATCH --job-name=gexcl_{analysis[:20]}_l{lib_num}
#SBATCH --output={log_dir}/gex_auto_clusters_lib{lib_num}_%j.out
#SBATCH --error={log_dir}/gex_auto_clusters_lib{lib_num}_%j.err
#SBATCH --time={SLURM_TIME}
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --partition={SLURM_PARTITION}
#SBATCH --nodes=1

set -euo pipefail
{module_block(AMBIENT_PYTHON_MODULES)}
mkdir -p "{log_dir}" "{os.path.dirname(output)}"
echo "GEX_AMBIENT {args.gex_cluster_source} clustering: {force_note} lib{lib_num} numbered groups"
python3 {_ambient_shell(os.path.abspath(__file__))} --_gex-auto-cluster-worker {_ambient_shell(spec_path)}
test -s {_ambient_shell(output)}
test -s {_ambient_shell(qc_output)}
test -s {_ambient_shell(contract_output)}
'''
    _ambient_write_text(script_path, script, executable=True)
    return script_path


def gex_run_cluster_intersection_worker(spec_path):
    """Write the exact cluster/MEX/valid-contamination intersection for a fit."""
    import gzip

    with open(spec_path, "r", encoding="utf-8") as handle:
        spec = json.load(handle)

    def canonical(value):
        return _gex_canonical_barcode(value)

    mex_barcode_list = []
    opener = gzip.open if str(spec["mex_barcodes"]).endswith(".gz") else open
    with opener(spec["mex_barcodes"], "rt", encoding="utf-8",
                errors="replace") as handle:
        for line in handle:
            if line.strip():
                mex_barcode_list.append(canonical(line))
    if len(set(mex_barcode_list)) != len(mex_barcode_list):
        raise RuntimeError(
            "filtered MEX has duplicate canonical barcodes; cohort subsetting "
            "would be ambiguous")
    mex_barcodes = set(mex_barcode_list)

    valid_rates = set()
    invalid_rate_rows = 0
    duplicate_rate_rows = 0
    with open(spec["contam_rate"], "r", encoding="utf-8",
              errors="replace") as handle:
        for line in handle:
            fields = line.split()
            if len(fields) < 2:
                if fields:
                    invalid_rate_rows += 1
                continue
            try:
                rate = float(fields[1])
            except ValueError:
                invalid_rate_rows += 1
                continue
            if not math.isfinite(rate) or not 0.0 <= rate < 1.0:
                invalid_rate_rows += 1
                continue
            barcode = canonical(fields[0])
            if barcode in valid_rates:
                duplicate_rate_rows += 1
            valid_rates.add(barcode)

    source_clusters = {}
    duplicate_cluster_rows = 0
    malformed_cluster_rows = 0
    with open(spec["source_clusters"], "r", encoding="utf-8",
              errors="replace") as handle:
        for line in handle:
            fields = line.split()
            if not fields or fields[0].startswith("#"):
                continue
            if len(fields) != 2:
                malformed_cluster_rows += 1
                continue
            barcode = canonical(fields[0])
            label = fields[1]
            if barcode in source_clusters:
                duplicate_cluster_rows += 1
                continue
            source_clusters[barcode] = label

    if malformed_cluster_rows or duplicate_cluster_rows:
        raise RuntimeError(
            "cluster source must be headerless unique-barcode two-column text; "
            f"malformed={malformed_cluster_rows}, duplicates={duplicate_cluster_rows}")
    if invalid_rate_rows or duplicate_rate_rows:
        raise RuntimeError(
            "contamination-rate input must contain unique barcodes and finite "
            f"rates in [0,1); invalid={invalid_rate_rows}, "
            f"duplicates={duplicate_rate_rows}")

    eligible = {
        barcode: label for barcode, label in source_clusters.items()
        if barcode in mex_barcodes and barcode in valid_rates}
    labels = sorted(set(eligible.values()))
    minimum = int(spec["minimum_cells"])
    if len(eligible) < minimum:
        raise RuntimeError(
            f"only {len(eligible)} clustered cells have expression and a valid "
            f"contamination rate; minimum is {minimum}")
    if len(labels) < 2:
        raise RuntimeError(
            "fewer than two clusters remain after intersecting expression and "
            "valid contamination-rate cells")

    label_counts = {
        label: sum(value == label for value in eligible.values())
        for label in labels}
    undersized = {
        label: count for label, count in label_counts.items()
        if count < minimum}
    if undersized:
        detail = ", ".join(
            f"{label}={count}" for label, count in sorted(undersized.items()))
        raise RuntimeError(
            "clusters below --gex-min-cluster-cells after intersecting "
            f"expression and valid contamination rates: {detail}; "
            f"minimum per cluster is {minimum}")

    output = Path(spec["output"])
    output.parent.mkdir(parents=True, exist_ok=True)
    tmp = output.with_name(output.name + f".tmp.{os.getpid()}")
    with tmp.open("w", encoding="utf-8") as handle:
        for barcode in sorted(eligible):
            handle.write(f"{barcode}\t{eligible[barcode]}\n")
    os.replace(tmp, output)

    # ambient_rna_gex iterates every matrix barcode and uses operator[] for
    # contamination rates, so a missing rate otherwise becomes an implicit
    # zero even when that barcode is absent from the cluster file.  Build the
    # exact eligible-column MEX passed to C++ to make that impossible.
    eligible_columns = {}
    eligible_order = []
    for old_column, barcode in enumerate(mex_barcode_list, start=1):
        if barcode in eligible:
            eligible_order.append(barcode)
            eligible_columns[old_column] = len(eligible_order)

    subset_barcodes = Path(spec["subset_barcodes"])
    subset_matrix = Path(spec["subset_matrix"])
    subset_barcodes.parent.mkdir(parents=True, exist_ok=True)
    barcode_tmp = subset_barcodes.with_name(
        subset_barcodes.name + f".tmp.{os.getpid()}")
    with gzip.open(barcode_tmp, "wt", encoding="utf-8") as handle:
        for barcode in eligible_order:
            handle.write(barcode + "\n")
    os.replace(barcode_tmp, subset_barcodes)

    matrix_opener = (
        gzip.open if str(spec["mex_matrix"]).endswith(".gz") else open)

    def read_matrix_header(handle):
        banner = handle.readline()
        if not banner.startswith("%%MatrixMarket matrix coordinate"):
            raise RuntimeError(
                f"unsupported MatrixMarket banner in {spec['mex_matrix']}: "
                f"{banner.strip()!r}")
        comments = []
        while True:
            line = handle.readline()
            if not line:
                raise RuntimeError("MatrixMarket file ends before dimensions")
            if line.startswith("%"):
                comments.append(line)
                continue
            fields = line.split()
            if len(fields) != 3:
                raise RuntimeError(
                    f"invalid MatrixMarket dimensions: {line.strip()!r}")
            try:
                dimensions = tuple(int(value) for value in fields)
            except ValueError as exc:
                raise RuntimeError(
                    f"invalid MatrixMarket dimensions: {line.strip()!r}") from exc
            return banner, comments, dimensions

    retained_nnz = 0
    with matrix_opener(spec["mex_matrix"], "rt", encoding="utf-8",
                       errors="strict") as handle:
        _, _, (n_rows, n_columns, declared_nnz) = read_matrix_header(handle)
        if n_columns != len(mex_barcode_list):
            raise RuntimeError(
                f"MEX matrix has {n_columns} columns but barcodes file has "
                f"{len(mex_barcode_list)} rows")
        observed_nnz = 0
        for line in handle:
            if not line.strip():
                continue
            fields = line.split()
            if len(fields) < 3:
                raise RuntimeError(
                    f"malformed MatrixMarket coordinate: {line.strip()!r}")
            try:
                old_column = int(fields[1])
            except ValueError as exc:
                raise RuntimeError(
                    f"invalid MatrixMarket column: {line.strip()!r}") from exc
            observed_nnz += 1
            retained_nnz += int(old_column in eligible_columns)
        if observed_nnz != declared_nnz:
            raise RuntimeError(
                f"MatrixMarket declares {declared_nnz} entries but contains "
                f"{observed_nnz}")

    matrix_tmp = subset_matrix.with_name(
        subset_matrix.name + f".tmp.{os.getpid()}")
    with matrix_opener(spec["mex_matrix"], "rt", encoding="utf-8",
                       errors="strict") as source, gzip.open(
                           matrix_tmp, "wt", encoding="utf-8") as destination:
        banner, comments, _ = read_matrix_header(source)
        destination.write(banner if banner.endswith("\n") else banner + "\n")
        for comment in comments:
            destination.write(
                comment if comment.endswith("\n") else comment + "\n")
        destination.write(f"{n_rows} {len(eligible_order)} {retained_nnz}\n")
        written_nnz = 0
        for line in source:
            if not line.strip():
                continue
            fields = line.split()
            old_column = int(fields[1])
            new_column = eligible_columns.get(old_column)
            if new_column is not None:
                destination.write(
                    " ".join([fields[0], str(new_column), *fields[2:]]) + "\n")
                written_nnz += 1
        if written_nnz != retained_nnz:
            raise RuntimeError(
                f"eligible MatrixMarket write mismatch: expected {retained_nnz}, "
                f"wrote {written_nnz}")
    os.replace(matrix_tmp, subset_matrix)

    _gex_atomic_tsv(
        spec["qc_output"],
        ["metric", "value"],
        [
            {"metric": "source_cluster_barcodes", "value": len(source_clusters)},
            {"metric": "mex_barcodes", "value": len(mex_barcodes)},
            {"metric": "valid_contam_rate_barcodes", "value": len(valid_rates)},
            {"metric": "eligible_cluster_barcodes", "value": len(eligible)},
            {"metric": "eligible_clusters", "value": len(labels)},
            {"metric": "eligible_matrix_nonzero_entries", "value": retained_nnz},
            {"metric": "source_clusters_missing_mex",
             "value": sum(barcode not in mex_barcodes for barcode in source_clusters)},
            {"metric": "source_clusters_missing_valid_contam_rate",
             "value": sum(barcode not in valid_rates for barcode in source_clusters)},
            {"metric": "invalid_contam_rate_rows", "value": invalid_rate_rows},
            {"metric": "duplicate_contam_rate_rows", "value": duplicate_rate_rows},
        ] + [
            {"metric": f"eligible_cluster_{label}_cells", "value": count}
            for label, count in label_counts.items()
        ])
    print(
        f"GEX_AMBIENT eligible cohort: {len(eligible)} cells across "
        f"{len(labels)} clusters -> {output}")


def _gex_atomic_tsv(path, fieldnames, rows):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + f".tmp.{os.getpid()}")
    with tmp.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t",
                                extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    os.replace(tmp, path)


def _gex_rank(values):
    order = sorted(range(len(values)), key=lambda i: values[i])
    ranks = [0.0] * len(values)
    start = 0
    while start < len(order):
        end = start + 1
        while end < len(order) and values[order[end]] == values[order[start]]:
            end += 1
        rank = (start + end - 1) / 2.0 + 1.0
        for pos in range(start, end):
            ranks[order[pos]] = rank
        start = end
    return ranks


def _gex_correlation(x, y):
    if len(x) < 2 or len(x) != len(y):
        return math.nan
    mx = sum(x) / len(x)
    my = sum(y) / len(y)
    dx = [v - mx for v in x]
    dy = [v - my for v in y]
    denom = math.sqrt(sum(v * v for v in dx) * sum(v * v for v in dy))
    return sum(a * b for a, b in zip(dx, dy)) / denom if denom else math.nan


def gex_run_summary_worker(spec_path):
    """Summarize ambient mass concentration without SNP ascertainment filters."""
    with open(spec_path, "r", encoding="utf-8") as handle:
        spec = json.load(handle)
    out_dir = Path(spec["summary_dir"])
    out_dir.mkdir(parents=True, exist_ok=True)
    metrics = []
    top_rows = []
    profiles = {}
    metadata = {}
    thresholds = (0.50, 0.80, 0.90, 0.95, 0.99)

    for run in spec["runs"]:
        path = Path(run["profile"])
        if not path.is_file() or path.stat().st_size == 0:
            raise RuntimeError(f"missing GEX ambient profile: {path}")
        values = {}
        with path.open("r", encoding="utf-8", errors="replace", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if not reader.fieldnames or "gene" not in reader.fieldnames or "ambient" not in reader.fieldnames:
                raise RuntimeError(f"invalid GEX profile schema: {path}")
            for row in reader:
                gene = str(row.get("gene", "")).strip()
                try:
                    value = float(row.get("ambient", "nan"))
                except (TypeError, ValueError):
                    value = math.nan
                if gene and math.isfinite(value) and value >= 0:
                    values[gene] = values.get(gene, 0.0) + value
        total = sum(values.values())
        if total <= 0:
            raise RuntimeError(f"GEX ambient profile has no positive mass: {path}")
        values = {gene: value / total for gene, value in values.items()}
        ordered = sorted(values.items(), key=lambda item: (-item[1], item[0]))
        cumulative = 0.0
        n_at = {}
        for rank, (gene, mass) in enumerate(ordered, start=1):
            cumulative += mass
            for threshold in thresholds:
                if threshold not in n_at and cumulative >= threshold:
                    n_at[threshold] = rank
            if rank <= int(spec.get("top_n", 200)):
                top_rows.append({
                    **run, "rank": rank, "gene": gene, "ambient_mass": mass,
                    "cumulative_mass": cumulative,
                    "gene_class": ("mitochondrial" if re.search(r"^MT-", gene, re.I)
                                   else "ribosomal_protein" if re.search(r"^RP[SL]", gene, re.I)
                                   else "other"),
                })
        all_masses = sorted(values.values())
        positive = [value for value in all_masses if value > 0]
        n_pos = len(positive)
        n_total = len(all_masses)
        gini = (sum((2 * i - n_total - 1) * value
                    for i, value in enumerate(all_masses, start=1)) /
                (n_total * sum(all_masses))) if n_total else math.nan
        entropy = -sum(value * math.log(value) for value in positive)
        row = {
            **run,
            "n_genes": len(values),
            "n_genes_nonzero": n_pos,
            "gini": gini,
            "shannon_effective_genes": math.exp(entropy),
            "simpson_effective_genes": 1.0 / sum(value * value for value in positive),
            "top10_mass": sum(value for _, value in ordered[:10]),
            "top25_mass": sum(value for _, value in ordered[:25]),
            "top50_mass": sum(value for _, value in ordered[:50]),
            "top100_mass": sum(value for _, value in ordered[:100]),
            "mitochondrial_mass": sum(value for gene, value in values.items()
                                      if re.search(r"^MT-", gene, re.I)),
            "ribosomal_protein_mass": sum(value for gene, value in values.items()
                                           if re.search(r"^RP[SL]", gene, re.I)),
        }
        for threshold in thresholds:
            row[f"n_genes_{int(threshold * 100)}pct"] = n_at.get(threshold, n_pos)
        metrics.append(row)
        key = run["run_key"]
        profiles[key] = values
        metadata[key] = run

    similarity = []
    keys = sorted(profiles)
    for i, key_a in enumerate(keys):
        for key_b in keys[i + 1:]:
            genes = sorted(set(profiles[key_a]) | set(profiles[key_b]))
            p = [profiles[key_a].get(gene, 0.0) for gene in genes]
            q = [profiles[key_b].get(gene, 0.0) for gene in genes]
            midpoint = [(a + b) / 2.0 for a, b in zip(p, q)]
            def kl_div(a, b):
                return sum(x * math.log(x / y) for x, y in zip(a, b)
                           if x > 0 and y > 0)
            denom = math.sqrt(sum(x * x for x in p) * sum(x * x for x in q))
            similarity.append({
                "run_a": key_a, "run_b": key_b,
                "library_a": metadata[key_a]["library"],
                "library_b": metadata[key_b]["library"],
                "condition_a": metadata[key_a]["condition"],
                "condition_b": metadata[key_b]["condition"],
                "assignment_source_a": metadata[key_a]["assignment_source"],
                "assignment_source_b": metadata[key_b]["assignment_source"],
                "cosine_similarity": (sum(a * b for a, b in zip(p, q)) / denom
                                      if denom else math.nan),
                "spearman_correlation": _gex_correlation(_gex_rank(p), _gex_rank(q)),
                "jensen_shannon_divergence": 0.5 * kl_div(p, midpoint) + 0.5 * kl_div(q, midpoint),
                "l1_distance": sum(abs(a - b) for a, b in zip(p, q)),
            })

    metric_fields = [
        "run_key", "library", "condition", "assignment_source", "profile",
        "cluster_source", "clusters", "n_genes", "n_genes_nonzero", "gini",
        "shannon_effective_genes", "simpson_effective_genes",
        "top10_mass", "top25_mass", "top50_mass", "top100_mass",
        "n_genes_50pct", "n_genes_80pct", "n_genes_90pct",
        "n_genes_95pct", "n_genes_99pct", "mitochondrial_mass",
        "ribosomal_protein_mass",
    ]
    _gex_atomic_tsv(out_dir / "gex_ambient_profile_metrics.tsv", metric_fields, metrics)
    _gex_atomic_tsv(
        out_dir / "gex_ambient_top_genes.tsv",
        ["run_key", "library", "condition", "assignment_source", "profile",
         "cluster_source", "clusters", "rank", "gene", "ambient_mass", "cumulative_mass",
         "gene_class"], top_rows)
    _gex_atomic_tsv(
        out_dir / "gex_ambient_profile_similarity.tsv",
        ["run_a", "run_b", "library_a", "library_b", "condition_a",
         "condition_b", "assignment_source_a", "assignment_source_b",
         "cosine_similarity", "spearman_correlation",
         "jensen_shannon_divergence", "l1_distance"], similarity)
    _ambient_write_json(out_dir / "gex_ambient_summary_contract.json", spec)
    print(f"GEX_AMBIENT summarized {len(metrics)} profiles under {out_dir}")


def generate_gex_ambient_script(lib_num, cond, assignment_source, args,
                                dep_job_ids=None, force=False):
    """Generate one isolated ambient_rna_gex run over filtered RNA counts."""
    abbrev = cond["abbrev"]
    analysis = _gex_safe_name(args.gex_analysis_name)
    source_prefix = get_contam_prefix(lib_num, abbrev, assignment_source)
    out_prefix = get_gex_ambient_prefix(lib_num, abbrev, assignment_source, args)
    out_dir = os.path.dirname(out_prefix)
    cluster_path = get_gex_cluster_path(lib_num, args)
    eligible_cluster_path = out_prefix + ".eligible_clusters.tsv"
    eligible_mex_dir = out_prefix + ".eligible_mtx"
    eligible_barcodes_path = eligible_mex_dir + "/barcodes.tsv.gz"
    eligible_matrix_path = eligible_mex_dir + "/matrix.mtx.gz"
    input_qc_path = out_prefix + ".gex_input_qc.tsv"
    intersection_spec_path = out_prefix + ".gex_input_spec.json"
    barcodes = get_filtered_barcodes(lib_num)
    features = get_filtered_features(lib_num)
    matrix = get_filtered_matrix(lib_num)
    binary = os.path.join(SOFTWARE_BIN, "tet_contam_estimate")
    source_suffix = "" if assignment_source == "demux" else f"_{assignment_source}"
    job_name = f"gex_{analysis[:18]}_{abbrev[:18]}{source_suffix}_l{lib_num}"
    log_dir = os.path.join(get_log_dir(), "GEX_AMBIENT", analysis, abbrev,
                           assignment_source)
    dep_line = _ambient_dependency_line(dep_job_ids)
    skip_regex = shlex.quote(args.gex_skip_genes_regex)
    round_flag = "" if args.gex_round_counts else "--noround"
    command = (
        f'"{binary}" -o "{out_prefix}" --interindividual '
        f'--barcodes "{eligible_barcodes_path}" --features "{features}" '
        f'--matrix "{eligible_matrix_path}" '
        f'--feature_type {shlex.quote(args.gex_feature_type)} '
        f'--clusts "{eligible_cluster_path}" --cellranger {round_flag} '
        f'--skip_genes_regex {skip_regex} -T "${{SLURM_CPUS_PER_TASK}}"')
    contract = {
        "schema_version": 1, "library": lib_num, "condition": abbrev,
        "assignment_source": assignment_source, "analysis_name": analysis,
        "source_prefix": source_prefix, "output_prefix": out_prefix,
        "barcodes": barcodes, "features": features, "matrix": matrix,
        "source_clusters": cluster_path,
        "eligible_clusters": eligible_cluster_path,
        "eligible_barcodes": eligible_barcodes_path,
        "eligible_matrix": eligible_matrix_path,
        "feature_type": args.gex_feature_type,
        "cluster_source": args.gex_cluster_source,
        "skip_genes_regex": args.gex_skip_genes_regex,
        "round_counts": bool(args.gex_round_counts),
        "scientific_scope": "gene-level ambient mass; not SNP opportunity or ASE",
    }
    os.makedirs(out_dir, exist_ok=True)
    planned_contract = out_prefix + ".gex_run_contract.planned.json"
    final_contract = out_prefix + ".gex_run_contract.json"
    _ambient_write_json(planned_contract, contract)
    _ambient_write_json(intersection_spec_path, {
        "schema_version": 1,
        "library": lib_num,
        "condition": abbrev,
        "assignment_source": assignment_source,
        "source_clusters": cluster_path,
        "mex_barcodes": barcodes,
        "mex_matrix": matrix,
        "contam_rate": source_prefix + ".contam_rate",
        "output": eligible_cluster_path,
        "subset_barcodes": eligible_barcodes_path,
        "subset_matrix": eligible_matrix_path,
        "qc_output": input_qc_path,
        "minimum_cells": args.gex_min_cluster_cells,
    })

    force_block = ""
    if force:
        force_block = f'''echo "Removing only prior GEX_AMBIENT products for this analysis namespace"
rm -f "{out_prefix}.gex_profile" \
      "{final_contract}" \
      "{out_prefix}_mtx/barcodes.tsv.gz" \
      "{out_prefix}_mtx/features.tsv.gz" \
      "{out_prefix}_mtx/matrix.mtx.gz"
'''
    script = f'''#!/bin/bash
#SBATCH --job-name={job_name}
#SBATCH --output={log_dir}/{job_name}_%j.out
#SBATCH --error={log_dir}/{job_name}_%j.err
#SBATCH --time={SLURM_TIME}
#SBATCH --cpus-per-task={args.gex_cpus}
#SBATCH --mem={args.gex_memory_gb}G
#SBATCH --partition={SLURM_PARTITION}
#SBATCH --nodes=1
#SBATCH --chdir={AGGREGATE_ROOT}
{dep_line}

set -euo pipefail
{module_block()}
mkdir -p "{out_dir}" "{log_dir}" "{out_prefix}_mtx" "{eligible_mex_dir}"

echo "================================================================"
echo "GEX_AMBIENT: ambient gene-expression profile"
echo "  Library: lib{lib_num}"
echo "  Condition: {abbrev}"
echo "  Assignment source: {assignment_source}"
echo "  Analysis namespace: {analysis}"
echo "  Source prefix: {source_prefix}"
echo "  Source clusters: {cluster_path}"
echo "  Eligible clusters: {eligible_cluster_path}"
echo "  Skip-gene regex: {args.gex_skip_genes_regex or '<empty; include all genes>'}"
echo "  Corrected counts: {'integer stochastic rounding' if args.gex_round_counts else 'fractional deterministic output'}"
echo "================================================================"

for path in "{binary}" "{source_prefix}.contam_rate" "{source_prefix}.contam_prof" \
            "{source_prefix}.assignments" "{source_prefix}.samples" \
            "{barcodes}" "{features}" "{matrix}" "{cluster_path}"; do
    [[ -s "$path" ]] || {{ echo "ERROR: required GEX_AMBIENT input missing or empty: $path" >&2; exit 1; }}
done

link_input() {{
    local source="$1" destination="$2"
    if [[ -e "$destination" && ! -L "$destination" ]]; then
        echo "ERROR: refusing to replace non-symlink input sidecar: $destination" >&2
        exit 1
    fi
    ln -sfn "$source" "$destination"
}}
link_input "{source_prefix}.contam_rate" "{out_prefix}.contam_rate"
link_input "{source_prefix}.contam_prof" "{out_prefix}.contam_prof"
link_input "{source_prefix}.assignments" "{out_prefix}.assignments"
link_input "{source_prefix}.samples" "{out_prefix}.samples"

python3 {_ambient_shell(os.path.abspath(__file__))} \
    --_gex-cluster-intersection-worker {_ambient_shell(intersection_spec_path)}
test -s "{eligible_cluster_path}"
test -s "{eligible_barcodes_path}"
test -s "{eligible_matrix_path}"
test -s "{input_qc_path}"

if [[ -s "{out_prefix}.gex_profile" ]] && \
   [[ -s "{out_prefix}_mtx/barcodes.tsv.gz" ]] && \
   [[ -s "{out_prefix}_mtx/features.tsv.gz" ]] && \
   [[ -s "{out_prefix}_mtx/matrix.mtx.gz" ]] && \
   [[ -s "{final_contract}" ]] && [[ "{1 if force else 0}" == "0" ]]; then
    echo "Complete GEX_AMBIENT outputs already exist; skipping"
    exit 0
fi
if [[ "{1 if force else 0}" == "0" ]] && \
   {{ [[ -e "{out_prefix}.gex_profile" ]] || [[ -e "{out_prefix}_mtx/matrix.mtx.gz" ]]; }}; then
    echo "ERROR: partial GEX_AMBIENT output exists; use --force or a new --gex-analysis-name" >&2
    exit 1
fi
{force_block}

echo "Running: {command}"
{command}

for path in "{out_prefix}.gex_profile" "{out_prefix}_mtx/barcodes.tsv.gz" \
            "{out_prefix}_mtx/features.tsv.gz" "{out_prefix}_mtx/matrix.mtx.gz"; do
    [[ -s "$path" ]] || {{ echo "ERROR: expected GEX_AMBIENT output missing or empty: $path" >&2; exit 1; }}
done
install -m 0644 "{planned_contract}" "{final_contract}.tmp.${{SLURM_JOB_ID:-$$}}"
mv -f "{final_contract}.tmp.${{SLURM_JOB_ID:-$$}}" "{final_contract}"
echo "GEX_AMBIENT complete: {out_prefix}"
'''
    os.makedirs(log_dir, exist_ok=True)
    script_path = os.path.join(
        get_script_dir(), f"{job_name}_{hashlib.sha1(out_prefix.encode()).hexdigest()[:8]}.sbatch")
    _ambient_write_text(script_path, script, executable=True)
    return script_path


def generate_gex_ambient_summary_script(args, runs, dep_job_ids=None):
    analysis = _gex_safe_name(args.gex_analysis_name)
    summary_dir = get_gex_ambient_summary_dir(args)
    os.makedirs(summary_dir, exist_ok=True)
    spec_path = os.path.join(summary_dir, "gex_ambient_summary_spec.json")
    _ambient_write_json(spec_path, {
        "schema_version": 1, "analysis_name": analysis,
        "summary_dir": summary_dir, "top_n": args.gex_top_n, "runs": runs,
        "interpretation_note": (
            "These tables describe inferred gene-level ambient mass. They do not "
            "measure SNP detectability, allelic imbalance, intronic/exonic origin, "
            "or prove that the corrected matrix is unbiased for DE/ASE."),
    })
    log_dir = os.path.join(get_log_dir(), "GEX_AMBIENT", analysis)
    dep_line = _ambient_dependency_line(dep_job_ids)
    script_path = os.path.join(get_script_dir(), f"gex_ambient_{analysis}_summary.sbatch")
    script = f'''#!/bin/bash
#SBATCH --job-name=gexsum_{analysis[:24]}
#SBATCH --output={log_dir}/gex_ambient_summary_%j.out
#SBATCH --error={log_dir}/gex_ambient_summary_%j.err
#SBATCH --time={SLURM_TIME}
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --partition={SLURM_PARTITION}
#SBATCH --nodes=1
{dep_line}

set -euo pipefail
module purge
module load miniforge/3
mkdir -p "{log_dir}" "{summary_dir}"
python3 {_ambient_shell(os.path.abspath(__file__))} --_gex-ambient-summary-worker {_ambient_shell(spec_path)}
test -s {_ambient_shell(os.path.join(summary_dir, 'gex_ambient_profile_metrics.tsv'))}
test -s {_ambient_shell(os.path.join(summary_dir, 'gex_ambient_top_genes.tsv'))}
echo "GEX_AMBIENT summary complete: {summary_dir}"
'''
    _ambient_write_text(script_path, script, executable=True)
    return script_path, spec_path


# =============================================================================
# AMBIENT_PLOTS: descriptive plots and deployed CellBouncer contam.R
# =============================================================================

def ambient_add_internal_worker_argument(parser):
    """Add hidden worker entry points to the orchestrator parser."""
    parser.add_argument(
        "--_ambient-plots-worker",
        metavar="SPEC_JSON",
        default=None,
        help=argparse.SUPPRESS,
    )
    parser.add_argument(
        "--_cleanup-results-worker",
        metavar="SPEC_JSON",
        default=None,
        help=argparse.SUPPRESS,
    )
    parser.add_argument(
        "--_publish-figure-shortcut",
        nargs=2,
        metavar=("ANALYSIS", "TARGET"),
        default=None,
        help=argparse.SUPPRESS,
    )


def ambient_maybe_run_internal_worker(args):
    """Run the hidden worker and return True when normal main must stop."""
    publish_request = getattr(args, "_publish_figure_shortcut", None)
    if publish_request:
        publish_figure_shortcut(*publish_request)
        return True
    cleanup_spec = getattr(args, "_cleanup_results_worker", None)
    if cleanup_spec:
        cleanup_results_worker(cleanup_spec)
        return True
    spec_path = getattr(args, "_ambient_plots_worker", None)
    if not spec_path:
        return False
    with open(spec_path, "r", encoding="utf-8") as handle:
        worker_kind = str(json.load(handle).get("worker_kind", "ambient_plots"))
    if worker_kind == "ambient_validation":
        ambient_validation_worker(spec_path)
    else:
        ambient_run_python_plot_worker(spec_path)
    return True


def _ambient_unique(items):
    return list(dict.fromkeys(items))


def _ambient_condition_slug(condition):
    slug = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(condition)).strip("._")
    if not slug:
        raise ValueError(f"condition has no safe filename characters: {condition!r}")
    return slug


def _ambient_shell(value):
    return shlex.quote(str(value))


def _ambient_dependency_line(job_ids):
    ids = [str(x).strip() for x in (job_ids or []) if str(x).strip()]
    return f"#SBATCH --dependency=afterok:{':'.join(ids)}" if ids else ""


def _ambient_write_text(path, text, executable=False):
    """Write an orchestrator-generated file and return its absolute path."""
    path = os.path.abspath(path)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    tmp = f"{path}.tmp.{os.getpid()}"
    with open(tmp, "w", encoding="utf-8") as handle:
        handle.write(text)
    os.replace(tmp, path)
    if executable:
        os.chmod(
            path,
            stat.S_IRWXU | stat.S_IRGRP | stat.S_IXGRP
            | stat.S_IROTH | stat.S_IXOTH,
        )
    return path


def _ambient_write_json(path, payload):
    text = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    return _ambient_write_text(path, text)


def ambient_condition_prefix_candidates(
    mapping_root,
    library,
    condition,
    assignment_source="demux",
    library_prefix="Tet_2025_Multiome-RNA_",
    demux_subdir="demux_nomito",
    contamination_subdir=AMBIENT_CONTAMINATION_SUBDIR,
    identity_ambient_candidate_set=None,
):
    """Return canonical then compatibility condition-result prefixes."""
    library = int(library)
    lib_name = f"lib{library}"
    base = os.path.join(
        os.path.abspath(mapping_root),
        f"{library_prefix}{library}",
        demux_subdir,
        contamination_subdir,
        str(condition),
    )
    if assignment_source == "reconciled":
        base = os.path.join(base, "reconciled")
    elif assignment_source in IDENTITY_AMBIENT_ARMS:
        base = os.path.join(
            base, IDENTITY_AMBIENT_DIRNAME,
            _identity_ambient_candidate_set(identity_ambient_candidate_set),
            assignment_source)
    elif assignment_source != "demux":
        raise ValueError(
            f"unknown contamination assignment source: {assignment_source}")
    return [
        os.path.join(base, f"{lib_name}_demuxed"),
        os.path.join(base, lib_name, f"{lib_name}_demuxed"),
    ]


def ambient_generate_plot_job_bundle(
    *,
    orchestrator_path,
    mapping_root,
    aggregate_root,
    plot_root=None,
    libraries,
    conditions=None,
    assignment_sources=("demux",),
    script_dir,
    log_dir,
    upstream_job_ids=None,
    partition="compute",
    library_prefix="Tet_2025_Multiome-RNA_",
    demux_subdir="demux_nomito",
    contamination_subdir=AMBIENT_CONTAMINATION_SUBDIR,
    deployed_contam_r=AMBIENT_DEPLOYED_CONTAM_R,
    r_array_parallelism=None,
    plot_formats=("pdf", "png"),
    reference_condition=None,
    identity_ambient_candidate_set=None,
):
    """Generate the R array, aggregate worker spec, and aggregate sbatch.

    Missing optional inputs are recorded as SKIP in per-task status files.
    Contract mismatches and contam.R execution/output-validation failures are
    recorded as FAIL and exit nonzero, so Slurm dependencies cannot silently
    accept a broken R array task. The aggregate output's
    ``data/contam_r_status.tsv`` retains the per-task audit trail.

    Returns a dictionary with ``r_manifest``, ``worker_spec``, ``r_sbatch``,
    and ``aggregate_sbatch``.  It does not submit anything.
    """
    libraries = _ambient_unique(int(x) for x in libraries)
    conditions = _ambient_unique(
        str(x) for x in (conditions or [AMBIENT_PLOT_DEFAULT_CONDITION])
    )
    assignment_sources = _ambient_unique(str(x) for x in assignment_sources)
    if not libraries:
        raise ValueError("AMBIENT_PLOTS needs at least one library")
    if not conditions:
        raise ValueError("AMBIENT_PLOTS needs at least one condition")
    valid_assignment_sources = set(CONTAM_ASSIGNMENT_SOURCES) | set(
        IDENTITY_AMBIENT_ARMS)
    if not assignment_sources or any(
            source not in valid_assignment_sources
            for source in assignment_sources):
        raise ValueError(
            "AMBIENT_PLOTS assignment sources must be demux, reconciled, or "
            "one of the reconciliation four-arm keys")
    four_arm_mode = all(
        source in IDENTITY_AMBIENT_ARMS for source in assignment_sources)
    identity_ambient_candidate_set = _identity_ambient_candidate_set(
        identity_ambient_candidate_set)
    if any(source in IDENTITY_AMBIENT_ARMS for source in assignment_sources) and not four_arm_mode:
        raise ValueError(
            "AMBIENT_PLOTS cannot mix legacy demux/reconciled series with "
            "reconciliation four-arm series in one report")
    if any(x < 1 for x in libraries):
        raise ValueError(f"invalid library list: {libraries}")
    if len(libraries) > 40:
        raise ValueError(
            f"AMBIENT_PLOTS supports at most 40 libraries, got {len(libraries)}"
        )
    if len(conditions) > 8:
        raise ValueError(
            f"AMBIENT_PLOTS supports at most 8 conditions, got {len(conditions)}"
        )
    if reference_condition is not None and reference_condition not in conditions:
        raise ValueError("ambient reference condition must be selected")
    if any("\t" in value or "\n" in value or "\r" in value for value in conditions):
        raise ValueError("condition names cannot contain tabs or newlines")
    series = []
    for condition in conditions:
        for assignment_source in assignment_sources:
            key = (condition if assignment_source == "demux" and
                   len(assignment_sources) == 1
                   else f"{condition}__{assignment_source}")
            arm_spec = IDENTITY_AMBIENT_ARMS.get(assignment_source, {})
            label = (
                arm_spec.get("short_label")
                if arm_spec else
                condition if len(assignment_sources) == 1 else
                f"{condition} [{assignment_source}]"
            )
            item = {
                "key": key,
                "label": label,
                "condition": condition,
                "assignment_source": assignment_source,
                "arm_key": assignment_source if arm_spec else "",
                "arm": arm_spec.get("arm", ""),
                "assignment_basis": arm_spec.get(
                    "assignment_basis", assignment_source),
                "roster_basis": arm_spec.get("roster_basis", "legacy"),
                "candidate_set": (
                    identity_ambient_candidate_set if arm_spec else "legacy"),
            }
            series.append(item)
    condition_slugs = [
        _ambient_condition_slug(item["key"]) for item in series]
    if len(condition_slugs) != len(set(condition_slugs)):
        raise ValueError("ambient plot series collapse to duplicate output slugs")
    reference_series = None
    if reference_condition is not None:
        candidates = [
            item["key"] for item in series
            if item["condition"] == reference_condition and
            item["assignment_source"] == "demux"
        ] or [
            item["key"] for item in series
            if item["condition"] == reference_condition
        ]
        reference_series = candidates[0]

    requested_formats = _ambient_unique(str(x).lower() for x in plot_formats)
    unsupported_formats = set(requested_formats) - {"pdf", "png"}
    if unsupported_formats:
        raise ValueError(
            "AMBIENT_PLOTS plot formats must be pdf and/or png, got "
            f"{sorted(unsupported_formats)}"
        )
    if not requested_formats:
        requested_formats = ["pdf"]
    validation_format = "pdf" if "pdf" in requested_formats else requested_formats[0]

    orchestrator_path = os.path.abspath(orchestrator_path)
    mapping_root = os.path.abspath(mapping_root)
    aggregate_root = os.path.abspath(aggregate_root)
    script_dir = os.path.abspath(script_dir)
    log_dir = os.path.abspath(log_dir)
    plot_root = os.path.abspath(
        plot_root or os.path.join(aggregate_root, "ambient_rna"))
    data_dir = os.path.join(plot_root, "data")
    status_dir = os.path.join(plot_root, "contam_r", "status")
    os.makedirs(script_dir, exist_ok=True)
    os.makedirs(log_dir, exist_ok=True)
    os.makedirs(data_dir, exist_ok=True)
    os.makedirs(status_dir, exist_ok=True)
    run_name = _ambient_condition_slug(os.path.basename(plot_root))
    run_token = (
        run_name[:36] + "_" + hashlib.sha1(plot_root.encode("utf-8")).hexdigest()[:8]
    )

    identity_arm_applicability = {}
    if four_arm_mode:
        for library in libraries:
            context_path = os.path.join(
                mapping_root, f"{library_prefix}{library}", demux_subdir,
                IDENTITY_AMBIENT_DIRNAME, identity_ambient_candidate_set,
                "plan", f"lib{library}.comparison_context.json",
            )
            try:
                with open(context_path, "r", encoding="utf-8") as handle:
                    context = json.load(handle)
            except (OSError, json.JSONDecodeError) as exc:
                raise ValueError(
                    f"AMBIENT_PLOTS cannot load current four-arm context "
                    f"for lib{library}: {context_path}: {exc}") from exc
            plan_fingerprint = str(
                context.get("plan_fingerprint", "")).strip()
            if (str(context.get("candidate_set", "")) !=
                    identity_ambient_candidate_set or not plan_fingerprint):
                raise ValueError(
                    f"AMBIENT_PLOTS current context candidate set/fingerprint "
                    f"is invalid for lib{library}: {context_path}")
            identity_arm_applicability[str(library)] = {
                "candidate_set": identity_ambient_candidate_set,
                "context_path": context_path,
                "plan_fingerprint": plan_fingerprint,
                "replacement_arm_eligible": bool(
                    context.get("replacement_arm_eligible", False)),
                "replacement_arm_skip_reason": str(
                    context.get("replacement_arm_skip_reason", "")),
            }

    manifest_path = os.path.join(data_dir, "contam_r_tasks.tsv")
    manifest_lines = [
        "library\tseries_key\tcondition\tassignment_source\tarm\tarm_key"
        "\tassignment_basis\troster_basis\tcandidate_set\tplan_fingerprint"
        "\tcondition_slug\tprimary_prefix"
        "\tcompat_prefix\tdemux_prefix\toutput_prefix"
    ]
    for library in libraries:
        for item in series:
            condition = item["condition"]
            assignment_source = item["assignment_source"]
            if (assignment_source == "reconciled_replacement" and
                    not identity_arm_applicability[str(library)][
                        "replacement_arm_eligible"]):
                # A stale D output from a previous plan must never be plotted.
                continue
            slug = _ambient_condition_slug(item["key"])
            prefixes = ambient_condition_prefix_candidates(
                mapping_root,
                library,
                condition,
                assignment_source=assignment_source,
                library_prefix=library_prefix,
                demux_subdir=demux_subdir,
                contamination_subdir=contamination_subdir,
                identity_ambient_candidate_set=(
                    identity_ambient_candidate_set),
            )
            output_prefix = os.path.join(
                plot_root, "contam_r", slug, f"lib{library}"
            )
            demux_prefix = os.path.join(
                mapping_root, f"{library_prefix}{library}", demux_subdir,
                f"lib{library}_demuxed")
            # Do not write empty tab fields: Bash treats tab as IFS whitespace
            # and would collapse adjacent empty columns in the array worker.
            values = [
                library, item["key"], condition, assignment_source,
                item.get("arm") or "NA", item.get("arm_key") or "NA",
                item.get("assignment_basis") or "NA",
                item.get("roster_basis") or "NA",
                (identity_ambient_candidate_set
                 if item.get("arm_key") else "NA"),
                (identity_arm_applicability[str(library)]["plan_fingerprint"]
                 if item.get("arm_key") else "NA"),
                slug, prefixes[0],
                prefixes[1], demux_prefix, output_prefix,
            ]
            manifest_lines.append("\t".join(str(x) for x in values))
    _ambient_write_text(manifest_path, "\n".join(manifest_lines) + "\n")

    worker_spec_path = os.path.join(data_dir, "ambient_plot_worker_spec.json")
    _ambient_write_json(
        worker_spec_path,
        {
            "schema_version": 3,
            "mapping_root": mapping_root,
            "aggregate_root": aggregate_root,
            "plot_root": plot_root,
            "libraries": libraries,
            "conditions": conditions,
            "assignment_sources": assignment_sources,
            "series": series,
            "library_prefix": library_prefix,
            "demux_subdir": demux_subdir,
            "contamination_subdir": contamination_subdir,
            "plot_formats": requested_formats,
            "reference_condition": reference_series,
            "comparison_mode": (
                IDENTITY_AMBIENT_SELECTOR if four_arm_mode else "legacy"),
            "identity_ambient_candidate_set": (
                identity_ambient_candidate_set),
            "identity_arm_applicability": identity_arm_applicability,
        },
    )

    r_script_path = os.path.join(
        script_dir, f"ambient_plots_{run_token}_contam_r.sbatch")
    n_tasks = len(manifest_lines) - 1
    if n_tasks < 1:
        raise ValueError("AMBIENT_PLOTS generated no applicable contam.R tasks")
    array_limit = ""
    if r_array_parallelism is not None:
        r_array_parallelism = max(1, min(int(r_array_parallelism), n_tasks))
        array_limit = f"%{r_array_parallelism}"
    r_dep_line = _ambient_dependency_line(upstream_job_ids)
    r_script = f'''#!/bin/bash
#SBATCH --job-name=ambR_{run_token}
#SBATCH --output={log_dir}/ambient_contam_r_%A_%a.out
#SBATCH --error={log_dir}/ambient_contam_r_%A_%a.err
#SBATCH --chdir={AGGREGATE_ROOT}
#SBATCH --time={SLURM_TIME}
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --partition={partition}
#SBATCH --nodes=1
#SBATCH --array=1-{n_tasks}{array_limit}
{r_dep_line}

set -euo pipefail
module purge
module load {AMBIENT_R_MODULE}
echo "Runtime host: $(hostname)"
echo "SLURM job: ${{SLURM_JOB_ID:-not-set}} array task ${{SLURM_ARRAY_TASK_ID:-not-set}}"
free -h || true

MANIFEST={_ambient_shell(manifest_path)}
CONTAM_R={_ambient_shell(os.path.abspath(deployed_contam_r))}
PLOT_ROOT={_ambient_shell(plot_root)}
STATUS_DIR={_ambient_shell(status_dir)}
FOUR_ARM_MODE={1 if four_arm_mode else 0}

row=$(awk -F '\t' -v n="$SLURM_ARRAY_TASK_ID" 'NR == n + 1 {{print; exit}}' "$MANIFEST")
if [ -z "$row" ]; then
    echo "No manifest row for array task $SLURM_ARRAY_TASK_ID" >&2
    exit 1
fi
IFS=$'\t' read -r library series_key condition assignment_source_key arm arm_key assignment_basis roster_basis candidate_set plan_fingerprint condition_slug primary_prefix compat_prefix demux_prefix output_prefix <<< "$row"
mkdir -p "$STATUS_DIR" "$(dirname "$output_prefix")" "$PLOT_ROOT/.staging"
status_file="$STATUS_DIR/${{condition_slug}}.lib${{library}}.tsv"

write_status() {{
    local state="$1"
    local reason="$2"
    local source="${{3:-}}"
    local tmp="${{status_file}}.tmp.${{SLURM_ARRAY_JOB_ID}}.${{SLURM_ARRAY_TASK_ID}}"
    printf 'library\tcondition\tassignment_source\tseries_key\tstatus\treason\tsource_prefix\n%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$library" "$condition" "$assignment_source_key" "$series_key" "$state" "$reason" "$source" > "$tmp"
    mv -f "$tmp" "$status_file"
}}

skip_task() {{
    echo "SKIP $condition lib$library: $1"
    write_status "SKIP" "$1" "${{2:-}}"
    exit 0
}}

fail_task() {{
    echo "FAIL $condition lib$library: $1" >&2
    write_status "FAIL" "$1" "${{2:-}}"
    exit 1
}}

source_prefix="$primary_prefix"
if [ ! -s "${{source_prefix}}.contam_rate" ] && [ -s "${{compat_prefix}}.contam_rate" ]; then
    source_prefix="$compat_prefix"
fi
[ -s "$CONTAM_R" ] || fail_task "missing_deployed_contam_R" "$source_prefix"
[ -s "${{source_prefix}}.contam_rate" ] || skip_task "missing_contam_rate" "$source_prefix"
[ -s "${{source_prefix}}.contam_prof" ] || skip_task "missing_individual_contam_prof" "$source_prefix"

assignment_suffix=""
assignment_source=""
if [ "$FOUR_ARM_MODE" -eq 1 ]; then
    # Four-arm estimators are run with --freeze_assignments.  contam.R must
    # therefore consume the exact estimator input, never an incidental legacy
    # fallback.  The Python worker separately verifies that the frozen final
    # assignment sidecar is byte-for-byte equivalent in identity calls.
    [ -s "${{source_prefix}}.assignments" ] || skip_task "missing_exact_arm_assignments" "$source_prefix"
    [ -s "${{source_prefix}}.decontam.assignments" ] || skip_task "missing_frozen_arm_assignments" "$source_prefix"
    [ -s "${{source_prefix}}.identity_ambient_arm.tsv" ] || skip_task "missing_identity_ambient_arm_contract" "$source_prefix"
    [ -s "${{source_prefix}}.cell_source_profile.tsv" ] || skip_task "missing_cell_source_profile" "$source_prefix"
    awk -F '\t' -v ec="$condition" -v ea="$arm" -v ek="$arm_key" \
        -v eb="$assignment_basis" -v er="$roster_basis" \
        -v ecs="$candidate_set" -v epf="$plan_fingerprint" '
        {{sub(/\r$/, "", $NF)}}
        NR == 1 {{
            for (i=1; i<=NF; i++) col[$i]=i
            header_ok=(col["condition"] && col["arm"] &&
                col["arm_key"] && col["assignment_basis"] &&
                col["roster_basis"] && col["candidate_set"] &&
                col["plan_fingerprint"] &&
                col["assignment_update_mode"] &&
                col["assignment_score_basis"])
            next
        }}
        {{
            rows++
            if ($(col["condition"]) != ec || $(col["arm"]) != ea ||
                $(col["arm_key"]) != ek ||
                $(col["assignment_basis"]) != eb ||
                $(col["roster_basis"]) != er ||
                $(col["candidate_set"]) != ecs ||
                $(col["plan_fingerprint"]) != epf ||
                $(col["assignment_update_mode"]) != "iterative_frozen" ||
                $(col["assignment_score_basis"]) != "original_demux_all_arms") bad=1
        }}
        END {{exit !(header_ok && rows == 1 && !bad)}}
    ' "${{source_prefix}}.identity_ambient_arm.tsv" || \
        fail_task "identity_ambient_arm_contract_mismatch" "$source_prefix"
    assignment_suffix=".assignments"
    assignment_source="${{source_prefix}}.assignments"
elif [ -s "${{source_prefix}}.decontam.assignments" ]; then
    assignment_suffix=".decontam.assignments"
    assignment_source="${{source_prefix}}.decontam.assignments"
elif [ -s "${{source_prefix}}.assignments" ]; then
    assignment_suffix=".assignments"
    assignment_source="${{source_prefix}}.assignments"
elif [ -s "${{demux_prefix}}.assignments" ]; then
    assignment_suffix=".assignments"
    assignment_source="${{demux_prefix}}.assignments"
else
    skip_task "missing_assignment_sidecar" "$source_prefix"
fi

awk 'NF {{seen=1; if (NF != 3) bad=1}} END {{exit !(seen && !bad)}}' \
    "${{source_prefix}}.contam_rate" || skip_task "contam_rate_not_three_columns" "$source_prefix"
awk 'NF {{seen=1; if (NF != 4) bad=1}} END {{exit !(seen && !bad)}}' \
    "$assignment_source" || skip_task "assignments_not_four_columns" "$source_prefix"
awk 'NF {{seen=1; if (NF < 2) bad=1}} END {{exit !(seen && !bad)}}' \
    "${{source_prefix}}.contam_prof" || skip_task "contam_prof_not_two_columns" "$source_prefix"

# Native species-only profiles are handled by the Python species plots.  The
# R script's ambient-vs-cell panel assumes individual profile keys.
if ! awk -F '\t' '$1 != "H" && $1 != "C" && $1 != "B" && $1 != "O" && $1 != "Hy" && $1 != "other_species" {{found=1}} END {{exit !found}}' \
    "${{source_prefix}}.contam_prof"; then
    skip_task "species_only_profile_not_valid_for_contam_R" "$source_prefix"
fi

if ! Rscript -e 'for (p in c("ggplot2", "ggsci", "cowplot")) if (!requireNamespace(p, quietly=TRUE)) stop(p); if (!capabilities("cairo")) stop("R lacks cairo capability")'; then
    fail_task "missing_R_plot_dependency" "$source_prefix"
fi

stage_dir="$PLOT_ROOT/.staging/${{condition_slug}}.lib${{library}}.${{SLURM_ARRAY_JOB_ID}}.${{SLURM_ARRAY_TASK_ID}}"
rm -rf "$stage_dir"
mkdir -p "$stage_dir"
trap 'rm -rf "$stage_dir"' EXIT
stage_prefix="$stage_dir/lib${{library}}"
ln -s "${{source_prefix}}.contam_rate" "${{stage_prefix}}.contam_rate"
ln -s "${{source_prefix}}.contam_prof" "${{stage_prefix}}.contam_prof"
ln -s "$assignment_source" "${{stage_prefix}}${{assignment_suffix}}"

if ! Rscript "$CONTAM_R" "$stage_prefix"; then
    fail_task "contam_R_nonzero_exit" "$source_prefix"
fi
for suffix in .contam.pdf .contam.png .contam.stats; do
    [ -s "${{stage_prefix}}${{suffix}}" ] || fail_task "missing_or_empty_output_${{suffix#.}}" "$source_prefix"
done

for suffix in .contam.pdf .contam.png .contam.stats; do
    publish_tmp="${{output_prefix}}${{suffix}}.tmp.${{SLURM_ARRAY_JOB_ID}}.${{SLURM_ARRAY_TASK_ID}}"
    install -m 0644 "${{stage_prefix}}${{suffix}}" "$publish_tmp"
    mv -f "$publish_tmp" "${{output_prefix}}${{suffix}}"
done
write_status "PASS" "validated" "$source_prefix"
echo "PASS $condition lib$library"
'''
    _ambient_write_text(r_script_path, r_script, executable=True)

    aggregate_script_path = os.path.join(
        script_dir, f"ambient_plots_{run_token}_aggregate.sbatch")
    python_module_lines = "\n".join(
        f"module load {module}" for module in AMBIENT_PYTHON_MODULES
    )
    aggregate_script = f'''#!/bin/bash
#SBATCH --job-name=ambP_{run_token}
#SBATCH --output={log_dir}/ambient_plots_%j.out
#SBATCH --error={log_dir}/ambient_plots_%j.err
#SBATCH --chdir={AGGREGATE_ROOT}
#SBATCH --time={SLURM_TIME}
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --partition={partition}
#SBATCH --nodes=1

set -euo pipefail
module purge
{python_module_lines}
echo "Runtime host: $(hostname)"
echo "SLURM job: ${{SLURM_JOB_ID:-not-set}}"
free -h || true
export MPLCONFIGDIR="${{SLURM_TMPDIR:-/tmp}}/matplotlib-${{SLURM_JOB_ID}}"
mkdir -p "$MPLCONFIGDIR"

python3 {_ambient_shell(orchestrator_path)} --_ambient-plots-worker {_ambient_shell(worker_spec_path)}
test -s {_ambient_shell(os.path.join(data_dir, 'ambient_rate_metrics.tsv'))}
test -s {_ambient_shell(os.path.join(plot_root, 'kde', f'kde_all_libraries.{validation_format}'))}
{chr(10).join('test -s ' + _ambient_shell(path) for path in ([
    os.path.join(data_dir, 'reconciliation_arm_contracts.tsv'),
    os.path.join(data_dir, 'reconciliation_rate_metrics.tsv'),
    os.path.join(data_dir, 'reconciliation_receiver_identity_metrics.tsv'),
    os.path.join(data_dir, 'reconciliation_normalized_profiles.tsv'),
    os.path.join(data_dir, 'reconciliation_exact_donor_burden.tsv'),
    os.path.join(data_dir, 'reconciliation_planned_contrasts.tsv'),
    os.path.join(data_dir, 'reconciliation_assignment_switch_summary.tsv'),
    os.path.join(data_dir, 'reconciliation_diagnostics.tsv'),
    os.path.join(plot_root, 'kde', f'reconciliation_arms_stratified.{validation_format}'),
] if four_arm_mode else []))}
test -s {_ambient_shell(os.path.join(data_dir, 'plot_manifest.tsv'))}
python3 {_ambient_shell(orchestrator_path)} --_publish-figure-shortcut \
  ambient_plots {_ambient_shell(plot_root)}
echo "AMBIENT_PLOTS aggregate outputs validated: {plot_root}"
'''
    _ambient_write_text(aggregate_script_path, aggregate_script, executable=True)

    return {
        "r_manifest": manifest_path,
        "worker_spec": worker_spec_path,
        "r_sbatch": r_script_path,
        "aggregate_sbatch": aggregate_script_path,
        "plot_root": plot_root,
        "n_r_tasks": n_tasks,
    }


def _ambient_submit_sbatch(script_path, dependency_job_ids=None):
    command = ["sbatch", "--parsable", f"--chdir={AGGREGATE_ROOT}"]
    deps = [str(x).strip() for x in (dependency_job_ids or []) if str(x).strip()]
    if deps:
        command.append(f"--dependency=afterok:{':'.join(deps)}")
    command.append(script_path)
    last_detail = "no response"
    for attempt in (1, 2):
        result = subprocess.run(command, check=False, capture_output=True, text=True)
        response = result.stdout.strip().split(";", 1)[0]
        if result.returncode == 0 and response.isdigit():
            return response
        last_detail = result.stderr.strip() or result.stdout.strip() or (
            f"exit {result.returncode}")
        if attempt == 1:
            print(f"  WARNING: sbatch attempt 1 failed for {script_path}: {last_detail}")
            print("  Retrying once...")
    raise RuntimeError(f"sbatch did not submit {script_path}: {last_detail}")


def ambient_submit_plot_job_bundle(bundle, upstream_job_ids=None, submit=False):
    """Submit the generated two-job chain, or describe it in dry-run mode."""
    if not submit:
        return {
            "r_job_id": None,
            "aggregate_job_id": None,
            "r_sbatch": bundle["r_sbatch"],
            "aggregate_sbatch": bundle["aggregate_sbatch"],
        }
    # Upstream dependencies can be supplied here instead of embedding them in
    # the generated R script.  Do not supply the same IDs both ways.
    r_job_id = _ambient_submit_sbatch(bundle["r_sbatch"], upstream_job_ids)
    aggregate_job_id = _ambient_submit_sbatch(
        bundle["aggregate_sbatch"], [r_job_id]
    )
    return {"r_job_id": r_job_id, "aggregate_job_id": aggregate_job_id}


def _ambient_atomic_tsv(path, fieldnames, rows):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(f".{path.name}.tmp.{os.getpid()}")
    with tmp.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
    os.replace(tmp, path)


def _ambient_load_profile_details(path):
    """Return raw and normalized ambient-profile values without hiding scale.

    The plotting code historically normalized column two of ``.contam_prof``
    to one.  That remains the display convention, but the four-arm diagnostic
    report also retains the pre-normalization sum and optional concentration
    column so an apparent donor-fraction change cannot be mistaken for an
    increase in absolute ambient burden.
    """
    profile = {}
    concentrations = {}
    path = Path(path)
    if not path.is_file() or path.stat().st_size == 0:
        return {
            "raw": {}, "normalized": {}, "concentration": {},
            "raw_total": 0.0,
        }
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            try:
                value = float(parts[1])
            except ValueError:
                continue
            if math.isfinite(value) and value >= 0:
                profile[parts[0]] = value
                if len(parts) >= 3:
                    try:
                        concentration = float(parts[2])
                    except ValueError:
                        concentration = math.nan
                    if math.isfinite(concentration):
                        concentrations[parts[0]] = concentration
    total = sum(profile.values())
    normalized = (
        {key: value / total for key, value in profile.items()}
        if total > 0 else {}
    )
    return {
        "raw": profile,
        "normalized": normalized,
        "concentration": concentrations,
        "raw_total": total,
    }


def _ambient_load_two_column_profile(path):
    """Load the normalized second column used by legacy profile figures."""
    return _ambient_load_profile_details(path)["normalized"]


def _ambient_species_code(identity):
    identity = str(identity)
    if identity.startswith("Chinobo"):
        return "Hy"
    if identity.startswith("Congo"):
        return "B"
    if identity.startswith("JOS"):
        return "O"
    if identity.startswith("C") and identity[1:].isdigit():
        return "C"
    if identity.startswith("H") or identity.startswith("KOLF"):
        return "H"
    if identity in {"H", "C", "B", "O", "Hy"}:
        return identity
    return None


def _ambient_profile_to_species(profile):
    out = {}
    for identity, value in (profile or {}).items():
        species = _ambient_species_code(identity)
        if species == "Hy":
            out["C"] = out.get("C", 0.0) + 0.5 * value
            out["B"] = out.get("B", 0.0) + 0.5 * value
        elif species:
            out[species] = out.get(species, 0.0) + value
    total = sum(out.values())
    return {key: value / total for key, value in out.items()} if total > 0 else {}


def _ambient_load_assignments(prefix, demux_fallback=None):
    candidates = [
        f"{prefix}.decontam.assignments",
        f"{prefix}.assignments",
    ]
    if demux_fallback:
        candidates.append(str(demux_fallback))
    for candidate in candidates:
        path = Path(candidate)
        if not path.is_file() or path.stat().st_size == 0:
            continue
        assignments = {}
        with path.open("r", encoding="utf-8", errors="replace") as handle:
            for line in handle:
                parts = line.rstrip("\n").split("\t")
                if len(parts) >= 2:
                    assignments[parts[0]] = parts[1]
        if assignments:
            return assignments, str(path)
    return {}, ""


def _ambient_load_assignment_file(path):
    """Load one exact assignment file without decontamination fallbacks."""
    path = Path(path)
    assignments = {}
    if not path.is_file() or path.stat().st_size == 0:
        return assignments
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2:
                assignments[parts[0]] = parts[1]
    return assignments


def _ambient_load_assignment_scores(path):
    """Load barcode -> finite column-4 LLR from one exact assignment file."""
    path = Path(path)
    scores = {}
    malformed = 0
    if not path.is_file() or path.stat().st_size == 0:
        return scores, 0
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4 or not parts[0]:
                malformed += 1
                continue
            try:
                score = float(parts[3])
            except ValueError:
                malformed += 1
                continue
            if not math.isfinite(score):
                malformed += 1
                continue
            scores[parts[0]] = score
    return scores, malformed


def _ambient_assignment_profile(assignments, equal_identities=False):
    weights = {}
    identities = (
        sorted(set(assignments.values())) if equal_identities
        else list(assignments.values())
    )
    for identity in identities:
        parts = [part for part in str(identity).split("+") if part]
        if not parts:
            continue
        weight = 1.0 / len(parts)
        for part in parts:
            weights[part] = weights.get(part, 0.0) + weight
    total = sum(weights.values())
    return {key: value / total for key, value in weights.items()} if total > 0 else {}


def _ambient_profile_l1(left, right):
    if not left or not right:
        return math.nan
    keys = set(left) | set(right)
    return sum(abs(left.get(key, 0.0) - right.get(key, 0.0)) for key in keys)


def _ambient_expand_species_profile(species_profile, model_b):
    by_species = {}
    for identity, value in (model_b or {}).items():
        species = _ambient_species_code(identity)
        if species == "Hy":
            for target in ("C", "B"):
                by_species.setdefault(target, {})[identity] = (
                    by_species.setdefault(target, {}).get(identity, 0.0)
                    + 0.5 * value
                )
        elif species:
            by_species.setdefault(species, {})[identity] = (
                by_species.setdefault(species, {}).get(identity, 0.0) + value
            )
    for species, values in by_species.items():
        total = sum(values.values())
        if total > 0:
            by_species[species] = {key: value / total for key, value in values.items()}
    expanded = {}
    for species, fraction in (species_profile or {}).items():
        if species == "Hy":
            targets = (("C", 0.5), ("B", 0.5))
        else:
            targets = ((species, 1.0),)
        for target, split in targets:
            for identity, within_weight in by_species.get(target, {}).items():
                expanded[identity] = expanded.get(identity, 0.0) + fraction * split * within_weight
    total = sum(expanded.values())
    return {key: value / total for key, value in expanded.items()} if total > 0 else {}


def _ambient_resolve_prefix(spec, library, condition,
                            assignment_source="demux"):
    candidates = ambient_condition_prefix_candidates(
        spec["mapping_root"],
        library,
        condition,
        assignment_source=assignment_source,
        library_prefix=spec.get("library_prefix", "Tet_2025_Multiome-RNA_"),
        demux_subdir=spec.get("demux_subdir", "demux_nomito"),
        contamination_subdir=spec.get(
            "contamination_subdir", AMBIENT_CONTAMINATION_SUBDIR
        ),
        identity_ambient_candidate_set=spec.get(
            "identity_ambient_candidate_set",
            IDENTITY_AMBIENT_CANDIDATE_SET),
    )
    for prefix in candidates:
        rate_path = Path(f"{prefix}.contam_rate")
        if rate_path.is_file() and rate_path.stat().st_size > 0:
            return prefix
    return candidates[0]


def _ambient_read_gate_summary(prefix):
    audit_candidates = [
        Path(f"{prefix}.geometry_gate_audit.tsv"),
        Path(f"{prefix}.contam_diagnostics.tsv"),
    ]
    for path in audit_candidates:
        if not path.is_file() or path.stat().st_size == 0:
            continue
        with path.open("r", encoding="utf-8", errors="replace", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            fields = set(reader.fieldnames or [])
            if "geometry_gate_triggered" not in fields:
                continue
            rows = list(reader)
        if not rows:
            continue
        triggered = 0
        movements = []
        triggered_movements = []
        for row in rows:
            trigger_text = str(row.get("geometry_gate_triggered", "")).strip().lower()
            is_triggered = trigger_text in {"1", "1.0", "true"}
            triggered += int(is_triggered)
            def finite(name):
                try:
                    value = float(row.get(name, "nan"))
                except (TypeError, ValueError):
                    return math.nan
                return value if math.isfinite(value) else math.nan
            base = finite("base_c_selected")
            selected = finite("selected_c")
            if not math.isfinite(selected):
                endpoint = str(row.get("selected_endpoint", "")).lower()
                if endpoint == "fallback" or is_triggered:
                    selected = finite("fallback_c_selected")
                else:
                    selected = base
            if math.isfinite(base) and math.isfinite(selected):
                movement = selected - base
                movements.append(movement)
                if is_triggered:
                    triggered_movements.append(movement)
        return {
            "n_rows": len(rows),
            "n_triggered": triggered,
            "fraction_triggered": triggered / len(rows),
            "mean_selected_minus_base_c": (
                sum(movements) / len(movements) if movements else math.nan
            ),
            "mean_selected_minus_base_c_triggered": (
                sum(triggered_movements) / len(triggered_movements)
                if triggered_movements else math.nan
            ),
            "audit_path": str(path),
        }
    return None


def ambient_run_python_plot_worker(spec_path):
    """Create descriptive, non-ranking aggregate contamination plots."""
    # Heavy dependencies remain worker-only by design.
    # Matplotlib otherwise tries to write below $HOME on some compute nodes.
    if not os.environ.get("MPLCONFIGDIR"):
        matplotlib_config = Path(
            os.environ.get("SLURM_TMPDIR", "/tmp")
        ) / f"matplotlib-{os.getpid()}"
        matplotlib_config.mkdir(parents=True, exist_ok=True)
        os.environ["MPLCONFIGDIR"] = str(matplotlib_config)
    import numpy as np
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from scipy.stats import gaussian_kde

    spec_path = os.path.abspath(spec_path)
    with open(spec_path, "r", encoding="utf-8") as handle:
        spec = json.load(handle)
    libraries = _ambient_unique(int(x) for x in spec["libraries"])
    series = spec.get("series") or [
        {
            "key": str(condition),
            "label": str(condition),
            "condition": str(condition),
            "assignment_source": "demux",
        }
        for condition in spec["conditions"]
    ]
    conditions = _ambient_unique(str(item["key"]) for item in series)
    series_by_key = {str(item["key"]): item for item in series}
    four_arm_mode = (
        spec.get("comparison_mode") == IDENTITY_AMBIENT_SELECTOR)
    identity_ambient_candidate_set = _identity_ambient_candidate_set(
        spec.get("identity_ambient_candidate_set"))
    if four_arm_mode:
        unknown_arms = sorted({
            str(item.get("assignment_source", "")) for item in series
        } - set(IDENTITY_AMBIENT_ARMS))
        if unknown_arms:
            raise ValueError(
                "four-arm worker spec contains non-arm series: "
                + ", ".join(unknown_arms))
        for item in series:
            arm_key = str(item["assignment_source"])
            arm_spec = IDENTITY_AMBIENT_ARMS[arm_key]
            for field in ("arm", "assignment_basis", "roster_basis"):
                observed = str(item.get(field, ""))
                expected = str(arm_spec[field])
                if observed != expected:
                    raise ValueError(
                        f"four-arm worker series {arm_key} has {field}="
                        f"{observed!r}, expected {expected!r}")
            if str(item.get("candidate_set", "")) != (
                    identity_ambient_candidate_set):
                raise ValueError(
                    f"four-arm worker series {arm_key} has candidate_set="
                    f"{item.get('candidate_set')!r}, expected "
                    f"{identity_ambient_candidate_set!r}")
    if not libraries or not conditions:
        raise ValueError("ambient plot worker spec has no libraries or conditions")
    max_plotted_series = 32 if four_arm_mode else 16
    if len(libraries) > 40 or len(conditions) > max_plotted_series:
        raise ValueError(
            "ambient plot worker supports at most 40 libraries and "
            f"{max_plotted_series} plotted series in this comparison mode"
        )
    slugs = [_ambient_condition_slug(condition) for condition in conditions]
    if len(slugs) != len(set(slugs)):
        raise ValueError("ambient plot worker condition slugs are not unique")

    plot_root = Path(spec.get("plot_root") or (
        Path(spec["aggregate_root"]) / AMBIENT_PLOT_DIRNAME
    ))
    data_dir = plot_root / "data"
    overview_dir = plot_root / "overview"
    kde_dir = plot_root / "kde"
    per_condition_dir = kde_dir / "by_condition"
    profile_dir = plot_root / "profiles"
    for directory in (data_dir, overview_dir, kde_dir, per_condition_dir, profile_dir):
        directory.mkdir(parents=True, exist_ok=True)

    requested_formats = {
        str(x).lower() for x in spec.get("plot_formats", ["pdf", "png"])
    }
    formats = [x for x in ("pdf", "png") if x in requested_formats]
    if not formats:
        formats = ["pdf"]
    generated = []
    plot_decisions = []

    plt.rcParams.update({
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.titleweight": "bold",
        "figure.facecolor": "white",
        "savefig.facecolor": "white",
    })

    def save_figure(fig, stem, dpi=160):
        stem = Path(stem)
        stem.parent.mkdir(parents=True, exist_ok=True)
        for suffix in formats:
            target = stem.with_suffix(f".{suffix}")
            tmp = target.with_name(f".{target.stem}.tmp.{os.getpid()}.{suffix}")
            kwargs = {"bbox_inches": "tight"}
            if suffix == "png":
                kwargs["dpi"] = dpi
            fig.savefig(tmp, format=suffix, **kwargs)
            if not tmp.is_file() or tmp.stat().st_size == 0:
                raise RuntimeError(f"empty plot output: {tmp}")
            os.replace(tmp, target)
            generated.append(str(target))
        plt.close(fig)
        plot_decisions.append({
            "figure": str(stem),
            "status": "GENERATED",
            "reason": "",
        })

    def skip_figure(stem, reason):
        """Record a deliberate omission and remove stale generated copies."""
        stem = Path(stem)
        for suffix in ("pdf", "png"):
            target = stem.with_suffix(f".{suffix}")
            if target.is_file() or target.is_symlink():
                target.unlink()
        plot_decisions.append({
            "figure": str(stem),
            "status": "SKIPPED",
            "reason": str(reason),
        })

    short_map = {
        "IND_CK_RF_SX0_RFREE_PFIT": "IR0",
        "IND_CK_RF_SX0_GATED_RFREE_PFIT": "SX0 gated",
        "IND_CK_RF_SX025_RFREE_PFIT": "IR25",
        "IND_CK_RF_SX025_GATED_RFREE_PFIT": "IRG25",
        "IND_CK_RF_SX050_RFREE_PFIT": "IR50",
        "IND_CK_RF_SX075_RFREE_PFIT": "IR75",
        "IND_CK_RF_SX1_RFREE_PFIT": "IR100",
    }
    biological_conditions = _ambient_unique(
        str(item["condition"]) for item in series
    )
    assignment_sources = _ambient_unique(
        str(item["assignment_source"]) for item in series
    )
    source_labels = {
        "demux": "Original demux identities",
        "reconciled": "Identity-reconciled identities",
    }
    source_short_labels = {
        "demux": "demux",
        "reconciled": "reconciled",
    }
    biological_labels = {
        condition: short_map.get(condition, condition[:24])
        for condition in biological_conditions
    }
    condition_labels = {}
    for condition in conditions:
        item = series_by_key[condition]
        base = biological_labels[item["condition"]]
        if four_arm_mode:
            label = IDENTITY_AMBIENT_ARMS[
                item["assignment_source"]]["short_label"]
            if len(biological_conditions) > 1:
                label = f"{base} — {label}"
        elif len(biological_conditions) == 1 and len(assignment_sources) > 1:
            label = source_labels.get(
                item["assignment_source"], item["assignment_source"])
        elif len(assignment_sources) > 1:
            source = source_short_labels.get(
                item["assignment_source"], item["assignment_source"])
            label = f"{base} — {source}"
        else:
            label = base
        condition_labels[condition] = label

    stable_colors = [
        "#4477AA", "#EE6677", "#228833", "#CCBB44",
        "#66CCEE", "#AA3377", "#BBBBBB", "#000000",
    ]
    biological_colors = {
        condition: stable_colors[index % len(stable_colors)]
        for index, condition in enumerate(biological_conditions)
    }
    source_colors = {
        "demux": "#4477AA",
        "reconciled": "#EE7733",
    }
    condition_colors = {}
    condition_linestyles = {}
    for condition in conditions:
        item = series_by_key[condition]
        if four_arm_mode:
            condition_colors[condition] = source_colors.get(
                item["assignment_basis"], "#666666")
        elif len(biological_conditions) == 1 and len(assignment_sources) > 1:
            condition_colors[condition] = source_colors.get(
                item["assignment_source"], "#666666")
        else:
            condition_colors[condition] = biological_colors[item["condition"]]
        if four_arm_mode:
            condition_linestyles[condition] = {
                "original": "-", "augmented": "--", "replacement": ":",
            }.get(item["roster_basis"], "-")
        else:
            condition_linestyles[condition] = (
                "--" if item["assignment_source"] == "reconciled" else "-"
            )

    rates = {}
    rate_nonfinite_barcodes = {}
    rate_malformed_rows = {}
    ses = {}
    assignments = {}
    input_assignments = {}
    final_assignments = {}
    input_assignment_scores = {}
    input_assignment_score_malformed = {}
    individual_profiles = {}
    profile_details = {}
    species_profiles = {}
    prefixes = {}
    arm_contracts = {}
    comparison_contexts = {}
    scrutiny_by_library = {}
    source_profile_paths = {}
    input_rows = []
    metrics_rows = []
    gate_rows = []

    def read_tsv_records(path):
        path = Path(path)
        if not path.is_file() or path.stat().st_size == 0:
            return [], []
        with path.open(
                "r", encoding="utf-8", errors="replace", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            return list(reader), list(reader.fieldnames or [])

    def read_identity_roster(path):
        path = Path(path)
        if not path.is_file() or path.stat().st_size == 0:
            return []
        values = []
        with path.open("r", encoding="utf-8", errors="replace") as handle:
            for raw in handle:
                value = raw.strip()
                if value and not value.startswith("#"):
                    values.append(value)
        return values

    identity_arm_applicability = spec.get(
        "identity_arm_applicability", {})
    if four_arm_mode:
        for library in libraries:
            applicability = identity_arm_applicability.get(str(library), {})
            context_path = str(applicability.get("context_path", ""))
            expected_fingerprint = str(
                applicability.get("plan_fingerprint", ""))
            try:
                with open(context_path, "r", encoding="utf-8") as handle:
                    current_context = json.load(handle)
            except (OSError, json.JSONDecodeError) as exc:
                raise RuntimeError(
                    f"four-arm worker cannot load current context for "
                    f"lib{library}: {context_path}: {exc}") from exc
            if (str(current_context.get("candidate_set", "")) !=
                    identity_ambient_candidate_set or
                    str(current_context.get("plan_fingerprint", "")) !=
                    expected_fingerprint or not expected_fingerprint or
                    bool(current_context.get("replacement_arm_eligible", False)) !=
                    bool(applicability.get("replacement_arm_eligible", False))):
                raise RuntimeError(
                    f"four-arm context changed after plot bundle generation "
                    f"for lib{library}: {context_path}")
            comparison_contexts[library] = current_context

    for library in libraries:
        demux_dir = (
            Path(spec["mapping_root"])
            / f"{spec.get('library_prefix', 'Tet_2025_Multiome-RNA_')}{library}"
            / spec.get("demux_subdir", "demux_nomito")
        )
        demux_assignment = demux_dir / f"lib{library}_demuxed.assignments"
        for condition in conditions:
            item = series_by_key[condition]
            if (four_arm_mode and
                    item["assignment_source"] == "reconciled_replacement" and
                    not identity_arm_applicability[str(library)][
                        "replacement_arm_eligible"]):
                input_rows.append({
                    "library": library,
                    "condition": condition,
                    "status": "NOT_APPLICABLE",
                    "reason": identity_arm_applicability[str(library)].get(
                        "replacement_arm_skip_reason") or
                        "replacement_arm_not_eligible_in_current_plan",
                    "prefix": "",
                    "n_cells": 0,
                    "assignment_path": "",
                    "profile_kind": "none",
                })
                continue
            prefix = _ambient_resolve_prefix(
                spec, library, item["condition"],
                item["assignment_source"])
            prefixes[(condition, library)] = prefix
            rate_path = Path(f"{prefix}.contam_rate")
            if not rate_path.is_file() or rate_path.stat().st_size == 0:
                input_rows.append({
                    "library": library,
                    "condition": condition,
                    "status": "MISSING",
                    "reason": "missing_contam_rate",
                    "prefix": prefix,
                    "n_cells": 0,
                    "assignment_path": "",
                    "profile_kind": "none",
                })
                continue
            barcode_values = {}
            barcode_ses = {}
            nonfinite_rate_barcodes = set()
            malformed_rate_rows = 0
            with rate_path.open("r", encoding="utf-8", errors="replace") as handle:
                for line in handle:
                    parts = line.rstrip("\n").split("\t")
                    if len(parts) < 2 or not parts[0]:
                        malformed_rate_rows += 1
                        continue
                    barcode = parts[0]
                    try:
                        value = float(parts[1])
                    except ValueError:
                        nonfinite_rate_barcodes.add(barcode)
                        continue
                    if not math.isfinite(value):
                        nonfinite_rate_barcodes.add(barcode)
                        continue
                    barcode_values[barcode] = value
                    try:
                        se_value = float(parts[2]) if len(parts) >= 3 else math.nan
                    except ValueError:
                        se_value = math.nan
                    barcode_ses[barcode] = (
                        se_value if math.isfinite(se_value) else math.nan)
            rate_nonfinite_barcodes[(condition, library)] = (
                nonfinite_rate_barcodes)
            rate_malformed_rows[(condition, library)] = malformed_rate_rows
            if not barcode_values:
                input_rows.append({
                    "library": library,
                    "condition": condition,
                    "status": "MISSING",
                    "reason": "no_finite_rates",
                    "prefix": prefix,
                    "n_cells": 0,
                    "assignment_path": "",
                    "profile_kind": "none",
                })
                continue
            rates[(condition, library)] = barcode_values
            ses[(condition, library)] = barcode_ses
            if four_arm_mode:
                input_map = _ambient_load_assignment_file(
                    f"{prefix}.assignments")
                final_map = _ambient_load_assignment_file(
                    f"{prefix}.decontam.assignments")
                assignment_map = input_map
                assignment_path = f"{prefix}.assignments" if input_map else ""
                input_assignments[(condition, library)] = input_map
                final_assignments[(condition, library)] = final_map
                score_map, malformed_scores = _ambient_load_assignment_scores(
                    f"{prefix}.assignments")
                input_assignment_scores[(condition, library)] = score_map
                input_assignment_score_malformed[
                    (condition, library)] = malformed_scores
            else:
                assignment_map, assignment_path = _ambient_load_assignments(
                    prefix, demux_fallback=demux_assignment
                )
            assignments[(condition, library)] = assignment_map

            details = _ambient_load_profile_details(f"{prefix}.contam_prof")
            profile_details[(condition, library)] = details
            contam_profile = details["normalized"]
            native_species_profile = _ambient_load_two_column_profile(f"{prefix}.species_prof")
            species_keys = {"H", "C", "B", "O", "Hy"}
            profile_kind = "none"
            if contam_profile:
                if set(contam_profile).issubset(species_keys | {"other_species"}):
                    species_profiles[(condition, library)] = {
                        key: value for key, value in contam_profile.items()
                        if key in species_keys
                    }
                    profile_kind = "species_in_contam_prof"
                else:
                    individual_profiles[(condition, library)] = contam_profile
                    profile_kind = "individual"
            if native_species_profile:
                species_profiles[(condition, library)] = native_species_profile
                profile_kind = (
                    "individual_and_species" if profile_kind == "individual" else "species"
                )

            if four_arm_mode:
                contract_path = Path(f"{prefix}.identity_ambient_arm.tsv")
                contract_records, contract_fields = read_tsv_records(
                    contract_path)
                required_contract = {
                    "library", "condition", "arm", "arm_key",
                    "assignment_basis", "roster_basis", "candidate_set",
                    "plan_fingerprint", "assignment_path",
                    "receiver_lines", "ambient_candidates", "scrutiny_cells",
                    "context_path", "assignment_update_mode",
                    "assignment_score_basis",
                }
                if len(contract_records) != 1 or not required_contract.issubset(
                        contract_fields):
                    raise RuntimeError(
                        f"invalid identity ambient arm contract: {contract_path}")
                contract = contract_records[0]
                expected_metadata = {
                    "condition": str(item["condition"]),
                    "arm": str(item["arm"]),
                    "arm_key": str(item["arm_key"]),
                    "assignment_basis": str(item["assignment_basis"]),
                    "roster_basis": str(item["roster_basis"]),
                    "candidate_set": identity_ambient_candidate_set,
                    "assignment_update_mode": "iterative_frozen",
                    "assignment_score_basis": "original_demux_all_arms",
                }
                mismatches = {
                    field: (contract.get(field, ""), expected)
                    for field, expected in expected_metadata.items()
                    if str(contract.get(field, "")) != expected
                }
                contract_library = str(contract.get("library", "")).strip()
                if contract_library.lower().startswith("lib"):
                    contract_library = contract_library[3:]
                if contract_library != str(library):
                    mismatches["library"] = (
                        contract.get("library", ""), str(library))
                if mismatches:
                    raise RuntimeError(
                        f"identity ambient arm contract metadata mismatch for "
                        f"{contract_path}: {mismatches}")
                if not input_map or not final_map:
                    raise RuntimeError(
                        f"four-arm output lacks exact/frozen assignments: {prefix}")
                arm_contracts[(condition, library)] = {
                    **contract, "contract_path": str(contract_path),
                }

                source_path = Path(f"{prefix}.cell_source_profile.tsv")
                with source_path.open(
                        "r", encoding="utf-8", errors="replace",
                        newline="") as handle:
                    source_reader = csv.DictReader(handle, delimiter="\t")
                    source_fields = list(source_reader.fieldnames or [])
                    has_source_row = next(source_reader, None) is not None
                required_source = {
                    "barcode", "source_label", "scoring_profile_mass",
                    "scoring_profile_mass_sum", "scoring_profile_status",
                }
                if not has_source_row or not required_source.issubset(source_fields):
                    raise RuntimeError(
                        f"invalid cell source-profile table: "
                        f"{prefix}.cell_source_profile.tsv")
                # Do not materialize this cell-by-source ledger.  Production
                # tables can contain millions of rows; the four-arm analysis
                # streams one arm at a time below and retains only O(cells +
                # sources * populations) aggregates.
                source_profile_paths[(condition, library)] = source_path

                context_path = contract.get("context_path", "")
                if library not in comparison_contexts:
                    try:
                        with open(context_path, "r", encoding="utf-8") as handle:
                            comparison_contexts[library] = json.load(handle)
                    except (OSError, json.JSONDecodeError) as exc:
                        raise RuntimeError(
                            f"invalid four-arm comparison context {context_path}: "
                            f"{exc}") from exc
                comparison_context = comparison_contexts[library]
                expected_context_path = identity_arm_applicability[
                    str(library)]["context_path"]
                if (os.path.abspath(str(contract.get("context_path", ""))) !=
                        os.path.abspath(expected_context_path) or
                        str(comparison_context.get("candidate_set", "")) !=
                        identity_ambient_candidate_set or
                        str(comparison_context.get("plan_fingerprint", "")) !=
                        str(contract.get("plan_fingerprint", "")) or
                        not str(contract.get("plan_fingerprint", ""))):
                    raise RuntimeError(
                        f"four-arm candidate-set/plan fingerprint mismatch "
                        f"for {contract_path}")
                scrutiny_path = contract.get("scrutiny_cells", "")
                if library not in scrutiny_by_library:
                    scrutiny_records, scrutiny_fields = read_tsv_records(
                        scrutiny_path)
                    required_scrutiny = {
                        "barcode", "demux_identity", "reconciled_identity",
                        "changed", "scrutinized", "stratum", "transition",
                    }
                    if (not scrutiny_records or
                            not required_scrutiny.issubset(scrutiny_fields)):
                        raise RuntimeError(
                            f"invalid four-arm scrutiny table: {scrutiny_path}")
                    scrutiny_by_library[library] = {
                        row["barcode"]: row for row in scrutiny_records
                    }

            values = np.asarray(list(barcode_values.values()), dtype=float)
            se_values = np.asarray(list(barcode_ses.values()), dtype=float)
            q25, q75 = np.percentile(values, [25, 75])
            finite_se = np.isfinite(se_values)
            positive_se = finite_se & (se_values > 0)
            metrics_rows.append({
                "condition": condition,
                "library": library,
                "n_cells": len(values),
                "mean": float(np.mean(values)),
                "median": float(np.median(values)),
                "q25": float(q25),
                "q75": float(q75),
                "iqr": float(q75 - q25),
                "frac_boundary": float(np.mean((values < 0.02) | (values > 0.90))),
                "frac_plausible": float(np.mean((values >= 0.05) & (values <= 0.50))),
                "frac_se_finite": float(np.mean(finite_se)),
                "frac_se_nonzero": (
                    float(np.mean(se_values[finite_se] > 0)) if np.any(finite_se) else math.nan
                ),
                "frac_se_zero": (
                    float(np.mean(se_values[finite_se] == 0)) if np.any(finite_se) else math.nan
                ),
                "median_se_nonzero": (
                    float(np.median(se_values[positive_se])) if np.any(positive_se) else math.nan
                ),
                "mean_se_nonzero": (
                    float(np.mean(se_values[positive_se])) if np.any(positive_se) else math.nan
                ),
            })
            input_rows.append({
                "library": library,
                "condition": condition,
                "status": "READY",
                "reason": "",
                "prefix": prefix,
                "n_cells": len(values),
                "assignment_path": assignment_path,
                "profile_kind": profile_kind,
            })
            gate = _ambient_read_gate_summary(prefix)
            if gate:
                gate_rows.append({"condition": condition, "library": library, **gate})

    if not rates:
        _ambient_atomic_tsv(
            data_dir / "input_status.tsv",
            ["library", "condition", "status", "reason", "prefix", "n_cells", "assignment_path", "profile_kind"],
            input_rows,
        )
        raise RuntimeError("AMBIENT_PLOTS found no finite contamination-rate data")

    metric_fields = [
        "condition", "library", "n_cells", "mean", "median", "q25", "q75", "iqr",
        "frac_boundary", "frac_plausible", "frac_se_finite",
        "frac_se_nonzero", "frac_se_zero", "median_se_nonzero",
        "mean_se_nonzero",
    ]
    _ambient_atomic_tsv(data_dir / "ambient_rate_metrics.tsv", metric_fields, metrics_rows)
    _ambient_atomic_tsv(
        data_dir / "input_status.tsv",
        ["library", "condition", "status", "reason", "prefix", "n_cells", "assignment_path", "profile_kind"],
        input_rows,
    )

    # The reconciliation comparison is a planned factorial experiment, not a
    # generic multi-condition ranking.  Build all summaries on fixed barcode
    # sets and retain both normalized profile composition and absolute
    # contamination burden.  In particular, a larger normalized Chinobo slice
    # is not called a larger Chinobo burden unless c * scoring_profile_mass also
    # rises on the same cells.
    four_arm_common = {}
    four_arm_hard_failures = []
    reconciliation_rate_rows = []
    receiver_identity_rows = []
    profile_component_rows = []
    donor_burden_rows = []
    contrast_rows = []
    switch_summary_rows = []
    contract_rows = []
    diagnostic_rows = []

    def finite_number(value):
        try:
            result = float(value)
        except (TypeError, ValueError):
            return math.nan
        return result if math.isfinite(result) else math.nan

    def diagnostic(library, biological_condition, item, check, status,
                   observed="", expected="", details=""):
        diagnostic_rows.append({
            "library": library,
            "condition": biological_condition,
            "arm": item.get("arm", "NA") if item else "NA",
            "arm_key": item.get("arm_key", "NA") if item else "NA",
            "check": check,
            "status": status,
            "observed": observed,
            "expected": expected,
            "details": details,
        })

    def describe_values(values):
        array = np.asarray(list(values), dtype=float)
        array = array[np.isfinite(array)]
        if not len(array):
            return {
                "n_cells": 0, "mean": math.nan, "median": math.nan,
                "q25": math.nan, "q75": math.nan, "iqr": math.nan,
            }
        q25, q75 = np.percentile(array, [25, 75])
        return {
            "n_cells": int(len(array)),
            "mean": float(np.mean(array)),
            "median": float(np.median(array)),
            "q25": float(q25),
            "q75": float(q75),
            "iqr": float(q75 - q25),
        }

    if four_arm_mode:
        arm_series = {
            (str(item["condition"]), str(item["assignment_source"])):
                str(item["key"])
            for item in series
        }
        mandatory_arms = IDENTITY_AMBIENT_ARM_ORDER[:3]
        contrast_plan = (
            ("roster_effect_demux", "demux_original", "demux_augmented"),
            ("assignment_effect_augmented", "demux_augmented",
             "reconciled_augmented"),
            ("replacement_sensitivity", "reconciled_augmented",
             "reconciled_replacement"),
            ("combined_production_change", "demux_original",
             "reconciled_augmented"),
        )
        for biological_condition in biological_conditions:
            missing_series = [
                arm_key for arm_key in mandatory_arms
                if (biological_condition, arm_key) not in arm_series
            ]
            if missing_series:
                raise ValueError(
                    f"four-arm plot spec for {biological_condition} lacks "
                    f"mandatory arms {missing_series}")
            for library in libraries:
                context = comparison_contexts.get(library, {})
                replacement_eligible = bool(
                    context.get("replacement_arm_eligible", False))
                required_arm_keys = list(mandatory_arms)
                if replacement_eligible:
                    required_arm_keys.append("reconciled_replacement")
                missing_outputs = [
                    arm_key for arm_key in required_arm_keys
                    if (arm_series.get((biological_condition, arm_key)), library)
                    not in rates
                ]
                if missing_outputs:
                    failure = (
                        f"lib{library} {biological_condition} missing required "
                        f"four-arm outputs: {','.join(missing_outputs)}")
                    four_arm_hard_failures.append(failure)
                    diagnostic(
                        library, biological_condition, None,
                        "mandatory_arm_outputs", "FAIL",
                        observed=",".join(missing_outputs), expected="A,B,C"
                        + (",D" if replacement_eligible else ""),
                        details=failure,
                    )
                    continue

                available_arm_keys = [
                    arm_key for arm_key in IDENTITY_AMBIENT_ARM_ORDER
                    if (arm_series.get((biological_condition, arm_key)), library)
                    in rates
                ]
                baseline_finite_key = arm_series[
                    (biological_condition, "demux_original")]
                baseline_finite_barcodes = set(
                    rates[(baseline_finite_key, library)])
                for finite_arm_key in available_arm_keys:
                    finite_key = arm_series[
                        (biological_condition, finite_arm_key)]
                    finite_item = series_by_key[finite_key]
                    finite_barcodes = set(rates[(finite_key, library)])
                    assignment_universe = set(input_assignments.get(
                        (finite_key, library), {}))
                    nonfinite_barcodes = rate_nonfinite_barcodes.get(
                        (finite_key, library), set())
                    seen_barcodes = finite_barcodes | nonfinite_barcodes
                    malformed = rate_malformed_rows.get(
                        (finite_key, library), 0)
                    missing_finite = sorted(
                        assignment_universe - finite_barcodes)
                    extra_rate_rows = sorted(
                        seen_barcodes - assignment_universe)
                    cross_arm_delta = sorted(
                        finite_barcodes ^ baseline_finite_barcodes)
                    finite_ok = (
                        not missing_finite and not extra_rate_rows and
                        not nonfinite_barcodes and malformed == 0 and
                        not cross_arm_delta and
                        finite_barcodes == assignment_universe)
                    diagnostic(
                        library, biological_condition, finite_item,
                        "finite_rate_barcode_contract",
                        "PASS" if finite_ok else "FAIL",
                        f"finite={len(finite_barcodes)};"
                        f"assignments={len(assignment_universe)};"
                        f"missing={len(missing_finite)};"
                        f"extra={len(extra_rate_rows)};"
                        f"nonfinite={len(nonfinite_barcodes)};"
                        f"malformed={malformed};"
                        f"cross_arm_delta={len(cross_arm_delta)}",
                        "finite=assignments;missing=extra=nonfinite="
                        "malformed=cross_arm_delta=0",
                        "examples=" + ",".join(_ambient_unique(
                            missing_finite[:3] + extra_rate_rows[:3] +
                            sorted(nonfinite_barcodes)[:3] +
                            cross_arm_delta[:3])),
                    )
                    if not finite_ok:
                        four_arm_hard_failures.append(
                            f"lib{library} {finite_key} finite-rate barcode "
                            "universe differs from assignments or another arm")
                common = None
                for arm_key in available_arm_keys:
                    key = arm_series[(biological_condition, arm_key)]
                    cells = set(rates[(key, library)])
                    common = cells if common is None else common & cells
                common = set(common or ())
                four_arm_common[(biological_condition, library)] = common
                if not common:
                    failure = (
                        f"lib{library} {biological_condition} has no barcode "
                        "intersection across available reconciliation arms")
                    four_arm_hard_failures.append(failure)
                diagnostic(
                    library, biological_condition, None,
                    "fixed_common_barcode_set",
                    "PASS" if common else "FAIL", len(common), ">0",
                    "intersection across " + ",".join(available_arm_keys),
                )

                scrutiny = scrutiny_by_library.get(library, {})
                population_sets = {"full_library": set(common)}
                for stratum in (
                        "changed_target", "scrutinized_other", "background"):
                    population_sets[stratum] = {
                        barcode for barcode in common
                        if scrutiny.get(barcode, {}).get("stratum") == stratum
                    }

                # The four arms may change identity and assignment type, but
                # must carry the same original demux evidence weight (column-4
                # LLR) for every barcode.  Otherwise B->C would confound an
                # identity change with a change in statistical weighting.
                baseline_key = arm_series[
                    (biological_condition, "demux_original")]
                baseline_scores = input_assignment_scores.get(
                    (baseline_key, library), {})
                for score_arm_key in available_arm_keys:
                    score_key = arm_series[(biological_condition, score_arm_key)]
                    score_item = series_by_key[score_key]
                    score_map = input_assignment_scores.get(
                        (score_key, library), {})
                    malformed = input_assignment_score_malformed.get(
                        (score_key, library), 0)
                    score_barcodes = set(baseline_scores) | set(score_map)
                    mismatched_scores = sorted(
                        barcode for barcode in score_barcodes
                        if barcode not in baseline_scores or
                        barcode not in score_map or
                        not math.isclose(
                            baseline_scores[barcode], score_map[barcode],
                            rel_tol=1e-12, abs_tol=1e-12)
                    )
                    score_ok = (
                        malformed == 0 and bool(baseline_scores) and
                        not mismatched_scores)
                    diagnostic(
                        library, biological_condition, score_item,
                        "assignment_column4_score_fixed_across_arms",
                        "PASS" if score_ok else "FAIL",
                        f"mismatches={len(mismatched_scores)};"
                        f"malformed={malformed};rows={len(score_map)}",
                        f"mismatches=0;malformed=0;rows={len(baseline_scores)}",
                        ("example_barcodes=" +
                         ",".join(mismatched_scores[:5]))
                        if mismatched_scores else "",
                    )
                    if not score_ok:
                        four_arm_hard_failures.append(
                            f"lib{library} {score_key} assignment column-4 "
                            "scores differ from Arm A")

                for arm_key in available_arm_keys:
                    key = arm_series[(biological_condition, arm_key)]
                    item = series_by_key[key]
                    prefix = prefixes[(key, library)]
                    contract = arm_contracts[(key, library)]
                    input_map = input_assignments[(key, library)]
                    final_map = final_assignments[(key, library)]
                    profile = profile_details[(key, library)]
                    source_path = source_profile_paths[(key, library)]
                    expected_sources = set(read_identity_roster(
                        contract.get("ambient_candidates", "")))

                    profile_status = (
                        "PASS" if profile["raw_total"] > 0 else "FAIL")
                    diagnostic(
                        library, biological_condition, item,
                        "contam_profile_raw_total", profile_status,
                        profile["raw_total"], ">0",
                        "normalized fractions are reported separately",
                    )
                    for source_label in sorted(profile["raw"]):
                        profile_component_rows.append({
                            "library": library,
                            "condition": biological_condition,
                            "series_key": key,
                            "arm": item["arm"],
                            "arm_key": arm_key,
                            "assignment_basis": item["assignment_basis"],
                            "roster_basis": item["roster_basis"],
                            "source_label": source_label,
                            "raw_profile_value": profile["raw"][source_label],
                            "normalized_profile_fraction": profile[
                                "normalized"].get(source_label, 0.0),
                            "profile_concentration": profile[
                                "concentration"].get(source_label, math.nan),
                            "raw_profile_total": profile["raw_total"],
                        })

                    all_assignment_barcodes = sorted(
                        set(input_map) | set(final_map))
                    n_final_switches = 0
                    n_planned_mismatches = 0
                    n_missing_planned = 0
                    for barcode in all_assignment_barcodes:
                        input_identity = input_map.get(barcode, "")
                        final_identity = final_map.get(barcode, "")
                        final_changed = bool(
                            input_identity and final_identity and
                            _canonical_identity(input_identity) !=
                            _canonical_identity(final_identity))
                        n_final_switches += int(final_changed)
                        scrutiny_row = scrutiny.get(barcode, {})
                        planned_identity = str(scrutiny_row.get(
                            "demux_identity" if item["assignment_basis"] == "demux"
                            else "reconciled_identity", "")).strip()
                        if planned_identity == "NA":
                            planned_identity = ""
                        if not planned_identity:
                            n_missing_planned += 1
                        planned_mismatch = bool(
                            input_identity and planned_identity and
                            _canonical_identity(input_identity) !=
                            _canonical_identity(planned_identity))
                        n_planned_mismatches += int(planned_mismatch)
                    switch_summary_rows.append({
                        "library": library,
                        "condition": biological_condition,
                        "series_key": key,
                        "arm": item["arm"],
                        "arm_key": arm_key,
                        "assignment_basis": item["assignment_basis"],
                        "n_input": len(input_map),
                        "n_final": len(final_map),
                        "n_union": len(all_assignment_barcodes),
                        "n_input_vs_planned_mismatch": n_planned_mismatches,
                        "n_missing_planned_identity": n_missing_planned,
                        "n_final_vs_input_switch": n_final_switches,
                        "frozen_assignment_invariant": (
                            "PASS" if n_final_switches == 0 and
                            set(input_map) == set(final_map) else "FAIL"),
                    })
                    frozen_ok = (
                        n_final_switches == 0 and set(input_map) == set(final_map))
                    if not frozen_ok:
                        four_arm_hard_failures.append(
                            f"lib{library} {key} violated frozen assignments")
                    diagnostic(
                        library, biological_condition, item,
                        "frozen_assignment_invariant",
                        "PASS" if frozen_ok else "FAIL",
                        f"switches={n_final_switches};input={len(input_map)};"
                        f"final={len(final_map)}",
                        "switches=0; identical barcode sets",
                    )
                    diagnostic(
                        library, biological_condition, item,
                        "assignment_matches_planned_basis",
                        "PASS" if n_planned_mismatches == 0 else "FAIL",
                        n_planned_mismatches, 0,
                        f"missing planned identities={n_missing_planned}",
                    )
                    if n_planned_mismatches:
                        four_arm_hard_failures.append(
                            f"lib{library} {key} assignment input does not "
                            "match its declared demux/reconciled basis")

                    # Stream the exact cell-specific scoring-profile ledger.
                    # Never retain its cell-by-source rows: four arms across 40
                    # libraries can otherwise exceed the worker's 32-GiB job.
                    observed_sources = set()
                    mass_sum_by_barcode = {}
                    source_mass_sums_by_population = {
                        population: {} for population in population_sets
                    }
                    source_burden_sums_by_population = {
                        population: {} for population in population_sets
                    }
                    invalid_mass_rows = 0
                    with source_path.open(
                            "r", encoding="utf-8", errors="replace",
                            newline="") as handle:
                        reader = csv.DictReader(handle, delimiter="\t")
                        for row in reader:
                            barcode = str(row.get("barcode", "")).strip()
                            source_label = str(row.get(
                                "source_label", "")).strip() or "NA"
                            mass = finite_number(row.get(
                                "scoring_profile_mass"))
                            if (not barcode or not math.isfinite(mass)
                                    or mass < 0):
                                invalid_mass_rows += 1
                                continue
                            if source_label != "other_species":
                                observed_sources.add(source_label)
                            if barcode not in common:
                                continue
                            mass_sum_by_barcode[barcode] = (
                                mass_sum_by_barcode.get(barcode, 0.0) + mass)
                            populations = ["full_library"]
                            stratum = scrutiny.get(
                                barcode, {}).get("stratum", "")
                            if stratum in population_sets:
                                populations.append(stratum)
                            c_value = rates[(key, library)][barcode]
                            for population in populations:
                                mass_sums = source_mass_sums_by_population[
                                    population]
                                burden_sums = source_burden_sums_by_population[
                                    population]
                                mass_sums[source_label] = (
                                    mass_sums.get(source_label, 0.0) + mass)
                                burden_sums[source_label] = (
                                    burden_sums.get(source_label, 0.0)
                                    + c_value * mass)

                    missing_sources = sorted(expected_sources - observed_sources)
                    unexpected_sources = sorted(observed_sources - expected_sources)
                    roster_status = (
                        "PASS" if not missing_sources and not unexpected_sources
                        else "FAIL")
                    if roster_status == "FAIL":
                        four_arm_hard_failures.append(
                            f"lib{library} {key} source roster mismatch")
                    diagnostic(
                        library, biological_condition, item,
                        "ambient_candidate_roster", roster_status,
                        observed=len(observed_sources),
                        expected=len(expected_sources),
                        details=(
                            f"missing={','.join(missing_sources) or 'none'};"
                            f"unexpected={','.join(unexpected_sources) or 'none'}"),
                    )
                    contract_rows.append({
                        "library": library,
                        "condition": biological_condition,
                        "series_key": key,
                        "arm": item["arm"],
                        "arm_key": arm_key,
                        "assignment_basis": item["assignment_basis"],
                        "roster_basis": item["roster_basis"],
                        "candidate_set": contract.get("candidate_set", ""),
                        "plan_fingerprint": contract.get(
                            "plan_fingerprint", ""),
                        "assignment_path": contract.get("assignment_path", ""),
                        "receiver_lines": contract.get("receiver_lines", ""),
                        "ambient_candidates": contract.get(
                            "ambient_candidates", ""),
                        "scrutiny_cells": contract.get("scrutiny_cells", ""),
                        "context_path": contract.get("context_path", ""),
                        "assignment_update_mode": contract.get(
                            "assignment_update_mode", ""),
                        "assignment_score_basis": contract.get(
                            "assignment_score_basis", ""),
                        "n_expected_sources": len(expected_sources),
                        "n_observed_sources": len(observed_sources),
                        "missing_sources": ",".join(missing_sources) or "NA",
                        "unexpected_sources": (
                            ",".join(unexpected_sources) or "NA"),
                        "profile_raw_total": profile["raw_total"],
                    })
                    mass_errors = [
                        abs(total - 1.0)
                        for total in mass_sum_by_barcode.values()
                    ]
                    max_mass_error = max(mass_errors) if mass_errors else math.inf
                    mass_status = (
                        "PASS" if invalid_mass_rows == 0 and
                        max_mass_error <= 1e-6 else "FAIL")
                    if mass_status == "FAIL":
                        four_arm_hard_failures.append(
                            f"lib{library} {key} invalid scoring-profile simplex")
                    diagnostic(
                        library, biological_condition, item,
                        "cell_scoring_profile_simplex", mass_status,
                        f"max_abs_error={max_mass_error:.6g};"
                        f"invalid_rows={invalid_mass_rows}",
                        "max_abs_error<=1e-6;invalid_rows=0",
                    )
                    coverage = (
                        len(common & set(mass_sum_by_barcode)) / len(common)
                        if common else 0.0)
                    coverage_status = "PASS" if coverage == 1.0 else "FAIL"
                    if coverage_status == "FAIL":
                        four_arm_hard_failures.append(
                            f"lib{library} {key} incomplete source-profile coverage")
                    diagnostic(
                        library, biological_condition, item,
                        "source_profile_common_barcode_coverage",
                        coverage_status, coverage, 1.0,
                    )

                    # Fixed-population rate metrics and exact donor burden.
                    for population, barcodes in population_sets.items():
                        summary = describe_values(
                            rates[(key, library)][barcode]
                            for barcode in barcodes)
                        reconciliation_rate_rows.append({
                            "library": library,
                            "condition": biological_condition,
                            "series_key": key,
                            "arm": item["arm"],
                            "arm_key": arm_key,
                            "assignment_basis": item["assignment_basis"],
                            "roster_basis": item["roster_basis"],
                            "population_kind": (
                                "full_library" if population == "full_library"
                                else "scrutiny_stratum"),
                            "population": population,
                            **summary,
                        })
                        if not barcodes:
                            continue
                        source_mass_sums = (
                            source_mass_sums_by_population[population])
                        source_burden_sums = (
                            source_burden_sums_by_population[population])
                        mean_c = summary["mean"]
                        summed_mean_burden = (
                            sum(source_burden_sums.values()) / len(barcodes))
                        burden_error = abs(summed_mean_burden - mean_c)
                        burden_status = (
                            "PASS" if burden_error <= 1e-6 else "FAIL")
                        diagnostic(
                            library, biological_condition, item,
                            f"exact_donor_burden_sum:{population}",
                            burden_status, summed_mean_burden, mean_c,
                            f"absolute_error={burden_error:.6g}",
                        )
                        if burden_status == "FAIL":
                            four_arm_hard_failures.append(
                                f"lib{library} {key} exact donor burden does "
                                f"not sum to mean c for {population}")
                        for source_label in sorted(source_burden_sums):
                            mean_burden = (
                                source_burden_sums[source_label] / len(barcodes))
                            donor_burden_rows.append({
                                "library": library,
                                "condition": biological_condition,
                                "series_key": key,
                                "arm": item["arm"],
                                "arm_key": arm_key,
                                "assignment_basis": item["assignment_basis"],
                                "roster_basis": item["roster_basis"],
                                "population": population,
                                "source_label": source_label,
                                "n_cells": len(barcodes),
                                "mean_scoring_profile_mass": (
                                    source_mass_sums[source_label] /
                                    len(barcodes)),
                                "mean_exact_contam_burden": mean_burden,
                                "fraction_of_mean_contam_burden": (
                                    mean_burden / mean_c
                                    if math.isfinite(mean_c) and mean_c > 0
                                    else math.nan),
                                "mean_total_contam_rate": mean_c,
                            })

                    # Per-receiver metrics make an unchanged identity such as
                    # Chinobo-mCherry directly auditable across A/B/C/D.
                    receiver_groups = {}
                    for barcode in common:
                        identity = input_map.get(barcode, "NA") or "NA"
                        receiver_groups.setdefault(identity, set()).add(barcode)
                    for identity, barcodes in sorted(receiver_groups.items()):
                        summary = describe_values(
                            rates[(key, library)][barcode]
                            for barcode in barcodes)
                        receiver_identity_rows.append({
                            "library": library,
                            "condition": biological_condition,
                            "series_key": key,
                            "arm": item["arm"],
                            "arm_key": arm_key,
                            "assignment_basis": item["assignment_basis"],
                            "roster_basis": item["roster_basis"],
                            "receiver_identity": identity,
                            **summary,
                        })

                    profile_diag_records, _ = read_tsv_records(
                        f"{prefix}.profile_fit_diagnostics.tsv")
                    if profile_diag_records:
                        profile_diag = profile_diag_records[0]
                        successful = finite_number(profile_diag.get(
                            "multistart_successful_starts"))
                        attempted = str(profile_diag.get(
                            "multistart_attempted", "")).strip()
                        nonunique = str(profile_diag.get(
                            "profile_nonunique_flag", "")).strip()
                        status = (
                            "PASS" if attempted in {"1", "true", "True"}
                            and math.isfinite(successful) and successful >= 1
                            else "WARN")
                        diagnostic(
                            library, biological_condition, item,
                            "profile_multistart", status,
                            f"attempted={attempted};successful={successful};"
                            f"nonunique={nonunique}",
                            "attempted=1;successful>=1",
                        )
                    else:
                        diagnostic(
                            library, biological_condition, item,
                            "profile_multistart", "WARN", "missing", "present",
                        )
                    model_records, model_fields = read_tsv_records(
                        f"{prefix}.model_fit.tsv")
                    model_metrics = {
                        row.get("metric", ""): row.get("value", "")
                        for row in model_records
                    } if {"metric", "value"}.issubset(model_fields) else {}
                    final_ll = finite_number(model_metrics.get("final_loglik"))
                    ll_valid = str(model_metrics.get(
                        "final_loglik_valid", "")).lower() == "true"
                    diagnostic(
                        library, biological_condition, item,
                        "model_final_log_likelihood",
                        "PASS" if ll_valid and math.isfinite(final_ll) else "WARN",
                        final_ll, "finite with final_loglik_valid=true",
                    )

                # Planned pairwise contrasts use each pair's own fixed barcode
                # intersection, then report the same predeclared strata.
                for contrast_name, left_arm, right_arm in contrast_plan:
                    left_key = arm_series.get((biological_condition, left_arm))
                    right_key = arm_series.get((biological_condition, right_arm))
                    if ((left_key, library) not in rates or
                            (right_key, library) not in rates):
                        continue
                    left_rates = rates[(left_key, library)]
                    right_rates = rates[(right_key, library)]
                    pair_common = set(left_rates) & set(right_rates)
                    contrast_populations = {
                        "full_library": pair_common,
                        **{
                            stratum: {
                                barcode for barcode in pair_common
                                if scrutiny.get(barcode, {}).get("stratum")
                                == stratum
                            }
                            for stratum in (
                                "changed_target", "scrutinized_other",
                                "background")
                        },
                    }
                    for population, barcodes in contrast_populations.items():
                        deltas = np.asarray([
                            right_rates[barcode] - left_rates[barcode]
                            for barcode in sorted(barcodes)
                        ], dtype=float)
                        contrast_rows.append({
                            "library": library,
                            "condition": biological_condition,
                            "contrast": contrast_name,
                            "left_arm": IDENTITY_AMBIENT_ARMS[left_arm]["arm"],
                            "right_arm": IDENTITY_AMBIENT_ARMS[right_arm]["arm"],
                            "left_arm_key": left_arm,
                            "right_arm_key": right_arm,
                            "population": population,
                            "n_common": len(deltas),
                            "mean_delta": (
                                float(np.mean(deltas)) if len(deltas)
                                else math.nan),
                            "median_delta": (
                                float(np.median(deltas)) if len(deltas)
                                else math.nan),
                            "median_absolute_delta": (
                                float(np.median(np.abs(deltas))) if len(deltas)
                                else math.nan),
                            "fraction_lower_in_right_arm": (
                                float(np.mean(deltas < 0)) if len(deltas)
                                else math.nan),
                        })

        def iter_assignment_switch_audit():
            """Regenerate cell rows lazily to keep 40-library memory bounded."""
            for biological_condition in biological_conditions:
                for library in libraries:
                    scrutiny = scrutiny_by_library.get(library, {})
                    for arm_key in IDENTITY_AMBIENT_ARM_ORDER:
                        key = arm_series.get((biological_condition, arm_key))
                        input_map = input_assignments.get((key, library), {})
                        final_map = final_assignments.get((key, library), {})
                        if not input_map and not final_map:
                            continue
                        item = series_by_key[key]
                        for barcode in sorted(set(input_map) | set(final_map)):
                            input_identity = input_map.get(barcode, "")
                            final_identity = final_map.get(barcode, "")
                            scrutiny_row = scrutiny.get(barcode, {})
                            planned_identity = str(scrutiny_row.get(
                                "demux_identity"
                                if item["assignment_basis"] == "demux"
                                else "reconciled_identity", "")).strip()
                            if planned_identity == "NA":
                                planned_identity = ""
                            planned_mismatch = bool(
                                input_identity and planned_identity and
                                _canonical_identity(input_identity) !=
                                _canonical_identity(planned_identity))
                            final_changed = bool(
                                input_identity and final_identity and
                                _canonical_identity(input_identity) !=
                                _canonical_identity(final_identity))
                            yield {
                                "library": library,
                                "condition": biological_condition,
                                "series_key": key,
                                "arm": item["arm"],
                                "arm_key": arm_key,
                                "barcode": barcode,
                                "stratum": scrutiny_row.get(
                                    "stratum", "unmapped"),
                                "planned_identity": planned_identity or "NA",
                                "input_identity": input_identity or "NA",
                                "final_identity": final_identity or "NA",
                                "input_vs_planned_changed": int(
                                    planned_mismatch),
                                "final_vs_input_changed": int(final_changed),
                                "input_missing": int(not input_identity),
                                "final_missing": int(not final_identity),
                            }

        def iter_planned_contrast_cells():
            """Yield paired-cell contrast rows without retaining row dicts."""
            for biological_condition in biological_conditions:
                for library in libraries:
                    scrutiny = scrutiny_by_library.get(library, {})
                    for contrast_name, left_arm, right_arm in contrast_plan:
                        left_key = arm_series.get(
                            (biological_condition, left_arm))
                        right_key = arm_series.get(
                            (biological_condition, right_arm))
                        left_rates = rates.get((left_key, library), {})
                        right_rates = rates.get((right_key, library), {})
                        for barcode in sorted(set(left_rates) & set(right_rates)):
                            yield {
                                "library": library,
                                "condition": biological_condition,
                                "contrast": contrast_name,
                                "barcode": barcode,
                                "stratum": scrutiny.get(
                                    barcode, {}).get("stratum", "unmapped"),
                                "left_rate": left_rates[barcode],
                                "right_rate": right_rates[barcode],
                                "right_minus_left": (
                                    right_rates[barcode] - left_rates[barcode]),
                            }

        _ambient_atomic_tsv(
            data_dir / "reconciliation_arm_contracts.tsv",
            ["library", "condition", "series_key", "arm", "arm_key",
             "assignment_basis", "roster_basis", "candidate_set",
             "plan_fingerprint", "assignment_path",
             "receiver_lines", "ambient_candidates", "scrutiny_cells",
             "context_path", "assignment_update_mode",
             "assignment_score_basis", "n_expected_sources",
             "n_observed_sources", "missing_sources", "unexpected_sources",
             "profile_raw_total"],
            contract_rows,
        )
        _ambient_atomic_tsv(
            data_dir / "reconciliation_rate_metrics.tsv",
            ["library", "condition", "series_key", "arm", "arm_key",
             "assignment_basis", "roster_basis", "population_kind",
             "population", "n_cells", "mean", "median", "q25", "q75",
             "iqr"],
            reconciliation_rate_rows,
        )
        _ambient_atomic_tsv(
            data_dir / "reconciliation_receiver_identity_metrics.tsv",
            ["library", "condition", "series_key", "arm", "arm_key",
             "assignment_basis", "roster_basis", "receiver_identity",
             "n_cells", "mean", "median", "q25", "q75", "iqr"],
            receiver_identity_rows,
        )
        _ambient_atomic_tsv(
            data_dir / "reconciliation_normalized_profiles.tsv",
            ["library", "condition", "series_key", "arm", "arm_key",
             "assignment_basis", "roster_basis", "source_label",
             "raw_profile_value", "normalized_profile_fraction",
             "profile_concentration", "raw_profile_total"],
            profile_component_rows,
        )
        _ambient_atomic_tsv(
            data_dir / "reconciliation_exact_donor_burden.tsv",
            ["library", "condition", "series_key", "arm", "arm_key",
             "assignment_basis", "roster_basis", "population",
             "source_label", "n_cells", "mean_scoring_profile_mass",
             "mean_exact_contam_burden", "fraction_of_mean_contam_burden",
             "mean_total_contam_rate"],
            donor_burden_rows,
        )
        _ambient_atomic_tsv(
            data_dir / "reconciliation_planned_contrasts.tsv",
            ["library", "condition", "contrast", "left_arm", "right_arm",
             "left_arm_key", "right_arm_key", "population", "n_common",
             "mean_delta", "median_delta", "median_absolute_delta",
             "fraction_lower_in_right_arm"],
            contrast_rows,
        )
        _ambient_atomic_tsv(
            data_dir / "reconciliation_planned_contrast_cells.tsv",
            ["library", "condition", "contrast", "barcode", "stratum",
             "left_rate", "right_rate", "right_minus_left"],
            iter_planned_contrast_cells(),
        )
        _ambient_atomic_tsv(
            data_dir / "reconciliation_assignment_switch_audit.tsv",
            ["library", "condition", "series_key", "arm", "arm_key",
             "barcode", "stratum", "planned_identity", "input_identity",
             "final_identity", "input_vs_planned_changed",
             "final_vs_input_changed", "input_missing", "final_missing"],
            iter_assignment_switch_audit(),
        )
        _ambient_atomic_tsv(
            data_dir / "reconciliation_assignment_switch_summary.tsv",
            ["library", "condition", "series_key", "arm", "arm_key",
             "assignment_basis", "n_input", "n_final", "n_union",
             "n_input_vs_planned_mismatch", "n_missing_planned_identity",
             "n_final_vs_input_switch", "frozen_assignment_invariant"],
            switch_summary_rows,
        )
        _ambient_atomic_tsv(
            data_dir / "reconciliation_diagnostics.tsv",
            ["library", "condition", "arm", "arm_key", "check", "status",
             "observed", "expected", "details"],
            diagnostic_rows,
        )

    # Paired descriptive deltas are produced only when the caller selected
    # multiple conditions.  No baseline is added implicitly.
    if len(conditions) > 1 and not four_arm_mode:
        reference = spec.get("reference_condition")
        if not reference:
            reference = (
                AMBIENT_PLOT_DEFAULT_CONDITION
                if AMBIENT_PLOT_DEFAULT_CONDITION in conditions else conditions[0]
            )
        comparators = [condition for condition in conditions if condition != reference]
        paired_rows = []
        scatter_data = {condition: [[], []] for condition in comparators}
        for condition in comparators:
            for library in libraries:
                reference_values = rates.get((reference, library), {})
                condition_values = rates.get((condition, library), {})
                common = sorted(set(reference_values) & set(condition_values))
                if not common:
                    continue
                x_values = np.asarray([reference_values[key] for key in common], dtype=float)
                y_values = np.asarray([condition_values[key] for key in common], dtype=float)
                delta = y_values - x_values
                correlation = (
                    float(np.corrcoef(x_values, y_values)[0, 1])
                    if len(common) > 1 and np.std(x_values) > 0 and np.std(y_values) > 0
                    else math.nan
                )
                paired_rows.append({
                    "reference_condition": reference,
                    "condition": condition,
                    "library": library,
                    "n_common": len(common),
                    "mean_delta": float(np.mean(delta)),
                    "median_delta": float(np.median(delta)),
                    "q25_delta": float(np.percentile(delta, 25)),
                    "q75_delta": float(np.percentile(delta, 75)),
                    "pearson_r": correlation,
                })
                # Bound figure size without changing the full-table statistics.
                stride = max(1, int(math.ceil(len(common) / 2500)))
                scatter_data[condition][0].extend(x_values[::stride].tolist())
                scatter_data[condition][1].extend(y_values[::stride].tolist())

        _ambient_atomic_tsv(
            data_dir / "paired_reference_deltas.tsv",
            ["reference_condition", "condition", "library", "n_common",
             "mean_delta", "median_delta", "q25_delta", "q75_delta", "pearson_r"],
            paired_rows,
        )
        if paired_rows:
            paired_lookup = {
                (row["condition"], int(row["library"])): float(row["median_delta"])
                for row in paired_rows
            }
            if len(comparators) == 1:
                comparator = comparators[0]
                rows = [
                    row for row in paired_rows
                    if row["condition"] == comparator
                ]
                rows.sort(key=lambda row: int(row["library"]))
                x = np.arange(len(rows))
                values = np.asarray(
                    [float(row["median_delta"]) for row in rows], dtype=float)
                limit = max(
                    0.05,
                    1.25 * float(np.max(np.abs(values))) if values.size else 0.05,
                )
                colors = [
                    "#CC6677" if value > 0 else
                    "#4477AA" if value < 0 else "#777777"
                    for value in values
                ]
                fig, ax = plt.subplots(
                    figsize=(max(8.5, 0.34 * len(rows) + 3.5), 5.3))
                ax.vlines(x, 0, values, color=colors, linewidth=2.0, alpha=0.8)
                ax.scatter(x, values, color=colors, s=38, zorder=3)
                ax.axhline(0, color="#333333", linewidth=0.9)
                annotation_offset = 0.035 * limit
                for position, value in zip(x, values):
                    ax.text(
                        position,
                        value + (annotation_offset if value >= 0 else -annotation_offset),
                        f"{value:+.3f}",
                        ha="center",
                        va="bottom" if value >= 0 else "top",
                        fontsize=7,
                    )
                ax.set_ylim(-limit, limit)
                ax.set_xticks(x)
                ax.set_xticklabels(
                    [f"lib{int(row['library'])}" for row in rows],
                    rotation=90 if len(rows) > 12 else 0,
                    fontsize=7,
                )
                ax.set_xlabel("Library")
                ax.set_ylabel("Median contamination-rate change")
                reference_item = series_by_key[reference]
                comparator_item = series_by_key[comparator]
                source_comparison = (
                    reference_item["condition"] == comparator_item["condition"]
                    and reference_item["assignment_source"]
                    != comparator_item["assignment_source"]
                )
                if source_comparison:
                    title = "Median contamination-rate change after identity reconciliation"
                else:
                    title = "Paired median contamination-rate change by library"
                fig.suptitle(title, fontsize=13, fontweight="bold", y=0.97)
                fig.text(
                    0.5, 0.915,
                    f"{condition_labels[comparator]} minus {condition_labels[reference]}",
                    ha="center", fontsize=9,
                )
                ax.grid(axis="y", alpha=0.2)
                fig.subplots_adjust(left=0.11, right=0.98, bottom=0.18, top=0.83)
                save_figure(fig, overview_dir / "paired_reference_median_deltas")
            else:
                matrix = np.full((len(comparators), len(libraries)), np.nan)
                for row_index, condition in enumerate(comparators):
                    for column_index, library in enumerate(libraries):
                        matrix[row_index, column_index] = paired_lookup.get(
                            (condition, library), np.nan)
                finite = np.abs(matrix[np.isfinite(matrix)])
                limit = (
                    max(0.05, float(np.percentile(finite, 95)))
                    if finite.size else 0.1
                )
                cmap = matplotlib.colormaps.get_cmap("RdBu_r").copy()
                cmap.set_bad("#E0E0E0")
                masked = np.ma.masked_invalid(matrix)
                fig, ax = plt.subplots(
                    figsize=(max(10, 0.34 * len(libraries) + 4),
                             max(4.5, 0.58 * len(comparators) + 2.5)))
                image = ax.imshow(
                    masked, aspect="auto", cmap=cmap,
                    vmin=-limit, vmax=limit, interpolation="nearest")
                if matrix.size <= 100:
                    for row_index in range(matrix.shape[0]):
                        for column_index in range(matrix.shape[1]):
                            value = matrix[row_index, column_index]
                            if np.isfinite(value):
                                red, green, blue, _ = image.cmap(image.norm(value))
                                luminance = 0.2126 * red + 0.7152 * green + 0.0722 * blue
                                ax.text(
                                    column_index, row_index, f"{value:+.3f}",
                                    ha="center", va="center", fontsize=7,
                                    color="black" if luminance > 0.55 else "white",
                                )
                ax.set_xticks(range(len(libraries)))
                ax.set_xticklabels(
                    [f"lib{x}" for x in libraries],
                    rotation=90 if len(libraries) > 12 else 0,
                    fontsize=6,
                )
                ax.set_yticks(range(len(comparators)))
                ax.set_yticklabels(
                    [condition_labels[x] for x in comparators], fontsize=7)
                ax.set_xlabel("Library")
                ax.set_title(
                    "Paired median contamination-rate change by library\n"
                    f"Each row minus {condition_labels[reference]}")
                fig.colorbar(
                    image, ax=ax, label="Median contamination-rate change")
                fig.subplots_adjust(left=0.17, right=0.94, bottom=0.18, top=0.86)
                save_figure(fig, overview_dir / "paired_reference_median_deltas")

            ncols = min(3, len(comparators))
            nrows = int(math.ceil(len(comparators) / ncols))
            fig, axes = plt.subplots(nrows, ncols,
                                     figsize=(5.2 * ncols, 4.8 * nrows),
                                     squeeze=False)
            for panel, condition in enumerate(comparators):
                ax = axes.ravel()[panel]
                x_values = np.asarray(scatter_data[condition][0], dtype=float)
                y_values = np.asarray(scatter_data[condition][1], dtype=float)
                ax.scatter(x_values, y_values, s=2, alpha=0.12,
                           color=condition_colors[condition], rasterized=True)
                ax.plot([0, 1], [0, 1], "k--", linewidth=0.8, alpha=0.6)
                ax.set_xlim(0, 1)
                ax.set_ylim(0, 1)
                ax.set_aspect("equal", adjustable="box")
                ax.set_xlabel(condition_labels[reference])
                ax.set_ylabel(condition_labels[condition])
                ax.set_title(
                    f"{condition_labels[condition]}\n"
                    "paired cells across selected libraries",
                    fontsize=8,
                )
            for panel in range(len(comparators), nrows * ncols):
                axes.ravel()[panel].axis("off")
            fig.suptitle(
                "Paired cell contamination rates\n"
                f"Reference: {condition_labels[reference]}"
            )
            fig.tight_layout(rect=[0, 0, 1, 0.96])
            save_figure(fig, overview_dir / "paired_reference_scatter")

    # When the same estimator condition was run from both original demux and
    # reconciled identities, write the cell-level audit that answers whether
    # corrected identities change the inferred ambient burden.  This is a
    # comparison, not a ranking of estimator conditions.
    identity_compare_rows = []
    identity_summary_rows = []
    changed_deltas = []
    unchanged_deltas = []
    for biological_condition in spec.get("conditions", []):
        by_source = {
            item["assignment_source"]: item["key"]
            for item in series
            if item["condition"] == biological_condition
        }
        if not {"demux", "reconciled"}.issubset(by_source):
            continue
        demux_key = by_source["demux"]
        reconciled_key = by_source["reconciled"]
        for library in libraries:
            demux_rates = rates.get((demux_key, library), {})
            reconciled_rates = rates.get((reconciled_key, library), {})
            common = sorted(set(demux_rates) & set(reconciled_rates))
            if not common:
                continue
            demux_prefix = prefixes[(demux_key, library)]
            reconciled_prefix = prefixes[(reconciled_key, library)]
            demux_ids = _ambient_load_assignment_file(
                f"{demux_prefix}.assignments")
            reconciled_ids = _ambient_load_assignment_file(
                f"{reconciled_prefix}.assignments")
            local_changed = []
            local_unchanged = []
            for barcode in common:
                demux_identity = demux_ids.get(barcode, "")
                reconciled_identity = reconciled_ids.get(barcode, "")
                identity_changed = bool(
                    demux_identity and reconciled_identity and
                    demux_identity != reconciled_identity)
                delta = reconciled_rates[barcode] - demux_rates[barcode]
                target = local_changed if identity_changed else local_unchanged
                target.append(delta)
                identity_compare_rows.append({
                    "library": library,
                    "condition": biological_condition,
                    "barcode": barcode,
                    "demux_identity": demux_identity or "NA",
                    "reconciled_identity": reconciled_identity or "NA",
                    "identity_changed": int(identity_changed),
                    "demux_contam_rate": demux_rates[barcode],
                    "reconciled_contam_rate": reconciled_rates[barcode],
                    "reconciled_minus_demux": delta,
                    "absolute_rate_change": abs(delta),
                })
            changed_deltas.extend(local_changed)
            unchanged_deltas.extend(local_unchanged)
            identity_summary_rows.append({
                "library": library,
                "condition": biological_condition,
                "n_common": len(common),
                "n_identity_changed": len(local_changed),
                "fraction_identity_changed": len(local_changed) / len(common),
                "median_delta_changed": (
                    float(np.median(local_changed)) if local_changed else math.nan),
                "median_abs_delta_changed": (
                    float(np.median(np.abs(local_changed)))
                    if local_changed else math.nan),
                "median_delta_unchanged": (
                    float(np.median(local_unchanged))
                    if local_unchanged else math.nan),
                "median_abs_delta_unchanged": (
                    float(np.median(np.abs(local_unchanged)))
                    if local_unchanged else math.nan),
            })
    if identity_compare_rows:
        _ambient_atomic_tsv(
            data_dir / "identity_reconciliation_ambient_cells.tsv",
            ["library", "condition", "barcode", "demux_identity",
             "reconciled_identity", "identity_changed", "demux_contam_rate",
             "reconciled_contam_rate", "reconciled_minus_demux",
             "absolute_rate_change"],
            identity_compare_rows,
        )
        _ambient_atomic_tsv(
            data_dir / "identity_reconciliation_ambient_summary.tsv",
            ["library", "condition", "n_common", "n_identity_changed",
             "fraction_identity_changed", "median_delta_changed",
             "median_abs_delta_changed", "median_delta_unchanged",
             "median_abs_delta_unchanged"],
            identity_summary_rows,
        )
        if changed_deltas or unchanged_deltas:
            fig, ax = plt.subplots(figsize=(8, 5.5))
            data = []
            labels = []
            if unchanged_deltas:
                data.append(np.asarray(unchanged_deltas, dtype=float))
                labels.append(f"Identity unchanged (n={len(unchanged_deltas):,})")
            if changed_deltas:
                data.append(np.asarray(changed_deltas, dtype=float))
                labels.append(f"Identity corrected (n={len(changed_deltas):,})")
            ax.violinplot(data, showmedians=True, showextrema=False)
            ax.set_xticks(range(1, len(labels) + 1))
            ax.set_xticklabels(labels)
            ax.axhline(0, color="black", linewidth=0.8)
            ax.set_ylabel("Reconciled - demux contamination rate")
            ax.set_title("Ambient-rate sensitivity to identity reconciliation")
            ax.grid(axis="y", alpha=0.2)
            fig.tight_layout()
            save_figure(
                fig, overview_dir / "identity_reconciliation_rate_changes")

    def draw_kde(
            ax, values, color, linewidth=1.4, label=None, linestyle="-"):
        # dict_values is iterable but NumPy does not coerce it directly.
        values = np.asarray(list(values), dtype=float)
        values = values[np.isfinite(values)]
        if len(values) < 20:
            return False
        grid = np.linspace(0, 1, 400)
        try:
            density = gaussian_kde(values, bw_method="scott")(grid)
            ax.plot(
                grid, density, color=color, linewidth=linewidth,
                label=label, linestyle=linestyle,
            )
        except Exception:
            hist, edges = np.histogram(values, bins=40, range=(0, 1), density=True)
            centers = 0.5 * (edges[:-1] + edges[1:])
            ax.plot(
                centers, hist, color=color, linewidth=linewidth,
                label=label, linestyle=linestyle,
            )
        return True

    # One readable small-multiple overlay shows every explicitly selected
    # series in every library.  Five columns keeps the 40-library report at a
    # useful panel size instead of shrinking labels into an eight-column grid.
    if len(libraries) <= 2:
        ncols = len(libraries)
    else:
        ncols = min(5, len(libraries))
    nrows = int(math.ceil(len(libraries) / ncols))
    fig, axes = plt.subplots(
        nrows, ncols,
        figsize=(max(8.5, 3.15 * ncols), max(5.5, 2.45 * nrows + 2.0)),
        squeeze=False,
    )
    for panel, library in enumerate(libraries):
        ax = axes.ravel()[panel]
        plotted = 0
        for condition in conditions:
            key = (condition, library)
            values = rates.get(key, {})
            if four_arm_mode and values:
                biological_condition = str(
                    series_by_key[condition]["condition"])
                common = four_arm_common.get(
                    (biological_condition, library), set())
                values = {
                    barcode: values[barcode] for barcode in common
                    if barcode in values
                }
            if values and draw_kde(
                ax, values.values(), condition_colors[condition],
                label=condition_labels[condition],
                linestyle=condition_linestyles[condition],
            ):
                plotted += 1
        ax.set_xlim(0, 1)
        ax.set_yticks([])
        ax.set_xticks([0, 0.25, 0.5, 0.75, 1.0])
        ax.tick_params(axis="x", labelsize=6)
        ax.grid(axis="x", alpha=0.15)
        ax.set_title(
            f"lib{library} | {plotted}/"
            f"{sum((condition, library) in rates for condition in conditions)} "
            "applicable series", fontsize=8)
        if panel % ncols == 0:
            ax.set_ylabel("Density", fontsize=7)
        if plotted == 0:
            ax.text(0.5, 0.5, "No data", transform=ax.transAxes, ha="center")
    for panel in range(len(libraries), nrows * ncols):
        axes.ravel()[panel].axis("off")
    legend = [
        Line2D(
            [0], [0], color=condition_colors[c], linewidth=2,
            linestyle=condition_linestyles[c], label=condition_labels[c],
        )
        for c in conditions
        if any((c, library) in rates for library in libraries)
    ]
    legend_ncol = min(4, len(conditions))
    legend_rows = int(math.ceil(len(conditions) / legend_ncol))
    if len(biological_conditions) == 1:
        biological_condition = biological_conditions[0]
        figure_title = (
            f"{biological_labels[biological_condition]} ambient-RNA "
            "contamination-rate KDEs"
        )
        scope_line = biological_condition
    else:
        figure_title = "Ambient-RNA contamination-rate KDEs"
        scope_line = (
            f"{len(biological_conditions)} selected estimator conditions"
        )
    scope_line += (
        f"  •  {len(libraries)} selected "
        f"{'library' if len(libraries) == 1 else 'libraries'}"
    )
    figure_width, figure_height = fig.get_size_inches()
    fig.suptitle(
        figure_title, fontsize=14, fontweight="bold",
        y=1.0 - 0.18 / figure_height,
    )
    fig.text(
        0.5, 1.0 - 0.48 / figure_height, scope_line,
        ha="center", va="top", fontsize=8.5,
    )
    if len(legend) > 1:
        fig.legend(
            handles=legend, loc="upper center",
            bbox_to_anchor=(0.5, 1.0 - 0.78 / figure_height),
            ncol=legend_ncol, fontsize=7.5, frameon=False,
        )
        top_band_inches = 1.35 + 0.18 * (legend_rows - 1)
    else:
        top_band_inches = 1.15
    axes_top = 1.0 - top_band_inches / figure_height
    axes_bottom = 0.62 / figure_height
    fig.supxlabel(
        "Estimated contamination rate", y=0.16 / figure_height, fontsize=9)
    fig.subplots_adjust(
        left=0.55 / figure_width,
        right=1.0 - 0.18 / figure_width,
        bottom=axes_bottom, top=axes_top,
        wspace=0.22, hspace=0.36,
    )
    save_figure(fig, kde_dir / "kde_all_libraries")

    if four_arm_mode:
        # Fixed-cell stratification: every curve in a panel uses exactly the
        # same barcodes.  This prevents a changed cell roster from masquerading
        # as a contamination-distribution shift.
        strata = (
            ("changed_target", "Changed / reassigned target cells"),
            ("scrutinized_other", "Scrutinized but unchanged cells"),
            ("background", "Background cells"),
        )
        library_columns = min(5, max(1, len(libraries)))
        library_blocks = int(math.ceil(len(libraries) / library_columns))
        nrows = len(biological_conditions) * len(strata) * library_blocks
        fig, axes = plt.subplots(
            nrows, library_columns,
            figsize=(max(9.0, 3.15 * library_columns),
                     max(7.0, 2.25 * nrows + 1.7)),
            squeeze=False,
        )
        panel = 0
        for biological_condition in biological_conditions:
            for stratum, stratum_label in strata:
                for block in range(library_blocks):
                    for column in range(library_columns):
                        library_index = block * library_columns + column
                        ax = axes.ravel()[panel]
                        panel += 1
                        if library_index >= len(libraries):
                            ax.axis("off")
                            continue
                        library = libraries[library_index]
                        common = four_arm_common.get(
                            (biological_condition, library), set())
                        scrutiny = scrutiny_by_library.get(library, {})
                        barcodes = sorted(
                            barcode for barcode in common
                            if scrutiny.get(barcode, {}).get("stratum") == stratum
                        )
                        plotted = 0
                        for arm_key in IDENTITY_AMBIENT_ARM_ORDER:
                            key = arm_series.get(
                                (biological_condition, arm_key))
                            values = rates.get((key, library), {})
                            arm_values = [
                                values[barcode] for barcode in barcodes
                                if barcode in values
                            ]
                            if draw_kde(
                                    ax, arm_values,
                                    condition_colors.get(key, "#666666"),
                                    label=condition_labels.get(key, arm_key),
                                    linestyle=condition_linestyles.get(key, "-")):
                                plotted += 1
                        ax.set_xlim(0, 1)
                        ax.set_yticks([])
                        ax.set_xticks([0, 0.25, 0.5, 0.75, 1])
                        ax.tick_params(axis="x", labelsize=6)
                        ax.grid(axis="x", alpha=0.15)
                        ax.set_title(
                            f"lib{library} | {stratum_label}\n"
                            f"fixed n={len(barcodes):,}", fontsize=7.5)
                        if plotted == 0:
                            message = (
                                "KDE omitted\n(n < 20)" if barcodes
                                else "No cells in stratum")
                            ax.text(
                                0.5, 0.5, message, transform=ax.transAxes,
                                ha="center", va="center", fontsize=7,
                            )
                        if column == 0:
                            ax.set_ylabel("Density", fontsize=7)
        for unused in range(panel, nrows * library_columns):
            axes.ravel()[unused].axis("off")
        handles = [
            Line2D(
                [0], [0],
                color=condition_colors[arm_series[(biological_conditions[0], arm_key)]],
                linestyle=condition_linestyles[
                    arm_series[(biological_conditions[0], arm_key)]],
                linewidth=2,
                label=IDENTITY_AMBIENT_ARMS[arm_key]["short_label"],
            )
            for arm_key in IDENTITY_AMBIENT_ARM_ORDER
            if (biological_conditions[0], arm_key) in arm_series
            and any(
                (arm_series[(biological_conditions[0], arm_key)], library)
                in rates for library in libraries)
        ]
        fig.suptitle(
            "Reconciliation four-arm contamination KDEs by fixed cell stratum",
            fontsize=14, fontweight="bold", y=0.995,
        )
        fig.legend(
            handles=handles, loc="upper center", bbox_to_anchor=(0.5, 0.958),
            ncol=min(4, len(handles)), fontsize=8, frameon=False,
        )
        fig.supxlabel("Estimated contamination rate", fontsize=9)
        fig.subplots_adjust(
            left=0.055, right=0.99, bottom=0.025, top=0.895,
            wspace=0.22, hspace=0.58,
        )
        save_figure(fig, kde_dir / "reconciliation_arms_stratified")

        # Planned contrast heatmap.  Rows are declared scientific contrasts
        # and fixed populations; columns are libraries.  It is intentionally
        # descriptive and does not rank conditions.
        contrast_order = [
            "roster_effect_demux", "assignment_effect_augmented",
            "replacement_sensitivity", "combined_production_change",
        ]
        population_order = [
            "full_library", "changed_target", "scrutinized_other", "background",
        ]
        contrast_lookup = {
            (str(row["condition"]), str(row["contrast"]),
             str(row["population"]), int(row["library"])):
                finite_number(row["median_delta"])
            for row in contrast_rows
        }
        row_keys = [
            (biological_condition, contrast, population)
            for biological_condition in biological_conditions
            for contrast in contrast_order
            for population in population_order
            if any(
                (biological_condition, contrast, population, library)
                in contrast_lookup for library in libraries)
        ]
        if row_keys:
            matrix = np.full((len(row_keys), len(libraries)), np.nan)
            for row_index, (bio, contrast, population) in enumerate(row_keys):
                for column_index, library in enumerate(libraries):
                    matrix[row_index, column_index] = contrast_lookup.get(
                        (bio, contrast, population, library), np.nan)
            finite = np.abs(matrix[np.isfinite(matrix)])
            limit = max(
                0.02, float(np.percentile(finite, 95))
                if finite.size else 0.1)
            cmap = matplotlib.colormaps.get_cmap("RdBu_r").copy()
            cmap.set_bad("#E0E0E0")
            fig, ax = plt.subplots(
                figsize=(max(9, 0.36 * len(libraries) + 5),
                         max(7, 0.36 * len(row_keys) + 2.5)))
            image = ax.imshow(
                np.ma.masked_invalid(matrix), aspect="auto", cmap=cmap,
                vmin=-limit, vmax=limit, interpolation="nearest")
            if matrix.size <= 140:
                for row_index in range(matrix.shape[0]):
                    for column_index in range(matrix.shape[1]):
                        value = matrix[row_index, column_index]
                        if np.isfinite(value):
                            ax.text(
                                column_index, row_index, f"{value:+.3f}",
                                ha="center", va="center", fontsize=6,
                            )
            ax.set_xticks(range(len(libraries)))
            ax.set_xticklabels(
                [f"lib{x}" for x in libraries],
                rotation=90 if len(libraries) > 12 else 0, fontsize=7)
            ax.set_yticks(range(len(row_keys)))
            ax.set_yticklabels([
                f"{biological_labels[bio]} | {contrast} | {population}"
                for bio, contrast, population in row_keys
            ], fontsize=6.5)
            ax.set_title(
                "Planned reconciliation contrasts\n"
                "Median paired contamination-rate delta (right arm − left arm)")
            fig.colorbar(image, ax=ax, label="Median paired Δc")
            fig.tight_layout()
            save_figure(
                fig, overview_dir / "reconciliation_planned_contrasts")

        # Freeze-assignment audit remains visible even when all expected bars
        # are zero; a zero-only panel is the desired all-clear result.
        if switch_summary_rows:
            rows = sorted(
                switch_summary_rows,
                key=lambda row: (
                    int(row["library"]), str(row["condition"]),
                    str(row["arm"])),
            )
            labels = [
                f"lib{row['library']} {row['arm']}" for row in rows]
            x = np.arange(len(rows))
            planned = np.asarray([
                int(row["n_input_vs_planned_mismatch"]) for row in rows])
            frozen = np.asarray([
                int(row["n_final_vs_input_switch"]) for row in rows])
            if int(planned.sum()) == 0 and int(frozen.sum()) == 0:
                fig, ax = plt.subplots(figsize=(8.5, 3.8))
                ax.axis("off")
                ax.text(
                    0.5, 0.62, "PASS", ha="center", va="center",
                    fontsize=30, fontweight="bold", color="#228833",
                    transform=ax.transAxes)
                ax.text(
                    0.5, 0.38,
                    f"0 assignment discrepancies across {len(rows)} arm/library checks",
                    ha="center", va="center", fontsize=12,
                    transform=ax.transAxes)
                ax.set_title(
                    "Reconciliation four-arm assignment-switch audit",
                    pad=14)
            else:
                fig, ax = plt.subplots(
                    figsize=(max(9, 0.42 * len(rows) + 4), 5.5))
                width = 0.38
                ax.bar(
                    x - width / 2, planned, width,
                    label="Input vs planned-basis mismatch", color="#CCBB44")
                ax.bar(
                    x + width / 2, frozen, width,
                    label="Estimator final vs frozen input switch", color="#CC6677")
                ax.axhline(0, color="#333333", linewidth=0.8)
                ax.set_xticks(x)
                ax.set_xticklabels(
                    labels, rotation=90 if len(labels) > 8 else 20,
                    ha="right", fontsize=7)
                ax.set_ylabel("Number of assignment discrepancies")
                ax.set_title("Reconciliation four-arm assignment-switch audit")
                ax.legend(frameon=False, fontsize=8)
                ax.grid(axis="y", alpha=0.2)
            fig.tight_layout()
            save_figure(
                fig, overview_dir / "reconciliation_assignment_switch_audit")

    # Transposed view: one condition per file, all libraries as panels.
    for condition in conditions:
        available = [library for library in libraries if (condition, library) in rates]
        if not available:
            continue
        ncols = min(5, max(1, len(available)))
        nrows = int(math.ceil(len(available) / ncols))
        fig, axes = plt.subplots(
            nrows, ncols,
            figsize=(max(8.0, 3.1 * ncols), max(5.0, 2.45 * nrows + 1.25)),
            squeeze=False,
        )
        for panel, library in enumerate(available):
            ax = axes.ravel()[panel]
            values = np.asarray(list(rates[(condition, library)].values()), dtype=float)
            draw_kde(
                ax, values, condition_colors[condition], linewidth=1.5,
                linestyle=condition_linestyles[condition],
            )
            ax.set_xlim(0, 1)
            ax.set_yticks([])
            ax.set_xticks([0, 0.25, 0.5, 0.75, 1.0])
            ax.tick_params(axis="x", labelsize=6)
            ax.grid(axis="x", alpha=0.15)
            ax.set_title(f"lib{library} | n={len(values):,} | med={np.median(values):.3f}", fontsize=7.5)
            if panel % ncols == 0:
                ax.set_ylabel("Density", fontsize=7)
        for panel in range(len(available), nrows * ncols):
            axes.ravel()[panel].axis("off")
        item = series_by_key[condition]
        series_title = condition_labels[condition]
        if len(biological_conditions) > 1:
            series_title = (
                f"{biological_labels[item['condition']]} | {series_title}"
            )
        figure_width, figure_height = fig.get_size_inches()
        fig.suptitle(
            f"{series_title}\nContamination-rate KDE by library",
            fontsize=13, fontweight="bold",
            y=1.0 - 0.14 / figure_height,
        )
        fig.supxlabel(
            "Estimated contamination rate", y=0.15 / figure_height,
            fontsize=9,
        )
        fig.subplots_adjust(
            left=0.55 / figure_width,
            right=1.0 - 0.18 / figure_width,
            bottom=0.60 / figure_height,
            top=1.0 - 1.15 / figure_height,
            wspace=0.22, hspace=0.36,
        )
        save_figure(fig, per_condition_dir / f"{_ambient_condition_slug(condition)}_kde")

    # Multi-library QC overview.  A heatmap is intentionally omitted for very
    # small selections, where one or two colored rectangles look quantitative
    # but communicate less than the metrics TSV and KDE/paired plots.
    metric_by_key = {(row["condition"], int(row["library"])): row for row in metrics_rows}
    ready_libraries = [
        library for library in libraries
        if any((condition, library) in metric_by_key for condition in conditions)
    ]
    qc_stem = overview_dir / "ambient_rate_qc_heatmaps"
    if len(ready_libraries) < 4:
        skip_figure(
            qc_stem,
            "requires at least four libraries; use ambient_rate_metrics.tsv "
            "and the KDE/paired plots for small selections",
        )
    else:
        def observed_upper(field, floor):
            values = [
                float(row[field]) for row in metrics_rows
                if math.isfinite(float(row.get(field, math.nan)))
            ]
            if not values:
                return floor
            return min(1.0, max(floor, 1.10 * float(np.percentile(values, 98))))

        panels = [
            ("Median contamination rate", "median", "viridis",
             0.0, observed_upper("median", 0.30), ".2f"),
            ("Interquartile range", "iqr", "YlGnBu",
             0.0, observed_upper("iqr", 0.20), ".2f"),
            ("Cells with finite standard error", "frac_se_finite", "RdYlGn",
             0.0, 1.0, ".2f"),
            ("Median positive standard error", "median_se_nonzero", "magma",
             0.0, observed_upper("median_se_nonzero", 0.15), ".2f"),
        ]
        fig, axes = plt.subplots(
            2, 2,
            figsize=(max(13, 0.34 * len(ready_libraries) + 6),
                     max(8, 0.48 * len(conditions) + 5)),
            squeeze=False,
        )
        for ax, (title, field, cmap_name, vmin, vmax, fmt) in zip(
                axes.ravel(), panels):
            matrix = np.full(
                (len(conditions), len(ready_libraries)), np.nan)
            for row_index, condition in enumerate(conditions):
                for column_index, library in enumerate(ready_libraries):
                    row = metric_by_key.get((condition, library))
                    if row:
                        value = row.get(field, math.nan)
                        matrix[row_index, column_index] = (
                            float(value) if value is not None else np.nan)
            cmap = matplotlib.colormaps.get_cmap(cmap_name).copy()
            cmap.set_bad("#E0E0E0")
            image = ax.imshow(
                np.ma.masked_invalid(matrix), aspect="auto", cmap=cmap,
                vmin=vmin, vmax=vmax, interpolation="nearest",
            )
            if matrix.size <= 120 and len(ready_libraries) <= 20:
                annotation_size = 7 if matrix.size <= 60 else 5.5
                for row_index in range(matrix.shape[0]):
                    for column_index in range(matrix.shape[1]):
                        value = matrix[row_index, column_index]
                        if np.isfinite(value):
                            red, green, blue, _ = image.cmap(image.norm(value))
                            luminance = (
                                0.2126 * red + 0.7152 * green + 0.0722 * blue)
                            ax.text(
                                column_index, row_index, format(value, fmt),
                                ha="center", va="center",
                                fontsize=annotation_size,
                                color="black" if luminance > 0.55 else "white",
                            )
            ax.set_xticks(range(len(ready_libraries)))
            ax.set_xticklabels(
                [f"lib{x}" for x in ready_libraries],
                rotation=90 if len(ready_libraries) > 12 else 0,
                fontsize=6,
            )
            ax.set_yticks(range(len(conditions)))
            ax.set_yticklabels(
                [condition_labels[x] for x in conditions], fontsize=7)
            ax.set_xlabel("Library", fontsize=8)
            ax.set_title(title, fontsize=10, fontweight="bold")
            fig.colorbar(image, ax=ax, fraction=0.03, pad=0.02)
        fig.suptitle(
            "Ambient-RNA contamination QC across selected libraries",
            fontsize=14, fontweight="bold", y=0.985,
        )
        fig.text(
            0.5, 0.025,
            "Gray = missing. SE coverage is the fraction of rate-bearing "
            "cells with a finite standard error.",
            ha="center", fontsize=8,
        )
        fig.subplots_adjust(
            left=0.15, right=0.95, bottom=0.12, top=0.91,
            wspace=0.32, hspace=0.42,
        )
        save_figure(fig, qc_stem)

    # The former standalone SE page duplicated information and concealed
    # missing-SE coverage.  Its useful metric now lives in the QC overview.
    skip_figure(
        overview_dir / "contamination_se_diagnostics",
        "retired: finite-SE coverage and median positive SE are in "
        "ambient_rate_qc_heatmaps",
    )

    # Per-library medians use a small-run paired-point fallback.  A one-value
    # boxplot is not a distribution; connected points make the actual
    # demux/reconciled comparison explicit instead.
    box_series = []
    for condition in conditions:
        values = [
            (library, float(metric_by_key[(condition, library)]["median"]))
            for library in libraries if (condition, library) in metric_by_key
        ]
        if values:
            box_series.append((condition, values))
    median_stem = overview_dir / "library_median_distributions"
    all_medians = [value for _, values in box_series for _, value in values]
    if not all_medians or (len(ready_libraries) < 3 and len(box_series) < 2):
        skip_figure(
            median_stem,
            "a distribution needs at least three libraries, or two plotted "
            "series for the small-run paired-point view",
        )
    elif len(ready_libraries) < 3:
        fig, ax = plt.subplots(
            figsize=(max(8.0, 1.65 * len(box_series) + 3.5), 5.6))
        x = np.arange(len(box_series))
        by_library = {}
        for series_index, (condition, values) in enumerate(box_series):
            for library, value in values:
                by_library.setdefault(library, []).append(
                    (series_index, value, condition))
        for library, points in sorted(by_library.items()):
            points.sort()
            ax.plot(
                [point[0] for point in points],
                [point[1] for point in points],
                color="#999999", linewidth=1.1, alpha=0.75, zorder=1,
            )
            for series_index, value, condition in points:
                ax.scatter(
                    series_index, value, s=55,
                    color=condition_colors[condition], zorder=3,
                )
                if sum(len(values) for _, values in box_series) <= 12:
                    ax.annotate(
                        f"lib{library}: {value:.3f}",
                        (series_index, value), xytext=(0, 7),
                        textcoords="offset points", ha="center", fontsize=7,
                    )
        ax.set_xticks(x)
        ax.set_xticklabels(
            [condition_labels[condition] for condition, _ in box_series],
            rotation=18 if len(box_series) > 3 else 0,
            ha="right" if len(box_series) > 3 else "center",
        )
        upper = min(1.0, max(0.30, 1.15 * max(all_medians)))
        ax.set_ylim(0, upper)
        ax.set_ylabel("Per-library median contamination rate")
        ax.set_title(
            "Per-library median contamination rate by plotted series",
            pad=18,
        )
        ax.text(
            0.5, 1.01, "Connected points represent the same library",
            transform=ax.transAxes, ha="center", va="bottom", fontsize=8,
        )
        ax.grid(axis="y", alpha=0.2)
        fig.subplots_adjust(left=0.13, right=0.98, bottom=0.18, top=0.84)
        save_figure(fig, median_stem)
    else:
        fig, ax = plt.subplots(
            figsize=(max(9.0, 0.25 * len(ready_libraries) + 6),
                     max(5.0, 0.62 * len(box_series) + 2.5)))
        box_values = [
            [value for _, value in values] for _, values in box_series
        ]
        box_labels = [
            f"{condition_labels[condition]} (n={len(values)})"
            for condition, values in box_series
        ]
        try:
            matplotlib_version = tuple(
                int(part) for part in matplotlib.__version__.split(".")[:2])
        except (TypeError, ValueError):
            matplotlib_version = (0, 0)
        label_argument = (
            {"tick_labels": box_labels}
            if matplotlib_version >= (3, 9)
            else {"labels": box_labels}
        )
        boxes = ax.boxplot(
            box_values, vert=False, patch_artist=True, showfliers=False,
            **label_argument,
        )
        for patch, (condition, _) in zip(boxes["boxes"], box_series):
            patch.set_facecolor(condition_colors[condition])
            patch.set_alpha(0.55)
        for index, (_, values) in enumerate(box_series, start=1):
            offsets = (
                np.linspace(-0.12, 0.12, len(values))
                if len(values) > 1 else np.array([0.0])
            )
            ax.scatter(
                [value for _, value in values], index + offsets,
                s=14, color="#222222", alpha=0.58, zorder=3,
            )
        upper = min(1.0, max(0.30, 1.15 * max(all_medians)))
        ax.set_xlim(0, upper)
        ax.invert_yaxis()
        ax.set_xlabel("Per-library median contamination rate")
        ax.set_title("Distribution of per-library median contamination rates")
        ax.grid(axis="x", alpha=0.2)
        fig.subplots_adjust(left=0.24, right=0.98, bottom=0.13, top=0.90)
        save_figure(fig, median_stem)

    # Profile figures and descriptive L1 distances.  The first selected
    # condition with assignments supplies Model A/B; there is no hidden
    # IND_BASE dependency.
    profile_distance_rows = []
    species_colors = {"H": "#4477AA", "C": "#EE6677", "B": "#228833", "O": "#CCBB44", "Hy": "#AA3377"}

    def stacked_profile_figure(
            rows, title, xlabel, colors=None, x_upper=1.0):
        components = sorted(set().union(*(profile.keys() for _, profile in rows)))
        fig, ax = plt.subplots(figsize=(11, max(4, 0.48 * len(rows) + 1.8)))
        left = np.zeros(len(rows))
        cmap = matplotlib.colormaps.get_cmap("tab20")
        for index, component in enumerate(components):
            values = np.array([profile.get(component, 0.0) for _, profile in rows])
            color = (colors or {}).get(component, cmap(index % 20))
            ax.barh(np.arange(len(rows)), values, left=left, height=0.72, color=color, label=component)
            left += values
        ax.set_yticks(np.arange(len(rows)))
        ax.set_yticklabels([label for label, _ in rows], fontsize=7)
        ax.set_xlim(0, x_upper)
        ax.invert_yaxis()
        ax.set_xlabel(xlabel)
        ax.set_title(title, fontweight="bold")
        ax.legend(loc="center left", bbox_to_anchor=(1.01, 0.5), fontsize=6, ncol=1 if len(components) < 16 else 2)
        ax.grid(axis="x", alpha=0.15)
        fig.tight_layout()
        return fig

    for library in libraries:
        demux_dir = (
            Path(spec["mapping_root"])
            / f"{spec.get('library_prefix', 'Tet_2025_Multiome-RNA_')}{library}"
            / spec.get("demux_subdir", "demux_nomito")
        )
        assignment_map = {}
        for condition in conditions:
            assignment_map = assignments.get((condition, library), {})
            if assignment_map:
                break
        if not assignment_map:
            assignment_map, _ = _ambient_load_assignments(
                str(demux_dir / f"lib{library}_demuxed")
            )
        model_a = _ambient_assignment_profile(assignment_map, equal_identities=True)
        model_b = _ambient_assignment_profile(assignment_map, equal_identities=False)
        empty_individual = _ambient_load_two_column_profile(
            demux_dir / f"lib{library}_raw.contam_prof_empty"
        )
        empty_species = _ambient_load_two_column_profile(
            demux_dir / f"lib{library}_raw.species_prof_empty"
        )

        individual_rows = []
        if model_a and not four_arm_mode:
            individual_rows.append(("Model A (equal identities)", model_a))
        if model_b and not four_arm_mode:
            individual_rows.append(("Model B (cell proportional)", model_b))
        if empty_individual and not four_arm_mode:
            individual_rows.append(("Empty drops", empty_individual))
        n_individual_conditions = 0
        for condition in conditions:
            profile = individual_profiles.get((condition, library), {})
            if not profile:
                profile = _ambient_expand_species_profile(
                    species_profiles.get((condition, library), {}), model_b
                )
            if profile:
                individual_rows.append((condition_labels[condition], profile))
                n_individual_conditions += 1
                if not four_arm_mode:
                    profile_distance_rows.append({
                        "condition": condition,
                        "library": library,
                        "resolution": "individual",
                        "l1_model_a": _ambient_profile_l1(profile, model_a),
                        "l1_model_b": _ambient_profile_l1(profile, model_b),
                        "l1_empty_drops": _ambient_profile_l1(
                            profile, empty_individual),
                    })
        if n_individual_conditions:
            fig = stacked_profile_figure(
                individual_rows,
                (f"lib{library}: reconciliation-arm normalized ambient "
                 "profile fractions" if four_arm_mode else
                 f"lib{library}: individual ambient profile"),
                ("Normalized profile fraction (composition, not absolute burden)"
                 if four_arm_mode else "Ambient fraction"),
            )
            save_figure(fig, profile_dir / f"lib{library}_ambient_profiles")

        if four_arm_mode:
            population_labels = {
                "full_library": "all fixed-common cells",
                "changed_target": "changed target cells",
                "scrutinized_other": "scrutinized unchanged cells",
                "background": "background cells",
            }
            burden_lookup = {}
            arm_label_lookup = {}
            for row in donor_burden_rows:
                if int(row["library"]) != library:
                    continue
                row_key = (str(row["series_key"]), str(row["population"]))
                burden_lookup.setdefault(row_key, {})[
                    str(row["source_label"])] = float(
                        row["mean_exact_contam_burden"])
                arm_label_lookup[str(row["series_key"])] = (
                    f"{row['arm']} {row['assignment_basis']}/"
                    f"{row['roster_basis']}")
            burden_rows_for_plot = []
            for condition in conditions:
                if series_by_key[condition].get("arm") not in {
                        "A", "B", "C", "D"}:
                    continue
                for population in (
                        "full_library", "changed_target",
                        "scrutinized_other", "background"):
                    profile = burden_lookup.get((condition, population), {})
                    if profile:
                        burden_rows_for_plot.append((
                            f"{arm_label_lookup.get(condition, condition)} | "
                            f"{population_labels[population]}",
                            profile,
                        ))
            if burden_rows_for_plot:
                burden_upper = min(
                    1.0,
                    max(
                        0.05,
                        1.10 * max(
                            sum(profile.values())
                            for _, profile in burden_rows_for_plot
                        ),
                    ),
                )
                fig = stacked_profile_figure(
                    burden_rows_for_plot,
                    f"lib{library}: exact donor-specific ambient burden",
                    "Mean c × cell-specific scoring-profile mass",
                    x_upper=burden_upper,
                )
                save_figure(
                    fig,
                    profile_dir /
                    f"lib{library}_reconciliation_exact_donor_burden",
                )

        model_a_species = _ambient_profile_to_species(model_a)
        model_b_species = _ambient_profile_to_species(model_b)
        species_rows = []
        if model_a_species and not four_arm_mode:
            species_rows.append(("Model A (equal identities)", model_a_species))
        if model_b_species and not four_arm_mode:
            species_rows.append(("Model B (cell proportional)", model_b_species))
        if empty_species and not four_arm_mode:
            species_rows.append(("Empty drops", empty_species))
        n_species_conditions = 0
        for condition in conditions:
            profile = species_profiles.get((condition, library), {})
            if not profile:
                profile = _ambient_profile_to_species(
                    individual_profiles.get((condition, library), {})
                )
            else:
                profile = _ambient_profile_to_species(profile)
            if profile:
                species_rows.append((condition_labels[condition], profile))
                n_species_conditions += 1
                if not four_arm_mode:
                    profile_distance_rows.append({
                        "condition": condition,
                        "library": library,
                        "resolution": "species",
                        "l1_model_a": _ambient_profile_l1(
                            profile, model_a_species),
                        "l1_model_b": _ambient_profile_l1(
                            profile, model_b_species),
                        "l1_empty_drops": _ambient_profile_l1(
                            profile, empty_species),
                    })
        if n_species_conditions:
            fig = stacked_profile_figure(
                species_rows,
                f"lib{library}: species-level ambient profile",
                "Ambient fraction",
                colors=species_colors,
            )
            save_figure(fig, profile_dir / f"lib{library}_species_profiles")

    _ambient_atomic_tsv(
        data_dir / "ambient_profile_distances.tsv",
        ["condition", "library", "resolution", "l1_model_a", "l1_model_b", "l1_empty_drops"],
        profile_distance_rows,
    )

    gate_stem = overview_dir / "geometry_gate_real_library_overview"
    gate_table = data_dir / "geometry_gate_by_library.tsv"
    if gate_rows:
        _ambient_atomic_tsv(
            gate_table,
            ["condition", "library", "n_rows", "n_triggered",
             "fraction_triggered", "mean_selected_minus_base_c",
             "mean_selected_minus_base_c_triggered", "audit_path"],
            gate_rows,
        )
        if sum(int(row["n_triggered"]) for row in gate_rows) == 0:
            skip_figure(
                gate_stem,
                "no geometry-gate activations in the selected libraries; "
                "see geometry_gate_by_library.tsv for the all-clear record",
            )
        else:
            fig, axes = plt.subplots(
                2, 1,
                figsize=(max(10, 0.38 * len(gate_rows) + 4), 8.2),
                sharex=True,
            )
            labels = [
                f"lib{row['library']}"
                if len(set(x["condition"] for x in gate_rows)) == 1
                else (
                    f"{condition_labels.get(row['condition'], row['condition'])} "
                    f"| lib{row['library']}"
                )
                for row in gate_rows
            ]
            colors = [
                condition_colors.get(row["condition"], "#4477AA")
                for row in gate_rows
            ]
            x = np.arange(len(gate_rows))
            fractions = np.asarray(
                [float(row["fraction_triggered"]) for row in gate_rows])
            fraction_bars = axes[0].bar(x, fractions, color=colors, alpha=0.82)
            axes[0].set_ylim(0, 1)
            axes[0].set_ylabel("Fraction gate-triggered")
            axes[0].set_title("Geometry-gate activation")
            for bar, row, value in zip(fraction_bars, gate_rows, fractions):
                axes[0].text(
                    bar.get_x() + bar.get_width() / 2,
                    min(0.98, value + 0.025),
                    f"{int(row['n_triggered'])}/{int(row['n_rows'])}",
                    ha="center", va="bottom", fontsize=7,
                )

            movements = np.asarray([
                float(row.get(
                    "mean_selected_minus_base_c_triggered", math.nan))
                for row in gate_rows
            ])
            finite_movements = np.abs(movements[np.isfinite(movements)])
            movement_limit = (
                max(0.05, 1.20 * float(finite_movements.max()))
                if finite_movements.size else 0.05
            )
            movement_values = np.nan_to_num(movements, nan=0.0)
            movement_bars = axes[1].bar(
                x, movement_values, color=colors, alpha=0.82)
            axes[1].axhline(0, color="#333333", linewidth=0.9)
            axes[1].set_ylim(-movement_limit, movement_limit)
            axes[1].set_ylabel("Mean selected c − base c\n(triggered cells)")
            axes[1].set_xticks(x)
            axes[1].set_xticklabels(
                labels, rotation=90 if len(labels) > 8 else 20,
                ha="right", fontsize=7,
            )
            for bar, value in zip(movement_bars, movements):
                if np.isfinite(value):
                    offset = 0.025 * movement_limit
                    axes[1].text(
                        bar.get_x() + bar.get_width() / 2,
                        value + (offset if value >= 0 else -offset),
                        f"{value:+.3f}", ha="center",
                        va="bottom" if value >= 0 else "top", fontsize=7,
                    )
            for ax in axes:
                ax.grid(axis="y", alpha=0.2)
            fig.suptitle(
                "Observed CK geometry-gate behavior",
                fontsize=14, fontweight="bold", y=0.98,
            )
            fig.subplots_adjust(
                left=0.12, right=0.98, bottom=0.20, top=0.88, hspace=0.38)
            save_figure(fig, gate_stem)
    else:
        if gate_table.is_file() or gate_table.is_symlink():
            gate_table.unlink()
        skip_figure(
            gate_stem,
            "no geometry-gate audit data were available for this selection",
        )

    # Combine per-task contam.R status without making R plots a gate for Python.
    status_rows = []
    for status_path in sorted((plot_root / "contam_r" / "status").glob("*.tsv")):
        with status_path.open("r", encoding="utf-8", errors="replace", newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                try:
                    selected_library = int(str(
                        row.get("library", "")).removeprefix("lib"))
                except ValueError:
                    continue
                if (selected_library in set(libraries) and
                        str(row.get("series_key", "")) in set(conditions)):
                    status_rows.append(row)
    _ambient_atomic_tsv(
        data_dir / "contam_r_status.tsv",
        ["library", "condition", "assignment_source", "series_key",
         "status", "reason", "source_prefix"],
        status_rows,
    )

    _ambient_atomic_tsv(
        data_dir / "plot_decisions.tsv",
        ["figure", "status", "reason"],
        plot_decisions,
    )
    _ambient_atomic_tsv(
        data_dir / "plot_manifest.tsv",
        ["path"],
        ({"path": path} for path in sorted(generated)),
    )
    if four_arm_mode and four_arm_hard_failures:
        unique_failures = _ambient_unique(four_arm_hard_failures)
        raise RuntimeError(
            "four-arm AMBIENT_PLOTS validation failed: "
            + "; ".join(unique_failures[:12])
            + (f"; plus {len(unique_failures) - 12} more"
               if len(unique_failures) > 12 else "")
        )
    print(f"AMBIENT_PLOTS wrote {len(generated)} figure files under {plot_root}")


# =============================================================================
# Diagnostics
# =============================================================================

def diagnose(libs, cond_list, condition_selection_label):
    """Print per-library status report for all requested conditions.

    Shows which upstream files exist and which conditions are ready/done/blocked.
    """
    print("=" * 72)
    print("PIPELINE DIAGNOSTICS")
    print(f"  Libraries: {len(libs)}")
    print(f"  Conditions: {len(cond_list)} ({condition_selection_label})")
    print(f"  Date: {datetime.now().strftime('%Y-%m-%d %H:%M')}")
    print("=" * 72)
    print()

    # CONDF stage: check condf files
    print("--- Stage 1 CONDF: .condf files ---")
    for key, path in CONDF_PATHS.items():
        status = "✅" if check_file_exists(path) else "❌"
        print(f"  {status} {key}: {path}")
    print()

    # Binary availability
    print("--- Binaries ---")
    for bname in ["demux_parallel", "tet_ambient_profile", "tet_contam_estimate"]:
        bpath = os.path.join(SOFTWARE_BIN, bname)
        status = "✅" if os.path.isfile(bpath) else "❌"
        print(f"  {status} {bname}: {bpath}")
    print()

    # Per-library status
    print("--- Per-library status ---")
    total_ready = 0
    total_done = 0
    total_blocked = 0

    for lib_num in libs:
        prefix = get_demux_prefix(lib_num)
        # Diagnosis is intentionally read-only: it reports the canonical file
        # if present but never regenerates expected-line metadata.
        el = get_expected_lines(lib_num, create=False)
        fb = get_filtered_barcodes(lib_num)

        # Check basic demux outputs
        has_counts = check_file_exists(prefix + ".counts")
        has_assignments = check_file_exists(prefix + ".assignments")
        has_samples = check_file_exists(prefix + ".samples")
        has_condf = check_file_exists(prefix + ".condf")
        has_el = el is not None
        has_fb = check_file_exists(fb)
        has_sp_counts = check_file_exists(prefix + ".species_counts")
        candidate_path = get_individual_ambient_candidates(lib_num, create=False)
        has_candidates = bool(candidate_path and check_file_exists(candidate_path))
        has_ed_ind = check_file_exists(get_empty_drops_indiv(lib_num))
        has_ed_sp = check_file_exists(get_empty_drops_species(lib_num))

        demux_ok = has_counts and has_assignments and has_samples and has_condf

        # Per-condition status
        n_done = 0
        n_ready = 0
        n_blocked = 0
        blocked_reasons = {}

        for cond in cond_list:
            abbrev = cond["abbrev"]

            # Check if output already exists
            if check_output_exists(lib_num, abbrev):
                n_done += 1
                continue

            # Check prerequisites
            ready, missing = check_upstream_ready(
                lib_num, cond, create_derived=False)
            if ready:
                n_ready += 1
            else:
                n_blocked += 1
                for m in missing:
                    # Categorize the blocker
                    if "empty_drops" in m:
                        blocked_reasons.setdefault("empty drops", set()).add(abbrev)
                    elif "ambient candidates" in m:
                        blocked_reasons.setdefault("ambient candidates", set()).add(abbrev)
                    elif "species_counts" in m or "species_condf" in m:
                        blocked_reasons.setdefault("species data", set()).add(abbrev)
                    elif "expected_lines" in m:
                        blocked_reasons.setdefault("expected_lines", set()).add(abbrev)
                    elif "panel_metadata" in m:
                        blocked_reasons.setdefault("panel_metadata", set()).add(abbrev)
                    else:
                        blocked_reasons.setdefault("demux outputs", set()).add(abbrev)

        total_ready += n_ready
        total_done += n_done
        total_blocked += n_blocked

        # Format upstream file status
        files_status = []
        files_status.append("D" if demux_ok else "d")
        files_status.append("E" if has_el else "e")
        files_status.append("F" if has_fb else "f")
        files_status.append("S" if has_sp_counts else "s")
        files_status.append("C" if has_candidates else "c")
        files_status.append("I" if has_ed_ind else "i")
        files_status.append("P" if has_ed_sp else "p")
        file_str = "".join(files_status)

        # Summary line
        parts = []
        if n_done > 0:
            parts.append(f"{n_done} done")
        if n_ready > 0:
            parts.append(f"{n_ready} ready")
        if n_blocked > 0:
            parts.append(f"{n_blocked} blocked")
        status_str = ", ".join(parts) if parts else "0 ready"

        print(f"  Lib {lib_num:>2} [{file_str}]: {status_str}")

        # Show blockers for small runs
        if blocked_reasons and len(libs) <= 10:
            for reason, conds in sorted(blocked_reasons.items()):
                print(f"         blocked by {reason}: {', '.join(sorted(conds))}")

    print()
    print(f"  File status key: D=demux E=expected_lines F=filtered_barcodes "
          f"S=species_counts C=ambient_candidates I=ed_indiv P=ed_species")
    print(f"  (lowercase = missing)")
    print()
    print(f"  Total: {total_done} done, {total_ready} ready, "
          f"{total_blocked} blocked")
    print(f"  across {len(libs)} libraries x {len(cond_list)} conditions "
          f"= {len(libs) * len(cond_list)} total cells")
    print("=" * 72)


# =============================================================================
# Orchestration
# =============================================================================

def parse_library_range(spec):
    """Parse a library specification like '1-40' or '1 2 8 19' into a list."""
    libs = []
    for item in spec:
        if "-" in item:
            parts = item.split("-", 1)
            start = int(parts[0])
            end = int(parts[1])
            libs.extend(range(start, end + 1))
        else:
            libs.append(int(item))
    return sorted(set(libs))


def _candidate_axis_preflight(lib_nums, args):
    failures = []
    targeted_event = bool(args.identity_candidate_axis_event_id)
    targeted_proposal = bool(args.identity_candidate_axis_proposal)
    if targeted_event != targeted_proposal:
        failures.append(
            "--identity-candidate-axis-event-id and proposal must be supplied together")
    if targeted_event and len(lib_nums) != 1:
        failures.append(
            "targeted candidate-axis mode requires exactly one requested library")
    raw_root = args.identity_candidate_axis_output_root or ""
    if not raw_root or not os.path.isabs(raw_root):
        failures.append(
            "--identity-candidate-axis-output-root must be an explicit absolute path")
    if (args.identity_candidate_axis_v6_3_root and
            not os.path.isabs(args.identity_candidate_axis_v6_3_root)):
        failures.append("--identity-candidate-axis-v6-3-root must be absolute")

    if raw_root and os.path.isabs(raw_root):
        paths = _candidate_axis_paths(args)
        for label in ("pairs", "scorer", "aggregate"):
            if os.path.lexists(paths[label]):
                failures.append(
                    "candidate-axis scientific output destination already exists; "
                    f"use a new output root: {paths[label]}")
        parent = Path(paths["root"]).parent
        while not parent.exists() and parent != parent.parent:
            parent = parent.parent
        if not parent.is_dir() or not os.access(parent, os.W_OK | os.X_OK):
            failures.append(
                f"candidate-axis output parent is not writable/searchable: {parent}")

    required = {
        "validation summary": os.path.join(
            get_identity_subdir(args, "validation"), "validation_summary.tsv"),
        "metadata expected genotypes": os.path.join(
            get_identity_subdir(args, "metadata"), "library_expected_genotypes.tsv"),
        "metadata biological lines": os.path.join(
            get_identity_subdir(args, "metadata"), "global_biological_lines.tsv"),
        "metadata donors": os.path.join(
            get_identity_subdir(args, "metadata"), "global_donors.tsv"),
        "pair builder": _identity_helper("identity_reconciliation.py"),
        "candidate-axis aggregate helper": _identity_helper(
            "identity_candidate_axis_aggregate.py"),
        "V6.3 source": _identity_helper(
            "identity_probability_aggregate.py"),
    }
    for label, path in required.items():
        if not check_file_exists(path):
            failures.append(f"candidate-axis {label} missing: {path}")

    event_presence = {}
    for lib_num in lib_nums:
        library = f"lib{lib_num}"
        prefix = get_demux_prefix(lib_num)
        library_required = {
            "decision table": os.path.join(
                get_identity_subdir(args, "decisions"),
                f"{library}.reconciled_cells.tsv.gz"),
            "nuclear samples": prefix + ".samples",
            "nuclear pileup sites": prefix + ".pileup_sites.tsv.gz",
            "nuclear pileup observations": prefix + ".pileup_obs.tsv.gz",
        }
        for label, path in library_required.items():
            if not check_file_exists(path):
                failures.append(f"candidate-axis {library} {label} missing: {path}")
        decision_path = library_required["decision table"]
        if check_file_exists(decision_path):
            try:
                event_presence[lib_num] = _candidate_axis_has_finalized_event(
                    lib_num, args
                )
            except Exception as exc:
                failures.append(
                    f"candidate-axis {library} finalized-event scan failed: {exc}")
                event_presence[lib_num] = False
    validation = required["validation summary"]
    if check_file_exists(validation):
        valid, detail = identity_validation_summary_passes(validation)
        if not valid:
            failures.append(
                "candidate-axis requires an all-PASS reconciliation validation "
                f"summary: {validation} ({detail})")

    scorer = os.path.join(SOFTWARE_BIN, "tetra_score_calls")
    if not os.path.isfile(scorer) or not os.access(scorer, os.X_OK):
        failures.append(f"candidate-axis scorer missing/not executable: {scorer}")
    else:
        try:
            capability = subprocess.run(
                [scorer, "--help"], capture_output=True, text=True,
                check=False, timeout=30)
            text = (capability.stdout or "") + (capability.stderr or "")
            options = (
                "--candidate-axis-output", "--candidate-axis-temp-dir",
                "--candidate-axis-self-test",
            )
            missing = [option for option in options if option not in text]
            if missing:
                failures.append(
                    "candidate-axis scorer help is missing: " + ", ".join(missing))
            selftest = subprocess.run(
                [scorer, "--candidate-axis-self-test"],
                capture_output=True, text=True, check=False, timeout=60)
            if selftest.returncode != 0:
                failures.append(
                    "candidate-axis scorer self-test failed: " +
                    (selftest.stderr.strip() or selftest.stdout.strip() or
                     f"exit {selftest.returncode}"))
        except (OSError, subprocess.SubprocessError) as exc:
            failures.append(
                f"could not verify candidate-axis scorer capabilities: {exc}")

    aggregate_helper = required["candidate-axis aggregate helper"]
    if check_file_exists(aggregate_helper):
        try:
            selftest = subprocess.run(
                [sys.executable, "-B", aggregate_helper, "--self-test"],
                capture_output=True, text=True, check=False, timeout=60)
            if selftest.returncode != 0:
                failures.append(
                    "candidate-axis aggregate self-test failed: " +
                    (selftest.stderr.strip() or selftest.stdout.strip() or
                     f"exit {selftest.returncode}"))
        except (OSError, subprocess.SubprocessError) as exc:
            failures.append(
                f"could not run candidate-axis aggregate self-test: {exc}")

    for lib_num in lib_nums:
        if lib_num not in event_presence:
            continue
        try:
            _candidate_axis_frozen_parameters(
                args, lib_num, event_presence[lib_num]
            )
        except Exception as exc:
            failures.append(
                f"candidate-axis lib{lib_num} frozen scoring parameters invalid: {exc}")
    temp_root = Path("/dev/shm")
    if (not temp_root.is_absolute() or not temp_root.is_dir() or
            not os.access(temp_root, os.W_OK | os.X_OK)):
        failures.append(
            "candidate-axis requires an absolute writable job-local /dev/shm; "
            "the score wrapper revalidates SLURM_TMPDIR or /dev/shm on the node")
    return len(failures) == 0, failures


def preflight_validate(lib_nums, cond_list, stages, args=None):
    """Check that all files needed by the requested stages exist on disk
    right now. Collects all failures before reporting.

    Note: shared memory segments (VCF daemon) are NOT checked here because
    they are node-local to pika/char and invisible from the head node.
    The sbatch scripts validate them at runtime on the compute node.

    Returns (ok, failures) where failures is a list of strings.
    """
    if set(stages) == {CANDIDATE_AXIS_STAGE}:
        return _candidate_axis_preflight(lib_nums, args)
    failures = []

    mt_manifest_libraries = None
    if (args is not None and
            ({"MT_FUSION", "MT_POPULATION"} & set(stages)) and
            args.mt_site_manifest and check_file_exists(args.mt_site_manifest)):
        try:
            mt_manifest_libraries = get_mt_manifest_library_numbers(args.mt_site_manifest)
        except Exception as e:
            failures.append(f"MT_FUSION site manifest could not be read: {e}")

    # 1. Binaries
    required_bins = []
    if "CONDF" in stages or "DEMUX" in stages:
        required_bins.append("demux_parallel")
        if args is not None and args.manage_vcf_daemons:
            required_bins.append("vcf_loader_daemon")
    if "EMPTY_DROPS" in stages:
        required_bins.append("tet_ambient_profile")
    if "CONTAM" in stages or "GEX_AMBIENT" in stages:
        required_bins.append("tet_contam_estimate")
    if {"AMBIENT_VALIDATE", "AMBIENT_SWAP_TEST"} & set(stages):
        required_bins.extend(["tet_ambient_profile", "tet_contam_estimate"])
    if "CONTAM" in stages and any(c["mode"] == LEGACY2C_MODE for c in cond_list):
        required_bins.append("legacy2c_contam_estimate")
    if "POSTHOC" in stages:
        required_bins.append("tetra_score_calls")
    if "TETRA_REFINE" in stages:
        required_bins.append("tetra_refine")
    if "MT_FUSION" in stages:
        required_bins.append("mt_fusion_ratio")
    if {"IDENTITY_RECONCILIATION", "IDENTITY_FINAL_EVIDENCE"} & set(stages):
        required_bins.extend([
            "tetra_score_calls", "mt_identity_score",
            "nuclear_panel_distinguishability", "tet_ambient_profile",
            "tet_contam_estimate"])
        if args.identity_evidence_mode == "rna-atac":
            required_bins.append("demux_parallel")
    if "IDENTITY_SCORE" in stages:
        required_bins.append("tetra_score_calls")
    for bname in required_bins:
        bpath = os.path.join(SOFTWARE_BIN, bname)
        if not os.path.isfile(bpath):
            failures.append(f"Binary missing: {bpath}")
        elif not os.access(bpath, os.X_OK):
            failures.append(f"Binary not executable: {bpath}")

    if "IDENTITY_SCORE" in stages:
        score_binary = os.path.join(SOFTWARE_BIN, "tetra_score_calls")
        if os.path.isfile(score_binary) and os.access(score_binary, os.X_OK):
            try:
                capability = subprocess.run(
                    [score_binary, "--help"], capture_output=True, text=True,
                    check=False, timeout=20)
                help_text = (
                    (capability.stdout or "") + (capability.stderr or ""))
                required_score_options = (
                    "--probability-output", "--pileup-molecules",
                    "--probability-resamples")
                missing = [
                    option for option in required_score_options
                    if option not in help_text]
                if missing:
                    failures.append(
                        "IDENTITY_SCORE requires the targeted nuclear-evidence "
                        "tetra_score_calls build; missing options: " +
                        ", ".join(missing))
                if "targeted original-vs-reconciliation-swap probability" not in help_text:
                    failures.append(
                        "IDENTITY_SCORE requires the reconciliation-targeted "
                        "tetra_score_calls build; the deployed binary still "
                        "accepts generic discovery comparisons")
            except (OSError, subprocess.SubprocessError) as exc:
                failures.append(
                    "could not verify tetra_score_calls probability "
                    f"capabilities: {exc}")
    if {"IDENTITY_RECONCILIATION", "IDENTITY_FINAL_EVIDENCE"} & set(stages):
        score_binary = os.path.join(SOFTWARE_BIN, "tetra_score_calls")
        if os.path.isfile(score_binary) and os.access(score_binary, os.X_OK):
            try:
                capability = subprocess.run(
                    [score_binary, "--help"], capture_output=True, text=True,
                    check=False, timeout=20)
                help_text = (
                    (capability.stdout or "") + (capability.stderr or ""))
                missing = [
                    option for option in (
                        "--candidate-axis-output",
                        "--candidate-axis-temp-dir",
                        "--candidate-axis-self-test")
                    if option not in help_text]
                if missing:
                    failures.append(
                        "IDENTITY_RECONCILIATION candidate-axis scorer is "
                        "missing options: " + ", ".join(missing))
            except (OSError, subprocess.SubprocessError) as exc:
                failures.append(
                    "could not verify tetra_score_calls candidate-axis "
                    f"capabilities: {exc}")
    if ("DEMUX" in stages and
            {"IDENTITY_RECONCILIATION", "IDENTITY_SCORE"} & set(stages)):
        demux_binary = os.path.join(SOFTWARE_BIN, "demux_parallel")
        if os.path.isfile(demux_binary) and os.access(demux_binary, os.X_OK):
            try:
                capability = subprocess.run(
                    [demux_binary, "--help"], capture_output=True, text=True,
                    check=False, timeout=20)
                help_text = (
                    (capability.stdout or "") + (capability.stderr or ""))
                if ".pileup_molecules.tsv.gz" not in help_text:
                    failures.append(
                        "Molecule-aware identity scoring after DEMUX requires "
                        "the supplied demux_parallel v2.15 build")
            except (OSError, subprocess.SubprocessError) as exc:
                failures.append(
                    "could not verify demux_parallel molecule-sidecar "
                    f"capability: {exc}")

    if "GEX_AMBIENT" in stages:
        estimator = os.path.join(SOFTWARE_BIN, "tet_contam_estimate")
        if os.path.isfile(estimator) and os.access(estimator, os.X_OK):
            try:
                help_result = subprocess.run(
                    [estimator, "--help"], capture_output=True, text=True,
                    check=False, timeout=15)
                help_text = (help_result.stdout or "") + (help_result.stderr or "")
                required_options = (
                    "--barcodes", "--features", "--matrix", "--clusts",
                    "--noround", "--skip_genes_regex")
                missing_options = [
                    option for option in required_options
                    if option not in help_text]
                if missing_options:
                    failures.append(
                        "GEX_AMBIENT requires a tet_contam_estimate build with "
                        "ambient_rna_gex support; missing help options: "
                        + ", ".join(missing_options))
            except (OSError, subprocess.SubprocessError) as exc:
                failures.append(
                    f"Could not verify tet_contam_estimate GEX capabilities: {exc}")

    if ("MT_FUSION" in stages and args is not None and
            args.mt_rna_ambient_condition):
        mt_binary = os.path.join(SOFTWARE_BIN, "mt_fusion_ratio")
        if os.path.isfile(mt_binary) and os.access(mt_binary, os.X_OK):
            try:
                capability = subprocess.run(
                    [mt_binary, "--help"], capture_output=True, text=True,
                    check=False, timeout=20)
                capability_text = (
                    (capability.stdout or "") + (capability.stderr or ""))
                if "--rna_ambient_fraction_file" not in capability_text:
                    failures.append(
                        "MT_FUSION RNA ambient integration requires the supplied "
                        "mt_fusion_ratio build with --rna_ambient_fraction_file")
            except (OSError, subprocess.SubprocessError) as exc:
                failures.append(
                    "could not verify mt_fusion_ratio RNA ambient capability: "
                    f"{exc}")

    four_arm_selected = (
        args is not None and
        args.contam_assignment_source == IDENTITY_AMBIENT_SELECTOR and
        bool({"CONTAM", "AMBIENT_PLOTS", "AMBIENT_VALIDATE"} & set(stages))
    )
    if four_arm_selected:
        unsupported = [
            condition["abbrev"] for condition in cond_list
            if condition["mode"] not in (1, 3, "1+SR")
        ]
        if unsupported:
            failures.append(
                "reconciliation-four-arm supports individual-native "
                "tet_contam_estimate conditions only; unsupported: "
                + ", ".join(unsupported))
    fixed_assignment_test_selected = bool(
        {"AMBIENT_VALIDATE", "AMBIENT_SWAP_TEST", "IDENTITY_RECONCILIATION"}
        & set(stages))
    # Plot-only reuse is validated from completed arm contracts below; it must
    # not depend on whichever estimator/helper is deployed after those outputs
    # were produced. Probe capabilities only when CONTAM will actually run.
    if ((four_arm_selected and "CONTAM" in stages) or
            fixed_assignment_test_selected):
        estimator = os.path.join(SOFTWARE_BIN, "tet_contam_estimate")
        if os.path.isfile(estimator) and os.access(estimator, os.X_OK):
            try:
                capability = subprocess.run(
                    [estimator, "--help"], capture_output=True, text=True,
                    check=False, timeout=20)
                capability_text = (
                    (capability.stdout or "") + (capability.stderr or ""))
                if "--freeze_assignments" not in capability_text:
                    failures.append(
                        "reconciliation-four-arm requires the supplied "
                        "tet_contam_estimate build with --freeze_assignments")
            except (OSError, subprocess.SubprocessError) as exc:
                failures.append(
                    "could not verify tet_contam_estimate fixed-assignment "
                    f"capability: {exc}")
        if fixed_assignment_test_selected:
            empty_binary = os.path.join(SOFTWARE_BIN, "tet_ambient_profile")
            if os.path.isfile(empty_binary) and os.access(empty_binary, os.X_OK):
                try:
                    capability = subprocess.run(
                        [empty_binary, "--help"], capture_output=True,
                        text=True, check=False, timeout=20)
                    capability_text = (
                        (capability.stdout or "") +
                        (capability.stderr or ""))
                    if "--profile_starts" not in capability_text:
                        failures.append(
                            "fixed-profile ambient validation requires the robust "
                            "bulk-profile "
                            "tet_ambient_profile build with --profile_starts")
                    version_probe = subprocess.run(
                        [empty_binary, "--version"], capture_output=True,
                        text=True, check=False, timeout=20)
                    version_text = (
                        (version_probe.stdout or "") +
                        (version_probe.stderr or ""))
                    if (version_probe.returncode != 0 or
                            f"v{AMBIENT_PROFILE_REQUIRED_VERSION}" not in
                            version_text):
                        failures.append(
                            "fixed-profile ambient validation requires "
                            "tet_ambient_profile "
                            f"v{AMBIENT_PROFILE_REQUIRED_VERSION}; deployed "
                            "binary is stale or incompatible")
                except (OSError, subprocess.SubprocessError) as exc:
                    failures.append(
                        "could not verify tet_ambient_profile bulk-profile "
                        f"capability: {exc}")
        if any(
                condition.get("runner") == "geometry_gated"
                for condition in cond_list):
            geometry_helper = resolve_geometry_gate_script()
            if os.path.isfile(geometry_helper):
                try:
                    capability = subprocess.run(
                        [sys.executable, geometry_helper, "--capabilities"],
                        capture_output=True, text=True, check=False,
                        timeout=20)
                    if (capability.returncode != 0 or
                            capability.stdout.strip() !=
                            "identity_ambient_fixed_assignments_fixed_profile_v2"):
                        detail = (
                            capability.stderr.strip() or
                            capability.stdout.strip() or
                            f"exit {capability.returncode}")
                        failures.append(
                            "reconciliation-four-arm geometry conditions require "
                            "the supplied endpoint-consistent geometry helper "
                            f"(capability probe failed: {detail})")
                except (OSError, subprocess.SubprocessError) as exc:
                    failures.append(
                        "could not verify geometry helper fixed-assignment "
                        f"capability: {exc}")

    if (args is not None and args.manage_vcf_daemons and
            ({"CONDF", "DEMUX"} & set(stages))):
        daemon_path = os.path.join(SOFTWARE_BIN, "vcf_loader_daemon")
        if os.path.isfile(daemon_path) and os.access(daemon_path, os.X_OK):
            try:
                help_result = subprocess.run(
                    [daemon_path, "--help"], capture_output=True, text=True,
                    check=False, timeout=15)
                help_text = (help_result.stdout or "") + (help_result.stderr or "")
                if "--ready-file" not in help_text:
                    failures.append(
                        "Managed VCF lifecycle requires the updated vcf_loader_daemon "
                        "with --ready-file support; rebuild/deploy the supplied source")
            except (OSError, subprocess.SubprocessError) as exc:
                failures.append(
                    f"Could not verify managed vcf_loader_daemon capabilities: {exc}")
        for key, path in VCF_SOURCE_PATHS.items():
            if not check_file_exists(path):
                failures.append(f"Managed VCF daemon source missing ({key}): {path}")
        if not get_vcf_daemon_nodes():
            failures.append(
                "Managed VCF lifecycle requires non-empty --daemon-nodes")
        if resolve_vcf_daemon_reference_bam(lib_nums) is None:
            failures.append(
                "Managed VCF lifecycle could not find a selected library BAM "
                "for chromosome/TID mapping")

    # 1b. Sidecar Python scripts used by POSTHOC/HYBRID.
    if "POSTHOC" in stages:
        for helper in required_posthoc_scripts():
            if not os.path.isfile(helper):
                failures.append(f"POSTHOC helper script missing: {helper}")
    if "POSTHOC_SUMMARY" in stages:
        for helper_name in ("swap_audit_summarize.py", "swap_audit_aggregate.py"):
            helper = resolve_process_script(helper_name)
            if not os.path.isfile(helper):
                failures.append(f"POSTHOC_SUMMARY helper script missing: {helper}")
    if "UNEXPECTED_COMPONENT_NN" in stages:
        helper = resolve_process_script("analyze_unexpected_component_nn.py")
        if not os.path.isfile(helper):
            failures.append(f"UNEXPECTED_COMPONENT_NN helper script missing: {helper}")
    if "HYBRID" in stages:
        for helper in required_hybrid_scripts():
            if not os.path.isfile(helper):
                failures.append(f"HYBRID helper script missing: {helper}")
    if {"AMBIENT_PLOTS", "IDENTITY_RECONCILIATION",
            "IDENTITY_FINAL_EVIDENCE"} & set(stages):
        contam_r = resolve_contam_r_script()
        if not os.path.isfile(contam_r):
            failures.append(f"AMBIENT_PLOTS deployed contam.R missing: {contam_r}")
    if {"AMBIENT_SWAP_TEST", "IDENTITY_RECONCILIATION",
            "IDENTITY_FINAL_EVIDENCE",
            "IDENTITY_FINAL_EVIDENCE_ONLY"} & set(stages):
        if not os.path.isfile(IDENTITY_RECONCILIATION_FIGURES):
            failures.append(
                "AMBIENT_SWAP_TEST figure generator missing: "
                + IDENTITY_RECONCILIATION_FIGURES)
    if ({"IDENTITY_RECONCILIATION", "IDENTITY_FINAL_EVIDENCE",
         "IDENTITY_RECONCILE_ONLY",
         "IDENTITY_SCORE", "IDENTITY_SCORE_AGGREGATE_ONLY"} & set(stages)):
        for helper_name in (
                "identity_reconciliation.py", "identity_reconciliation_common.py",
                "validate_identity_reconciliation.py",
                "identity_probability_aggregate.py"):
            helper = resolve_process_script(helper_name)
            if not os.path.isfile(helper):
                failures.append(
                    f"Identity reconciliation helper script missing: {helper}")
    if {"IDENTITY_FINAL_EVIDENCE_ONLY",
            "IDENTITY_FINALIZE_ONLY"} & set(stages):
        stage_label = (
            "IDENTITY_FINALIZE_ONLY"
            if "IDENTITY_FINALIZE_ONLY" in stages
            else "IDENTITY_FINAL_EVIDENCE_ONLY")
        for helper_name in (
                "identity_reconciliation.py",
                "identity_reconciliation_common.py"):
            helper = resolve_process_script(helper_name)
            if not os.path.isfile(helper):
                failures.append(
                    f"{stage_label} helper missing: {helper}")
    if {"IDENTITY_RECONCILIATION", "IDENTITY_FINAL_EVIDENCE"} & set(stages):
        if not check_file_exists(IDENTITY_LIBRARY_CONVERSIONS_XLSX):
            failures.append(
                "IDENTITY_RECONCILIATION global/library metadata workbook missing: "
                + IDENTITY_LIBRARY_CONVERSIONS_XLSX)
        if not check_file_exists(VCF_SOURCE_PATHS["interindiv_20M"]):
            failures.append(
                "IDENTITY_RECONCILIATION NoMito nuclear panel missing: "
                + VCF_SOURCE_PATHS["interindiv_20M"])
        for helper_name in (
                "identity_candidate_axis_aggregate.py",
                "identity_reconciliation.py"):
            helper = resolve_process_script(helper_name)
            if not os.path.isfile(helper):
                failures.append(
                    f"IDENTITY_RECONCILIATION helper missing: {helper}")
        if args.identity_review_input and not check_file_exists(
                args.identity_review_input):
            failures.append(
                "IDENTITY_RECONCILIATION review input missing: "
                + args.identity_review_input)
    if {"IDENTITY_SCORE", "IDENTITY_SCORE_AGGREGATE_ONLY"} & set(stages):
        metadata_root = get_identity_subdir(args, "metadata")
        for path in (
                os.path.join(metadata_root, "library_expected_genotypes.tsv"),
                os.path.join(metadata_root, "global_biological_lines.tsv"),
                os.path.join(metadata_root, "global_donors.tsv")):
            if not check_file_exists(path):
                failures.append(
                    f"identity score existing metadata input missing: {path}")
        validation = os.path.join(
            get_identity_subdir(args, "validation"),
            "validation_summary.tsv")
        validation_ok, validation_detail = identity_validation_summary_passes(
            validation
        )
        if not validation_ok:
            failures.append(
                "identity scoring requires an all-PASS reconciliation "
                f"validation summary: {validation} ({validation_detail})")
    if "IDENTITY_RECONCILE_ONLY" in stages:
        metadata_root = get_identity_subdir(args, "metadata")
        for path in (
            os.path.join(metadata_root, "library_expected_genotypes.tsv"),
            os.path.join(metadata_root, "global_biological_lines.tsv"),
            os.path.join(metadata_root, "global_donors.tsv"),
            os.path.join(metadata_root, "nuclear_panel_distinguishability.tsv"),
        ):
            if not check_file_exists(path):
                failures.append(f"IDENTITY_RECONCILE_ONLY existing metadata/evidence input missing: {path}")
    if {"IDENTITY_FINAL_EVIDENCE", "IDENTITY_FINAL_EVIDENCE_ONLY",
            "IDENTITY_FINALIZE_ONLY"} & set(stages):
        stage_label = (
            "IDENTITY_FINALIZE_ONLY"
            if "IDENTITY_FINALIZE_ONLY" in stages else
            "IDENTITY_FINAL_EVIDENCE_ONLY"
            if "IDENTITY_FINAL_EVIDENCE_ONLY" in stages
            else "IDENTITY_FINAL_EVIDENCE")
        validation = os.path.join(
            get_identity_subdir(args, "validation"),
            "validation_summary.tsv")
        validation_ok, validation_detail = identity_validation_summary_passes(
            validation)
        if not validation_ok:
            failures.append(
                f"{stage_label} requires an all-PASS reconciliation "
                f"validation summary: {validation} ({validation_detail})")
        for path in (
                get_identity_event_path(args),
                get_identity_reconciliation_manifest_path(args),
                os.path.join(
                    get_identity_subdir(args, "decisions"),
                    "all_libraries.reconciled_cells.tsv.gz")):
            if not check_file_exists(path):
                failures.append(
                    f"{stage_label} input missing: {path}")
        if (args.identity_review_input and
                not check_file_exists(args.identity_review_input)):
            failures.append(
                f"{stage_label} review input missing: "
                + args.identity_review_input)
        if "IDENTITY_FINALIZE_ONLY" in stages:
            try:
                resolve_identity_finalize_only_evidence_roots(
                    lib_nums, args, require_completed=True)
            except (FileNotFoundError, OSError, ValueError) as exc:
                failures.append(
                    "IDENTITY_FINALIZE_ONLY evidence input invalid: "
                    + str(exc))

    if "MT_FUSION" in stages:
        for path, label in [
            (args.mt_vcf, "MT_FUSION mitochondrial VCF/BCF"),
            (args.mt_site_manifest, "MT_FUSION site manifest"),
        ]:
            if not path or not check_file_exists(path):
                failures.append(f"{label} missing: {path}")
        if args.mt_site_calibration and not check_file_exists(args.mt_site_calibration):
            failures.append(f"MT_FUSION site calibration missing: {args.mt_site_calibration}")
        if args.mt_rho_reference and not check_file_exists(args.mt_rho_reference):
            failures.append(f"MT_FUSION rho reference missing: {args.mt_rho_reference}")
        if args.mt_mask_bed and not check_file_exists(args.mt_mask_bed):
            failures.append(f"MT_FUSION mt mask BED missing: {args.mt_mask_bed}")
        if (
                args.mt_rho_mode in ("fixed", "low_information_fixed") and
                args.mt_pooled_rho is None and not args.mt_rho_reference):
            failures.append(
                "MT_FUSION fixed/low_information_fixed rho mode requires "
                "--mt-pooled-rho or --mt-rho-reference")
        ambient_choices = sum(bool(x) for x in (
            args.mt_ambient_profile_template, args.mt_empty_barcodes_template,
            args.mt_ambient_none))
        if ambient_choices > 1:
            failures.append(
                "MT_FUSION accepts only one ambient mode: --mt-ambient-profile-template, "
                "--mt-empty-barcodes-template, or --mt-ambient-none")
        if args.mt_assay_mode == "ATAC" and not args.mt_bam_template:
            failures.append("MT_FUSION ATAC mode requires --mt-bam-template so the RNA BAM is not used accidentally")
    if "MT_POPULATION" in stages:
        for helper in required_mt_population_scripts():
            if not os.path.isfile(helper):
                failures.append(f"MT_POPULATION helper script missing: {helper}")
        for lib_num in lib_nums:
            if (mt_manifest_libraries is not None and
                    lib_num not in mt_manifest_libraries):
                continue
            reconciled = get_reconciled_cells_path(lib_num)
            reconciliation_generated_upstream = (
                "IDENTITY_RECONCILIATION" in stages or
                "IDENTITY_RECONCILE_ONLY" in stages)
            if (not reconciliation_generated_upstream and
                    not check_file_exists(reconciled)):
                failures.append(
                    f"lib{lib_num}: MT_POPULATION reconciled-cell input missing: {reconciled}")
            if "MT_FUSION" not in stages:
                if (args.mt_rna_ambient_condition and
                        not mt_outputs_complete(lib_num, args)):
                    failures.append(
                        f"lib{lib_num}: MT_POPULATION inputs are incomplete or "
                        "have stale RNA ambient provenance")
                else:
                    for path in (
                            get_mt_ratio_path(lib_num, args),
                            get_mt_profile_path(lib_num, args)):
                        if not check_file_exists(path):
                            failures.append(
                                f"lib{lib_num}: MT_POPULATION input missing: {path}")
    if "CONTAM" in stages and any(
        c.get("runner") == "geometry_gated" for c in cond_list
    ):
        for helper in required_geometry_gate_scripts():
            if not os.path.isfile(helper):
                failures.append(f"Geometry-gate helper script missing: {helper}")
    if "PLOIDY_NN" in stages:
        for helper in required_ploidy_nn_scripts():
            if not os.path.isfile(helper):
                failures.append(f"PLOIDY_NN helper script missing: {helper}")
        for path, label in [
            (args.ploidy_nn_h5ad, "PLOIDY_NN h5ad"),
            (args.ploidy_nn_weights, "PLOIDY_NN weights"),
        ]:
            if not path or not os.path.isfile(path):
                failures.append(f"{label} missing: {path}")
        if args.ploidy_nn_weights:
            scaler = args.ploidy_nn_weights[:-3] + "_scaler.npz" if args.ploidy_nn_weights.endswith(".pt") else args.ploidy_nn_weights + "_scaler.npz"
            if not os.path.isfile(scaler):
                failures.append(f"PLOIDY_NN scaler missing: {scaler}")
    if "TETRA_REFINE" in stages:
        for helper in required_tetra_refine_scripts():
            if not os.path.isfile(helper):
                failures.append(f"TETRA_REFINE helper script missing: {helper}")

    if "GEX_AMBIENT" in stages:
        if args.gex_cluster_source == "manual" and not args.gex_clusters_template:
            failures.append(
                "manual GEX_AMBIENT clustering requires --gex-clusters-template")
        if (args.gex_cluster_source == "h5ad" and
                not check_file_exists(args.ploidy_nn_h5ad)):
            failures.append(
                "GEX_AMBIENT NN H5AD missing: "
                f"{args.ploidy_nn_h5ad}")
        if args.gex_marker_genes and not check_file_exists(args.gex_marker_genes):
            failures.append(
                f"GEX_AMBIENT marker-gene file missing: {args.gex_marker_genes}")
        if args.gex_min_cluster_cells < 2:
            failures.append("--gex-min-cluster-cells must be at least 2")

    # 2. Panel metadata is producer/orchestrator-side metadata.
    # DEMUX needs it to write native .species_* artifacts; EMPTY_DROPS and
    # native-species CONTAM need it to generate/consume species-native
    # expected_lines and run species-native profile estimation.
    needs_panel = False
    if ("DEMUX" in stages and not args.individual_only_demux) or "EMPTY_DROPS" in stages:
        needs_panel = True
    if "CONTAM" in stages:
        for cond in cond_list:
            if cond["mode"] in (2, 4, 5, "1+SR"):
                needs_panel = True
                break
    if ("POSTHOC" in stages or "POSTHOC_SUMMARY" in stages or
            "IDENTITY_RECONCILIATION" in stages or
            "IDENTITY_FINAL_EVIDENCE" in stages or
            "IDENTITY_RECONCILE_ONLY" in stages):
        needs_panel = True
    if needs_panel and not check_file_exists(PANEL_METADATA):
        failures.append(f"Panel metadata missing: {PANEL_METADATA}")

    needs_expected_lines = bool(
        {"DEMUX", "EMPTY_DROPS", "CONTAM", "TETRA_REFINE",
         "IDENTITY_RECONCILIATION", "IDENTITY_FINAL_EVIDENCE",
         "AMBIENT_SWAP_TEST"} & set(stages))
    if needs_expected_lines and not check_file_exists(LIBRARY_CONVERSIONS_XLSX):
        failures.append(
            f"Library conversions workbook missing: {LIBRARY_CONVERSIONS_XLSX}")

    ambient_swap_discovery = None
    if "AMBIENT_SWAP_TEST" in stages:
        try:
            ambient_swap_discovery = discover_ambient_swap_events(args)
        except Exception as exc:
            failures.append(f"AMBIENT_SWAP_TEST discovery failed: {exc}")

    # 3. Per-library files
    for lib_num in lib_nums:
        # BAM (needed by DEMUX)
        if "DEMUX" in stages:
            bam = get_bam_path(lib_num)
            if bam is None:
                failures.append(f"lib{lib_num}: BAM file missing in {get_lib_dir(lib_num)}")

        # Filtered barcodes (needed by DEMUX and EMPTY_DROPS)
        if "DEMUX" in stages or "EMPTY_DROPS" in stages:
            fb = get_filtered_barcodes(lib_num)
            if not check_file_exists(fb):
                failures.append(f"lib{lib_num}: Filtered barcodes missing: {fb}")

        # Expected lines are needed by DEMUX (-I), EMPTY_DROPS (--ids/-i),
        # and CONTAM (-X). Native species EMPTY_DROPS and CONTAM receive a
        # generated species-resolution expected-lines file.
        if ({"DEMUX", "EMPTY_DROPS", "CONTAM", "AMBIENT_SWAP_TEST",
             "IDENTITY_FINAL_EVIDENCE"}
                & set(stages)):
            try:
                el = get_expected_lines(lib_num)
            except Exception as e:
                el = None
                failures.append(
                    f"lib{lib_num}: Expected lines generation failed: {e}")
            if el is None:
                failures.append(f"lib{lib_num}: Expected lines file missing")
            needs_species_el = "EMPTY_DROPS" in stages or (
                "CONTAM" in stages and any(c["mode"] in (2, 4, 5) for c in cond_list)
            )
            if needs_species_el and el is not None:
                try:
                    get_species_expected_lines(lib_num)
                except Exception as e:
                    failures.append(f"lib{lib_num}: Species expected lines generation failed: {e}")
            needs_individual_candidates = (
                ("CONTAM" in stages or "AMBIENT_SWAP_TEST" in stages)
                and any(c["mode"] in (1, 3, "1+SR", LEGACY2C_MODE) for c in cond_list)
            )
            if needs_individual_candidates and el is not None:
                try:
                    candidate_path = get_individual_ambient_candidates(lib_num)
                    if not check_file_exists(candidate_path):
                        failures.append(
                            f"lib{lib_num}: Individual ambient candidates missing: {candidate_path}")
                except Exception as e:
                    failures.append(
                        f"lib{lib_num}: Individual ambient candidate generation failed: {e}")

        if four_arm_selected:
            try:
                comparison_context = prepare_identity_ambient_comparison(
                    lib_num, candidate_set=args.reconciliation_candidate_set)
                four_arm_analysis_requested = bool(
                    {"CONTAM", "AMBIENT_PLOTS"} & set(stages))
                if (four_arm_analysis_requested and
                        int(comparison_context.get("n_scrutinized", 0)) == 0):
                    failures.append(
                        f"lib{lib_num}: four-arm plan has no scrutinized cells")
                if (four_arm_analysis_requested and
                        int(comparison_context.get("n_background", 0)) == 0):
                    failures.append(
                        f"lib{lib_num}: four-arm plan has no background cells")
            except Exception as exc:
                failures.append(
                    f"lib{lib_num}: four-arm reconciliation plan failed: {exc}")

        if "AMBIENT_VALIDATE" in stages:
            try:
                validation_plan = prepare_ambient_validation_plan(
                    lib_num, args.reconciliation_candidate_set)
            except Exception as exc:
                failures.append(
                    f"lib{lib_num}: AMBIENT_VALIDATE plan failed: {exc}")
                validation_plan = None
            if validation_plan and validation_plan.get("applicable"):
                raw_prefix = os.path.join(
                    get_demux_dir(lib_num), f"lib{lib_num}_raw")
                for path in (
                        raw_prefix + ".counts", raw_prefix + ".condf",
                        get_demux_prefix(lib_num) + ".samples",
                        get_filtered_barcodes(lib_num)):
                    if not check_file_exists(path):
                        failures.append(
                            f"lib{lib_num}: AMBIENT_VALIDATE input missing: {path}")

        if ("AMBIENT_SWAP_TEST" in stages and
                ambient_swap_discovery is not None):
            try:
                swap_plan = prepare_ambient_swap_test_plan(
                    lib_num, args, discovery=ambient_swap_discovery)
            except Exception as exc:
                failures.append(
                    f"lib{lib_num}: AMBIENT_SWAP_TEST plan failed: {exc}")
                swap_plan = None
            if swap_plan and swap_plan.get("applicable"):
                for path in (
                        swap_plan["raw_prefix"] + ".counts",
                        swap_plan["raw_prefix"] + ".condf",
                        swap_plan["demux_prefix"] + ".samples",
                        get_filtered_barcodes(lib_num)):
                    if not check_file_exists(path):
                        failures.append(
                            f"lib{lib_num}: AMBIENT_SWAP_TEST input missing: {path}")

        if "CONTAM" in stages:
            for assignment_source in resolve_contam_assignment_sources(
                    args.contam_assignment_source):
                if assignment_source == "reconciled_replacement":
                    try:
                        applicable, _ = contam_source_applicable(
                            lib_num, assignment_source)
                    except Exception:
                        applicable = False
                    if not applicable:
                        continue
                assignments = get_contam_assignments_path(
                    lib_num, assignment_source)
                generated_upstream = (
                    assignment_source == "demux" and "DEMUX" in stages)
                if not generated_upstream and not check_file_exists(assignments):
                    failures.append(
                        f"lib{lib_num}: CONTAM {assignment_source} assignments missing: "
                        f"{assignments}")

        if "GEX_AMBIENT" in stages:
            for path, label in (
                    (get_filtered_barcodes(lib_num), "barcodes"),
                    (get_filtered_features(lib_num), "features"),
                    (get_filtered_matrix(lib_num), "matrix")):
                if not check_file_exists(path):
                    failures.append(
                        f"lib{lib_num}: GEX_AMBIENT filtered {label} missing: {path}")
            try:
                clusters = get_gex_cluster_path(lib_num, args)
            except Exception as exc:
                failures.append(
                    f"lib{lib_num}: GEX cluster configuration invalid: {exc}")
                clusters = None
            if (args.gex_cluster_source == "manual" and clusters and
                    not check_file_exists(clusters)):
                failures.append(
                    f"lib{lib_num}: GEX_AMBIENT cluster file missing: {clusters}")
            if "CONTAM" not in stages:
                for cond in cond_list:
                    if cond["mode"] != 1:
                        continue
                    for assignment_source in resolve_contam_assignment_sources(
                            args.contam_assignment_source):
                        try:
                            applicable, _ = contam_source_applicable(
                                lib_num, assignment_source)
                        except Exception:
                            applicable = False
                        if not applicable:
                            continue
                        source_prefix = get_contam_prefix(
                            lib_num, cond["abbrev"], assignment_source)
                        for suffix in (".contam_rate", ".contam_prof",
                                       ".assignments", ".samples"):
                            path = source_prefix + suffix
                            if not check_file_exists(path):
                                failures.append(
                                    f"lib{lib_num}: GEX_AMBIENT source missing: {path}")

        if "TETRA_REFINE" in stages:
            prefix = get_demux_prefix(lib_num)
            refine_needed = [prefix + ext for ext in [".assignments", ".diagnostics.gz", ".runner_ups.gz"]]
            if "DEMUX" not in stages:
                for f in refine_needed:
                    if not check_file_exists(f):
                        failures.append(f"lib{lib_num}: TETRA_REFINE input missing: {f}")
            if get_expected_lines(lib_num) is None:
                failures.append(f"lib{lib_num}: TETRA_REFINE expected lines missing")
            if args.require_refine_external_ploidy and "PLOIDY_NN" not in stages:
                pc = get_ploidy_calls_path(lib_num)
                if not check_file_exists(pc):
                    failures.append(f"lib{lib_num}: required external ploidy calls missing: {pc}")
            refine_contam_planned = (
                "CONTAM" in stages and args.refine_contam_condition in {
                    condition["abbrev"] for condition in cond_list})
            if args.require_refine_contam_rate and not refine_contam_planned:
                cr = get_contam_prefix(lib_num, args.refine_contam_condition) + ".contam_rate"
                if not check_file_exists(cr):
                    failures.append(f"lib{lib_num}: required refine contam_rate missing: {cr}")

        if "AMBIENT_PLOTS" in stages and "CONTAM" not in stages:
            for cond in cond_list:
                for assignment_source in resolve_contam_assignment_sources(
                        args.contam_assignment_source):
                    if assignment_source == "reconciled_replacement":
                        try:
                            applicable, _ = contam_source_applicable(
                                lib_num, assignment_source)
                        except Exception:
                            applicable = False
                        if not applicable:
                            continue
                    plot_prefix = get_contam_prefix(
                        lib_num, cond["abbrev"], assignment_source)
                    if assignment_source in IDENTITY_AMBIENT_ARMS:
                        if not check_output_exists(
                                lib_num, cond["abbrev"], assignment_source):
                            failures.append(
                                f"lib{lib_num}: AMBIENT_PLOTS {assignment_source} "
                                "output is missing, incomplete, or belongs to a "
                                f"stale comparison plan: {plot_prefix}")
                        continue
                    for suffix in (".contam_rate", ".contam_prof"):
                        path = plot_prefix + suffix
                        if not check_file_exists(path):
                            failures.append(
                                f"lib{lib_num}: AMBIENT_PLOTS input missing: {path}")

        if "POSTHOC" in stages:
            if not check_file_exists(get_expected_pool_metadata()):
                failures.append(f"Expected pool metadata missing: {get_expected_pool_metadata()}")
            if get_allowed_identities() and not check_file_exists(get_allowed_identities()):
                failures.append(f"Allowed identities file missing: {get_allowed_identities()}")
            prefix = get_demux_prefix(lib_num)
            posthoc_needed = [prefix + ext for ext in [
                ".counts", ".condf", ".samples", ".assignments", ".diagnostics.gz", ".runner_ups.gz",
                ".species_counts", ".species_condf", ".species_samples"]]
            if "DEMUX" not in stages:
                for f in posthoc_needed:
                    if not check_file_exists(f):
                        failures.append(f"lib{lib_num}: POSTHOC input missing: {f}")

        if "POSTHOC_SUMMARY" in stages:
            if not check_file_exists(get_expected_pool_metadata()):
                failures.append(f"Expected pool metadata missing: {get_expected_pool_metadata()}")
            manifest = os.path.join(get_audit_lib_dir(lib_num), f"lib{lib_num}.capabilities.json")
            if not check_file_exists(manifest):
                failures.append(f"lib{lib_num}: POSTHOC_SUMMARY capabilities manifest missing: {manifest}")

        if "UNEXPECTED_COMPONENT_NN" in stages:
            report = get_swap_report_path(lib_num)
            if not check_file_exists(report):
                failures.append(f"lib{lib_num}: UNEXPECTED_COMPONENT_NN swap report missing: {report}")
            elif has_unexpected_component_signal(lib_num):
                uc_needed = [
                    os.path.join(get_audit_lib_dir(lib_num), f"lib{lib_num}.swap_scores.tsv"),
                    get_call_qc_path(lib_num),
                    get_ploidy_calls_path(lib_num),
                ]
                for f in uc_needed:
                    if not check_file_exists(f):
                        failures.append(f"lib{lib_num}: UNEXPECTED_COMPONENT_NN input missing: {f}")

        if "HYBRID" in stages and "CONTAM" not in stages and "POSTHOC" not in stages:
            hybrid_needed = [
                get_swap_report_path(lib_num), get_call_qc_path(lib_num), get_species_qc_path(lib_num),
                get_contam_prefix(lib_num, args.hybrid_individual_condition) + ".contam_rate",
                get_contam_prefix(lib_num, args.hybrid_species_condition) + ".contam_rate",
                get_contam_prefix(lib_num, args.hybrid_species_condition) + ".species_prof",
            ]
            if args.hybrid_fixed_species_condition:
                hybrid_needed.extend([
                    get_contam_prefix(lib_num, args.hybrid_fixed_species_condition) + ".contam_rate",
                    get_contam_prefix(lib_num, args.hybrid_fixed_species_condition) + ".species_prof",
                ])
            for f in hybrid_needed:
                if not check_file_exists(f):
                    failures.append(f"lib{lib_num}: HYBRID input missing: {f}")

        if ({"IDENTITY_RECONCILIATION", "IDENTITY_FINAL_EVIDENCE",
             "IDENTITY_SCORE"} & set(stages)):
            identity_stage_label = (
                "IDENTITY_RECONCILIATION"
                if "IDENTITY_RECONCILIATION" in stages else
                "IDENTITY_FINAL_EVIDENCE"
                if "IDENTITY_FINAL_EVIDENCE" in stages else
                "IDENTITY_SCORE")
            prefix = get_demux_prefix(lib_num)
            for ext in (".counts", ".samples", ".assignments"):
                path = prefix + ext
                if not check_file_exists(path):
                    failures.append(
                        f"lib{lib_num}: {identity_stage_label} core demux input missing: {path}")
            for ext in (".pileup_sites.tsv.gz", ".pileup_obs.tsv.gz"):
                path = prefix + ext
                produced_by_forced_demux = (
                    "DEMUX" in stages and args.force and
                    not args.individual_only_demux)
                if not produced_by_forced_demux and not check_file_exists(path):
                    failures.append(
                        f"lib{lib_num}: {identity_stage_label} nuclear pileup missing: "
                        f"{path} (generate with forced full DEMUX)")

        if "IDENTITY_FINAL_EVIDENCE_ONLY" in stages:
            selected_inputs = {
                "candidate decisions": get_identity_candidate_path(
                    lib_num, args),
                "reconciled decisions": get_reconciled_cells_path(lib_num),
                "reconciled assignments": get_reconciled_assignments_path(
                    lib_num),
            }
            missing_selected = [
                label for label, path in selected_inputs.items()
                if not check_file_exists(path)]
            if missing_selected:
                failures.append(
                    f"lib{lib_num}: IDENTITY_FINAL_EVIDENCE_ONLY missing "
                    + ", ".join(missing_selected))

        if "IDENTITY_FINALIZE_ONLY" in stages:
            selected_inputs = {
                "demux assignments": get_demux_prefix(lib_num)
                + ".assignments",
                "candidate decisions": get_identity_candidate_path(
                    lib_num, args),
                "reconciled decisions": get_reconciled_cells_path(lib_num),
                "reconciled assignments": get_reconciled_assignments_path(
                    lib_num),
            }
            missing_selected = [
                label for label, path in selected_inputs.items()
                if not check_file_exists(path)]
            if missing_selected:
                failures.append(
                    f"lib{lib_num}: IDENTITY_FINALIZE_ONLY missing "
                    + ", ".join(missing_selected))

        if "IDENTITY_RECONCILIATION" in stages:
            prefix = get_demux_prefix(lib_num)
            diagnostics = prefix + ".diagnostics.gz"
            if not check_file_exists(diagnostics):
                failures.append(
                    f"lib{lib_num}: IDENTITY_RECONCILIATION demux diagnostics missing: "
                    f"{diagnostics}")
            nn_path = get_ploidy_calls_path(lib_num)
            if "PLOIDY_NN" not in stages and not check_file_exists(nn_path):
                failures.append(
                    f"lib{lib_num}: IDENTITY_RECONCILIATION NN input missing: {nn_path}")
            identity_audit_lib = os.path.join(
                args.identity_audit_root, f"lib{lib_num}")
            identity_audit_paths = (
                os.path.join(identity_audit_lib, f"lib{lib_num}.call_qc.tsv.gz"),
                os.path.join(identity_audit_lib, f"lib{lib_num}.swap_report.tsv"),
            )
            for apath in identity_audit_paths:
                if "POSTHOC" not in stages and not check_file_exists(apath):
                    failures.append(
                        f"lib{lib_num}: IDENTITY_RECONCILIATION audit input missing: {apath}")
            # RNA mode is deliberately ATAC-blind: no ATAC path is expanded,
            # validated, globbed, or opened here.

        if {"IDENTITY_SCORE", "IDENTITY_SCORE_AGGREGATE_ONLY"} & set(stages):
            decisions = os.path.join(
                get_identity_subdir(args, "decisions"),
                f"lib{lib_num}.reconciled_cells.tsv.gz")
            if not check_file_exists(decisions):
                failures.append(
                    f"lib{lib_num}: identity scoring requires an existing "
                    f"reconciliation decision table: {decisions}")
            if "IDENTITY_SCORE_AGGREGATE_ONLY" in stages:
                frozen_inputs = (
                    get_identity_score_pair_path(lib_num, args),
                    get_identity_score_pair_summary_path(lib_num, args),
                    get_identity_probability_path(lib_num, args),
                )
                for path in frozen_inputs:
                    if not check_file_exists(path):
                        failures.append(
                            f"lib{lib_num}: aggregate-only frozen score input missing: {path}")

        if "IDENTITY_RECONCILE_ONLY" in stages:
            prefix = get_demux_prefix(lib_num)
            required_existing = [
                prefix + ".assignments",
                get_identity_candidate_path(lib_num, args),
                get_identity_nuclear_score_path(lib_num, args),
                os.path.join(get_identity_subdir(args, "nuclear"), f"lib{lib_num}.identity_site_fold_scores.tsv.gz"),
                get_ploidy_calls_path(lib_num),
                os.path.join(args.identity_audit_root, f"lib{lib_num}", f"lib{lib_num}.call_qc.tsv.gz"),
                os.path.join(get_identity_subdir(args, "doublet_context"), f"lib{lib_num}.doublet_dragon_summary.tsv"),
            ]
            if args.identity_require_mt:
                required_existing.append(get_mt_identity_score_path(lib_num, args))
            if args.identity_evidence_mode == "rna-atac":
                required_existing.append(get_atac_identity_score_path(lib_num, args))
            for path in required_existing:
                if not check_file_exists(path):
                    failures.append(
                        f"lib{lib_num}: IDENTITY_RECONCILE_ONLY existing evidence input missing: {path}")

        if "MT_FUSION" in stages:
            if (mt_manifest_libraries is not None and
                    lib_num not in mt_manifest_libraries):
                continue
            bam = get_mt_bam_path(lib_num, args)
            if not bam or not check_file_exists(bam):
                failures.append(f"lib{lib_num}: MT_FUSION BAM missing: {bam}")
            auto_empty_barcodes = (
                args.mt_ambient_qc_max is not None and
                not args.mt_ambient_none and
                not args.mt_ambient_profile_template and
                not args.mt_empty_barcodes_template
            )
            if auto_empty_barcodes:
                filtered_barcodes = get_filtered_barcodes(lib_num)
                if not check_file_exists(filtered_barcodes):
                    failures.append(
                        f"lib{lib_num}: MT_FUSION filtered barcodes missing for automatic "
                        f"MT non-cell barcode derivation: {filtered_barcodes}")
            try:
                library_id = get_mt_library_id(lib_num, args)
            except Exception as e:
                failures.append(f"lib{lib_num}: invalid --mt-library-id-template: {e}")
                library_id = None
            if not library_id:
                failures.append(f"lib{lib_num}: MT_FUSION library_id resolved empty")
            assignments = get_mt_assignments_path(lib_num, args)
            assignment_generated_upstream = (
                (args.mt_assignment_source in ("auto", "reconciled") and
                 ("IDENTITY_RECONCILIATION" in stages or "IDENTITY_RECONCILE_ONLY" in stages)) or
                (args.mt_assignment_source == "refined" and
                 "TETRA_REFINE" in stages) or
                (args.mt_assignment_source == "demux" and
                 "DEMUX" in stages))
            if (not assignment_generated_upstream and
                    not check_file_exists(assignments)):
                failures.append(
                    f"lib{lib_num}: MT_FUSION assignments missing: {assignments}")
            for template, label in [
                (args.mt_ambient_profile_template, "ambient profile"),
                (args.mt_empty_barcodes_template, "empty barcodes"),
            ]:
                if template:
                    try:
                        path = _expand_mt_template(template, lib_num)
                    except Exception as e:
                        failures.append(f"lib{lib_num}: invalid MT_FUSION {label} template: {e}")
                        continue
                    if not check_file_exists(path):
                        failures.append(f"lib{lib_num}: MT_FUSION {label} missing: {path}")
            rna_ambient_path = get_mt_rna_ambient_path(lib_num, args)
            if rna_ambient_path:
                ambient_source = args.mt_rna_ambient_assignment_source
                ambient_planned = (
                    "CONTAM" in stages and
                    args.mt_rna_ambient_condition in {
                        condition["abbrev"] for condition in cond_list} and
                    ambient_source in resolve_contam_assignment_sources(
                        args.contam_assignment_source))
                if ambient_source in IDENTITY_AMBIENT_ARMS:
                    try:
                        applicable, reason = contam_source_applicable(
                            lib_num, ambient_source)
                        if not applicable:
                            failures.append(
                                f"lib{lib_num}: MT_FUSION RNA ambient arm "
                                f"{ambient_source} is not applicable: {reason}")
                            continue
                    except Exception as exc:
                        failures.append(
                            f"lib{lib_num}: MT_FUSION RNA ambient comparison "
                            f"plan failed: {exc}")
                        continue
                if ambient_planned:
                    # The generated MT job is dependency-bound to this exact
                    # CONTAM job below; no pre-existing rate is required.
                    continue
                try:
                    if ambient_source in IDENTITY_AMBIENT_ARMS:
                        get_mt_rna_ambient_provenance(lib_num, args)
                    elif not check_file_exists(rna_ambient_path):
                        raise FileNotFoundError(rna_ambient_path)
                except (OSError, ValueError) as exc:
                    failures.append(
                        f"lib{lib_num}: MT_FUSION RNA ambient input is missing, "
                        f"incomplete, or stale: {exc}")

    # Shared memory segments (used by CONDF and DEMUX) are node-local to
    # pika/char where the VCF daemon runs. Pre-flight runs on the head node
    # (ash) and cannot see them. The sbatch scripts validate shared memory
    # at runtime on the compute node instead.

    # 5. Shared CONDF source is only needed for DEMUX generation when the
    # selected output directory does not already contain its self-contained
    # demux_input_20M.condf copy. EMPTY_DROPS and CONTAM consume the local
    # per-library .condf links and therefore do not depend on CONDF_DIR.
    if "CONDF" not in stages and "DEMUX" in stages and not args.individual_only_demux:
        for lib_num in lib_nums:
            local_condf = get_local_condf_path(lib_num, "interindiv_20M")
            if not check_file_exists(local_condf):
                shared_condf = CONDF_PATHS["interindiv_20M"]
                if not check_file_exists(shared_condf):
                    failures.append(
                        f"lib{lib_num}: local CONDF missing and shared DEMUX source unavailable: "
                        f"{local_condf} ; {shared_condf}")

    return (len(failures) == 0, failures)


def run(args):
    """Main orchestration logic."""
    global AUDIT_ROOT, HYBRID_ROOT, MT_FUSION_ROOT, IDENTITY_RECONCILIATION_ROOT, TETRA_REFINE_ROOT, PLOIDY_CALLS_ROOT
    global EXPECTED_POOL_METADATA, ALLOWED_IDENTITIES, SHARED_VCF, DAEMON_NODELIST
    global MANAGE_VCF_DAEMONS, MANAGED_VCF_RUN_ID, MANAGED_VCF_READY_FILE
    global CONDF_DIR, CONDF_PATHS, PANEL_METADATA, DEMUX_OUTPUT_ROOT
    global IDENTITY_AMBIENT_CANDIDATE_SET
    AUDIT_ROOT = args.audit_root
    HYBRID_ROOT = args.hybrid_root
    MT_FUSION_ROOT = os.path.abspath(args.mt_output_root)
    IDENTITY_RECONCILIATION_ROOT = os.path.abspath(args.identity_reconciliation_root)
    TETRA_REFINE_ROOT = args.refined_assignments_root
    PLOIDY_CALLS_ROOT = args.ploidy_calls_root
    EXPECTED_POOL_METADATA = args.expected_pool_metadata
    if (not os.path.isfile(EXPECTED_POOL_METADATA)) and os.path.isfile(os.path.join(PROCESS_SCRIPTS_DIR, "pool_combinations.tsv")):
        EXPECTED_POOL_METADATA = os.path.join(PROCESS_SCRIPTS_DIR, "pool_combinations.tsv")
    ALLOWED_IDENTITIES = args.allowed_identities
    DEMUX_OUTPUT_ROOT = (
        os.path.abspath(args.demux_output_root)
        if args.demux_output_root else None
    )
    CONDF_DIR = os.path.abspath(args.condf_dir)
    CONDF_PATHS = {
        "interindiv_20M": os.path.join(CONDF_DIR, "demux_input_20M.condf"),
        "interindiv_het_10M": os.path.join(CONDF_DIR, "demux_input_HET_10M.condf"),
        "species_20M": os.path.join(CONDF_DIR, "demux_input_species_20M.condf"),
    }
    PANEL_METADATA = os.path.abspath(args.panel_metadata)
    IDENTITY_AMBIENT_CANDIDATE_SET = args.reconciliation_candidate_set
    SHARED_VCF["interindiv_20M"] = args.shared_main_segment
    SHARED_VCF["species_20M"] = args.shared_species_segment
    SHARED_VCF["interindiv_het_10M"] = args.shared_het_segment
    DAEMON_NODELIST = args.daemon_nodes.strip()
    MANAGE_VCF_DAEMONS = bool(args.manage_vcf_daemons)
    if args.refine_contam_condition == "":
        args.refine_contam_condition = None
    if args.ploidy_nn_module == "":
        args.ploidy_nn_module = None
    if args.identity_with_atac:
        args.identity_evidence_mode = "rna-atac"
    if (args.identity_probability_resamples < 0 or
            not 0.0 <= args.identity_poor_fit_residual <= 1.0):
        print(
            "❌ Identity scoring requires non-negative resamples and "
            "--identity-poor-fit-residual within [0,1].")
        sys.exit(1)

    if args.lib19_full_test:
        args.with_posthoc = True
        args.with_hybrid = True
        if not args.libraries:
            args.libraries = ["19"]
        if not args.stage:
            args.stage = "CONDF,DEMUX,EMPTY_DROPS,CONTAM,PLOIDY_NN,TETRA_REFINE,POSTHOC,HYBRID"

    # Resolve libraries
    if args.libraries:
        lib_nums = parse_library_range(args.libraries)
    else:
        lib_nums = list(ALL_LIBRARIES)

    # Validate library numbers
    bad_libs = [l for l in lib_nums if l not in LIB_INFO]
    if bad_libs:
        print(f"❌ Unknown library numbers: {bad_libs}")
        print(f"   Valid range: 1-40")
        sys.exit(1)

    # Resolve conditions. Explicit --conditions wins. Otherwise use the named
    # roster; the focused production condition is the default and ``all`` is
    # retained as an explicit compatibility set.
    if args.conditions:
        bad_conds = [c for c in args.conditions if c not in COND_BY_ABBREV]
        if bad_conds:
            print(f"❌ Unknown conditions: {', '.join(bad_conds)}")
            print(f"   Available: {', '.join(ALL_CONDITION_ABBREVS)}")
            sys.exit(1)
        cond_list = [COND_BY_ABBREV[c] for c in args.conditions]
        condition_selection_label = "explicit --conditions"
    else:
        cond_list = resolve_condition_set(args.condition_set)
        condition_selection_label = f"--condition-set {args.condition_set}"

    if (args.mt_rna_ambient_condition and
            args.mt_rna_ambient_condition not in COND_BY_ABBREV):
        print(
            "❌ Unknown --mt-rna-ambient-condition: "
            f"{args.mt_rna_ambient_condition}")
        sys.exit(1)
    if (args.mt_rna_ambient_condition and
            args.mt_rna_ambient_assignment_source in IDENTITY_AMBIENT_ARMS and
            COND_BY_ABBREV[args.mt_rna_ambient_condition]["mode"] not in
            (1, 3, "1+SR")):
        print(
            "❌ Four-arm RNA ambient inputs for MT_FUSION require an "
            "individual-native contamination condition.")
        sys.exit(1)
    if (args.mt_rna_ambient_max is not None and
            not 0.0 <= args.mt_rna_ambient_max <= 0.99):
        print("❌ --mt-rna-ambient-max must be between 0 and 0.99")
        sys.exit(1)
    if args.mt_rna_ambient_max is not None and not args.mt_rna_ambient_condition:
        print("❌ --mt-rna-ambient-max requires --mt-rna-ambient-condition")
        sys.exit(1)

    print(
        f"Condition selection: {condition_selection_label}; "
        f"{len(cond_list)} condition(s)"
    )
    if args.condition_set != "all" or args.conditions:
        print("  " + " ".join(c["abbrev"] for c in cond_list))
    print()

    legacy_requested = [c["abbrev"] for c in cond_list if c["mode"] == "1+SR"]
    if legacy_requested:
        print("⚠️  Including parked legacy 1+SR comparator condition(s): "
              + ", ".join(legacy_requested))
        print("   These remain legacy comparison paths, not native species or production joint-evidence modes.")

    # Resolve stages
    if args.stage:
        stages = {stage.strip() for stage in args.stage.upper().split(",") if stage.strip()}
        valid_stages = {"CLEANUP_RESULTS", "CONDF", "DEMUX", "EMPTY_DROPS", "CONTAM", "GEX_AMBIENT", "AMBIENT_PLOTS", "AMBIENT_VALIDATE", "AMBIENT_SWAP_TEST", "PLOIDY_NN", "TETRA_REFINE", "POSTHOC", "POSTHOC_SUMMARY", "UNEXPECTED_COMPONENT_NN", "HYBRID", "IDENTITY_SCORE", "IDENTITY_SCORE_AGGREGATE_ONLY", CANDIDATE_AXIS_STAGE, "IDENTITY_RECONCILIATION", "IDENTITY_FINAL_EVIDENCE", "IDENTITY_FINAL_EVIDENCE_ONLY", "IDENTITY_FINALIZE_ONLY", "IDENTITY_RECONCILE_ONLY", "MT_FUSION", "MT_POPULATION"}
        bad_stages = stages - valid_stages
        if bad_stages:
            print(f"❌ Unknown stages: {bad_stages}. Valid: {', '.join(sorted(valid_stages))}")
            sys.exit(1)
    else:
        stages = {"CONDF", "DEMUX", "CONTAM"}
        if any(c.get("needs_empty_drops") for c in cond_list):
            stages.add("EMPTY_DROPS")

    if args.skip_demux:
        stages.discard("DEMUX")

    # A dedicated CONDF refresh must not force downstream stages.  It simply
    # ensures Stage 1 is included and rebuilds all central CONDF artifacts.
    # Any DEMUX jobs submitted in the same invocation inherit the normal
    # afterok dependency on the newly submitted CONDF jobs.
    if args.regenerate_condf:
        stages.add("CONDF")

    if args.with_refine:
        stages.add("TETRA_REFINE")
    if args.with_ploidy_nn:
        stages.add("PLOIDY_NN")
        stages.add("TETRA_REFINE")
    if args.with_posthoc:
        stages.add("POSTHOC")
    if args.with_ambient_plots:
        stages.add("AMBIENT_PLOTS")
    if args.with_hybrid:
        stages.add("HYBRID")
        stages.add("POSTHOC")
    if args.with_mt_fusion:
        stages.add("MT_FUSION")
    if args.with_mt_population:
        stages.add("MT_FUSION")
        stages.add("MT_POPULATION")
    if args.with_identity_reconciliation:
        stages.add("IDENTITY_RECONCILIATION")
    if args.with_identity_score:
        stages.add("IDENTITY_SCORE")
    if args.identity_reconcile_only:
        # Standalone shortcut: preserve the historical meaning of
        # --identity-reconcile-only as reconcile + validate only.  When the
        # caller explicitly requested another stage (for example MT_FUSION),
        # compose the lightweight reconciliation refresh with that stage
        # instead of silently replacing the entire requested stage set.
        if args.stage or args.with_mt_fusion or args.with_mt_population:
            stages.add("IDENTITY_RECONCILE_ONLY")
        else:
            stages = {"IDENTITY_RECONCILE_ONLY"}
    if "IDENTITY_RECONCILE_ONLY" in stages:
        stages.discard("IDENTITY_RECONCILIATION")
    if "IDENTITY_RECONCILIATION" in stages:
        # The routine reconciliation graph owns its frozen candidate-axis and
        # ambient evidence tail. Legacy probability scoring remains explicit.
        stages.discard("IDENTITY_SCORE")
        stages.discard("IDENTITY_FINAL_EVIDENCE")
    if ("IDENTITY_FINAL_EVIDENCE" in stages and
            stages != {"IDENTITY_FINAL_EVIDENCE"}):
        print("❌ IDENTITY_FINAL_EVIDENCE is a standalone resume stage.")
        sys.exit(1)
    if ("IDENTITY_FINAL_EVIDENCE_ONLY" in stages and
            stages != {"IDENTITY_FINAL_EVIDENCE_ONLY"}):
        print("❌ IDENTITY_FINAL_EVIDENCE_ONLY is a standalone resume stage.")
        sys.exit(1)
    if ("IDENTITY_FINALIZE_ONLY" in stages and
            stages != {"IDENTITY_FINALIZE_ONLY"}):
        print("❌ IDENTITY_FINALIZE_ONLY is a standalone resume stage.")
        sys.exit(1)
    if "IDENTITY_SCORE_AGGREGATE_ONLY" in stages and len(stages) != 1:
        print("❌ IDENTITY_SCORE_AGGREGATE_ONLY is a standalone frozen-input stage.")
        print("   Finish reconciliation/scoring first, then run aggregate-only separately.")
        sys.exit(1)
    if CANDIDATE_AXIS_STAGE in stages and stages != {CANDIDATE_AXIS_STAGE}:
        print("❌ IDENTITY_CANDIDATE_AXIS is a standalone frozen-input stage.")
        print("   Run it alone after the requested reconciliations are validated.")
        sys.exit(1)
    if CANDIDATE_AXIS_STAGE in stages:
        targeted_event = bool(args.identity_candidate_axis_event_id)
        targeted_proposal = bool(args.identity_candidate_axis_proposal)
        if targeted_event != targeted_proposal:
            print(
                "❌ Candidate-axis event ID and proposal must be supplied together.")
            sys.exit(1)
        if targeted_event and len(lib_nums) != 1:
            print("❌ Targeted candidate-axis mode requires exactly one library.")
            sys.exit(1)
        if not args.identity_candidate_axis_output_root:
            print("❌ IDENTITY_CANDIDATE_AXIS requires an explicit output root.")
            sys.exit(1)
        if not os.path.isabs(args.identity_candidate_axis_output_root):
            print("❌ --identity-candidate-axis-output-root must be absolute.")
            sys.exit(1)
        if (args.identity_candidate_axis_v6_3_root and
                not os.path.isabs(args.identity_candidate_axis_v6_3_root)):
            print("❌ --identity-candidate-axis-v6-3-root must be absolute.")
            sys.exit(1)
    if {"IDENTITY_SCORE", "IDENTITY_RECONCILE_ONLY"} <= stages:
        print("❌ IDENTITY_SCORE cannot run with IDENTITY_RECONCILE_ONLY in one command.")
        print("   Refresh and validate reconciliation first, then launch IDENTITY_SCORE.")
        sys.exit(1)
    if "MT_FUSION" in stages and not args.mt_ratios_only:
        stages.add("MT_POPULATION")

    if "CLEANUP_RESULTS" in stages and stages != {"CLEANUP_RESULTS"}:
        print("❌ CLEANUP_RESULTS is a standalone maintenance stage.")
        print("   Run it by itself before submitting the next pipeline run.")
        sys.exit(1)

    if "AMBIENT_VALIDATE" in stages:
        if (args.contam_assignment_source != IDENTITY_AMBIENT_SELECTOR or
                args.reconciliation_candidate_set != "applied"):
            print(
                "AMBIENT_VALIDATE automatically uses applied reconciliation "
                "assignments and rosters.")
        args.contam_assignment_source = IDENTITY_AMBIENT_SELECTOR
        args.reconciliation_candidate_set = "applied"
        IDENTITY_AMBIENT_CANDIDATE_SET = "applied"

    if "AMBIENT_SWAP_TEST" in stages and len(cond_list) > 8:
        print("❌ AMBIENT_SWAP_TEST accepts at most eight conditions per run.")
        print("   Use --conditions or a smaller --condition-set.")
        sys.exit(1)
    if "AMBIENT_PLOTS" in stages and len(cond_list) > 8:
        print("❌ AMBIENT_PLOTS accepts at most eight selected conditions per run.")
        print("   Use --conditions or a smaller --condition-set; no ranking requires all conditions together.")
        sys.exit(1)
    if {"AMBIENT_VALIDATE", "AMBIENT_SWAP_TEST"} & stages:
        unsupported = [
            condition["abbrev"] for condition in cond_list
            if condition["mode"] not in (1, 3)]
        if unsupported:
            selected = ", ".join(sorted(
                {"AMBIENT_VALIDATE", "AMBIENT_SWAP_TEST"} & stages))
            print(f"❌ {selected} supports individual modes 1/3 only: "
                  + ", ".join(unsupported))
            sys.exit(1)
    if (args.contam_assignment_source == IDENTITY_AMBIENT_SELECTOR and
            ({"CONTAM", "AMBIENT_PLOTS", "AMBIENT_VALIDATE", "GEX_AMBIENT"}
             & stages) and
            ({"IDENTITY_RECONCILIATION", "IDENTITY_RECONCILE_ONLY"} & stages)):
        print("❌ reconciliation-four-arm consumes a completed, validated identity reconciliation.")
        print("   Run identity reconciliation first, then launch CONTAM/AMBIENT_PLOTS in a second command.")
        sys.exit(1)
    if ("AMBIENT_SWAP_TEST" in stages and
            ({"IDENTITY_RECONCILIATION", "IDENTITY_RECONCILE_ONLY"} & stages)):
        print("❌ AMBIENT_SWAP_TEST consumes a completed identity reconciliation ledger.")
        print("   Finish identity reconciliation first, then launch the swap test.")
        sys.exit(1)
    if (args.mt_rna_ambient_condition and
            args.mt_rna_ambient_assignment_source in IDENTITY_AMBIENT_ARMS and
            "MT_FUSION" in stages and
            ({"IDENTITY_RECONCILIATION", "IDENTITY_RECONCILE_ONLY"} & stages)):
        print("❌ Four-arm RNA ambient MT_FUSION consumes a completed identity reconciliation plan.")
        print("   Finish identity reconciliation first, then launch CONTAM/MT_FUSION in a second command.")
        sys.exit(1)
    if ("AMBIENT_PLOTS" in stages and args.ambient_reference_condition and
            args.ambient_reference_condition not in {
                condition["abbrev"] for condition in cond_list}):
        print("❌ --ambient-reference-condition must be one of the selected conditions.")
        sys.exit(1)
    if "GEX_AMBIENT" in stages:
        if (args.gex_cluster_source == "manual" and
                not args.gex_clusters_template):
            print("❌ --gex-cluster-source manual requires --gex-clusters-template.")
            sys.exit(1)
        unsupported = [
            condition["abbrev"] for condition in cond_list
            if condition["mode"] != 1
        ]
        if unsupported:
            print("❌ GEX_AMBIENT currently accepts individual-native mode-1 conditions only.")
            print("   Unsupported selected conditions: " + ", ".join(unsupported))
            print("   For the planned pilot use --condition-set ck-minimal.")
            sys.exit(1)
        if (args.gex_cpus < 1 or args.gex_memory_gb < 1 or args.gex_top_n < 1 or
                args.gex_auto_data_features < 0 or args.gex_auto_min_genes < 1 or
                args.gex_auto_neighbors < 2 or args.gex_auto_pca_components < 2 or
                args.gex_auto_stability_repeats < 2):
            print("❌ GEX resources/top-N must be positive; auto feature count may be 0, and min genes/neighbors/PCA/repeats must be valid")
            sys.exit(1)
        try:
            re.compile(args.gex_skip_genes_regex)
            re.compile(args.gex_auto_exclude_genes_regex)
            resolutions = [
                float(value) for value in args.gex_auto_resolution_grid.split(",")
                if value.strip()]
            if not resolutions or any(value <= 0 for value in resolutions):
                raise ValueError("--gex-auto-resolution-grid requires positive values")
            _gex_safe_name(args.gex_analysis_name)
        except (re.error, ValueError) as exc:
            print(f"❌ Invalid GEX_AMBIENT analysis option: {exc}")
            sys.exit(1)

    # Reconciliation-owned files are products of IDENTITY_RECONCILIATION, not
    # prerequisites for a first MT run. Add that composite stage automatically
    # only when the selected MT path actually needs those products and they are
    # not already present. Generic --force applies to the explicitly selected
    # downstream stage; identity regeneration must be requested explicitly.
    mt_needs_identity = False
    if "MT_POPULATION" in stages:
        mt_needs_identity = any(
            not check_file_exists(get_reconciled_cells_path(lib_num))
            for lib_num in lib_nums)
    if ("MT_FUSION" in stages and
            args.mt_assignment_source in ("auto", "reconciled")):
        mt_needs_identity = mt_needs_identity or any(
            not check_file_exists(get_reconciled_assignments_path(lib_num))
            for lib_num in lib_nums)
    if mt_needs_identity and "IDENTITY_RECONCILE_ONLY" not in stages:
        if (args.mt_rna_ambient_condition and
                args.mt_rna_ambient_assignment_source in IDENTITY_AMBIENT_ARMS):
            print("❌ Four-arm RNA ambient MT_FUSION needs completed identity reconciliation outputs.")
            print("   Run identity reconciliation first, then launch the controlled ambient comparison.")
            sys.exit(1)
        stages.add("IDENTITY_RECONCILIATION")

    daemon_stages_selected = bool({"CONDF", "DEMUX"} & stages)
    if MANAGE_VCF_DAEMONS and daemon_stages_selected:
        try:
            daemon_nodes = get_vcf_daemon_nodes()
        except ValueError as exc:
            print(f"❌ {exc}")
            sys.exit(1)
        if not daemon_nodes:
            print("❌ Managed VCF daemons require at least one --daemon-nodes host.")
            sys.exit(1)
        configure_managed_vcf_segments()
    else:
        MANAGED_VCF_RUN_ID = None
        MANAGED_VCF_READY_FILE = None

    mt_analysis_libs = list(lib_nums)
    mt_nonapplicable_libs = []
    if ({"MT_FUSION", "MT_POPULATION"} & stages and
            args.mt_site_manifest and check_file_exists(args.mt_site_manifest)):
        try:
            mt_manifest_libraries = get_mt_manifest_library_numbers(args.mt_site_manifest)
        except Exception:
            # Preflight below reports the concrete manifest error. Keep the full
            # library set here so no unrelated exception masks that diagnosis.
            mt_manifest_libraries = None
        if mt_manifest_libraries is not None:
            mt_analysis_libs = [lib for lib in lib_nums if lib in mt_manifest_libraries]
            mt_nonapplicable_libs = [lib for lib in lib_nums if lib not in mt_manifest_libraries]

    # Diagnose-only mode
    if args.diagnose_only:
        diagnose(lib_nums, cond_list, condition_selection_label)
        return

    # Pre-flight validation (runs by default, skip with --skip-validation)
    if not args.skip_validation:
        ok, failures = preflight_validate(lib_nums, cond_list, stages, args=args)
        if not ok:
            print("=" * 72)
            print("PRE-FLIGHT VALIDATION FAILED")
            print(f"  {len(failures)} problem(s) found:")
            print()
            for f in failures:
                print(f"  ❌ {f}")
            print()
            print("  Fix the issues above, or use --skip-validation to bypass")
            print("  (e.g., when upstream stages will create missing files)")
            print("=" * 72)
            sys.exit(1)
        else:
            print("✅ Pre-flight validation passed")
            print()

    # Create authoritative per-library directories and the small helper-facing
    # symlink index only after validation/diagnosis has finished.
    if stages != {CANDIDATE_AXIS_STAGE}:
        prepare_output_layout(lib_nums, cond_list, stages, args)

    # Arm A is numerically independent of candidate-set selection. Reuse a
    # validated pre-applied calculation through lightweight per-library links,
    # while writing a fresh current-plan contract for plotting/provenance.
    reused_arm_a = []
    if (args.contam_assignment_source == IDENTITY_AMBIENT_SELECTOR and
            args.reconciliation_candidate_set == "applied" and
            ({"CONTAM", "AMBIENT_PLOTS"} & stages) and not args.force):
        for lib_num in lib_nums:
            for condition in cond_list:
                source_prefix = materialize_reusable_identity_ambient_arm_a(
                    lib_num, condition["abbrev"])
                if source_prefix:
                    reused_arm_a.append((
                        lib_num, condition["abbrev"], source_prefix))
        if reused_arm_a:
            print("Reused candidate-independent Arm A output(s):")
            for lib_num, condition, source_prefix in reused_arm_a:
                print(
                    f"  lib{lib_num} {condition}: linked from {source_prefix}")
            print()

    # Print plan header
    print("=" * 72)
    print("PIPELINE ORCHESTRATOR")
    print(f"  Release: {ORCHESTRATOR_RELEASE}")
    print(f"  Libraries: {len(lib_nums)} ({min(lib_nums)}-{max(lib_nums)})")
    print(f"  Conditions: {len(cond_list)}")
    if {"CONTAM", "AMBIENT_PLOTS", "GEX_AMBIENT"} & stages:
        print(
            "  Contamination assignment source(s): " + ", ".join(
                resolve_contam_assignment_sources(
                    args.contam_assignment_source)))
        if args.contam_assignment_source == IDENTITY_AMBIENT_SELECTOR:
            print(
                "  Reconciliation candidate set: "
                + args.reconciliation_candidate_set)
    elif "AMBIENT_VALIDATE" in stages:
        print("  Validation arms: E demux/fixed-profile, F reconciled/fixed-profile")
        print("  Reconciliation candidate set: applied (automatic)")
    elif "AMBIENT_SWAP_TEST" in stages:
        print(
            "  Swap-test arms: G original/fixed-profile, "
            "H proposed/fixed-profile, J original/joint-profile, "
            "K proposed/joint-profile")
        print(f"  Candidate discovery: {get_identity_event_path(args)}")
    stage_order = ["CONDF", "DEMUX", "EMPTY_DROPS", "CONTAM", "GEX_AMBIENT", "AMBIENT_PLOTS", "AMBIENT_VALIDATE", "AMBIENT_SWAP_TEST", "PLOIDY_NN", "TETRA_REFINE", "POSTHOC", "POSTHOC_SUMMARY", "UNEXPECTED_COMPONENT_NN", "HYBRID", "IDENTITY_SCORE", "IDENTITY_SCORE_AGGREGATE_ONLY", CANDIDATE_AXIS_STAGE, "IDENTITY_RECONCILIATION", "IDENTITY_FINAL_EVIDENCE", "IDENTITY_FINAL_EVIDENCE_ONLY", "IDENTITY_FINALIZE_ONLY", "IDENTITY_RECONCILE_ONLY", "MT_FUSION", "MT_POPULATION"]
    stages_display = [s for s in stage_order if s in stages]
    print(f"  Stages: {', '.join(stages_display)}")
    print(f"  Submit: {'YES' if args.submit else 'DRY RUN'}")
    print(f"  Force all selected stages: {'YES' if args.force else 'no'}")
    print(f"  Regenerate CONDF: {'YES' if args.regenerate_condf else 'no'}")
    print(f"  Individual-only DEMUX: {'enabled' if args.individual_only_demux else 'disabled'}")
    print(
        "  HET VCF panel in DEMUX: "
        + ("disabled" if (args.skip_het_vcf or args.individual_only_demux) else "enabled")
    )
    print(
        "  Species/raw/pileup DEMUX extras: "
        + ("disabled" if args.individual_only_demux else "normal production behavior")
    )
    print(
        "  DEMUX output root: "
        + (DEMUX_OUTPUT_ROOT if DEMUX_OUTPUT_ROOT else
           f"in-place under {BEEGFS_ROOT}/<library>/{DEMUX_SUBDIR}")
    )
    print(f"  CONDF directory: {CONDF_DIR}")
    print(f"  Panel metadata: {PANEL_METADATA}")
    print(f"  NoMito individual VCF source: {VCF_SOURCE_PATHS['interindiv_20M']}")
    print(f"  NoMito HET VCF source: {VCF_SOURCE_PATHS['interindiv_het_10M']}")
    print(f"  NoMito species VCF source: {VCF_SOURCE_PATHS['species_20M']}")
    if daemon_stages_selected:
        if MANAGE_VCF_DAEMONS:
            print(f"  VCF daemon lifecycle: orchestrator-managed run {MANAGED_VCF_RUN_ID}")
            print(f"  Managed main segment: {SHARED_VCF['interindiv_20M']}")
            print(f"  Managed HET segment: {SHARED_VCF['interindiv_het_10M']}")
            print(f"  Managed species segment: {SHARED_VCF['species_20M']}")
        else:
            print("  VCF daemon lifecycle: externally managed")
    if ({"IDENTITY_SCORE", "IDENTITY_SCORE_AGGREGATE_ONLY", "IDENTITY_RECONCILIATION",
         "IDENTITY_FINAL_EVIDENCE", "IDENTITY_FINAL_EVIDENCE_ONLY",
         "IDENTITY_FINALIZE_ONLY",
         "IDENTITY_RECONCILE_ONLY"} & stages):
        print(f"  IDENTITY_RECONCILIATION root: {args.identity_reconciliation_root}")
        print(f"  IDENTITY evidence mode: {args.identity_evidence_mode}")
        print(f"  IDENTITY audit root: {args.identity_audit_root}")
        if {"IDENTITY_SCORE", "IDENTITY_SCORE_AGGREGATE_ONLY"} & stages:
            print(f"  Legacy score summary: {get_identity_score_output_root(args)}")
            print(
                "  Legacy molecule-aware probability: preferred when "
                "pileup_molecules exists; explicit site fallback otherwise")
            print(
                "  Legacy ambient annotation root: "
                f"{args.identity_score_ambient_root or 'disabled'}")
        if "IDENTITY_RECONCILIATION" in stages:
            print(
                "  Final evidence: frozen candidate axis plus controlled and "
                "four-arm ambient comparisons")
        if "IDENTITY_FINAL_EVIDENCE" in stages:
            print(
                "  Final-evidence resume: reuse completed reconciliation and "
                "submit only missing evidence plus canonical finalization")
        if "IDENTITY_FINAL_EVIDENCE_ONLY" in stages:
            print(
                "  Final-evidence-only resume: reuse all completed per-library "
                "evidence; submit selected-library aggregates and finalization")
        if "IDENTITY_FINALIZE_ONLY" in stages:
            print(
                "  Finalize-only resume: reuse completed evidence aggregates; "
                "submit only canonical finalization")
        if "IDENTITY_SCORE" in stages:
            print("  IDENTITY_SCORE: existing reconciled proposals -> targeted original/swap probabilities only")
        if "IDENTITY_SCORE_AGGREGATE_ONLY" in stages:
            print("  IDENTITY_SCORE_AGGREGATE_ONLY: reuse frozen pair/probability files; no scoring jobs")
        if "IDENTITY_RECONCILE_ONLY" in stages:
            print("  IDENTITY reconcile-only: reuse existing evidence; submit reconcile + validate only")
    if CANDIDATE_AXIS_STAGE in stages:
        candidate_paths = _candidate_axis_paths(args)
        print(f"  IDENTITY_CANDIDATE_AXIS root: {candidate_paths['root']}")
        if args.identity_candidate_axis_event_id:
            print(
                "  Fixed event: " + args.identity_candidate_axis_event_id +
                " -> " + args.identity_candidate_axis_proposal)
        else:
            print(
                "  Finalized reconciliation events: enumerated automatically")
        print(
            "  Candidate-axis interpretation: geometric fixed-pair position; "
            "not confidence or probability an identity is correct")
        print(
            "  Jobs: PAIRS -> SCORE -> AGGREGATE; aggregate depends on both upstream jobs")
    if "MT_FUSION" in stages:
        print(f"  MT_FUSION VCF: {args.mt_vcf}")
        print(f"  MT_FUSION site manifest: {args.mt_site_manifest}")
        print(
            "  MT_FUSION RNA ambient covariate: " +
            (f"{args.mt_rna_ambient_condition} "
             f"[{args.mt_rna_ambient_assignment_source}]"
             if args.mt_rna_ambient_condition else "disabled"))
        if mt_nonapplicable_libs:
            print(
                "  MT_FUSION not applicable (zero ratio-manifest rows): " +
                ", ".join(f"lib{x}" for x in mt_nonapplicable_libs))
    if "MT_POPULATION" in stages:
        print(f"  MT_POPULATION reconciliation root: {get_identity_reconciliation_root()}")
    if "AMBIENT_PLOTS" in stages:
        print(f"  AMBIENT_PLOTS output: {get_ambient_plot_run_dir(args, cond_list)}")
    if "GEX_AMBIENT" in stages:
        print(f"  GEX_AMBIENT analysis: {_gex_safe_name(args.gex_analysis_name)}")
        print(f"  GEX_AMBIENT summary: {get_gex_ambient_summary_dir(args)}")
        if args.gex_cluster_source == "auto":
            print(
                "  GEX clusters: automatic RNA-only stable Leiden partition "
                f"from all filtered-MEX cells ({args.gex_auto_resolution_grid})")
        elif args.gex_cluster_source == "h5ad":
            print(
                f"  GEX clusters: H5AD obs[{args.gex_h5ad_cluster_column!r}] "
                f"from {args.ploidy_nn_h5ad}")
        else:
            print(f"  GEX clusters: manual template {args.gex_clusters_template}")
        if args.gex_marker_genes:
            print(f"  GEX auto marker union: {args.gex_marker_genes}")
        print("  GEX gene scope: " +
              (f"exclude /{args.gex_skip_genes_regex}/"
               if args.gex_skip_genes_regex else
               "all parsed Gene Expression features"))
    if "AMBIENT_VALIDATE" in stages:
        print(
            "  AMBIENT_VALIDATE output: "
            f"{get_ambient_plot_run_dir(args, cond_list)}/fixed_profile_validation")
        print(
            "  Fixed empty-drop profile: one deterministic sample of "
            f"up to {AMBIENT_VALIDATION_PROFILE_MAX_EMPTY} droplets per "
            f"applicable library (seed {AMBIENT_VALIDATION_PROFILE_SEED})")
    if "AMBIENT_SWAP_TEST" in stages:
        swap_discovery = discover_ambient_swap_events(args)
        print(
            "  AMBIENT_SWAP_TEST output: "
            + get_ambient_swap_test_output_root(
                args, cond_list,
                ambient_swap_discovery_id(swap_discovery)))
        print(
            "  Candidate selection: STRONG/DECISIVE supported unexpected "
            "exact-identity events; diploid-to-heterotypic cells require "
            f"QC-passing NN P(tet) >= "
            f"{AMBIENT_SWAP_HETEROTYPIC_NN_MIN_PROB:.2f}")
        print(
            "  Fixed empty-drop profile: one deterministic sample of "
            f"up to {AMBIENT_VALIDATION_PROFILE_MAX_EMPTY} droplets per "
            f"applicable library (seed {AMBIENT_VALIDATION_PROFILE_SEED})")
    if "CLEANUP_RESULTS" in stages:
        print(f"  Cleanup root: {AGGREGATE_ROOT}")
        print(
            "  Retention: newest completed result per generated workflow; "
            "newest logical log/script attempt")
        print(
            "  Protected: CONDF, indexes, reconciliation, posthoc, hybrid, "
            "ploidy, tetra-refine, mitochondrial, GEX ambient, migration journals")
    print(f"  Daemon-bound nodes: {DAEMON_NODELIST if DAEMON_NODELIST else 'UNPINNED BY EXPLICIT OVERRIDE'}")
    print(f"  Date: {datetime.now().strftime('%Y-%m-%d %H:%M')}")
    print("=" * 72)
    print()

    # Track generated scripts and job IDs for dependency chaining
    generated = []
    condf_job_ids = []
    demux_job_ids = {}   # lib_num -> job_id
    empty_job_ids = {}   # lib_num -> job_id
    contam_job_ids = {}  # (condition_abbrev, assignment_source, lib_num) -> job_id
    contam_generated_keys = set()  # populated for submitted and dry-run jobs
    gex_cluster_job_ids = {}  # lib_num -> automatic/H5AD cluster job id
    gex_ambient_job_ids = {}  # (condition, assignment source, lib) -> job id
    ambient_validation_profile_job_ids = {}
    ambient_validation_fixed_job_ids = []
    ambient_swap_profile_job_ids = {}
    ambient_swap_arm_job_ids = []
    ploidy_nn_job_id = None
    tetra_refine_job_ids = {}  # lib_num -> job_id
    posthoc_job_ids = {} # lib_num -> job_id
    posthoc_summary_job_ids = {}  # lib_num -> job_id
    unexpected_component_nn_job_ids = {}  # lib_num -> job_id
    hybrid_job_ids = {}  # lib_num -> job_id
    mt_fusion_job_ids = {}  # lib_num -> job_id
    mt_fusion_changed = False
    identity_metadata_job_id = None
    identity_candidates_job_ids = {}  # lib_num -> job_id
    identity_candidates_aggregate_job_id = None
    identity_doublet_context_job_id = None
    identity_score_job_ids = {}
    identity_score_pair_manifest_job_id = None
    identity_probability_job_ids = {}
    identity_probability_aggregate_job_id = None
    identity_reconcile_job_id = None
    identity_validate_job_id = None
    identity_final_evidence_job_id = None
    vcf_holder_job_ids = []
    vcf_consumer_job_ids = []
    vcf_abort_guard = None
    vcf_cleanup_armed = {"value": False}

    if CANDIDATE_AXIS_STAGE in stages:
        print("--- IDENTITY_CANDIDATE_AXIS: finalized-event fixed-pair scoring ---")
        pair_job_ids = {}
        score_job_ids = {}
        for lib_num in lib_nums:
            pair_script = generate_identity_candidate_axis_pairs_script(
                lib_num, args
            )
            generated.append((
                f"IDENTITY_CANDIDATE_AXIS_PAIRS_lib{lib_num}",
                f"lib{lib_num}", pair_script,
            ))
            pair_job_id = None
            if args.submit:
                pair_job_id = submit_job(pair_script)
                pair_job_ids[lib_num] = pair_job_id
                print(
                    f"  ✅ IDENTITY_CANDIDATE_AXIS_PAIRS_lib{lib_num}: "
                    f"job {pair_job_id}"
                )
            else:
                print(f"  Generated PAIRS lib{lib_num}: {pair_script}")

            score_script = generate_identity_candidate_axis_score_script(
                lib_num, args, [pair_job_id] if pair_job_id else None
            )
            generated.append((
                f"IDENTITY_CANDIDATE_AXIS_SCORE_lib{lib_num}",
                f"lib{lib_num}", score_script,
            ))
            score_job_id = None
            if args.submit:
                score_job_id = submit_job(score_script)
                score_job_ids[lib_num] = score_job_id
                print(
                    f"  ✅ IDENTITY_CANDIDATE_AXIS_SCORE_lib{lib_num}: "
                    f"job {score_job_id} afterok:{pair_job_id}"
                )
            else:
                print(
                    f"  Generated SCORE lib{lib_num}: {score_script} "
                    f"[afterok: IDENTITY_CANDIDATE_AXIS_PAIRS_lib{lib_num}]"
                )

        aggregate_dependencies = [
            pair_job_ids[lib_num] for lib_num in lib_nums
            if lib_num in pair_job_ids
        ] + [
            score_job_ids[lib_num] for lib_num in lib_nums
            if lib_num in score_job_ids
        ]
        aggregate_script = generate_identity_candidate_axis_aggregate_script(
            lib_nums, args, aggregate_dependencies or None)
        generated.append((
            "IDENTITY_CANDIDATE_AXIS_AGGREGATE", "all_requested", aggregate_script))
        if args.submit:
            aggregate_job_id = submit_job(aggregate_script)
            print(
                "  ✅ IDENTITY_CANDIDATE_AXIS_AGGREGATE: job "
                f"{aggregate_job_id} after all pair and score jobs")
        else:
            print(
                f"  Generated AGGREGATE: {aggregate_script} "
                "[afterok: all requested pair and score jobs]")
        print()
        print(f"Generated exactly {2 * len(lib_nums) + 1} candidate-axis jobs.")
        print(
            "No assignment, decision, resample, molecule, ambient, MT, ATAC, "
            "plot, or upstream demux job is part of this stage.")
        return

    if "CLEANUP_RESULTS" in stages:
        print("--- CLEANUP_RESULTS: remove obsolete generated results ---")
        cleanup_bundle = generate_cleanup_results_script(lib_nums)
        generated.append((
            "CLEANUP_RESULTS", "all_stage_results", cleanup_bundle["script"]))
        print(f"  Generated: {cleanup_bundle['script']}")
        print(
            "  Ledger: "
            f"{os.path.join(AGGREGATE_ROOT, 'cleanup_results_summary.tsv')}")
        if args.submit:
            cleanup_job_id = submit_job(cleanup_bundle["script"])
            print(f"  ✅ all-stage result cleanup: job {cleanup_job_id}")
        print()

    condf_force = args.force or args.regenerate_condf
    condf_work_planned = (
        "CONDF" in stages and any(
            condf_force or not check_file_exists(path)
            for path in CONDF_PATHS.values()))
    demux_work_planned = (
        "DEMUX" in stages and any(
            (args.force or not demux_outputs_complete(
                lib_num, individual_only=args.individual_only_demux)) and
            check_file_exists(get_bam_path(lib_num))
            for lib_num in lib_nums))
    vcf_service_needed = bool(
        MANAGE_VCF_DAEMONS and (condf_work_planned or demux_work_planned))

    if vcf_service_needed:
        reference_bam = resolve_vcf_daemon_reference_bam(lib_nums)
        if reference_bam is None:
            raise RuntimeError(
                "managed VCF daemon service has work but no selected reference BAM")
        print("--- Managed VCF daemons: run-scoped holder launch ---")
        state_dir = get_vcf_daemon_state_dir(MANAGED_VCF_RUN_ID)
        os.makedirs(state_dir, exist_ok=True)
        sentinel = os.path.join(state_dir, "TEARDOWN")
        if os.path.exists(sentinel):
            os.unlink(sentinel)
        for node in get_vcf_daemon_nodes():
            holder_script = generate_vcf_daemon_holder_script(
                node, reference_bam, MANAGED_VCF_RUN_ID)
            generated.append(("VCF_DAEMON_HOLDER", node, holder_script))
            print(f"  Generated holder: {node}")
            if args.submit:
                holder_id = submit_job(holder_script)
                vcf_holder_job_ids.append(holder_id)
                print(f"  ✅ {node}: holder job {holder_id}")
        if args.submit and vcf_holder_job_ids:
            def _signal_vcf_holders_on_submit_abort():
                if vcf_cleanup_armed["value"]:
                    return
                try:
                    os.makedirs(state_dir, exist_ok=True)
                    Path(sentinel).touch()
                    print(
                        "WARNING: orchestration ended before the VCF cleanup "
                        "watcher was armed; teardown sentinel planted for "
                        "submitted holder jobs", file=sys.stderr)
                except OSError as exc:
                    print(
                        f"ERROR: could not plant emergency VCF teardown "
                        f"sentinel {sentinel}: {exc}", file=sys.stderr)
            vcf_abort_guard = _signal_vcf_holders_on_submit_abort
            atexit.register(vcf_abort_guard)
        print()

    # ---- Stage 1 CONDF: .condf generation ----
    if "CONDF" in stages:
        print("--- Stage 1 CONDF: .condf generation ---")
        for vcf_key in CONDF_PATHS:
            if check_file_exists(CONDF_PATHS[vcf_key]) and not condf_force:
                print(f"  ✅ {vcf_key}: already exists, skipping")
                continue

            script_path, _ = generate_condf_script(vcf_key, force=condf_force)
            generated.append(("CONDF", vcf_key, script_path))
            print(f"  Generated: {script_path}")

            if args.submit:
                jid = submit_job(script_path)
                if jid:
                    condf_job_ids.append(jid)
                    vcf_consumer_job_ids.append(jid)
                    print(f"  ✅ Submitted: job {jid}")
        print()

    # ---- Stage 2 DEMUX: Demux ----
    if "DEMUX" in stages:
        print("--- Stage 2 DEMUX: Demux ---")
        for lib_num in lib_nums:
            # Match the generated job's complete bundle contract so a partial
            # run cannot be skipped on the login node.
            demux_done = demux_outputs_complete(
                lib_num, individual_only=args.individual_only_demux)
            if demux_done and not args.force:
                print(f"  ✅ lib{lib_num}: demux outputs exist, skipping")
                continue

            # Check BAM exists
            bam = get_bam_path(lib_num)
            if not check_file_exists(bam):
                print(f"  ⚠️  lib{lib_num}: BAM missing ({bam}), skipping")
                continue

            dep_ids = condf_job_ids if condf_job_ids else None
            script_path = generate_demux_script(
                lib_num, condf_job_ids=dep_ids, force=args.force,
                use_het_vcf=not args.skip_het_vcf,
                individual_only=args.individual_only_demux)
            generated.append(("DEMUX", f"lib{lib_num}", script_path))
            print(f"  Generated: lib{lib_num}")

            if args.submit:
                jid = submit_job(script_path)
                if jid:
                    demux_job_ids[lib_num] = jid
                    vcf_consumer_job_ids.append(jid)
                    print(f"  ✅ lib{lib_num}: job {jid}")
        print()

    if vcf_service_needed:
        cleanup_script = generate_vcf_daemon_cleanup_script(
            MANAGED_VCF_RUN_ID, vcf_holder_job_ids,
            vcf_consumer_job_ids)
        generated.append(("VCF_DAEMON_CLEANUP", MANAGED_VCF_RUN_ID,
                          cleanup_script))
        print("--- Managed VCF daemons: cleanup watcher ---")
        print(f"  Generated: {cleanup_script}")
        if args.submit:
            cleanup_id = submit_job(cleanup_script)
            print(f"  ✅ Cleanup watcher: job {cleanup_id}")
            vcf_cleanup_armed["value"] = True
            if vcf_abort_guard is not None:
                atexit.unregister(vcf_abort_guard)
        print()

    # ---- Stage 3 EMPTY_DROPS: Empty drops ----
    if "EMPTY_DROPS" in stages:
        print("--- Stage 3 EMPTY_DROPS: Empty drops ---")
        for lib_num in lib_nums:
            # Check if already done
            if (check_file_exists(get_empty_drops_indiv(lib_num)) and
                check_file_exists(get_empty_drops_species(lib_num)) and
                not args.force):
                print(f"  ✅ lib{lib_num}: empty drops outputs exist, skipping")
                continue

            # Determine dependency
            dep_id = demux_job_ids.get(lib_num)

            # If not submitting DEMUX in this run, verify demux outputs exist
            if dep_id is None and "DEMUX" not in stages:
                prefix = get_demux_prefix(lib_num)
                if not check_file_exists(prefix + ".counts"):
                    print(f"  ⚠️  lib{lib_num}: demux outputs missing, skipping")
                    continue

            # EMPTY_DROPS now runs both individual and species expected-lines
            # restrictions. Keep this guard here too so --skip-validation does
            # not produce a generated sbatch with a missing --ids input.
            try:
                el = get_expected_lines(lib_num)
                species_el = get_species_expected_lines(lib_num) if el is not None else None
            except Exception as e:
                print(f"  ⚠️  lib{lib_num}: expected_lines generation failed, skipping ({e})")
                continue
            if el is None or species_el is None:
                print(f"  ⚠️  lib{lib_num}: expected_lines missing, skipping EMPTY_DROPS")
                continue

            script_path = generate_empty_drops_script(
                lib_num, demux_job_id=dep_id, force=args.force)
            generated.append(("EMPTY_DROPS", f"lib{lib_num}", script_path))
            print(f"  Generated: lib{lib_num}")

            if args.submit:
                jid = submit_job(script_path)
                if jid:
                    empty_job_ids[lib_num] = jid
                    print(f"  ✅ lib{lib_num}: job {jid}")
        print()

    # ---- Stage 4 CONTAM: Contamination estimation ----
    if "CONTAM" in stages:
        print("--- Stage 4 CONTAM: Contamination estimation ---")
        q_count = 0
        q_skipped = 0
        q_blocked = 0
        q_not_applicable = 0

        contamination_sources = resolve_contam_assignment_sources(
            args.contam_assignment_source)
        for cond in cond_list:
            abbrev = cond["abbrev"]

            for assignment_source in contamination_sources:
                cond_generated = 0
                for lib_num in lib_nums:
                    try:
                        applicable, reason = contam_source_applicable(
                            lib_num, assignment_source)
                    except Exception as exc:
                        q_blocked += 1
                        print(
                            f"  ⚠️  {abbrev}/{assignment_source}/lib{lib_num}: "
                            f"comparison plan unavailable ({exc})")
                        continue
                    if not applicable:
                        q_not_applicable += 1
                        print(
                            f"  ↪ {abbrev}/{assignment_source}/lib{lib_num}: "
                            f"not applicable ({reason})")
                        continue

                    # Check if output exists
                    if check_output_exists(
                            lib_num, abbrev, assignment_source) and not args.force:
                        q_skipped += 1
                        continue

                    # Determine which dependency to use
                    if cond["needs_empty_drops"]:
                        dep_id = empty_job_ids.get(lib_num)
                        if dep_id is None and "EMPTY_DROPS" not in stages:
                            if cond["needs_empty_drops"] == "individual":
                                if not check_file_exists(get_empty_drops_indiv(lib_num)):
                                    q_blocked += 1
                                    continue
                            elif cond["needs_empty_drops"] in ("species", "species_fixed"):
                                if not check_file_exists(get_empty_drops_species(lib_num)):
                                    q_blocked += 1
                                    continue
                    else:
                        dep_id = demux_job_ids.get(lib_num)
                        if dep_id is None and "DEMUX" not in stages:
                            prefix = get_demux_prefix(lib_num)
                            if not check_file_exists(prefix + ".counts"):
                                q_blocked += 1
                                continue

                    # Check expected_lines at the mode's native resolution and
                    # materialize the explicit individual donor roster when needed.
                    try:
                        arm_inputs = (
                            identity_ambient_arm_inputs(lib_num, assignment_source)
                            if assignment_source in IDENTITY_AMBIENT_ARMS else None)
                        el = (
                            arm_inputs["receiver_lines"] if arm_inputs else
                            get_species_expected_lines(lib_num)
                            if cond["mode"] in (2, 4, 5) else
                            get_expected_lines(lib_num))
                        candidates = (
                            arm_inputs["ambient_candidates"] if arm_inputs else
                            get_individual_ambient_candidates(lib_num)
                            if cond["mode"] in (1, 3, "1+SR", LEGACY2C_MODE)
                            else None)
                    except Exception:
                        el = None
                        candidates = None
                    if el is None or (
                            cond["mode"] in (1, 3, "1+SR", LEGACY2C_MODE)
                            and not (candidates and check_file_exists(candidates))):
                        q_blocked += 1
                        continue

                    script_path = generate_contam_script(
                        lib_num, cond, dep_job_id=dep_id, force=args.force,
                        assignment_source=assignment_source)
                    generated.append((
                        "CONTAM",
                        f"{abbrev}/{assignment_source}/lib{lib_num}",
                        script_path))
                    contam_generated_keys.add(
                        (abbrev, assignment_source, lib_num))
                    cond_generated += 1
                    q_count += 1

                    if args.submit:
                        jid = submit_job(script_path)
                        if jid:
                            contam_job_ids[
                                (abbrev, assignment_source, lib_num)] = jid
                            print(
                                f"  ✅ {abbrev}/{assignment_source}/lib{lib_num}: "
                                f"job {jid}")

                if cond_generated > 0:
                    print(
                        f"  {abbrev}/{assignment_source}: "
                        f"{cond_generated} jobs generated")

        print()
        print(f"  CONTAM totals: {q_count} generated, "
              f"{q_skipped} skipped (output exists), "
              f"{q_blocked} blocked (missing prereqs), "
              f"{q_not_applicable} not applicable")
        print()

    # ---- GEX_AMBIENT: gene-level ambient profile and concentration tables ----
    if "GEX_AMBIENT" in stages:
        print("--- GEX_AMBIENT: ambient gene-expression profiles ---")
        gex_runs = []
        generated_count = 0
        skipped_count = 0
        blocked_count = 0
        not_applicable_count = 0
        gex_blocked_keys = set()
        assignment_sources = resolve_contam_assignment_sources(
            args.contam_assignment_source)
        if args.gex_cluster_source in {"auto", "h5ad"}:
            if args.gex_cluster_source == "auto":
                print(
                    "  auto clusters: raw filtered MEX for every library; "
                    "writing deterministic numbered groups")
            else:
                print(
                    "  h5ad clusters: reusing --ploidy-nn-h5ad and writing "
                    "deterministic numbered groups")
            for lib_num in lib_nums:
                if (gex_auto_cluster_outputs_complete(lib_num, args) and
                        not args.force):
                    print(f"    lib{lib_num}: generated clusters already exist")
                    continue
                cluster_script = generate_gex_auto_cluster_script(
                    lib_num, args, force=args.force)
                generated.append((
                    "GEX_AUTO_CLUSTERS", f"lib{lib_num}", cluster_script))
                print(f"    Generated {args.gex_cluster_source} clusters: lib{lib_num}")
                if args.submit:
                    jid = submit_job(cluster_script)
                    if jid:
                        gex_cluster_job_ids[lib_num] = jid
                        print(f"      ✅ job {jid}")
        for cond in cond_list:
            abbrev = cond["abbrev"]
            for assignment_source in assignment_sources:
                for lib_num in lib_nums:
                    try:
                        applicable, reason = contam_source_applicable(
                            lib_num, assignment_source)
                    except Exception as exc:
                        blocked_count += 1
                        print(
                            f"  ⚠️  {abbrev}/{assignment_source}/lib{lib_num}: "
                            f"comparison plan unavailable ({exc})")
                        continue
                    if not applicable:
                        not_applicable_count += 1
                        print(
                            f"  ↪ {abbrev}/{assignment_source}/lib{lib_num}: "
                            f"not applicable ({reason})")
                        continue
                    out_prefix = get_gex_ambient_prefix(
                        lib_num, abbrev, assignment_source, args)
                    cluster_path = get_gex_cluster_path(lib_num, args)
                    run_key = f"lib{lib_num}__{abbrev}__{assignment_source}"
                    gex_runs.append({
                        "run_key": run_key,
                        "library": lib_num,
                        "condition": abbrev,
                        "assignment_source": assignment_source,
                        "cluster_source": args.gex_cluster_source,
                        "profile": out_prefix + ".gex_profile",
                        "clusters": cluster_path,
                    })

                    if gex_ambient_outputs_complete(
                            lib_num, abbrev, assignment_source, args) and not args.force:
                        skipped_count += 1
                        continue

                    dep_id = contam_job_ids.get(
                        (abbrev, assignment_source, lib_num))
                    dep_ids = [
                        job_id for job_id in (
                            dep_id, gex_cluster_job_ids.get(lib_num))
                        if job_id is not None]
                    source_prefix = get_contam_prefix(
                        lib_num, abbrev, assignment_source)
                    source_ready = all(check_file_exists(
                            source_prefix + suffix) for suffix in (
                                ".contam_rate", ".contam_prof",
                                ".assignments", ".samples"))
                    source_planned = (
                        (abbrev, assignment_source, lib_num)
                        in contam_generated_keys)
                    if dep_id is None and not source_ready and not source_planned:
                        print(
                            f"  ⚠️  {abbrev}/{assignment_source}/lib{lib_num}: "
                            "source contamination bundle is unavailable")
                        blocked_count += 1
                        gex_blocked_keys.add((abbrev, assignment_source, lib_num))
                        continue
                    script_path = generate_gex_ambient_script(
                        lib_num, cond, assignment_source, args,
                        dep_job_ids=dep_ids, force=args.force)
                    generated.append((
                        "GEX_AMBIENT",
                        f"{abbrev}/{assignment_source}/lib{lib_num}",
                        script_path))
                    generated_count += 1
                    print(
                        f"  Generated: {abbrev}/{assignment_source}/lib{lib_num}")
                    if args.submit:
                        jid = submit_job(script_path)
                        if jid:
                            gex_ambient_job_ids[
                                (abbrev, assignment_source, lib_num)] = jid
                            print(f"    ✅ job {jid}")

        if gex_blocked_keys:
            gex_runs = [
                run for run in gex_runs
                if (run["condition"], run["assignment_source"], run["library"])
                not in gex_blocked_keys
            ]

        if gex_runs:
            summary_script, _ = generate_gex_ambient_summary_script(
                args, gex_runs,
                dep_job_ids=list(gex_ambient_job_ids.values()))
            generated.append(("GEX_AMBIENT_SUMMARY", args.gex_analysis_name,
                              summary_script))
            print(f"  Summary: {get_gex_ambient_summary_dir(args)}")
            if args.submit:
                summary_jid = submit_job(summary_script)
                print(f"  ✅ GEX_AMBIENT summary job {summary_jid}")
        else:
            print("  ⚠️  No runnable or completed GEX_AMBIENT profiles to summarize")
        print(
            f"  GEX_AMBIENT totals: {generated_count} generated, "
            f"{skipped_count} complete, {blocked_count} blocked, "
            f"{not_applicable_count} not applicable")
        print()


    # ---- AMBIENT_PLOTS: selected-condition descriptive contamination QC ----
    if "AMBIENT_PLOTS" in stages:
        print("--- AMBIENT_PLOTS: descriptive contamination plots (no ranking) ---")
        plot_run_dir = get_ambient_plot_run_dir(args, cond_list)
        upstream_ids = list(contam_job_ids.values())
        bundle = ambient_generate_plot_job_bundle(
            orchestrator_path=os.path.abspath(__file__),
            mapping_root=get_demux_mapping_root(),
            aggregate_root=AGGREGATE_ROOT,
            plot_root=plot_run_dir,
            libraries=lib_nums,
            conditions=[condition["abbrev"] for condition in cond_list],
            assignment_sources=resolve_contam_assignment_sources(
                args.contam_assignment_source),
            script_dir=get_script_dir(),
            log_dir=os.path.join(get_log_dir(), "AMBIENT_PLOTS"),
            upstream_job_ids=None,
            partition=SLURM_PARTITION,
            library_prefix=LIB_PREFIX,
            demux_subdir=DEMUX_SUBDIR,
            contamination_subdir="contamination",
            deployed_contam_r=resolve_contam_r_script(),
            r_array_parallelism=None,
            plot_formats=("pdf", "png"),
            reference_condition=args.ambient_reference_condition,
            identity_ambient_candidate_set=(
                args.reconciliation_candidate_set),
        )
        generated.append(("AMBIENT_PLOTS_R", "selected_conditions", bundle["r_sbatch"]))
        generated.append(("AMBIENT_PLOTS", "selected_conditions", bundle["aggregate_sbatch"]))
        print(f"  Output: {bundle['plot_root']}")
        print(f"  contam.R tasks: {bundle['n_r_tasks']}")
        print(f"  Generated: {bundle['r_sbatch']}")
        print(f"  Generated: {bundle['aggregate_sbatch']}")
        if args.submit:
            plot_jobs = ambient_submit_plot_job_bundle(
                bundle, upstream_job_ids=upstream_ids, submit=True)
            print(f"  ✅ contam.R array: job {plot_jobs['r_job_id']}")
            print(f"  ✅ aggregate plots: job {plot_jobs['aggregate_job_id']}")
        print()

    # ---- AMBIENT_VALIDATE: production fixed-profile E/F comparison ----
    if "AMBIENT_VALIDATE" in stages:
        print("--- AMBIENT_VALIDATE: fixed-profile original-vs-reconciled assignments ---")
        validation_plans = {}
        for lib_num in lib_nums:
            plan = prepare_ambient_validation_plan(
                lib_num, args.reconciliation_candidate_set)
            validation_plans[lib_num] = plan
            if not plan.get("applicable"):
                print(
                    f"  ↪ lib{lib_num}: not applicable "
                    f"({plan['skip_reason']})")
                continue
            profile_script = generate_ambient_validation_profile_script(
                lib_num, plan, force=args.force)
            generated.append((
                "AMBIENT_VALIDATE_PROFILE", f"lib{lib_num}", profile_script))
            print(
                f"  lib{lib_num}: targets={','.join(plan['targets_list'])}; "
                f"production donors={len(plan['production_roster_list'])}")
            if args.submit:
                jid = submit_job(profile_script)
                ambient_validation_profile_job_ids[lib_num] = jid
                print(f"  ✅ fixed profile lib{lib_num}: job {jid}")

        for lib_num, plan in validation_plans.items():
            if not plan.get("applicable"):
                continue
            for cond in cond_list:
                fixed_script = generate_ambient_validation_fixed_script(
                    lib_num, cond, args, plan,
                    dep_job_id=ambient_validation_profile_job_ids.get(lib_num),
                    force=args.force)
                generated.append((
                    "AMBIENT_VALIDATE_FIXED",
                    f"{cond['abbrev']}/lib{lib_num}", fixed_script))
                if args.submit:
                    jid = submit_job(fixed_script)
                    ambient_validation_fixed_job_ids.append(jid)
                    print(
                        f"  ✅ fixed-profile E/F {cond['abbrev']}/lib{lib_num}: "
                        f"job {jid}")

        aggregate_bundle = generate_ambient_validation_aggregate_script(
            args, cond_list, validation_plans,
            fixed_job_ids=(
                ambient_validation_fixed_job_ids
                if args.submit else None))
        generated.append((
            "AMBIENT_VALIDATE", "selected_libraries",
            aggregate_bundle["script"]))
        print(f"  Output: {aggregate_bundle['output_root']}")
        print(f"  Generated: {aggregate_bundle['script']}")
        if args.submit:
            jid = submit_job(aggregate_bundle["script"])
            print(f"  ✅ validation aggregate: job {jid}")
        print()

    # ---- AMBIENT_SWAP_TEST: auto-discovered fixed G/H and joint J/K comparison ----
    if "AMBIENT_SWAP_TEST" in stages:
        print("--- AMBIENT_SWAP_TEST: auto-discovered original-vs-proposed assignments ---")
        swap_discovery = discover_ambient_swap_events(args)
        print(
            f"  Discovery: {swap_discovery['selected_event_count']} supported "
            f"unexpected events across "
            f"{swap_discovery['selected_library_count']} libraries; "
            f"event_min_cells={swap_discovery['event_min_cells']}")
        print(
            "  Cell gate: diploid-to-heterotypic proposals require "
            f"NN QC PASS and P(tet)>="
            f"{AMBIENT_SWAP_HETEROTYPIC_NN_MIN_PROB:.2f}; "
            "event mass is rechecked afterward")
        swap_plans = {}
        for lib_num in lib_nums:
            plan = prepare_ambient_swap_test_plan(
                lib_num, args, discovery=swap_discovery)
            swap_plans[lib_num] = plan
            if not plan.get("applicable"):
                print(
                    f"  ↪ lib{lib_num}: not applicable "
                    f"({plan['skip_reason']})")
                continue
            profile_script = generate_ambient_swap_profile_script(
                lib_num, plan, force=args.force)
            generated.append((
                "AMBIENT_SWAP_PROFILE", f"lib{lib_num}", profile_script))
            print(
                f"  lib{lib_num}: candidates="
                f"{','.join(plan['candidate_identities'])}; "
                f"selected cells={plan['n_selected_cells']}; "
                f"source controls={plan['n_source_controls']}; "
                f"profile donors={len(plan['production_roster_list'])}")
            if args.submit:
                jid = submit_job(profile_script)
                ambient_swap_profile_job_ids[lib_num] = jid
                print(f"  ✅ discovered-candidate profile lib{lib_num}: job {jid}")

        for lib_num, plan in swap_plans.items():
            if not plan.get("applicable"):
                continue
            for cond in cond_list:
                arm_script = generate_ambient_swap_arm_script(
                    lib_num, cond, plan,
                    dep_job_id=ambient_swap_profile_job_ids.get(lib_num),
                    force=args.force)
                generated.append((
                    "AMBIENT_SWAP_ARMS",
                    f"{cond['abbrev']}/lib{lib_num}", arm_script))
                if args.submit:
                    jid = submit_job(arm_script)
                    ambient_swap_arm_job_ids.append(jid)
                    print(
                        f"  ✅ fixed-profile G/H + joint-profile J/K "
                        f"{cond['abbrev']}/lib{lib_num}: job {jid}")

        aggregate_bundle = generate_ambient_swap_aggregate_script(
            args, cond_list, swap_plans, swap_discovery,
            arm_job_ids=(
                ambient_swap_arm_job_ids if args.submit else None))
        generated.append((
            "AMBIENT_SWAP_TEST", "selected_libraries",
            aggregate_bundle["script"]))
        print(f"  Output: {aggregate_bundle['output_root']}")
        print(f"  Generated: {aggregate_bundle['script']}")
        if args.submit:
            ambient_swap_aggregate_job_id = submit_job(
                aggregate_bundle["script"])
            print(
                f"  ✅ swap-test aggregate: job "
                f"{ambient_swap_aggregate_job_id}")
        print()

    # ---- Stage 5 PLOIDY_NN: optional NN ploidy inference ----
    if "PLOIDY_NN" in stages:
        print("--- Stage 5 PLOIDY_NN: optional NN ploidy inference ---")
        missing = [lib for lib in lib_nums if not check_file_exists(get_ploidy_calls_path(lib))]
        if not missing and not args.force:
            print(f"  ✅ ploidy NN outputs already exist for {len(lib_nums)} libraries, skipping")
        else:
            script_path = generate_ploidy_nn_script(lib_nums, force=args.force, args=args)
            generated.append(("PLOIDY_NN", f"libs{min(lib_nums)}-{max(lib_nums)}", script_path))
            print(f"  Generated: libraries {', '.join('lib' + str(x) for x in lib_nums)}")
            if args.submit:
                jid = submit_job(script_path)
                if jid:
                    ploidy_nn_job_id = jid
                    print(f"  ✅ PLOIDY_NN: job {jid}")
        print()

    # ---- Stage 6 TETRA_REFINE: optional biological assignment refinement ----
    if "TETRA_REFINE" in stages:
        print("--- Stage 6 TETRA_REFINE: optional biological assignment refinement ---")
        for lib_num in lib_nums:
            if tetra_refine_outputs_complete(lib_num) and not args.force:
                print(f"  ✅ lib{lib_num}: tetra_refine outputs exist, skipping")
                continue

            deps = []
            if ploidy_nn_job_id:
                deps.append(ploidy_nn_job_id)
            if demux_job_ids.get(lib_num):
                deps.append(demux_job_ids[lib_num])
            if args.refine_contam_condition:
                jid = contam_job_ids.get(
                    (args.refine_contam_condition, "demux", lib_num))
                if jid:
                    deps.append(jid)

            if "DEMUX" not in stages:
                prefix = get_demux_prefix(lib_num)
                needed = [prefix + x for x in [".assignments", ".diagnostics.gz", ".runner_ups.gz"]]
                if any(not check_file_exists(x) for x in needed):
                    print(f"  ⚠️  lib{lib_num}: TETRA_REFINE prerequisites missing, skipping")
                    continue

            script_path = generate_tetra_refine_script(
                lib_num, dep_job_ids=deps, force=args.force, args=args)
            generated.append(("TETRA_REFINE", f"lib{lib_num}", script_path))
            print(f"  Generated: lib{lib_num}")
            if args.submit:
                jid = submit_job(script_path)
                if jid:
                    tetra_refine_job_ids[lib_num] = jid
                    print(f"  ✅ lib{lib_num}: job {jid}")
        print()

    # ---- Stage 7 POSTHOC: Per-cell scoring + swap audit ----
    if "POSTHOC" in stages:
        print("--- Stage 7 POSTHOC: tetra_score_calls + sample-swap audit ---")
        for lib_num in lib_nums:
            if (check_file_exists(get_call_qc_path(lib_num)) and
                check_file_exists(get_species_qc_path(lib_num)) and
                check_file_exists(get_swap_report_path(lib_num)) and
                not args.force):
                print(f"  ✅ lib{lib_num}: POSTHOC outputs exist, skipping")
                continue

            dep_id = tetra_refine_job_ids.get(lib_num) or demux_job_ids.get(lib_num)
            if dep_id is None and "DEMUX" not in stages:
                prefix = get_demux_prefix(lib_num)
                needed = [prefix + x for x in [
                    ".counts", ".condf", ".samples", ".assignments", ".diagnostics.gz", ".runner_ups.gz",
                    ".species_counts", ".species_condf", ".species_samples"]]
                if any(not check_file_exists(x) for x in needed):
                    print(f"  ⚠️  lib{lib_num}: POSTHOC prerequisites missing, skipping")
                    continue

            script_path = generate_posthoc_script(lib_num, dep_job_id=dep_id, force=args.force)
            generated.append(("POSTHOC", f"lib{lib_num}", script_path))
            print(f"  Generated: lib{lib_num}")
            if args.submit:
                jid = submit_job(script_path)
                if jid:
                    posthoc_job_ids[lib_num] = jid
                    print(f"  ✅ lib{lib_num}: job {jid}")

        # Shared aggregate files are written exactly once, after every POSTHOC
        # job submitted by this invocation has completed successfully.
        if posthoc_job_ids:
            aggregate_script = generate_posthoc_aggregate_script(
                dep_job_ids=list(posthoc_job_ids.values()) if args.submit else None)
            generated.append(("POSTHOC_AGG", "all_libraries", aggregate_script))
            print("  Generated: final all-library POSTHOC aggregate")
            if args.submit:
                jid = submit_job(aggregate_script)
                print(f"  ✅ POSTHOC aggregate: job {jid}")
        print()

    # ---- Stage 7b POSTHOC_SUMMARY: refresh verdicts from existing audit outputs ----
    if "POSTHOC_SUMMARY" in stages:
        print("--- Stage 7b POSTHOC_SUMMARY: refresh swap-audit summaries only ---")
        for lib_num in lib_nums:
            manifest = os.path.join(get_audit_lib_dir(lib_num), f"lib{lib_num}.capabilities.json")
            if not check_file_exists(manifest):
                print(f"  WARNING lib{lib_num}: capabilities manifest missing, skipping")
                continue

            script_path = generate_posthoc_summary_script(lib_num)
            generated.append(("POSTHOC_SUMMARY", f"lib{lib_num}", script_path))
            print(f"  Generated: lib{lib_num}")
            if args.submit:
                jid = submit_job(script_path)
                if jid:
                    posthoc_summary_job_ids[lib_num] = jid
                    print(f"  OK lib{lib_num}: job {jid}")

        if posthoc_summary_job_ids or not args.submit:
            aggregate_script = generate_posthoc_aggregate_script(
                dep_job_ids=list(posthoc_summary_job_ids.values()) if args.submit else None)
            generated.append(("POSTHOC_SUMMARY_AGG", "all_libraries", aggregate_script))
            print("  Generated: final all-library aggregate after refreshed summaries")
            if args.submit:
                jid = submit_job(aggregate_script)
                print(f"  OK POSTHOC_SUMMARY aggregate: job {jid}")
        print()

    # ---- Stage 7c UNEXPECTED_COMPONENT_NN: genotype alerts x NN P(tet) ----
    if "UNEXPECTED_COMPONENT_NN" in stages:
        print("--- Stage 7c UNEXPECTED_COMPONENT_NN: supported unexpected genotype cells x NN ploidy ---")
        uc_job_ids = unexpected_component_nn_job_ids
        for lib_num in lib_nums:
            if not has_unexpected_component_signal(lib_num):
                print(f"  lib{lib_num}: no UNEXPECTED_COMPONENT_SIGNAL, skipping")
                continue
            dep_id = posthoc_summary_job_ids.get(lib_num) or posthoc_job_ids.get(lib_num)
            script_path = generate_unexpected_component_nn_script(lib_num, dep_job_id=dep_id)
            generated.append(("UNEXPECTED_COMPONENT_NN", f"lib{lib_num}", script_path))
            print(f"  Generated: lib{lib_num}")
            if args.submit:
                jid = submit_job(script_path)
                if jid:
                    uc_job_ids[lib_num] = jid
                    print(f"  OK lib{lib_num}: job {jid}")

        if uc_job_ids or not args.submit:
            agg_script = generate_unexpected_component_nn_aggregate_script(
                dep_job_ids=list(uc_job_ids.values()) if args.submit else None,
                libraries=lib_nums)
            generated.append(("UNEXPECTED_COMPONENT_NN_AGG", "all_libraries", agg_script))
            print("  Generated: unexpected-component NN aggregate")
            if args.submit:
                jid = submit_job(agg_script)
                print(f"  OK UNEXPECTED_COMPONENT_NN aggregate: job {jid}")
        print()

    # ---- Stage 8 HYBRID: Post-hoc reconciliation ----
    if "HYBRID" in stages:
        print("--- Stage 8 HYBRID: post-hoc individual/species reconciliation ---")
        for lib_num in lib_nums:
            if check_file_exists(get_hybrid_summary_path(lib_num)) and not args.force:
                print(f"  ✅ lib{lib_num}: HYBRID output exists, skipping")
                continue

            deps = []
            if posthoc_job_ids.get(lib_num):
                deps.append(posthoc_job_ids[lib_num])
            for cond_name in (args.hybrid_individual_condition, args.hybrid_species_condition, args.hybrid_fixed_species_condition):
                if cond_name:
                    jid = contam_job_ids.get((cond_name, "demux", lib_num))
                    if jid:
                        deps.append(jid)

            # Validate every prerequisite independently. One available producer
            # must not mask an unrelated missing condition or POSTHOC product.
            missing_hybrid = []
            if not posthoc_job_ids.get(lib_num):
                for path in (get_swap_report_path(lib_num), get_call_qc_path(lib_num),
                             get_species_qc_path(lib_num)):
                    if not check_file_exists(path):
                        missing_hybrid.append(path)
            for cond_name, needs_profile in (
                    (args.hybrid_individual_condition, False),
                    (args.hybrid_species_condition, True),
                    (args.hybrid_fixed_species_condition, True)):
                if not cond_name or contam_job_ids.get(
                        (cond_name, "demux", lib_num)):
                    continue
                prefix = get_contam_prefix(lib_num, cond_name)
                if not check_file_exists(prefix + ".contam_rate"):
                    missing_hybrid.append(prefix + ".contam_rate")
                if (needs_profile and
                        not (check_file_exists(prefix + ".species_prof") or
                             check_file_exists(prefix + ".contam_prof"))):
                    missing_hybrid.append(prefix + ".species_prof")
            if missing_hybrid:
                print(f"  ⚠️  lib{lib_num}: HYBRID prerequisites missing, skipping")
                for path in missing_hybrid:
                    print(f"       {path}")
                continue

            script_path = generate_hybrid_script(
                lib_num, dep_job_ids=deps, force=args.force,
                individual_condition=args.hybrid_individual_condition,
                species_condition=args.hybrid_species_condition,
                fixed_species_condition=args.hybrid_fixed_species_condition)
            generated.append(("HYBRID", f"lib{lib_num}", script_path))
            print(f"  Generated: lib{lib_num}")
            if args.submit:
                jid = submit_job(script_path)
                if jid:
                    hybrid_job_ids[lib_num] = jid
                    print(f"  ✅ lib{lib_num}: job {jid}")
        print()

    # ---- Stage 8a IDENTITY_SCORE_AGGREGATE_ONLY: frozen-score v6.1 reporting ----
    if "IDENTITY_SCORE_AGGREGATE_ONLY" in stages:
        print("--- Stage 8a IDENTITY_SCORE_AGGREGATE_ONLY: frozen pairs/probabilities -> QC-scoped aggregate ---")
        aggregate_script = generate_identity_probability_aggregate_script(
            lib_nums, args, reuse_frozen=True)
        generated.append((
            "IDENTITY_PROBABILITY_AGGREGATE", "selected_libraries",
            aggregate_script))
        if args.submit:
            identity_probability_aggregate_job_id = submit_job(
                aggregate_script)
            print(
                "  IDENTITY_PROBABILITY_AGGREGATE: job "
                f"{identity_probability_aggregate_job_id}")
        print()

    # ---- Stage 8b IDENTITY_SCORE: post-reconciliation targeted probabilities ----
    if "IDENTITY_SCORE" in stages:
        print("--- Stage 8b IDENTITY_SCORE: reconciled proposals -> exact original/swap pairs -> probabilities -> aggregate ---")
        pair_script = generate_identity_score_pair_manifest_script(
            lib_nums, args)
        generated.append((
            "IDENTITY_SCORE_PAIRS", "selected_libraries", pair_script))
        if args.submit:
            identity_score_pair_manifest_job_id = submit_job(pair_script)
            print(
                "  IDENTITY_SCORE_PAIRS: job "
                f"{identity_score_pair_manifest_job_id}")
        for lib_num in lib_nums:
            pair_dependency = (
                [identity_score_pair_manifest_job_id]
                if identity_score_pair_manifest_job_id else None)
            score_script = generate_identity_probability_score_script(
                lib_num, args,
                dep_job_ids=pair_dependency)
            generated.append((
                "IDENTITY_PROBABILITY_SCORE", f"lib{lib_num}",
                score_script))
            if args.submit:
                job_id = submit_job(score_script)
                if job_id:
                    identity_probability_job_ids[lib_num] = job_id
                print(f"  lib{lib_num}: targeted probability job {job_id}")

        aggregate_deps = list(identity_probability_job_ids.values())
        if identity_score_pair_manifest_job_id:
            aggregate_deps.append(identity_score_pair_manifest_job_id)
        aggregate_script = generate_identity_probability_aggregate_script(
            lib_nums, args,
            dep_job_ids=aggregate_deps if args.submit else None)
        generated.append((
            "IDENTITY_PROBABILITY_AGGREGATE", "selected_libraries",
            aggregate_script))
        if args.submit:
            identity_probability_aggregate_job_id = submit_job(
                aggregate_script)
            print(
                "  IDENTITY_PROBABILITY_AGGREGATE: job "
                f"{identity_probability_aggregate_job_id}")
        print()

    # ---- Stage 9 IDENTITY_RECONCILIATION: canonical identity framework ----
    if "IDENTITY_RECONCILIATION" in stages:
        print("--- Stage 9 IDENTITY_RECONCILIATION: preliminary reconciliation -> frozen candidate axis + ambient evidence -> canonical finalization ---")
        metadata_script = generate_identity_metadata_script(lib_nums, args)
        generated.append(("IDENTITY_METADATA", "selected_libraries", metadata_script))
        if args.submit:
            identity_metadata_job_id = submit_job(metadata_script)
            print(f"  IDENTITY_METADATA: job {identity_metadata_job_id}")

        # Candidate construction is independent by library once the shared
        # metadata contract exists.  Each candidate job depends only on shared
        # prerequisites plus that library's own POSTHOC job (when POSTHOC is
        # being regenerated in this invocation).  This removes the previous
        # all-library candidate barrier before nuclear/MT scoring.
        for lib_num in lib_nums:
            candidate_deps = []
            if demux_job_ids.get(lib_num):
                candidate_deps.append(demux_job_ids[lib_num])
            if identity_metadata_job_id:
                candidate_deps.append(identity_metadata_job_id)
            if ploidy_nn_job_id:
                candidate_deps.append(ploidy_nn_job_id)
            if posthoc_job_ids.get(lib_num):
                candidate_deps.append(posthoc_job_ids[lib_num])

            candidate_script = generate_identity_candidates_script(
                lib_num, args, dep_job_ids=candidate_deps)
            generated.append((
                "IDENTITY_CANDIDATES", f"lib{lib_num}", candidate_script))
            if args.submit:
                jid = submit_job(candidate_script)
                identity_candidates_job_ids[lib_num] = jid
                print(f"  lib{lib_num}: identity candidates job {jid}")

        candidate_aggregate_script = generate_identity_candidates_aggregate_script(
            lib_nums, args,
            dep_job_ids=list(identity_candidates_job_ids.values())
            if args.submit else None)
        generated.append((
            "IDENTITY_CANDIDATES_AGGREGATE", "selected_libraries",
            candidate_aggregate_script))
        if args.submit:
            identity_candidates_aggregate_job_id = submit_job(
                candidate_aggregate_script)
            print(
                f"  IDENTITY_CANDIDATES_AGGREGATE: job "
                f"{identity_candidates_aggregate_job_id}")

        doublet_deps = [identity_metadata_job_id] if identity_metadata_job_id else []
        doublet_script = generate_identity_doublet_context_script(
            lib_nums, args, dep_job_ids=doublet_deps)
        generated.append(("IDENTITY_DOUBLET_CONTEXT", "selected_libraries", doublet_script))
        if args.submit:
            identity_doublet_context_job_id = submit_job(doublet_script)
            print(f"  IDENTITY_DOUBLET_CONTEXT: job {identity_doublet_context_job_id}")

        for lib_num in lib_nums:
            candidate_job_id = identity_candidates_job_ids.get(lib_num)
            deps = [candidate_job_id] if candidate_job_id else []
            nuclear_script = generate_identity_nuclear_script(
                lib_num, args, dep_job_ids=deps)
            mt_identity_script = generate_identity_mt_script(
                lib_num, args, dep_job_ids=deps)
            generated.append((
                "NUCLEAR_IDENTITY_SCORE", f"lib{lib_num}", nuclear_script))
            generated.append((
                "MT_IDENTITY_SCORE", f"lib{lib_num}", mt_identity_script))
            lib_score_jobs = []
            if args.submit:
                lib_score_jobs.append(submit_job(nuclear_script))
                lib_score_jobs.append(submit_job(mt_identity_script))
            if args.identity_evidence_mode == "rna-atac":
                atac_script = generate_identity_atac_script(
                    lib_num, args, dep_job_ids=deps)
                generated.append((
                    "ATAC_IDENTITY_SCORE", f"lib{lib_num}", atac_script))
                if args.submit:
                    lib_score_jobs.append(submit_job(atac_script))
            if args.submit:
                identity_score_job_ids[lib_num] = [
                    j for j in lib_score_jobs if j]
                print(
                    f"  lib{lib_num}: identity score jobs "
                    + ",".join(identity_score_job_ids[lib_num]))

        all_score_jobs = [
            j for jobs in identity_score_job_ids.values() for j in jobs]
        reconcile_deps = list(all_score_jobs)
        if identity_candidates_aggregate_job_id:
            reconcile_deps.append(identity_candidates_aggregate_job_id)
        if identity_doublet_context_job_id:
            reconcile_deps.append(identity_doublet_context_job_id)
        reconcile_script = generate_identity_reconcile_script(
            lib_nums, args, dep_job_ids=reconcile_deps)
        generated.append((
            "IDENTITY_RECONCILE", "selected_libraries", reconcile_script))
        if args.submit:
            identity_reconcile_job_id = submit_job(reconcile_script)
            print(f"  IDENTITY_RECONCILE: job {identity_reconcile_job_id}")
        validate_script = generate_identity_validate_script(
            lib_nums, args,
            dep_job_ids=[identity_reconcile_job_id]
            if identity_reconcile_job_id else None)
        generated.append((
            "IDENTITY_VALIDATE", "selected_libraries", validate_script))
        if args.submit:
            identity_validate_job_id = submit_job(validate_script)
            print(f"  IDENTITY_VALIDATE: job {identity_validate_job_id}")

        final_evidence_script = generate_identity_final_evidence_planner_script(
            lib_nums, args,
            dep_job_ids=[identity_validate_job_id]
            if identity_validate_job_id else None)
        generated.append((
            "IDENTITY_FINAL_EVIDENCE", "selected_libraries",
            final_evidence_script))
        if args.submit:
            identity_final_evidence_job_id = submit_job(final_evidence_script)
            print(
                "  IDENTITY_FINAL_EVIDENCE (event-aware continuation): job "
                f"{identity_final_evidence_job_id}")
        print()

    # ---- Stage 9a IDENTITY_FINAL_EVIDENCE: resume validated evidence tail ----
    if "IDENTITY_FINAL_EVIDENCE" in stages:
        print("--- Stage 9a IDENTITY_FINAL_EVIDENCE: reuse validated reconciliation -> complete evidence and canonical finalization ---")
        final_evidence_script = generate_identity_final_evidence_planner_script(
            lib_nums, args)
        generated.append((
            "IDENTITY_FINAL_EVIDENCE", "selected_libraries",
            final_evidence_script))
        if args.submit:
            identity_final_evidence_job_id = submit_job(final_evidence_script)
            print(
                "  IDENTITY_FINAL_EVIDENCE (event-aware continuation): job "
                f"{identity_final_evidence_job_id}")
        print()

    # ---- Stage 9b IDENTITY_FINAL_EVIDENCE_ONLY: subset-only aggregation ----
    if "IDENTITY_FINAL_EVIDENCE_ONLY" in stages:
        print("--- Stage 9b IDENTITY_FINAL_EVIDENCE_ONLY: reuse completed evidence -> selected-library aggregates and finalization ---")
        final_evidence_script = generate_identity_final_evidence_planner_script(
            lib_nums, args, reuse_only=True)
        generated.append((
            "IDENTITY_FINAL_EVIDENCE_ONLY", "selected_libraries",
            final_evidence_script))
        if args.submit:
            identity_final_evidence_job_id = submit_job(final_evidence_script)
            print(
                "  IDENTITY_FINAL_EVIDENCE_ONLY (reuse-only continuation): job "
                f"{identity_final_evidence_job_id}")
        print()

    # ---- Stage 9c IDENTITY_FINALIZE_ONLY: completed evidence -> finalizer ----
    if "IDENTITY_FINALIZE_ONLY" in stages:
        print("--- Stage 9c IDENTITY_FINALIZE_ONLY: completed evidence -> canonical finalization only ---")
        evidence_roots = resolve_identity_finalize_only_evidence_roots(
            lib_nums, args, require_completed=True)
        finalizer_script = generate_identity_finalize_script(
            lib_nums, args,
            candidate_axis_root=evidence_roots["candidate_axis_root"],
            frozen_ambient_root=evidence_roots["frozen_ambient_root"],
            four_arm_root=evidence_roots["four_arm_root"],
            dep_job_ids=None)
        generated.append((
            "IDENTITY_FINALIZE", "selected_libraries", finalizer_script))
        if args.submit:
            finalizer_job_id = submit_job(finalizer_script)
            print(f"  IDENTITY_FINALIZE: job {finalizer_job_id}")
        print()

    # ---- Stage 9d IDENTITY_RECONCILE_ONLY: reuse existing evidence ----
    if "IDENTITY_RECONCILE_ONLY" in stages:
        print("--- Stage 9b IDENTITY_RECONCILE_ONLY: existing evidence -> decisions -> validation ---")
        reconcile_script = generate_identity_reconcile_script(lib_nums, args)
        generated.append((
            "IDENTITY_RECONCILE", "selected_libraries", reconcile_script))
        if args.submit:
            identity_reconcile_job_id = submit_job(reconcile_script)
            print(f"  IDENTITY_RECONCILE: job {identity_reconcile_job_id}")
        validate_script = generate_identity_validate_script(
            lib_nums, args,
            dep_job_ids=[identity_reconcile_job_id]
            if identity_reconcile_job_id else None)
        generated.append((
            "IDENTITY_VALIDATE", "selected_libraries", validate_script))
        if args.submit:
            identity_validate_job_id = submit_job(validate_script)
            print(f"  IDENTITY_VALIDATE: job {identity_validate_job_id}")
        print()

    # ---- Stage 9 MT_FUSION: per-cell mitochondrial ratios ----
    if "MT_FUSION" in stages:
        print("--- Stage 9 MT_FUSION: per-cell mitochondrial ratios ---")
        for lib_num in mt_nonapplicable_libs:
            print(
                f"  lib{lib_num}: no rows in frozen MT ratio site manifest; "
                "skipping MT_FUSION as not applicable")
        for lib_num in mt_analysis_libs:
            mt_reconciliation_refresh = (
                args.mt_assignment_source in ("auto", "reconciled") and
                identity_validate_job_id is not None)
            mt_ambient_job_id = None
            if args.mt_rna_ambient_condition:
                mt_ambient_job_id = get_mt_rna_ambient_job_id(
                    lib_num, args, contam_job_ids)
            mt_ambient_refresh = mt_ambient_job_id is not None
            if (mt_outputs_complete(lib_num, args) and
                    not args.force and not mt_reconciliation_refresh and
                    not mt_ambient_refresh):
                print(f"  ✅ lib{lib_num}: MT_FUSION outputs exist, skipping")
                continue

            deps = []
            if mt_reconciliation_refresh:
                deps.append(identity_validate_job_id)
            if mt_ambient_refresh:
                deps.append(mt_ambient_job_id)
            if not deps:
                for jid in (
                    hybrid_job_ids.get(lib_num),
                    unexpected_component_nn_job_ids.get(lib_num),
                    posthoc_summary_job_ids.get(lib_num),
                    posthoc_job_ids.get(lib_num),
                    tetra_refine_job_ids.get(lib_num),
                    demux_job_ids.get(lib_num),
                ):
                    if jid:
                        deps.append(jid)
                        break

            script_path = generate_mt_fusion_script(
                lib_num, dep_job_ids=deps,
                force=(args.force or mt_reconciliation_refresh or
                       mt_ambient_refresh), args=args)
            generated.append(("MT_FUSION", f"lib{lib_num}", script_path))
            mt_fusion_changed = True
            print(f"  Generated: lib{lib_num}")
            if args.submit:
                jid = submit_job(script_path)
                if jid:
                    mt_fusion_job_ids[lib_num] = jid
                    print(f"  ✅ lib{lib_num}: job {jid}")
        print()

    # ---- Stage 10 MT_POPULATION: optional metadata-driven population structure ----
    if "MT_POPULATION" in stages:
        print("--- Stage 10 MT_POPULATION: mitochondrial population structure ---")
        if not mt_analysis_libs:
            print("  No selected libraries have rows in the frozen MT ratio site manifest; skipping MT_POPULATION")
        else:
            population_deps = list(mt_fusion_job_ids.values())
            if identity_validate_job_id:
                population_deps.append(identity_validate_job_id)
            population_script = generate_mt_population_script(
                mt_analysis_libs,
                dep_job_ids=population_deps if args.submit else None,
                force=(args.force or mt_fusion_changed or
                       identity_validate_job_id is not None),
                args=args)
            generated.append(("MT_POPULATION", "selected_libraries", population_script))
            print(
                "  Generated: mitochondrial population-structure job over "
                f"{len(mt_analysis_libs)} manifest-applicable libraries")
            if args.submit:
                jid = submit_job(population_script)
                print(f"  ✅ MT population: job {jid}")
        print()

    # ---- Summary ----
    print("=" * 72)
    print("SUMMARY")
    stage_counts = {}
    for stage, label, _ in generated:
        stage_counts[stage] = stage_counts.get(stage, 0) + 1

    for stage in ["CLEANUP_RESULTS", "VCF_DAEMON_HOLDER", "CONDF", "DEMUX", "VCF_DAEMON_CLEANUP", "EMPTY_DROPS", "CONTAM", "GEX_AUTO_CLUSTERS", "GEX_AMBIENT", "GEX_AMBIENT_SUMMARY", "AMBIENT_PLOTS_R", "AMBIENT_PLOTS", "AMBIENT_VALIDATE_PROFILE", "AMBIENT_VALIDATE_FIXED", "AMBIENT_VALIDATE", "AMBIENT_SWAP_PROFILE", "AMBIENT_SWAP_ARMS", "AMBIENT_SWAP_TEST", "PLOIDY_NN", "TETRA_REFINE", "POSTHOC", "POSTHOC_AGG", "POSTHOC_SUMMARY", "POSTHOC_SUMMARY_AGG", "UNEXPECTED_COMPONENT_NN", "UNEXPECTED_COMPONENT_NN_AGG", "HYBRID", "IDENTITY_METADATA", "IDENTITY_CANDIDATES", "IDENTITY_CANDIDATES_AGGREGATE", "IDENTITY_DOUBLET_CONTEXT", "IDENTITY_SCORE_PAIRS", "IDENTITY_PROBABILITY_SCORE", "NUCLEAR_IDENTITY_SCORE", "MT_IDENTITY_SCORE", "ATAC_IDENTITY_SCORE", "IDENTITY_RECONCILE", "IDENTITY_VALIDATE", "IDENTITY_FINAL_EVIDENCE", "IDENTITY_FINAL_EVIDENCE_ONLY", "IDENTITY_FINALIZE", "IDENTITY_PROBABILITY_AGGREGATE", "MT_FUSION", "MT_POPULATION"]:
        if stage in stage_counts:
            script_label = (
                "script" if stage_counts[stage] == 1 else "scripts")
            print(
                f"  {stage}: {stage_counts[stage]} {script_label} generated")

    total_label = "script" if len(generated) == 1 else "scripts"
    print(f"  Total: {len(generated)} {total_label}")

    if args.submit:
        # Count generated/submitted scripts by stage. Dependency-held jobs are
        # still valid submissions. This is a bookkeeping summary only; sbatch
        # job IDs printed above remain the source of truth.
        n_submitted = sum(1 for _, _, _ in generated)
        job_label = "job" if n_submitted == 1 else "jobs"
        print(f"  Submitted: {n_submitted} {job_label}")
        print(f"  Monitor with: squeue -u $USER")
    else:
        print(f"  Dry run complete. Re-run with --submit to launch jobs.")

    print(f"  Scripts in: {get_script_dir()}/")
    print(f"  Logs in: {get_log_dir()}/")
    print("=" * 72)


# =============================================================================
# CLI
# =============================================================================

def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Orchestrate the production tetraploid CellBouncer workflow.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  orchestrate_tetraploid.py --libraries 13
  orchestrate_tetraploid.py --diagnose-only
  orchestrate_tetraploid.py --libraries 1-40 --stage CONTAM --submit
  orchestrate_tetraploid.py --libraries 1-40 --stage AMBIENT_PLOTS --submit
  orchestrate_tetraploid.py --libraries 19 --stage CONTAM,AMBIENT_PLOTS --contam-assignment-source reconciliation-four-arm --submit
  orchestrate_tetraploid.py --libraries 1-40 --stage AMBIENT_VALIDATE --submit
  orchestrate_tetraploid.py --libraries 1-40 --stage AMBIENT_SWAP_TEST --submit
  orchestrate_tetraploid.py --libraries 19 --stage GEX_AMBIENT --condition-set ck-minimal
  orchestrate_tetraploid.py --stage CLEANUP_RESULTS --submit
  orchestrate_tetraploid.py --condition-set ck-minimal --stage CONTAM,AMBIENT_PLOTS --submit
  orchestrate_tetraploid.py --condition-set all --stage CONTAM --submit
  orchestrate_tetraploid.py --stage DEMUX --libraries 1 2 --submit
  orchestrate_tetraploid.py --stage MT_FUSION --libraries 19 --submit
  orchestrate_tetraploid.py --stage IDENTITY_SCORE --libraries 1-40 --submit
  orchestrate_tetraploid.py --stage IDENTITY_SCORE_AGGREGATE_ONLY --libraries 19 --identity-score-output-root /path/to/new_v6_1_summary --submit
  orchestrate_tetraploid.py --stage IDENTITY_CANDIDATE_AXIS --libraries 1 8 19 --identity-candidate-axis-output-root /absolute/candidate_axis_round2 --submit
  orchestrate_tetraploid.py --stage IDENTITY_CANDIDATE_AXIS --libraries 19 --identity-candidate-axis-event-id EVENT --identity-candidate-axis-proposal C40210+H27322 --identity-candidate-axis-output-root /absolute/candidate_axis_targeted --submit
  orchestrate_tetraploid.py --stage IDENTITY_FINAL_EVIDENCE --libraries 7 11 12 19 36 --identity-reconciliation-root /absolute/identity_run --submit
  orchestrate_tetraploid.py --stage IDENTITY_FINAL_EVIDENCE_ONLY --libraries 1-37 39 --identity-reconciliation-root /absolute/identity_run --submit
  orchestrate_tetraploid.py --stage IDENTITY_FINALIZE_ONLY --libraries 1-37 39 --identity-reconciliation-root /absolute/identity_run --submit
  orchestrate_tetraploid.py --stage MT_FUSION --libraries 19 --mt-rna-ambient-condition IND_CK_RF_SX0_GATED_RFREE_PFIT --submit

Stages (run in order):
  CLEANUP_RESULTS = standalone cleanup retaining newest completed generated results
  CONDF       = .condf generation (3 jobs total, library-independent)
  DEMUX       = demux (1 job per library; managed VCF holders coexist by reservation)
  EMPTY_DROPS = empty drops ambient profile (1 job per library)
  CONTAM      = contamination estimation (1 job per condition x library)
  GEX_AMBIENT = infer ambient gene profiles with RNA-Leiden, H5AD-column, or manual clusters
  AMBIENT_PLOTS = selected-condition descriptive plots plus CellBouncer contam.R
  AMBIENT_VALIDATE = fixed-profile original-vs-reconciled assignment validation
  AMBIENT_SWAP_TEST = auto-discovered fixed and joint-profile swap test; figures by identity_reconciliation_figures.py
  PLOIDY_NN   = optional NN ploidy inference for homotypic/tetraploid evidence
  TETRA_REFINE= optional biological assignment refinement using expected lines and NN calls
  POSTHOC     = tetra_score_calls + swap audit from demux/refined outputs
  POSTHOC_SUMMARY = refresh swap-audit verdicts from existing POSTHOC outputs only
  UNEXPECTED_COMPONENT_NN = join supported unexpected-component cells to existing NN P(tet) calls
  HYBRID      = post-hoc reconciliation over independent native outputs
  IDENTITY_SCORE = post-reconciliation original-vs-nominated-swap probabilities and aggregation
  IDENTITY_SCORE_AGGREGATE_ONLY = reuse frozen pair/probability files; run QC-scoped aggregation only
  IDENTITY_CANDIDATE_AXIS = standalone finalized-event fixed-pair geometric diagnostic
  IDENTITY_RECONCILIATION = preliminary reconciliation, frozen axis/ambient evidence, and canonical finalization
  IDENTITY_FINAL_EVIDENCE = resume validated reconciliation at its evidence/finalization tail
  IDENTITY_FINAL_EVIDENCE_ONLY = reuse completed evidence; rerun selected-library aggregates and finalization only
  IDENTITY_FINALIZE_ONLY = reuse completed evidence aggregates; rerun canonical finalization only
  IDENTITY_RECONCILE_ONLY = reuse existing identity evidence; rerun only decisions + validation
  MT_FUSION   = per-cell mitochondrial ratio profiling; followed by MT_POPULATION unless --mt-ratios-only
  MT_POPULATION = population structure over MT profiles joined to reconciled-cell metadata

Conditions:
  Mode 1 (interindiv):   IND_BASE IND_LOO IND_AD IND_AD_LOO
                         IND_RF IND_RF_LOO IND_RA IND_RA_LOO
  Mode 1+SR legacy:      IND_AD_LOO_SR IND_RA_LOO_SR
  Legacy 2C comparators: LEGACY2C_CLASSIC LEGACY2C_TET_AWARE
  Mode 2 (species):      SP_BASE SP_LOO SP_AD SP_AD_LOO SP_RF SP_RF_LOO
  Mode 5 (species fixed): SP_FIXED_EMPTY
  Mode 3 (ind+WS):       IND_WS IND_WS_LOO IND_WS_AD_LOO IND_WS_RA_LOO
  CK factorial (16):     IND_CK IND_LOO_CK IND_AD_CK IND_AD_LOO_CK
                         IND_RF_CK IND_RF_LOO_CK IND_RA_CK IND_RA_LOO_CK
                         IND_WS_CK IND_WS_LOO_CK IND_WS_AD_CK IND_WS_AD_LOO_CK
                         IND_WS_RF_CK IND_WS_RF_LOO_CK IND_WS_RA_CK IND_WS_RA_LOO_CK
  CK gated validation:   IND_CK_RF_SX025_GATED_RFREE_PFIT
  Mode 4 (sp+WS):        SP_WS SP_WS_LOO SP_WS_AD_LOO
  Final CK additions:    IND_CK_RF_SX0_RFREE_PFIT
                         IND_CK_RF_SX0_GATED_RFREE_PFIT
                         IND_CK_RF_SX025_RFREE_PFIT
                         IND_CK_RF_SX050_RFREE_PFIT
                         IND_CK_RF_SX075_RFREE_PFIT
                         IND_CK_RF_SX1_RFREE_PFIT

Named condition sets:
  focused                 IND_CK_RF_SX0_GATED_RFREE_PFIT (default)
  all                     every historical and final condition
  real-library-comparison final CK roster plus essential historical controls
  ck-finalists            seven frozen CK finalists/components
  ck-rate                 SX0, SX0_GATED, SX025, SX025_GATED, SX1
  ck-profile              SX0, SX050, SX075, SX1
  ck-minimal              SX0, SX0_GATED, SX025_GATED
""")

    parser.add_argument("--libraries", nargs="+", default=None,
                        help="Library numbers or ranges (e.g., 1 2 8 or 1-40)")
    parser.add_argument("--conditions", nargs="+", default=None,
                        help=("Explicit condition abbreviations to run. Overrides "
                              "--condition-set."))
    parser.add_argument(
        "--contam-assignment-source",
        choices=["demux", "reconciled", "both", IDENTITY_AMBIENT_SELECTOR],
        default="demux",
        help=("Assignments used by CONTAM, GEX_AMBIENT, and AMBIENT_PLOTS. AMBIENT_VALIDATE "
              "automatically uses the applied E/F comparison. demux preserves "
              "the historical output location; reconciled writes below each "
              "condition's reconciled/ folder; both enables the legacy paired "
              "comparison; reconciliation-four-arm runs the controlled A-D "
              "identity/ambient-roster experiment with frozen assignments."))
    parser.add_argument(
        "--reconciliation-candidate-set",
        choices=["applied", "exploratory"], default="applied",
        help=("Candidate identities propagated into reconciliation-four-arm "
              "ambient rosters. applied is the production-safe default and "
              "admits only applied reconciliations/validated final receivers. "
              "exploratory additionally admits strong but unapplied alternative "
              "candidates and must not be interpreted as presence evidence. "
              "Default: applied."))
    parser.add_argument("--condition-set", choices=sorted(CONDITION_SETS),
                        default=DEFAULT_CONDITION_SET,
                        help=("Named condition roster. Default: focused. The full "
                              "historical registry remains available with all."))
    parser.add_argument("--stage", default=None,
                        help=("Run comma-separated stages; use AMBIENT_PLOTS for "
                              "non-ranking plots or AMBIENT_VALIDATE for the "
                              "standard applied-assignment control; use "
                              "AMBIENT_SWAP_TEST for an automatically discovered proposed-"
                              "assignment test; use GEX_AMBIENT for genome-wide "
                              "ambient-expression profiles; use CLEANUP_RESULTS alone to retain "
                              "only the newest generated results across all stages"))
    parser.add_argument("--submit", action="store_true",
                        help="Actually submit jobs to SLURM (default is dry run)")
    parser.add_argument("--skip-demux", action="store_true",
                        help="Skip DEMUX stage; validate demux outputs at sbatch runtime")
    parser.add_argument("--force", action="store_true",
                        help="Ignore existing outputs, rerun every selected stage")
    parser.add_argument("--regenerate-condf", action="store_true",
                        help="Rebuild all central CONDF files and include the CONDF stage without forcing downstream stages")
    parser.add_argument("--skip-het-vcf", action="store_true",
                        help="Disable --shared_het_vcf for DEMUX Pass 1; run only the individual and species VCF panels")
    parser.add_argument(
        "--individual-only-demux", action="store_true",
        help=("Run only the filtered-barcode interindividual DEMUX pass. "
              "Disables HET, species-panel work, the raw-barcode pass, pileup "
              "sidecars, and DEMUX CONDF symlinks while retaining the normal "
              "expected-lines restriction, selection audit, resources, daemon "
              "node placement, force/reuse behavior, and output-root layout."))
    parser.add_argument(
        "--demux-output-root",
        dest="demux_output_root", default=None,
        help=("Optional alternate root for DEMUX outputs. Inputs still come from "
              "the normal mapping_output tree; outputs are written under "
              "<root>/Tet_2025_Multiome-RNA_<N>/demux_nomito/. "
              "Use a different root for each A/B/C arm."))
    parser.add_argument("--condf-dir", default=CONDF_DIR,
                        help=("Optional shared CONDF generation/staging directory used by CONDF/DEMUX. "
                              "Current production EMPTY_DROPS/CONTAM read the self-contained "
                              "CONDF copies in each demux_nomito directory."))
    parser.add_argument("--panel-metadata", default=PANEL_METADATA,
                        help=("Panel metadata TSV matched to the selected main/species panels; "
                              "forwarded to both filtered and raw species counting and downstream stages."))
    parser.add_argument("--shared-main-segment", default=SHARED_VCF["interindiv_20M"],
                        help="Shared-memory segment for the interindividual panel")
    parser.add_argument("--shared-species-segment", default=SHARED_VCF["species_20M"],
                        help="Shared-memory segment for the species panel")
    parser.add_argument("--shared-het-segment", default=SHARED_VCF["interindiv_het_10M"],
                        help="Shared-memory segment for the HET panel; ignored with --skip-het-vcf")
    parser.add_argument("--daemon-nodes", default=DAEMON_NODELIST,
                        help=("Comma-separated SLURM node list for CONDF and DEMUX, "
                              "which consume node-local shared-memory panels. "
                              "Pass an empty string only when the requested segments "
                              "are provisioned on every eligible node."))
    parser.set_defaults(manage_vcf_daemons=True)
    parser.add_argument(
        "--manage-vcf-daemons", dest="manage_vcf_daemons",
        action="store_true",
        help=("Automatically launch one run-scoped VCF shared-memory daemon on "
              "each --daemon-nodes host and tear it down after all CONDF/DEMUX "
              "consumers settle (default)."))
    parser.add_argument(
        "--no-manage-vcf-daemons", dest="manage_vcf_daemons",
        action="store_false",
        help=("Use externally provisioned --shared-*-segment names and retain "
              "the historical manual daemon lifecycle."))
    parser.add_argument("--diagnose-only", action="store_true",
                        help="Print per-library status report, then exit")
    parser.add_argument("--skip-validation", action="store_true",
                        help="Skip pre-flight validation (for when upstream stages "
                             "will create missing files via dependency chains)")
    parser.add_argument("--with-ploidy-nn", action="store_true",
                        help="Also run optional PLOIDY_NN and TETRA_REFINE stages")
    parser.add_argument("--with-refine", action="store_true",
                        help="Also run optional TETRA_REFINE stage")
    parser.add_argument("--with-posthoc", action="store_true",
                        help="Also run tetra_score_calls and the swap-audit POSTHOC stage")
    parser.add_argument("--with-ambient-plots", action="store_true",
                        help="Also run selected-condition aggregate contamination plotting")
    parser.add_argument("--ambient-plot-root", default=AMBIENT_PLOT_ROOT,
                        help="Aggregate root for AMBIENT_PLOTS outputs")
    parser.add_argument("--ambient-plot-label", default=None,
                        help="Optional output label; default is deterministic from selected conditions")
    parser.add_argument("--ambient-reference-condition", default=None,
                        help="Selected reference condition for descriptive paired deltas")
    parser.add_argument(
        "--gex-cluster-source", choices=["auto", "h5ad", "manual"],
        default="auto",
        help=("Cluster source for GEX_AMBIENT. auto (default) creates stable "
              "RNA-only numbered Leiden groups; h5ad reuses an obs column; "
              "manual consumes --gex-clusters-template."))
    parser.add_argument(
        "--gex-clusters-template", default=None,
        help=("Manual-mode headerless barcode<tab>cluster file template supporting "
              "{lib}, {lib_num}, or {library}. Labels may be biological names or "
              "arbitrary group numbers."))
    parser.add_argument(
        "--gex-h5ad-cluster-column", default="leiden",
        help="H5AD obs column reused by --gex-cluster-source h5ad")
    parser.add_argument(
        "--gex-marker-genes", default=None,
        help=("Optional one-gene-per-line file forced into the automatic RNA "
              "feature union; it is not a cell-type assignment file."))
    parser.add_argument(
        "--gex-auto-data-features", type=int, default=5000,
        help=("Top raw-UMI Poisson-deviance genes used by automatic clustering "
              "before adding optional marker genes; 0 uses all detected genes."))
    parser.add_argument(
        "--gex-auto-min-genes", type=int, default=200,
        help=("Permissive detected-gene floor for cells that define the Leiden "
              "graph. Cells below it are retained by nearest-cluster transfer."))
    parser.add_argument(
        "--gex-auto-neighbors", type=int, default=40,
        help="Neighbors in the automatic RNA graph (default follows the final notebook workflow)")
    parser.add_argument(
        "--gex-auto-pca-components", type=int, default=20,
        help="Raw-MEX log-normalized truncated-SVD dimensions for automatic clustering")
    parser.add_argument(
        "--gex-auto-resolution-grid", default="0.15,0.25,0.35,0.5",
        help=("Comma-separated Leiden resolutions tested for stability and "
              "separation; no target cluster count is imposed."))
    parser.add_argument(
        "--gex-auto-stability-repeats", type=int, default=3,
        help="Leiden repeats per resolution used to measure ARI stability")
    parser.add_argument(
        "--gex-auto-exclude-genes-regex", default=r"^(MT-|RPS|RPL)",
        help=("Genes excluded only from automatic clustering. They remain in "
              "the ambient_rna_gex fit and genome-wide output."))
    parser.add_argument("--gex-auto-random-seed", type=int, default=1729,
                        help="Deterministic SVD/Leiden/transfer seed")
    parser.add_argument(
        "--gex-analysis-name", default="full_gene_rna_leiden_v1",
        help=("Versioned GEX_AMBIENT namespace. Change this when clusters, feature "
              "selection, or other analysis-defining inputs change."))
    parser.add_argument("--gex-ambient-root", default=GEX_AMBIENT_ROOT,
                        help="Aggregate root for GEX_AMBIENT comparison tables")
    parser.add_argument("--gex-feature-type", default="Gene Expression",
                        help="Feature type passed to tet_contam_estimate --feature_type")
    parser.add_argument(
        "--gex-skip-genes-regex", default="",
        help=("Regex passed to --skip_genes_regex. Default is empty so the "
              "genome-wide diagnostic does not silently discard MT genes."))
    parser.add_argument(
        "--gex-round-counts", action="store_true",
        help=("Stochastically round the corrected MEX to integers. Default keeps "
              "deterministic fractional counts via --noround for diagnostics."))
    parser.add_argument("--gex-min-cluster-cells", type=int, default=25,
                        help=("Minimum cells required in every nuisance cluster "
                              "after intersecting MEX and valid contam_rate barcodes"))
    parser.add_argument("--gex-cpus", type=int, default=8)
    parser.add_argument("--gex-memory-gb", type=int, default=128)
    parser.add_argument("--gex-top-n", type=int, default=200,
                        help="Number of ranked ambient genes retained per run in the summary")

    parser.add_argument("--with-hybrid", action="store_true",
                        help="Also run HYBRID post-hoc individual/species reconciliation; implies --with-posthoc")
    parser.add_argument("--with-mt-fusion", action="store_true",
                        help="Append the complete final mitochondrial analysis: MT_FUSION then automatic MT_POPULATION")
    parser.add_argument("--with-mt-population", action="store_true",
                        help="Compatibility alias for the complete MT_FUSION + MT_POPULATION analysis")
    parser.add_argument("--mt-ratios-only", action="store_true",
                        help="Run MT_FUSION per-cell ratios/profiles without the final population analysis")
    parser.add_argument("--mt-output-root", default=MT_FUSION_ROOT,
                        help="Root for per-library mitochondrial ratio/profile outputs and population results")
    parser.add_argument("--mt-vcf", default=MT_PANEL_FILES["vcf"],
                        help=("Dedicated mitochondrial panel BCF/VCF. Default: production NoMito "
                              "tet.vars.mt_fusion_ratio.bcf; override only for a deliberate panel revision."))
    parser.add_argument("--mt-site-manifest", default=MT_PANEL_FILES["site_manifest"],
                        help=("Library-aware mitochondrial site manifest. Default: matching production "
                              "NoMito tet.vars.mt_fusion_ratio.site_manifest.tsv."))
    parser.add_argument(
        "--identity-reconciliation-root", default=IDENTITY_RECONCILIATION_ROOT,
        help=("Root for identity metadata, candidates, evidence, decisions, reports, "
              "and validation. MT_FUSION/MT_POPULATION consume the decisions here."))
    parser.add_argument("--mt-library-id-template", default="lib{lib}",
                        help="Manifest library_id template; supports {lib}, {lib_num}, and {library}")
    parser.add_argument("--mt-bam-template", default=None,
                        help="Optional BAM path template; supports {lib}, {lib_num}, and {library}; default uses get_bam_path()")
    parser.add_argument(
        "--mt-assignment-source", choices=["reconciled", "demux", "refined", "auto"],
        default="reconciled",
        help=("Assignment source for mt_fusion_ratio. Default: reconciled single-cell assignments; "
              "auto is retained as an alias for reconciled and does not silently fall back."))
    parser.add_argument(
        "--mt-rna-ambient-condition", default=None,
        help=("Optional contamination condition whose per-cell RNA contam_rate "
              "is annotated in MT_FUSION and analyzed by MT_POPULATION."))
    parser.add_argument(
        "--mt-rna-ambient-assignment-source",
        choices=["demux", "reconciled", *IDENTITY_AMBIENT_ARM_ORDER],
        default="reconciled",
        help=("Assignment/roster variant of --mt-rna-ambient-condition. "
              "The four explicit arm keys let MT_FUSION consume the matched "
              "controlled per-barcode RNA covariate; reconciled_augmented is "
              "arm C. This does not change MT cell identities, which remain "
              "controlled independently by --mt-assignment-source."))
    parser.add_argument(
        "--mt-rna-ambient-max", type=float, default=None,
        help=("Optional RNA contamination cutoff. MT_FUSION records the flag; "
              "MT_POPULATION excludes missing or above-cutoff cells. Omit to "
              "annotate and report associations without filtering."))
    parser.add_argument("--mt-assay-mode", choices=["RNA", "ATAC", "GENERIC"], default="RNA")
    parser.add_argument("--mt-likelihood", choices=["beta_binomial", "binomial"], default="beta_binomial")
    parser.add_argument("--mt-ambient-profile-template", default=None,
                        help="Optional per-library mt ambient-profile template; explicit source overrides automatic MT non-cell barcode generation")
    parser.add_argument("--mt-empty-barcodes-template", default=None,
                        help="Optional per-library non-cell barcode template for mt ambient estimation; explicit source overrides automatic generation")
    parser.add_argument(
        "--mt-ambient-qc-max",
        type=float,
        default=None,
        help=("Maximum independently estimated MT ambient fraction allowed for parental-ratio inference; "
              "with no explicit ambient source, derive same-library chrM non-cell barcodes automatically")
    )
    parser.add_argument("--mt-ambient-none", action="store_true",
                        help=("Explicit no-ambient mode. This is already the default when neither an "
                              "mt ambient profile nor an empty-barcode file is supplied."))
    parser.add_argument("--mt-site-calibration", default=None,
                        help="Optional learned site-calibration TSV passed to mt_fusion_ratio")
    parser.add_argument("--mt-site-calibration-stratum-template", default=None,
                        help="Optional site-calibration stratum template; supports {lib}, {lib_num}, and {library}")
    parser.add_argument("--mt-rho-mode", choices=["free", "fixed", "low_information_fixed", "shrink"], default="free")
    parser.add_argument("--mt-pooled-rho", type=float, default=None)
    parser.add_argument("--mt-rho-reference", default=None)
    parser.add_argument("--mt-rho-low-information-molecules", type=int, default=50)
    parser.add_argument("--mt-rho-prior-strength", type=float, default=10.0)
    parser.add_argument(
        "--mt-site-influence-mode", choices=["full", "none"], default="none",
        help=("Per-site leave-one-out influence diagnostics for MT_FUSION. Default: none "
              "for production speed; site-count rows and schema are still written. Use full "
              "when LOO diagnostic values are specifically required."))
    parser.add_argument("--mt-mask-bed", default=None)
    parser.add_argument("--mt-atac-include-singletons", action="store_true")
    parser.add_argument("--mt-single-parent-epsilon", type=float, default=0.02)
    parser.add_argument("--mt-profile-grid-step", type=float, default=0.005)
    parser.add_argument("--mt-group-by", default="assay_mode,parent_pair_id,fusion_replicate_id")
    parser.add_argument("--mt-population-prefix", default=None,
                        help="Optional output prefix for final population files; default is <mt-output-root>/population/mt_population")
    parser.add_argument("--mt-meaningful-sd", type=float, default=0.05)
    parser.add_argument("--mt-min-cells", type=int, default=20)
    parser.add_argument("--mt-delta-bic", type=float, default=10.0)
    parser.add_argument("--mt-min-component-fraction", type=float, default=0.10)
    parser.add_argument("--mt-min-component-cells", type=int, default=10)
    parser.add_argument("--mt-min-component-separation", type=float, default=0.10)
    parser.add_argument("--mt-membership-threshold", type=float, default=0.80)
    parser.add_argument("--mt-min-confident-fraction", type=float, default=0.70)
    parser.add_argument("--with-identity-reconciliation", action="store_true",
                        help=("Append the RNA-centered IDENTITY_RECONCILIATION "
                              "composite stage with frozen nuclear and ambient evidence"))
    parser.add_argument(
        "--with-identity-score", action="store_true",
        help=("Append non-mutating IDENTITY_SCORE using existing validated "
              "reconciliation proposals; generic constrained/unconstrained "
              "runners are never scored"))
    parser.add_argument("--identity-reconcile-only", action="store_true",
                        help=("Run only identity reconciliation + validation using existing metadata, "
                              "candidates, nuclear/MT[/ATAC] scores, NN calls, audit QC, and doublet "
                              "context; does not rerun candidate generation or BAM scoring"))
    parser.add_argument("--identity-audit-root", default=IDENTITY_AUDIT_ROOT,
                        help="Authoritative NoMito swap-audit root consumed by identity reconciliation")
    parser.add_argument("--mt-panel", default=MT_PANEL_FILES["vcf"],
                        help="Mitochondrial identity panel BCF")
    parser.add_argument("--mt-haplotype-groups", default=MT_PANEL_FILES["haplotype_groups"],
                        help="MT haplotype-group sidecar")
    parser.add_argument("--mt-haplotype-pairwise", default=MT_PANEL_FILES["haplotype_pairwise"],
                        help="MT pairwise-resolution sidecar")
    parser.add_argument("--identity-evidence-mode", choices=["rna", "rna-atac"], default="rna",
                        help="rna never opens/validates ATAC inputs; rna-atac adds optional ATAC evidence")
    parser.add_argument("--identity-with-atac", action="store_true",
                        help="Alias for --identity-evidence-mode rna-atac")
    parser.add_argument("--atac-panel", default=None,
                        help="Explicit ATAC genotype panel for rna-atac scoring")
    parser.add_argument("--atac-het-panel", default=None,
                        help="Optional ATAC het-balance panel")
    parser.add_argument("--atac-barcode-map", default=None,
                        help="Optional ATAC->RNA barcode map path/template")
    parser.add_argument("--rna-barcode-whitelist", default=None,
                        help="RNA whitelist path/template used to construct an ATAC map")
    parser.add_argument("--atac-barcode-whitelist", default=None,
                        help="ATAC whitelist path/template used to construct an ATAC map")
    parser.add_argument("--identity-atac-root", default=ATAC_MAPPING_ROOT,
                        help="Root containing Tet_2025_Multiome-ATAC_N/atac.bam")
    parser.add_argument("--identity-atac-bam-template", default=None,
                        help="Optional explicit ATAC BAM template supporting {lib}")
    parser.add_argument("--identity-require-mt", action="store_true",
                        help="Make missing/failed mt identity evidence fatal")
    parser.add_argument("--identity-auto-apply", dest="identity_auto_apply", action="store_true", default=True,
                        help="Apply DECISIVE policy-supported identity changes (default)")
    parser.add_argument("--identity-no-auto-apply", dest="identity_auto_apply", action="store_false",
                        help="Review-only override: write reconciliation decisions without applying changes")
    parser.add_argument("--identity-candidate-top-k", type=int, default=5)
    parser.add_argument("--identity-max-candidates", type=int, default=12)
    parser.add_argument("--identity-site-folds", type=int, default=5)
    parser.add_argument(
        "--identity-probability-resamples", type=int, default=100,
        help="Deterministic bootstrap/downsample replicates per identity pair")
    parser.add_argument(
        "--identity-probability-seed", type=int, default=1729,
        help="Base seed for reproducible identity-score stability audits")
    parser.add_argument(
        "--identity-poor-fit-residual", type=float, default=0.30,
        help=("Absolute allele-fraction residual above which neither compared "
              "candidate is treated as fitting (default 0.30)"))
    parser.add_argument(
        "--identity-score-output-root", default=None,
        help=("Aggregate score output; default "
              "<identity-reconciliation-root>/"
              "reconciliation_swap_score_summary_v6_1 (or the separate "
              "v6_1_frozen_reuse directory for aggregate-only)"))
    parser.add_argument(
        "--identity-score-ambient-root", default=AMBIENT_PLOT_ROOT,
        help=("Optional completed controlled ambient-analysis root used only "
              "for separate c/multimodality annotations; empty disables"))
    parser.add_argument(
        "--identity-score-ambient-condition", default=None,
        help=("Optional controlled ambient condition; default prefers the "
              "production CK gated condition when present"))
    parser.add_argument(
        "--identity-candidate-axis-event-id", default=None,
        help=("Optional exact frozen supported event ID; omit with proposal for "
              "routine automatic enumeration, or supply both for one-library diagnostics"))
    parser.add_argument(
        "--identity-candidate-axis-proposal", default=None,
        help=("Optional exact fixed candidate-B identity paired with the event-ID "
              "diagnostic override; never expands to a roster-wide search"))
    parser.add_argument(
        "--identity-candidate-axis-output-root", default=None,
        help=("Absolute isolated output root containing score_pairs, scorer, "
              "aggregate, and logs; required for standalone candidate-axis "
              "runs and optional as a compatible frozen root for routine reconciliation"))
    parser.add_argument(
        "--identity-candidate-axis-v6-3-root", default=None,
        help=("Optional absolute frozen V6.3 aggregate root containing cell and "
              "review-only score tables for verbatim baseline priority"))
    parser.add_argument(
        "--identity-review-input", default=None,
        help=("Optional final review TSV with library, barcode_or_event_id, "
              "review_disposition, and rationale"))
    parser.add_argument("--lib19-full-test", action="store_true",
                        help=("Preset for lib19 stages/libraries only; condition selection "
                              "still uses the selected --condition-set unless --conditions is supplied"))
    parser.add_argument("--audit-root", default=AUDIT_ROOT,
                        help="Root for POSTHOC swap-audit/QC outputs")
    parser.add_argument("--hybrid-root", default=HYBRID_ROOT,
                        help="Root for HYBRID post-hoc reconciliation outputs")
    parser.add_argument("--refined-assignments-root", default=TETRA_REFINE_ROOT,
                        help="Optional tetra_refine output root shaped as <root>/lib<N>/lib<N>.refined_assignments")
    parser.add_argument("--ploidy-calls-root", default=PLOIDY_CALLS_ROOT,
                        help="Root containing lib<N>.ploidy_calls_nn.tsv and all_libs.ploidy_calls_nn.tsv")
    parser.add_argument("--ploidy-nn-h5ad", default=PLOIDY_NN_H5AD,
                        help="h5ad input for optional PLOIDY_NN stage")
    parser.add_argument("--ploidy-nn-weights", default=PLOIDY_NN_WEIGHTS,
                        help="ploidy NN weights .pt file; sibling _scaler.npz must exist")
    parser.add_argument("--ploidy-nn-module", default=PLOIDY_NN_MODULE,
                        help="Environment module to load for PLOIDY_NN; empty string disables module loading")
    parser.add_argument("--ploidy-nn-qc-only", action="store_true",
                        help="Pass --qc_only to run_ploidy_nn_inference.py")
    parser.add_argument("--refine-contam-condition", default=REFINE_CONTAM_CONDITION,
                        help="Contamination condition whose contam_rate is passed to tetra_refine; empty string disables")
    parser.add_argument("--refine-external-ploidy-min-prob", type=float, default=REFINE_EXTERNAL_PLOIDY_MIN_PROB,
                        help="Minimum external ploidy confidence required for tetra_refine to relabel A -> A+A")
    parser.add_argument("--require-refine-external-ploidy", action="store_true",
                        help="Make TETRA_REFINE fail if lib<N>.ploidy_calls_nn.tsv is absent")
    parser.add_argument("--require-refine-contam-rate", action="store_true",
                        help="Make TETRA_REFINE fail if the selected contam_rate is absent")
    parser.add_argument("--refine-scoring-only", action="store_true",
                        help="Pass --scoring_only to run_tetra_refine_for_library.py")
    parser.add_argument("--expected-pool-metadata", default=EXPECTED_POOL_METADATA,
                        help="Expected pool metadata TSV, e.g. pool_combinations.tsv")
    parser.add_argument("--allowed-identities", default=None,
                        help="Optional explicit constrained identity universe for swap audit")
    parser.add_argument("--hybrid-individual-condition", default="IND_RF",
                        help="Individual-native condition used by HYBRID reconciliation")
    parser.add_argument("--hybrid-species-condition", default="SP_RF",
                        help="Species-native condition used by HYBRID reconciliation")
    parser.add_argument("--hybrid-fixed-species-condition", default="SP_FIXED_EMPTY",
                        help="Fixed-empty species comparator used by HYBRID reconciliation")

    gex_add_internal_worker_argument(parser)
    ambient_add_internal_worker_argument(parser)
    parser.add_argument(
        "--_identity-final-evidence-worker", default=None,
        help=argparse.SUPPRESS)

    return parser.parse_args()


def main():
    args = parse_args()
    if gex_maybe_run_internal_worker(args):
        return
    if ambient_maybe_run_internal_worker(args):
        return
    if args._identity_final_evidence_worker is not None:
        identity_final_evidence_worker(
            args._identity_final_evidence_worker)
        return
    run(args)


if __name__ == "__main__":
    main()
