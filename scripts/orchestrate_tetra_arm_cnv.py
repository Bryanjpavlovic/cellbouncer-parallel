#!/usr/bin/env python3
"""Submit the downstream Tetraploid chromosome-arm ASE/CNV workflow.

This orchestrator deliberately starts at the canonical outputs of
``orchestrate_tetraploid.py``.  It does not modify or extend that upstream
orchestrator.  The fixed job graph is::

    REFERENCE -----------------> ASE[] ---\
    LEDGER -> PREPARE[] -------> ASE[] ----+-> CALL -> REPORT
                         EXPRESSION[] -----/

Planning and script generation do not submit jobs.  ``--submit``
submits the selected stages with SLURM ``afterok`` dependencies.  Array task
manifests are ordinary, deterministic TSV files so every command and input can
be inspected before submission.
"""

from __future__ import annotations

import argparse
import csv
import glob
import gzip
import hashlib
import itertools
import json
import math
import os
import re
import shlex
import sqlite3
import stat
import statistics
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from typing import Iterable, Mapping, Sequence


RELEASE = "2.4.0"
STAGES = ("REFERENCE", "LEDGER", "PREPARE", "ASE", "EXPRESSION", "CALL", "REPORT")
STAGE_PARENTS = {
    "REFERENCE": (),
    "LEDGER": (),
    "PREPARE": ("LEDGER",),
    "ASE": ("REFERENCE", "PREPARE"),
    "EXPRESSION": ("PREPARE",),
    "CALL": ("ASE", "EXPRESSION"),
    "REPORT": ("CALL",),
}

# Production paths inherited from orchestrate_tetraploid.py.  Production keeps
# its historical combined mapping/analysis tree.  Alternate remaps use a
# separate mapping-input tree and upstream-analysis tree, matching
# orchestrate_tetraploid.py's isolation contract.
PROJECT_ROOT = "/mnt/beegfs/tetmultiome_rna_mapped"
MAPPING_ROOT = os.path.join(PROJECT_ROOT, "mapping_output")
DEFAULT_PANEL_METADATA = os.path.join(
    PROJECT_ROOT, "Misc_Metadata", "panel_metadata.tsv")
DEFAULT_DEMUX_POOL_WORKBOOK = os.path.join(
    PROJECT_ROOT, "Misc_Metadata", "Library_conversions.xlsx")
DEFAULT_EXPECTED_POOL_METADATA = os.path.join(
    PROJECT_ROOT, "Misc_Metadata", "pool_combinations.tsv")
DEFAULT_IDENTITY_WORKBOOK = os.path.join(
    MAPPING_ROOT, "Library_conversions.xlsx")
DEFAULT_CONDITION = "IND_CK_RF_SX0_GATED_RFREE_PFIT"
DEFAULT_GEX_AMBIENT_ANALYSIS = "full_gene_rna_leiden_v1"
LIBRARY_PREFIX = "Tet_2025_Multiome-RNA_"
DEMUX_SUBDIR = "demux_nomito"
NOMITO_PANEL_ROOT = (
    "/mnt/beegfs/home/b/vcfdownsample/"
    "Downsample_ATAC_Species_poolInformer/NoMito")
DEFAULT_INTERINDIVIDUAL_PANEL = os.path.join(
    NOMITO_PANEL_ROOT, "tet.vars.downsampled_20M.bcf")
DEFAULT_HET_PANEL = os.path.join(NOMITO_PANEL_ROOT, "tet.vars.het_10M.bcf")
DEFAULT_SPECIES_PANEL = os.path.join(NOMITO_PANEL_ROOT, "tet.vars.species_20M.bcf")
DEFAULT_PLOIDY_NN_WEIGHTS = os.path.join(
    PROJECT_ROOT, "ploidy_classifier", "retrain_nomito_20260814", "model",
    "ploidy_nn_weights.pt")

CELLBOUNCER_ROOT = "/nvme/software/packages/cellbouncer/dev"
DEPLOYED_BIN = os.path.join(CELLBOUNCER_ROOT, "bin")
DEPLOYED_SCRIPTS = DEPLOYED_BIN
FUSEBOX_DATA = "/nvme/software/packages/fusebox/latest/data"

DEFAULT_GENE_ARMS = os.path.join(
    FUSEBOX_DATA, "hg38_gene_arms_noXCI.txt")
DEFAULT_SOURCE_ARMS = os.path.join(FUSEBOX_DATA, "hg38_arms.bed")
DEFAULT_ARM_BUILDER_SCRIPT = os.path.join(
    DEPLOYED_SCRIPTS, "tetra_arm_gene_synteny.py")
DEFAULT_HAL_REFERENCE_SCRIPT = os.path.join(
    DEPLOYED_SCRIPTS, "tetra_arm_hal_liftover.py")
DEFAULT_PREPARE_SCRIPT = os.path.join(DEPLOYED_SCRIPTS, "tetra_arm_prepare.py")
DEFAULT_EXPRESSION_SCRIPT = os.path.join(
    DEPLOYED_SCRIPTS, "tetra_arm_expression.py")
DEFAULT_CALL_SCRIPT = os.path.join(DEPLOYED_SCRIPTS, "tetra_arm_call.py")
DEFAULT_REPORT_SCRIPT = os.path.join(DEPLOYED_SCRIPTS, "tetra_arm_report.py")
DEFAULT_ASE_BINARY = os.path.join(DEPLOYED_BIN, "tetra_arm_ase")
DEFAULT_PLOIDY_NN_HELPER = os.path.join(
    DEPLOYED_SCRIPTS, "run_ploidy_nn_inference.py")

BASE_PYTHON_MODULES = ("miniforge/3",)
SCIENTIFIC_PYTHON_MODULES = ("miniforge/3", "genomics-base/latest")
ASE_MODULES = ("htslib/1.20", "cellbouncer/dev")
DEFAULT_PARTITION = "compute"
DEFAULT_TIME = "7-00:00:00"

PREPARE_HEADER = (
    "task_index", "library", "ledger", "final_assignments",
    "ambient_standard_prefix", "ambient_arm_a_prefix",
    "ambient_arm_c_prefix", "panel_metadata", "ploidy_nn", "cell_groups",
    "output_dir",
)
ASE_HEADER = (
    "task_index", "library", "samples", "pileup_sites",
    "pileup_molecules", "pileup_observations", "cell_manifest",
    "ambient_sources", "arms_bed", "output", "qc",
)
ASE_V2_HEADER = (
    "library", "barcode", "donor_a", "donor_b", "donor_pair", "arm",
    "chromosome", "arm_start", "arm_end", "ambient_c", "ambient_c_se",
    "n_sites", "n_molecules", "n_informative_molecules",
    "n_molecules_ref", "n_molecules_alt", "n_molecules_mixed",
    "a_ref", "b_ref", "a_alt", "b_alt", "a_mixed", "b_mixed",
    "n_ambiguous", "soft_a", "soft_b", "soft_a_ref", "soft_b_ref",
    "soft_a_alt", "soft_b_alt", "soft_a_mixed", "soft_b_mixed",
    "soft_a_sumsq_ref", "soft_a_sumsq_alt", "soft_a_sumsq_mixed",
    "effective_a_ref", "effective_b_ref", "effective_a_alt", "effective_b_alt",
    "effective_a_mixed", "effective_b_mixed",
    "effective_weight_ref", "effective_weight_alt", "effective_weight_mixed",
    "ambient_a_ref", "ambient_a_alt", "ambient_a_mixed",
    "ambient_genotyped_mass", "qname_fallback_fraction",
    "mean_sites_per_molecule", "evidence_basis", "model_eligible",
    "evidence_status", "schema_version",
)

# Exact schemas emitted by the attached production POSTHOC implementation.
# These are intentionally strict: the upstream summarizer otherwise accepts a
# truncated intersection of barcode sets and silently falls back to reduced
# feature modes when a sidecar is merely present but empty.
POSTHOC_CALL_QC_HEADER = (
    "barcode", "libname", "assignment", "assignment_type", "ploidy_status",
    "llr_vs_runner_up", "runnerup_comparison_state", "margin_softmax_score",
    "total_depth", "n_informative_bins", "n_informative_depth", "n_close",
    "depth_normalized_llr_vs_runner_up", "dosage_concordance",
    "dosage_runnerup_identity", "dosage_runnerup_comparison_state",
    "runnerup_dosage_concordance", "dosage_gap_constrained",
    "residual_mismatch", "expected_species_set", "species_support_expected",
    "species_conflict_flag", "species_relation",
    "species_missing_expected_component", "species_has_unexpected_component",
    "species_disjoint_wrong_species", "species_best_identity",
    "species_best_support", "species_gap", "call_qc_flags", "warnings",
)
POSTHOC_SPECIES_QC_HEADER = (
    "library", "n_cells_with_species_evidence",
    "median_species_support_expected", "frac_cells_species_conflict",
    "frac_cells_species_exact_match", "frac_cells_species_component_missing",
    "frac_cells_species_unexpected_extra",
    "frac_cells_species_partial_overlap_extra",
    "frac_cells_species_disjoint_wrong",
    "frac_cells_species_unexpected_or_disjoint", "expected_species_set",
    "observed_species_evidence", "warnings",
)
POSTHOC_RUNNER_UP_HEADER = (
    "barcode", "rank", "identity", "llr_vs_winner", "min_margin",
    "comparison_state", "winner_candidate", "selection_resolved",
    "schema_version",
)
POSTHOC_REPORT_HEADER = (
    "library", "preflight_status", "audit_verdict", "audit_flags",
    "feature_mode", "has_call_qc", "has_species_qc", "has_atac_qc",
    "has_refined_assignments", "has_calibrated_thresholds",
    "missing_optional_features", "warnings", "thresholds_source", "n_cells",
    "expected_identity", "audit_best_identity_unconstrained",
    "audit_best_identity_constrained", "audit_best_fraction",
    "expected_rank_median", "expected_rank_p90",
    "delta_ll_best_vs_expected_median", "median_llr_vs_runner_up",
    "frac_cells_runner_up_comparison_complete", "frac_cells_high_n_close",
    "frac_cells_component_overlap_2", "frac_cells_component_overlap_1",
    "frac_cells_component_overlap_0", "median_dosage_concordance_assigned",
    "median_dosage_concordance_unconstrained", "median_dosage_gap_constrained",
    "median_dosage_gap_unconstrained", "frac_cells_neg_gap_constrained",
    "frac_cells_neg_gap_unconstrained", "expected_vs_observed_l1",
    "expected_vs_observed_cosine", "jensen_shannon_distance",
    "unexpected_identity_fraction", "missing_expected_identity_fraction",
    "frac_cells_unconstrained_unexpected_identity",
    "frac_cells_unconstrained_unexpected_identity_supported",
    "frac_cells_unconstrained_unexpected_component",
    "frac_cells_unconstrained_unexpected_component_supported",
    "top_unexpected_component_supported",
    "top_unexpected_component_fraction_supported", "refined_n_cells",
    "refined_n_total_available", "refined_n_overlap_with_audit",
    "refined_overlap_fraction", "refined_n_changed", "refined_n_singlet",
    "refined_n_heterotypic", "refined_n_homotypic", "refined_n_plus_identity",
    "refined_n_biological_fusion", "refined_n_doublet",
    "refined_n_droplet_flagged", "refined_top_identity",
    "refined_top_fraction", "refined_biological_unexpected_fraction",
    "refined_biological_missing_expected_fraction", "species_expected",
    "species_audit_best", "median_species_support_expected",
    "frac_cells_species_conflict", "frac_cells_species_exact_match",
    "frac_cells_species_component_missing", "frac_cells_species_unexpected_extra",
    "frac_cells_species_partial_overlap_extra", "frac_cells_species_disjoint_wrong",
    "frac_cells_species_unexpected_or_disjoint", "unresolved_homotypic_ratio",
    "atac_qc_mode", "median_atac_dosage_concordance", "atac_best_identity",
    "rna_atac_same_identity_fraction", "rna_atac_same_species_fraction",
    "frac_cells_rna_atac_discordant",
)
POSTHOC_SCORE_HEADER = (
    "barcode", "libname", "raw_demux_assignment", "snp_resolvable_assignment",
    "refined_biological_assignment", "refined_assignment_source",
    "refined_assignment_confidence", "droplet_doublet_flag", "quad_pattern_score",
    "original_assignment", "audit_unconstrained_best", "audit_constrained_best",
    "expected_components", "audit_best_components", "shared_components",
    "missing_expected_components", "unexpected_components",
    "unexpected_pool_components", "component_overlap",
    "expected_rank_unconstrained", "expected_rank_constrained",
    "delta_ll_best_vs_expected", "dosage_concordance_assigned",
    "dosage_concordance_unconstrained", "dosage_gap_constrained",
    "dosage_gap_unconstrained", "species_support_expected", "species_relation",
    "species_missing_expected_component", "species_has_unexpected_component",
    "species_disjoint_wrong_species", "combined_qc_flags", "swap_cell_verdict",
    "feature_mode", "warnings",
)
EXPRESSION_HEADER = (
    "task_index", "library", "barcodes", "features", "matrix",
    "cell_manifest", "gene_arms", "output_dir",
)

CALL_INTEGER_DEFAULTS = (
    ("min_call_sites", 3),
    ("min_calibration_sites", 2),
    ("min_fallback_calibration_cells", 12),
    ("calibration_crossfit_folds", 5),
    ("min_cell_baseline_arms", 4),
    ("min_empirical_null_cells", 20),
    ("absolute_min_null_cells", 5),
    ("min_expression_genes", 3),
    ("min_uid_cells_per_arm", 1),
    ("min_pair_recurrence_cells", 2),
    ("max_pair_test_cells", 500),
    ("max_pair_cells_per_uid", 100),
)
CALL_FLOAT_DEFAULTS = (
    ("min_call_effective_weight", 8.0),
    ("min_calibration_effective_weight", 4.0),
    ("min_orientation_calibration_effective_weight", 20.0),
    ("calibration_shrinkage_cells", 20.0),
    ("calibration_max_cell_weight", 50.0),
    ("calibration_huber_z", 2.5),
    ("calibration_min_robust_weight", 0.25),
    ("calibration_expression_mad_cutoff", 4.0),
    ("cell_baseline_shrinkage_arms", 8.0),
    ("default_rho", 0.02),
    ("min_rho", 0.001),
    ("max_rho", 0.25),
    ("event_prior", 0.04),
    ("min_event_posterior", 0.90),
    ("min_balanced_posterior", 0.80),
    ("max_event_q", 0.05),
    ("empirical_depth_fold", 2.0),
    ("empirical_ambient_window", 0.10),
    ("min_expression_counts", 100.0),
    ("min_expression_sigma", 0.15),
    ("max_expression_sigma", 2.0),
    ("expression_weight", 0.50),
    ("max_expression_log_bf", 4.0),
    ("site_fallback_likelihood_weight", 0.50),
    ("max_qname_fallback_fraction", 0.50),
    ("min_ambient_genotyped_mass", 0.50),
    ("min_uid_arm_posterior", 0.80),
    ("min_uid_state_concordance", 0.67),
    ("max_aggregate_cell_log_bf", 8.0),
    ("max_whole_chromosome_q", 0.05),
    ("min_pair_state_concordance", 0.60),
    ("min_pair_cell_posterior", 0.50),
    ("min_pair_arm_posterior", 0.90),
    ("max_pair_arm_q", 0.05),
)


@dataclass(frozen=True)
class LibraryPaths:
    library: int
    mapping_bam: str
    mapping_bam_index: str
    final_assignments: str
    demux_prefix: str
    samples: str
    pileup_sites: str
    pileup_molecules: str
    pileup_observations: str
    expression_barcodes: str
    expression_features: str
    expression_matrix: str
    ambient_standard_prefix: str
    ambient_arm_a_prefix: str
    ambient_arm_c_prefix: str
    ploidy_nn: str
    cell_groups: str
    cell_groups_source: str
    split_ledger: str
    cell_manifest: str
    ambient_sources: str
    prepare_qc: str
    prepare_contract: str
    ase: str
    ase_qc: str
    expression: str
    expression_qc: str
    expression_contract: str


@dataclass(frozen=True)
class RunPaths:
    root: str
    logs: str
    scripts: str
    manifests: str
    reference: str
    ledger: str
    prepare: str
    ase: str
    expression: str
    call: str
    report: str
    call_prefix: str

    @property
    def generated_arms_bed(self) -> str:
        return os.path.join(self.reference, "ancestral_arms.bed")

    @property
    def reference_qc(self) -> str:
        return os.path.join(self.reference, "reference_qc.tsv")

    @property
    def reference_contract(self) -> str:
        return os.path.join(self.reference, "reference_contract.json")

    @property
    def call_outputs(self) -> tuple[str, ...]:
        return tuple(self.call_prefix + suffix for suffix in (
            ".arm_calls.tsv.gz",
            ".uid_chromosome_flags.tsv.gz",
            ".calibration.tsv.gz",
            ".qc.tsv",
            ".contract.json",
            ".donor_pair_arm_summary.tsv.gz",
        ))

    @property
    def report_outputs(self) -> tuple[str, ...]:
        return (
            os.path.join(self.report, "tetra_arm_cnv_summary.tsv"),
            os.path.join(self.report, "tetra_arm_cnv_summary.json"),
            os.path.join(self.report, "tetra_arm_cnv_report.html"),
        )

    @property
    def donor_pair_summary(self) -> str:
        return self.call_prefix + ".donor_pair_arm_summary.tsv.gz"

    @property
    def input_provenance(self) -> str:
        return os.path.join(self.root, "input_provenance.json")


def parse_libraries(values: Sequence[str]) -> list[int]:
    result: set[int] = set()
    for raw in values:
        for token in str(raw).split(","):
            token = token.strip()
            if not token:
                continue
            match = re.fullmatch(r"(?:lib)?(\d+)(?:-(\d+))?", token, re.I)
            if not match:
                raise ValueError(f"invalid library selection: {token!r}")
            first = int(match.group(1))
            last = int(match.group(2) or first)
            if last < first:
                raise ValueError(f"descending library range is not allowed: {token}")
            if first < 1 or last > 40:
                raise ValueError(f"library selection is outside 1-40: {token}")
            result.update(range(first, last + 1))
    if not result:
        raise ValueError("at least one library must be selected")
    return sorted(result)


def parse_stages(values: Sequence[str] | None) -> tuple[str, ...]:
    if not values:
        return STAGES
    selected: set[str] = set()
    for raw in values:
        for token in str(raw).split(","):
            stage = token.strip().upper()
            if not stage:
                continue
            if stage == "ALL":
                selected.update(STAGES)
            elif stage in STAGES:
                selected.add(stage)
            else:
                raise ValueError(
                    f"unknown stage {token!r}; choose from {', '.join(STAGES)}")
    if not selected:
        raise ValueError("at least one stage must be selected")
    return tuple(stage for stage in STAGES if stage in selected)


def absolute(value: str) -> str:
    return os.path.abspath(os.path.expanduser(value))


def validate_text(value: str, label: str) -> str:
    if "\n" in value or "\r" in value or "\t" in value:
        raise ValueError(f"{label} contains a tab or newline")
    return value


def call_parameters(args: argparse.Namespace) -> dict[str, int | float]:
    names = [name for name, _default in (
        *CALL_INTEGER_DEFAULTS, *CALL_FLOAT_DEFAULTS)]
    return {name: getattr(args, name) for name in names}


def expand_template(template: str, library: int, label: str) -> str:
    if not template:
        return ""
    try:
        rendered = template.format(
            lib=library, lib_num=library, library=f"lib{library}")
    except (KeyError, IndexError, ValueError) as exc:
        raise ValueError(
            f"{label} supports only {{lib}}, {{lib_num}}, and {{library}}: {exc}") \
            from exc
    return absolute(validate_text(rendered, label))


def regular_nonempty(path: str) -> bool:
    try:
        return os.path.isfile(path) and os.path.getsize(path) > 0
    except OSError:
        return False


def lexical_path_contains(parent: str, child: str) -> bool:
    """Check the declared namespace without resolving scientific symlinks."""
    parent_path = absolute(parent)
    child_path = absolute(child)
    try:
        return os.path.commonpath((parent_path, child_path)) == parent_path
    except ValueError:
        return False


def sha256_file(path: str) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def provenance_file_record(path: str, label: str,
                           include_sha256: bool = False) -> dict[str, object]:
    path = absolute(path)
    if not regular_nonempty(path):
        raise ValueError(f"{label} is missing or empty: {path}")
    stat_result = os.stat(path)
    record: dict[str, object] = {
        "label": label,
        "path": path,
        "realpath": os.path.realpath(path),
        "size": stat_result.st_size,
        "mtime_ns": stat_result.st_mtime_ns,
        "device": stat_result.st_dev,
        "inode": stat_result.st_ino,
    }
    if include_sha256:
        record["sha256"] = sha256_file(path)
    return record


def record_matches_path(record: object, path: str,
                        require_sha256: bool = False) -> tuple[bool, str]:
    """Compare a producer-recorded file identity with the current file.

    ``st_dev`` and ``st_ino`` are retained in records for local diagnostics but
    are deliberately not used as cross-host identity fields.  Parallel-file-
    system clients may synthesize different values for the same path on the
    submit and compute nodes.  Path, resolved path, size, mtime, and an optional
    digest are portable across those clients.
    """
    if not isinstance(record, dict):
        return False, "record is not an object"
    path = absolute(path)
    try:
        current = os.stat(path)
    except OSError as exc:
        return False, f"cannot stat current file: {exc}"
    expected = {
        "path": path,
        "realpath": os.path.realpath(path),
        "size": current.st_size,
        "mtime_ns": current.st_mtime_ns,
    }
    differences = [
        key for key, value in expected.items()
        if record.get(key) != value
    ]
    if differences:
        return False, "mismatched field(s): " + ",".join(differences)
    recorded_sha = record.get("sha256")
    if require_sha256 and (
            not isinstance(recorded_sha, str) or
            re.fullmatch(r"[0-9a-f]{64}", recorded_sha) is None):
        return False, "required sha256 is absent or malformed"
    if recorded_sha and recorded_sha != sha256_file(path):
        return False, "sha256 mismatch"
    return True, "PASS"


def active_shell_payload(path: str) -> str:
    """Return non-comment, non-echo shell text for provenance inspection."""
    active = []
    with open(path, "r", encoding="utf-8") as handle:
        for raw in handle:
            stripped = raw.strip()
            if not stripped or stripped.startswith("#"):
                continue
            candidate = stripped[:-1].rstrip() if stripped.endswith("\\") else stripped
            try:
                tokens = shlex.split(candidate, comments=True, posix=True)
            except ValueError:
                continue
            if not tokens:
                continue
            # Generated jobs do not place scientific commands after an
            # unconditional top-level exit.  Do not let unreachable text serve
            # as provenance evidence; indented exits inside normal shell gates
            # remain valid and do not terminate this static scan.
            if raw == raw.lstrip() and tokens[0] in {"exit", "return"}:
                break
            if tokens and tokens[0] in {"echo", "printf"}:
                continue
            active.append(" ".join(tokens))
    return "\n".join(active)


def shell_command_argvs(path: str, executable: str) -> list[list[str]]:
    """Parse one-line generated commands whose executable basename matches."""
    commands: list[list[str]] = []
    with open(path, "r", encoding="utf-8") as handle:
        for raw in handle:
            stripped = raw.strip()
            if not stripped or stripped.startswith("#"):
                continue
            try:
                tokens = shlex.split(stripped, comments=True, posix=True)
            except ValueError:
                continue
            if not tokens:
                continue
            if raw == raw.lstrip() and tokens[0] in {"exit", "return"}:
                break
            first = os.path.basename(tokens[0])
            if first == executable:
                commands.append(tokens)
            elif (first.startswith("python") and len(tokens) > 1 and
                  os.path.basename(tokens[1]) == executable):
                commands.append(tokens[1:])
    return commands


def logical_program_argvs(path: str, executable: str) -> list[list[str]]:
    """Return every logical invocation, including control-separated commands."""
    commands: list[list[str]] = []
    for tokens in logical_shell_argvs(path):
        first = os.path.basename(tokens[0])
        if first == executable:
            commands.append(tokens)
        elif (first.startswith("python") and len(tokens) > 1 and
              os.path.basename(tokens[1]) == executable):
            commands.append(tokens[1:])
    return commands


def split_shell_control(text: str) -> list[str]:
    """Split shell control operators outside quotes for provenance parsing."""
    result: list[str] = []
    start = 0
    quote = ""
    escaped = False
    index = 0
    while index < len(text):
        character = text[index]
        if escaped:
            escaped = False
            index += 1
            continue
        if quote == "'":
            if character == "'":
                quote = ""
            index += 1
            continue
        if quote == '"':
            if character == "\\":
                escaped = True
            elif character == '"':
                quote = ""
            index += 1
            continue
        if character in {"'", '"'}:
            quote = character
            index += 1
            continue
        if character == "\\":
            escaped = True
            index += 1
            continue
        if character in ";&|":
            candidate = text[start:index].strip()
            if candidate:
                result.append(candidate)
            while index < len(text) and text[index] in ";&|":
                index += 1
            start = index
            continue
        index += 1
    candidate = text[start:].strip()
    if candidate:
        result.append(candidate)
    return result


def logical_shell_argvs(path: str) -> list[list[str]]:
    """Parse shell commands across continuations and multiline quotations."""
    commands: list[list[str]] = []
    pending = ""
    with open(path, "r", encoding="utf-8") as handle:
        for raw in handle:
            stripped = raw.strip()
            if not pending and (not stripped or stripped.startswith("#")):
                continue
            pending += raw
            candidate = pending.replace("\\\n", " ")
            try:
                shlex.split(candidate, comments=True, posix=True)
            except ValueError as exc:
                if ("No closing quotation" in str(exc) or
                        "No escaped character" in str(exc)):
                    continue
                raise ValueError(
                    f"cannot parse generated shell command in {path}: {exc}") \
                    from exc
            if raw.rstrip().endswith("\\"):
                continue
            for command in split_shell_control(candidate):
                try:
                    tokens = shlex.split(command, comments=True, posix=True)
                except ValueError as exc:
                    raise ValueError(
                        f"cannot parse generated shell command in {path}: {exc}") \
                        from exc
                if tokens:
                    commands.append(tokens)
            pending = ""
    if pending:
        raise ValueError(f"unterminated shell command in: {path}")
    return commands


def vcf_holder_contract(path: str) -> dict[str, str]:
    """Extract the one foreground loader invocation from a generated holder."""
    logical = logical_shell_argvs(path)
    assignments: dict[str, str] = {}
    for argv in logical:
        if len(argv) != 1:
            continue
        match = re.fullmatch(r"([A-Z][A-Z0-9_]*)=(.*)", argv[0])
        if match:
            assignments[match.group(1)] = match.group(2)
    commands: list[list[str]] = []
    for argv in logical:
        executable = argv[0]
        if executable.startswith("$"):
            executable = assignments.get(executable[1:], "")
        if (os.path.basename(executable) == "vcf_loader_daemon" and
                "--foreground" in argv and "--destroy" not in argv):
            commands.append(argv)
    if len(commands) != 1:
        raise ValueError(
            "holder does not contain exactly one foreground "
            f"vcf_loader_daemon command: {path}")
    argv = commands[0]
    daemon = argv[0]
    if daemon.startswith("$"):
        daemon = assignments.get(daemon[1:], "")
    base = assignments.get("BASE", "")
    ready = assignments.get("READY", "")
    if (not base or not ready or argv_value(argv, "--name") != "$BASE" or
            argv_value(argv, "--qual") != "50" or
            argv.count("--foreground") != 1 or
            argv.count("--ready-file") != 1 or
            argv_value(argv, "--ready-file") != "$READY"):
        raise ValueError(f"holder loader contract is malformed: {path}")
    required = {
        "daemon": daemon,
        "base": base,
        "vcf": argv_value(argv, "--vcf"),
        "het_vcf": argv_value(argv, "--het_vcf"),
        "species_vcf": argv_value(argv, "--species_vcf"),
        "bam": argv_value(argv, "--bam"),
    }
    if not all(isinstance(value, str) and value for value in required.values()):
        raise ValueError(f"holder loader contract lacks a required value: {path}")
    return required


def sbatch_nodelist(path: str) -> set[str]:
    """Read the simple comma-separated nodelist used by this orchestrator."""
    with open(path, "r", encoding="utf-8") as handle:
        for raw in handle:
            match = re.match(
                r"^\s*#SBATCH\s+--nodelist(?:=|\s+)(\S+)\s*$", raw)
            if match:
                return {value for value in match.group(1).split(",") if value}
    return set()


def argv_value(argv: Sequence[str], option: str) -> str | None:
    """Return a unique option value, rejecting absent or repeated options."""
    positions = [index for index, value in enumerate(argv) if value == option]
    if len(positions) != 1 or positions[0] + 1 >= len(argv):
        return None
    return argv[positions[0] + 1]


def canonical_barcode(value: object) -> str:
    text = str(value).strip()
    text = re.sub(r"^(?:lib)?\d+_", "", text, flags=re.I)
    return re.sub(r"(?:-\d+)+$", "", text)


def validate_ploidy_calls(path: str, library: int,
                           barcodes_path: str) -> dict[str, object]:
    """Require one valid NN call for every selected filtered-MEX barcode."""
    with gzip.open(barcodes_path, "rt", encoding="utf-8") as handle:
        expected = [canonical_barcode(line) for line in handle if line.strip()]
    if not expected or len(expected) != len(set(expected)) or not all(expected):
        raise ValueError(
            f"lib{library} filtered MEX has empty/duplicate canonical barcodes")
    expected_set = set(expected)
    seen: set[str] = set()
    required = {
        "barcode", "library", "ploidy_call", "ploidy_probability",
        "classification_group", "prob_tetraploid", "qc_pass",
    }
    with open(path, "r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = set(reader.fieldnames or [])
        if not required <= fields:
            raise ValueError(
                f"lib{library} PLOIDY_NN calls lack fields: "
                + ",".join(sorted(required - fields)))
        for line_number, row in enumerate(reader, start=2):
            if None in row:
                raise ValueError(
                    f"lib{library} malformed PLOIDY_NN row {line_number}")
            barcode = canonical_barcode(row.get("barcode", ""))
            if not barcode or barcode in seen:
                raise ValueError(
                    f"lib{library} empty/duplicate PLOIDY_NN barcode at "
                    f"line {line_number}: {barcode!r}")
            try:
                row_library = int(str(row.get("library", "")).strip())
                confidence = float(row.get("ploidy_probability", ""))
                probability = float(row.get("prob_tetraploid", ""))
                qc_pass = int(str(row.get("qc_pass", "")).strip())
            except (TypeError, ValueError) as exc:
                raise ValueError(
                    f"lib{library} invalid PLOIDY_NN numeric value at "
                    f"line {line_number}") from exc
            call = str(row.get("ploidy_call", "")).strip().lower()
            expected_call = "tetraploid" if probability >= 0.5 else "diploid"
            expected_confidence = max(probability, 1.0 - probability)
            if (row_library != library or call not in {"diploid", "tetraploid"} or
                    call != expected_call or
                    not math.isfinite(confidence) or
                    not math.isfinite(probability) or
                    not 0.0 <= probability <= 1.0 or
                    not 0.5 <= confidence <= 1.0 or
                    abs(confidence - expected_confidence) > 1e-5 or
                    qc_pass not in {0, 1}):
                raise ValueError(
                    f"lib{library} inconsistent PLOIDY_NN call at "
                    f"line {line_number}")
            seen.add(barcode)
    if seen != expected_set:
        missing = sorted(expected_set - seen)[:5]
        extra = sorted(seen - expected_set)[:5]
        raise ValueError(
            f"lib{library} PLOIDY_NN/MEX barcode mismatch: "
            f"missing={missing}, extra={extra}")
    return {"library": library, "calls": len(seen), "status": "PASS"}


def mapping_validation_passes(
        path: str,
        required_rows: Sequence[tuple[str, str, str]]) -> tuple[bool, str]:
    """Verify exact modality, library, status, and path mapping rows."""
    try:
        with open(path, "r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if list(reader.fieldnames or []) != [
                    "modality", "library", "status", "path"]:
                return False, "unexpected schema"
            rows: dict[str, tuple[str, str, str]] = {}
            for row in reader:
                if None in row:
                    return False, "malformed row"
                row_path = absolute(row["path"])
                if row_path in rows:
                    return False, f"duplicate path row: {row_path}"
                rows[row_path] = (
                    row["modality"].strip(), row["library"].strip(),
                    row["status"].strip().upper())
        mismatches = []
        for required_path, required_modality, required_library in required_rows:
            observed = rows.get(absolute(required_path))
            expected = (required_modality, required_library, "OK")
            if observed != expected:
                mismatches.append(
                    f"{required_path} expected={expected!r} observed={observed!r}")
        if mismatches:
            return False, "missing/mismatched row(s): " + "; ".join(mismatches[:3])
        return True, f"PASS ({len(required_rows)} exact selected rows)"
    except (OSError, UnicodeError, csv.Error) as exc:
        return False, f"unreadable: {exc}"


def ledger_event_libraries(path: str,
                           libraries: Sequence[int]) -> set[int]:
    requested = set(libraries)
    opener = gzip.open if path.endswith(".gz") else open
    result: set[int] = set()
    with opener(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = set(reader.fieldnames or [])
        required = {
            "library", "barcode", "production_assignment",
            "production_assignment_source", "review_required",
            "downstream_release_status", "ambient_production_arm",
            "ambient_production_c", "event_id", "final_schema_version",
        }
        if not required <= fields:
            raise ValueError(
                "canonical final reconciliation ledger lacks required fields: "
                + ",".join(sorted(required - fields)))
        seen = set()
        for row in reader:
            raw_library = str(row.get("library", "")).strip().lower()
            raw_library = raw_library.removeprefix("lib")
            try:
                library = int(raw_library)
            except ValueError as exc:
                raise ValueError(
                    f"invalid library in final reconciliation ledger: "
                    f"{row.get('library')!r}") from exc
            if library not in requested:
                continue
            barcode = str(row.get("barcode", "")).strip()
            key = (library, barcode)
            if not barcode or key in seen:
                raise ValueError(
                    f"empty/duplicate selected cell in final reconciliation ledger: "
                    f"lib{library}/{barcode}")
            seen.add(key)
            schema = str(row.get("final_schema_version", "")).strip()
            if schema != "identity_reconciliation_final_v2_phase3_dispositions":
                raise ValueError(
                    "unsupported final reconciliation ledger schema: " + schema)
            event_id = str(row.get("event_id", "")).strip()
            if event_id and event_id.upper() not in {"NA", "N/A", "NONE", "."}:
                result.add(library)
    missing_libraries = requested - {library for library, _barcode in seen}
    if missing_libraries:
        raise ValueError(
            "final reconciliation ledger has no cells for: " +
            ",".join(f"lib{value}" for value in sorted(missing_libraries)))
    return result


def validate_final_assignment_bundle(
        ledger_path: str, paths: Sequence[LibraryPaths]) -> dict[int, int]:
    """Require every audited final assignment to equal the canonical ledger."""
    requested = {item.library for item in paths}
    ledger: dict[int, dict[str, str]] = {
        library: {} for library in requested}
    opener = gzip.open if ledger_path.endswith(".gz") else open
    with opener(ledger_path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            "library", "barcode", "production_assignment",
            "final_schema_version",
        }
        fields = set(reader.fieldnames or [])
        if not required <= fields:
            raise ValueError(
                "canonical final reconciliation ledger lacks assignment "
                "fields: " + ",".join(sorted(required - fields)))
        for row in reader:
            raw_library = str(row.get("library", "")).strip().lower()
            raw_library = raw_library.removeprefix("lib")
            try:
                library = int(raw_library)
            except ValueError:
                continue
            if library not in requested:
                continue
            if (str(row.get("final_schema_version", "")).strip() !=
                    "identity_reconciliation_final_v2_phase3_dispositions"):
                raise ValueError(
                    f"lib{library} canonical ledger has an unsupported schema")
            barcode = str(row.get("barcode", "")).strip()
            assignment = str(row.get("production_assignment", "")).strip()
            if not barcode or not assignment or barcode in ledger[library]:
                raise ValueError(
                    f"lib{library} canonical ledger has an empty/duplicate cell")
            ledger[library][barcode] = assignment

    counts: dict[int, int] = {}
    for item in paths:
        expected = ledger[item.library]
        if not expected:
            raise ValueError(
                f"canonical final ledger has no lib{item.library} cells")
        with open(item.samples, "r", encoding="utf-8") as handle:
            sample_ids = {
                line.strip() for line in handle if line.strip()}
        if not sample_ids:
            raise ValueError(f"lib{item.library} DEMUX sample order is empty")
        observed: dict[str, str] = {}
        with open(item.final_assignments, "r", encoding="utf-8") as handle:
            for line_number, raw in enumerate(handle, start=1):
                if not raw.strip():
                    continue
                fields = raw.rstrip("\r\n").split("\t")
                if len(fields) != 4:
                    raise ValueError(
                        f"lib{item.library} final assignments line "
                        f"{line_number} is not canonical four-column output")
                barcode, assignment, assignment_type, _score = (
                    value.strip() for value in fields)
                if not barcode or not assignment or barcode in observed:
                    raise ValueError(
                        f"lib{item.library} final assignments have an "
                        f"empty/duplicate cell at line {line_number}")
                donor_count = len([
                    token for token in assignment.split("+") if token.strip()])
                expected_type = (
                    "D" if assignment.startswith("M{") or donor_count >= 2
                    else "S")
                if assignment_type != expected_type:
                    raise ValueError(
                        f"lib{item.library}/{barcode} final assignment type "
                        "does not match its production assignment")
                if not assignment.startswith("M{"):
                    missing_donors = [
                        token.strip() for token in assignment.split("+")
                        if token.strip() not in sample_ids]
                    if missing_donors:
                        raise ValueError(
                            f"lib{item.library}/{barcode} final assignment "
                            "contains donor(s) absent from the main-panel "
                            f"sample order: {missing_donors}")
                observed[barcode] = assignment
        if observed != expected:
            missing = sorted(set(expected) - set(observed))[:5]
            extra = sorted(set(observed) - set(expected))[:5]
            changed = sorted(
                barcode for barcode in set(expected) & set(observed)
                if expected[barcode] != observed[barcode])[:5]
            raise ValueError(
                f"lib{item.library} final assignment/ledger mismatch: "
                f"missing={missing}, extra={extra}, changed={changed}")
        counts[item.library] = len(observed)
    return counts


def validate_finalization_summary(
        path: str, libraries: Sequence[int],
        assignment_counts: Mapping[int, int]) -> dict[str, object]:
    """Validate the Phase-3 finalizer's per-library and ALL accounting rows."""
    required = {
        "library", "input_barcodes", "output_ledger_rows",
        "accounting_status", "final_schema_version",
    }
    with open(path, "r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = set(reader.fieldnames or [])
        if not required <= fields:
            raise ValueError(
                "identity finalization run summary lacks fields: "
                + ",".join(sorted(required - fields)))
        rows = list(reader)
    by_key: dict[str, dict[str, str]] = {}
    for row in rows:
        raw = str(row.get("library", "")).strip()
        key = "ALL" if raw.upper() == "ALL" else raw.lower().removeprefix("lib")
        if key != "ALL":
            try:
                key = str(int(key))
            except ValueError as exc:
                raise ValueError(
                    f"invalid library in identity finalization summary: {raw!r}") \
                    from exc
        if key in by_key:
            raise ValueError(
                f"duplicate identity finalization summary row: {raw}")
        by_key[key] = row
    expected_keys = {str(int(value)) for value in libraries} | {"ALL"}
    if set(by_key) != expected_keys:
        raise ValueError(
            "identity finalization summary library coverage mismatch: "
            f"expected={sorted(expected_keys)} observed={sorted(by_key)}")

    total = 0
    for library in libraries:
        row = by_key[str(library)]
        try:
            input_count = int(str(row["input_barcodes"]).strip())
            output_count = int(str(row["output_ledger_rows"]).strip())
        except (TypeError, ValueError) as exc:
            raise ValueError(
                f"lib{library} identity finalization counts are invalid") from exc
        expected_count = int(assignment_counts.get(library, -1))
        if (input_count < 1 or input_count != output_count or
                output_count != expected_count or
                str(row["accounting_status"]).strip().upper() != "PASS" or
                str(row["final_schema_version"]).strip() !=
                "identity_reconciliation_final_v2_phase3_dispositions"):
            raise ValueError(
                f"lib{library} identity finalization accounting is not PASS")
        total += output_count
    overall = by_key["ALL"]
    try:
        overall_input = int(str(overall["input_barcodes"]).strip())
        overall_output = int(str(overall["output_ledger_rows"]).strip())
    except (TypeError, ValueError) as exc:
        raise ValueError("ALL identity finalization counts are invalid") from exc
    if (overall_input != total or overall_output != total or
            str(overall["accounting_status"]).strip().upper() != "PASS" or
            str(overall["final_schema_version"]).strip() !=
            "identity_reconciliation_final_v2_phase3_dispositions"):
        raise ValueError("ALL identity finalization accounting is not PASS")
    return {"status": "PASS", "libraries": len(libraries), "cells": total}


def validate_identity_metadata_manifest(
        path: str, workbook: str, panel_metadata: str,
        libraries: Sequence[int]) -> dict[str, object]:
    """Bind canonical reconciliation to its workbook and panel metadata."""
    try:
        with open(path, "r", encoding="utf-8") as handle:
            manifest = json.load(handle)
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise ValueError(f"invalid identity metadata manifest: {path}: {exc}") \
            from exc
    try:
        manifest_libraries = {
            int(value) for value in manifest.get("libraries", [])}
    except (TypeError, ValueError) as exc:
        raise ValueError("identity metadata manifest has invalid libraries") from exc
    if (manifest.get("schema_version") != "identity_reconciliation_v1" or
            absolute(str(manifest.get("workbook", ""))) != workbook or
            manifest.get("workbook_sha256") != sha256_file(workbook) or
            absolute(str(manifest.get("panel_metadata", ""))) != panel_metadata or
            manifest.get("panel_metadata_sha256") !=
            sha256_file(panel_metadata) or
            manifest.get("global_biological_line_source") != "2025_LineMeta" or
            manifest.get("ambient_rna_evaluated") is not False or
            not isinstance(manifest.get("n_global_biological_lines"), int) or
            int(manifest.get("n_global_biological_lines", 0)) < 1 or
            not isinstance(manifest.get("n_global_donors"), int) or
            int(manifest.get("n_global_donors", 0)) < 1 or
            not set(libraries) <= manifest_libraries):
        raise ValueError(
            "identity metadata manifest does not match the declared workbook, "
            "panel metadata, and upstream cohort")
    return {
        "status": "PASS", "libraries": len(libraries),
        "workbook_sha256": manifest["workbook_sha256"],
        "panel_metadata_sha256": manifest["panel_metadata_sha256"],
    }


def normalize_library_number(value: object) -> int | None:
    """Parse the library spellings accepted by the upstream metadata tools."""
    text = str(value).strip()
    match = re.search(r"(?:lib|RNA[_-]?)(\d+)$", text, re.I)
    if match:
        return int(match.group(1))
    if re.fullmatch(r"\d+(?:\.0+)?", text):
        return int(float(text))
    return None


def canonical_pool_identity(value: object, sample_ids: set[str] | None,
                            context: str, allow_x: bool = True,
                            max_components: int | None = 2) -> str:
    """Canonicalize a panel identity and reject unknown donor components."""
    raw = re.sub(r"\s+", "", str(value).strip())
    if allow_x:
        raw = raw.replace("x", "+")
    parts = [part for part in raw.split("+") if part]
    if (not raw or not parts or "+".join(parts) != raw or
            (max_components is not None and len(parts) > max_components)):
        raise ValueError(f"{context} has a malformed donor identity: {value!r}")
    missing = ([part for part in parts if part not in sample_ids]
               if sample_ids is not None else [])
    if missing:
        raise ValueError(
            f"{context} contains donor(s) absent from the main-panel sample "
            f"order: {missing}")
    return "+".join(sorted(parts))


def sample_orders(paths: Sequence[LibraryPaths]) -> dict[int, list[str]]:
    """Read unique, ordered main-panel sample identifiers for each library."""
    result: dict[int, list[str]] = {}
    for item in paths:
        with open(item.samples, "r", encoding="utf-8") as handle:
            samples = [line.strip().split()[0] for line in handle if line.strip()]
        if not samples or len(samples) != len(set(samples)):
            raise ValueError(
                f"lib{item.library} DEMUX sample order is empty or duplicated")
        result[item.library] = samples
    return result


def workbook_expected_pools(
        path: str, libraries: Sequence[int], sample_by_library: Mapping[int, list[str]],
        identity_semantics: bool) -> dict[int, set[str]]:
    """Reproduce the upstream workbook-to-expected-genotype conversion."""
    try:
        import pandas as pd
    except ImportError as exc:
        raise ValueError(
            "pandas/openpyxl are required to validate donor-pool workbooks") \
            from exc
    try:
        excel = pd.ExcelFile(path)
        if not {"convert", "libs"} <= set(excel.sheet_names):
            raise ValueError("required convert/libs sheets are absent")
        read_options = {"dtype": str} if identity_semantics else {}
        convert = pd.read_excel(
            path, sheet_name="convert", **read_options)
        libs = pd.read_excel(path, sheet_name="libs", **read_options)
        convert_extra = (
            pd.read_excel(path, sheet_name="convert_extra", dtype=str)
            if identity_semantics and "convert_extra" in excel.sheet_names
            else None)
    except ValueError:
        raise
    except Exception as exc:
        raise ValueError(f"cannot read donor-pool workbook {path}: {exc}") from exc

    def normalized_column(frame: object, *names: str) -> object | None:
        columns = getattr(frame, "columns")
        by_name = {
            re.sub(r"[^a-z0-9]", "", str(column).lower()): column
            for column in columns}
        for name in names:
            key = re.sub(r"[^a-z0-9]", "", name.lower())
            if key in by_name:
                return by_name[key]
        return None

    if identity_semantics:
        label_column = normalized_column(
            convert, "Library", "library_label", "line")
        vcf_column = normalized_column(
            convert, "VCF_ID", "vcf_id", "canonical_vcf_id")
        library_column = normalized_column(libs, "lib", "library", "LibNum")
        line_column = normalized_column(libs, "line", "Line", "WGS_Key")
    else:
        label_column = "Library" if "Library" in convert.columns else None
        vcf_column = "VCF_ID" if "VCF_ID" in convert.columns else None
        library_column = "lib" if "lib" in libs.columns else None
        line_column = "line" if "line" in libs.columns else None
    if (label_column is None or vcf_column is None or
            library_column is None or line_column is None):
        raise ValueError(
            f"donor-pool workbook lacks required conversion columns: {path}")

    aliases: dict[str, str] = {}

    def add_alias(label: str, vcf_id: str, source: str,
                  only_if_absent: bool = False) -> None:
        if not label or not vcf_id:
            return
        if label in aliases and aliases[label] != vcf_id:
            if only_if_absent:
                return
            raise ValueError(
                f"conflicting donor alias {label!r} in {source}: "
                f"{aliases[label]!r} versus {vcf_id!r}")
        aliases.setdefault(label, vcf_id)
        if identity_semantics:
            aliases.setdefault(vcf_id, vcf_id)

    for _, row in convert.iterrows():
        if pd.isna(row[label_column]) or pd.isna(row[vcf_column]):
            continue
        label = str(row[label_column])
        vcf_id = str(row[vcf_column])
        if identity_semantics:
            label = label.strip()
            vcf_id = vcf_id.strip()
        add_alias(label, vcf_id, path)
    if convert_extra is not None:
        extra_label = normalized_column(
            convert_extra, "Library", "library_label", "line")
        extra_vcf = normalized_column(
            convert_extra, "VCF_ID", "vcf_id", "canonical_vcf_id")
        if extra_label is not None and extra_vcf is not None:
            for _, row in convert_extra.iterrows():
                if pd.isna(row[extra_label]) or pd.isna(row[extra_vcf]):
                    continue
                label = str(row[extra_label]).strip()
                vcf_id = str(row[extra_vcf]).strip()
                add_alias(label, vcf_id, path, True)
                if label and vcf_id:
                    aliases.setdefault(vcf_id, aliases.get(label, vcf_id))
    if not aliases:
        raise ValueError(f"donor-pool workbook has no usable aliases: {path}")

    requested = set(int(value) for value in libraries)
    result: dict[int, set[str]] = {library: set() for library in requested}
    for row_number, (_, row) in enumerate(libs.iterrows(), start=2):
        if pd.isna(row[library_column]) or pd.isna(row[line_column]):
            continue
        if identity_semantics:
            library = normalize_library_number(row[library_column])
        else:
            library = next((value for value in requested
                            if row[library_column] == value), None)
        if library not in requested:
            continue
        raw = str(row[line_column]).strip()
        delimiter = r"[+x]" if identity_semantics else "x"
        parts = [part.strip() for part in re.split(delimiter, raw) if part.strip()]
        if not parts or len(parts) not in {1, 2}:
            raise ValueError(
                f"{path} libs row {row_number} has malformed identity {raw!r}")
        missing = [part for part in parts if part not in aliases]
        if missing:
            raise ValueError(
                f"{path} libs row {row_number} has unmapped line label(s): "
                f"{missing}")
        identity = canonical_pool_identity(
            "+".join(aliases[part] for part in parts),
            set(sample_by_library[library]),
            f"{path} libs row {row_number}")
        result[library].add(identity)
    empty = [library for library, values in result.items() if not values]
    if empty:
        raise ValueError(
            f"donor-pool workbook has no expected identities for: "
            + ",".join(f"lib{value}" for value in sorted(empty)))
    return result


def generated_expected_pools(
        paths: Sequence[LibraryPaths],
        sample_by_library: Mapping[int, list[str]]) -> dict[int, set[str]]:
    """Parse the exact expected-lines files consumed by forced DEMUX."""
    result: dict[int, set[str]] = {}
    for item in paths:
        path = os.path.join(
            os.path.dirname(item.mapping_bam),
            f"lib{item.library}_expected_lines.txt")
        values: list[str] = []
        with open(path, "r", encoding="utf-8") as handle:
            for line_number, raw in enumerate(handle, start=1):
                value = raw.strip()
                if not value:
                    continue
                values.append(canonical_pool_identity(
                    value, set(sample_by_library[item.library]),
                    f"{path} line {line_number}"))
        if not values or len(values) != len(set(values)):
            raise ValueError(
                f"lib{item.library} expected-lines file is empty or duplicated")
        result[item.library] = set(values)
    return result


def tsv_expected_pools(
        path: str, libraries: Sequence[int],
        sample_by_library: Mapping[int, list[str]]) -> dict[int, set[str]]:
    """Parse the expected-pool TSV using the upstream accepted columns."""
    with open(path, "r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = set(reader.fieldnames or [])
        library_column = next((name for name in (
            "library_id", "library", "lib", "Lib") if name in fields), None)
        identity_column = next((name for name in (
            "cell_identities", "expected_identities", "expected_identity",
            "identity", "identities") if name in fields), None)
        if library_column is None or identity_column is None:
            raise ValueError(
                "expected-pool metadata lacks library and identity columns")
        requested = set(int(value) for value in libraries)
        result: dict[int, set[str]] = {}
        for row_number, row in enumerate(reader, start=2):
            raw_library = str(row.get(library_column, "")).strip()
            match = re.fullmatch(r"(?i:lib)(\d+)", raw_library)
            if match:
                library = int(match.group(1))
            elif raw_library.isdigit():
                library = int(raw_library)
            else:
                library = None
            if library not in requested:
                continue
            if library in result:
                raise ValueError(
                    f"expected-pool metadata has duplicate lib{library} rows")
            raw_values = [value.strip() for value in
                          str(row.get(identity_column, "")).split(",")
                          if value.strip()]
            values = {
                canonical_pool_identity(
                    value, set(sample_by_library[library]),
                    f"{path} row {row_number}", False)
                for value in raw_values}
            if not values:
                raise ValueError(
                    f"expected-pool metadata has no identities for lib{library}")
            result[library] = values
    missing = set(int(value) for value in libraries) - set(result)
    if missing:
        raise ValueError(
            "expected-pool metadata lacks: " +
            ",".join(f"lib{value}" for value in sorted(missing)))
    return result


def validate_identity_metadata_tables(
        expected_genotypes_path: str, resolution_audit_path: str,
        warnings_path: str, uid_members_path: str,
        global_lines_path: str, global_donors_path: str,
        libraries: Sequence[int], sample_by_library: Mapping[int, list[str]],
        core_pools: Mapping[int, set[str]],
        manifest: Mapping[str, object]) -> dict[str, object]:
    """Validate reconciliation's expanded roster and critical metadata warnings."""
    with open(expected_genotypes_path, "rb") as handle:
        expected_bytes = handle.read()
    with open(resolution_audit_path, "rb") as handle:
        audit_bytes = handle.read()
    if expected_bytes != audit_bytes:
        raise ValueError(
            "identity expected-genotype and resolution-audit tables differ")
    expanded: dict[int, set[str]] = {int(value): set() for value in libraries}
    expected_uids: dict[tuple[int, str], set[str]] = {}
    with open(expected_genotypes_path, "r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            "library", "canonical_genotype", "uid_candidate_count",
            "uid_candidates", "reconciled_uid", "uid_resolution_status",
            "uid_resolution_scope", "donor_components",
            "expected_ploidy_class"}
        fields = set(reader.fieldnames or [])
        if not required <= fields:
            raise ValueError(
                "identity expected-genotype table lacks fields: " +
                ",".join(sorted(required - fields)))
        for row_number, row in enumerate(reader, start=2):
            library = normalize_library_number(row.get("library", ""))
            if library not in expanded:
                continue
            genotype = canonical_pool_identity(
                row.get("canonical_genotype", ""),
                set(sample_by_library[library]),
                f"{expected_genotypes_path} row {row_number}")
            if genotype in expanded[library]:
                raise ValueError(
                    f"identity expected-genotype table duplicates "
                    f"lib{library}/{genotype}")
            try:
                uid_count = int(str(row.get("uid_candidate_count", "")).strip())
            except ValueError as exc:
                raise ValueError(
                    f"identity expected-genotype table has an invalid UID "
                    f"count for lib{library}/{genotype}") from exc
            status = str(row.get("uid_resolution_status", "")).strip()
            raw_uids = str(row.get("uid_candidates", "")).strip()
            raw_reconciled = str(row.get("reconciled_uid", "")).strip()
            uid_values = {
                value.strip() for value in raw_uids.split("|")
                if value.strip() and value.strip().upper() not in {"NA", "."}}
            components = genotype.split("+")
            expected_ploidy = (
                "DIPLOID" if len(components) == 1 else
                "HOMOTYPIC_TETRAPLOID" if len(set(components)) == 1 else
                "HETEROTYPIC_TETRAPLOID")
            expected_status = (
                "EXACT_LIBRARY_METADATA_MATCH" if uid_count == 1 else
                "MULTIPLE_EXPECTED_UIDS_SAME_GENOTYPE")
            if (uid_count < 1 or len(uid_values) != uid_count or
                    raw_reconciled != raw_uids or status != expected_status or
                    str(row.get("uid_resolution_scope", "")).strip() !=
                    f"library:{library}" or
                    str(row.get("donor_components", "")).strip() !=
                    ",".join(components) or
                    str(row.get("expected_ploidy_class", "")).strip() !=
                    expected_ploidy or status not in {
                    "EXACT_LIBRARY_METADATA_MATCH",
                    "MULTIPLE_EXPECTED_UIDS_SAME_GENOTYPE"}):
                raise ValueError(
                    f"lib{library}/{genotype} has unresolved identity metadata: "
                    f"status={status!r}, uid_candidate_count={uid_count}")
            expanded[library].add(genotype)
            expected_uids[(library, genotype)] = uid_values
    for library in libraries:
        if not core_pools[library] <= expanded[library]:
            raise ValueError(
                f"lib{library} identity metadata omits DEMUX expected "
                f"genotype(s): {sorted(core_pools[library] - expanded[library])}")

    observed_uids: dict[tuple[int, str], set[str]] = {}
    with open(uid_members_path, "r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            "library", "canonical_genotype", "uid",
            "uid_resolution_scope"}
        fields = set(reader.fieldnames or [])
        if not required <= fields:
            raise ValueError(
                "identity UID-member table lacks fields: " +
                ",".join(sorted(required - fields)))
        for row_number, row in enumerate(reader, start=2):
            library = normalize_library_number(row.get("library", ""))
            if library not in expanded:
                continue
            genotype = canonical_pool_identity(
                row.get("canonical_genotype", ""),
                set(sample_by_library[library]),
                f"{uid_members_path} row {row_number}")
            uid = str(row.get("uid", "")).strip()
            uid_values = {
                value.strip() for value in re.split(r"[|,;]", uid)
                if value.strip() and value.strip().upper() not in {"NA", "."}}
            if (genotype not in expanded[library] or not uid_values or
                    str(row.get("uid_resolution_scope", "")).strip() !=
                    f"library:{library}"):
                raise ValueError(
                    f"invalid selected UID-member row in {uid_members_path}:"
                    f"{row_number}")
            observed_uids.setdefault((library, genotype), set()).update(uid_values)
    if observed_uids != expected_uids:
        raise ValueError(
            "identity UID-member rows do not reproduce selected expected UID sets")

    global_lines: dict[str, set[str]] = {}
    with open(global_lines_path, "r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            "canonical_genotype", "donor_components",
            "biological_ploidy_class", "n_line_meta_rows", "source_sheet",
            "source_workbook_sha256"}
        fields = set(reader.fieldnames or [])
        if not required <= fields:
            raise ValueError(
                "global biological-line table lacks fields: " +
                ",".join(sorted(required - fields)))
        for row_number, row in enumerate(reader, start=2):
            genotype = canonical_pool_identity(
                row.get("canonical_genotype", ""), None,
                f"{global_lines_path} row {row_number}", True, None)
            components = genotype.split("+")
            expected_ploidy = (
                "DIPLOID" if len(components) == 1 else
                "HOMOTYPIC_TETRAPLOID" if len(set(components)) == 1 else
                "HETEROTYPIC_TETRAPLOID" if len(components) == 2 else
                "OTHER")
            try:
                row_count = int(str(row.get("n_line_meta_rows", "")).strip())
            except ValueError as exc:
                raise ValueError(
                    f"invalid global line count at row {row_number}") from exc
            if (genotype in global_lines or row_count < 1 or
                    str(row.get("donor_components", "")).strip() !=
                    ",".join(components) or
                    str(row.get("biological_ploidy_class", "")).strip() !=
                    expected_ploidy or row.get("source_sheet") != "2025_LineMeta" or
                    row.get("source_workbook_sha256") !=
                    manifest.get("workbook_sha256")):
                raise ValueError(
                    f"invalid global biological line at row {row_number}")
            global_lines[genotype] = set(components)
    try:
        manifest_line_count = int(manifest.get("n_global_biological_lines", -1))
    except (TypeError, ValueError) as exc:
        raise ValueError("identity manifest has invalid global line count") from exc
    if manifest_line_count != len(global_lines):
        raise ValueError("identity manifest/global biological-line count mismatch")

    expected_donor_lines: dict[str, set[str]] = {}
    for genotype, components in global_lines.items():
        for donor in components:
            expected_donor_lines.setdefault(donor, set()).add(genotype)
    observed_donor_lines: dict[str, set[str]] = {}
    with open(global_donors_path, "r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            "donor_id", "n_global_biological_lines",
            "global_biological_lines", "source_sheet",
            "source_workbook_sha256"}
        fields = set(reader.fieldnames or [])
        if not required <= fields:
            raise ValueError(
                "global donor table lacks fields: " +
                ",".join(sorted(required - fields)))
        for row_number, row in enumerate(reader, start=2):
            donor = str(row.get("donor_id", "")).strip()
            lines = {
                value.strip() for value in
                str(row.get("global_biological_lines", "")).split("|")
                if value.strip()}
            try:
                line_count = int(str(
                    row.get("n_global_biological_lines", "")).strip())
            except ValueError as exc:
                raise ValueError(
                    f"invalid global donor count at row {row_number}") from exc
            if (not donor or donor in observed_donor_lines or
                    line_count != len(lines) or
                    row.get("source_sheet") != "2025_LineMeta" or
                    row.get("source_workbook_sha256") !=
                    manifest.get("workbook_sha256")):
                raise ValueError(f"invalid global donor at row {row_number}")
            observed_donor_lines[donor] = lines
    try:
        manifest_donor_count = int(manifest.get("n_global_donors", -1))
    except (TypeError, ValueError) as exc:
        raise ValueError("identity manifest has invalid global donor count") from exc
    if (observed_donor_lines != expected_donor_lines or
            manifest_donor_count != len(observed_donor_lines)):
        raise ValueError(
            "identity global donor table does not reproduce global lines")

    critical = {
        "MISSING_DONOR_ALIAS", "DONOR_NOT_IN_NUCLEAR_PANEL",
        "MISSING_UID", "NO_LIBRARY_METADATA_MATCH",
        "DUPLICATE_METADATA_ROW",
    }
    failures: list[str] = []
    with open(warnings_path, "r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"library", "canonical_genotype", "warning", "detail"}
        fields = set(reader.fieldnames or [])
        if not required <= fields:
            raise ValueError(
                "identity metadata warnings table lacks fields: " +
                ",".join(sorted(required - fields)))
        requested = set(int(value) for value in libraries)
        for row in reader:
            library = normalize_library_number(row.get("library", ""))
            warning = str(row.get("warning", "")).strip()
            if library in requested and warning in critical:
                failures.append(
                    f"lib{library}:{warning}:{row.get('detail', '')}")
            elif (str(row.get("library", "")).strip().upper() == "GLOBAL" and
                  warning == "MISSING_DONOR_ALIAS"):
                failures.append(
                    f"GLOBAL:{warning}:{row.get('detail', '')}")
    if failures:
        raise ValueError(
            "critical identity metadata warnings remain: " +
            "; ".join(failures[:10]))
    return {
        "status": "PASS", "libraries": len(libraries),
        "expanded_genotypes": sum(len(values) for values in expanded.values()),
        "global_biological_lines": len(global_lines),
        "global_donors": len(observed_donor_lines),
        "critical_warnings": 0,
    }


def panel_metadata_sample_ids(path: str) -> set[str]:
    """Read donor identifiers using identity reconciliation's header rules."""
    with open(path, "r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = list(reader.fieldnames or [])
        if not fields:
            raise ValueError("panel metadata has no header")
        id_column = next((column for column in fields if
                          re.sub(r"[^a-z0-9]", "", column.lower()) in
                          {"individ", "vcfid", "sample", "id"}), fields[0])
        result = {
            str(row.get(id_column, "")).strip() for row in reader
            if str(row.get(id_column, "")).strip() not in {
                "", ".", "NA", "na", "NaN", "nan", "None", "none",
                "null", "NULL", "unavailable", "UNAVAILABLE"}}
    if not result:
        raise ValueError("panel metadata contains no donor identifiers")
    return result


def clean_metadata_value(value: object) -> str:
    text = str(value).strip()
    if text.lower() in {
            "", ".", "na", "nan", "none", "null", "unavailable"}:
        return ""
    return text


def validate_final_ledger_uids(
        args: argparse.Namespace, paths: Sequence[LibraryPaths],
        libraries: Sequence[int]) -> dict[str, object]:
    """Bind every released production genotype to its recomputed UID set.

    The attached upstream finalizer can retain a preliminary UID after an
    explicit review changes ``production_assignment``.  Refuse that state so
    arm evidence cannot be pooled under the wrong biological UID.
    """
    sample_by_library = sample_orders(paths)
    local: dict[tuple[int, str], tuple[str, str]] = {}
    with open(args.identity_expected_genotypes, "r", encoding="utf-8",
              newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row_number, row in enumerate(reader, start=2):
            library = normalize_library_number(row.get("library", ""))
            if library not in sample_by_library:
                continue
            genotype = canonical_pool_identity(
                row.get("canonical_genotype", ""),
                set(sample_by_library[library]),
                f"{args.identity_expected_genotypes} row {row_number}")
            uid_set = clean_metadata_value(row.get("uid_candidates", ""))
            status = clean_metadata_value(row.get("uid_resolution_status", ""))
            local[(library, genotype)] = (uid_set, status)

    global_lines: dict[str, str] = {}
    with open(args.identity_global_lines, "r", encoding="utf-8",
              newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = set(reader.fieldnames or [])
        if not {"canonical_genotype", "uid_candidates"} <= fields:
            raise ValueError(
                "global biological-line table lacks UID candidate fields")
        for row_number, row in enumerate(reader, start=2):
            genotype = canonical_pool_identity(
                row.get("canonical_genotype", ""), None,
                f"{args.identity_global_lines} row {row_number}", True, None)
            uid_set = clean_metadata_value(row.get("uid_candidates", ""))
            if genotype in global_lines:
                raise ValueError(
                    f"duplicate global genotype at row {row_number}")
            global_lines[genotype] = uid_set

    released = 0
    global_resolutions = 0
    opener = gzip.open if args.ledger_input.endswith(".gz") else open
    with opener(args.ledger_input, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            "library", "barcode", "production_assignment",
            "uid_or_uid_set", "uid_resolution_status",
            "downstream_release_status"}
        fields = set(reader.fieldnames or [])
        if not required <= fields:
            raise ValueError(
                "final reconciliation ledger lacks UID fields: " +
                ",".join(sorted(required - fields)))
        requested = set(int(value) for value in libraries)
        for row_number, row in enumerate(reader, start=2):
            library = normalize_library_number(row.get("library", ""))
            if (library not in requested or
                    clean_metadata_value(
                        row.get("downstream_release_status", "")).upper() !=
                    "READY"):
                continue
            genotype = canonical_pool_identity(
                row.get("production_assignment", ""),
                set(sample_by_library[library]),
                f"{args.ledger_input} row {row_number}")
            if (library, genotype) in local:
                expected_uid, expected_status = local[(library, genotype)]
            elif genotype in global_lines:
                expected_uid = global_lines[genotype]
                uid_values = [
                    value for value in expected_uid.split("|") if value]
                expected_status = (
                    "EXACT_GLOBAL_METADATA_MATCH" if len(uid_values) == 1 else
                    "MULTIPLE_GLOBAL_UIDS_SAME_GENOTYPE"
                    if len(uid_values) > 1 else "GLOBAL_LINE_MISSING_UID")
                global_resolutions += 1
            else:
                raise ValueError(
                    f"released lib{library}/{row.get('barcode')} production "
                    f"genotype {genotype!r} has no local/global line metadata")
            observed_uid = clean_metadata_value(row.get("uid_or_uid_set", ""))
            observed_status = clean_metadata_value(
                row.get("uid_resolution_status", ""))
            if (not expected_uid or observed_uid != expected_uid or
                    observed_status != expected_status or
                    "MISSING" in expected_status or
                    "CONFLICT" in expected_status):
                raise ValueError(
                    f"released lib{library}/{row.get('barcode')} UID does not "
                    f"match production genotype {genotype}: expected "
                    f"{expected_uid or 'NA'}/{expected_status}, observed "
                    f"{observed_uid or 'NA'}/{observed_status}. This is the "
                    "upstream explicit-review UID carryover bug; do not run "
                    "arm CNV until reconciliation is corrected.")
            released += 1
    if released < 1:
        raise ValueError("final ledger has no released cells with validated UIDs")
    return {
        "status": "PASS", "released_cells": released,
        "global_scope_resolutions": global_resolutions,
    }


def validate_donor_pool_chain(
        args: argparse.Namespace, paths: Sequence[LibraryPaths],
        libraries: Sequence[int]) -> dict[str, object]:
    """Require all donor-roster sources used upstream to agree exactly."""
    sample_by_library = sample_orders(paths)
    demux_workbook = workbook_expected_pools(
        args.demux_pool_workbook, libraries, sample_by_library, False)
    generated = generated_expected_pools(paths, sample_by_library)
    pool_metadata = tsv_expected_pools(
        args.expected_pool_metadata, libraries, sample_by_library)
    identity_workbook = workbook_expected_pools(
        args.identity_metadata_workbook, libraries, sample_by_library, True)
    panel_ids = panel_metadata_sample_ids(args.panel_metadata)
    for library in libraries:
        sources = {
            "DEMUX workbook": demux_workbook[library],
            "generated expected-lines": generated[library],
            "expected-pool metadata": pool_metadata[library],
            "identity workbook libs sheet": identity_workbook[library],
        }
        reference = demux_workbook[library]
        mismatches = {
            name: sorted(values) for name, values in sources.items()
            if values != reference}
        if mismatches:
            raise ValueError(
                f"lib{library} donor-pool sources disagree: {mismatches}")
        expected_path = os.path.join(
            os.path.dirname(next(
                item.mapping_bam for item in paths
                if item.library == library)),
            f"lib{library}_expected_lines.txt")
        expected_bytes = "".join(
            f"{identity}\n" for identity in sorted(reference)).encode("utf-8")
        with open(expected_path, "rb") as handle:
            observed_bytes = handle.read()
        if observed_bytes != expected_bytes:
            raise ValueError(
                f"lib{library} generated expected-lines bytes do not match "
                "the DEMUX workbook conversion")
        missing_panel_ids = sorted({
            component for identity in reference
            for component in identity.split("+") if component not in panel_ids})
        if missing_panel_ids:
            raise ValueError(
                f"lib{library} expected donor(s) are absent from panel metadata: "
                f"{missing_panel_ids}")
    try:
        with open(args.identity_metadata_manifest, "r", encoding="utf-8") as handle:
            identity_manifest = json.load(handle)
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise ValueError(
            f"cannot load identity metadata manifest: {exc}") from exc
    metadata_detail = validate_identity_metadata_tables(
        args.identity_expected_genotypes, args.identity_resolution_audit,
        args.identity_metadata_warnings, args.identity_uid_members,
        args.identity_global_lines, args.identity_global_donors,
        libraries, sample_by_library, demux_workbook,
        identity_manifest)
    return {
        "status": "PASS", "libraries": len(libraries),
        "expected_identities": sum(
            len(values) for values in demux_workbook.values()),
        "identity_metadata": metadata_detail,
    }


def iter_strict_posthoc_tsv(
        path: str, expected_header: Sequence[str], label: str
) -> Iterable[dict[str, str]]:
    """Stream a TSV without DictReader's silent truncation behavior."""
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        try:
            header = next(reader)
        except StopIteration as exc:
            raise ValueError(f"{label} is empty: {path}") from exc
        if tuple(header) != tuple(expected_header):
            raise ValueError(
                f"{label} has a noncanonical header: {path}")
        for line_number, fields in enumerate(reader, start=2):
            if not fields or (len(fields) == 1 and fields[0] == ""):
                continue
            if len(fields) != len(header):
                raise ValueError(
                    f"{label} line {line_number} has {len(fields)} fields; "
                    f"expected {len(header)}: {path}")
            yield dict(zip(header, fields))


def strict_posthoc_tsv(
        path: str, expected_header: Sequence[str], label: str
) -> list[dict[str, str]]:
    """Materialize a small strict TSV; large ledgers use the iterator above."""
    return list(iter_strict_posthoc_tsv(path, expected_header, label))


def strict_assignment_rows(
        path: str, label: str
) -> tuple[list[str], dict[str, tuple[str, str, str]]]:
    """Read exact CellBouncer barcode/identity/type/LLR rows and order."""
    order: list[str] = []
    assignments: dict[str, tuple[str, str, str]] = {}
    with open(path, "r", encoding="utf-8", newline="") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) != 4:
                raise ValueError(
                    f"{label} line {line_number} does not have four fields: "
                    f"{path}")
            if any(field != field.strip() for field in fields):
                raise ValueError(
                    f"{label} line {line_number} has surrounding whitespace: "
                    f"{path}")
            barcode, identity, assignment_type, llr_text = fields
            if (not barcode or not identity or assignment_type not in {"S", "D"}
                    or barcode in assignments):
                raise ValueError(
                    f"{label} has an invalid or duplicate assignment at line "
                    f"{line_number}: {path}")
            canonical_identity = model_identity(identity)
            expected_type = "D" if "+" in canonical_identity else "S"
            if not canonical_identity or assignment_type != expected_type:
                raise ValueError(
                    f"{label} line {line_number} type does not match its "
                    f"identity: {path}")
            try:
                llr = float(llr_text)
            except ValueError as exc:
                raise ValueError(
                    f"{label} line {line_number} has invalid LLR: {path}") \
                    from exc
            if not math.isfinite(llr) or llr <= 0:
                raise ValueError(
                    f"{label} line {line_number} has nonpositive/nonfinite "
                    f"LLR: {path}")
            order.append(barcode)
            assignments[barcode] = (identity, assignment_type, llr_text)
    if not assignments:
        raise ValueError(f"{label} has no assignments: {path}")
    return order, assignments


def strict_assignment_map(path: str, label: str) -> dict[str, str]:
    """Read CellBouncer's headerless barcode/identity/S-D/LLR table."""
    _order, rows = strict_assignment_rows(path, label)
    assignments = {
        barcode: values[0] for barcode, values in rows.items()}
    return assignments


def unique_posthoc_barcodes(
        rows: Sequence[Mapping[str, str]], label: str, path: str
) -> dict[str, Mapping[str, str]]:
    indexed: dict[str, Mapping[str, str]] = {}
    for row_number, row in enumerate(rows, start=2):
        barcode = str(row.get("barcode", "")).strip()
        if not barcode or barcode in indexed:
            raise ValueError(
                f"{label} has an empty or duplicate barcode at line "
                f"{row_number}: {path}")
        indexed[barcode] = row
    return indexed


def validate_runner_up_coverage(
        path: str, required_barcodes: set[str], label: str
) -> None:
    rows = strict_posthoc_tsv(
        path, POSTHOC_RUNNER_UP_HEADER, label)
    observed: set[str] = set()
    barcode_ranks: set[tuple[str, int]] = set()
    for line_number, row in enumerate(rows, start=2):
        barcode = row["barcode"].strip()
        try:
            rank = int(row["rank"])
        except ValueError as exc:
            raise ValueError(
                f"{label} line {line_number} has invalid rank: {path}") from exc
        key = (barcode, rank)
        if (not barcode or rank < 1 or key in barcode_ranks or
                row["schema_version"] != "demux_parallel_runner_ups_v3"):
            raise ValueError(
                f"{label} line {line_number} has invalid barcode/rank/schema: "
                f"{path}")
        barcode_ranks.add(key)
        observed.add(barcode)
    missing = required_barcodes - observed
    if missing:
        raise ValueError(
            f"{label} lacks runner-up coverage for {len(missing)} core "
            f"barcode(s): {path}")


def verify_posthoc_generation(
        args: argparse.Namespace,
        paths: Sequence[LibraryPaths]) -> dict[str, object]:
    """Prove POSTHOC used the declared pool table and main-panel samples."""
    scripts_root = os.path.join(
        args.upstream_analysis_root, "aggregate_library_analysis",
        "slurm_scripts")
    audit_root = os.path.join(
        args.upstream_analysis_root, "aggregate_library_analysis", "posthoc")
    refined_root = os.path.join(
        args.upstream_analysis_root, "aggregate_library_analysis", "tetra_refine")
    sample_by_library = sample_orders(paths)
    expected_by_library = tsv_expected_pools(
        args.expected_pool_metadata, [item.library for item in paths],
        sample_by_library)
    distinct_sample_orders = {
        tuple(values) for values in sample_by_library.values()}
    if len(distinct_sample_orders) != 1:
        raise ValueError(
            "POSTHOC uses shared allowed-identity filenames but audited "
            "libraries have different main-panel sample orders")
    common_samples = next(iter(distinct_sample_orders))

    def posthoc_identity_sort_key(identity: str) -> tuple[int, str]:
        parts = identity.split("+")
        rank = 0 if len(parts) == 1 else 1 if len(set(parts)) == 1 else 2
        return rank, identity

    all_singlets_payload = "".join(
        f"{sample}\n" for sample in common_samples).encode("utf-8")
    allowed_identities = set(common_samples)
    for values in expected_by_library.values():
        allowed_identities.update(values)
    constrained_payload = "".join(
        f"{identity}\n" for identity in
        sorted(allowed_identities, key=posthoc_identity_sort_key)).encode("utf-8")
    records: list[dict[str, object]] = []
    for item in paths:
        path = os.path.join(scripts_root, f"posthoc_lib{item.library}.sbatch")
        records.append(provenance_file_record(
            path, f"lib{item.library} POSTHOC generation script", True))
        with open(path, "r", encoding="utf-8") as handle:
            script_payload = handle.read()
        if "Skipping (use --force to rerun)" in script_payload:
            raise ValueError(
                f"lib{item.library} POSTHOC script permits stale-output reuse; "
                "the upstream run was not generated with --force")
        audit_dir = os.path.join(audit_root, f"lib{item.library}")
        manifest_path = os.path.join(
            audit_dir, f"lib{item.library}.capabilities.json")
        out_prefix = os.path.join(audit_dir, f"lib{item.library}")
        prepare_commands = logical_program_argvs(
            path, "swap_audit_prepare.py")
        expected_prepare = [
            os.path.join(DEPLOYED_SCRIPTS, "swap_audit_prepare.py"),
            "--vcf_samples", item.samples,
            "--expected_pool_metadata", args.expected_pool_metadata,
            "--demux_root", os.path.dirname(item.demux_prefix),
            "--audit_root", audit_root,
            "--refined_assignments_root", refined_root,
            "--libraries", f"lib{item.library}",
            "--panel_metadata", args.panel_metadata,
            "--overwrite",
        ]
        matched_prepare = [argv for argv in prepare_commands
                           if list(argv) == expected_prepare]
        if len(prepare_commands) != 1 or len(matched_prepare) != 1:
            raise ValueError(
                f"lib{item.library} POSTHOC preparation is not tied to the "
                "declared samples, pool metadata, and panel metadata")

        summarize_commands = logical_program_argvs(
            path, "swap_audit_summarize.py")
        expected_summary = [
            os.path.join(DEPLOYED_SCRIPTS, "swap_audit_summarize.py"),
            "--lib", f"lib{item.library}",
            "--capabilities_manifest", manifest_path,
            "--audit_root", audit_root,
            "--expected_metadata", args.expected_pool_metadata,
            "--panel_metadata", args.panel_metadata,
            "--out_prefix", out_prefix,
            "--overwrite",
        ]
        matched_summary = [argv for argv in summarize_commands
                           if list(argv) == expected_summary]
        if len(summarize_commands) != 1 or len(matched_summary) != 1:
            raise ValueError(
                f"lib{item.library} POSTHOC summary is not tied to the "
                "declared pool and panel metadata")

        manifest_record = provenance_file_record(
            manifest_path, f"lib{item.library} POSTHOC capabilities", True)
        try:
            with open(manifest_path, "r", encoding="utf-8") as handle:
                manifest = json.load(handle)
        except (OSError, UnicodeError, json.JSONDecodeError) as exc:
            raise ValueError(
                f"invalid lib{item.library} POSTHOC capabilities: {exc}") from exc
        expected_core = {
            key: os.path.realpath(item.demux_prefix + suffix)
            for key, suffix in {
                "counts": ".counts", "samples": ".samples",
                "condf": ".condf", "assignments": ".assignments",
                "diagnostics": ".diagnostics.gz",
                "runner_ups": ".runner_ups.gz",
            }.items()}
        refined_candidates = [
            item.demux_prefix + ".refined_assignments",
            os.path.join(refined_root, f"lib{item.library}",
                         f"lib{item.library}.refined_assignments"),
            os.path.join(refined_root,
                         f"lib{item.library}.refined_assignments"),
            os.path.join(refined_root,
                         f"lib{item.library}_demuxed.refined_assignments"),
        ]
        expected_refined = next((
            os.path.realpath(candidate) for candidate in refined_candidates
            if os.path.exists(candidate)), "")
        if not expected_refined:
            raise ValueError(
                f"lib{item.library} has no current tetra_refine rich sidecar")
        optional_suffixes = {
            "species_counts": ".species_counts",
            "species_condf": ".species_condf",
            "species_samples": ".species_samples",
            "species_assignments": ".species_assignments",
            "atac_counts": ".atac.counts",
        }
        expected_optional: dict[str, str | None] = {
            key: (os.path.realpath(item.demux_prefix + suffix)
                  if os.path.exists(item.demux_prefix + suffix) else None)
            for key, suffix in optional_suffixes.items()}
        expected_optional.update({
            "call_qc": os.path.realpath(out_prefix + ".call_qc.tsv.gz"),
            "species_qc": os.path.realpath(out_prefix + ".species_qc.tsv"),
            "atac_call_qc": next((os.path.realpath(candidate) for candidate in (
                out_prefix + ".atac.call_qc.tsv.gz",
                item.demux_prefix + ".atac.call_qc.tsv.gz")
                if os.path.exists(candidate)), None),
            "refined_assignments": expected_refined,
            "thresholds": None,
            "panel_metadata": os.path.realpath(args.panel_metadata),
        })
        if any(expected_optional[key] is None for key in (
                "species_counts", "species_condf", "species_samples")):
            raise ValueError(
                f"lib{item.library} lacks the current species POSTHOC bundle")
        expected_audit_files = {
            "all_singlets": os.path.join(
                audit_root,
                f"all_{len(sample_by_library[item.library])}_individuals.txt"),
            "constrained_identities": os.path.join(
                audit_root,
                "global_allowed_singlets_homotypics_fusions.txt"),
            "unconstrained_prefix": out_prefix + "_audit_unconstrained",
            "constrained_prefix": out_prefix + "_audit_constrained",
        }
        for list_path, expected_payload in (
                (expected_audit_files["all_singlets"], all_singlets_payload),
                (expected_audit_files["constrained_identities"],
                 constrained_payload)):
            with open(list_path, "rb") as handle:
                observed_payload = handle.read()
            if observed_payload != expected_payload:
                raise ValueError(
                    f"lib{item.library} POSTHOC shared identity list has "
                    f"unexpected bytes: {list_path}")
        allowed = manifest.get("allowed_identities", {})
        optional = manifest.get("optional_files", {})
        if (manifest.get("schema_version") != "swap_audit_v1" or
                manifest.get("library") != f"lib{item.library}" or
                manifest.get("demux_prefix") !=
                os.path.realpath(item.demux_prefix) or
                manifest.get("core_files") != expected_core or
                not isinstance(optional, dict) or
                optional != expected_optional or
                manifest.get("global_files", {}).get("panel_metadata") !=
                os.path.realpath(args.panel_metadata) or
                manifest.get("audit_files") != expected_audit_files or
                not isinstance(allowed, dict) or
                allowed.get("allowed_identities_source") != "expected_pool" or
                allowed.get("allowed_identities_input") != "NA"):
            raise ValueError(
                f"lib{item.library} POSTHOC capabilities are not bound to "
                "the declared new-run inputs")
        try:
            manifest_identities = {
                canonical_pool_identity(
                    value, set(sample_by_library[item.library]),
                    f"{manifest_path} expected_identities", False)
                for value in manifest.get("expected_identities", [])}
        except TypeError as exc:
            raise ValueError(
                f"lib{item.library} POSTHOC expected identities are invalid") \
                from exc
        if manifest_identities != expected_by_library[item.library]:
            raise ValueError(
                f"lib{item.library} POSTHOC capabilities use a different "
                "expected donor pool")

        audit_command_path = os.path.join(
            audit_dir, f"lib{item.library}.run_audit_demux_parallel_commands.sh")
        score_command_path = os.path.join(
            audit_dir, f"lib{item.library}.run_tetra_score_calls_commands.sh")
        audit_commands = logical_program_argvs(
            audit_command_path, "demux_parallel")
        expected_audit_commands = [
            ["demux_parallel", "-o",
             expected_audit_files["unconstrained_prefix"], "--reuse_counts",
             "-i", expected_audit_files["all_singlets"],
             "-D", "0.5", "-N", "0", "--n_runner_ups", "250",
             "--close_threshold", "20", "-t", "8"],
            ["demux_parallel", "-o",
             expected_audit_files["constrained_prefix"], "--reuse_counts",
             "-I", expected_audit_files["constrained_identities"],
             "-D", "0.5", "-N", "0", "--n_runner_ups", "250",
             "--close_threshold", "20", "-t", "8"],
        ]
        if audit_commands != expected_audit_commands:
            raise ValueError(
                f"lib{item.library} POSTHOC audit command file is not canonical")
        score_commands = logical_program_argvs(
            score_command_path, "tetra_score_calls")
        expected_score_command = [
            "tetra_score_calls",
            "--counts", expected_core["counts"],
            "--samples", expected_core["samples"],
            "--assignments", expected_core["assignments"],
            "--diagnostics", expected_core["diagnostics"],
            "--runner_ups", expected_core["runner_ups"],
            "--output", expected_optional["call_qc"],
            "--libname", f"lib{item.library}",
            "--panel_metadata", expected_optional["panel_metadata"],
            "--species_counts", expected_optional["species_counts"],
            "--species_condf", expected_optional["species_condf"],
            "--species_samples", expected_optional["species_samples"],
        ]
        if score_commands != [expected_score_command]:
            raise ValueError(
                f"lib{item.library} POSTHOC score command file is not canonical")
        for generated in (
                audit_command_path, score_command_path,
                expected_audit_files["all_singlets"],
                expected_audit_files["constrained_identities"]):
            records.append(provenance_file_record(
                generated, f"lib{item.library} POSTHOC execution input", True))

        optional_inputs = [
            value for key, value in expected_optional.items()
            if value and key not in {"call_qc", "species_qc"}]
        input_paths = list(expected_core.values()) + [
            path, args.expected_pool_metadata, args.panel_metadata,
            *optional_inputs]
        newest_input = max(os.stat(value).st_mtime_ns for value in input_paths)
        validation_mtime = os.stat(args.identity_validation).st_mtime_ns
        manifest_mtime = os.stat(manifest_path).st_mtime_ns
        if not newest_input <= manifest_mtime <= validation_mtime:
            raise ValueError(
                f"lib{item.library} POSTHOC capabilities are outside the "
                "new-run input-to-identity-validation chronology")

        intermediate_paths = (
            out_prefix + ".call_qc.tsv.gz",
            out_prefix + ".species_qc.tsv",
            expected_audit_files["unconstrained_prefix"] + ".assignments",
            expected_audit_files["unconstrained_prefix"] + ".runner_ups.gz",
            expected_audit_files["constrained_prefix"] + ".assignments",
            expected_audit_files["constrained_prefix"] + ".runner_ups.gz",
        )
        execution_mtime = max(
            manifest_mtime, os.stat(audit_command_path).st_mtime_ns,
            os.stat(score_command_path).st_mtime_ns)
        for output in intermediate_paths:
            record = provenance_file_record(
                output, f"lib{item.library} POSTHOC intermediate")
            output_mtime = os.stat(output).st_mtime_ns
            if output_mtime < execution_mtime or output_mtime > validation_mtime:
                raise ValueError(
                    f"lib{item.library} POSTHOC intermediate is outside the "
                    f"new-run chronology: {output}")
            records.append(record)
        newest_intermediate = max(
            os.stat(output).st_mtime_ns for output in intermediate_paths)
        summary_paths = (
            out_prefix + ".swap_report.tsv",
            out_prefix + ".swap_scores.tsv",
            out_prefix + ".switch_matrix.tsv",
            out_prefix + ".identity_fractions.tsv",
        )
        for output in summary_paths:
            output_record = provenance_file_record(
                output, f"lib{item.library} POSTHOC output")
            output_mtime = os.stat(output).st_mtime_ns
            if output_mtime < newest_intermediate or output_mtime > validation_mtime:
                raise ValueError(
                    f"lib{item.library} POSTHOC output is outside the "
                    f"new-run chronology: {output}")
            records.append(output_record)

        # The attached summarizer intentionally intersects barcode universes
        # and supports reduced modes.  That is useful interactively but unsafe
        # as a production boundary: a truncated audit or header-only QC file can
        # still produce a small PASS report.  Bind every per-cell product to the
        # complete core DEMUX assignment universe before accepting the summary.
        core_assignments = strict_assignment_map(
            expected_core["assignments"],
            f"lib{item.library} core DEMUX assignments")
        core_barcodes = set(core_assignments)
        unconstrained_path = (
            expected_audit_files["unconstrained_prefix"] + ".assignments")
        constrained_path = (
            expected_audit_files["constrained_prefix"] + ".assignments")
        unconstrained = strict_assignment_map(
            unconstrained_path,
            f"lib{item.library} unconstrained POSTHOC assignments")
        constrained = strict_assignment_map(
            constrained_path,
            f"lib{item.library} constrained POSTHOC assignments")
        if set(unconstrained) != core_barcodes:
            raise ValueError(
                f"lib{item.library} unconstrained POSTHOC barcode universe "
                "does not exactly match core DEMUX assignments")
        if set(constrained) != core_barcodes:
            raise ValueError(
                f"lib{item.library} constrained POSTHOC barcode universe "
                "does not exactly match core DEMUX assignments")

        call_qc_path = out_prefix + ".call_qc.tsv.gz"
        call_qc_rows = strict_posthoc_tsv(
            call_qc_path, POSTHOC_CALL_QC_HEADER,
            f"lib{item.library} POSTHOC call QC")
        call_qc = unique_posthoc_barcodes(
            call_qc_rows, f"lib{item.library} POSTHOC call QC", call_qc_path)
        if set(call_qc) != core_barcodes:
            raise ValueError(
                f"lib{item.library} POSTHOC call-QC barcode universe does "
                "not exactly match core DEMUX assignments")
        for barcode, row in call_qc.items():
            if (row["libname"] != f"lib{item.library}" or
                    row["assignment"] != core_assignments[barcode]):
                raise ValueError(
                    f"lib{item.library} POSTHOC call QC disagrees with core "
                    f"DEMUX for barcode {barcode}")

        species_qc_path = out_prefix + ".species_qc.tsv"
        species_qc_rows = strict_posthoc_tsv(
            species_qc_path, POSTHOC_SPECIES_QC_HEADER,
            f"lib{item.library} POSTHOC species QC")
        if len(species_qc_rows) != 1:
            raise ValueError(
                f"lib{item.library} POSTHOC species QC must have one row")
        species_qc = species_qc_rows[0]
        try:
            species_evidence_cells = int(
                species_qc["n_cells_with_species_evidence"])
        except ValueError as exc:
            raise ValueError(
                f"lib{item.library} POSTHOC species QC has invalid cell count") \
                from exc
        disabling_species_warnings = {
            "NO_SPECIES_INPUTS", "SPECIES_SCORING_DISABLED",
            "SPECIES_COUNTS_MISSING", "SPECIES_SAMPLES_MISSING",
            "PANEL_METADATA_MISSING",
        }
        observed_species_warnings = {
            value for value in species_qc["warnings"].split(",") if value}
        if (species_qc["library"] != f"lib{item.library}" or
                species_evidence_cells < 1 or
                disabling_species_warnings & observed_species_warnings):
            raise ValueError(
                f"lib{item.library} POSTHOC species QC is not full mode")

        validate_runner_up_coverage(
            expected_core["runner_ups"], core_barcodes,
            f"lib{item.library} core DEMUX runner ups")
        validate_runner_up_coverage(
            expected_audit_files["unconstrained_prefix"] + ".runner_ups.gz",
            core_barcodes,
            f"lib{item.library} unconstrained POSTHOC runner ups")
        validate_runner_up_coverage(
            expected_audit_files["constrained_prefix"] + ".runner_ups.gz",
            core_barcodes,
            f"lib{item.library} constrained POSTHOC runner ups")

        score_rows = strict_posthoc_tsv(
            summary_paths[1], POSTHOC_SCORE_HEADER,
            f"lib{item.library} POSTHOC swap scores")
        score_by_barcode = unique_posthoc_barcodes(
            score_rows, f"lib{item.library} POSTHOC swap scores",
            summary_paths[1])
        if set(score_by_barcode) != core_barcodes:
            raise ValueError(
                f"lib{item.library} POSTHOC swap-score barcode universe does "
                "not exactly match core DEMUX assignments")
        for barcode, row in score_by_barcode.items():
            if (row["libname"] != f"lib{item.library}" or
                    row["raw_demux_assignment"] != core_assignments[barcode] or
                    row["original_assignment"] != core_assignments[barcode] or
                    row["audit_unconstrained_best"] !=
                    unconstrained[barcode] or
                    row["audit_constrained_best"] != constrained[barcode]):
                raise ValueError(
                    f"lib{item.library} POSTHOC swap score disagrees with its "
                    f"per-cell inputs for barcode {barcode}")

        report_rows = strict_posthoc_tsv(
            summary_paths[0], POSTHOC_REPORT_HEADER,
            f"lib{item.library} POSTHOC report")
        if len(report_rows) != 1:
            raise ValueError(
                f"lib{item.library} POSTHOC report does not have one row")
        report = report_rows[0]
        try:
            report_cells = int(str(report.get("n_cells", "")).strip())
            report_identities = {
                canonical_pool_identity(
                    value, set(sample_by_library[item.library]),
                    f"lib{item.library} POSTHOC report", False)
                for value in str(report.get("expected_identity", "")).split(",")
                if value.strip()}
        except ValueError as exc:
            raise ValueError(
                f"lib{item.library} POSTHOC report has invalid values") from exc
        full_feature_modes = {
            "CORE_PLUS_CALL_QC_AND_SPECIES",
            "CORE_PLUS_CALL_QC_AND_SPECIES_AND_ATAC",
        }
        if (report.get("library") != f"lib{item.library}" or
                report.get("preflight_status") != "PASS" or
                report_cells != len(core_barcodes) or
                report.get("has_call_qc") != "true" or
                report.get("has_species_qc") != "true" or
                report.get("has_refined_assignments") != "true" or
                report.get("feature_mode") not in full_feature_modes or
                not str(report.get("audit_verdict", "")).strip() or
                report.get("audit_verdict") == "NOT_RUN" or
                report_identities != expected_by_library[item.library]):
            raise ValueError(
                f"lib{item.library} POSTHOC report is not a completed "
                "new-run donor audit")
        records.append(manifest_record)
    return {"status": "PASS", "records": records}


def validate_panel_distinguishability(
        path: str, samples: Sequence[str]) -> dict[str, object]:
    """Validate the complete deterministic donor-pair table and its algebra."""
    rows = strict_posthoc_tsv(
        path, ("donor1", "donor2", "joint_callable_sites",
               "discordant_genotype_sites", "expected_information",
               "confusability_status"),
        "nuclear panel distinguishability output")
    expected_pairs = {
        (samples[left], samples[right])
        for left in range(len(samples))
        for right in range(left + 1, len(samples))}
    observed: set[tuple[str, str]] = set()
    for row in rows:
        pair = (row["donor1"].strip(), row["donor2"].strip())
        if pair not in expected_pairs or pair in observed:
            raise ValueError(
                "nuclear panel distinguishability has an unexpected or "
                f"duplicate donor pair: {pair}")
        observed.add(pair)
        try:
            joint = int(row["joint_callable_sites"])
            discordant = int(row["discordant_genotype_sites"])
            information = float(row["expected_information"])
        except ValueError as exc:
            raise ValueError(
                f"nuclear panel distinguishability has nonnumeric values for "
                f"{pair}") from exc
        expected_status = (
            "NO_JOINT_CALLABLE_SITES" if joint == 0 else
            "INDISTINGUISHABLE_ON_PANEL" if discordant == 0 else
            "DISTINGUISHABLE_ON_PANEL")
        if (joint < 0 or discordant < 0 or discordant > joint or
                not math.isfinite(information) or information < 0 or
                information < 0.25 * discordant - 1e-8 or
                information > discordant + 1e-8 or
                ((discordant > 0) != (information > 0)) or
                row["confusability_status"].strip() != expected_status):
            raise ValueError(
                f"nuclear panel distinguishability accounting is invalid for "
                f"{pair}")
    if observed != expected_pairs:
        raise ValueError(
            "nuclear panel distinguishability does not contain every ordered "
            "main-panel donor pair")
    return {"status": "PASS", "samples": len(samples), "pairs": len(rows)}


def verify_identity_metadata_generation(
        args: argparse.Namespace,
        paths: Sequence[LibraryPaths]) -> dict[str, object]:
    """Bind reconciliation metadata outputs to the deployed producer code."""
    scripts_root = os.path.join(
        args.upstream_analysis_root, "aggregate_library_analysis",
        "slurm_scripts")
    script_path = os.path.join(scripts_root, "identity_metadata.sbatch")
    metadata_root = os.path.dirname(args.identity_metadata_manifest)
    helper = args.identity_metadata_helper
    panel_utility = args.panel_distinguishability_binary
    metadata_commands = logical_program_argvs(
        script_path, "identity_reconciliation.py")
    expected_metadata_command = [
        helper, "metadata",
        "--xlsx", args.identity_metadata_workbook,
        "--panel-metadata", args.panel_metadata,
        "--output-root", metadata_root,
    ]
    matched_metadata = [argv for argv in metadata_commands
                        if list(argv) == expected_metadata_command]
    if len(metadata_commands) != 1 or len(matched_metadata) != 1:
        raise ValueError(
            "identity metadata script is not tied to the declared workbook, "
            "panel metadata, and output root")
    panel_commands = logical_program_argvs(
        script_path, "nuclear_panel_distinguishability")
    distinguishability = os.path.join(
        metadata_root, "nuclear_panel_distinguishability.tsv")
    expected_panel_command = [
        panel_utility, "--vcf", args.interindividual_panel,
        "--output", distinguishability,
    ]
    matched_panel = [argv for argv in panel_commands
                     if list(argv) == expected_panel_command]
    if len(panel_commands) != 1 or len(matched_panel) != 1:
        raise ValueError(
            "identity metadata script is not tied to the declared main panel")
    records = [
        provenance_file_record(
            script_path, "identity metadata generation script", True),
        provenance_file_record(
            helper, "identity reconciliation metadata helper", True),
        provenance_file_record(
            panel_utility, "nuclear panel distinguishability binary", True),
        provenance_file_record(
            distinguishability, "nuclear panel distinguishability output", True),
    ]
    producer_mtime = max(
        os.stat(script_path).st_mtime_ns, os.stat(helper).st_mtime_ns,
        os.stat(panel_utility).st_mtime_ns,
        os.stat(args.identity_metadata_workbook).st_mtime_ns,
        os.stat(args.panel_metadata).st_mtime_ns)
    manifest_mtime = os.stat(args.identity_metadata_manifest).st_mtime_ns
    if producer_mtime > manifest_mtime:
        raise ValueError(
            "identity metadata manifest predates its generation script, helper, "
            "panel utility, workbook, or panel metadata")
    if not manifest_mtime <= os.stat(distinguishability).st_mtime_ns <= \
            os.stat(args.identity_validation).st_mtime_ns:
        raise ValueError(
            "nuclear panel distinguishability output is outside the current "
            "identity-generation chronology")
    orders = {tuple(order) for order in sample_orders(paths).values()}
    if len(orders) != 1:
        raise ValueError(
            "selected libraries do not share one main-panel sample order")
    distinguishability_detail = validate_panel_distinguishability(
        distinguishability, next(iter(orders)))
    deterministic_outputs = (
        "donor_aliases.tsv",
        "global_biological_lines.tsv",
        "global_donors.tsv",
        "library_uid_members.tsv",
        "library_expected_genotypes.tsv",
        "library_resolution_audit.tsv",
        "genotype_to_physical_lines.tsv",
        "metadata_warnings.tsv",
        "preparation_relationships.tsv",
        "metadata_manifest.json",
    )
    with tempfile.TemporaryDirectory(
            prefix="tetra_arm_identity_metadata.") as temporary:
        command = [
            sys.executable, helper, "metadata",
            "--xlsx", args.identity_metadata_workbook,
            "--panel-metadata", args.panel_metadata,
            "--output-root", temporary,
        ]
        result = subprocess.run(
            command, capture_output=True, text=True, check=False)
        if result.returncode != 0:
            detail = result.stderr.strip() or result.stdout.strip()
            raise ValueError(
                "cannot reproduce identity metadata with the bound helper: "
                + detail[:1000])
        for name in deterministic_outputs:
            reproduced = os.path.join(temporary, name)
            selected = os.path.join(metadata_root, name)
            if (not regular_nonempty(reproduced) or
                    not regular_nonempty(selected) or
                    sha256_file(reproduced) != sha256_file(selected)):
                raise ValueError(
                    f"selected identity metadata does not reproduce from the "
                    f"bound workbook/helper: {name}")
            records.append(provenance_file_record(
                selected, f"reproduced identity metadata {name}", True))
    return {
        "status": "PASS",
        "panel_distinguishability": distinguishability_detail,
        "records": records,
    }


def verify_demux_generation(args: argparse.Namespace,
                            paths: Sequence[LibraryPaths]) -> dict[str, object]:
    """Tie every sidecar to a forced DEMUX script and the selected BAM."""
    scripts_root = os.path.join(
        args.upstream_analysis_root, "aggregate_library_analysis",
        "slurm_scripts")
    panel_paths = (
        args.interindividual_panel, args.het_panel, args.species_panel)
    holder_candidates = sorted(glob.glob(os.path.join(
        scripts_root, "vcf_daemon_holder_*.sbatch")))
    allowed_reference_bams = {item.mapping_bam for item in paths}
    matching_holders: dict[str, dict[str, str]] = {}
    for candidate in holder_candidates:
        try:
            contract = vcf_holder_contract(candidate)
        except (OSError, UnicodeError, ValueError):
            continue
        if (contract["vcf"] == args.interindividual_panel and
                contract["het_vcf"] == args.het_panel and
                contract["species_vcf"] == args.species_panel and
                contract["daemon"] == os.path.join(
                    DEPLOYED_BIN, "vcf_loader_daemon") and
                contract["bam"] in allowed_reference_bams):
            matching_holders[candidate] = contract
    if not matching_holders:
        raise ValueError(
            "no upstream managed-VCF holder script ties the selected mapping "
            "root to the declared interindividual/HET/species panels under: "
            + scripts_root)

    demux_records = []
    used_holders: set[str] = set()
    # The daemon consumes panel BCF bytes directly and DEMUX consumes the
    # filtered barcode roster and panel metadata.  All resulting sidecars must
    # be at least as new as every one of those scientific inputs.
    newest_shared_input_mtime = max(
        os.stat(path).st_mtime_ns
        for path in (*panel_paths, args.panel_metadata))
    for item in paths:
        script_path = os.path.join(scripts_root, f"demux_lib{item.library}.sbatch")
        record = provenance_file_record(
            script_path, f"lib{item.library} DEMUX generation script", True)
        commands = shell_command_argvs(script_path, "demux_parallel")
        raw_prefix = os.path.join(
            os.path.dirname(item.demux_prefix), f"lib{item.library}_raw")
        expected_lines = os.path.join(
            os.path.dirname(item.mapping_bam),
            f"lib{item.library}_expected_lines.txt")
        demux_binary = os.path.join(DEPLOYED_BIN, "demux_parallel")

        def exact_filtered_command(argv: Sequence[str]) -> bool:
            shared = argv_value(argv, "--shared_vcf")
            if not shared:
                return False
            expected = [
                demux_binary, "-b", item.mapping_bam,
                "-o", item.demux_prefix,
                "--shared_vcf", shared,
                "--shared_het_vcf", shared + "_het",
                "-f", "--barcodes", item.expression_barcodes,
                "-I", expected_lines,
                "--species_shared_vcf", shared + "_species",
                "--species_counts_output", "--species_assignment_output",
                "--panel_metadata", args.panel_metadata,
                "--species_panel_mode", "count_only",
                "--dump_selection_audit",
                "--dump_pileup", item.demux_prefix,
                "--force_recount", "-t", "80",
            ]
            return list(argv) == expected

        def exact_raw_command(argv: Sequence[str]) -> bool:
            shared = argv_value(argv, "--shared_vcf")
            if not shared:
                return False
            expected = [
                demux_binary, "-b", item.mapping_bam,
                "-o", raw_prefix,
                "--shared_vcf", shared, "-f", "-I", expected_lines,
                "--species_shared_vcf", shared + "_species",
                "--species_counts_output",
                "--panel_metadata", args.panel_metadata,
                "--species_panel_mode", "count_only",
                "--skip_assignment", "--force_recount", "-t", "8",
            ]
            return list(argv) == expected

        matched_commands = [argv for argv in commands
                            if exact_filtered_command(argv)]
        matched_raw_commands = [argv for argv in commands
                                if exact_raw_command(argv)]
        if (len(commands) != 2 or len(matched_commands) != 1 or
                len(matched_raw_commands) != 1):
            raise ValueError(
                f"lib{item.library} DEMUX script does not contain the exact "
                "forced filtered pileup pass and raw empty-drop counting pass "
                f"for the selected BAM: {script_path}")
        command = matched_commands[0]
        raw_command = matched_raw_commands[0]
        shared_segment = argv_value(command, "--shared_vcf")
        het_segment = argv_value(command, "--shared_het_vcf")
        species_segment = argv_value(command, "--species_shared_vcf")
        if (not shared_segment or het_segment != shared_segment + "_het" or
                species_segment != shared_segment + "_species" or
                argv_value(raw_command, "--shared_vcf") != shared_segment or
                argv_value(raw_command, "--species_shared_vcf") !=
                species_segment):
            raise ValueError(
                f"lib{item.library} DEMUX shared-memory roles do not match "
                "the managed main/HET/species segment contract")
        compatible_holders = [
            path for path, contract in matching_holders.items()
            if contract["base"] == shared_segment]
        demux_nodes = sbatch_nodelist(script_path)
        holder_by_node: dict[str, str] = {}
        ambiguous_nodes = set()
        for holder in compatible_holders:
            nodes = sbatch_nodelist(holder)
            if len(nodes) != 1:
                continue
            node = next(iter(nodes))
            if node in holder_by_node:
                ambiguous_nodes.add(node)
            holder_by_node[node] = holder
        if (not demux_nodes or ambiguous_nodes or
                demux_nodes != set(holder_by_node)):
            raise ValueError(
                f"lib{item.library} DEMUX main-panel shared-memory segment is "
                "not tied one-to-one to its allowed daemon node holders")
        selected_holders = [holder_by_node[node] for node in sorted(demux_nodes)]
        holder_mtime = max(os.stat(path).st_mtime_ns for path in selected_holders)
        bam_mtime = os.stat(item.mapping_bam).st_mtime_ns
        barcodes_mtime = os.stat(item.expression_barcodes).st_mtime_ns
        expected_lines_mtime = os.stat(expected_lines).st_mtime_ns
        script_mtime = os.stat(script_path).st_mtime_ns
        if holder_mtime > script_mtime:
            raise ValueError(
                f"lib{item.library} matching VCF holder was generated after "
                "the DEMUX script; the generation cannot be tied uniquely")
        newest_demux_mtime = 0
        filtered_generated = (
            item.demux_prefix + ".counts",
            item.demux_prefix + ".assignments",
            item.demux_prefix + ".summary",
            item.demux_prefix + ".diagnostics.gz",
            item.demux_prefix + ".runner_ups.gz",
            item.samples,
            item.demux_prefix + ".species_counts",
            item.demux_prefix + ".species_condf",
            item.demux_prefix + ".species_samples",
            item.demux_prefix + ".species_assignments",
            item.pileup_sites,
            item.pileup_molecules,
            item.pileup_observations,
        )
        for sidecar in filtered_generated:
            sidecar_mtime = os.stat(sidecar).st_mtime_ns
            newest_demux_mtime = max(newest_demux_mtime, sidecar_mtime)
            if (sidecar_mtime < bam_mtime or
                    sidecar_mtime < barcodes_mtime or
                    sidecar_mtime < expected_lines_mtime or
                    sidecar_mtime < script_mtime or
                    sidecar_mtime < holder_mtime or
                    sidecar_mtime < newest_shared_input_mtime):
                raise ValueError(
                    f"lib{item.library} DEMUX artifact predates its selected "
                    "BAM, barcode roster, expected-pool file, panel, panel "
                    "metadata, VCF holder, or generation script: " + sidecar)
            demux_records.append(provenance_file_record(
                sidecar, f"lib{item.library} forced filtered DEMUX artifact"))
        newest_filtered_mtime = max(
            os.stat(sidecar).st_mtime_ns for sidecar in filtered_generated)

        raw_generated = (
            raw_prefix + ".counts",
            raw_prefix + ".species_counts",
            raw_prefix + ".species_condf",
        )
        for sidecar in raw_generated:
            sidecar_mtime = os.stat(sidecar).st_mtime_ns
            newest_demux_mtime = max(newest_demux_mtime, sidecar_mtime)
            if (sidecar_mtime < bam_mtime or
                    sidecar_mtime < expected_lines_mtime or
                    sidecar_mtime < script_mtime or
                    sidecar_mtime < holder_mtime or
                    sidecar_mtime < newest_shared_input_mtime):
                raise ValueError(
                    f"lib{item.library} raw DEMUX artifact predates its "
                    "selected BAM, expected-pool file, panel, panel metadata, "
                    "VCF holder, or forced two-pass script: " + sidecar)
            demux_records.append(provenance_file_record(
                sidecar, f"lib{item.library} forced raw DEMUX artifact"))
        if min(os.stat(sidecar).st_mtime_ns for sidecar in raw_generated) < \
                newest_filtered_mtime:
            raise ValueError(
                f"lib{item.library} raw DEMUX Pass 2 predates a filtered "
                "Pass-1 artifact; the two-pass bundle is mixed")

        raw_metadata = (
            raw_prefix + ".samples",
            raw_prefix + ".condf",
            raw_prefix + ".species_samples",
        )
        for sidecar in raw_metadata:
            demux_records.append(provenance_file_record(
                sidecar, f"lib{item.library} raw DEMUX metadata", True))
        if (os.path.realpath(raw_prefix + ".samples") !=
                os.path.realpath(item.samples) or
                sha256_file(raw_prefix + ".species_samples") !=
                sha256_file(item.demux_prefix + ".species_samples") or
                os.path.realpath(raw_prefix + ".condf") !=
                os.path.realpath(item.demux_prefix + ".condf")):
            raise ValueError(
                f"lib{item.library} raw/filtered DEMUX metadata bundles do not "
                "share the same main/species sample order and CONDF")
        for downstream in (item.final_assignments,):
            if os.stat(downstream).st_mtime_ns < newest_demux_mtime:
                raise ValueError(
                    f"lib{item.library} reconciled/ploidy artifact predates "
                    f"the newest DEMUX bundle: {downstream}")
        demux_records.append(record)
        demux_records.append(provenance_file_record(
            expected_lines, f"lib{item.library} DEMUX expected donor pool", True))
        used_holders.update(selected_holders)
    return {
        "status": "PASS",
        "demux_scripts": demux_records,
        "vcf_holder_scripts": [
            provenance_file_record(path, "managed VCF holder script", True)
            for path in sorted(used_holders)],
    }


def assignment_identity_map(path: str, label: str) -> dict[str, str]:
    """Read a four-column assignment export while allowing NA final scores."""
    result: dict[str, str] = {}
    with open(path, "r", encoding="utf-8", newline="") as handle:
        for line_number, raw in enumerate(handle, start=1):
            if not raw.strip():
                continue
            fields = [value.strip() for value in raw.rstrip("\r\n").split("\t")]
            if len(fields) != 4 or not fields[0] or not fields[1]:
                raise ValueError(
                    f"{label} has a malformed row at line {line_number}: {path}")
            if fields[0] in result:
                raise ValueError(
                    f"{label} has duplicate barcode {fields[0]}: {path}")
            result[fields[0]] = fields[1]
    if not result:
        raise ValueError(f"{label} has no rows: {path}")
    return result


def model_identity(value: str) -> str:
    """Collapse biological A+A to the SNP-resolvable CellBouncer identity A."""
    parts = [part.strip() for part in value.split("+") if part.strip()]
    if not parts:
        return ""
    if len(set(parts)) == 1:
        return parts[0]
    return "+".join(sorted(parts))


def identity_components(value: str) -> tuple[str, ...]:
    return tuple(part.strip() for part in str(value or "").split("+")
                 if part.strip())


def upstream_truthy(value: object) -> bool:
    return str(value or "").strip().lower() in {
        "1", "true", "t", "yes", "y"}


def read_identity_sequence(path: str, label: str) -> list[str]:
    """Reproduce upstream identity-line parsing and first-seen deduplication."""
    values: list[str] = []
    seen: set[str] = set()
    with open(path, "r", encoding="utf-8") as handle:
        for raw in handle:
            line = raw.split("#", 1)[0].strip()
            if not line:
                continue
            value = line.split("\t", 1)[0].split(None, 1)[0].strip()
            if value and value not in seen:
                seen.add(value)
                values.append(value)
    if not values:
        raise ValueError(f"{label} is empty: {path}")
    return values


def selected_reconciliation_candidate(
        row: Mapping[str, str], candidate_set: str) -> bool:
    """Exact candidate selection used by upstream comparison-plan V4."""
    if upstream_truthy(row.get("reassignment_applied")):
        return True
    if candidate_set == "applied":
        return False
    if candidate_set != "exploratory":
        raise ValueError(
            f"unknown reconciliation ambient candidate set: {candidate_set}")
    proposed = str(row.get("proposed_donor_genotype", "")).strip()
    original = str(row.get("original_demux_assignment", "")).strip()
    current = (
        str(row.get("current_refined_assignment", "")).strip() or original)
    if (not proposed or model_identity(proposed) in {
            model_identity(original), model_identity(current)}):
        return False
    if str(row.get("proposed_droplet_state", "")).strip() != "SINGLE_CELL":
        return False
    if str(row.get("proposed_biological_admissibility", "")).strip() not in {
            "BIOLOGICAL_SINGLE_CELL_ALLOWED", "SINGLET_IDENTITY_CANDIDATE"}:
        return False
    confidence = str(row.get("decision_confidence", "")).strip()
    action = str(row.get("final_action", "")).strip()
    return (confidence in {"DECISIVE", "STRONG_NOT_AUTOAPPLIED"} and
            action not in {"", "KEEP"})


@dataclass
class AmbientPlanExpectation:
    demux_order: list[str]
    demux_identities: Mapping[str, str]
    validated_reconciled: Mapping[str, str]
    selected_rows: dict[str, dict[str, str]]
    selected_source_identities: set[str]
    selected_candidates: set[str]
    candidate_keys: set[str]
    original_receivers: list[str]
    original_donors: set[str]
    final_receivers: list[str]
    augmented_receivers: list[str]
    augmented_donors: set[str]
    displaced: dict[str, str]
    replacement_receivers: list[str]
    replacement_donors: set[str]
    removable_displaced_components: set[str]

    def reconciled_identity(self, barcode: str) -> str:
        return self.validated_reconciled.get(
            barcode, self.demux_identities[barcode])

    def is_scrutinized(self, barcode: str) -> bool:
        return (
            barcode in self.selected_rows or
            model_identity(self.demux_identities[barcode]) in
            self.selected_source_identities or
            model_identity(self.reconciled_identity(barcode)) in
            self.candidate_keys)


def reconstruct_ambient_plan_expectation(
        decisions_path: str, candidate_set: str, demux_order: list[str],
        demux_identities: Mapping[str, str],
        validated_reconciled: Mapping[str, str],
        original_receiver_path: str, original_donor_path: str,
        label: str,
) -> AmbientPlanExpectation:
    """Recompute the V4 plan from its hashed DEMUX and decision inputs."""
    if not set(validated_reconciled) <= set(demux_identities):
        raise ValueError(
            f"{label} validated identities extend beyond the DEMUX universe")
    required_fields = {
        "barcode", "original_demux_assignment", "proposed_donor_genotype",
        "reconciled_donor_genotype", "proposed_droplet_state",
        "proposed_biological_admissibility", "final_action",
        "decision_confidence", "reassignment_applied",
        "explicit_multiplet_evidence", "occupancy_resolution_status",
        "library_exchange_evidence_eligible", "event_id",
    }
    selected_rows: dict[str, dict[str, str]] = {}
    selected_source_identities: set[str] = set()
    selected_candidates: set[str] = set()
    missing_validated = set(validated_reconciled)
    seen_decisions: set[str] = set()
    opener = gzip.open if decisions_path.endswith(".gz") else open
    with opener(decisions_path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        header = list(reader.fieldnames or [])
        if (len(header) != len(set(header)) or
                not required_fields <= set(header)):
            raise ValueError(
                f"{label} decision table has an incompatible schema: "
                f"{decisions_path}")
        for line_number, source_row in enumerate(reader, start=2):
            if (None in source_row or
                    any(value is None for value in source_row.values())):
                raise ValueError(
                    f"{label} decision table has a malformed row at line "
                    f"{line_number}")
            row = {key: str(value).strip()
                   for key, value in source_row.items()}
            barcode = row["barcode"]
            if (not barcode or barcode in seen_decisions or
                    barcode not in demux_identities):
                raise ValueError(
                    f"{label} decision table has an empty, duplicate, or "
                    f"non-DEMUX barcode at line {line_number}")
            seen_decisions.add(barcode)
            recorded_original = row["original_demux_assignment"]
            if (not recorded_original or
                    model_identity(recorded_original) !=
                    model_identity(demux_identities[barcode])):
                raise ValueError(
                    f"{label}/{barcode} decision original identity does not "
                    "match DEMUX")
            if barcode in validated_reconciled:
                missing_validated.discard(barcode)
                if (model_identity(row["reconciled_donor_genotype"]) !=
                        model_identity(validated_reconciled[barcode])):
                    raise ValueError(
                        f"{label}/{barcode} decision final identity does not "
                        "match the validated assignment export")
            if (upstream_truthy(row["reassignment_applied"]) and
                    barcode not in validated_reconciled):
                raise ValueError(
                    f"{label}/{barcode} applied decision is absent from the "
                    "validated assignment export")
            if not selected_reconciliation_candidate(row, candidate_set):
                continue
            selected_rows[barcode] = row
            selected_source_identities.add(model_identity(
                row["original_demux_assignment"] or
                demux_identities[barcode]))
            reconciled_candidate = (
                row["reconciled_donor_genotype"] or
                validated_reconciled.get(barcode, demux_identities[barcode]))
            candidate_identities = (
                (reconciled_candidate,)
                if candidate_set == "applied" else
                (row["proposed_donor_genotype"], reconciled_candidate))
            selected_candidates.update(
                identity for identity in candidate_identities if identity)
    if missing_validated:
        raise ValueError(
            f"{label} validated assignments lack decision rows: "
            + ",".join(sorted(missing_validated)[:10]))

    original_receivers = read_identity_sequence(
        original_receiver_path, f"{label} original receiver roster")
    original_donors = set(read_identity_sequence(
        original_donor_path, f"{label} original ambient donor roster"))
    final_receivers: list[str] = []
    seen_final_receivers: set[str] = set()
    for barcode in demux_order:
        identity = validated_reconciled.get(
            barcode, demux_identities[barcode])
        if identity not in seen_final_receivers:
            seen_final_receivers.add(identity)
            final_receivers.append(identity)
    original_receiver_keys = {
        model_identity(identity) for identity in original_receivers}
    selected_candidates.update(
        identity for identity in final_receivers
        if model_identity(identity) not in original_receiver_keys)
    candidate_keys = {
        model_identity(identity) for identity in selected_candidates}

    augmented_receivers = list(original_receivers)
    augmented_receiver_keys = set(original_receiver_keys)
    for identity in final_receivers:
        key = model_identity(identity)
        if key not in augmented_receiver_keys:
            augmented_receivers.append(identity)
            augmented_receiver_keys.add(key)
    augmented_donors = set(original_donors)
    for identity in itertools.chain(final_receivers, selected_candidates):
        augmented_donors.update(identity_components(identity))

    applied_added_donors: set[str] = set()
    for row in selected_rows.values():
        if upstream_truthy(row["reassignment_applied"]):
            applied_added_donors.update(identity_components(
                row["reconciled_donor_genotype"]))
    final_receiver_donors: set[str] = set()
    for identity in final_receivers:
        final_receiver_donors.update(identity_components(identity))
    if candidate_set == "applied":
        ungrounded = sorted(
            donor for donor in augmented_donors - original_donors
            if donor not in applied_added_donors and
            donor not in final_receiver_donors)
        if ungrounded:
            raise ValueError(
                f"{label} applied roster has ungrounded donors: "
                + ",".join(ungrounded))

    displaced: dict[str, str] = {}
    for source_key in sorted(selected_source_identities):
        source_barcodes = [
            barcode for barcode in demux_order
            if model_identity(demux_identities[barcode]) == source_key]
        replacements = {
            model_identity(validated_reconciled.get(
                barcode, demux_identities[barcode]))
            for barcode in source_barcodes}
        if (source_barcodes and source_key not in replacements and
                len(replacements) == 1):
            displaced[source_key] = next(iter(replacements))
    replacement_receivers = [
        identity for identity in augmented_receivers
        if model_identity(identity) not in displaced]
    retained_receiver_donors: set[str] = set()
    for identity in replacement_receivers:
        retained_receiver_donors.update(identity_components(identity))
    displaced_components: set[str] = set()
    for identity in displaced:
        displaced_components.update(identity_components(identity))
    removable = displaced_components - retained_receiver_donors
    replacement_donors = augmented_donors - removable
    return AmbientPlanExpectation(
        demux_order=demux_order,
        demux_identities=demux_identities,
        validated_reconciled=validated_reconciled,
        selected_rows=selected_rows,
        selected_source_identities=selected_source_identities,
        selected_candidates=selected_candidates,
        candidate_keys=candidate_keys,
        original_receivers=original_receivers,
        original_donors=original_donors,
        final_receivers=final_receivers,
        augmented_receivers=augmented_receivers,
        augmented_donors=augmented_donors,
        displaced=displaced,
        replacement_receivers=replacement_receivers,
        replacement_donors=replacement_donors,
        removable_displaced_components=removable,
    )


def read_exact_generated_roster(path: str, label: str) -> list[str]:
    values: list[str] = []
    with open(path, "r", encoding="utf-8", newline="") as handle:
        for line_number, raw in enumerate(handle, start=1):
            value = raw.rstrip("\r\n")
            if (not value or value != value.strip() or
                    any(character.isspace() for character in value)):
                raise ValueError(
                    f"{label} has a malformed identity at line {line_number}: "
                    f"{path}")
            values.append(value)
    if not values:
        raise ValueError(f"{label} is empty: {path}")
    return values


def validate_ambient_plan_artifacts(
        expectation: AmbientPlanExpectation,
        context: Mapping[str, object],
        plan_paths: Mapping[str, str], sample_path: str,
        candidate_set: str, label: str,
) -> dict[str, int]:
    """Validate every deterministic V4 roster/provenance/context product."""
    replacement_eligible = bool(
        expectation.removable_displaced_components)
    expected_rosters: dict[str, list[str]] = {
        "augmented_receiver_lines": expectation.augmented_receivers,
        "augmented_ambient_candidates": sorted(
            expectation.augmented_donors),
    }
    if replacement_eligible:
        expected_rosters.update({
            "replacement_receiver_lines": expectation.replacement_receivers,
            "replacement_ambient_candidates": sorted(
                expectation.replacement_donors),
        })
    for name, expected in expected_rosters.items():
        path = plan_paths.get(name, "")
        observed = read_exact_generated_roster(
            path, f"{label} {name}")
        if observed != expected:
            raise ValueError(
                f"{label} {name} does not reproduce from the hashed "
                "reconciliation inputs")

    samples = read_identity_sequence(sample_path, f"{label} DEMUX samples")
    unknown_donors = sorted(expectation.augmented_donors - set(samples))
    if unknown_donors:
        raise ValueError(
            f"{label} augmented donors are absent from the DEMUX/main-panel "
            "sample universe: " + ",".join(unknown_donors))

    selected_context = context.get("selected_candidate_identities")
    if (not isinstance(selected_context, list) or
            any(not isinstance(value, str) or not value
                for value in selected_context) or
            len(selected_context) != len(set(selected_context)) or
            set(selected_context) != expectation.selected_candidates or
            [model_identity(value) for value in selected_context] !=
            sorted(model_identity(value) for value in selected_context)):
        raise ValueError(
            f"{label} selected_candidate_identities does not reproduce from "
            "the reconciliation decisions")
    expected_context_values: dict[str, object] = {
        "n_original_receivers": len(expectation.original_receivers),
        "n_original_candidates": len(expectation.original_donors),
        "n_augmented_receivers": len(expectation.augmented_receivers),
        "n_augmented_candidates": len(expectation.augmented_donors),
        "n_replacement_receivers": (
            len(expectation.replacement_receivers)
            if replacement_eligible else 0),
        "n_replacement_candidates": (
            len(expectation.replacement_donors)
            if replacement_eligible else 0),
        "replacement_arm_eligible": replacement_eligible,
        "replacement_arm_skip_reason": (
            "" if replacement_eligible else
            "no fully displaced source changed the ambient donor roster"),
        "added_candidate_components": sorted(
            expectation.augmented_donors - expectation.original_donors),
        "fully_displaced_identity_map": expectation.displaced,
        "removed_candidate_components": sorted(
            expectation.removable_displaced_components),
    }
    for name, expected in expected_context_values.items():
        observed = context.get(name)
        same_type = (
            isinstance(observed, bool) if isinstance(expected, bool) else
            type(observed) is int if type(expected) is int else
            isinstance(observed, type(expected)))
        if not same_type or observed != expected:
            raise ValueError(
                f"{label} context {name} does not reproduce from the hashed "
                "reconciliation inputs")

    evidence: dict[str, dict[str, set[str]]] = {}
    for barcode, row in expectation.selected_rows.items():
        fields = (
            ("reconciled_donor_genotype",)
            if candidate_set == "applied" else
            ("proposed_donor_genotype", "reconciled_donor_genotype"))
        identities = {
            row[field] for field in fields if row.get(field, "")}
        components: set[str] = set()
        for identity in identities:
            components.update(identity_components(identity))
        for donor in components:
            donor_evidence = evidence.setdefault(donor, {
                "barcodes": set(), "applied": set(), "unapplied": set(),
                "confidences": set(), "events": set(), "identities": set(),
            })
            donor_evidence["barcodes"].add(barcode)
            target = ("applied" if upstream_truthy(
                row.get("reassignment_applied")) else "unapplied")
            donor_evidence[target].add(barcode)
            if row.get("decision_confidence", ""):
                donor_evidence["confidences"].add(
                    row["decision_confidence"])
            if row.get("event_id", ""):
                donor_evidence["events"].add(row["event_id"])
            donor_evidence["identities"].update(identities)
    final_receiver_counts = {
        donor: 0 for donor in expectation.augmented_donors}
    for barcode in expectation.demux_order:
        for donor in identity_components(
                expectation.reconciled_identity(barcode)):
            if donor in final_receiver_counts:
                final_receiver_counts[donor] += 1

    provenance_path = plan_paths["ambient_candidate_provenance"]
    provenance_header = (
        "donor", "roster_status", "inclusion_reason", "selected_rows",
        "applied_rows", "unapplied_exploratory_rows",
        "final_receiver_cells", "decision_confidences", "event_ids",
        "source_identities")
    observed_rows = strict_posthoc_tsv(
        provenance_path, provenance_header,
        f"{label} ambient candidate provenance")
    if len(observed_rows) != len(expectation.augmented_donors):
        raise ValueError(
            f"{label} ambient candidate provenance donor count is wrong")
    for donor, observed in zip(
            sorted(expectation.augmented_donors), observed_rows):
        donor_evidence = evidence.get(donor, {})
        applied_count = len(donor_evidence.get("applied", ()))
        unapplied_count = len(donor_evidence.get("unapplied", ()))
        final_count = final_receiver_counts[donor]
        if donor in expectation.original_donors:
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
        expected = {
            "donor": donor,
            "roster_status": roster_status,
            "inclusion_reason": inclusion_reason,
            "selected_rows": str(len(donor_evidence.get("barcodes", ()))),
            "applied_rows": str(applied_count),
            "unapplied_exploratory_rows": str(unapplied_count),
            "final_receiver_cells": str(final_count),
            "decision_confidences": ";".join(sorted(
                donor_evidence.get("confidences", ()))) or "NA",
            "event_ids": ";".join(sorted(
                donor_evidence.get("events", ()))) or "NA",
            "source_identities": ";".join(sorted(
                donor_evidence.get("identities", ()))) or "NA",
        }
        if observed != expected:
            raise ValueError(
                f"{label} ambient candidate provenance does not reproduce "
                f"for donor {donor}")
    return {
        "comparison_barcodes": len(expectation.demux_identities),
        "validated_barcodes": len(expectation.validated_reconciled),
        "preserved_barcodes": (
            len(expectation.demux_identities) -
            len(expectation.validated_reconciled)),
    }


def gzip_payload_sha256(path: str) -> str:
    digest = hashlib.sha256()
    with gzip.open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def identity_roster(path: str, label: str) -> set[str]:
    values: set[str] = set()
    with open(path, "r", encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, start=1):
            value = raw.split("#", 1)[0].strip()
            if not value:
                continue
            value = value.split("\t", 1)[0].split(None, 1)[0].strip()
            if not value or value in values:
                raise ValueError(
                    f"{label} has an empty/duplicate identity at line "
                    f"{line_number}: {path}")
            values.add(value)
    if not values:
        raise ValueError(f"{label} is empty: {path}")
    return values


def identity_component_roster(path: str, label: str) -> set[str]:
    """Expand receiver identities such as A+B into ambient donor components."""
    values: set[str] = set()
    with open(path, "r", encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, start=1):
            value = raw.split("#", 1)[0].strip()
            if not value:
                continue
            value = value.split("\t", 1)[0].split(None, 1)[0].strip()
            components = [part.strip() for part in value.split("+")
                          if part.strip()]
            if not components:
                raise ValueError(
                    f"{label} has an empty identity at line {line_number}: "
                    f"{path}")
            values.update(components)
    if not values:
        raise ValueError(f"{label} is empty: {path}")
    return values


@dataclass
class ArmProfileSummary:
    """Bounded sufficient statistics from one direct cell/source profile."""

    cell_count: int
    profile_rows: int
    source_labels: set[str]
    population_counts: dict[str, int]
    population_rate_sums: dict[str, float]
    source_mass_sums: dict[tuple[str, str], float]
    source_burden_sums: dict[tuple[str, str], float]


class AmbientValidationStore:
    """Disk-backed exact per-cell index used only during validation."""

    def __init__(self, _run_root: str) -> None:
        # SQLite performs random I/O, so keep this transient validation index
        # on node-local scratch rather than on the shared run filesystem.
        local_tmp_root = (
            os.environ.get("SLURM_TMPDIR") or
            os.environ.get("TMPDIR") or
            tempfile.gettempdir())
        fd, path = tempfile.mkstemp(
            prefix="ambient_validation_", suffix=".sqlite3",
            dir=local_tmp_root)
        os.close(fd)
        self.path = path
        self.connection = sqlite3.connect(path)
        self.connection.execute("PRAGMA journal_mode=OFF")
        self.connection.execute("PRAGMA synchronous=OFF")
        self.connection.execute("PRAGMA temp_store=MEMORY")
        self.connection.execute("PRAGMA cache_size=-131072")
        self.connection.executescript("""
            CREATE TABLE cell_rate (
                library INTEGER NOT NULL,
                arm TEXT NOT NULL,
                barcode TEXT NOT NULL,
                identity TEXT NOT NULL,
                rate REAL NOT NULL,
                rate_text TEXT NOT NULL,
                se_text TEXT NOT NULL,
                PRIMARY KEY (library, arm, barcode)
            ) WITHOUT ROWID;
            CREATE TABLE profile_seen (
                library INTEGER NOT NULL,
                arm TEXT NOT NULL,
                barcode TEXT NOT NULL,
                PRIMARY KEY (library, arm, barcode)
            ) WITHOUT ROWID;
            CREATE TABLE cell_stratum (
                library INTEGER NOT NULL,
                barcode TEXT NOT NULL,
                stratum TEXT NOT NULL,
                PRIMARY KEY (library, barcode)
            ) WITHOUT ROWID;
            CREATE TABLE contrast_cell (
                library INTEGER NOT NULL,
                contrast TEXT NOT NULL,
                barcode TEXT NOT NULL,
                stratum TEXT NOT NULL,
                left_text TEXT NOT NULL,
                right_text TEXT NOT NULL,
                delta_text TEXT NOT NULL,
                delta REAL NOT NULL,
                abs_delta REAL NOT NULL,
                PRIMARY KEY (library, contrast, barcode)
            ) WITHOUT ROWID;
            CREATE INDEX contrast_delta_full
                ON contrast_cell (library, contrast, delta);
            CREATE INDEX contrast_abs_delta_full
                ON contrast_cell (library, contrast, abs_delta);
            CREATE INDEX contrast_delta_stratum
                ON contrast_cell (library, contrast, stratum, delta);
            CREATE INDEX contrast_abs_delta_stratum
                ON contrast_cell (library, contrast, stratum, abs_delta);
        """)
        self._closed = False

    def close(self) -> None:
        if self._closed:
            return
        self._closed = True
        try:
            self.connection.close()
        finally:
            for suffix in ("", "-journal", "-wal", "-shm"):
                try:
                    os.unlink(self.path + suffix)
                except FileNotFoundError:
                    pass

    def add_scrutiny(
            self, library: int, path: str, label: str,
            expectation: AmbientPlanExpectation,
    ) -> dict[str, int]:
        expected_header = (
            "barcode", "demux_identity", "reconciled_identity",
            "proposed_identity", "selected_candidate_row",
            "reassignment_applied", "changed", "scrutinized", "stratum",
            "transition", "event_id", "decision_confidence")
        count = 0
        changed_count = 0
        scrutinized_count = 0
        background_count = 0
        expected_barcodes = iter(sorted(expectation.demux_identities))
        with open(path, "r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            header = list(reader.fieldnames or [])
            if tuple(header) != expected_header:
                raise ValueError(
                    f"{label} has a noncanonical header: {path}")
            for line_number, row in enumerate(reader, start=2):
                if (None in row or
                        any(value is None for value in row.values())):
                    raise ValueError(
                        f"{label} has a malformed row at line {line_number}")
                if any(value != value.strip() for value in row.values()):
                    raise ValueError(
                        f"{label} has surrounding whitespace at line "
                        f"{line_number}")
                try:
                    expected_barcode = next(expected_barcodes)
                except StopIteration as exc:
                    raise ValueError(
                        f"{label} has more rows than the exact DEMUX "
                        "universe") from exc
                barcode = row["barcode"]
                if barcode != expected_barcode:
                    raise ValueError(
                        f"{label} does not preserve the upstream sorted "
                        f"barcode universe at line {line_number}")
                demux_identity = expectation.demux_identities[barcode]
                reconciled_identity = expectation.reconciled_identity(barcode)
                selected_row = expectation.selected_rows.get(barcode)
                selected = selected_row is not None
                changed = (
                    model_identity(demux_identity) !=
                    model_identity(reconciled_identity))
                scrutinized = expectation.is_scrutinized(barcode)
                stratum = (
                    "changed_target" if changed else
                    "scrutinized_other" if scrutinized else
                    "background")
                expected_row = {
                    "barcode": barcode,
                    "demux_identity": demux_identity or "NA",
                    "reconciled_identity": reconciled_identity or "NA",
                    "proposed_identity": (
                        selected_row.get("proposed_donor_genotype", "")
                        if selected_row else "") or "NA",
                    "selected_candidate_row": str(int(selected)),
                    "reassignment_applied": str(int(upstream_truthy(
                        selected_row.get("reassignment_applied", "")
                        if selected_row else ""))),
                    "changed": str(int(changed)),
                    "scrutinized": str(int(scrutinized)),
                    "stratum": stratum,
                    "transition": (
                        f"{demux_identity or 'NA'} -> "
                        f"{reconciled_identity or 'NA'}"),
                    "event_id": (
                        selected_row.get("event_id", "")
                        if selected_row else "") or "NA",
                    "decision_confidence": (
                        selected_row.get("decision_confidence", "")
                        if selected_row else "") or "NA",
                }
                if any(row[name] != expected_row[name]
                       for name in expected_header):
                    raise ValueError(
                        f"{label} row for {barcode} does not reproduce from "
                        "the hashed DEMUX/reconciliation decision inputs")
                try:
                    self.connection.execute(
                        "INSERT INTO cell_stratum VALUES (?, ?, ?)",
                        (library, barcode, stratum))
                except sqlite3.IntegrityError as exc:
                    raise ValueError(
                        f"{label} has duplicate barcode {barcode!r}") from exc
                count += 1
                changed_count += int(changed)
                scrutinized_count += int(scrutinized)
                background_count += int(stratum == "background")
        if count < 1:
            raise ValueError(f"{label} has no rows: {path}")
        try:
            next(expected_barcodes)
        except StopIteration:
            pass
        else:
            raise ValueError(
                f"{label} has fewer rows than the exact DEMUX universe")
        return {
            "rows": count,
            "changed": changed_count,
            "scrutinized": scrutinized_count,
            "background": background_count,
        }

    def stratum(self, library: int, barcode: str) -> str:
        row = self.connection.execute(
            "SELECT stratum FROM cell_stratum WHERE library=? AND barcode=?",
            (library, barcode)).fetchone()
        return str(row[0]) if row is not None else "unmapped"

    def add_rate(
            self, library: int, arm: str, barcode: str, identity: str,
            rate: float, rate_text: str, se_text: str, label: str) -> None:
        try:
            self.connection.execute(
                "INSERT INTO cell_rate VALUES (?, ?, ?, ?, ?, ?, ?)",
                (library, arm, barcode, identity, rate, rate_text, se_text))
        except sqlite3.IntegrityError as exc:
            raise ValueError(
                f"{label} has duplicate direct rate for {barcode}") from exc

    def rate(
            self, library: int, arm: str, barcode: str
    ) -> tuple[float, str, str, str] | None:
        row = self.connection.execute(
            "SELECT rate, rate_text, se_text, identity FROM cell_rate "
            "WHERE library=? AND arm=? AND barcode=?",
            (library, arm, barcode)).fetchone()
        if row is None:
            return None
        return float(row[0]), str(row[1]), str(row[2]), str(row[3])

    def rate_count(self, library: int, arm: str) -> int:
        row = self.connection.execute(
            "SELECT COUNT(*) FROM cell_rate WHERE library=? AND arm=?",
            (library, arm)).fetchone()
        return int(row[0])

    def scrutiny_count(self, library: int) -> int:
        row = self.connection.execute(
            "SELECT COUNT(*) FROM cell_stratum WHERE library=?",
            (library,)).fetchone()
        return int(row[0])

    def same_scrutiny_rate_universe(self, library: int, arm: str) -> bool:
        scrutiny_only = self.connection.execute("""
            SELECT EXISTS(
                SELECT barcode FROM cell_stratum WHERE library=?
                EXCEPT
                SELECT barcode FROM cell_rate WHERE library=? AND arm=?
            )
        """, (library, library, arm)).fetchone()[0]
        rate_only = self.connection.execute("""
            SELECT EXISTS(
                SELECT barcode FROM cell_rate WHERE library=? AND arm=?
                EXCEPT
                SELECT barcode FROM cell_stratum WHERE library=?
            )
        """, (library, arm, library)).fetchone()[0]
        return not scrutiny_only and not rate_only

    def arm_names(self, library: int) -> set[str]:
        return {str(row[0]) for row in self.connection.execute(
            "SELECT DISTINCT arm FROM cell_rate WHERE library=? AND arm!='S'",
            (library,))}

    def same_rate_universe(self, library: int, left: str, right: str) -> bool:
        query = """
            SELECT EXISTS(
                SELECT barcode FROM cell_rate WHERE library=? AND arm=?
                EXCEPT
                SELECT barcode FROM cell_rate WHERE library=? AND arm=?
            )
        """
        left_only = self.connection.execute(
            query, (library, left, library, right)).fetchone()[0]
        right_only = self.connection.execute(
            query, (library, right, library, left)).fetchone()[0]
        return not left_only and not right_only

    def add_profile_seen(
            self, library: int, arm: str, barcode: str, label: str) -> None:
        try:
            self.connection.execute(
                "INSERT INTO profile_seen VALUES (?, ?, ?)",
                (library, arm, barcode))
        except sqlite3.IntegrityError as exc:
            raise ValueError(
                f"{label} profile barcode block is repeated for {barcode}") \
                from exc

    def profile_count(self, library: int, arm: str) -> int:
        row = self.connection.execute(
            "SELECT COUNT(*) FROM profile_seen WHERE library=? AND arm=?",
            (library, arm)).fetchone()
        return int(row[0])

    def add_contrast(
            self, library: int, contrast: str, barcode: str, stratum: str,
            left_text: str, right_text: str, delta_text: str,
            delta: float, label: str) -> None:
        try:
            self.connection.execute(
                "INSERT INTO contrast_cell VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)",
                (library, contrast, barcode, stratum, left_text, right_text,
                 delta_text, delta, abs(delta)))
        except sqlite3.IntegrityError as exc:
            raise ValueError(
                f"{label} has duplicate {contrast}/{barcode} row") from exc

    def contrast_count(self, library: int, contrast: str) -> int:
        row = self.connection.execute(
            "SELECT COUNT(*) FROM contrast_cell "
            "WHERE library=? AND contrast=?", (library, contrast)).fetchone()
        return int(row[0])

    def contrast_summary(
            self, library: int, contrast: str, population: str
    ) -> tuple[int, float, float, float, float]:
        where = "library=? AND contrast=?"
        parameters: list[object] = [library, contrast]
        if population != "full_library":
            where += " AND stratum=?"
            parameters.append(population)
        row = self.connection.execute(
            f"SELECT COUNT(*), SUM(delta), "
            f"SUM(CASE WHEN delta < 0 THEN 1 ELSE 0 END) "
            f"FROM contrast_cell WHERE {where}", parameters).fetchone()
        count = int(row[0])
        if count == 0:
            return 0, math.nan, math.nan, math.nan, math.nan

        def middle(column: str) -> float:
            take = 1 if count % 2 else 2
            offset = (count - 1) // 2
            values = [float(value[0]) for value in self.connection.execute(
                f"SELECT {column} FROM contrast_cell WHERE {where} "
                f"ORDER BY {column} LIMIT ? OFFSET ?",
                [*parameters, take, offset])]
            return sum(values) / len(values)

        mean = float(row[1]) / count
        fraction = float(row[2]) / count
        return count, mean, middle("delta"), middle("abs_delta"), fraction

    def contrast_rows(
            self, library: int, barcode: str
    ) -> dict[str, tuple[str, str, str, str]]:
        return {
            str(row[0]): (str(row[1]), str(row[2]), str(row[3]), str(row[4]))
            for row in self.connection.execute(
                "SELECT contrast, stratum, left_text, right_text, delta_text "
                "FROM contrast_cell WHERE library=? AND barcode=?",
                (library, barcode))
        }


class FourArmFinalizerLookup:
    """Bounded lookup that reconstructs one final-ledger row on demand."""

    def __init__(
            self, store: AmbientValidationStore,
            required_contrasts: Mapping[int, set[str]],
            burden_text: Mapping[tuple[int, str], str],
            background_text: Mapping[int, str]) -> None:
        self.store = store
        self.required_contrasts = required_contrasts
        self.burden_text = burden_text
        self.background_text = background_text

    def get(self, library: int, barcode: str) -> dict[str, object] | None:
        required = self.required_contrasts.get(library)
        if required is None:
            return None
        rows = self.store.contrast_rows(library, barcode)
        if set(rows) != required:
            return None
        roster = rows["roster_effect_demux"]
        assignment = rows["assignment_effect_augmented"]
        combined = rows["combined_production_change"]
        replacement = rows.get("replacement_sensitivity")
        stratum = roster[0]
        if any(row[0] != stratum for row in rows.values()):
            return None
        burden = self.burden_text.get((library, stratum))
        if burden is None:
            burden = self.burden_text.get((library, "full_library"), "NA")
        return {
            "ambient_arm_a_c": roster[1],
            "ambient_arm_b_c": roster[2],
            "ambient_arm_c_c": assignment[2],
            "ambient_arm_d_c": replacement[2] if replacement else None,
            "ambient_roster_effect_b_minus_a": roster[3],
            "ambient_assignment_effect_c_minus_b": assignment[3],
            "ambient_replacement_effect_d_minus_c": (
                replacement[3] if replacement else None),
            "ambient_combined_augmented_c_minus_a": combined[3],
            "ambient_production_arm": "C",
            "ambient_production_c": assignment[2],
            "ambient_production_minus_original_c": combined[3],
            "ambient_evaluation_status": (
                "PAIRED_A_B_C_D" if replacement else
                "PAIRED_A_B_C_D_NOT_APPLICABLE"),
            "ambient_exact_donor_burden_fields": burden,
            "ambient_background_shift_fields": self.background_text[library],
        }


def validate_geometry_selected_bundle(
        prefix: str, expected_barcodes: set[str], label: str,
        expected_identities: Mapping[str, str] | None = None,
        expected_sources: set[str] | None = None,
        require_cell_profile: bool = True,
        *, validation_store: AmbientValidationStore,
        library: int, arm: str,
) -> tuple[dict[str, object], ArmProfileSummary | None]:
    """Stream-check the selected rate, diagnostic, audit, and profile ledgers."""
    rate_path = prefix + ".contam_rate"
    audit_fields = (
        "geometry_gate_version", "geometry_gate_triggered",
        "geometry_gate_parent_axis_alpha",
        "geometry_gate_ambient_orthogonal_norm",
        "geometry_gate_max_raw_excluded_parent_mass",
        "geometry_gate_parent_axis_alpha_threshold",
        "geometry_gate_ambient_orthogonal_norm_threshold",
        "geometry_gate_parent_mass_threshold",
        "source_exclusion_strength_base",
        "source_exclusion_strength_fallback",
        "source_exclusion_strength_selected",
        "geometry_gate_selection_reason", "base_c_selected",
        "fallback_c_selected", "selected_c", "geometry_gate_evaluable",
        "base_rate_status", "fallback_rate_status", "selected_endpoint",
        "selected_rate_status", "base_optimizer_status",
        "fallback_optimizer_status", "base_profile_validation_status",
        "fallback_profile_validation_status", "base_fit_failure_reason",
        "fallback_fit_failure_reason", "selected_c_se",
    )

    diagnostic_path = prefix + ".contam_diagnostics.tsv"
    audit_path = prefix + ".geometry_gate_audit.tsv"
    diagnostic_required = {"barcode", "identity", *audit_fields}
    row_count = 0
    sentinel = object()

    def rate_rows(handle: object) -> Iterable[tuple[int, str]]:
        for line_number, raw in enumerate(handle, start=1):
            if raw.strip():
                yield line_number, raw

    with (open(rate_path, "r", encoding="utf-8") as rate_handle,
          open(diagnostic_path, "r", encoding="utf-8", newline="")
          as diagnostic_handle,
          open(audit_path, "r", encoding="utf-8", newline="")
          as audit_handle):
        diagnostic_reader = csv.DictReader(
            diagnostic_handle, delimiter="\t")
        audit_reader = csv.DictReader(audit_handle, delimiter="\t")
        diagnostic_header = list(diagnostic_reader.fieldnames or [])
        audit_header = list(audit_reader.fieldnames or [])
        if (len(diagnostic_header) != len(set(diagnostic_header)) or
                not diagnostic_required <= set(diagnostic_header)):
            raise ValueError(
                f"{label} has an incompatible table: {diagnostic_path}")
        if (len(audit_header) != len(set(audit_header)) or
                not {"barcode", *audit_fields} <= set(audit_header)):
            raise ValueError(
                f"{label} has an incompatible table: {audit_path}")
        for rate_item, diagnostic, audit in itertools.zip_longest(
                rate_rows(rate_handle), diagnostic_reader, audit_reader,
                fillvalue=sentinel):
            if sentinel in (rate_item, diagnostic, audit):
                raise ValueError(
                    f"{label} selected rate/diagnostic/audit row counts differ")
            line_number, raw = rate_item
            rate_fields = raw.split()
            if (len(rate_fields) != 3 or not rate_fields[0] or
                    rate_fields[0] not in expected_barcodes):
                raise ValueError(
                    f"{label} selected rate row is malformed/duplicate at "
                    f"{rate_path}:{line_number}")
            try:
                rate = float(rate_fields[1])
                standard_error = float(rate_fields[2])
            except ValueError as exc:
                raise ValueError(
                    f"{label} has a nonnumeric selected rate at "
                    f"{rate_path}:{line_number}") from exc
            if (not math.isfinite(rate) or not 0 <= rate < 1 or
                    math.isinf(standard_error) or
                    (math.isfinite(standard_error) and standard_error < 0)):
                raise ValueError(
                    f"{label} has an invalid selected rate at "
                    f"{rate_path}:{line_number}")
            if (None in diagnostic or None in audit or
                    any(value is None for value in diagnostic.values()) or
                    any(value is None for value in audit.values())):
                raise ValueError(
                    f"{label} has a malformed selected endpoint row")
            barcode = rate_fields[0]
            if (barcode != diagnostic["barcode"] or
                    barcode != audit["barcode"] or
                    any(diagnostic[field] != audit[field]
                        for field in audit_fields) or
                    diagnostic["selected_c"] != rate_fields[1] or
                    diagnostic["selected_c_se"] != rate_fields[2] or
                    diagnostic["selected_rate_status"] != "finite" or
                    diagnostic["selected_endpoint"] not in {
                        "base", "fallback"} or
                    not diagnostic["geometry_gate_selection_reason"]):
                raise ValueError(
                    f"{label} selected endpoint evidence disagrees for "
                    f"{barcode}")
            if expected_identities is not None:
                expected = model_identity(
                    expected_identities.get(barcode, ""))
                if model_identity(diagnostic["identity"]) != expected:
                    raise ValueError(
                        f"{label} selected identity disagrees for {barcode}")
            validation_store.add_rate(
                library, arm, barcode, diagnostic["identity"], rate,
                rate_fields[1], rate_fields[2], label)
            row_count += 1
    if row_count != len(expected_barcodes):
        raise ValueError(f"{label} selected rate roster is incomplete")

    if not require_cell_profile:
        return ({"rows": row_count, "profile_rows": 0}, None)

    profile_path = prefix + ".cell_source_profile.tsv"
    profile_required = {
        "barcode", "identity", "source_label", "scoring_profile_mass",
        "scoring_profile_mass_sum", "scoring_profile_status"}
    profile_rows = 0
    profile_cells = 0
    source_labels: set[str] = set()
    population_counts: dict[str, int] = {}
    population_rate_sums: dict[str, float] = {}
    source_mass_sums: dict[tuple[str, str], float] = {}
    source_burden_sums: dict[tuple[str, str], float] = {}
    current_barcode = ""
    current_identity = ""
    current_masses: dict[str, float] = {}
    populations = {"changed_target", "scrutinized_other", "background"}

    def finish_profile_block() -> None:
        nonlocal profile_cells
        if not current_barcode:
            return
        direct = validation_store.rate(library, arm, current_barcode)
        if direct is None:
            raise ValueError(
                f"{label} profile barcode is absent from direct rates: "
                f"{current_barcode}")
        rate, _rate_text, _se_text, diagnostic_identity = direct
        if model_identity(current_identity) != model_identity(
                diagnostic_identity):
            raise ValueError(
                f"{label} profile/diagnostic identity differs for "
                f"{current_barcode}")
        observed_sources = set(current_masses)
        if (expected_sources is not None and
                observed_sources - {"other_species"} != expected_sources):
            raise ValueError(
                f"{label} cell/source roster differs for {current_barcode}")
        if not math.isclose(
                sum(current_masses.values()), 1.0,
                rel_tol=0.0, abs_tol=1e-6):
            raise ValueError(
                f"{label} cell/source profile is not a simplex for "
                f"{current_barcode}")
        validation_store.add_profile_seen(
            library, arm, current_barcode, label)
        stratum = validation_store.stratum(library, current_barcode)
        selected_populations = ["full_library"]
        if stratum in populations:
            selected_populations.append(stratum)
        for population in selected_populations:
            population_counts[population] = (
                population_counts.get(population, 0) + 1)
            population_rate_sums[population] = (
                population_rate_sums.get(population, 0.0) + rate)
            for source, mass in current_masses.items():
                key = (population, source)
                source_mass_sums[key] = source_mass_sums.get(key, 0.0) + mass
                source_burden_sums[key] = (
                    source_burden_sums.get(key, 0.0) + rate * mass)
        source_labels.update(observed_sources)
        profile_cells += 1

    with open(profile_path, "r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        header = list(reader.fieldnames or [])
        if (len(header) != len(set(header)) or
                not profile_required <= set(header)):
            raise ValueError(
                f"{label} has an incompatible table: {profile_path}")
        for row in reader:
            if (None in row or any(value is None for value in row.values())):
                raise ValueError(
                    f"{label} has a malformed cell/source profile row")
            barcode = row["barcode"]
            source = row["source_label"]
            expected_identity = (
                model_identity(expected_identities.get(barcode, ""))
                if expected_identities is not None else "")
            if (barcode not in expected_barcodes or not source or
                    (expected_identities is not None and
                     model_identity(row["identity"]) != expected_identity) or
                    not row["scoring_profile_status"]):
                raise ValueError(
                    f"{label} cell/source profile is inconsistent for "
                    f"{barcode}")
            if barcode != current_barcode:
                finish_profile_block()
                current_barcode = barcode
                current_identity = row["identity"]
                current_masses = {}
            elif model_identity(row["identity"]) != model_identity(
                    current_identity):
                raise ValueError(
                    f"{label} profile identity changes within {barcode}")
            if source in current_masses:
                raise ValueError(
                    f"{label} has duplicate profile source {source!r} for "
                    f"{barcode}")
            try:
                mass = float(row["scoring_profile_mass"])
                reported_sum = float(row["scoring_profile_mass_sum"])
            except ValueError as exc:
                raise ValueError(
                    f"{label} cell/source profile has nonnumeric mass") \
                    from exc
            if (not math.isfinite(mass) or mass < 0 or
                    not math.isfinite(reported_sum) or
                    not math.isclose(reported_sum, 1.0, rel_tol=0.0,
                                     abs_tol=1e-6)):
                raise ValueError(
                    f"{label} cell/source profile has invalid mass")
            current_masses[source] = mass
            profile_rows += 1
    finish_profile_block()
    if (profile_cells != row_count or
            validation_store.profile_count(library, arm) != row_count):
        raise ValueError(
            f"{label} cell/source profile does not cover every direct rate")
    summary = ArmProfileSummary(
        cell_count=profile_cells,
        profile_rows=profile_rows,
        source_labels=source_labels,
        population_counts=population_counts,
        population_rate_sums=population_rate_sums,
        source_mass_sums=source_mass_sums,
        source_burden_sums=source_burden_sums,
    )
    return ({"rows": row_count, "profile_rows": profile_rows}, summary)


def validate_ambient_profile(
        path: str, expected_donors: set[str], label: str) -> dict[str, object]:
    observed: set[str] = set()
    total = 0.0
    with open(path, "r", encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, start=1):
            if not raw.strip():
                continue
            fields = raw.split()
            if len(fields) not in {2, 3} or not fields[0] or fields[0] in observed:
                raise ValueError(
                    f"{label} has a malformed/duplicate row at line "
                    f"{line_number}: {path}")
            try:
                numeric = [float(value) for value in fields[1:]]
            except ValueError as exc:
                raise ValueError(
                    f"{label} has a nonnumeric row at line {line_number}: "
                    f"{path}") from exc
            if any(not math.isfinite(value) or value < 0 for value in numeric):
                raise ValueError(
                    f"{label} has an invalid value at line {line_number}: "
                    f"{path}")
            total += numeric[0]
            observed.add(fields[0])
    if observed != expected_donors or not math.isclose(
            total, 1.0, rel_tol=0.0, abs_tol=1e-3):
        raise ValueError(
            f"{label} donor universe/normalization mismatch: "
            f"expected={sorted(expected_donors)}, observed={sorted(observed)}, "
            f"total={total}")
    return {"rows": len(observed), "profile_total": total}


def validate_sized_path_record(
        record: object, expected_path: str, label: str,
        require_sha256: bool = False) -> None:
    if not isinstance(record, dict):
        raise ValueError(f"{label} path record is not an object")
    expected_path = absolute(expected_path)
    try:
        current_size = os.stat(expected_path).st_size
    except OSError as exc:
        raise ValueError(f"{label} cannot be statted: {expected_path}: {exc}") \
            from exc
    observed_size = record.get("size_bytes", record.get("size"))
    if (absolute(str(record.get("path", ""))) != expected_path or
            observed_size != current_size):
        raise ValueError(
            f"{label} path/size record no longer matches {expected_path}")
    if require_sha256 and record.get("sha256") != sha256_file(expected_path):
        raise ValueError(f"{label} SHA256 no longer matches {expected_path}")


def option_float_matches(
        argv: Sequence[str], option: str, expected: float) -> bool:
    raw = argv_value(argv, option)
    try:
        value = float(raw) if raw is not None else math.nan
    except ValueError:
        return False
    return math.isfinite(value) and math.isclose(
        value, expected, rel_tol=0.0, abs_tol=1e-12)


def validate_geometry_command(
        args: argparse.Namespace, script_path: str, output_prefix: str,
        receiver_lines: str, ambient_candidates: str, condf: str,
        frozen: bool, label: str) -> list[str]:
    with open(script_path, "r", encoding="utf-8") as handle:
        script_payload = handle.read()
    if "Skipping (use --force to rerun)" in script_payload:
        raise ValueError(f"{label} was not generated with --force")
    commands = logical_program_argvs(
        script_path, os.path.basename(args.geometry_gate_helper))
    expected = [
        args.geometry_gate_helper,
        "--estimator-binary", args.contam_binary,
        "--condition-key", args.ambient_condition,
        "--gate-version", "CK_GEOMETRY_GATE_V1",
        "--base-strength", "0",
        "--fallback-strength", "1",
        "--parent-axis-alpha-threshold", "0.80000000000000004",
        "--ambient-orthogonal-norm-threshold", "0.10000000000000001",
        "--parent-mass-threshold", "0.90000000000000002",
        "--", "-o", output_prefix,
        "--interindividual", "-X", receiver_lines,
        "--ambient_candidates", ambient_candidates,
        "--condf", condf,
        "--strict_condf", "--run_class", "production",
        "--assignments_basis", "library",
        "--expected_lines_basis", "library",
        "--ambient_candidates_basis", "library",
        "--condition_key", args.ambient_condition,
        "--run_contract", output_prefix + ".run_contract.json",
        "--profile_restarts", "8",
        "--r_feedback", "--candidate_keyed_rows",
        "--candidate_keyed_split", "0.05",
    ]
    if frozen:
        expected.append("--freeze_assignments")
    expected.extend(("-T", "${SLURM_CPUS_PER_TASK}"))
    matched = [argv for argv in commands if list(argv) == expected]
    if len(commands) != 1 or len(matched) != 1:
        raise ValueError(
            f"{label} does not contain the exact selected geometry-gated "
            "ambient command")
    return matched[0]


def validate_geometry_contract(
        path: str, output_prefix: str, condition: str,
        expected_inputs: Mapping[str, str], frozen: bool,
        label: str) -> dict[str, object]:
    try:
        with open(path, "r", encoding="utf-8") as handle:
            contract = json.load(handle)
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise ValueError(f"cannot load {label}: {path}: {exc}") from exc
    exact = {
        "contract_version": "geometry_gated_contam_estimate_run_contract_V3",
        "tool": "geometry_gated_contam_estimate",
        "condition_key": condition,
        "output_prefix": output_prefix,
        "run_class": "validation_candidate",
        "panel_mode": "interindividual",
        "geometry_gate_version": "CK_GEOMETRY_GATE_V1",
        "fixed_identity_ambient_comparison": frozen,
        "freeze_assignments": frozen,
        "run_once": False,
        "assignment_update_mode": (
            "iterative_frozen" if frozen else "iterative_reclassification"),
        "assignments_basis": "library",
        "expected_lines_basis": "library",
        "ambient_candidates_basis": "library",
        "fixed_ambient_enabled": False,
        "endpoint_input_contracts_validated": True,
        "production_contract_pass": True,
    }
    mismatches = {
        key: contract.get(key) for key, expected in exact.items()
        if contract.get(key) != expected}
    numeric = {
        "base_source_exclusion_strength": 0.0,
        "fallback_source_exclusion_strength": 1.0,
        "geometry_gate_parent_axis_alpha_threshold": 0.8,
        "geometry_gate_ambient_orthogonal_norm_threshold": 0.1,
        "geometry_gate_parent_mass_threshold": 0.9,
    }
    for key, expected in numeric.items():
        try:
            observed = float(contract.get(key))
        except (TypeError, ValueError):
            mismatches[key] = contract.get(key)
            continue
        if not math.isclose(observed, expected, rel_tol=0.0, abs_tol=1e-12):
            mismatches[key] = observed
    if mismatches:
        raise ValueError(f"{label} has contract mismatches: {mismatches}")
    for field, expected_path in expected_inputs.items():
        validate_sized_path_record(
            contract.get(field), expected_path, f"{label} {field}")
    return contract


def verify_four_arm_aggregate(
        args: argparse.Namespace, event_libraries: set[int],
        newest_arm_output_mtime: int,
        final_boundary_mtime: int,
        contexts: Mapping[int, Mapping[str, object]],
        validation_store: AmbientValidationStore,
        arm_summaries: Mapping[tuple[int, str], ArmProfileSummary],
        ) -> tuple[dict[str, object], FourArmFinalizerLookup]:
    """Validate the applied A/B/C/D tables consumed by the Phase-3 finalizer."""
    if not event_libraries:
        return ({"status": "NOT_APPLICABLE_ZERO_EVENT", "records": []},
                FourArmFinalizerLookup(validation_store, {}, {}, {}))
    root = os.path.join(
        args.upstream_analysis_root, "aggregate_library_analysis",
        "ambient_rna",
        f"{args.ambient_condition}__reconciliation-four-arm-"
        f"{args.ambient_candidate_set}")
    data = os.path.join(root, "data")
    schemas = {
        "reconciliation_arm_contracts.tsv": (
            "library", "condition", "series_key", "arm", "arm_key",
            "assignment_basis", "roster_basis", "candidate_set",
            "plan_fingerprint", "assignment_path", "receiver_lines",
            "ambient_candidates", "scrutiny_cells", "context_path",
            "assignment_update_mode", "assignment_score_basis",
            "n_expected_sources", "n_observed_sources", "missing_sources",
            "unexpected_sources", "profile_raw_total"),
        "reconciliation_planned_contrasts.tsv": (
            "library", "condition", "contrast", "left_arm", "right_arm",
            "left_arm_key", "right_arm_key", "population", "n_common",
            "mean_delta", "median_delta", "median_absolute_delta",
            "fraction_lower_in_right_arm"),
        "reconciliation_planned_contrast_cells.tsv": (
            "library", "condition", "contrast", "barcode", "stratum",
            "left_rate", "right_rate", "right_minus_left"),
        "reconciliation_exact_donor_burden.tsv": (
            "library", "condition", "series_key", "arm", "arm_key",
            "assignment_basis", "roster_basis", "population",
            "source_label", "n_cells", "mean_scoring_profile_mass",
            "mean_exact_contam_burden", "fraction_of_mean_contam_burden",
            "mean_total_contam_rate"),
        "reconciliation_assignment_switch_summary.tsv": (
            "library", "condition", "series_key", "arm", "arm_key",
            "assignment_basis", "n_input", "n_final", "n_union",
            "n_input_vs_planned_mismatch", "n_missing_planned_identity",
            "n_final_vs_input_switch", "frozen_assignment_invariant"),
        "reconciliation_diagnostics.tsv": (
            "library", "condition", "arm", "arm_key", "check", "status",
            "observed", "expected", "details"),
    }
    tables: dict[str, list[dict[str, str]]] = {}
    contrast_cells_path = os.path.join(
        data, "reconciliation_planned_contrast_cells.tsv")
    records: list[dict[str, object]] = []
    for name, header in schemas.items():
        path = os.path.join(data, name)
        if name != "reconciliation_planned_contrast_cells.tsv":
            rows = strict_posthoc_tsv(
                path, header, f"applied four-arm aggregate {name}")
            if not rows:
                raise ValueError(
                    f"applied four-arm aggregate table is empty: {path}")
            tables[name] = rows
        output_mtime = os.stat(path).st_mtime_ns
        if not newest_arm_output_mtime <= output_mtime <= final_boundary_mtime:
            raise ValueError(
                "applied four-arm aggregate is outside the arm-output-to-"
                f"final-ledger chronology: {path}")
        records.append(provenance_file_record(
            path, f"applied four-arm aggregate {name}", True))

    def selected_library(row: Mapping[str, str], label: str) -> int:
        library = normalize_library_number(row.get("library"))
        if library is None or library not in event_libraries:
            raise ValueError(
                f"{label} has an unexpected event-library value: "
                f"{row.get('library')!r}")
        if row.get("condition", "").strip() != args.ambient_condition:
            raise ValueError(
                f"{label} has an unexpected ambient condition for lib{library}")
        return library

    arm_contract_specs = {
        "demux_original": ("A", "demux", "original"),
        "demux_augmented": ("B", "demux", "augmented"),
        "reconciled_augmented": ("C", "reconciled", "augmented"),
        "reconciled_replacement": ("D", "reconciled", "replacement"),
    }
    required_arms_by_library = {
        library: ({"demux_original", "demux_augmented",
                   "reconciled_augmented", "reconciled_replacement"}
                  if contexts.get(library, {}).get(
                      "replacement_arm_eligible") is True else
                  {"demux_original", "demux_augmented",
                   "reconciled_augmented"})
        for library in event_libraries}
    contracts_by_library: dict[int, set[str]] = {
        library: set() for library in event_libraries}
    for row in tables["reconciliation_arm_contracts.tsv"]:
        library = selected_library(row, "four-arm contract")
        arm_key = row["arm_key"].strip()
        if arm_key in contracts_by_library[library]:
            raise ValueError(
                f"lib{library} has a duplicate four-arm contract for {arm_key}")
        contracts_by_library[library].add(arm_key)
        context = contexts.get(library, {})
        expected_basis = arm_contract_specs.get(arm_key)
        assignment_path = (
            context.get("demux_assignments", "")
            if arm_key in {"demux_original", "demux_augmented"} else
            context.get("reconciled_assignments", ""))
        receiver_path = {
            "demux_original": context.get("original_receiver_lines", ""),
            "demux_augmented": context.get("original_receiver_lines", ""),
            "reconciled_augmented": context.get(
                "augmented_receiver_lines", ""),
            "reconciled_replacement": context.get(
                "replacement_receiver_lines", ""),
        }.get(arm_key, "")
        candidate_path = {
            "demux_original": context.get(
                "original_ambient_candidates", ""),
            "demux_augmented": context.get(
                "augmented_ambient_candidates", ""),
            "reconciled_augmented": context.get(
                "augmented_ambient_candidates", ""),
            "reconciled_replacement": context.get(
                "replacement_ambient_candidates", ""),
        }.get(arm_key, "")
        if (expected_basis is None or
                (row["arm"].strip(), row["assignment_basis"].strip(),
                 row["roster_basis"].strip()) != expected_basis or
                row["series_key"].strip() !=
                f"{args.ambient_condition}__{arm_key}" or
                row["candidate_set"].strip() != args.ambient_candidate_set or
                row["plan_fingerprint"].strip() !=
                str(context.get("plan_fingerprint", "")) or
                absolute(row["assignment_path"]) !=
                absolute(str(assignment_path)) or
                absolute(row["receiver_lines"]) !=
                absolute(str(receiver_path)) or
                absolute(row["ambient_candidates"]) !=
                absolute(str(candidate_path)) or
                absolute(row["scrutiny_cells"]) !=
                absolute(str(context.get("scrutiny_cells", ""))) or
                absolute(row["context_path"]) != absolute(os.path.join(
                    os.path.dirname(str(context.get("scrutiny_cells", ""))),
                    f"lib{library}.comparison_context.json")) or
                row["assignment_update_mode"].strip() != "iterative_frozen" or
                row["assignment_score_basis"].strip() !=
                "original_demux_all_arms" or
                row["missing_sources"].strip().upper() not in {"NA", "NONE"} or
                row["unexpected_sources"].strip().upper() not in {"NA", "NONE"}):
            raise ValueError(
                f"lib{library} has a malformed applied four-arm contract row")
        try:
            expected_sources = int(row["n_expected_sources"])
            observed_sources = int(row["n_observed_sources"])
            profile_total = float(row["profile_raw_total"])
        except ValueError as exc:
            raise ValueError(
                f"lib{library} four-arm contract has nonnumeric accounting") \
                from exc
        if (expected_sources < 1 or observed_sources != expected_sources or
                not math.isfinite(profile_total) or profile_total <= 0):
            raise ValueError(
                f"lib{library} four-arm profile/source accounting is invalid")
    for library, observed in contracts_by_library.items():
        expected = required_arms_by_library[library]
        if observed != expected:
            raise ValueError(
                f"lib{library} applied four-arm contract set is wrong: "
                f"expected={sorted(expected)} observed={sorted(observed)}")

    contrast_specs = {
        "roster_effect_demux": ("A", "B", "demux_original",
                                "demux_augmented"),
        "assignment_effect_augmented": ("B", "C", "demux_augmented",
                                        "reconciled_augmented"),
        "combined_production_change": ("A", "C", "demux_original",
                                       "reconciled_augmented"),
        "replacement_sensitivity": ("C", "D", "reconciled_augmented",
                                    "reconciled_replacement"),
    }
    required_contrasts_by_library = {
        library: ({"roster_effect_demux", "assignment_effect_augmented",
                   "combined_production_change", "replacement_sensitivity"}
                  if "reconciled_replacement" in
                  required_arms_by_library[library] else
                  {"roster_effect_demux", "assignment_effect_augmented",
                   "combined_production_change"})
        for library in event_libraries}
    expected_populations = {
        "full_library", "changed_target", "scrutinized_other", "background"}
    summary_rows: dict[tuple[int, str, str], dict[str, str]] = {}
    for row in tables["reconciliation_planned_contrasts.tsv"]:
        library = selected_library(row, "four-arm contrast summary")
        contrast = row["contrast"].strip()
        population = row["population"].strip()
        if (contrast not in required_contrasts_by_library[library] or
                population not in expected_populations):
            raise ValueError(
                f"lib{library} has an unexpected four-arm contrast summary "
                f"row: {contrast}/{population}")
        left_arm, right_arm, left_key, right_key = contrast_specs[contrast]
        if (row["left_arm"].strip(), row["right_arm"].strip(),
                row["left_arm_key"].strip(),
                row["right_arm_key"].strip()) != (
                    left_arm, right_arm, left_key, right_key):
            raise ValueError(
                f"lib{library} contrast {contrast} has the wrong arm mapping")
        key = (library, contrast, population)
        if key in summary_rows:
            raise ValueError(
                f"lib{library} has duplicate {contrast}/{population} "
                "contrast summaries")
        summary_rows[key] = row

    for row in iter_strict_posthoc_tsv(
            contrast_cells_path,
            schemas["reconciliation_planned_contrast_cells.tsv"],
            "applied four-arm aggregate "
            "reconciliation_planned_contrast_cells.tsv"):
        library = selected_library(row, "four-arm cell contrast")
        if not row["barcode"].strip():
            raise ValueError(
                f"lib{library} four-arm cell contrast has an empty barcode")
        contrast = row["contrast"].strip()
        if contrast not in required_contrasts_by_library[library]:
            raise ValueError(
                f"lib{library} has an unexpected cell contrast: {contrast}")
        barcode = row["barcode"].strip()
        stratum = row["stratum"].strip() or "unmapped"
        authoritative_stratum = validation_store.stratum(library, barcode)
        if stratum != authoritative_stratum:
            raise ValueError(
                f"lib{library}/{barcode} aggregate stratum {stratum!r} does "
                f"not match hashed scrutiny stratum "
                f"{authoritative_stratum!r}")
        try:
            left = float(row["left_rate"])
            right = float(row["right_rate"])
            delta = float(row["right_minus_left"])
        except ValueError as exc:
            raise ValueError(
                f"lib{library} four-arm cell contrast has nonnumeric rates") \
                from exc
        if (not all(math.isfinite(value) and 0 <= value < 1
                    for value in (left, right)) or
                not math.isfinite(delta) or
                not math.isclose(delta, right - left, rel_tol=0.0,
                                 abs_tol=1e-8)):
            raise ValueError(
                f"lib{library} four-arm cell contrast rate arithmetic is invalid")
        expected_left_arm, expected_right_arm = contrast_specs[contrast][:2]
        direct_left = validation_store.rate(
            library, expected_left_arm, barcode)
        direct_right = validation_store.rate(
            library, expected_right_arm, barcode)
        if (direct_left is None or
                not math.isclose(left, direct_left[0], rel_tol=0.0,
                                 abs_tol=1e-8)):
            raise ValueError(
                f"lib{library}/{barcode} aggregate {contrast} left rate does "
                f"not match direct Arm {expected_left_arm}")
        if (direct_right is None or
                not math.isclose(right, direct_right[0], rel_tol=0.0,
                                 abs_tol=1e-8)):
            raise ValueError(
                f"lib{library}/{barcode} aggregate {contrast} right rate does "
                f"not match direct Arm {expected_right_arm}")
        validation_store.add_contrast(
            library, contrast, barcode, stratum,
            row["left_rate"].strip(), row["right_rate"].strip(),
            row["right_minus_left"].strip(), delta,
            f"lib{library} four-arm cell contrast")

    def numeric_matches(observed: str, expected: float,
                        tolerance: float = 1e-8) -> bool:
        try:
            value = float(observed)
        except (TypeError, ValueError):
            return False
        if math.isnan(expected):
            return math.isnan(value)
        return (math.isfinite(value) and
                math.isclose(value, expected, rel_tol=0.0,
                             abs_tol=tolerance))

    for library in event_libraries:
        required_contrasts = required_contrasts_by_library[library]
        expected_count = validation_store.rate_count(library, "C")
        if (validation_store.scrutiny_count(library) != expected_count or
                not validation_store.same_scrutiny_rate_universe(
                    library, "C")):
            raise ValueError(
                f"lib{library} hashed scrutiny table does not cover the exact "
                "production Arm C cell universe")
        expected_arm_names = {
            arm_contract_specs[key][0]
            for key in required_arms_by_library[library]}
        if validation_store.arm_names(library) != expected_arm_names:
            raise ValueError(
                f"lib{library} direct ambient arm set does not match its plan")
        for arm_name in expected_arm_names:
            if not validation_store.same_rate_universe(
                    library, arm_name, "C"):
                raise ValueError(
                    f"lib{library} Arm {arm_name} rate universe differs from "
                    "production Arm C")
        for contrast in required_contrasts:
            if validation_store.contrast_count(
                    library, contrast) != expected_count:
                raise ValueError(
                    f"lib{library} aggregate contrast {contrast} does not "
                    "cover the complete production Arm C cell universe")
            summary_keys = {
                population for (summary_library, summary_contrast, population)
                in summary_rows
                if summary_library == library and
                summary_contrast == contrast}
            if summary_keys != expected_populations:
                raise ValueError(
                    f"lib{library} contrast {contrast} summary population set "
                    "is incomplete")
            for population in expected_populations:
                row = summary_rows[(library, contrast, population)]
                (count, expected_mean, expected_median,
                 expected_abs_median, expected_fraction) = (
                    validation_store.contrast_summary(
                        library, contrast, population))
                try:
                    observed_count = int(row["n_common"])
                except ValueError as exc:
                    raise ValueError(
                        f"lib{library} contrast {contrast}/{population} has "
                        "invalid n_common") from exc
                if (observed_count != count or
                        not numeric_matches(row["mean_delta"], expected_mean) or
                        not numeric_matches(
                            row["median_delta"], expected_median) or
                        not numeric_matches(
                            row["median_absolute_delta"],
                            expected_abs_median) or
                        not numeric_matches(
                            row["fraction_lower_in_right_arm"],
                            expected_fraction)):
                    raise ValueError(
                        f"lib{library} contrast {contrast}/{population} "
                        "summary does not reproduce from direct arm rates")

    expected_burden: dict[
        tuple[int, str, str, str], dict[str, float | int]] = {}
    for library in event_libraries:
        expected_arm_names = {
            arm_contract_specs[key][0]
            for key in required_arms_by_library[library]}
        for arm in expected_arm_names:
            profile = arm_summaries.get((library, arm))
            if (profile is None or profile.cell_count !=
                    validation_store.rate_count(library, arm)):
                raise ValueError(
                    f"lib{library} Arm {arm} profile universe differs from "
                    "its rate universe")
            for population in expected_populations:
                count = profile.population_counts.get(population, 0)
                if count == 0:
                    continue
                mean_total = profile.population_rate_sums[population] / count
                population_sources = {
                    source for (source_population, source) in
                    profile.source_mass_sums
                    if source_population == population}
                for source in population_sources:
                    mean_mass = profile.source_mass_sums[
                        (population, source)] / count
                    mean_burden = profile.source_burden_sums[
                        (population, source)] / count
                    expected_burden[(library, arm, population, source)] = {
                        "n_cells": count,
                        "mean_scoring_profile_mass": mean_mass,
                        "mean_exact_contam_burden": mean_burden,
                        "fraction_of_mean_contam_burden": (
                            mean_burden / mean_total
                            if mean_total > 0 else math.nan),
                        "mean_total_contam_rate": mean_total,
                    }

    burden_groups: dict[
        tuple[int, str, str], list[dict[str, str]]] = {}
    observed_burden: set[tuple[int, str, str, str]] = set()
    for row in tables["reconciliation_exact_donor_burden.tsv"]:
        library = selected_library(row, "four-arm exact donor burden")
        arm = row["arm"].strip()
        arm_key = row["arm_key"].strip()
        expected = arm_contract_specs.get(arm_key)
        population = row["population"].strip()
        if (expected is None or arm_key not in
                required_arms_by_library[library] or
                (arm, row["assignment_basis"].strip(),
                 row["roster_basis"].strip()) != expected or
                row["series_key"].strip() !=
                f"{args.ambient_condition}__{arm_key}" or
                population not in expected_populations):
            raise ValueError(
                f"lib{library} has a malformed exact-donor burden row")
        source = row["source_label"].strip()
        key = (library, arm, population, source)
        if not source or key in observed_burden or key not in expected_burden:
            raise ValueError(
                f"lib{library} Arm {arm}/{population} has a duplicate or "
                "unexpected burden source")
        observed_burden.add(key)
        expected_values = expected_burden[key]
        try:
            n_cells = int(row["n_cells"])
        except ValueError as exc:
            raise ValueError(
                f"lib{library} Arm {arm}/{population} burden row is "
                "nonnumeric") from exc
        if (n_cells != expected_values["n_cells"] or
                any(not numeric_matches(row[field], float(expected_values[field]))
                    for field in (
                        "mean_scoring_profile_mass",
                        "mean_exact_contam_burden",
                        "fraction_of_mean_contam_burden",
                        "mean_total_contam_rate"))):
            raise ValueError(
                f"lib{library} Arm {arm}/{population} burden row does not "
                "reproduce from the direct per-cell source profile")
        burden_groups.setdefault((library, arm, population), []).append(row)
    if observed_burden != set(expected_burden):
        missing = sorted(set(expected_burden) - observed_burden)[:5]
        raise ValueError(
            "applied four-arm exact-donor burden table is incomplete; "
            f"missing examples={missing}")

    def natural_sort_key(value: str) -> tuple[object, ...]:
        return tuple(
            int(part) if part.isdigit() else part.lower()
            for part in re.split(r"(\d+)", value))

    def final_text(value: object) -> str:
        text = str(value).strip()
        return text if text and text.upper() != "NAN" else "NA"

    burden_text_by_population: dict[tuple[int, str], str] = {}
    background_text_by_library: dict[int, str] = {}
    for library in event_libraries:

        background_rows = [
            summary_rows[(library, contrast, "background")]
            for contrast in required_contrasts_by_library[library]]
        background_text_by_library[library] = ";".join(
            f"{row['contrast'].strip()}:n={final_text(row['n_common'])},"
            f"mean={final_text(row['mean_delta'])},"
            f"median={final_text(row['median_delta'])}"
            for row in sorted(
                background_rows,
                key=lambda value: natural_sort_key(value["contrast"]))) or "NA"
        for population in expected_populations:
            burden_rows = burden_groups.get((library, "C", population), [])
            if not burden_rows:
                continue
            burden_text_by_population[(library, population)] = ";".join(
                f"{row['source_label'].strip()}="
                f"{final_text(row['mean_exact_contam_burden'])}"
                for row in sorted(
                    burden_rows,
                    key=lambda value: natural_sort_key(
                        value["source_label"]))) or "NA"

    for row in tables["reconciliation_assignment_switch_summary.tsv"]:
        library = selected_library(row, "four-arm switch summary")
        try:
            switches = int(row["n_final_vs_input_switch"])
        except ValueError as exc:
            raise ValueError(
                f"lib{library} four-arm switch summary is nonnumeric") from exc
        if (switches != 0 or
                row["frozen_assignment_invariant"].strip().upper() != "PASS"):
            raise ValueError(
                f"lib{library} applied four-arm fit changed a frozen assignment")
    diagnostic_libraries = set()
    documented_warnings = 0
    for row in tables["reconciliation_diagnostics.tsv"]:
        library = selected_library(row, "four-arm diagnostic")
        diagnostic_libraries.add(library)
        status = row["status"].strip().upper()
        check = row["check"].strip()
        if status == "WARN" and check in {
                "profile_multistart", "model_final_log_likelihood"}:
            documented_warnings += 1
        elif status != "PASS":
            raise ValueError(
                f"lib{library} applied four-arm diagnostic has an unsupported "
                f"status {status!r}: {check}")
    if diagnostic_libraries != event_libraries:
        raise ValueError(
            "applied four-arm diagnostics do not cover every event library")
    lookup = FourArmFinalizerLookup(
        validation_store, required_contrasts_by_library,
        burden_text_by_population, background_text_by_library)
    return ({
        "status": "PASS", "root": root,
        "event_libraries": sorted(event_libraries),
        "documented_warning_rows": documented_warnings, "records": records,
    }, lookup)


def verify_empty_drop_generation(
        args: argparse.Namespace, paths: Sequence[LibraryPaths],
        final_boundary_mtime: int) -> dict[str, object]:
    scripts_root = os.path.join(
        args.upstream_analysis_root, "aggregate_library_analysis",
        "slurm_scripts")
    records: list[dict[str, object]] = [provenance_file_record(
        args.ambient_profile_binary, "EMPTY_DROPS ambient-profile binary", True)]
    for item in paths:
        demux_dir = os.path.dirname(item.demux_prefix)
        raw_prefix = os.path.join(demux_dir, f"lib{item.library}_raw")
        expected_lines = os.path.join(
            os.path.dirname(item.mapping_bam),
            f"lib{item.library}_expected_lines.txt")
        species_lines = os.path.join(
            demux_dir, f"lib{item.library}_expected_species_lines.txt")
        decompressed_barcodes = os.path.join(
            demux_dir, "filtered_barcodes.tsv")
        script_path = os.path.join(
            scripts_root, f"empty_lib{item.library}.sbatch")
        with open(script_path, "r", encoding="utf-8") as handle:
            script_payload = handle.read()
        if "Skipping (use --force to rerun)" in script_payload:
            raise ValueError(
                f"lib{item.library} EMPTY_DROPS was not generated with --force")
        commands = logical_program_argvs(
            script_path, os.path.basename(args.ambient_profile_binary))
        common = [
            args.ambient_profile_binary, "-o", raw_prefix,
            "--filtered_barcodes", decompressed_barcodes,
        ]
        expected_individual = common + [
            "--interindividual", "--condf", raw_prefix + ".condf",
            "--ids", expected_lines,
            "-T", "${SLURM_CPUS_PER_TASK}",
        ]
        expected_species = common + [
            "--interspecies", "--condf", raw_prefix + ".species_condf",
            "--ids", species_lines,
            "--panel_metadata", args.panel_metadata,
            "-T", "${SLURM_CPUS_PER_TASK}",
        ]
        individual = [argv for argv in commands
                      if list(argv) == expected_individual]
        species = [argv for argv in commands
                   if list(argv) == expected_species]
        if len(commands) != 2 or len(individual) != 1 or len(species) != 1:
            raise ValueError(
                f"lib{item.library} EMPTY_DROPS script does not contain the "
                "exact individual/species raw-count commands")
        records.append(provenance_file_record(
            script_path, f"lib{item.library} EMPTY_DROPS script", True))
        records.append(provenance_file_record(
            decompressed_barcodes,
            f"lib{item.library} decompressed filtered barcodes", True))
        if sha256_file(decompressed_barcodes) != gzip_payload_sha256(
                item.expression_barcodes):
            raise ValueError(
                f"lib{item.library} EMPTY_DROPS filtered barcode payload does "
                "not match the newest filtered MEX")
        newest_input = max(os.stat(path).st_mtime_ns for path in (
            script_path, raw_prefix + ".counts", raw_prefix + ".samples",
            raw_prefix + ".condf", raw_prefix + ".species_counts",
            raw_prefix + ".species_condf", raw_prefix + ".species_samples",
            decompressed_barcodes, expected_lines, species_lines,
            args.panel_metadata, args.ambient_profile_binary))
        for suffix, label in (
                (".contam_prof_empty", "individual empty-drop profile"),
                (".species_prof_empty", "species empty-drop profile")):
            output = raw_prefix + suffix
            output_mtime = os.stat(output).st_mtime_ns
            if not newest_input <= output_mtime <= final_boundary_mtime:
                raise ValueError(
                    f"lib{item.library} {label} is outside the forced "
                    "raw-input-to-final-reconciliation chronology")
            records.append(provenance_file_record(
                output, f"lib{item.library} {label}", True))
        validate_ambient_profile(
            raw_prefix + ".contam_prof_empty",
            identity_component_roster(
                expected_lines,
                f"lib{item.library} empty-drop individual donor roster"),
            f"lib{item.library} individual empty-drop profile")
        validate_ambient_profile(
            raw_prefix + ".species_prof_empty",
            identity_component_roster(
                species_lines,
                f"lib{item.library} empty-drop species roster"),
            f"lib{item.library} species empty-drop profile")
    return {"status": "PASS", "records": records}


def verify_ambient_generation(
        args: argparse.Namespace, paths: Sequence[LibraryPaths],
        event_libraries: set[int], final_boundary_mtime: int,
        validation_boundary_mtime: int
) -> dict[str, object]:
    """Run ambient validation with bounded disk-backed per-cell state."""
    validation_store = AmbientValidationStore(args.run_root)
    try:
        return _verify_ambient_generation(
            args, paths, event_libraries, final_boundary_mtime,
            validation_boundary_mtime, validation_store)
    finally:
        validation_store.close()


def _verify_ambient_generation(
        args: argparse.Namespace, paths: Sequence[LibraryPaths],
        event_libraries: set[int], final_boundary_mtime: int,
        validation_boundary_mtime: int,
        validation_store: AmbientValidationStore,
) -> dict[str, object]:
    """Bind standard and applied reconciliation ambient estimates to inputs."""
    scripts_root = os.path.join(
        args.upstream_analysis_root, "aggregate_library_analysis",
        "slurm_scripts")
    records: list[dict[str, object]] = [
        provenance_file_record(
            args.geometry_gate_helper, "geometry-gated ambient helper", True),
        provenance_file_record(
            args.contam_binary, "ambient contamination estimator", True),
    ]
    standard_details: dict[int, dict[str, object]] = {}
    arm_details: dict[int, dict[str, object]] = {}
    arm_summaries: dict[tuple[int, str], ArmProfileSummary] = {}
    comparison_contexts: dict[int, dict[str, object]] = {}
    newest_arm_output_mtime = 0
    for item in paths:
        core_assignments = strict_assignment_map(
            item.demux_prefix + ".assignments",
            f"lib{item.library} core DEMUX assignments")
        core_barcodes = set(core_assignments)
        demux_dir = os.path.dirname(item.demux_prefix)
        expected_lines = os.path.join(
            os.path.dirname(item.mapping_bam),
            f"lib{item.library}_expected_lines.txt")
        ambient_candidates = os.path.join(
            demux_dir, f"lib{item.library}_ambient_candidates.txt")
        expected_donors = identity_roster(
            ambient_candidates,
            f"lib{item.library} standard ambient donor roster")
        prefix = item.ambient_standard_prefix
        script_path = os.path.join(
            scripts_root,
            f"qc_{args.ambient_condition}_lib{item.library}.sbatch")
        validate_geometry_command(
            args, script_path, prefix, expected_lines, ambient_candidates,
            item.demux_prefix + ".condf", False,
            f"lib{item.library} standard ambient script")
        for suffix in (".counts", ".assignments", ".samples", ".condf"):
            staged = prefix + suffix
            source = item.demux_prefix + suffix
            if os.path.realpath(staged) != os.path.realpath(source):
                raise ValueError(
                    f"lib{item.library} standard ambient staged {suffix} does "
                    "not resolve to the selected forced DEMUX bundle")
        contract_path = prefix + ".run_contract.json"
        validate_geometry_contract(
            contract_path, prefix, args.ambient_condition, {
                "counts": prefix + ".counts",
                "condf": item.demux_prefix + ".condf",
                "assignments": prefix + ".assignments",
                "samples": prefix + ".samples",
                "expected_lines": expected_lines,
                "ambient_candidates": ambient_candidates,
            }, False, f"lib{item.library} standard ambient run contract")
        required_suffixes = (
            ".contam_rate", ".contam_prof", ".allele_ratio",
            ".contam_diagnostics.tsv", ".profile_fit_diagnostics.tsv",
            ".condf_coverage.tsv", ".run_contract.json",
            ".geometry_gate_audit.tsv",
        )
        newest_input = max(os.stat(path).st_mtime_ns for path in (
            script_path, item.demux_prefix + ".counts",
            item.demux_prefix + ".assignments", item.demux_prefix + ".samples",
            item.demux_prefix + ".condf", expected_lines,
            ambient_candidates, args.geometry_gate_helper,
            args.contam_binary))
        for suffix in required_suffixes:
            output = prefix + suffix
            output_mtime = os.stat(output).st_mtime_ns
            if not newest_input <= output_mtime <= validation_boundary_mtime:
                raise ValueError(
                    f"lib{item.library} standard ambient artifact is outside "
                    f"the forced DEMUX-to-identity-validation chronology: "
                    f"{output}")
            records.append(provenance_file_record(
                output, f"lib{item.library} standard ambient {suffix}", True))
        standard_selected, _standard_profile_masses = (
            validate_geometry_selected_bundle(
                prefix, core_barcodes,
                f"lib{item.library} standard ambient", core_assignments,
                expected_donors, require_cell_profile=False,
                validation_store=validation_store,
                library=item.library, arm="S"))
        standard_details[item.library] = {
            "rates": {
                "rows": standard_selected["rows"],
                "finite_rates": standard_selected["rows"],
            },
            "profile": validate_ambient_profile(
                prefix + ".contam_prof", expected_donors,
                f"lib{item.library} standard ambient profile"),
            "selected_bundle": standard_selected,
        }
        records.append(provenance_file_record(
            script_path, f"lib{item.library} standard ambient script", True))

        if item.library not in event_libraries:
            continue
        plan_dir = os.path.join(
            demux_dir, "reconciliation_four_arm", args.ambient_candidate_set,
            "plan")
        context_path = os.path.join(
            plan_dir, f"lib{item.library}.comparison_context.json")
        try:
            with open(context_path, "r", encoding="utf-8") as handle:
                context = json.load(handle)
        except (OSError, UnicodeError, json.JSONDecodeError) as exc:
            raise ValueError(
                f"lib{item.library} invalid reconciliation ambient context: "
                f"{exc}") from exc
        if (context.get("schema_version") != 2 or
                context.get("plan_version") !=
                "identity_ambient_comparison_plan_V4" or
                context.get("candidate_set") != args.ambient_candidate_set or
                context.get("library") != item.library or
                context.get("assignment_score_basis") !=
                "original_demux_all_arms"):
            raise ValueError(
                f"lib{item.library} reconciliation ambient context is not the "
                "applied V4 final-call comparison")
        plan_dir = os.path.dirname(context_path)
        expected_input_paths = {
            "original_receiver_lines": os.path.join(
                os.path.dirname(item.mapping_bam),
                f"lib{item.library}_expected_lines.txt"),
            "original_ambient_candidates": os.path.join(
                demux_dir, f"lib{item.library}_ambient_candidates.txt"),
            "demux_assignments": item.demux_prefix + ".assignments",
            "validated_reconciled_assignments": os.path.join(
                args.identity_root, "decisions",
                f"lib{item.library}.reconciled_single_cell.assignments"),
            "reconciled_cells": os.path.join(
                args.identity_root, "decisions",
                f"lib{item.library}.reconciled_cells.tsv.gz"),
            "demux_samples": item.samples,
        }
        expected_plan_paths = {
            "reconciled_assignments": os.path.join(
                plan_dir,
                f"lib{item.library}.comparison_reconciled.assignments"),
            "scrutiny_cells": os.path.join(
                plan_dir, f"lib{item.library}.scrutiny_cells.tsv"),
            "ambient_candidate_provenance": os.path.join(
                plan_dir,
                f"lib{item.library}.ambient_candidate_provenance.tsv"),
            "augmented_receiver_lines": os.path.join(
                plan_dir,
                f"lib{item.library}.augmented_receiver_lines.txt"),
            "augmented_ambient_candidates": os.path.join(
                plan_dir,
                f"lib{item.library}.augmented_ambient_candidates.txt"),
        }
        replacement_eligible = context.get("replacement_arm_eligible") is True
        if replacement_eligible:
            expected_plan_paths.update({
                "replacement_receiver_lines": os.path.join(
                    plan_dir,
                    f"lib{item.library}.replacement_receiver_lines.txt"),
                "replacement_ambient_candidates": os.path.join(
                    plan_dir,
                    f"lib{item.library}.replacement_ambient_candidates.txt"),
            })
        for group_name, expected_paths in (
                ("input_records", expected_input_paths),
                ("plan_artifact_records", expected_plan_paths)):
            group = context.get(group_name)
            if not isinstance(group, dict) or set(group) != set(expected_paths):
                raise ValueError(
                    f"lib{item.library} ambient context {group_name} does not "
                    "have the exact V4 producer record set")
            for name, expected_path in expected_paths.items():
                validate_sized_path_record(
                    group[name], expected_path,
                    f"lib{item.library} ambient context {group_name}/{name}",
                    True)
        fingerprint_payload = {
            "plan_version": "identity_ambient_comparison_plan_V4",
            "candidate_set": args.ambient_candidate_set,
            "inputs": context["input_records"],
        }
        expected_fingerprint = hashlib.sha256(json.dumps(
            fingerprint_payload, sort_keys=True,
            separators=(",", ":")).encode("utf-8")).hexdigest()
        if context.get("plan_fingerprint") != expected_fingerprint:
            raise ValueError(
                f"lib{item.library} reconciliation ambient plan fingerprint "
                "does not reproduce from its exact input records")
        context_inputs = context["input_records"]
        context_samples_path = str(
            context_inputs.get("demux_samples", {}).get("path", ""))
        if (absolute(str(context.get("demux_assignments", ""))) !=
                item.demux_prefix + ".assignments" or
                absolute(context_samples_path) != item.samples):
            raise ValueError(
                f"lib{item.library} ambient context is not tied to the selected "
                "forced DEMUX assignments/sample order")
        comparison_assignments_path = absolute(str(
            context.get("reconciled_assignments", "")))
        exact_context_paths = {
            "demux_assignments": expected_input_paths["demux_assignments"],
            "validated_reconciled_assignments": expected_input_paths[
                "validated_reconciled_assignments"],
            "reconciled_cells": expected_input_paths["reconciled_cells"],
            "scrutiny_cells": expected_plan_paths["scrutiny_cells"],
            "ambient_candidate_provenance": expected_plan_paths[
                "ambient_candidate_provenance"],
            "original_receiver_lines": expected_input_paths[
                "original_receiver_lines"],
            "original_ambient_candidates": expected_input_paths[
                "original_ambient_candidates"],
            "augmented_receiver_lines": expected_plan_paths[
                "augmented_receiver_lines"],
            "augmented_ambient_candidates": expected_plan_paths[
                "augmented_ambient_candidates"],
            "replacement_receiver_lines": expected_plan_paths.get(
                "replacement_receiver_lines", ""),
            "replacement_ambient_candidates": expected_plan_paths.get(
                "replacement_ambient_candidates", ""),
        }
        if (comparison_assignments_path !=
                expected_plan_paths["reconciled_assignments"] or
                any(absolute(str(context.get(key, ""))) != absolute(value)
                    if value else str(context.get(key, "")) != ""
                    for key, value in exact_context_paths.items())):
            raise ValueError(
                f"lib{item.library} reconciliation ambient context paths do "
                "not match the exact V4 input/plan artifacts")
        demux_order, demux_assignment_rows = strict_assignment_rows(
            item.demux_prefix + ".assignments",
            f"lib{item.library} exact DEMUX assignments")
        demux_identities = {
            barcode: values[0]
            for barcode, values in demux_assignment_rows.items()}
        validated_reconciled = assignment_identity_map(
            expected_input_paths["validated_reconciled_assignments"],
            f"lib{item.library} validated reconciliation assignments")
        if not set(validated_reconciled) <= set(demux_assignment_rows):
            raise ValueError(
                f"lib{item.library} validated reconciliation contains cells "
                "outside the exact DEMUX universe")
        comparison_order, comparison_assignment_rows = strict_assignment_rows(
            comparison_assignments_path,
            f"lib{item.library} ambient comparison assignments")
        if comparison_order != demux_order:
            raise ValueError(
                f"lib{item.library} ambient comparison does not preserve the "
                "exact DEMUX barcode order and universe")
        for barcode in demux_order:
            demux_identity, _demux_type, demux_llr = (
                demux_assignment_rows[barcode])
            expected_identity = model_identity(
                validated_reconciled.get(barcode, demux_identity))
            expected_type = "D" if "+" in expected_identity else "S"
            observed_identity, observed_type, observed_llr = (
                comparison_assignment_rows[barcode])
            if ((observed_identity, observed_type, observed_llr) !=
                    (expected_identity, expected_type, demux_llr)):
                raise ValueError(
                    f"lib{item.library}/{barcode} ambient comparison must "
                    "overlay only the validated identity, recompute S/D, and "
                    "preserve the exact original DEMUX LLR token")
        comparison_assignments = {
            barcode: values[0]
            for barcode, values in comparison_assignment_rows.items()}
        plan_expectation = reconstruct_ambient_plan_expectation(
            expected_input_paths["reconciled_cells"],
            args.ambient_candidate_set, demux_order, demux_identities,
            validated_reconciled,
            expected_input_paths["original_receiver_lines"],
            expected_input_paths["original_ambient_candidates"],
            f"lib{item.library} reconciliation ambient plan")
        plan_counts = validate_ambient_plan_artifacts(
            plan_expectation, context, expected_plan_paths, item.samples,
            args.ambient_candidate_set,
            f"lib{item.library} reconciliation ambient plan")
        final_assignments = assignment_identity_map(
            item.final_assignments,
            f"lib{item.library} Phase-3 final assignments")
        if set(comparison_assignments) != set(final_assignments):
            raise ValueError(
                f"lib{item.library} applied ambient comparison and Phase-3 "
                "final assignments do not have the same cell universe")
        checked_final = 0
        for barcode, final_identity in final_assignments.items():
            if final_identity.startswith("M{"):
                continue
            if (barcode not in comparison_assignments or
                    model_identity(comparison_assignments[barcode]) !=
                    model_identity(final_identity)):
                raise ValueError(
                    f"lib{item.library}/{barcode} applied reconciliation "
                    "ambient plan does not use the Phase-3 final assignment")
            checked_final += 1
        if checked_final < 1:
            raise ValueError(
                f"lib{item.library} has no final single-cell assignments "
                "represented in the applied reconciliation ambient plan")
        comparison_contexts[item.library] = context
        scrutiny_counts = validation_store.add_scrutiny(
            item.library, expected_plan_paths["scrutiny_cells"],
            f"lib{item.library} hashed reconciliation scrutiny cells",
            plan_expectation)
        expected_context_counts = {
            "n_comparison_barcodes": scrutiny_counts["rows"],
            "n_validated_single_cell_barcodes": plan_counts[
                "validated_barcodes"],
            "n_preserved_non_single_cell_barcodes": plan_counts[
                "preserved_barcodes"],
            "n_changed": scrutiny_counts["changed"],
            "n_scrutinized": scrutiny_counts["scrutinized"],
            "n_background": scrutiny_counts["background"],
        }
        for count_name, expected_count in expected_context_counts.items():
            observed_count = context.get(count_name)
            if (not isinstance(observed_count, int) or
                    isinstance(observed_count, bool) or
                    observed_count != expected_count):
                raise ValueError(
                    f"lib{item.library} ambient context {count_name} does not "
                    "reproduce from its hashed scrutiny/assignment inputs")

        arm_details[item.library] = {}
        arm_root = os.path.join(
            demux_dir, "contamination", args.ambient_condition,
            "reconciliation_four_arm", args.ambient_candidate_set)
        arm_specs = [
            ("demux_original", "A", item.ambient_arm_a_prefix,
             "original_receiver_lines", "original_ambient_candidates",
             item.demux_prefix + ".assignments", "demux", "original"),
            ("demux_augmented", "B", os.path.join(
                arm_root, "demux_augmented",
                f"lib{item.library}_demuxed"),
             "original_receiver_lines", "augmented_ambient_candidates",
             item.demux_prefix + ".assignments", "demux", "augmented"),
            ("reconciled_augmented", "C", item.ambient_arm_c_prefix,
             "augmented_receiver_lines", "augmented_ambient_candidates",
             comparison_assignments_path, "reconciled", "augmented"),
        ]
        if replacement_eligible:
            arm_specs.append((
                "reconciled_replacement", "D", os.path.join(
                    arm_root, "reconciled_replacement",
                    f"lib{item.library}_demuxed"),
                "replacement_receiver_lines",
                "replacement_ambient_candidates",
                comparison_assignments_path, "reconciled", "replacement"))
        for (arm_key, arm_name, arm_prefix, receiver_key, candidate_key,
             expected_assignment_path, expected_basis,
             expected_roster) in arm_specs:
            receiver = absolute(str(context.get(receiver_key, "")))
            candidates = absolute(str(context.get(candidate_key, "")))
            arm_script = os.path.join(
                scripts_root,
                f"qc_{args.ambient_condition}_arm{arm_name}_lib"
                f"{item.library}.sbatch")
            validate_geometry_command(
                args, arm_script, arm_prefix, receiver, candidates,
                item.demux_prefix + ".condf", True,
                f"lib{item.library} reconciliation ambient Arm {arm_name}")
            arm_order, arm_assignment_rows = strict_assignment_rows(
                arm_prefix + ".assignments",
                f"lib{item.library} ambient Arm {arm_name} assignments")
            expected_order, expected_assignment_rows = strict_assignment_rows(
                expected_assignment_path,
                f"lib{item.library} ambient Arm {arm_name} source assignments")
            if (arm_order != expected_order or
                    arm_assignment_rows != expected_assignment_rows):
                raise ValueError(
                    f"lib{item.library} ambient Arm {arm_name} assignments do "
                    "not exactly preserve the declared barcode order, "
                    "identity, S/D type, and LLR")
            arm_assignments = {
                barcode: values[0]
                for barcode, values in arm_assignment_rows.items()}
            contract_path = arm_prefix + ".run_contract.json"
            validate_geometry_contract(
                contract_path, arm_prefix, args.ambient_condition, {
                    "counts": arm_prefix + ".counts",
                    "condf": item.demux_prefix + ".condf",
                    "assignments": arm_prefix + ".assignments",
                    "samples": arm_prefix + ".samples",
                    "expected_lines": receiver,
                    "ambient_candidates": candidates,
                }, True,
                f"lib{item.library} ambient Arm {arm_name} run contract")
            arm_contract_path = arm_prefix + ".identity_ambient_arm.tsv"
            arm_contract_rows = strict_posthoc_tsv(
                arm_contract_path, (
                    "library", "condition", "arm", "arm_key",
                    "assignment_basis", "roster_basis", "candidate_set",
                    "plan_fingerprint", "assignment_path", "receiver_lines",
                    "ambient_candidates", "scrutiny_cells", "context_path",
                    "assignment_update_mode", "assignment_score_basis"),
                f"lib{item.library} ambient Arm {arm_name} contract")
            if len(arm_contract_rows) != 1:
                raise ValueError(
                    f"lib{item.library} ambient Arm {arm_name} contract must "
                    "have exactly one row")
            arm_row = arm_contract_rows[0]
            if (arm_row["library"] != f"lib{item.library}" or
                    arm_row["condition"] != args.ambient_condition or
                    arm_row["arm"] != arm_name or
                    arm_row["arm_key"] != arm_key or
                    arm_row["assignment_basis"] != expected_basis or
                    arm_row["roster_basis"] != expected_roster or
                    arm_row["candidate_set"] != args.ambient_candidate_set or
                    arm_row["plan_fingerprint"] !=
                    str(context.get("plan_fingerprint", "")) or
                    absolute(arm_row["assignment_path"]) !=
                    expected_assignment_path or
                    absolute(arm_row["receiver_lines"]) != receiver or
                    absolute(arm_row["ambient_candidates"]) != candidates or
                    absolute(arm_row["scrutiny_cells"]) !=
                    absolute(str(context.get("scrutiny_cells", ""))) or
                    absolute(arm_row["context_path"]) != context_path or
                    arm_row["assignment_update_mode"] != "iterative_frozen" or
                    arm_row["assignment_score_basis"] !=
                    "original_demux_all_arms"):
                raise ValueError(
                    f"lib{item.library} ambient Arm {arm_name} contract does "
                    "not match the applied reconciliation plan")
            arm_candidates = identity_roster(
                candidates,
                f"lib{item.library} ambient Arm {arm_name} donor roster")
            required_arm_suffixes = (
                ".contam_rate", ".contam_prof", ".allele_ratio",
                ".contam_diagnostics.tsv", ".cell_source_profile.tsv",
                ".profile_fit_diagnostics.tsv", ".condf_coverage.tsv",
                ".run_contract.json", ".decontam.assignments",
                ".identity_ambient_arm.tsv", ".geometry_gate_audit.tsv",
            )
            newest_arm_input = max(os.stat(path).st_mtime_ns for path in (
                arm_script, context_path, expected_assignment_path, receiver,
                candidates, item.demux_prefix + ".counts",
                item.demux_prefix + ".condf", item.samples,
                args.geometry_gate_helper, args.contam_binary))
            for suffix in required_arm_suffixes:
                output = arm_prefix + suffix
                output_mtime = os.stat(output).st_mtime_ns
                newest_arm_output_mtime = max(
                    newest_arm_output_mtime, output_mtime)
                if not newest_arm_input <= output_mtime <= final_boundary_mtime:
                    raise ValueError(
                        f"lib{item.library} ambient Arm {arm_name} artifact is "
                        f"outside the plan-to-final chronology: {output}")
                records.append(provenance_file_record(
                    output,
                    f"lib{item.library} reconciliation ambient Arm "
                    f"{arm_name} {suffix}", True))
            selected_detail, profile_masses = validate_geometry_selected_bundle(
                arm_prefix, set(arm_assignments),
                f"lib{item.library} ambient Arm {arm_name}",
                arm_assignments, arm_candidates,
                validation_store=validation_store,
                library=item.library, arm=arm_name)
            if profile_masses is None:
                raise ValueError(
                    f"lib{item.library} ambient Arm {arm_name} did not yield "
                    "a direct cell/source profile summary")
            arm_details[item.library][arm_name] = {
                "rates": {
                    "rows": selected_detail["rows"],
                    "finite_rates": selected_detail["rows"],
                },
                "profile": validate_ambient_profile(
                    arm_prefix + ".contam_prof", arm_candidates,
                    f"lib{item.library} ambient Arm {arm_name} profile"),
                "selected_bundle": selected_detail,
                "final_assignments_checked": checked_final,
            }
            arm_summaries[(item.library, arm_name)] = profile_masses
            records.append(provenance_file_record(
                arm_script,
                f"lib{item.library} reconciliation ambient Arm {arm_name} "
                "script", True))
        records.append(provenance_file_record(
            context_path,
            f"lib{item.library} applied reconciliation ambient context", True))
    aggregate_detail, aggregate_finalizer_fields = verify_four_arm_aggregate(
        args, event_libraries, newest_arm_output_mtime,
        final_boundary_mtime, comparison_contexts, validation_store,
        arm_summaries)
    records.extend(aggregate_detail.get("records", []))

    # Phase 3 selects Arm C uniformly for event-bearing libraries and Arm A
    # for zero-event libraries.  Verify that each finalized rate is taken from
    # the selected fit.  The upstream finalizer records NA for zero-event c;
    # in that documented case the standard demux fit is the Arm-A equivalent.
    opener = gzip.open if args.ledger_input.endswith(".gz") else open
    ledger_counts: dict[int, int] = {item.library: 0 for item in paths}
    selected_libraries = set(ledger_counts)
    with opener(args.ledger_input, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        exact_ambient_fields = {
            "ambient_arm_a_c", "ambient_arm_b_c", "ambient_arm_c_c",
            "ambient_arm_d_c", "ambient_roster_effect_b_minus_a",
            "ambient_assignment_effect_c_minus_b",
            "ambient_replacement_effect_d_minus_c",
            "ambient_combined_augmented_c_minus_a",
            "ambient_production_arm", "ambient_production_c",
            "ambient_production_minus_original_c",
            "ambient_exact_donor_burden_fields",
            "ambient_background_shift_fields", "ambient_evaluation_status",
        }
        fields = set(reader.fieldnames or [])
        if not exact_ambient_fields <= fields:
            raise ValueError(
                "final reconciliation ledger lacks exact four-arm ambient "
                "fields: " + ",".join(sorted(exact_ambient_fields - fields)))

        def final_numeric_text(value: object) -> str:
            text = "" if value is None else str(value).strip()
            return "NA" if not text or text.upper() == "NAN" else text

        numeric_fields = (
            "ambient_arm_a_c", "ambient_arm_b_c", "ambient_arm_c_c",
            "ambient_arm_d_c", "ambient_roster_effect_b_minus_a",
            "ambient_assignment_effect_c_minus_b",
            "ambient_replacement_effect_d_minus_c",
            "ambient_combined_augmented_c_minus_a",
            "ambient_production_c",
            "ambient_production_minus_original_c",
        )
        for row in reader:
            library = normalize_library_number(row.get("library"))
            if library not in selected_libraries:
                continue
            barcode = str(row.get("barcode", "")).strip()
            expected_arm = "C" if library in event_libraries else "A"
            observed_arm = str(
                row.get("ambient_production_arm", "")).strip().upper()
            rate_record = validation_store.rate(
                library, "C" if library in event_libraries else "S", barcode)
            if (not barcode or rate_record is None or
                    observed_arm != expected_arm):
                raise ValueError(
                    f"lib{library}/{barcode or 'EMPTY'} is not bound to its "
                    f"finalized Arm {expected_arm} ambient-rate source")
            if library in event_libraries:
                expected = aggregate_finalizer_fields.get(library, barcode)
                if expected is None:
                    raise ValueError(
                        f"lib{library}/{barcode} has no exact four-arm "
                        "aggregate record")
                if (any(final_numeric_text(row.get(field)) !=
                        final_numeric_text(expected.get(field))
                        for field in numeric_fields) or
                        str(row.get("ambient_production_arm", "")).strip() !=
                        expected["ambient_production_arm"] or
                        str(row.get("ambient_evaluation_status", "")).strip() !=
                        expected["ambient_evaluation_status"] or
                        str(row.get(
                            "ambient_exact_donor_burden_fields", "")).strip() !=
                        expected["ambient_exact_donor_burden_fields"] or
                        str(row.get(
                            "ambient_background_shift_fields", "")).strip() !=
                        expected["ambient_background_shift_fields"]):
                    raise ValueError(
                        f"lib{library}/{barcode} Phase-3 four-arm ambient "
                        "fields do not reproduce from the direct A/B/C/D fits")
            else:
                zero_event_numeric = {
                    field: ("0" if field ==
                            "ambient_production_minus_original_c" else "NA")
                    for field in numeric_fields}
                if (str(row.get(
                            "ambient_production_arm", "")).strip() != "A" or
                        any(final_numeric_text(row.get(field)) !=
                        final_numeric_text(zero_event_numeric[field])
                        for field in numeric_fields) or
                        str(row.get("ambient_evaluation_status", "")).strip() !=
                        "NOT_APPLICABLE_ZERO_EVENT" or
                        str(row.get(
                            "ambient_exact_donor_burden_fields", "")).strip() !=
                        "NA" or
                        str(row.get(
                            "ambient_background_shift_fields", "")).strip() !=
                        "NA"):
                    raise ValueError(
                        f"lib{library}/{barcode} zero-event Phase-3 ambient "
                        "fields do not match the upstream finalizer contract")
            ledger_counts[library] += 1
    if any(count < 1 for count in ledger_counts.values()):
        raise ValueError(
            "final ambient production validation did not cover every upstream "
            "cohort library")
    return {
        "status": "PASS", "standard": standard_details,
        "reconciliation_arms": arm_details,
        "four_arm_aggregate": aggregate_detail,
        "final_ambient_rows": ledger_counts, "records": records,
    }


def verify_ploidy_generation(
        args: argparse.Namespace, paths: Sequence[LibraryPaths],
        libraries: Sequence[int]) -> dict[str, object]:
    """Prove that newest-MEX H5AD, NN job, and per-library calls agree."""
    if not args.ploidy_input_h5ad:
        raise ValueError(
            "a non-production provenance check requires "
            "--ploidy-input-h5ad so old production NN inputs cannot be reused")
    h5ad = absolute(args.ploidy_input_h5ad)
    if (not args.allow_external_input_templates and
            not lexical_path_contains(args.upstream_analysis_root, h5ad)):
        raise ValueError(
            "--ploidy-input-h5ad must be inside --upstream-analysis-root "
            "unless --allow-external-input-templates is explicitly set")
    provenance_file_record(h5ad, "newest-remap ploidy H5AD")
    with open(h5ad, "rb") as handle:
        if handle.read(8) != b"\x89HDF\r\n\x1a\n":
            raise ValueError(f"ploidy input is not an HDF5/H5AD file: {h5ad}")
    contract_path = h5ad + ".contract.json"
    contract_record = provenance_file_record(
        contract_path, "newest-remap ploidy H5AD contract", True)
    try:
        with open(contract_path, "r", encoding="utf-8") as handle:
            contract = json.load(handle)
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise ValueError(f"invalid ploidy H5AD contract: {contract_path}: {exc}") \
            from exc
    try:
        contract_libraries = {
            int(value) for value in contract.get("libraries", [])}
    except (TypeError, ValueError) as exc:
        raise ValueError("ploidy H5AD contract has an invalid library set") from exc
    selected_libraries = {int(value) for value in libraries}
    if (contract.get("schema_version") != "tetra_arm_ploidy_h5ad_v3" or
            str(contract.get("status", "")).upper() != "PASS" or
            absolute(str(contract.get("mapping_input_root", ""))) !=
            absolute(args.mapping_input_root) or
            not selected_libraries <= contract_libraries):
        raise ValueError(
            "ploidy H5AD contract does not match the selected mapping root and "
            f"library set: {contract_path}")
    output_ok, output_detail = record_matches_path(
        contract.get("output"), h5ad, True)
    if not output_ok:
        raise ValueError(
            f"ploidy H5AD no longer matches its producer contract: "
            f"{output_detail}")
    h5ad_record = dict(contract["output"])
    h5ad_record["label"] = "newest-remap ploidy H5AD"
    normalization = contract.get("normalization")
    if (not isinstance(normalization, dict) or
            normalization.get("transform") !=
            "X=ln(1+counts*target_sum_by_library[library]/cell_total_counts)" or
            normalization.get("normalization_scope") != "per_library"):
        raise ValueError(
            "ploidy H5AD contract lacks the validated normalize_total + "
            "natural-log1p transform")
    try:
        tolerance = float(normalization.get("tolerance"))
        maximum_error = float(normalization.get("maximum_absolute_X_error"))
        reference_total = int(normalization.get("reference_total_cells"))
        reference_checked = int(normalization.get("reference_cells_checked"))
        reference_nonempty = int(
            normalization.get("reference_nonempty_cells"))
        reference_zero = int(
            normalization.get("reference_zero_count_cells"))
        targets = {
            int(key): float(value) for key, value in
            dict(normalization.get("target_sum_by_library", {})).items()}
        spreads = {
            int(key): float(value) for key, value in dict(
                normalization.get(
                    "maximum_target_relative_spread_by_library", {})).items()}
        checked_by_library = {
            int(key): int(value) for key, value in dict(
                normalization.get(
                    "reference_cells_checked_by_library", {})).items()}
        nonempty_by_library = {
            int(key): int(value) for key, value in dict(
                normalization.get(
                    "reference_nonempty_cells_by_library", {})).items()}
    except (TypeError, ValueError) as exc:
        raise ValueError("ploidy H5AD normalization proof is malformed") from exc
    if (not all(math.isfinite(value) for value in (
            tolerance, maximum_error)) or not 0 < tolerance <= 0.01 or
            not 0 <= maximum_error <= tolerance or
            normalization.get("reference_scan_mode") !=
            "all_cells_bounded_chunks" or
            reference_total < 2 or reference_checked != reference_total or
            not 2 <= reference_nonempty <= reference_total or
            not 0 <= reference_zero <= reference_total or
            reference_nonempty + reference_zero != reference_total):
        raise ValueError("ploidy H5AD normalization proof did not pass")
    try:
        reference_chunk_cells = int(normalization.get("reference_chunk_cells"))
    except (TypeError, ValueError) as exc:
        raise ValueError(
            "ploidy H5AD normalization proof has an invalid chunk size") from exc
    if not 1 <= reference_chunk_cells <= 4096:
        raise ValueError(
            "ploidy H5AD normalization proof has an invalid chunk size")
    for library in selected_libraries:
        if (library not in targets or library not in spreads or
                library not in checked_by_library or
                library not in nonempty_by_library or
                not math.isfinite(targets[library]) or targets[library] <= 0 or
                not math.isfinite(spreads[library]) or
                not 0 <= spreads[library] <= tolerance or
                not 0 <= checked_by_library[library] <= reference_total or
                not 0 <= nonempty_by_library[library] <= reference_total or
                checked_by_library[library] < nonempty_by_library[library] or
                nonempty_by_library[library] < 2):
            raise ValueError(
                f"ploidy H5AD normalization proof is invalid for lib{library}")
    if contract.get("expression_source") != "newest filtered MEX only":
        raise ValueError("ploidy H5AD contract has an unsupported expression source")
    contract_inputs = contract.get("inputs")
    if not isinstance(contract_inputs, list):
        raise ValueError("ploidy H5AD contract has no MEX input records")
    by_library = {
        int(record.get("library")): record for record in contract_inputs
        if isinstance(record, dict) and str(record.get("library", "")).isdigit()
    }
    if not selected_libraries <= set(by_library):
        raise ValueError(
            "ploidy H5AD contract lacks a selected library's MEX records")
    validated_mex_records: list[dict[str, object]] = []
    for item in paths:
        input_record = by_library[item.library]
        for key, expected_path in (
                ("barcodes", item.expression_barcodes),
                ("features", item.expression_features),
                ("matrix", item.expression_matrix)):
            matched, detail = record_matches_path(
                input_record.get(key), expected_path, True)
            if not matched:
                raise ValueError(
                    f"lib{item.library} {key} no longer matches the ploidy "
                    f"H5AD contract: {detail}")
            frozen_record = dict(input_record[key])
            frozen_record["label"] = (
                f"lib{item.library} ploidy H5AD source {key}")
            validated_mex_records.append(frozen_record)

    scripts_root = os.path.join(
        args.upstream_analysis_root, "aggregate_library_analysis",
        "slurm_scripts")
    output_roots = {os.path.dirname(item.ploidy_nn) for item in paths}
    if len(output_roots) != 1:
        raise ValueError(
            "selected PLOIDY_NN call templates do not share one output directory")
    output_root = next(iter(output_roots))
    weights = absolute(args.ploidy_nn_weights)
    scaler = (weights[:-3] + "_scaler.npz"
              if weights.endswith(".pt") else weights + "_scaler.npz")
    helper = absolute(args.ploidy_nn_helper)
    weight_record = provenance_file_record(
        weights, "PLOIDY_NN weights", True)
    scaler_record = provenance_file_record(
        scaler, "PLOIDY_NN scaler", True)
    helper_record = provenance_file_record(helper, "PLOIDY_NN helper", True)
    candidates = sorted(glob.glob(os.path.join(scripts_root, "ploidy_nn_*.sbatch")))
    invocations: list[tuple[str, list[str]]] = []
    for candidate in candidates:
        try:
            commands = logical_program_argvs(
                candidate, "run_ploidy_nn_inference.py")
        except (OSError, UnicodeError, ValueError):
            continue
        for argv in commands:
            invocations.append((candidate, argv))
    if len(invocations) != 1:
        raise ValueError(
            "expected exactly one upstream PLOIDY_NN inference invocation; "
            f"observed {len(invocations)} under: {scripts_root}")
    producer_script, argv = invocations[0]
    lib_range = argv_value(argv, "--lib_range")
    try:
        command_libraries = (
            set(parse_libraries([lib_range])) if lib_range else set())
    except ValueError:
        command_libraries = set()
    expected_prefix = [
        helper, "--h5ad", h5ad, "--weights", weights,
        "--output_dir", output_root, "--lib_range"]
    if (len(argv) != 10 or argv[:8] != expected_prefix or
            argv[8] != lib_range or argv[9] != "--force" or
            command_libraries != selected_libraries):
        raise ValueError(
            "the sole upstream PLOIDY_NN invocation is not the exact forced "
            "all-cohort H5AD/weights/output command")
    newest_script_mtime = os.stat(producer_script).st_mtime_ns
    h5ad_mtime = max(
        os.stat(h5ad).st_mtime_ns, os.stat(contract_path).st_mtime_ns,
        os.stat(weights).st_mtime_ns, os.stat(scaler).st_mtime_ns,
        os.stat(helper).st_mtime_ns)
    calls_validation = []
    for item in paths:
        calls_mtime = os.stat(item.ploidy_nn).st_mtime_ns
        if calls_mtime < h5ad_mtime or calls_mtime < newest_script_mtime:
            raise ValueError(
                f"lib{item.library} PLOIDY_NN calls predate the selected H5AD "
                f"or generation script: {item.ploidy_nn}")
        calls_validation.append(validate_ploidy_calls(
            item.ploidy_nn, item.library, item.expression_barcodes))
    args.runtime_ploidy_output_root = output_root
    args.runtime_ploidy_libraries = sorted(selected_libraries)
    return {
        "status": "PASS",
        "h5ad": h5ad_record,
        "h5ad_contract": contract_record,
        "mex_inputs": validated_mex_records,
        "ploidy_scripts": [provenance_file_record(
            producer_script, "PLOIDY_NN generation script", True)],
        "weights": weight_record,
        "scaler": scaler_record,
        "helper": helper_record,
        "calls_validation": calls_validation,
    }


def validate_and_record_input_bundle(
        args: argparse.Namespace, run: RunPaths,
        paths: Sequence[LibraryPaths], libraries: Sequence[int],
        selected: Sequence[str]) -> str:
    """Fail closed on mixed roots and freeze the matched upstream inputs."""
    mapping_root = args.mapping_input_root
    analysis_root = args.upstream_analysis_root
    identity_root = args.identity_root
    strict = (
        (absolute(mapping_root) != absolute(MAPPING_ROOT) or
         absolute(analysis_root) != absolute(MAPPING_ROOT)) and
        not args.skip_upstream_provenance_check)

    for root, label in (
            (mapping_root, "mapping input root"),
            (analysis_root, "upstream analysis root"),
            (identity_root, "identity reconciliation root")):
        if not os.path.isdir(root):
            raise ValueError(f"{label} is not an existing directory: {root}")

    declared_paths: list[tuple[str, str, str]] = []
    for item in paths:
        for path, label in (
                (item.mapping_bam, "mapping BAM"),
                (item.mapping_bam_index, "mapping BAM index"),
                (item.expression_barcodes, "expression barcodes"),
                (item.expression_features, "expression features"),
                (item.expression_matrix, "expression matrix")):
            declared_paths.append((mapping_root, path, f"lib{item.library} {label}"))
        for path, label in (
                (item.final_assignments, "final reconciled assignments"),
                (item.samples, "DEMUX samples"),
                (item.pileup_sites, "interindividual pileup sites"),
                (item.pileup_molecules, "pileup molecules"),
                (item.pileup_observations, "pileup observations"),
                (item.ambient_standard_prefix + ".contam_rate", "standard ambient rates"),
                (item.ambient_standard_prefix + ".contam_prof", "standard ambient profile"),
                (item.ploidy_nn, "PLOIDY_NN calls")):
            declared_paths.append((analysis_root, path, f"lib{item.library} {label}"))
        if item.cell_groups:
            declared_paths.append((
                analysis_root, item.cell_groups,
                f"lib{item.library} GEX calibration groups"))
    declared_paths.extend((
        (identity_root, args.ledger_input, "final reconciliation ledger"),
        (identity_root, args.identity_validation, "identity validation summary"),
        (identity_root, args.identity_validation_failures,
         "identity validation failures"),
        (identity_root, args.identity_run_summary,
         "identity finalization run summary"),
        (identity_root, args.identity_metadata_manifest,
         "identity metadata manifest"),
        (identity_root, args.identity_expected_genotypes,
         "identity expected genotypes"),
        (identity_root, args.identity_resolution_audit,
         "identity resolution audit"),
        (identity_root, args.identity_metadata_warnings,
         "identity metadata warnings"),
        (identity_root, args.identity_uid_members,
         "identity UID members"),
        (identity_root, args.identity_global_lines,
         "global biological lines"),
        (identity_root, args.identity_global_donors,
         "global donors"),
    ))
    if not args.allow_external_input_templates:
        failures = [f"{label}: {path}"
                    for root, path, label in declared_paths
                    if not lexical_path_contains(root, path)]
        if failures:
            raise ValueError(
                "resolved scientific input escaped its declared root; use "
                "--allow-external-input-templates only for an audited advanced "
                "override: " + "; ".join(failures[:3]))
        if not lexical_path_contains(analysis_root, identity_root):
            raise ValueError(
                "--identity-root is outside --upstream-analysis-root; use "
                "--allow-external-input-templates only for an audited override")

    records: list[dict[str, object]] = []
    # These are source inputs even though BAM itself is not rescanned by the
    # downstream workflow. Its stat lineage proves which remap generated the
    # selected pileup sidecars.
    for item in paths:
        records.extend((
            provenance_file_record(item.mapping_bam,
                                   f"lib{item.library} mapping BAM"),
            provenance_file_record(item.mapping_bam_index,
                                   f"lib{item.library} mapping BAM index"),
            provenance_file_record(item.samples,
                                   f"lib{item.library} DEMUX samples", True),
            provenance_file_record(item.pileup_sites,
                                   f"lib{item.library} pileup sites"),
            provenance_file_record(item.pileup_molecules,
                                   f"lib{item.library} pileup molecules"),
            provenance_file_record(item.pileup_observations,
                                   f"lib{item.library} pileup observations"),
            provenance_file_record(item.demux_prefix + ".assignments",
                                   f"lib{item.library} DEMUX assignments"),
            provenance_file_record(item.final_assignments,
                                   f"lib{item.library} final assignments", True),
            provenance_file_record(item.ploidy_nn,
                                   f"lib{item.library} PLOIDY_NN calls", True),
        ))
        for path, label in (
                (item.ambient_standard_prefix + ".contam_rate", "standard ambient rates"),
                (item.ambient_standard_prefix + ".contam_prof", "standard ambient profile"),
                (item.expression_barcodes, "filtered MEX barcodes"),
                (item.expression_features, "filtered MEX features"),
                (item.expression_matrix, "filtered MEX matrix")):
            records.append(provenance_file_record(
                path, f"lib{item.library} {label}"))
        output_libraries = set(getattr(
            args, "downstream_output_libraries", libraries))
        if ("PREPARE" in selected and item.library in output_libraries and
                not item.cell_groups):
            if not args.allow_library_calibration_fallback:
                raise ValueError(
                    f"lib{item.library} has no unique newest-run GEX cluster "
                    "table; run upstream GEX_AMBIENT or explicitly use "
                    "--allow-library-calibration-fallback")
        elif item.cell_groups:
            records.append(provenance_file_record(
                item.cell_groups, f"lib{item.library} GEX calibration groups", True))

    records.extend((
        provenance_file_record(args.ledger_input,
                               "final reconciliation ledger", True),
        provenance_file_record(args.identity_validation,
                               "identity validation summary", True),
        provenance_file_record(args.identity_validation_failures,
                               "identity validation failures", True),
        provenance_file_record(args.identity_run_summary,
                               "identity finalization run summary", True),
        provenance_file_record(args.identity_metadata_manifest,
                               "identity metadata manifest", True),
        provenance_file_record(args.identity_expected_genotypes,
                               "identity expected genotypes", True),
        provenance_file_record(args.identity_resolution_audit,
                               "identity resolution audit", True),
        provenance_file_record(args.identity_metadata_warnings,
                               "identity metadata warnings", True),
        provenance_file_record(args.identity_uid_members,
                               "identity UID members", True),
        provenance_file_record(args.identity_global_lines,
                               "global biological lines", True),
        provenance_file_record(args.identity_global_donors,
                               "global donors", True),
        provenance_file_record(args.identity_metadata_workbook,
                               "identity metadata workbook", True),
        provenance_file_record(args.demux_pool_workbook,
                               "DEMUX donor-pool workbook", True),
        provenance_file_record(args.expected_pool_metadata,
                               "expected-pool metadata", True),
        provenance_file_record(args.panel_metadata, "panel metadata", True),
        provenance_file_record(args.interindividual_panel,
                               "main interindividual panel", True),
        provenance_file_record(args.het_panel, "HET diagnostic panel", True),
        provenance_file_record(args.species_panel, "species count panel", True),
    ))

    event_libraries = ledger_event_libraries(args.ledger_input, libraries)
    assignment_counts = validate_final_assignment_bundle(
        args.ledger_input, paths)
    finalization_detail = validate_finalization_summary(
        args.identity_run_summary, libraries, assignment_counts)
    identity_metadata_detail = validate_identity_metadata_manifest(
        args.identity_metadata_manifest, args.identity_metadata_workbook,
        args.panel_metadata, libraries)
    donor_pool_detail = validate_donor_pool_chain(
        args, paths, libraries)
    final_uid_detail = validate_final_ledger_uids(
        args, paths, libraries)
    newest_selected_demux = max(
        os.stat(path).st_mtime_ns for item in paths for path in (
            item.pileup_sites, item.pileup_molecules,
            item.pileup_observations, item.demux_prefix + ".assignments"))
    if os.stat(args.ledger_input).st_mtime_ns < newest_selected_demux:
        raise ValueError(
            "canonical final reconciliation ledger predates a selected DEMUX "
            "artifact; the roots do not prove one completed generation")
    ledger_mtime = os.stat(args.ledger_input).st_mtime_ns
    summary_mtime = os.stat(args.identity_run_summary).st_mtime_ns
    metadata_mtime = os.stat(args.identity_metadata_manifest).st_mtime_ns
    validation_mtime = os.stat(args.identity_validation).st_mtime_ns
    validation_failures_mtime = os.stat(
        args.identity_validation_failures).st_mtime_ns
    validation_boundary_mtime = max(
        validation_mtime, validation_failures_mtime)
    newest_identity_metadata_input = max(
        os.stat(args.identity_metadata_workbook).st_mtime_ns,
        os.stat(args.panel_metadata).st_mtime_ns,
        os.stat(args.identity_expected_genotypes).st_mtime_ns,
        os.stat(args.identity_resolution_audit).st_mtime_ns,
        os.stat(args.identity_metadata_warnings).st_mtime_ns,
        os.stat(args.identity_uid_members).st_mtime_ns,
        os.stat(args.identity_global_lines).st_mtime_ns,
        os.stat(args.identity_global_donors).st_mtime_ns,
    )
    if metadata_mtime < newest_identity_metadata_input:
        raise ValueError(
            "identity metadata manifest predates its workbook, panel metadata, "
            "or generated metadata tables")
    if newest_selected_demux < os.stat(args.demux_pool_workbook).st_mtime_ns:
        raise ValueError(
            "selected DEMUX artifacts predate the declared donor-pool workbook")
    newest_posthoc_script = max(
        os.stat(os.path.join(
            args.upstream_analysis_root, "aggregate_library_analysis",
            "slurm_scripts", f"posthoc_lib{item.library}.sbatch")).st_mtime_ns
        for item in paths)
    if validation_mtime < max(
            os.stat(args.expected_pool_metadata).st_mtime_ns,
            newest_posthoc_script):
        raise ValueError(
            "identity validation predates its expected-pool metadata or "
            "POSTHOC generation scripts")
    if not (metadata_mtime <= validation_mtime and
            newest_selected_demux <= validation_mtime and
            validation_boundary_mtime <= ledger_mtime <=
            summary_mtime):
        raise ValueError(
            "DEMUX and identity metadata do not both precede validation, "
            "the final ledger, and the Phase-3 run summary")
    for item in paths:
        if os.stat(item.final_assignments).st_mtime_ns < summary_mtime:
            raise ValueError(
                f"lib{item.library} final assignment predates the Phase-3 "
                "finalization run summary")
    for item in paths:
        newest_library_demux = max(os.stat(path).st_mtime_ns for path in (
            item.pileup_sites, item.pileup_molecules,
            item.pileup_observations, item.demux_prefix + ".assignments"))
        for ambient_path in (
                item.ambient_standard_prefix + ".contam_rate",
                item.ambient_standard_prefix + ".contam_prof"):
            if os.stat(ambient_path).st_mtime_ns < newest_library_demux:
                raise ValueError(
                    f"lib{item.library} standard ambient artifact predates "
                    f"the forced DEMUX bundle: {ambient_path}")
        if item.cell_groups:
            newest_mex = max(os.stat(path).st_mtime_ns for path in (
                item.expression_barcodes, item.expression_features,
                item.expression_matrix))
            if os.stat(item.cell_groups).st_mtime_ns < newest_mex:
                raise ValueError(
                    f"lib{item.library} GEX calibration groups predate the "
                    "selected filtered MEX")
    identity_failure_rows = strict_posthoc_tsv(
        args.identity_validation_failures,
        ("check", "library", "barcode", "detail"),
        "identity reconciliation validation failures")
    if not args.skip_identity_validation and identity_failure_rows:
        raise ValueError(
            "identity validation_failures.tsv is not empty: "
            f"{args.identity_validation_failures}")
    identity_ok, identity_detail = identity_validation_summary_passes(
        args.identity_validation, libraries)
    if not args.skip_identity_validation and not identity_ok:
        raise ValueError(
            "identity validation boundary is not PASS: "
            f"{args.identity_validation} ({identity_detail})")
    for item in paths:
        if item.library not in event_libraries:
            continue
        for prefix, arm in (
                (item.ambient_arm_a_prefix, "A"),
                (item.ambient_arm_c_prefix, "C")):
            arm_paths = (prefix + ".contam_rate", prefix + ".contam_prof")
            records.extend((
                provenance_file_record(
                    arm_paths[0],
                    f"lib{item.library} reconciliation ambient Arm {arm} rates"),
                provenance_file_record(
                    arm_paths[1],
                    f"lib{item.library} reconciliation ambient Arm {arm} profile"),
            ))
            newest_library_demux = max(os.stat(path).st_mtime_ns for path in (
                item.pileup_sites, item.pileup_molecules,
                item.pileup_observations, item.demux_prefix + ".assignments"))
            if any(os.stat(path).st_mtime_ns < newest_library_demux
                   for path in arm_paths):
                raise ValueError(
                    f"lib{item.library} reconciliation ambient Arm {arm} "
                    "predates the forced DEMUX bundle")

    mapping_records = []
    mapping_detail = "not required for historical production layout"
    demux_generation: dict[str, object] = {"status": "NOT_REQUIRED"}
    ploidy_generation: dict[str, object] = {"status": "NOT_REQUIRED"}
    posthoc_generation: dict[str, object] = {"status": "NOT_REQUIRED"}
    identity_metadata_generation: dict[str, object] = {"status": "NOT_REQUIRED"}
    empty_drop_generation: dict[str, object] = {"status": "NOT_REQUIRED"}
    ambient_generation: dict[str, object] = {"status": "NOT_REQUIRED"}
    if strict:
        if not args.mapping_run_root:
            raise ValueError(
                "a non-production mapping input requires --mapping-run-root "
                "unless it has the canonical <run>/rna3/mapping_output layout")
        expected_mapping_root = os.path.join(
            args.mapping_run_root, "rna3", "mapping_output")
        if absolute(expected_mapping_root) != absolute(mapping_root):
            raise ValueError(
                "--mapping-input-root is not the rna3/mapping_output tree below "
                f"--mapping-run-root: {expected_mapping_root}")
        required_mapping_records = (
            ("control/run_config.json", "mapping run configuration"),
            ("control/job_plan.json", "mapping job plan"),
            ("rna3/libs.txt", "mapping library manifest"),
            ("rna3/mapping_project/MAPPING_COMPLETE.ok", "RNA mapping completion marker"),
            ("rna3/mapping_project/rg_metadata.tsv", "RNA read-group metadata"),
            ("rna3/mapping_project/run_assignments.tsv", "RNA run assignments"),
            ("validation/output_validation.tsv", "mapping output validation"),
            ("validation/RUN_COMPLETE.ok", "mapping run completion marker"),
        )
        for relative, label in required_mapping_records:
            mapping_records.append(provenance_file_record(
                os.path.join(args.mapping_run_root, relative), label,
                relative.endswith((".json", ".ok"))))
        newest_mapping_input = max(
            os.stat(value).st_mtime_ns for item in paths for value in (
                item.mapping_bam, item.mapping_bam_index,
                item.expression_barcodes, item.expression_features,
                item.expression_matrix))
        for relative in (
                "rna3/mapping_project/MAPPING_COMPLETE.ok",
                "validation/output_validation.tsv",
                "validation/RUN_COMPLETE.ok"):
            marker = os.path.join(args.mapping_run_root, relative)
            if os.stat(marker).st_mtime_ns < newest_mapping_input:
                raise ValueError(
                    f"mapping completion/validation record predates a selected "
                    f"BAM or MEX input: {marker}")
        validation_path = os.path.join(
            args.mapping_run_root, "validation", "output_validation.tsv")
        required_validated_rows = [
            (value, "rna3", f"{LIBRARY_PREFIX}{item.library}")
            for item in paths for value in (
                item.mapping_bam, item.mapping_bam_index,
                item.expression_barcodes, item.expression_features,
                item.expression_matrix)]
        mapping_ok, mapping_detail = mapping_validation_passes(
            validation_path, required_validated_rows)
        if not mapping_ok:
            raise ValueError(
                f"newest mapping validation does not cover selected inputs: "
                f"{validation_path} ({mapping_detail})")
        demux_generation = verify_demux_generation(args, paths)
        empty_drop_generation = verify_empty_drop_generation(
            args, paths, validation_boundary_mtime)
        ambient_generation = verify_ambient_generation(
            args, paths, event_libraries, ledger_mtime,
            validation_boundary_mtime)
        ploidy_generation = verify_ploidy_generation(args, paths, libraries)
        posthoc_generation = verify_posthoc_generation(args, paths)
        identity_metadata_generation = verify_identity_metadata_generation(
            args, paths)

    payload = {
        "schema_version": "tetra_arm_input_provenance_v1",
        "release": RELEASE,
        "status": "PASS" if not args.skip_upstream_provenance_check
                  else "PASS_WITH_AUDIT_OVERRIDE",
        "libraries": list(libraries),
        "downstream_output_libraries": list(
            getattr(args, "downstream_output_libraries", libraries)),
        "upstream_cohort_libraries": list(libraries),
        "roots": {
            "mapping_run_root": args.mapping_run_root or None,
            "mapping_input_root": mapping_root,
            "upstream_analysis_root": analysis_root,
            "identity_root": identity_root,
            "downstream_run_root": run.root,
        },
        "panel_roles": {
            "ASE_DONOR_GENOTYPES": args.interindividual_panel,
            "PLOIDY_DIAGNOSTICS_ONLY": args.het_panel,
            "SPECIES_COUNTING_ONLY": args.species_panel,
        },
        "mapping_validation": {
            "status": "PASS" if strict else "NOT_REQUIRED",
            "detail": mapping_detail,
            "records": mapping_records,
        },
        "demux_generation": demux_generation,
        "empty_drop_generation": empty_drop_generation,
        "ambient_generation": ambient_generation,
        "posthoc_generation": posthoc_generation,
        "ploidy_generation": ploidy_generation,
        "identity_validation": {
            "status": ("AUDIT_OVERRIDE" if args.skip_identity_validation
                       else "PASS"),
            "detail": identity_detail,
        },
        "identity_finalization": finalization_detail,
        "identity_metadata": identity_metadata_detail,
        "identity_metadata_generation": identity_metadata_generation,
        "donor_pool_chain": donor_pool_detail,
        "final_uid_binding": final_uid_detail,
        "event_libraries_requiring_four_arm_ambient": sorted(event_libraries),
        "input_records": records,
        "audit_overrides": {
            "skip_upstream_provenance_check": bool(
                args.skip_upstream_provenance_check),
            "allow_external_input_templates": bool(
                args.allow_external_input_templates),
            "allow_library_calibration_fallback": bool(
                args.allow_library_calibration_fallback),
        },
    }
    return publish_new_or_identical(
        run.input_provenance,
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        "input provenance contract")


def identity_validation_summary_passes(
        path: str, libraries: Sequence[int]) -> tuple[bool, str]:
    """Require the complete attached-validator checklist to pass."""
    if not regular_nonempty(path):
        return False, "missing or empty"
    try:
        with open(path, "r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            rows = list(reader)
            fields = list(reader.fieldnames or [])
        if not rows:
            return False, "contains no validation rows"
        if fields != ["check", "status", "n_failures", "detail"]:
            return False, "header is not the canonical validation-summary schema"
        per_library_checks = {
            "all_input_barcodes_once",
            "original_assignments_preserved",
            "candidate_global_biological_universe",
            "technical_multiplet_candidate_invariants",
            "applied_change_invariants",
            "component_singlet_classification_invariants",
            "occupancy_ambiguity_invariants",
            "multiplets_excluded_from_single_cell_assignments",
            "reconciled_assignments_invariants",
            "doublet_dragon_context_invariants",
            "uid_resolution_invariants",
            "mt_guardrail",
            "atac_mode_isolation",
        }
        global_checks = {
            "unique_library_barcode",
            "event_mass_threshold",
            "component_singlet_event_invariants",
            "cell_exchange_evidence_invariants",
            "library_exchange_invariants",
            "library_exchange_donor_evidence_invariants",
            "aggregate_row_count",
            "manifest_invariants",
        }
        expected = {
            (check, f"lib{int(library)}")
            for library in libraries for check in per_library_checks}
        expected.update((check, "") for check in global_checks)
        observed: set[tuple[str, str]] = set()
        failed = []
        for index, row in enumerate(rows, start=1):
            check = str(row.get("check", "")).strip()
            detail = str(row.get("detail", "")).strip()
            key = (check, detail)
            if key in observed:
                return False, f"duplicate validation row: {check}/{detail}"
            observed.add(key)
            try:
                failures = int(str(row.get("n_failures", "")).strip())
            except (TypeError, ValueError):
                failures = -1
            if (str(row.get("status", "")).strip() != "PASS"
                    or failures != 0):
                failed.append(check or f"row {index}")
        if failed:
            return False, "failed checks: " + ",".join(
                value or "UNKNOWN" for value in failed[:5])
        if observed != expected:
            missing = sorted(expected - observed)[:5]
            extra = sorted(observed - expected)[:5]
            return False, f"checklist mismatch: missing={missing}, extra={extra}"
        return True, f"PASS ({len(rows)} complete checks)"
    except (OSError, csv.Error, UnicodeError) as exc:
        return False, f"unreadable: {exc}"


def manifest_payload(header: Sequence[str],
                     rows: Iterable[Mapping[str, object]]) -> str:
    lines = ["\t".join(header)]
    for row in rows:
        values = []
        for field in header:
            value = validate_text(str(row.get(field, "")), f"manifest field {field}")
            values.append(value)
        lines.append("\t".join(values))
    return "\n".join(lines) + "\n"


def publish_new_or_identical(path: str, payload: str, label: str) -> str:
    """Atomically create an execution input; never replace different content."""
    destination = absolute(path)
    os.makedirs(os.path.dirname(destination), exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(
        prefix=os.path.basename(destination) + ".tmp.",
        dir=os.path.dirname(destination), text=True)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8", newline="") as handle:
            handle.write(payload)
            handle.flush()
            os.fsync(handle.fileno())
        try:
            os.link(temporary, destination)
        except FileExistsError:
            try:
                with open(destination, "r", encoding="utf-8") as handle:
                    existing = handle.read()
            except OSError as exc:
                raise ValueError(
                    f"existing {label} cannot be verified: "
                    f"{destination}: {exc}") from exc
            if existing != payload:
                raise ValueError(
                    f"refusing to replace a differing {label}: {destination}; "
                    "use a new --run-root so a queued job cannot be retargeted")
        return destination
    finally:
        try:
            os.unlink(temporary)
        except FileNotFoundError:
            pass


def write_task_manifest(path: str, header: Sequence[str],
                        rows: Iterable[Mapping[str, object]]) -> str:
    """Write a minimal array-index map without retargeting a queued job."""
    return publish_new_or_identical(
        path, manifest_payload(header, rows), "array task map")


def run_paths(root: str) -> RunPaths:
    root = absolute(root)
    return RunPaths(
        root=root,
        logs=os.path.join(root, "logs"),
        scripts=os.path.join(root, "generated_scripts"),
        manifests=os.path.join(root, "task_manifests"),
        reference=os.path.join(root, "reference"),
        ledger=os.path.join(root, "ledger"),
        prepare=os.path.join(root, "prepare"),
        ase=os.path.join(root, "ase"),
        expression=os.path.join(root, "expression"),
        call=os.path.join(root, "call"),
        report=os.path.join(root, "report"),
        call_prefix=os.path.join(root, "call", "tetra_arm_cnv"),
    )


def path_contains(parent: str, child: str) -> bool:
    """Return whether child is the same path as parent or lies below it."""
    parent = os.path.realpath(absolute(parent))
    child = os.path.realpath(absolute(child))
    try:
        return os.path.commonpath((parent, child)) == parent
    except ValueError:
        return False


def configure_input_roots(args: argparse.Namespace) -> None:
    """Resolve the mapping-input and upstream-analysis namespaces.

    The production layout historically stores both namespaces below
    ``MAPPING_ROOT``.  A non-production remap must name the isolated analysis
    root produced by orchestrate_tetraploid.py so MEX files cannot be silently
    mixed with demux, ambient, or reconciled-identity products from another run.
    """
    mapping_input_root = absolute(args.mapping_input_root)
    production_mapping_root = absolute(MAPPING_ROOT)
    supplied_analysis_root = args.upstream_analysis_root
    if mapping_input_root != production_mapping_root and not supplied_analysis_root:
        raise ValueError(
            "a non-production --mapping-input-root requires "
            "--upstream-analysis-root so remapped MEX inputs cannot be mixed "
            "with production demux/reconciliation products")

    upstream_analysis_root = absolute(
        supplied_analysis_root or mapping_input_root)
    if mapping_input_root != production_mapping_root and (
            path_contains(mapping_input_root, upstream_analysis_root) or
            path_contains(upstream_analysis_root, mapping_input_root)):
        raise ValueError(
            "--mapping-input-root and --upstream-analysis-root must be "
            "separate, non-nested directories for a non-production remap")

    args.mapping_input_root = mapping_input_root
    args.upstream_analysis_root = upstream_analysis_root
    if args.mapping_run_root:
        args.mapping_run_root = absolute(args.mapping_run_root)
    elif (os.path.basename(mapping_input_root) == "mapping_output" and
          os.path.basename(os.path.dirname(mapping_input_root)) == "rna3"):
        args.mapping_run_root = os.path.dirname(os.path.dirname(mapping_input_root))
    else:
        args.mapping_run_root = ""
    args.run_root = absolute(
        args.run_root or os.path.join(
            upstream_analysis_root, "aggregate_library_analysis",
            "tetra_arm_cnv"))


def default_templates(args: argparse.Namespace) -> None:
    mapping_input_root = absolute(args.mapping_input_root)
    upstream_analysis_root = absolute(args.upstream_analysis_root)
    identity_root = absolute(
        args.identity_root or os.path.join(
            upstream_analysis_root, "aggregate_library_analysis",
            "identity_reconciliation"))
    args.identity_root = identity_root
    mapping_library_root = os.path.join(
        mapping_input_root, f"{LIBRARY_PREFIX}{{lib}}")
    analysis_library_root = os.path.join(
        upstream_analysis_root, f"{LIBRARY_PREFIX}{{lib}}")
    demux_prefix = os.path.join(
        analysis_library_root, DEMUX_SUBDIR, "lib{lib}_demuxed")
    contam_root = os.path.join(
        analysis_library_root, DEMUX_SUBDIR, "contamination",
        args.ambient_condition)
    four_arm_root = os.path.join(
        contam_root, "reconciliation_four_arm", args.ambient_candidate_set)

    args.ledger_input = absolute(args.ledger_input or os.path.join(
        identity_root, "aggregate", "identity_reconciliation_final_cells.tsv.gz"))
    args.identity_validation = absolute(
        args.identity_validation or os.path.join(
            identity_root, "validation", "validation_summary.tsv"))
    args.identity_validation_failures = absolute(
        args.identity_validation_failures or os.path.join(
            os.path.dirname(args.identity_validation),
            "validation_failures.tsv"))
    args.identity_run_summary = absolute(
        args.identity_run_summary or os.path.join(
            identity_root, "aggregate",
            "identity_reconciliation_run_summary.tsv"))
    args.identity_metadata_manifest = absolute(
        args.identity_metadata_manifest or os.path.join(
            identity_root, "metadata", "metadata_manifest.json"))
    args.identity_expected_genotypes = absolute(
        args.identity_expected_genotypes or os.path.join(
            identity_root, "metadata", "library_expected_genotypes.tsv"))
    args.identity_resolution_audit = absolute(
        args.identity_resolution_audit or os.path.join(
            identity_root, "metadata", "library_resolution_audit.tsv"))
    args.identity_metadata_warnings = absolute(
        args.identity_metadata_warnings or os.path.join(
            identity_root, "metadata", "metadata_warnings.tsv"))
    args.identity_uid_members = absolute(
        args.identity_uid_members or os.path.join(
            identity_root, "metadata", "library_uid_members.tsv"))
    args.identity_global_lines = absolute(
        args.identity_global_lines or os.path.join(
            identity_root, "metadata", "global_biological_lines.tsv"))
    args.identity_global_donors = absolute(
        args.identity_global_donors or os.path.join(
            identity_root, "metadata", "global_donors.tsv"))
    args.final_assignments_template = args.final_assignments_template or os.path.join(
        identity_root, "final_assignments", "lib{lib}.reconciled.assignments")
    args.mapping_bam_template = args.mapping_bam_template or os.path.join(
        mapping_library_root, "gex.bam")
    args.demux_prefix_template = args.demux_prefix_template or demux_prefix
    args.expression_barcodes_template = (
        args.expression_barcodes_template
        or os.path.join(
            mapping_library_root, "filtered", "barcodes.tsv.gz"))
    args.expression_features_template = (
        args.expression_features_template
        or os.path.join(
            mapping_library_root, "filtered", "features.tsv.gz"))
    args.expression_matrix_template = (
        args.expression_matrix_template
        or os.path.join(
            mapping_library_root, "filtered", "matrix.mtx.gz"))
    args.ambient_standard_template = (
        args.ambient_standard_template
        or os.path.join(contam_root, "lib{lib}_demuxed"))
    args.ambient_arm_a_template = (
        args.ambient_arm_a_template
        or os.path.join(four_arm_root, "demux_original", "lib{lib}_demuxed"))
    args.ambient_arm_c_template = (
        args.ambient_arm_c_template
        or os.path.join(four_arm_root, "reconciled_augmented", "lib{lib}_demuxed"))
    args.ploidy_nn_template = (
        args.ploidy_nn_template
        or os.path.join(
            upstream_analysis_root, "aggregate_library_analysis", "ploidy",
            "lib{lib}.ploidy_calls_nn.tsv"))


def resolve_cell_groups(args: argparse.Namespace, upstream_analysis_root: str,
                        library: int) -> tuple[str, str]:
    if args.cell_groups_template:
        return (expand_template(
            args.cell_groups_template, library, "cell groups template"),
                "EXPLICIT_TEMPLATE")
    cluster_dir = os.path.join(
        upstream_analysis_root, "aggregate_library_analysis", "gex_ambient",
        args.gex_ambient_analysis, "clusters")
    pattern = os.path.join(cluster_dir, f"lib{library}.*.tsv")
    matches = sorted(
        absolute(path) for path in glob.glob(pattern)
        if os.path.isfile(path) and not path.endswith(".qc.tsv"))
    if len(matches) > 1:
        raise ValueError(
            f"lib{library} has multiple GEX cluster matches under {cluster_dir}: "
            f"{', '.join(matches)}; set --cell-groups-template explicitly")
    if matches:
        return matches[0], "AUTO_DISCOVERED"
    return "", "LIBRARY_FALLBACK"


def make_library_paths(args: argparse.Namespace, run: RunPaths,
                       libraries: Sequence[int],
                       selected: Sequence[str]) -> list[LibraryPaths]:
    result = []
    upstream_analysis_root = absolute(args.upstream_analysis_root)
    for library in libraries:
        demux = expand_template(
            args.demux_prefix_template, library, "demux prefix template")
        if "PREPARE" in selected:
            cell_groups, cell_groups_source = resolve_cell_groups(
                args, upstream_analysis_root, library)
        else:
            cell_groups, cell_groups_source = "", "NOT_SELECTED"
        result.append(LibraryPaths(
            library=library,
            mapping_bam=expand_template(
                args.mapping_bam_template, library, "mapping BAM template"),
            mapping_bam_index=expand_template(
                args.mapping_bam_template, library, "mapping BAM template") + ".bai",
            final_assignments=expand_template(
                args.final_assignments_template, library,
                "final assignments template"),
            demux_prefix=demux,
            samples=demux + ".samples",
            pileup_sites=demux + ".pileup_sites.tsv.gz",
            pileup_molecules=demux + ".pileup_molecules.tsv.gz",
            pileup_observations=demux + ".pileup_obs.tsv.gz",
            expression_barcodes=expand_template(
                args.expression_barcodes_template, library,
                "expression barcodes template"),
            expression_features=expand_template(
                args.expression_features_template, library,
                "expression features template"),
            expression_matrix=expand_template(
                args.expression_matrix_template, library,
                "expression matrix template"),
            ambient_standard_prefix=expand_template(
                args.ambient_standard_template, library,
                "ambient standard template"),
            ambient_arm_a_prefix=expand_template(
                args.ambient_arm_a_template, library,
                "ambient Arm A template"),
            ambient_arm_c_prefix=expand_template(
                args.ambient_arm_c_template, library,
                "ambient Arm C template"),
            ploidy_nn=expand_template(
                args.ploidy_nn_template, library, "ploidy NN template"),
            cell_groups=cell_groups,
            cell_groups_source=cell_groups_source,
            split_ledger=os.path.join(
                run.ledger, f"lib{library}.final_cells.tsv.gz"),
            cell_manifest=os.path.join(
                run.prepare, f"lib{library}.cell_manifest.tsv.gz"),
            ambient_sources=os.path.join(
                run.prepare, f"lib{library}.ambient_sources.tsv.gz"),
            prepare_qc=os.path.join(
                run.prepare, f"lib{library}.prepare_qc.tsv"),
            prepare_contract=os.path.join(
                run.prepare, f"lib{library}.prepare_contract.json"),
            ase=os.path.join(run.ase, f"lib{library}.arm_ase.tsv.gz"),
            ase_qc=os.path.join(run.ase, f"lib{library}.ase_qc.tsv"),
            expression=os.path.join(
                run.expression, f"lib{library}.arm_expression.tsv.gz"),
            expression_qc=os.path.join(
                run.expression, f"lib{library}.expression_qc.tsv"),
            expression_contract=os.path.join(
                run.expression, f"lib{library}.expression_contract.json"),
        ))
    return result


def module_block(modules: Sequence[str], imports: Sequence[str] = ()) -> str:
    loads = "\n".join(f"module load {module}" for module in modules)
    import_block = ""
    if imports:
        import_block = (
            "\ncommand -v python3 >/dev/null 2>&1"
            "\npython3 - <<'PY'\n" +
            "\n".join(f"import {name}" for name in imports) +
            "\nPY")
    return f"""module purge
{loads}
module list 2>&1
command -v date >/dev/null 2>&1
command -v hostname >/dev/null 2>&1{import_block}

date
hostname"""


def input_provenance_guard(run: RunPaths) -> str:
    """Re-stat frozen upstream inputs inside every submitted stage."""
    if not regular_nonempty(run.input_provenance):
        return ""
    return f"""
INPUT_PROVENANCE={shlex.quote(run.input_provenance)}
if [[ -s "$INPUT_PROVENANCE" ]]; then
    python3 - "$INPUT_PROVENANCE" <<'PY'
import hashlib
import json
import os
import sys

path = sys.argv[1]
with open(path, "r", encoding="utf-8") as handle:
    payload = json.load(handle)
if payload.get("schema_version") != "tetra_arm_input_provenance_v1":
    raise SystemExit(f"ERROR: invalid input provenance schema: {{path}}")

records = []
def visit(value):
    if isinstance(value, dict):
        required = {{"path", "realpath", "size", "mtime_ns"}}
        if required <= set(value):
            records.append(value)
        for child in value.values():
            visit(child)
    elif isinstance(value, list):
        for child in value:
            visit(child)
visit(payload)
if not records:
    raise SystemExit(f"ERROR: input provenance has no file records: {{path}}")

checked = set()
for record in records:
    current_path = os.path.abspath(str(record["path"]))
    signature = (
        current_path, str(record["realpath"]), int(record["size"]),
        int(record["mtime_ns"]), str(record.get("sha256", "")))
    if signature in checked:
        continue
    checked.add(signature)
    try:
        stat_result = os.stat(current_path)
    except OSError as exc:
        raise SystemExit(
            f"ERROR: frozen upstream input is unavailable: {{current_path}}: {{exc}}")
    actual = (
        os.path.realpath(current_path), stat_result.st_size,
        stat_result.st_mtime_ns)
    expected = signature[1:4]
    if actual != expected:
        raise SystemExit(
            f"ERROR: frozen upstream input changed after planning: {{current_path}}")
    expected_sha = signature[4]
    if expected_sha:
        digest = hashlib.sha256()
        with open(current_path, "rb") as source:
            for block in iter(lambda: source.read(1024 * 1024), b""):
                digest.update(block)
        if digest.hexdigest() != expected_sha:
            raise SystemExit(
                f"ERROR: frozen upstream input hash changed: {{current_path}}")
print(f"PASS: {{len(checked)}} frozen upstream input records")
PY
else
    echo "ERROR: missing input provenance contract: $INPUT_PROVENANCE" >&2
    exit 1
fi
"""


def sbatch_header(stage: str, run: RunPaths, args: argparse.Namespace,
                  cpus: int, memory: str, modules: Sequence[str],
                  imports: Sequence[str] = (), tasks: int | None = None) -> str:
    stage_lower = stage.lower()
    log_token = "%A_%a" if tasks is not None else "%j"
    array = ""
    if tasks is not None:
        if tasks < 1:
            raise ValueError(f"cannot render empty {stage} array")
        suffix = (f"%{args.array_throttle}"
                  if args.array_throttle is not None else "")
        array = f"#SBATCH --array=0-{tasks - 1}{suffix}\n"
    return f"""#!/bin/bash
# Generated deterministically by orchestrate_tetra_arm_cnv.py {RELEASE}.
# Generated script and task map are immutable execution inputs.
# Change orchestrator options and use a new --run-root instead of editing them.
#SBATCH --job-name=tetarm_{stage_lower}
#SBATCH --output={run.logs}/tetarm_{stage_lower}_{log_token}.out
#SBATCH --error={run.logs}/tetarm_{stage_lower}_{log_token}.err
#SBATCH --partition={args.partition}
#SBATCH --nodes=1
#SBATCH --cpus-per-task={cpus}
#SBATCH --mem={memory}
#SBATCH --time={args.time}
{array}
set -euo pipefail

{module_block(modules, imports)}
{input_provenance_guard(run)}
"""


def reference_script(args: argparse.Namespace, run: RunPaths) -> str:
    modules = ("miniforge/3", "hal/latest") \
        if args.reference_mode == "HAL_LIFTOVER" else BASE_PYTHON_MODULES
    script = sbatch_header(
        "REFERENCE", run, args, 2, args.reference_memory,
        modules, ("json",))
    output_checks = [run.reference_qc, run.reference_contract]
    if args.reference_mode != "EXPLICIT_BED":
        output_checks.insert(0, args.arms_bed)
    checks = "\n".join(
        (f"if [[ -e {shlex.quote(path)} || -L {shlex.quote(path)} ]]; then\n"
         f"    echo {shlex.quote('ERROR: refusing to overwrite REFERENCE output: ' + path + '; use a new --run-root')} >&2\n"
         "    exit 1\nfi")
        for path in output_checks)
    include_sex = "\ncommand+=(--include-sex-chromosomes)" \
        if args.include_sex_chromosomes else ""
    if args.reference_mode == "GENE_SYNTENY":
        body = f"""
ARM_BUILDER={shlex.quote(args.arm_builder_script)}
ANNOTATION={shlex.quote(args.gene_annotation)}
GENE_ARMS={shlex.quote(args.gene_arms)}
REFERENCE_FAI={shlex.quote(args.reference_fai)}
ARMS_BED={shlex.quote(args.arms_bed)}
REFERENCE_QC={shlex.quote(run.reference_qc)}
REFERENCE_CONTRACT={shlex.quote(run.reference_contract)}

test -s "$ARM_BUILDER"
test -s "$ANNOTATION"
test -s "$GENE_ARMS"
if [[ {shlex.quote(args.gene_projection_mode)} == anchored-contigs ]]; then
    test -s "$REFERENCE_FAI"
fi
python3 "$ARM_BUILDER" --version
{checks}
command=(
    python3 "$ARM_BUILDER"
    --annotation "$ANNOTATION"
    --gene-arms "$GENE_ARMS"
    --projection-mode {shlex.quote(args.gene_projection_mode)}
    --output "$ARMS_BED"
    --qc "$REFERENCE_QC"
    --contract "$REFERENCE_CONTRACT"
)
if [[ {shlex.quote(args.gene_projection_mode)} == anchored-contigs ]]; then
    command+=(--reference-fai "$REFERENCE_FAI")
fi
{include_sex}
"${{command[@]}}"
test -s "$ARMS_BED"
test -s "$REFERENCE_QC"
test -s "$REFERENCE_CONTRACT"
echo "COMPLETE: REFERENCE gene-synteny arm projection"
date
"""
    elif args.reference_mode == "HAL_LIFTOVER":
        body = f"""
HAL_REFERENCE={shlex.quote(args.hal_reference_script)}
HAL_FILE={shlex.quote(args.hal_file)}
SOURCE_ARMS={shlex.quote(args.source_arms_bed)}
ARMS_BED={shlex.quote(args.arms_bed)}
REFERENCE_QC={shlex.quote(run.reference_qc)}
REFERENCE_CONTRACT={shlex.quote(run.reference_contract)}

test -s "$HAL_REFERENCE"
test -s "$HAL_FILE"
test -s "$SOURCE_ARMS"
command -v halStats >/dev/null 2>&1
command -v halLiftover >/dev/null 2>&1
python3 "$HAL_REFERENCE" --version
{checks}
command=(
    python3 "$HAL_REFERENCE"
    --hal "$HAL_FILE"
    --source-genome {shlex.quote(args.hal_source_genome)}
    --target-genome {shlex.quote(args.hal_target_genome)}
    --source-arms "$SOURCE_ARMS"
    --output "$ARMS_BED"
    --qc "$REFERENCE_QC"
    --contract "$REFERENCE_CONTRACT"
)
{include_sex}
"${{command[@]}}"
test -s "$ARMS_BED"
test -s "$REFERENCE_QC"
test -s "$REFERENCE_CONTRACT"
echo "COMPLETE: REFERENCE HAL chromosome-arm projection"
date
"""
    else:
        body = f"""
ARMS_BED={shlex.quote(args.arms_bed)}
REFERENCE_QC={shlex.quote(run.reference_qc)}
REFERENCE_CONTRACT={shlex.quote(run.reference_contract)}

test -s "$ARMS_BED"
{checks}
python3 - "$ARMS_BED" "$REFERENCE_QC" "$REFERENCE_CONTRACT" <<'PY'
import json
import os
import sys

arms_path, qc_path, contract_path = sys.argv[1:]
intervals = 0
names = set()
contigs = set()
with open(arms_path, "r", encoding="utf-8") as handle:
    for line_number, line in enumerate(handle, start=1):
        if not line.strip() or line.startswith("#"):
            continue
        fields = line.rstrip("\\r\\n").split("\\t")
        if len(fields) < 4:
            raise SystemExit(f"{{arms_path}}:{{line_number}}: expected BED4")
        try:
            start, end = int(fields[1]), int(fields[2])
        except ValueError:
            raise SystemExit(f"{{arms_path}}:{{line_number}}: noninteger coordinates")
        if not fields[0] or not fields[3] or start < 0 or end <= start:
            raise SystemExit(f"{{arms_path}}:{{line_number}}: invalid BED4 interval")
        intervals += 1
        contigs.add(fields[0])
        names.add(fields[3])
if intervals == 0:
    raise SystemExit(f"{{arms_path}}: no BED intervals")
with open(qc_path, "x", encoding="utf-8", newline="") as handle:
    handle.write("mode\\tintervals\\tcontigs\\tlogical_arms\\tstatus\\tschema_version\\n")
    handle.write(f"EXPLICIT_BED\\t{{intervals}}\\t{{len(contigs)}}\\t{{len(names)}}\\tPASS\\ttetra_arm_reference_qc_v1\\n")
with open(contract_path, "x", encoding="utf-8") as handle:
    json.dump({{
        "schema_version": "tetra_arm_reference_contract_v1",
        "mode": "EXPLICIT_BED", "arms_bed": os.path.abspath(arms_path),
        "intervals": intervals, "contigs": len(contigs),
        "logical_arms": sorted(names), "status": "PASS",
    }}, handle, sort_keys=True, indent=2)
    handle.write("\\n")
PY
test -s "$REFERENCE_QC"
test -s "$REFERENCE_CONTRACT"
echo "COMPLETE: REFERENCE explicit BED validation"
date
"""
    return script + body


def ledger_script(args: argparse.Namespace, run: RunPaths,
                  libraries: Sequence[int]) -> str:
    runtime_ploidy_root = str(getattr(
        args, "runtime_ploidy_output_root", ""))
    runtime_ploidy_libraries = list(getattr(
        args, "runtime_ploidy_libraries", []))
    runtime_ploidy = bool(runtime_ploidy_root and runtime_ploidy_libraries)
    modules = BASE_PYTHON_MODULES + ASE_MODULES
    if runtime_ploidy:
        modules = (BASE_PYTHON_MODULES + (args.ploidy_nn_module,) +
                   ASE_MODULES)
    script = sbatch_header(
        "LEDGER", run, args, args.ledger_cpus, args.ledger_memory,
        modules, ("csv", "gzip", "json"))
    library_words = " ".join(shlex.quote(str(value)) for value in libraries)
    ploidy_runtime_block = ""
    if runtime_ploidy:
        ploidy_runtime_block = f"""
PLOIDY_HELPER={shlex.quote(args.ploidy_nn_helper)}
PLOIDY_H5AD={shlex.quote(args.ploidy_input_h5ad)}
PLOIDY_WEIGHTS={shlex.quote(args.ploidy_nn_weights)}
PLOIDY_SCALER={shlex.quote(args.ploidy_nn_weights[:-3] + '_scaler.npz' if args.ploidy_nn_weights.endswith('.pt') else args.ploidy_nn_weights + '_scaler.npz')}
PLOIDY_SELECTED_ROOT={shlex.quote(runtime_ploidy_root)}
PLOIDY_LIBRARIES_CSV={shlex.quote(','.join(str(value) for value in runtime_ploidy_libraries))}
test -s "$PLOIDY_HELPER"
test -s "$PLOIDY_H5AD"
test -s "$PLOIDY_WEIGHTS"
test -s "$PLOIDY_SCALER"
mkdir -p "$panel_tmp_root/ploidy"
python3 "$PLOIDY_HELPER" \
    --h5ad "$PLOIDY_H5AD" \
    --weights "$PLOIDY_WEIGHTS" \
    --output_dir "$panel_tmp_root/ploidy" \
    --lib_range "$PLOIDY_LIBRARIES_CSV" \
    --force
python3 - "$PLOIDY_SELECTED_ROOT" "$panel_tmp_root/ploidy" \
    "$PLOIDY_LIBRARIES_CSV" <<'PY'
import csv
import math
import os
import sys

selected_root, reproduced_root, library_csv = sys.argv[1:]
libraries = [int(value) for value in library_csv.split(",") if value]
header = [
    "barcode", "library", "ploidy_call", "ploidy_probability",
    "classification_group", "prob_tetraploid", "qc_pass",
]

def load(path, library):
    if not os.path.isfile(path) or os.path.getsize(path) <= 0:
        raise ValueError(f"missing PLOIDY_NN output: {{path}}")
    with open(path, "r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if list(reader.fieldnames or []) != header:
            raise ValueError(f"PLOIDY_NN header mismatch: {{path}}")
        rows = {{}}
        for row in reader:
            barcode = str(row.get("barcode", "")).strip()
            key = (str(row.get("library", "")).strip(), barcode)
            if (key[0] != str(library) or not barcode or key in rows):
                raise ValueError(f"PLOIDY_NN key mismatch: {{path}}")
            rows[key] = row
    if not rows:
        raise ValueError(f"header-only PLOIDY_NN output: {{path}}")
    return rows


for library in libraries:
    selected = load(os.path.join(
        selected_root, f"lib{{library}}.ploidy_calls_nn.tsv"), library)
    reproduced = load(os.path.join(
        reproduced_root, f"lib{{library}}.ploidy_calls_nn.tsv"), library)
    if set(selected) != set(reproduced):
        raise ValueError(f"lib{{library}} PLOIDY_NN cell universe changed")
    for key in selected:
        left, right = selected[key], reproduced[key]
        for field in ("barcode", "library", "ploidy_call",
                      "classification_group", "qc_pass"):
            if left[field] != right[field]:
                raise ValueError(
                    f"lib{{library}}/{{key[1]}} PLOIDY_NN {{field}} changed")
        for field in ("ploidy_probability", "prob_tetraploid"):
            try:
                observed = float(left[field])
                expected = float(right[field])
            except ValueError as exc:
                raise ValueError(
                    f"lib{{library}}/{{key[1]}} invalid PLOIDY_NN {{field}}") \
                    from exc
            if (not math.isfinite(observed) or not math.isfinite(expected) or
                    not math.isclose(observed, expected, rel_tol=0.0,
                                     abs_tol=1e-7)):
                raise ValueError(
                    f"lib{{library}}/{{key[1]}} PLOIDY_NN {{field}} does not "
                    "reproduce from the declared H5AD/model")
PY
"""
    output_checks = "\n".join(
        (f"if [[ -e {shlex.quote(path)} || -L {shlex.quote(path)} ]]; then\n"
         f"    echo {shlex.quote('ERROR: refusing to overwrite LEDGER output: ' + path + '; use a new --run-root')} >&2\n"
         "    exit 1\nfi")
        for path in (
            *(os.path.join(run.ledger, f"lib{library}.final_cells.tsv.gz")
              for library in libraries),
            os.path.join(run.ledger, "split_ledger_summary.tsv"),
            os.path.join(run.ledger, "split_ledger_contract.json"),
        ))
    return script + f"""
command -v python3 >/dev/null 2>&1
command -v cmp >/dev/null 2>&1
command -v mktemp >/dev/null 2>&1
PREPARE_SCRIPT={shlex.quote(args.prepare_script)}
LEDGER_INPUT={shlex.quote(args.ledger_input)}
IDENTITY_VALIDATION={shlex.quote(args.identity_validation)}
PANEL_UTILITY={shlex.quote(args.panel_distinguishability_binary)}
INTERINDIVIDUAL_PANEL={shlex.quote(args.interindividual_panel)}
DISTINGUISHABILITY={shlex.quote(os.path.join(os.path.dirname(args.identity_metadata_manifest), 'nuclear_panel_distinguishability.tsv'))}
SKIP_IDENTITY_VALIDATION={"1" if args.skip_identity_validation else "0"}
LEDGER_ROOT={shlex.quote(run.ledger)}
LIBRARY_CSV={shlex.quote(','.join(str(value) for value in libraries))}

test -s "$PREPARE_SCRIPT"
python3 "$PREPARE_SCRIPT" --version
test -x "$PANEL_UTILITY"
test -s "$INTERINDIVIDUAL_PANEL"
test -s "$DISTINGUISHABILITY"

panel_tmp_root="$(mktemp -d "${{SLURM_TMPDIR:-/tmp}}/tetra_arm_panel.XXXXXX")"
cleanup_panel_tmp() {{
    rm -rf -- "$panel_tmp_root"
}}
trap cleanup_panel_tmp EXIT
"$PANEL_UTILITY" \
    --vcf "$INTERINDIVIDUAL_PANEL" \
    --output "$panel_tmp_root/nuclear_panel_distinguishability.tsv"
if ! cmp -s \
        "$panel_tmp_root/nuclear_panel_distinguishability.tsv" \
        "$DISTINGUISHABILITY"; then
    echo "ERROR: selected donor-pair distinguishability table does not reproduce exactly from the declared interindividual panel" >&2
    exit 1
fi
{ploidy_runtime_block}

identity_validation_valid() {{
    python3 - "$IDENTITY_VALIDATION" "$LIBRARY_CSV" <<'PY'
import csv
import os
import re
import sys

path, library_csv = sys.argv[1:]
requested = {{int(value) for value in library_csv.split(",") if value}}
try:
    if not os.path.isfile(path) or os.path.getsize(path) <= 0:
        raise ValueError("missing or empty")
    with open(path, "r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
        fields = list(reader.fieldnames or [])
    normalized = {{field.strip().lower(): field for field in fields}}
    status_field = normalized.get("status")
    failure_field = next((normalized[name] for name in
        ("n_failures", "failure_count", "failures") if name in normalized), None)
    library_field = next((normalized[name] for name in
        ("library", "lib", "lib_num", "library_number") if name in normalized), None)
    if not rows or status_field is None or failure_field is None:
        raise ValueError("invalid validation summary schema")
    relevant = []
    for row in rows:
        raw = str(row.get(library_field, "")).strip() if library_field else ""
        match = re.fullmatch(r"(?:lib)?(\\d+)", raw, re.I)
        if match and int(match.group(1)) not in requested:
            continue
        relevant.append(row)
    if not relevant:
        raise ValueError("no relevant validation rows")
    for row in relevant:
        if (str(row.get(status_field, "")).strip().upper() != "PASS"
                or int(float(str(row.get(failure_field, "")))) != 0):
            raise ValueError("identity validation failure")
except (OSError, ValueError, TypeError, csv.Error):
    raise SystemExit(1)
PY
}}

if [[ "$SKIP_IDENTITY_VALIDATION" == "1" ]]; then
    echo "AUDIT OVERRIDE: --skip-identity-validation; upstream boundary was not validated" >&2
else
    test -s "$IDENTITY_VALIDATION"
    identity_validation_valid || {{
        echo "ERROR: identity validation boundary is not PASS: $IDENTITY_VALIDATION" >&2
        exit 1
    }}
fi

ledger_bundle_valid() {{
    python3 - "$LEDGER_ROOT" "$LIBRARY_CSV" <<'PY'
import csv
import gzip
import json
import os
import sys

root, library_csv = sys.argv[1:]
libraries = [int(value) for value in library_csv.split(",") if value]
try:
    counts = {{}}
    for library in libraries:
        path = os.path.join(root, f"lib{{library}}.final_cells.tsv.gz")
        if not os.path.isfile(path) or os.path.getsize(path) <= 0:
            raise ValueError(path)
        with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
            reader = csv.reader(handle, delimiter="\t")
            header = next(reader)
            if (len(header) != len(set(header)) or not
                    {{"library", "barcode", "production_assignment"}} <= set(header)):
                raise ValueError(path)
            library_col = header.index("library")
            count = 0
            for row in reader:
                if len(row) != len(header) or row[library_col] != str(library):
                    raise ValueError(path)
                count += 1
        if count < 1:
            raise ValueError(path)
        counts[library] = count
    summary = os.path.join(root, "split_ledger_summary.tsv")
    contract_path = os.path.join(root, "split_ledger_contract.json")
    if not os.path.isfile(summary) or os.path.getsize(summary) <= 0:
        raise ValueError(summary)
    with open(summary, "r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
        if list(reader.fieldnames or []) != [
                "library", "cells", "output", "schema_version"]:
            raise ValueError(summary)
    if len(rows) != len(libraries):
        raise ValueError(summary)
    seen = set()
    for row in rows:
        library = int(row["library"])
        expected = os.path.abspath(os.path.join(
            root, f"lib{{library}}.final_cells.tsv.gz"))
        if (library not in counts or library in seen
                or int(row["cells"]) != counts[library]
                or os.path.abspath(row["output"]) != expected
                or row["schema_version"] != "tetra_arm_split_ledger_v1"):
            raise ValueError(summary)
        seen.add(library)
    with open(contract_path, "r", encoding="utf-8") as handle:
        contract = json.load(handle)
    if (contract.get("schema_version") != "tetra_arm_split_ledger_contract_v1"
            or str(contract.get("status", "")).upper() != "PASS"
            or [int(value) for value in contract.get("libraries", [])] != libraries
            or int(contract.get("cells", -1)) != sum(counts.values())):
        raise ValueError(contract_path)
except (OSError, ValueError, TypeError, EOFError):
    raise SystemExit(1)
PY
}}

{output_checks}
test -s "$LEDGER_INPUT"
python3 "$PREPARE_SCRIPT" split-ledger \
    --input "$LEDGER_INPUT" \
    --libraries {library_words} \
    --output-root "$LEDGER_ROOT"

ledger_bundle_valid || {{
    echo "ERROR: LEDGER did not produce its complete validated bundle" >&2
    exit 1
}}
echo "COMPLETE: LEDGER $LEDGER_ROOT"
date
"""


def task_row_block(manifest: str, variables: Sequence[str]) -> str:
    names = " ".join(variables)
    expected_header = "\t".join(variables)
    return f"""TASK_MANIFEST={shlex.quote(manifest)}
test -s "$TASK_MANIFEST"
row="$(awk -F '\t' -v task="$SLURM_ARRAY_TASK_ID" \
    -v expected_header={shlex.quote(expected_header)} \
    -v expected_fields={len(variables)} '
    NR == 1 {{ if ($0 != expected_header) exit 2; next }}
    NF != expected_fields {{ exit 2 }}
    NR == task + 2 {{ print; found=1 }}
    END {{ if (!found) exit 3 }}
' "$TASK_MANIFEST")"
IFS=$'\t' read -r {names} <<< "$row"
test "$task_index" = "$SLURM_ARRAY_TASK_ID"
"""


def prepare_script(args: argparse.Namespace, run: RunPaths,
                   manifest: str, tasks: int) -> str:
    script = sbatch_header(
        "PREPARE", run, args, 2, args.prepare_memory,
        BASE_PYTHON_MODULES, ("csv", "gzip", "json"), tasks)
    script += "\ncommand -v python3 >/dev/null 2>&1\ncommand -v awk >/dev/null 2>&1\n"
    script += task_row_block(manifest, PREPARE_HEADER)
    return script + f"""
PREPARE_SCRIPT={shlex.quote(args.prepare_script)}
test -s "$PREPARE_SCRIPT"
python3 "$PREPARE_SCRIPT" --version

cell_manifest="$output_dir/lib${{library}}.cell_manifest.tsv.gz"
ambient_sources="$output_dir/lib${{library}}.ambient_sources.tsv.gz"
prepare_qc="$output_dir/lib${{library}}.prepare_qc.tsv"
prepare_contract="$output_dir/lib${{library}}.prepare_contract.json"

prepare_bundle_valid() {{
    python3 - "$cell_manifest" "$ambient_sources" "$prepare_qc" \
        "$prepare_contract" "$library" <<'PY'
import csv
import gzip
import json
import os
import sys

(cells, ambient, qc, contract_path, library) = sys.argv[1:]

try:
    for path in (cells, ambient, qc, contract_path):
        if not os.path.isfile(path) or os.path.getsize(path) <= 0:
            raise ValueError(path)
    with gzip.open(cells, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        cell_fields = list(reader.fieldnames or [])
        if len(cell_fields) != len(set(cell_fields)):
            raise ValueError(cells)
        if not {{"library", "barcode", "donor_a", "donor_b", "ambient_c",
                "model_eligible", "schema_version"}} <= set(cell_fields):
            raise ValueError(cells)
        cell_count = 0
        barcodes = set()
        for row in reader:
            barcode = row.get("barcode", "")
            if (None in row or any(value is None for value in row.values())
                    or row.get("library") != str(library) or not barcode
                    or barcode in barcodes
                    or row.get("schema_version")
                       != "tetra_arm_cell_manifest_v1"):
                raise ValueError(cells)
            barcodes.add(barcode)
            cell_count += 1
    with gzip.open(ambient, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        ambient_fields = list(reader.fieldnames or [])
        if len(ambient_fields) != len(set(ambient_fields)):
            raise ValueError(ambient)
        if not {{"library", "barcode", "source_label", "scoring_profile_mass",
                "schema_version"}} <= set(ambient_fields):
            raise ValueError(ambient)
        ambient_count = 0
        for row in reader:
            if (None in row or any(value is None for value in row.values())
                    or row.get("library") != str(library)
                    or row.get("schema_version")
                       != "tetra_arm_ambient_sources_v1"):
                raise ValueError(ambient)
            ambient_count += 1
    if cell_count < 1 or ambient_count < 1:
        raise ValueError("header-only PREPARE output")
    with open(qc, "r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        qc_fields = list(reader.fieldnames or [])
        rows = list(reader)
    expected_qc = ["library", "cells", "model_eligible",
                   "calibration_eligible", "ambient_rows", "profile_modes",
                   "eligibility_reasons", "status", "schema_version"]
    if (qc_fields != expected_qc or len(rows) != 1
            or rows[0].get("status", "").upper() != "PASS"
            or rows[0].get("schema_version") != "tetra_arm_prepare_qc_v1"
            or rows[0].get("library") != str(library)
            or int(rows[0].get("cells", -1)) != cell_count
            or int(rows[0].get("ambient_rows", -1)) != ambient_count):
        raise ValueError(qc)
    with open(contract_path, "r", encoding="utf-8") as handle:
        contract = json.load(handle)
    if (contract.get("schema_version") != "tetra_arm_prepare_contract_v1"
            or str(contract.get("status", "")).upper() != "PASS"
            or int(contract.get("library")) != int(library)
            or int(contract.get("cells", -1)) != cell_count):
        raise ValueError(contract_path)
except (OSError, ValueError, TypeError, EOFError):
    raise SystemExit(1)
PY
}}

for target in "$cell_manifest" "$ambient_sources" "$prepare_qc" "$prepare_contract"; do
    if [[ -e "$target" || -L "$target" ]]; then
        echo "ERROR: refusing to overwrite PREPARE output: $target; use a new --run-root" >&2
        exit 1
    fi
done
test -s "$ledger"
test -s "$final_assignments"
test -s "${{ambient_standard_prefix}}.contam_rate"
test -s "${{ambient_standard_prefix}}.contam_prof"

command=(
    python3 "$PREPARE_SCRIPT" prepare-library
    --library "$library"
    --ledger "$ledger"
    --final-assignments "$final_assignments"
    --ambient-standard-prefix "$ambient_standard_prefix"
    --ambient-arm-a-prefix "$ambient_arm_a_prefix"
    --ambient-arm-c-prefix "$ambient_arm_c_prefix"
    --min-calibration-tet-probability {args.min_calibration_tet_probability:.17g}
    --output-dir "$output_dir"
)
if [[ "$panel_metadata" != "NA" && -s "$panel_metadata" ]]; then
    command+=(--panel-metadata "$panel_metadata")
fi
if [[ -s "$ploidy_nn" ]]; then
    command+=(--ploidy-nn "$ploidy_nn")
fi
if [[ "$cell_groups" != "NA" ]]; then
    test -s "$cell_groups"
    command+=(--cell-groups "$cell_groups")
fi
"${{command[@]}}"

prepare_bundle_valid || {{
    echo "ERROR: PREPARE lib${{library}} did not produce its complete bundle" >&2
    exit 1
}}
echo "COMPLETE: PREPARE lib${{library}}"
date
"""


def ase_script(args: argparse.Namespace, run: RunPaths,
               manifest: str, tasks: int) -> str:
    script = sbatch_header(
        "ASE", run, args, args.ase_threads, args.ase_memory,
        ASE_MODULES, (), tasks)
    script += (
        "\ncommand -v awk >/dev/null 2>&1\n"
        "command -v gzip >/dev/null 2>&1\n"
        "command -v ln >/dev/null 2>&1\n"
        "command -v rm >/dev/null 2>&1\n")
    script += task_row_block(manifest, ASE_HEADER)
    sex_option = ("\ncommand+=(--include-sex-chromosomes)"
                  if args.include_sex_chromosomes else "")
    return script + f"""
ASE_BINARY={shlex.quote(args.ase_binary)}
EXPECTED_ASE_HEADER={shlex.quote(chr(9).join(ASE_V2_HEADER))}
test -x "$ASE_BINARY"
command -v "$ASE_BINARY" >/dev/null 2>&1
"$ASE_BINARY" --version

ase_paths_valid() {{
    local data_path="$1"
    local qc_path="$2"
    local status
    local threshold
    local sex_state
    local evidence_basis
    local qc_schema
    local allele_contract
    local likelihood_contract
    local output_rows
    local passing_rows
    [[ -s "$data_path" && -s "$qc_path" ]] || return 1
    awk -F '\t' '
        NR == 1 {{ if (NF != 2 || $1 != "metric" || $2 != "value") exit 1; next }}
        {{ if (NF != 2 || seen[$1]++) exit 1 }}
        END {{ if (NR < 2) exit 1 }}
    ' "$qc_path" || return 1
    status="$(awk -F '\t' '
        $1 == "status" {{ count++; value=$2 }}
        END {{
            if (count != 1) exit 1
            print value
        }}
    ' "$qc_path")" || return 1
    [[ "$status" == "PASS" || "$status" == "PASS_NO_HETEROTYPIC_TARGETS" \
        || "$status" == "PASS_NO_INFORMATIVE_EVIDENCE" ]] \
        || return 1
    threshold="$(awk -F '\t' '
        $1 == "hard_allele_threshold" {{ count++; value=$2 }}
        END {{
            if (count != 1) exit 1
            print value
        }}
    ' "$qc_path")" || return 1
    awk -v observed="$threshold" \
        -v expected={args.hard_allele_threshold:.17g} \
        'BEGIN {{ delta=observed-expected; if (delta<0) delta=-delta; exit !(delta<=1e-12) }}' \
        || return 1
    sex_state="$(awk -F '\t' '
        $1 == "sex_chromosomes" {{ count++; value=$2 }}
        END {{
            if (count != 1) exit 1
            print value
        }}
    ' "$qc_path")" || return 1
    [[ "$sex_state" == {shlex.quote("INCLUDED_BY_OPT_IN" if args.include_sex_chromosomes else "EXCLUDED_DEFAULT")} ]] \
        || return 1
    evidence_basis="$(awk -F '\t' '
        $1 == "evidence_basis" {{ count++; value=$2 }}
        END {{
            if (count != 1) exit 1
            print value
        }}
    ' "$qc_path")" || return 1
    [[ "$evidence_basis" == "MOLECULE" \
        || "$evidence_basis" == "SITE_FALLBACK" ]] || return 1
    qc_schema="$(awk -F '\t' '
        $1 == "schema_version" {{ count++; value=$2 }}
        END {{ if (count != 1) exit 1; print value }}
    ' "$qc_path")" || return 1
    [[ "$qc_schema" == "tetra_arm_ase_qc_v2" ]] || return 1
    allele_contract="$(awk -F '\t' '
        $1 == "allele_evidence_contract" {{ count++; value=$2 }}
        END {{ if (count != 1) exit 1; print value }}
    ' "$qc_path")" || return 1
    [[ "$allele_contract" == "ORIENTATION_CONFIDENCE_WEIGHTED_SOFT_PRIMARY_HARD_CALLS_QC_ONLY" ]] \
        || return 1
    likelihood_contract="$(awk -F '\t' '
        $1 == "soft_likelihood_contract" {{ count++; value=$2 }}
        END {{ if (count != 1) exit 1; print value }}
    ' "$qc_path")" || return 1
    [[ "$likelihood_contract" == "EFFECTIVE_FRACTIONAL_QUASI_LIKELIHOOD_NOT_EXACT_BINOMIAL_PMF" ]] \
        || return 1
    output_rows="$(awk -F '\t' '
        $1 == "output_rows" {{ count++; value=$2 }}
        END {{ if (count != 1 || value !~ /^[0-9]+$/) exit 1; print value }}
    ' "$qc_path")" || return 1
    passing_rows="$(awk -F '\t' '
        $1 == "passing_rows" {{ count++; value=$2 }}
        END {{ if (count != 1 || value !~ /^[0-9]+$/) exit 1; print value }}
    ' "$qc_path")" || return 1
    gzip -t "$data_path" || return 1
    gzip -cd "$data_path" | awk -F '\t' -v expected_status="$status" \
        -v expected_header="$EXPECTED_ASE_HEADER" \
        -v expected_library="$library" \
        -v expected_rows="$output_rows" -v passing_rows="$passing_rows" '
        NR == 1 {{
            if ($0 != expected_header) exit 1
            header_fields = NF
            next
        }}
        {{
            if (NF != header_fields || $1 != expected_library ||
                $NF != "tetra_arm_ase_evidence_v2") exit 1
            if ($(NF - 1) == "PASS" || $(NF - 1) == "PASS_SITE_FALLBACK")
                observed_passing++
            else if ($(NF - 1) != "NO_DIRECTIONAL_SOFT_EVIDENCE") exit 1
        }}
        END {{
            if (NR - 1 != expected_rows) exit 1
            if (passing_rows > expected_rows || observed_passing != passing_rows)
                exit 1
            if (expected_status == "PASS" && passing_rows < 1) exit 1
            if (expected_status == "PASS_NO_HETEROTYPIC_TARGETS" &&
                (NR != 1 || passing_rows != 0)) exit 1
            if (expected_status == "PASS_NO_INFORMATIVE_EVIDENCE" &&
                passing_rows != 0) exit 1
        }}
    '
}}

for target in "$output" "$qc"; do
    if [[ -e "$target" || -L "$target" ]]; then
        echo "ERROR: refusing to overwrite ASE output: $target; use a new --run-root" >&2
        exit 1
    fi
done
test -s "$samples"
test -s "$pileup_sites"
test -s "$cell_manifest"
test -s "$ambient_sources"
test -s "$arms_bed"
if [[ ! -s "$pileup_molecules" && ! -s "$pileup_observations" ]]; then
    echo "ERROR: lib${{library}} has neither molecule nor observation pileup" >&2
    exit 1
fi
for evidence_input in "$pileup_molecules" "$pileup_observations"; do
    if [[ -s "$evidence_input" ]]; then
        gzip -t "$evidence_input"
    fi
done

tmp_output="${{output%.gz}}.work.${{SLURM_JOB_ID}}.${{SLURM_ARRAY_TASK_ID}}.gz"
tmp_qc="${{qc}}.work.${{SLURM_JOB_ID}}.${{SLURM_ARRAY_TASK_ID}}"
cleanup() {{
    rm -f -- "$tmp_output" "$tmp_qc"
}}
trap cleanup EXIT
cleanup

command=(
    "$ASE_BINARY"
    --samples "$samples"
    --panel-bcf {shlex.quote(args.interindividual_panel)}
    --pileup-sites "$pileup_sites"
    --cells "$cell_manifest"
    --ambient-sources "$ambient_sources"
    --arms "$arms_bed"
    --library "$library"
    --output "$tmp_output"
    --qc "$tmp_qc"
    --threads "$SLURM_CPUS_PER_TASK"
    --hard-allele-threshold {args.hard_allele_threshold:.17g}
){sex_option}
if [[ -s "$pileup_molecules" ]]; then
    command+=(--pileup-molecules "$pileup_molecules")
fi
if [[ -s "$pileup_observations" ]]; then
    command+=(--pileup-observations "$pileup_observations")
fi
"${{command[@]}}"

ase_paths_valid "$tmp_output" "$tmp_qc"
ln -- "$tmp_qc" "$qc" || {{
    echo "ERROR: could not publish ASE QC without overwriting $qc" >&2
    exit 1
}}
if ! ln -- "$tmp_output" "$output"; then
    rm -f -- "$qc"
    echo "ERROR: could not publish ASE data without overwriting $output" >&2
    exit 1
fi
ase_paths_valid "$output" "$qc" || {{
    rm -f -- "$output" "$qc"
    echo "ERROR: ASE lib${{library}} did not produce its complete bundle" >&2
    exit 1
}}
echo "COMPLETE: ASE lib${{library}}"
date
"""


def expression_script(args: argparse.Namespace, run: RunPaths,
                      manifest: str, tasks: int) -> str:
    script = sbatch_header(
        "EXPRESSION", run, args, args.expression_cpus,
        args.expression_memory, SCIENTIFIC_PYTHON_MODULES,
        ("numpy", "scipy", "scipy.io", "scipy.sparse"), tasks)
    script += "\ncommand -v python3 >/dev/null 2>&1\ncommand -v awk >/dev/null 2>&1\n"
    script += task_row_block(manifest, EXPRESSION_HEADER)
    optional = ""
    if args.include_sex_chromosomes:
        optional += "\ncommand+=(--include-sex-chromosomes)"
    if args.allow_missing_expression_barcodes:
        optional += "\ncommand+=(--allow-missing-barcodes)"
    if args.include_nonheterotypic_expression:
        optional += "\ncommand+=(--include-nonheterotypic)"
    if args.expression_ambient_corrected:
        optional += "\ncommand+=(--ambient-corrected)"
    return script + f"""
EXPRESSION_SCRIPT={shlex.quote(args.expression_script)}
test -s "$EXPRESSION_SCRIPT"
python3 "$EXPRESSION_SCRIPT" --version

expression="$output_dir/lib${{library}}.arm_expression.tsv.gz"
expression_qc="$output_dir/lib${{library}}.expression_qc.tsv"
expression_contract="$output_dir/lib${{library}}.expression_contract.json"

expression_bundle_valid() {{
    python3 - "$expression" "$expression_qc" "$expression_contract" "$library" \
        {shlex.quote("UPSTREAM_AMBIENT_CORRECTED" if args.expression_ambient_corrected else "OBSERVED_FILTERED_COUNTS")} <<'PY'
import csv
import gzip
import json
import os
import sys

(evidence, qc, contract_path, library, expected_state) = sys.argv[1:]
expected_evidence = [
    "library", "barcode", "arm", "chromosome", "arm_counts",
    "other_autosomal_counts", "reference_autosomal_counts",
    "total_autosomal_counts", "arm_fraction", "log2_arm_to_other",
    "log2_arm_to_reference", "mapped_genes_on_arm", "nonzero_genes_on_arm",
    "matrix_value_type", "expression_input_state", "schema_version",
]
expected_qc = [
    "library", "manifest_cells", "selected_manifest_cells", "matrix_cells",
    "overlap_cells", "missing_manifest_cells", "matrix_genes",
    "mapped_genes", "arms", "rows", "status", "schema_version",
]

try:
    for path in (evidence, qc, contract_path):
        if not os.path.isfile(path) or os.path.getsize(path) <= 0:
            raise ValueError(path)
    with gzip.open(evidence, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = list(reader.fieldnames or [])
        if fields != expected_evidence:
            raise ValueError(evidence)
        row_count = 0
        for row in reader:
            if (None in row or any(value is None for value in row.values())
                    or row.get("library") != str(library)
                    or row.get("schema_version")
                       != "tetra_arm_expression_evidence_v1"
                    or row.get("expression_input_state") != expected_state):
                raise ValueError(evidence)
            row_count += 1
    with open(qc, "r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        qc_fields = list(reader.fieldnames or [])
        rows = list(reader)
    if (qc_fields != expected_qc or len(rows) != 1
            or rows[0].get("library") != str(library)
            or rows[0].get("schema_version") != "tetra_arm_expression_qc_v1"
            or int(rows[0].get("rows", -1)) != row_count):
        raise ValueError(qc)
    qc_status = rows[0].get("status", "").upper()
    with open(contract_path, "r", encoding="utf-8") as handle:
        contract = json.load(handle)
    terminal_state = str(contract.get("terminal_state", "")).upper()
    row_state_ok = (
        row_count > 0 and qc_status == "PASS" and not terminal_state
        or (row_count == 0 and qc_status == "PASS_NO_HETEROTYPIC_TARGETS"
            and terminal_state == "PASS_NO_HETEROTYPIC_TARGETS"))
    if (contract.get("schema_version") != "tetra_arm_expression_contract_v1"
            or str(contract.get("status", "")).upper() != "PASS"
            or int(contract.get("library")) != int(library)
            or int(contract.get("rows", -1)) != row_count
            or contract.get("evidence_fields") != expected_evidence
            or not row_state_ok
            or (not terminal_state
                and contract.get("expression_input_state") != expected_state)):
        raise ValueError(contract_path)
except (OSError, ValueError, TypeError, EOFError):
    raise SystemExit(1)
PY
}}

for target in "$expression" "$expression_qc" "$expression_contract"; do
    if [[ -e "$target" || -L "$target" ]]; then
        echo "ERROR: refusing to overwrite EXPRESSION output: $target; use a new --run-root" >&2
        exit 1
    fi
done
test -s "$barcodes"
test -s "$features"
test -s "$matrix"
test -s "$cell_manifest"
test -s "$gene_arms"

command=(
    python3 "$EXPRESSION_SCRIPT"
    --library "$library"
    --barcodes "$barcodes"
    --features "$features"
    --matrix "$matrix"
    --cell-manifest "$cell_manifest"
    --gene-arms "$gene_arms"
    --output-dir "$output_dir"
    --pseudocount {args.expression_pseudocount:.17g}
){optional}
"${{command[@]}}"

expression_bundle_valid || {{
    echo "ERROR: EXPRESSION lib${{library}} did not produce its complete bundle" >&2
    exit 1
}}
echo "COMPLETE: EXPRESSION lib${{library}}"
date
"""


def call_script(args: argparse.Namespace, run: RunPaths,
                paths: Sequence[LibraryPaths]) -> str:
    script = sbatch_header(
        "CALL", run, args, args.call_cpus, args.call_memory,
        BASE_PYTHON_MODULES, ("csv", "gzip", "json"))
    ase_inputs = [path.ase for path in paths]
    manifests = [path.cell_manifest for path in paths]
    expression_inputs = [path.expression for path in paths]
    input_checks = "\n".join(
        f"test -s {shlex.quote(value)}"
        for value in (*ase_inputs, *manifests, *expression_inputs))
    ase_words = "\n".join(f"    {shlex.quote(value)}" for value in ase_inputs)
    manifest_words = "\n".join(
        f"    {shlex.quote(value)}" for value in manifests)
    expression_words = "\n".join(
        f"    {shlex.quote(value)}" for value in expression_inputs)
    parameter_words = "\n".join(
        f"    --{name.replace('_', '-')} " +
        (str(value) if isinstance(value, int) else f"{value:.17g}")
        for name, value in call_parameters(args).items())
    return script + f"""
command -v python3 >/dev/null 2>&1
CALL_SCRIPT={shlex.quote(args.call_script)}
OUTPUT_PREFIX={shlex.quote(run.call_prefix)}
test -s "$CALL_SCRIPT"
python3 "$CALL_SCRIPT" --version

call_bundle_valid() {{
    python3 - "$OUTPUT_PREFIX" <<'PY'
import csv
import gzip
import json
import os
import sys

prefix = sys.argv[1]
suffixes = (
    ".arm_calls.tsv.gz", ".uid_chromosome_flags.tsv.gz",
    ".calibration.tsv.gz", ".qc.tsv", ".contract.json",
    ".donor_pair_arm_summary.tsv.gz")
try:
    for suffix in suffixes:
        path = prefix + suffix
        if not os.path.isfile(path) or os.path.getsize(path) <= 0:
            raise ValueError(path)
    gzip_headers = (
        (
            prefix + ".arm_calls.tsv.gz",
            {{"library", "barcode", "donor_pair", "arm", "chromosome",
              "best_state", "best_state_posterior",
              "empirical_q_resolution_floor", "call_state", "call_status",
              "call_schema_version"}},
            179, "call_schema_version", "tetra_arm_cnv_calls_v2",
        ),
        (
            prefix + ".calibration.tsv.gz",
            {{"library", "barcode", "calibration_group", "donor_pair", "arm",
              "crossfit_fold", "excluded_chromosome",
              "cell_loo_baseline_status", "ref_source_level",
              "ref_quasi_overdispersion_rho", "expression_source_level",
              "calibration_status", "schema_version"}},
            52, "schema_version", "tetra_arm_calibration_v2",
        ),
        (
            prefix + ".uid_chromosome_flags.tsv.gz",
            {{"uid", "donor_pair", "chromosome", "p_state", "q_state",
              "whole_chromosome_flag", "summary_status", "schema_version"}},
            57, "schema_version", "tetra_arm_uid_chromosome_flags_v2",
        ),
        (
            prefix + ".donor_pair_arm_summary.tsv.gz",
            {{"library", "calibration_group", "donor_pair", "arm", "chromosome",
              "individual_evaluable_cells", "aggregate_eligible_cells",
              "tested_cells", "tested_uid_blocks", "supporting_cells",
              "concordant_cells", "best_state", "partial_conjunction_method",
              "dependence_assumption", "fdr_interpretation",
              "partial_conjunction_q_value",
              "partial_conjunction_q_resolution_floor", "recurrence_flag",
              "summary_status", "schema_version"}},
            45, "schema_version", "tetra_arm_donor_pair_arm_summary_v2",
        ),
    )
    row_counts = {{}}
    for path, required, expected_columns, schema_field, expected_schema in gzip_headers:
        with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            fields = list(reader.fieldnames or [])
            if (len(fields) != expected_columns or len(fields) != len(set(fields))
                    or not required <= set(fields)):
                raise ValueError(path)
            count = 0
            for row in reader:
                if (None in row or any(value is None for value in row.values())
                        or row.get(schema_field) != expected_schema):
                    raise ValueError(path)
                count += 1
            row_counts[path] = count
    with open(prefix + ".qc.tsv", "r", encoding="utf-8", newline="") as handle:
        rows = list(csv.reader(handle, delimiter="\t"))
    if not rows or rows[0] != ["metric", "value"]:
        raise ValueError("call QC")
    metrics = {{}}
    for row in rows[1:]:
        if len(row) != 2 or row[0] in metrics:
            raise ValueError("call QC")
        metrics[row[0]] = row[1]
    with open(prefix + ".contract.json", "r", encoding="utf-8") as handle:
        contract = json.load(handle)
    output_schemas = contract.get("output_schemas", {{}})
    output_suffixes = contract.get("output_suffixes", {{}})
    terminal_state = str(contract.get("terminal_state", "NONE")).upper()
    status = str(contract.get("status", "")).upper()
    qc_status = str(metrics.get("status", "")).upper()
    qc_terminal = str(metrics.get("terminal_state", "NONE")).upper()
    allowed_terminal = {{
        "PASS_NO_HETEROTYPIC_TARGETS",
        "PASS_NO_OBSERVED_ASE",
        "PASS_NO_CALLABLE_ASE",
    }}
    calls_rows = row_counts[prefix + ".arm_calls.tsv.gz"]
    calibration_rows = row_counts[prefix + ".calibration.tsv.gz"]
    uid_rows = row_counts[prefix + ".uid_chromosome_flags.tsv.gz"]
    pair_rows = row_counts[prefix + ".donor_pair_arm_summary.tsv.gz"]
    state_ok = (
        terminal_state == "NONE" and qc_terminal == "NONE"
        and status == "PASS" and qc_status == "PASS" and calls_rows > 0
        or terminal_state in allowed_terminal
        and qc_terminal == terminal_state and status == terminal_state
        and qc_status == terminal_state
        and (terminal_state != "PASS_NO_HETEROTYPIC_TARGETS"
             or calls_rows == 0)
        and (terminal_state != "PASS_NO_CALLABLE_ASE" or calls_rows > 0)
    )
    if (contract.get("schema_version") != "tetra_arm_call_contract_v2"
            or metrics.get("schema_version") != "tetra_arm_call_qc_v2"
            or int(metrics.get("expanded_output_rows", -1)) != calls_rows
            or int(metrics.get("calibration_rows", -1)) != calibration_rows
            or int(metrics.get("uid_chromosome_rows", -1)) != uid_rows
            or int(metrics.get("donor_pair_arm_summary_rows", -1)) != pair_rows
            or not state_ok
            or not isinstance(output_schemas, dict)
            or output_schemas.get("calls") != "tetra_arm_cnv_calls_v2"
            or output_schemas.get("calibration") != "tetra_arm_calibration_v2"
            or output_schemas.get("uid_chromosome_flags")
               != "tetra_arm_uid_chromosome_flags_v2"
            or output_schemas.get("donor_pair_arm_summary")
               != "tetra_arm_donor_pair_arm_summary_v2"
            or output_schemas.get("qc") != "tetra_arm_call_qc_v2"
            or output_schemas.get("contract")
               != "tetra_arm_call_contract_v2"
            or not isinstance(output_suffixes, dict)
            or output_suffixes.get("calls") != ".arm_calls.tsv.gz"
            or output_suffixes.get("calibration") != ".calibration.tsv.gz"
            or output_suffixes.get("uid_chromosome_flags")
               != ".uid_chromosome_flags.tsv.gz"
            or output_suffixes.get("donor_pair_arm_summary")
               != ".donor_pair_arm_summary.tsv.gz"
            or output_suffixes.get("qc") != ".qc.tsv"
            or output_suffixes.get("contract") != ".contract.json"):
        raise ValueError("call contract")
except (OSError, ValueError, TypeError, EOFError):
    raise SystemExit(1)
PY
}}

for suffix in \
    .arm_calls.tsv.gz .uid_chromosome_flags.tsv.gz .calibration.tsv.gz \
    .qc.tsv .contract.json .donor_pair_arm_summary.tsv.gz; do
    target="${{OUTPUT_PREFIX}}${{suffix}}"
    if [[ -e "$target" || -L "$target" ]]; then
        echo "ERROR: refusing to overwrite CALL output: $target; use a new --run-root" >&2
        exit 1
    fi
done
{input_checks}
command=(
    python3 "$CALL_SCRIPT"
    --ase
{ase_words}
    --cell-manifest
{manifest_words}
    --expression
{expression_words}
    --output-prefix "$OUTPUT_PREFIX"
{parameter_words}
)
"${{command[@]}}"

call_bundle_valid || {{
    echo "ERROR: CALL did not produce its complete validated bundle" >&2
    exit 1
}}
echo "COMPLETE: CALL $OUTPUT_PREFIX"
date
"""


def report_script(args: argparse.Namespace, run: RunPaths,
                  qc_files: Sequence[str]) -> str:
    script = sbatch_header(
        "REPORT", run, args, args.report_cpus, args.report_memory,
        BASE_PYTHON_MODULES, ("csv", "gzip", "json"))
    call_qc = run.call_prefix + ".qc.tsv"
    supplemental_qc = [value for value in qc_files if value != call_qc]
    optional_qc_words = "\n".join(
        f"    {shlex.quote(value)}" for value in supplemental_qc)
    optional_qc_array = (
        f"OPTIONAL_REPORT_QC_FILES=(\n{optional_qc_words}\n)"
        if supplemental_qc else "OPTIONAL_REPORT_QC_FILES=()")
    return script + f"""
command -v python3 >/dev/null 2>&1
REPORT_SCRIPT={shlex.quote(args.report_script)}
CALLS={shlex.quote(run.call_prefix + '.arm_calls.tsv.gz')}
CALIBRATION={shlex.quote(run.call_prefix + '.calibration.tsv.gz')}
UID_SUMMARY={shlex.quote(run.call_prefix + '.uid_chromosome_flags.tsv.gz')}
DONOR_PAIR_SUMMARY={shlex.quote(run.donor_pair_summary)}
CALL_CONTRACT={shlex.quote(run.call_prefix + '.contract.json')}
CALL_QC={shlex.quote(call_qc)}
REPORT_DIR={shlex.quote(run.report)}
SUMMARY="$REPORT_DIR/tetra_arm_cnv_summary.tsv"
SUMMARY_JSON="$REPORT_DIR/tetra_arm_cnv_summary.json"
REPORT_HTML="$REPORT_DIR/tetra_arm_cnv_report.html"
{optional_qc_array}
REPORT_QC_FILES=("$CALL_QC")
for candidate in "${{OPTIONAL_REPORT_QC_FILES[@]}}"; do
    if [[ -s "$candidate" ]]; then
        REPORT_QC_FILES+=("$candidate")
    fi
done

test -s "$REPORT_SCRIPT"
python3 "$REPORT_SCRIPT" --version

report_bundle_valid() {{
    python3 - "$SUMMARY" "$SUMMARY_JSON" "$REPORT_HTML" <<'PY'
import csv
import json
import os
import sys
summary, json_path, html_path = sys.argv[1:]

try:
    for path in (summary, json_path, html_path):
        if not os.path.isfile(path) or os.path.getsize(path) <= 0:
            raise ValueError(path)
    with open(summary, "r", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        if next(reader) != ["metric", "value", "definition"]:
            raise ValueError(summary)
        metrics = {{}}
        for row in reader:
            if len(row) != 3 or not row[0] or row[0] in metrics:
                raise ValueError(summary)
            metrics[row[0]] = row[1]
    with open(json_path, "r", encoding="utf-8") as handle:
        payload = json.load(handle)
    if (payload.get("schema_version") != "tetra_arm_cnv_report_v2"
            or str(payload.get("status", "")).upper()
               not in {{"PASS", "REVIEW", "PASS_NO_HETEROTYPIC_TARGETS",
                       "PASS_NO_OBSERVED_ASE", "PASS_NO_CALLABLE_ASE"}}
            or metrics.get("schema_version") != "tetra_arm_cnv_report_v2"
            or metrics.get("status", "").upper()
               != str(payload.get("status", "")).upper()):
        raise ValueError(json_path)
except (OSError, ValueError, TypeError, EOFError):
    raise SystemExit(1)
PY
}}

for target in "$SUMMARY" "$SUMMARY_JSON" "$REPORT_HTML"; do
    if [[ -e "$target" || -L "$target" ]]; then
        echo "ERROR: refusing to overwrite REPORT output: $target; use a new --run-root" >&2
        exit 1
    fi
done
test -s "$CALLS"
test -s "$CALIBRATION"
test -s "$UID_SUMMARY"
test -s "$DONOR_PAIR_SUMMARY"
test -s "$CALL_CONTRACT"
test -s "$CALL_QC"

command=(
    python3 "$REPORT_SCRIPT"
    --calls "$CALLS"
    --calibration "$CALIBRATION"
    --uid-summary "$UID_SUMMARY"
    --donor-pair-summary "$DONOR_PAIR_SUMMARY"
    --call-contract "$CALL_CONTRACT"
    --output-dir "$REPORT_DIR"
    --title {shlex.quote(args.report_title)}
    --top-calls {args.report_top_calls}
    --qc-files "${{REPORT_QC_FILES[@]}}"
)
"${{command[@]}}"

report_bundle_valid || {{
    echo "ERROR: REPORT did not produce its complete validated bundle" >&2
    exit 1
}}
echo "COMPLETE: REPORT $REPORT_HTML"
date
"""


def render_script(path: str, payload: str) -> str:
    if not payload.startswith("#!/bin/bash\n"):
        raise RuntimeError("generated sbatch lacks the required bash shebang")
    required = (
        "set -euo pipefail", "module purge", "module list 2>&1", "command -v",
    )
    missing = [token.strip() for token in required if token not in payload]
    if missing:
        raise RuntimeError(
            f"generated sbatch lacks required invariant(s): {', '.join(missing)}")
    module_loads = re.findall(r"(?m)^module load ([^\s]+)$", payload)
    if not module_loads or any("/" not in module for module in module_loads):
        raise RuntimeError(
            "generated sbatch must load at least one concrete versioned module")
    result = publish_new_or_identical(path, payload, "generated sbatch")
    os.chmod(result, os.stat(result).st_mode | stat.S_IXUSR)
    checked = subprocess.run(
        ["bash", "-n", result], capture_output=True, text=True, check=False)
    if checked.returncode != 0:
        detail = checked.stderr.strip() or checked.stdout.strip()
        raise RuntimeError(f"generated sbatch failed bash -n: {result}: {detail}")
    return result


def submit_job(script: str, dependencies: Sequence[str], run_root: str) -> str:
    command = ["sbatch", "--parsable", "--chdir", absolute(run_root)]
    if dependencies:
        command.append("--dependency=afterok:" + ":".join(dependencies))
    command.append(script)
    failures = []
    transient = re.compile(
        r"temporar|timed?\s*out|unable\s+to\s+contact.*controller|"
        r"connection\s+(?:refused|reset)", re.I)
    for attempt in range(2):
        result = subprocess.run(
            command, capture_output=True, text=True, check=False)
        token = result.stdout.strip().split(";", 1)[0]
        if result.returncode == 0 and re.fullmatch(r"\d+", token):
            return token
        detail = result.stderr.strip() or result.stdout.strip() or "no scheduler response"
        failures.append(f"attempt {attempt + 1}: {detail}")
        if result.returncode == 0 or attempt == 1 or not transient.search(detail):
            break
    raise RuntimeError(
        f"sbatch --parsable failed for {script}: {'; '.join(failures)}")


def resolved_dependencies(stage: str, selected: Sequence[str],
                          job_ids: Mapping[str, str]) -> list[str]:
    """Return the nearest submitted ancestors on every upstream DAG branch."""
    selected_set = set(selected)
    frontier: set[str] = set()

    def visit(node: str) -> None:
        for parent in STAGE_PARENTS[node]:
            if parent in selected_set and parent in job_ids:
                frontier.add(parent)
            else:
                visit(parent)

    visit(stage)
    return [job_ids[parent] for parent in STAGES if parent in frontier]


def prepare_manifests(
        args: argparse.Namespace, run: RunPaths,
        paths: Sequence[LibraryPaths],
        selected: Sequence[str]) -> dict[str, str]:
    manifests: dict[str, str] = {}
    library_key = "-".join(str(path.library) for path in paths)

    prepare_rows = [{
        "task_index": index,
        "library": path.library,
        "ledger": path.split_ledger,
        "final_assignments": path.final_assignments,
        "ambient_standard_prefix": path.ambient_standard_prefix,
        "ambient_arm_a_prefix": path.ambient_arm_a_prefix,
        "ambient_arm_c_prefix": path.ambient_arm_c_prefix,
        "panel_metadata": args.panel_metadata or "NA",
        "ploidy_nn": path.ploidy_nn,
        "cell_groups": path.cell_groups or "NA",
        "output_dir": run.prepare,
    } for index, path in enumerate(paths)]
    if "PREPARE" in selected:
        manifests["PREPARE"] = write_task_manifest(
            os.path.join(
                run.manifests, f"prepare_tasks.libs_{library_key}.tsv"),
            PREPARE_HEADER, prepare_rows)

    ase_rows = [{
        "task_index": index,
        "library": path.library,
        "samples": path.samples,
        "pileup_sites": path.pileup_sites,
        "pileup_molecules": path.pileup_molecules,
        "pileup_observations": path.pileup_observations,
        "cell_manifest": path.cell_manifest,
        "ambient_sources": path.ambient_sources,
        "arms_bed": args.arms_bed,
        "output": path.ase,
        "qc": path.ase_qc,
    } for index, path in enumerate(paths)]
    if "ASE" in selected:
        manifests["ASE"] = write_task_manifest(
            os.path.join(run.manifests, f"ase_tasks.libs_{library_key}.tsv"),
            ASE_HEADER, ase_rows)

    expression_rows = [{
        "task_index": index,
        "library": path.library,
        "barcodes": path.expression_barcodes,
        "features": path.expression_features,
        "matrix": path.expression_matrix,
        "cell_manifest": path.cell_manifest,
        "gene_arms": args.gene_arms,
        "output_dir": run.expression,
    } for index, path in enumerate(paths)]
    if "EXPRESSION" in selected:
        manifests["EXPRESSION"] = write_task_manifest(
            os.path.join(
                run.manifests, f"expression_tasks.libs_{library_key}.tsv"),
            EXPRESSION_HEADER, expression_rows)
    return manifests


def report_qc_files(run: RunPaths,
                    paths: Sequence[LibraryPaths]) -> list[str]:
    return [
        run.reference_qc,
        run.call_prefix + ".qc.tsv",
        *(path.prepare_qc for path in paths),
        *(path.ase_qc for path in paths),
        *(path.expression_qc for path in paths),
    ]


def stage_primary_outputs(stage: str, run: RunPaths,
                          paths: Sequence[LibraryPaths],
                          args: argparse.Namespace) -> list[str]:
    if stage == "REFERENCE":
        outputs = [run.reference_qc, run.reference_contract]
        if args.reference_mode != "EXPLICIT_BED":
            outputs.insert(0, args.arms_bed)
        return outputs
    if stage == "LEDGER":
        return [
            *(path.split_ledger for path in paths),
            os.path.join(run.ledger, "split_ledger_summary.tsv"),
            os.path.join(run.ledger, "split_ledger_contract.json"),
        ]
    if stage == "PREPARE":
        return [value for path in paths for value in (
            path.cell_manifest, path.ambient_sources,
            path.prepare_qc, path.prepare_contract)]
    if stage == "ASE":
        return [value for path in paths for value in (path.ase, path.ase_qc)]
    if stage == "EXPRESSION":
        return [value for path in paths for value in (
            path.expression, path.expression_qc, path.expression_contract)]
    if stage == "CALL":
        return list(run.call_outputs)
    if stage == "REPORT":
        return list(run.report_outputs)
    raise ValueError(f"unknown stage: {stage}")


def reject_output_conflicts(selected: Sequence[str], run: RunPaths,
                            paths: Sequence[LibraryPaths],
                            args: argparse.Namespace) -> None:
    for stage in selected:
        conflicts = [path for path in stage_primary_outputs(stage, run, paths, args)
                     if os.path.lexists(path)]
        if conflicts:
            preview = ", ".join(conflicts[:3])
            suffix = " ..." if len(conflicts) > 3 else ""
            raise ValueError(
                f"selected stage {stage} already has {len(conflicts)} output "
                f"path(s): {preview}{suffix}; use a new --run-root rather than "
                "overwriting primary scientific outputs")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Generate and optionally submit the downstream chromosome-arm "
            "ASE/CNV SLURM DAG."))
    parser.add_argument("--version", action="version", version=f"%(prog)s {RELEASE}")
    parser.add_argument("--submit", action="store_true",
                        help="Submit generated jobs; otherwise print the plan")
    parser.add_argument(
        "--stage", "--stages", dest="stages", action="append",
        help="Selected stage(s), comma-separated or repeated; default ALL")
    parser.add_argument("--libraries", nargs="+", default=["1-40"],
                        help="Libraries as N, libN, N-M, or comma-separated values")
    parser.add_argument(
        "--upstream-cohort-libraries", nargs="+", default=None,
        help=("Optional larger upstream cohort to validate while producing "
              "outputs only for --libraries. Use 1-40 for a lib19 pilot whose "
              "reconciliation must include all cross-library evidence."))
    parser.add_argument("--array-throttle", type=int, default=None,
                        help="Optional maximum concurrent tasks per array")

    parser.add_argument(
        "--run-root", default=None,
        help=("Output root for this downstream DAG; defaults to "
              "<upstream-analysis-root>/aggregate_library_analysis/"
              "tetra_arm_cnv"))
    parser.add_argument(
        "--mapping-input-root", "--mapping-root",
        dest="mapping_input_root", default=MAPPING_ROOT,
        help=("Mapping-output root containing per-library filtered MEX files; "
              "--mapping-root is retained as a legacy alias"))
    parser.add_argument(
        "--upstream-analysis-root", default=None,
        help=("Root produced by orchestrate_tetraploid.py containing per-library "
              "DEMUX/ambient products and aggregate reconciliation, ploidy, and "
              "GEX-calibration products. Required with a non-production "
              "--mapping-input-root; production defaults to the historical "
              "combined mapping tree."))
    parser.add_argument(
        "--mapping-run-root", default="",
        help=("Top-level mapping run containing control/, validation/, and "
              "rna3/mapping_output. Auto-derived for the canonical layout and "
              "required by the non-production provenance check."))
    parser.add_argument(
        "--identity-root", default=None,
        help=("Advanced override; defaults to <upstream-analysis-root>/"
              "aggregate_library_analysis/identity_reconciliation"))
    parser.add_argument(
        "--identity-validation", default=None,
        help=("Identity reconciliation validation_summary.tsv; defaults below "
              "--identity-root and must be all-PASS before LEDGER submission"))
    parser.add_argument(
        "--identity-validation-failures", default=None,
        help=("Identity reconciliation validation_failures.tsv; defaults "
              "beside --identity-validation and must contain only its header"))
    parser.add_argument(
        "--identity-run-summary", default=None,
        help=("Phase-3 identity_reconciliation_run_summary.tsv; defaults "
              "below --identity-root and must contain complete PASS accounting"))
    parser.add_argument(
        "--identity-metadata-manifest", default=None,
        help=("Identity metadata_manifest.json; defaults below --identity-root"))
    parser.add_argument(
        "--identity-metadata-workbook", default=DEFAULT_IDENTITY_WORKBOOK,
        help=("Fixed Library_conversions.xlsx consumed by the upstream "
              "identity metadata stage"))
    parser.add_argument(
        "--identity-expected-genotypes", default=None,
        help=("Identity library_expected_genotypes.tsv; defaults below "
              "--identity-root"))
    parser.add_argument(
        "--identity-resolution-audit", default=None,
        help=("Identity library_resolution_audit.tsv; defaults below "
              "--identity-root"))
    parser.add_argument(
        "--identity-metadata-warnings", default=None,
        help=("Identity metadata_warnings.tsv; defaults below --identity-root"))
    parser.add_argument(
        "--identity-uid-members", default=None,
        help=("Identity library_uid_members.tsv; defaults below --identity-root"))
    parser.add_argument(
        "--identity-global-lines", default=None,
        help=("Identity global_biological_lines.tsv; defaults below "
              "--identity-root"))
    parser.add_argument(
        "--identity-global-donors", default=None,
        help=("Identity global_donors.tsv; defaults below --identity-root"))
    parser.add_argument(
        "--identity-metadata-helper",
        default=os.path.join(DEPLOYED_SCRIPTS, "identity_reconciliation.py"),
        help="Exact deployed identity_reconciliation.py metadata producer")
    parser.add_argument(
        "--panel-distinguishability-binary",
        default=os.path.join(DEPLOYED_BIN, "nuclear_panel_distinguishability"),
        help="Exact nuclear-panel distinguishability producer")
    parser.add_argument(
        "--ambient-profile-binary",
        default=os.path.join(DEPLOYED_BIN, "tet_ambient_profile"),
        help="Exact EMPTY_DROPS ambient-profile producer")
    parser.add_argument(
        "--geometry-gate-helper",
        default=os.path.join(DEPLOYED_SCRIPTS,
                             "geometry_gated_contam_estimate.py"),
        help="Exact geometry-gated ambient estimator wrapper")
    parser.add_argument(
        "--contam-binary",
        default=os.path.join(DEPLOYED_BIN, "tet_contam_estimate"),
        help="Exact contamination estimator invoked by the geometry wrapper")
    parser.add_argument(
        "--demux-pool-workbook", default=DEFAULT_DEMUX_POOL_WORKBOOK,
        help=("Fixed Library_conversions.xlsx used to generate DEMUX "
              "expected-lines files"))
    parser.add_argument(
        "--expected-pool-metadata", default=DEFAULT_EXPECTED_POOL_METADATA,
        help=("pool_combinations.tsv consumed by upstream POSTHOC identity "
              "auditing"))
    parser.add_argument(
        "--skip-identity-validation", action="store_true",
        help=("AUDIT OVERRIDE: bypass the required upstream identity boundary; "
              "recorded in the LEDGER job log"))
    parser.add_argument("--ledger-input", default=None)
    parser.add_argument("--final-assignments-template", default=None)
    parser.add_argument(
        "--mapping-bam-template", default=None,
        help=("Newest mapping BAM template used only to prove DEMUX lineage; "
              "defaults to <mapping-input-root>/Tet_..._{lib}/gex.bam"))
    parser.add_argument("--demux-prefix-template", default=None)
    parser.add_argument("--expression-barcodes-template", default=None,
                        help="MEX barcode path template; supports {lib} or {library}")
    parser.add_argument("--expression-features-template", default=None,
                        help="MEX feature path template; supports {lib} or {library}")
    parser.add_argument("--expression-matrix-template", default=None,
                        help="MEX matrix path template; supports {lib} or {library}")
    parser.add_argument("--ambient-standard-template", default=None)
    parser.add_argument("--ambient-arm-a-template", default=None)
    parser.add_argument("--ambient-arm-c-template", default=None)
    parser.add_argument(
        "--interindividual-panel", default=DEFAULT_INTERINDIVIDUAL_PANEL,
        help="Main donor-genotype BCF represented by the ASE pileup sidecars")
    parser.add_argument(
        "--het-panel", default=DEFAULT_HET_PANEL,
        help="HET BCF used upstream for ploidy diagnostics, not donor ASE")
    parser.add_argument(
        "--species-panel", default=DEFAULT_SPECIES_PANEL,
        help="Species BCF used upstream for species counting, not donor ASE")
    parser.add_argument(
        "--skip-upstream-provenance-check", action="store_true",
        help=("AUDIT OVERRIDE: skip mapping-run marker and forced-DEMUX-script "
              "lineage checks; input files and reconciliation schema are still "
              "validated and recorded"))
    parser.add_argument(
        "--allow-external-input-templates", action="store_true",
        help=("AUDIT OVERRIDE: permit explicit MEX/DEMUX/ambient/identity "
              "templates outside their declared roots"))
    parser.add_argument(
        "--allow-library-calibration-fallback", action="store_true",
        help=("AUDIT OVERRIDE: permit calibration grouping by library when no "
              "unique newest-run GEX cluster table exists"))
    parser.add_argument(
        "--ploidy-nn-template", default=None,
        help=("Optional PLOIDY_NN calls template; passed only when its expanded "
              "file exists and is nonempty"))
    parser.add_argument(
        "--ploidy-input-h5ad", default="",
        help=("Normalized newest-remap H5AD supplied to upstream PLOIDY_NN; "
              "required for non-production provenance validation"))
    parser.add_argument(
        "--ploidy-nn-weights", default=DEFAULT_PLOIDY_NN_WEIGHTS,
        help="Exact trained weights supplied to upstream PLOIDY_NN")
    parser.add_argument(
        "--ploidy-nn-helper", default=DEFAULT_PLOIDY_NN_HELPER,
        help="Exact run_ploidy_nn_inference.py used by upstream PLOIDY_NN")
    parser.add_argument(
        "--ploidy-nn-module", default="ploidy-inference/latest",
        help="Cluster module used for the runtime PLOIDY_NN reproduction gate")
    parser.add_argument("--ambient-condition", default=DEFAULT_CONDITION)
    parser.add_argument("--ambient-candidate-set", default="applied",
                        choices=("applied", "exploratory"))
    parser.add_argument(
        "--cell-groups-template", default="",
        help=("Explicit per-library cell-group template; otherwise discover one "
              "libN.*.tsv cluster artifact under the selected GEX analysis"))
    parser.add_argument(
        "--gex-ambient-analysis", default=DEFAULT_GEX_AMBIENT_ANALYSIS,
        help="GEX_AMBIENT analysis namespace used for cluster auto-discovery")
    parser.add_argument("--panel-metadata", default=DEFAULT_PANEL_METADATA)
    reference_input = parser.add_mutually_exclusive_group(required=True)
    reference_input.add_argument(
        "--arms-bed", default="",
        help="Existing BED4 in the same coordinates as demux pileup sites")
    reference_input.add_argument(
        "--gene-annotation", default="",
        help=("Ancestral-coordinate GTF/GFF; gene names are joined to "
              "--gene-arms and a BED is generated in the run root"))
    reference_input.add_argument(
        "--hal-file", default="",
        help=("Multispecies HAL containing the source genome and the exact "
              "ancestral target genome used by the pileup"))
    parser.add_argument("--gene-arms", default=DEFAULT_GENE_ARMS)
    parser.add_argument(
        "--reference-fai", default="",
        help=("Ancestral FASTA .fai with true contig lengths; required by "
              "the default anchored-contigs gene projection"))
    parser.add_argument(
        "--gene-projection-mode", choices=("anchored-contigs", "gene-spans"),
        default="anchored-contigs",
        help=("anchored-contigs partitions complete mapped contigs using ordered "
              "gene anchors; gene-spans includes only annotated gene spans"))
    parser.add_argument(
        "--arm-builder-script", default=DEFAULT_ARM_BUILDER_SCRIPT)
    parser.add_argument("--source-arms-bed", default=DEFAULT_SOURCE_ARMS)
    parser.add_argument("--hal-source-genome", default="Human")
    parser.add_argument("--hal-target-genome", default="human_chimp_bonobo")
    parser.add_argument(
        "--hal-reference-script", default=DEFAULT_HAL_REFERENCE_SCRIPT)

    parser.add_argument("--prepare-script", default=DEFAULT_PREPARE_SCRIPT)
    parser.add_argument("--expression-script", default=DEFAULT_EXPRESSION_SCRIPT)
    parser.add_argument("--ase-binary", default=DEFAULT_ASE_BINARY)
    parser.add_argument("--call-script", default=DEFAULT_CALL_SCRIPT)
    parser.add_argument("--report-script", default=DEFAULT_REPORT_SCRIPT)

    parser.add_argument("--partition", default=DEFAULT_PARTITION)
    parser.add_argument("--time", default=DEFAULT_TIME)
    parser.add_argument("--reference-memory", default="16G")
    parser.add_argument("--ledger-memory", default="256G")
    parser.add_argument("--prepare-memory", default="16G")
    parser.add_argument("--ase-memory", default="256G")
    parser.add_argument("--expression-memory", default="128G")
    parser.add_argument("--call-memory", default="256G")
    parser.add_argument("--report-memory", default="32G")
    parser.add_argument("--ase-threads", type=int, default=16)
    parser.add_argument("--ledger-cpus", type=int, default=16)
    parser.add_argument("--expression-cpus", type=int, default=8)
    parser.add_argument("--call-cpus", type=int, default=2)
    parser.add_argument("--report-cpus", type=int, default=2)

    call_model = parser.add_argument_group(
        "CALL model", "Arguments are passed explicitly to the global caller")
    for name, default in CALL_INTEGER_DEFAULTS:
        call_model.add_argument(
            "--" + name.replace("_", "-"), type=int, default=default)
    for name, default in CALL_FLOAT_DEFAULTS:
        call_model.add_argument(
            "--" + name.replace("_", "-"), type=float, default=default)

    parser.add_argument("--min-calibration-tet-probability", type=float, default=0.90)
    hard_threshold = parser.add_mutually_exclusive_group()
    hard_threshold.add_argument(
        "--hard-allele-threshold", dest="hard_allele_threshold", type=float,
        default=argparse.SUPPRESS,
        help="Intra-molecule weighted A/B support-fraction threshold [0.80]")
    hard_threshold.add_argument(
        "--hard-posterior", dest="hard_allele_threshold", type=float,
        default=argparse.SUPPRESS,
        help="Deprecated alias for --hard-allele-threshold")
    parser.set_defaults(hard_allele_threshold=0.80)
    parser.add_argument("--expression-pseudocount", type=float, default=0.5)
    parser.add_argument("--include-sex-chromosomes", action="store_true")
    parser.add_argument("--allow-missing-expression-barcodes", action="store_true")
    parser.add_argument("--include-nonheterotypic-expression", action="store_true")
    parser.add_argument(
        "--expression-ambient-corrected", action="store_true",
        help=("Mark the selected expression MEX as already ambient-corrected; "
              "does not transform counts"))
    parser.add_argument("--report-title",
                        default="Tetraploid chromosome-arm ASE/CNV report")
    parser.add_argument("--report-top-calls", type=int, default=100)
    return parser


def validate_args(args: argparse.Namespace) -> tuple[list[int], tuple[str, ...]]:
    args.runtime_ploidy_output_root = ""
    args.runtime_ploidy_libraries = []
    libraries = parse_libraries(args.libraries)
    upstream_cohort = (
        parse_libraries(args.upstream_cohort_libraries)
        if args.upstream_cohort_libraries else list(libraries))
    if not set(libraries) <= set(upstream_cohort):
        raise ValueError(
            "--upstream-cohort-libraries must contain every output library")
    args.upstream_cohort_libraries_resolved = upstream_cohort
    stages = parse_stages(args.stages)
    if args.array_throttle is not None and args.array_throttle < 1:
        raise ValueError("--array-throttle must be positive")
    for name in (
            "ledger_cpus", "ase_threads", "expression_cpus", "call_cpus",
            "report_cpus"):
        if not 1 <= int(getattr(args, name)) <= 256:
            raise ValueError(f"--{name.replace('_', '-')} must be between 1 and 256")
    if not 0.5 < args.hard_allele_threshold < 1.0:
        raise ValueError(
            "--hard-allele-threshold must be strictly between 0.5 and 1")
    for name, _default in CALL_INTEGER_DEFAULTS:
        if getattr(args, name) < 1:
            raise ValueError(f"--{name.replace('_', '-')} must be positive")
    probability_names = (
        "event_prior", "max_event_q", "max_whole_chromosome_q",
        "min_event_posterior", "min_balanced_posterior",
        "min_uid_arm_posterior", "max_qname_fallback_fraction",
        "min_ambient_genotyped_mass", "min_rho", "max_rho", "default_rho",
        "min_uid_state_concordance", "min_pair_state_concordance",
        "min_pair_cell_posterior", "min_pair_arm_posterior",
        "max_pair_arm_q", "empirical_ambient_window",
    )
    for name, _default in CALL_FLOAT_DEFAULTS:
        if not math.isfinite(getattr(args, name)):
            raise ValueError(f"--{name.replace('_', '-')} must be finite")
    for name in probability_names:
        if not 0.0 <= getattr(args, name) <= 1.0:
            raise ValueError(f"--{name.replace('_', '-')} must be in [0,1]")
    if not 0.0 < args.event_prior < 1.0:
        raise ValueError("--event-prior must be strictly between zero and one")
    if args.min_rho <= 0.0 or args.max_rho <= args.min_rho:
        raise ValueError("invalid CALL quasi-likelihood rho bounds")
    if not args.min_rho <= args.default_rho <= args.max_rho:
        raise ValueError("--default-rho must be within the CALL rho bounds")
    if args.absolute_min_null_cells > args.min_empirical_null_cells:
        raise ValueError(
            "--absolute-min-null-cells cannot exceed --min-empirical-null-cells")
    if args.calibration_crossfit_folds < 2:
        raise ValueError("--calibration-crossfit-folds must be at least two")
    if args.min_pair_recurrence_cells > args.max_pair_test_cells:
        raise ValueError(
            "--min-pair-recurrence-cells cannot exceed --max-pair-test-cells")
    if args.min_expression_counts < 0.0:
        raise ValueError("--min-expression-counts cannot be negative")
    if not 0.0 < args.min_expression_sigma <= args.max_expression_sigma:
        raise ValueError("invalid CALL expression sigma bounds")
    if args.expression_weight < 0.0 or args.max_expression_log_bf < 0.0:
        raise ValueError("CALL expression weights must be nonnegative")
    if not 0.0 < args.site_fallback_likelihood_weight < 1.0:
        raise ValueError("--site-fallback-likelihood-weight must be in (0,1)")
    if args.empirical_depth_fold <= 1.0:
        raise ValueError("--empirical-depth-fold must exceed one")
    if (args.min_call_effective_weight <= 0.0
            or args.min_calibration_effective_weight <= 0.0
            or args.min_orientation_calibration_effective_weight <= 0.0):
        raise ValueError("CALL effective-weight thresholds must be positive")
    if not 0.0 < args.calibration_min_robust_weight <= 1.0:
        raise ValueError("--calibration-min-robust-weight must be in (0,1]")
    if (args.calibration_shrinkage_cells < 0.0
            or args.calibration_max_cell_weight <= 0.0
            or args.calibration_huber_z <= 0.0
            or args.calibration_expression_mad_cutoff <= 0.0
            or args.cell_baseline_shrinkage_arms <= 0.0
            or args.max_aggregate_cell_log_bf <= 0.0):
        raise ValueError("invalid CALL calibration tuning value")
    if not 0.0 <= args.min_calibration_tet_probability <= 1.0:
        raise ValueError("--min-calibration-tet-probability must be in [0,1]")
    if args.expression_pseudocount <= 0:
        raise ValueError("--expression-pseudocount must be positive")
    if not 1 <= args.report_top_calls <= 10000:
        raise ValueError("--report-top-calls must be between 1 and 10000")
    if not re.fullmatch(r"[A-Za-z0-9_.-]+", args.partition):
        raise ValueError("--partition contains unsupported characters")
    if not re.fullmatch(r"[A-Za-z0-9_.-]+", args.ambient_condition):
        raise ValueError("--ambient-condition contains unsupported characters")
    if (not args.skip_upstream_provenance_check and
            args.ambient_condition != DEFAULT_CONDITION):
        raise ValueError(
            "production provenance currently requires the audited ambient "
            f"condition {DEFAULT_CONDITION}")
    if (not args.skip_upstream_provenance_check and
            args.ambient_candidate_set != "applied"):
        raise ValueError(
            "production provenance requires --ambient-candidate-set applied; "
            "exploratory reconciliation cannot drive final arm-CNV calls")
    if not re.fullmatch(r"[A-Za-z0-9_.-]+", args.gex_ambient_analysis):
        raise ValueError("--gex-ambient-analysis contains unsupported characters")
    if not re.fullmatch(r"[A-Za-z0-9_./+-]+", args.ploidy_nn_module):
        raise ValueError("--ploidy-nn-module contains unsupported characters")
    if (not args.skip_upstream_provenance_check and
            set(stages) & {"PREPARE", "ASE", "EXPRESSION"} and
            "LEDGER" not in stages):
        raise ValueError(
            "strict scientific stage selections must include LEDGER so the "
            "runtime panel and PLOIDY_NN reproduction gates cannot be bypassed")
    if not re.fullmatch(r"[0-9]+-[0-9]{2}:[0-9]{2}:[0-9]{2}|[0-9]{1,3}:[0-9]{2}:[0-9]{2}", args.time):
        raise ValueError("--time must use D-HH:MM:SS or HH:MM:SS")
    for name in (
            "reference_memory", "ledger_memory", "prepare_memory", "ase_memory",
            "expression_memory", "call_memory", "report_memory"):
        if not re.fullmatch(r"[1-9][0-9]*(?:[KMGTP])?", getattr(args, name)):
            raise ValueError(
                f"--{name.replace('_', '-')} must be a positive SLURM memory token")

    configure_input_roots(args)
    default_templates(args)
    if args.gene_annotation:
        args.reference_mode = "GENE_SYNTENY"
        args.gene_annotation = absolute(args.gene_annotation)
        args.hal_file = ""
        args.reference_fai = (
            absolute(args.reference_fai) if args.reference_fai else "")
        if args.gene_projection_mode == "anchored-contigs" and not args.reference_fai:
            raise ValueError(
                "--reference-fai is required with --gene-projection-mode "
                "anchored-contigs")
        args.arms_bed = run_paths(args.run_root).generated_arms_bed
    elif args.hal_file:
        args.reference_mode = "HAL_LIFTOVER"
        args.gene_annotation = ""
        args.reference_fai = ""
        args.hal_file = absolute(args.hal_file)
        args.arms_bed = run_paths(args.run_root).generated_arms_bed
    else:
        args.reference_mode = "EXPLICIT_BED"
        args.gene_annotation = ""
        args.reference_fai = ""
        args.hal_file = ""
        args.arms_bed = absolute(args.arms_bed)
    for field in (
        "run_root", "mapping_input_root", "upstream_analysis_root",
        "identity_root", "identity_validation",
        "identity_validation_failures", "identity_run_summary",
        "identity_metadata_manifest", "identity_metadata_workbook",
        "identity_expected_genotypes", "identity_resolution_audit",
        "identity_metadata_warnings", "identity_uid_members",
        "identity_global_lines", "identity_global_donors",
        "identity_metadata_helper", "panel_distinguishability_binary",
        "ambient_profile_binary", "geometry_gate_helper", "contam_binary",
        "demux_pool_workbook", "expected_pool_metadata",
        "ledger_input",
        "arms_bed", "gene_arms", "arm_builder_script", "source_arms_bed",
        "hal_reference_script", "prepare_script",
        "expression_script", "ase_binary", "call_script", "report_script",
        "interindividual_panel", "het_panel", "species_panel",
        "ploidy_nn_weights", "ploidy_nn_helper"):
        setattr(args, field, absolute(getattr(args, field)))
    args.mapping_run_root = (
        absolute(args.mapping_run_root) if args.mapping_run_root else "")
    args.ploidy_input_h5ad = (
        absolute(args.ploidy_input_h5ad) if args.ploidy_input_h5ad else "")
    args.panel_metadata = (
        absolute(args.panel_metadata) if args.panel_metadata else "")
    validate_text(args.report_title, "--report-title")
    for field in (
            "run_root", "mapping_run_root", "mapping_input_root",
            "upstream_analysis_root",
            "identity_root", "identity_validation",
            "identity_validation_failures", "identity_run_summary",
            "identity_metadata_manifest", "identity_metadata_workbook",
            "identity_expected_genotypes", "identity_resolution_audit",
            "identity_metadata_warnings", "identity_uid_members",
            "identity_global_lines", "identity_global_donors",
            "identity_metadata_helper", "panel_distinguishability_binary",
            "ambient_profile_binary", "geometry_gate_helper", "contam_binary",
            "demux_pool_workbook", "expected_pool_metadata",
            "ledger_input",
            "panel_metadata", "ploidy_input_h5ad", "arms_bed",
            "gene_annotation", "gene_arms",
            "reference_fai", "hal_file", "source_arms_bed",
            "hal_source_genome", "hal_target_genome",
            "arm_builder_script", "hal_reference_script", "prepare_script",
            "expression_script", "ase_binary", "call_script", "report_script",
            "interindividual_panel", "het_panel", "species_panel",
            "ploidy_nn_weights", "ploidy_nn_helper"):
        validate_text(getattr(args, field), f"--{field.replace('_', '-')}")
    if any(character.isspace() for character in args.run_root):
        raise ValueError("--run-root cannot contain whitespace in SLURM log paths")
    for field in (
        "final_assignments_template", "mapping_bam_template",
        "demux_prefix_template",
        "expression_barcodes_template", "expression_features_template",
        "expression_matrix_template", "ambient_standard_template",
        "ambient_arm_a_template", "ambient_arm_c_template",
        "ploidy_nn_template", "cell_groups_template"):
        validate_text(getattr(args, field), f"--{field.replace('_', '-')}")
    return libraries, stages


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        libraries, selected = validate_args(args)
        run = run_paths(args.run_root)
        stage_directories = {
            "REFERENCE": run.reference,
            "LEDGER": run.ledger,
            "PREPARE": run.prepare,
            "ASE": run.ase,
            "EXPRESSION": run.expression,
            "CALL": run.call,
            "REPORT": run.report,
        }
        directories = [run.root, run.logs, run.scripts]
        if any(stage in selected for stage in ("PREPARE", "ASE", "EXPRESSION")):
            directories.append(run.manifests)
        directories.extend(stage_directories[stage] for stage in selected)
        for directory in dict.fromkeys(directories):
            os.makedirs(directory, exist_ok=True)
        paths = make_library_paths(args, run, libraries, selected)
        args.downstream_output_libraries = list(libraries)
        upstream_cohort = list(args.upstream_cohort_libraries_resolved)
        if upstream_cohort == list(libraries):
            audit_paths = paths
        else:
            output_paths_by_library = {item.library: item for item in paths}
            cohort_paths = make_library_paths(
                args, run, upstream_cohort,
                tuple(stage for stage in selected if stage != "PREPARE"))
            audit_paths = [
                output_paths_by_library.get(item.library, item)
                for item in cohort_paths]
        if (args.reference_mode != "EXPLICIT_BED" and
                "REFERENCE" not in selected and not regular_nonempty(args.arms_bed)):
            raise ValueError(
                "the selected reference mode requires the REFERENCE stage unless the generated "
                f"BED already exists: {args.arms_bed}")
        reject_output_conflicts(selected, run, paths, args)
        input_provenance = ""
        if set(selected) & {"LEDGER", "PREPARE", "ASE", "EXPRESSION"}:
            input_provenance = validate_and_record_input_bundle(
                args, run, audit_paths, upstream_cohort, selected)
        qc_files = report_qc_files(run, paths)
        manifests = prepare_manifests(args, run, paths, selected)

        identity_ok = True
        identity_detail = "not checked (LEDGER not selected)"
        if "LEDGER" in selected:
            identity_ok, identity_detail = identity_validation_summary_passes(
                args.identity_validation, upstream_cohort)
            if args.submit and not args.skip_identity_validation and not identity_ok:
                raise ValueError(
                    "identity validation boundary is not PASS: "
                    f"{args.identity_validation} ({identity_detail}); use "
                    "--skip-identity-validation only as an explicit audited override")

        scripts: dict[str, str] = {}
        if "REFERENCE" in selected:
            scripts["REFERENCE"] = render_script(
                os.path.join(run.scripts, "00_reference.sbatch"),
                reference_script(args, run))
        if "LEDGER" in selected:
            scripts["LEDGER"] = render_script(
                os.path.join(run.scripts, "01_ledger.sbatch"),
                ledger_script(args, run, libraries))
        if "PREPARE" in selected:
            scripts["PREPARE"] = render_script(
                os.path.join(run.scripts, "02_prepare_array.sbatch"),
                prepare_script(args, run, manifests["PREPARE"], len(paths)))
        if "ASE" in selected:
            scripts["ASE"] = render_script(
                os.path.join(run.scripts, "03_ase_array.sbatch"),
                ase_script(args, run, manifests["ASE"], len(paths)))
        if "EXPRESSION" in selected:
            scripts["EXPRESSION"] = render_script(
                os.path.join(run.scripts, "04_expression_array.sbatch"),
                expression_script(
                    args, run, manifests["EXPRESSION"], len(paths)))
        if "CALL" in selected:
            scripts["CALL"] = render_script(
                os.path.join(run.scripts, "05_call.sbatch"),
                call_script(args, run, paths))
        if "REPORT" in selected:
            scripts["REPORT"] = render_script(
                os.path.join(run.scripts, "06_report.sbatch"),
                report_script(args, run, qc_files))

        print(f"Tetraploid Arm CNV orchestrator {RELEASE}")
        print(f"Run root: {run.root}")
        print(f"Mapping input root: {args.mapping_input_root}")
        print(f"Upstream analysis root: {args.upstream_analysis_root}")
        if input_provenance:
            print(f"Matched input provenance: {input_provenance}")
        print("Output libraries: " + ",".join(str(value) for value in libraries))
        print("Validated upstream cohort: " + ",".join(
            str(value) for value in upstream_cohort))
        print("Selected stages: " + ",".join(selected))
        print(f"Arm reference mode: {args.reference_mode}")
        if args.reference_mode == "GENE_SYNTENY":
            print(f"Gene projection mode: {args.gene_projection_mode}")
        elif args.reference_mode == "HAL_LIFTOVER":
            print(
                f"HAL projection: {args.hal_source_genome} -> "
                f"{args.hal_target_genome}")
        print(f"Arm BED consumed by ASE: {args.arms_bed}")
        if manifests:
            print("Array task maps:")
            for stage in ("PREPARE", "ASE", "EXPRESSION"):
                if stage in manifests:
                    print(f"  {stage:<10} {manifests[stage]}")
        else:
            print("Array task maps: none (no array stage selected)")
        if "LEDGER" in selected:
            print(
                "Identity boundary: " +
                ("SKIPPED (AUDIT OVERRIDE)"
                 if args.skip_identity_validation else
                 (identity_detail if identity_ok
                  else f"NOT PASS ({identity_detail})")) +
                f" [{args.identity_validation}]")
        else:
            print("Identity boundary: not checked (LEDGER not selected)")
        if "PREPARE" in selected:
            fallback = [path.library for path in paths if not path.cell_groups]
            for library in fallback:
                print(
                    f"WARNING: lib{library} has no unique auto-discovered GEX "
                    f"cluster artifact; calibration group falls back to lib{library}")
        for stage in STAGES:
            state = "planned" if stage in selected else "not selected"
            print(f"  {stage:<10} {state}")
            if stage in scripts:
                print(f"              script: {scripts[stage]}")

        if not args.submit:
            print("Planning only: no jobs submitted. Re-run with --submit to launch.")
            return 0

        job_ids: dict[str, str] = {}
        for stage in STAGES:
            if stage not in scripts:
                continue
            dependencies = resolved_dependencies(stage, selected, job_ids)
            job_id = submit_job(scripts[stage], dependencies, run.root)
            job_ids[stage] = job_id
            dependency_text = ":".join(dependencies) if dependencies else "none"
            print(
                f"Submitted {stage:<10} job {job_id}; afterok={dependency_text}; "
                f"logs={run.logs}")
        return 0
    except (OSError, ValueError, RuntimeError, subprocess.SubprocessError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
