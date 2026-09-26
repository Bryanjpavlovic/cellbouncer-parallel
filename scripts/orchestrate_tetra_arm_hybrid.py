#!/usr/bin/env python3
"""Plan and optionally submit the complete hybrid arm-CNV workflow.

This is the authoritative, standalone orchestrator for the hybrid expression
and CellBouncer ASE workflow.  It owns input discovery, the nine-stage DAG,
job rendering, immutable planning, recovery validation, and safe scheduler
submission.  Scientific calculations remain in their dedicated workers.

Planning is the default and never calls ``sbatch``.  Re-run the same command
with ``--submit`` only after inspecting the plan and generated scripts.
"""

from __future__ import annotations

import argparse
import csv
import fcntl
import glob
import gzip
import json
import math
import os
import re
import shlex
import stat
import subprocess
import sys
import tempfile
from dataclasses import dataclass, replace
from typing import Iterable, Mapping, Sequence


RELEASE = "1.0.2"
ORCHESTRATOR_NAME = "orchestrate_tetra_arm_hybrid.py"
WORKFLOW = "arm-cnv-hybrid"


HYBRID_STAGES = (
    "REFERENCE", "LEDGER", "PREPARE", "ASE", "EXPRESSION",
    "HYBRID_SHARD", "MODEL", "CALL", "REPORT",
)

STAGES = HYBRID_STAGES


HYBRID_STAGE_PARENTS = {
    "REFERENCE": (),
    "LEDGER": (),
    "PREPARE": ("LEDGER",),
    "ASE": ("REFERENCE", "PREPARE"),
    "EXPRESSION": ("PREPARE",),
    "HYBRID_SHARD": ("ASE", "EXPRESSION"),
    "MODEL": ("HYBRID_SHARD",),
    "CALL": ("MODEL",),
    "REPORT": ("CALL",),
}


PROJECT_ROOT = "/mnt/beegfs/tetraploid_multiome_cis_trans"


CURRENT_RNA_RUN_ROOT = os.path.join(PROJECT_ROOT, "3P")


MAPPING_ROOT = os.path.join(CURRENT_RNA_RUN_ROOT, "mapping_output")


PRODUCTION_ANALYSIS_ROOT = os.path.join(CURRENT_RNA_RUN_ROOT, "analysis")


DEFAULT_PANEL_METADATA = os.path.join(
    PROJECT_ROOT, "Misc_Metadata", "panel_metadata.tsv")


DEFAULT_CONDITION = "IND_CK_RF_SX0_GATED_RFREE_PFIT"


DEFAULT_GEX_AMBIENT_ANALYSIS = "full_gene_rna_leiden_v1"


LIBRARY_PREFIX = "Tet_2025_Multiome-RNA_"


DEMUX_SUBDIR = "demux_nomito"


NOMITO_PANEL_ROOT = (
    "/mnt/beegfs/home/b/vcfdownsample/"
    "Downsample_ATAC_Species_poolInformer/NoMito")


DEFAULT_INTERINDIVIDUAL_PANEL = os.path.join(
    NOMITO_PANEL_ROOT, "tet.vars.downsampled_20M.bcf")


DEFAULT_PLOIDY_NN_WEIGHTS = (
    "/mnt/beegfs/tetmultiome_rna_mapped/ploidy_classifier/"
    "ploidy_nn_weights.pt"
)


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


DEFAULT_REPORT_SCRIPT = os.path.join(DEPLOYED_SCRIPTS, "tetra_arm_report.py")


DEFAULT_HYBRID_SHARD_SCRIPT = os.path.join(
    DEPLOYED_SCRIPTS, "tetra_arm_hybrid_shard.py")


DEFAULT_HYBRID_MODEL_SCRIPT = os.path.join(
    DEPLOYED_SCRIPTS, "tetra_arm_hybrid_model.py")


DEFAULT_HYBRID_AGGREGATE_SCRIPT = os.path.join(
    DEPLOYED_SCRIPTS, "tetra_arm_hybrid_aggregate.py")


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


EXPRESSION_HEADER = (
    "task_index", "library", "barcodes", "features", "matrix",
    "cell_manifest", "gene_arms", "output_dir",
)


HYBRID_INPUT_HEADER = (
    "library", "hybrid_cell_manifest", "ase", "expression",
    "expression_model",
)


HYBRID_SHARD_TASK_HEADER = (
    "task_index", "library", "input_manifest", "output_directory", "qc",
    "contract",
)


HYBRID_MODEL_TASK_HEADER = (
    "task_index", "chromosome", "logical_arms", "shard_manifest",
    "output_directory", "qc", "contract",
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
    hybrid_cell_manifest: str
    expression_model: str


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


@dataclass(frozen=True)
class HybridRunPaths:
    """Isolated run layout for the opt-in hybrid-v1 workflow."""

    root: str
    logs: str
    scripts: str
    manifests: str
    reference: str
    ledger: str
    prepare: str
    ase: str
    expression: str
    shard: str
    model: str
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
    def input_provenance(self) -> str:
        # Hybrid reuse is deliberately guarded by a non-hash source manifest.
        # An empty value disables the legacy SHA/stat execution guard in the
        # shared SLURM header even if an unrelated file exists in the run root.
        return ""

    @property
    def source_reuse_manifest(self) -> str:
        return os.path.join(self.manifests, "hybrid_source_reuse.json")

    @property
    def hybrid_inputs(self) -> str:
        return os.path.join(self.manifests, "hybrid_inputs.tsv")

    @property
    def shard_tasks(self) -> str:
        return os.path.join(self.manifests, "hybrid_shard_tasks.tsv")

    @property
    def model_tasks(self) -> str:
        return os.path.join(self.manifests, "hybrid_model_tasks.tsv")

    @property
    def call_outputs(self) -> tuple[str, ...]:
        return tuple(self.call_prefix + suffix for suffix in (
            ".arm_calls.tsv.gz",
            ".calibration.tsv.gz",
            ".expression_components.tsv.gz",
            ".uid_chromosome_flags.tsv.gz",
            ".donor_pair_arm_summary.tsv.gz",
            ".qc.tsv",
            ".contract.json",
        ))

    @property
    def report_outputs(self) -> tuple[str, ...]:
        return (
            os.path.join(self.report, "tetra_arm_cnv_hybrid_summary.tsv"),
            os.path.join(self.report, "tetra_arm_cnv_hybrid_summary.json"),
            os.path.join(self.report, "tetra_arm_cnv_hybrid_report.html"),
        )


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


def parse_hybrid_stages(values: Sequence[str] | None) -> tuple[str, ...]:
    if not values:
        return HYBRID_STAGES
    selected: set[str] = set()
    for raw in values:
        for token in str(raw).split(","):
            stage = token.strip().upper()
            if not stage:
                continue
            if stage == "ALL":
                selected.update(HYBRID_STAGES)
            elif stage in HYBRID_STAGES:
                selected.add(stage)
            else:
                raise ValueError(
                    f"unknown hybrid stage {token!r}; choose from "
                    f"{', '.join(HYBRID_STAGES)}")
    if not selected:
        raise ValueError("at least one hybrid stage must be selected")
    return tuple(stage for stage in HYBRID_STAGES if stage in selected)


def absolute(value: str) -> str:
    return os.path.abspath(os.path.expanduser(value))


def validate_text(value: str, label: str) -> str:
    if "\n" in value or "\r" in value or "\t" in value:
        raise ValueError(f"{label} contains a tab or newline")
    return value


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


def hybrid_run_paths(root: str) -> HybridRunPaths:
    root = absolute(root)
    return HybridRunPaths(
        root=root,
        logs=os.path.join(root, "logs"),
        scripts=os.path.join(root, "generated_scripts"),
        manifests=os.path.join(root, "task_manifests"),
        reference=os.path.join(root, "reference"),
        ledger=os.path.join(root, "ledger"),
        prepare=os.path.join(root, "prepare"),
        ase=os.path.join(root, "ase"),
        expression=os.path.join(root, "expression"),
        shard=os.path.join(root, "hybrid_shard"),
        model=os.path.join(root, "model"),
        call=os.path.join(root, "call"),
        report=os.path.join(root, "report"),
        call_prefix=os.path.join(root, "call", "tetra_arm_cnv_hybrid"),
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

    The canonical layout keeps mapping and analysis as separate siblings under
    the physical 3P root.  A non-production remap must name its matching
    isolated analysis root so MEX files cannot be mixed with demux, ambient, or
    reconciled-identity products from another run.
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
    args.identity_metadata_manifest = absolute(
        args.identity_metadata_manifest or os.path.join(
            identity_root, "metadata", "metadata_manifest.json"))
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
            hybrid_cell_manifest=os.path.join(
                run.prepare, f"lib{library}.hybrid_cell_manifest.tsv.gz"),
            expression_model=os.path.join(
                run.expression, f"lib{library}.arm_expression_model.tsv.gz"),
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
    nodelist = (f"#SBATCH --nodelist={args.nodelist}\n"
                if args.nodelist else "")
    return f"""#!/bin/bash
# Generated deterministically by {ORCHESTRATOR_NAME} {RELEASE}.
# Generated script and task map are immutable execution inputs.
# Change orchestrator options and use a new --run-root instead of editing them.
#SBATCH --job-name=tetarm_{stage_lower}
#SBATCH --output={run.logs}/tetarm_{stage_lower}_{log_token}.out
#SBATCH --error={run.logs}/tetarm_{stage_lower}_{log_token}.err
#SBATCH --partition={args.partition}
{nodelist}#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task={cpus}
#SBATCH --mem={memory}
#SBATCH --time={args.time}
{array}
set -eo pipefail

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
import math
import os
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
    if fields != ["check", "status", "n_failures", "detail"] or not rows:
        raise ValueError("invalid validation summary schema")
    observed = set()
    for row in rows:
        check = str(row.get("check", "")).strip()
        detail = str(row.get("detail", "")).strip()
        if not check:
            raise ValueError("empty validation check name")
        key = (check, detail)
        if key in observed:
            raise ValueError("duplicate identity validation row")
        observed.add(key)
        failures_value = float(str(row.get("n_failures", "")).strip())
        failures = (
            int(failures_value)
            if math.isfinite(failures_value) and
            failures_value.is_integer() else -1)
        if (not math.isfinite(failures_value)
                or failures_value != failures
                or str(row.get("status", "")).strip().upper() != "PASS"
                or failures != 0):
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
                    {{"library", "barcode", "assignment_status",
                      "final_assignment"}} <= set(header)):
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


def read_json_object(path: str, label: str) -> dict[str, object]:
    if not regular_nonempty(path):
        raise ValueError(f"completed {label} is missing or empty: {path}")
    try:
        with open(path, "r", encoding="utf-8") as handle:
            payload = json.load(handle)
    except (OSError, ValueError, TypeError, json.JSONDecodeError) as exc:
        raise ValueError(f"completed {label} is not valid JSON: {path}: {exc}") \
            from exc
    if not isinstance(payload, dict):
        raise ValueError(f"completed {label} is not a JSON object: {path}")
    return payload


def read_metric_map(path: str, label: str) -> dict[str, str]:
    if not regular_nonempty(path):
        raise ValueError(f"completed {label} is missing or empty: {path}")
    try:
        with open(path, "r", encoding="utf-8", newline="") as handle:
            reader = csv.reader(handle, delimiter="\t")
            header = next(reader)
            if header != ["metric", "value"]:
                raise ValueError("expected metric/value header")
            metrics: dict[str, str] = {}
            for row in reader:
                if len(row) != 2 or not row[0] or row[0] in metrics:
                    raise ValueError("invalid or duplicate metric row")
                metrics[row[0]] = row[1]
    except (OSError, ValueError, StopIteration, csv.Error) as exc:
        raise ValueError(f"completed {label} is invalid: {path}: {exc}") from exc
    return metrics


def read_tsv_header(path: str, label: str, compressed: bool = False) -> list[str]:
    if not regular_nonempty(path):
        raise ValueError(f"completed {label} is missing or empty: {path}")
    opener = gzip.open if compressed else open
    try:
        with opener(path, "rt", encoding="utf-8", newline="") as handle:
            header = next(csv.reader(handle, delimiter="\t"))
    except (OSError, EOFError, ValueError, StopIteration, csv.Error) as exc:
        raise ValueError(f"completed {label} has no valid TSV header: {path}: {exc}") \
            from exc
    if not header or len(header) != len(set(header)):
        raise ValueError(f"completed {label} has an empty or duplicate TSV header: {path}")
    return header


def require_contract(
        path: str, label: str, schema: str,
        allowed_statuses: Sequence[str] = ("PASS",)) -> dict[str, object]:
    payload = read_json_object(path, label)
    status = str(payload.get("status", "")).upper()
    if payload.get("schema_version") != schema or status not in allowed_statuses:
        raise ValueError(
            f"completed {label} contract is not valid: {path}; "
            f"schema={payload.get('schema_version')!r}, status={status!r}")
    return payload


def validate_completed_stage(
        stage: str, run: RunPaths, paths: Sequence[LibraryPaths],
        args: argparse.Namespace) -> None:
    """Validate a completed ancestor bundle before it is reused by resume."""
    missing = [
        path for path in stage_primary_outputs(stage, run, paths, args)
        if not regular_nonempty(path)
    ]
    if missing:
        preview = ", ".join(missing[:3])
        suffix = " ..." if len(missing) > 3 else ""
        raise ValueError(
            f"cannot resume: completed {stage} bundle is missing or empty: "
            f"{preview}{suffix}")

    if stage == "REFERENCE":
        schema_by_mode = {
            "EXPLICIT_BED": "tetra_arm_reference_contract_v1",
            "GENE_SYNTENY": "tetra_arm_gene_synteny_contract_v1",
            "HAL_LIFTOVER": "tetra_arm_hal_liftover_contract_v1",
        }
        contract = require_contract(
            run.reference_contract, "REFERENCE",
            schema_by_mode[args.reference_mode])
        output_field = "arms_bed" if args.reference_mode == "EXPLICIT_BED" else "output"
        if absolute(str(contract.get(output_field, ""))) != absolute(args.arms_bed):
            raise ValueError(
                "cannot resume: REFERENCE contract does not name the selected "
                f"arm BED: {run.reference_contract}")
        return

    if stage == "LEDGER":
        contract_path = os.path.join(run.ledger, "split_ledger_contract.json")
        contract = require_contract(
            contract_path, "LEDGER", "tetra_arm_split_ledger_contract_v1")
        expected_libraries = [path.library for path in paths]
        observed = [int(value) for value in contract.get("libraries", [])]
        if observed != expected_libraries or int(contract.get("cells", 0)) < 1:
            raise ValueError(
                f"cannot resume: LEDGER contract does not match output "
                f"libraries {expected_libraries}: {contract_path}")
        for item in paths:
            fields = set(read_tsv_header(
                item.split_ledger, f"LEDGER lib{item.library}", compressed=True))
            if not {"library", "barcode", "assignment_status",
                    "final_assignment"} <= fields:
                raise ValueError(
                    f"cannot resume: LEDGER header is incomplete: {item.split_ledger}")
        return

    if stage == "PREPARE":
        for item in paths:
            contract = require_contract(
                item.prepare_contract, f"PREPARE lib{item.library}",
                "tetra_arm_prepare_contract_v1")
            if (int(contract.get("library", -1)) != item.library
                    or int(contract.get("cells", 0)) < 1):
                raise ValueError(
                    f"cannot resume: PREPARE contract is inconsistent: "
                    f"{item.prepare_contract}")
            cell_fields = set(read_tsv_header(
                item.cell_manifest, f"PREPARE lib{item.library} cell manifest",
                compressed=True))
            ambient_fields = set(read_tsv_header(
                item.ambient_sources,
                f"PREPARE lib{item.library} ambient sources", compressed=True))
            if not {"library", "barcode", "donor_a", "donor_b", "donor_pair",
                    "model_eligible", "schema_version"} <= cell_fields:
                raise ValueError(
                    f"cannot resume: PREPARE cell-manifest header is incomplete: "
                    f"{item.cell_manifest}")
            if not {"library", "barcode", "source_label",
                    "scoring_profile_mass", "schema_version"} <= ambient_fields:
                raise ValueError(
                    f"cannot resume: PREPARE ambient-source header is incomplete: "
                    f"{item.ambient_sources}")
        return

    if stage == "ASE":
        allowed = {
            "PASS", "PASS_NO_HETEROTYPIC_TARGETS",
            "PASS_NO_INFORMATIVE_EVIDENCE",
        }
        for item in paths:
            header = tuple(read_tsv_header(
                item.ase, f"ASE lib{item.library}", compressed=True))
            metrics = read_metric_map(item.ase_qc, f"ASE lib{item.library} QC")
            if (header != ASE_V2_HEADER
                    or metrics.get("schema_version") != "tetra_arm_ase_qc_v2"
                    or metrics.get("status", "").upper() not in allowed):
                raise ValueError(
                    f"cannot resume: ASE bundle is not schema-valid: {item.ase}")
        return

    if stage == "EXPRESSION":
        expected = [
            "library", "barcode", "arm", "chromosome", "arm_counts",
            "other_autosomal_counts", "reference_autosomal_counts",
            "total_autosomal_counts", "arm_fraction", "log2_arm_to_other",
            "log2_arm_to_reference", "mapped_genes_on_arm",
            "nonzero_genes_on_arm", "matrix_value_type",
            "expression_input_state", "schema_version",
        ]
        for item in paths:
            contract = require_contract(
                item.expression_contract, f"EXPRESSION lib{item.library}",
                "tetra_arm_expression_contract_v1")
            if int(contract.get("library", -1)) != item.library:
                raise ValueError(
                    f"cannot resume: EXPRESSION contract has the wrong library: "
                    f"{item.expression_contract}")
            if read_tsv_header(
                    item.expression, f"EXPRESSION lib{item.library}",
                    compressed=True) != expected:
                raise ValueError(
                    f"cannot resume: EXPRESSION header is invalid: {item.expression}")
        return

    if stage == "CALL":
        allowed = (
            "PASS", "PASS_NO_HETEROTYPIC_TARGETS", "PASS_NO_OBSERVED_ASE",
            "PASS_NO_CALLABLE_ASE")
        contract = require_contract(
            run.call_prefix + ".contract.json", "CALL",
            "tetra_arm_call_contract_v2", allowed)
        schemas = contract.get("output_schemas", {})
        if (not isinstance(schemas, dict)
                or schemas.get("calls") != "tetra_arm_cnv_calls_v2"
                or schemas.get("calibration") != "tetra_arm_calibration_v2"
                or schemas.get("uid_chromosome_flags")
                   != "tetra_arm_uid_chromosome_flags_v3"
                or schemas.get("donor_pair_arm_summary")
                   != "tetra_arm_donor_pair_arm_summary_v3"):
            raise ValueError(
                f"cannot resume: CALL output schemas are invalid: "
                f"{run.call_prefix}.contract.json")
        return

    raise ValueError(f"unsupported resume prerequisite stage: {stage}")


def _hybrid_natural_key(value: object) -> tuple[object, ...]:
    return tuple(
        int(part) if part.isdigit() else part.lower()
        for part in re.split(r"(\d+)", str(value)))


def _hybrid_chromosome(value: object) -> str:
    return str(value).strip().removeprefix("chr").removeprefix("Chr")


def hybrid_chromosome_arms(
        gene_arms_path: str, include_sex_chromosomes: bool
        ) -> dict[str, tuple[str, ...]]:
    """Discover logical p/q arms from the configured reference map."""
    path = absolute(gene_arms_path)
    if not regular_nonempty(path):
        raise ValueError(f"gene-arm map is missing or empty: {path}")
    opener = gzip.open if path.endswith(".gz") else open
    result: dict[str, set[str]] = {}
    try:
        with opener(path, "rt", encoding="utf-8", newline="") as handle:
            for line_number, raw in enumerate(handle, start=1):
                fields = raw.rstrip("\r\n").split("\t")
                if not fields or not fields[0]:
                    continue
                if len(fields) < 2 or not fields[1]:
                    raise ValueError(
                        f"malformed gene-arm row {path}:{line_number}")
                arm = fields[1].strip()
                if len(arm) < 2 or arm[-1].lower() not in {"p", "q"}:
                    continue
                chromosome = _hybrid_chromosome(arm[:-1])
                autosomal = chromosome.isdigit() and 1 <= int(chromosome) <= 22
                sex = chromosome.upper() in {"X", "Y"}
                if not autosomal and not (include_sex_chromosomes and sex):
                    continue
                canonical_arm = chromosome + arm[-1].lower()
                result.setdefault(chromosome, set()).add(canonical_arm)
    except (OSError, EOFError, UnicodeError) as exc:
        raise ValueError(f"cannot read gene-arm map {path}: {exc}") from exc
    if not result:
        raise ValueError(f"gene-arm map contains no selected p/q arms: {path}")
    return {
        chromosome: tuple(sorted(arms, key=_hybrid_natural_key))
        for chromosome, arms in sorted(
            result.items(), key=lambda item: _hybrid_natural_key(item[0]))
    }


def _hybrid_file_record(path: str, schema: str, stage: str,
                        library: int | None = None) -> dict[str, object]:
    target = absolute(path)
    if not regular_nonempty(target):
        raise ValueError(f"source {stage} file is missing or empty: {target}")
    record: dict[str, object] = {
        "stage": stage, "path": target, "size": os.path.getsize(target),
        "schema_version": schema,
    }
    if library is not None:
        record["library"] = library
    return record


def validate_hybrid_source_run(
        source_root: str, destination_root: str, libraries: Sequence[int]
        ) -> tuple[dict[str, object], dict[int, dict[str, str]]]:
    """Validate a native run without modifying or hashing its files."""
    source = absolute(source_root)
    destination = absolute(destination_root)
    if not os.path.isdir(source):
        raise ValueError(f"--hybrid-source-run-root is not a directory: {source}")
    if (path_contains(source, destination) or
            path_contains(destination, source)):
        raise ValueError(
            "--hybrid-source-run-root and --run-root must be separate, "
            "non-nested directories")
    source_run = run_paths(source)
    files: list[dict[str, object]] = []

    reference_contract = read_json_object(
        source_run.reference_contract, "source REFERENCE")
    reference_schema = str(reference_contract.get("schema_version", ""))
    if (reference_schema not in {
            "tetra_arm_reference_contract_v1",
            "tetra_arm_gene_synteny_contract_v1",
            "tetra_arm_hal_liftover_contract_v1",
            } or str(reference_contract.get("status", "")).upper() != "PASS"):
        raise ValueError(
            f"source REFERENCE contract is incompatible: "
            f"{source_run.reference_contract}")
    arms_bed = source_run.generated_arms_bed
    if not regular_nonempty(arms_bed):
        named = reference_contract.get("output") or reference_contract.get("arms_bed")
        if named and regular_nonempty(absolute(str(named))):
            arms_bed = absolute(str(named))
        else:
            raise ValueError(
                f"source REFERENCE arm BED is unavailable: {arms_bed}")
    files.extend([
        _hybrid_file_record(
            source_run.reference_contract, reference_schema, "REFERENCE"),
        _hybrid_file_record(arms_bed, "BED4", "REFERENCE"),
    ])
    if regular_nonempty(source_run.reference_qc):
        files.append(_hybrid_file_record(
            source_run.reference_qc, "REFERENCE_QC", "REFERENCE"))

    ledger_contract_path = os.path.join(
        source_run.ledger, "split_ledger_contract.json")
    ledger_contract = require_contract(
        ledger_contract_path, "source LEDGER",
        "tetra_arm_split_ledger_contract_v1")
    try:
        observed_libraries = [int(value)
                              for value in ledger_contract.get("libraries", [])]
    except (TypeError, ValueError) as exc:
        raise ValueError(
            f"source LEDGER libraries are malformed: {ledger_contract_path}") \
            from exc
    if observed_libraries != list(libraries):
        raise ValueError(
            "source LEDGER libraries do not exactly match the resolved "
            "upstream calibration cohort: "
            f"observed={observed_libraries}, expected={list(libraries)}. "
            "--hybrid-source-run-root must cover every library in "
            "--upstream-cohort-libraries; otherwise omit the source root, "
            "provide a reference input, and select REFERENCE, LEDGER, "
            "PREPARE, and ASE so the full cohort is built in the new run")
    files.append(_hybrid_file_record(
        ledger_contract_path, "tetra_arm_split_ledger_contract_v1", "LEDGER"))

    reused: dict[int, dict[str, str]] = {}
    allowed_ase_status = {
        "PASS", "PASS_NO_HETEROTYPIC_TARGETS",
        "PASS_NO_INFORMATIVE_EVIDENCE",
    }
    for library in libraries:
        split_ledger = os.path.join(
            source_run.ledger, f"lib{library}.final_cells.tsv.gz")
        split_fields = set(read_tsv_header(
            split_ledger, f"source LEDGER lib{library}", compressed=True))
        if not {"library", "barcode", "assignment_status",
                "final_assignment"} <= split_fields:
            raise ValueError(
                f"source LEDGER header is incomplete: {split_ledger}")

        cell_manifest = os.path.join(
            source_run.prepare, f"lib{library}.cell_manifest.tsv.gz")
        cell_header = read_tsv_header(
            cell_manifest, f"source cell manifest lib{library}", compressed=True)
        if not {"library", "barcode", "donor_a", "donor_b",
                "schema_version"} <= set(cell_header):
            raise ValueError(
                f"source cell manifest header is incomplete: {cell_manifest}")
        manifest_cells: dict[str, tuple[str, str]] = {}
        try:
            with gzip.open(
                    cell_manifest, "rt", encoding="utf-8", newline="") as handle:
                for row in csv.DictReader(handle, delimiter="\t"):
                    barcode = str(row.get("barcode", "")).strip()
                    if (str(row.get("library", "")).strip() != str(library)
                            or row.get("schema_version")
                            != "tetra_arm_cell_manifest_v1"
                            or not barcode or barcode in manifest_cells):
                        raise ValueError(
                            f"invalid source cell-manifest key: "
                            f"lib{library}/{barcode}")
                    manifest_cells[barcode] = (
                        str(row.get("donor_a", "")).strip(),
                        str(row.get("donor_b", "")).strip())
        except (OSError, EOFError, csv.Error, UnicodeError) as exc:
            raise ValueError(
                f"source cell manifest is unreadable: {cell_manifest}: {exc}") \
                from exc

        ase = os.path.join(source_run.ase, f"lib{library}.arm_ase.tsv.gz")
        if tuple(read_tsv_header(
                ase, f"source ASE lib{library}", compressed=True)) != ASE_V2_HEADER:
            raise ValueError(f"source ASE header is incompatible: {ase}")
        ase_keys: set[tuple[str, str]] = set()
        try:
            with gzip.open(ase, "rt", encoding="utf-8", newline="") as handle:
                for row in csv.DictReader(handle, delimiter="\t"):
                    barcode = str(row.get("barcode", "")).strip()
                    arm = str(row.get("arm", "")).strip()
                    key = (barcode, arm)
                    if (str(row.get("library", "")).strip() != str(library)
                            or row.get("schema_version")
                            != "tetra_arm_ase_evidence_v2"
                            or not barcode or not arm or key in ase_keys
                            or barcode not in manifest_cells):
                        raise ValueError(
                            f"invalid source ASE key: lib{library}/{barcode}/{arm}")
                    if manifest_cells[barcode] != (
                            str(row.get("donor_a", "")).strip(),
                            str(row.get("donor_b", "")).strip()):
                        raise ValueError(
                            f"source ASE donor order disagrees with its cell "
                            f"manifest: lib{library}/{barcode}/{arm}")
                    ase_keys.add(key)
        except (OSError, EOFError, csv.Error, UnicodeError) as exc:
            raise ValueError(f"source ASE is unreadable: {ase}: {exc}") from exc
        ase_qc = os.path.join(source_run.ase, f"lib{library}.ase_qc.tsv")
        metrics = read_metric_map(ase_qc, f"source ASE lib{library} QC")
        if (metrics.get("schema_version") != "tetra_arm_ase_qc_v2" or
                metrics.get("status", "").upper() not in allowed_ase_status):
            raise ValueError(f"source ASE QC is incompatible: {ase_qc}")

        files.extend([
            _hybrid_file_record(
                split_ledger, "tetra_arm_split_ledger_v1", "LEDGER", library),
            _hybrid_file_record(
                cell_manifest, "tetra_arm_cell_manifest_v1", "PREPARE_CONTEXT",
                library),
            _hybrid_file_record(
                ase, "tetra_arm_ase_evidence_v2", "ASE", library),
            _hybrid_file_record(
                ase_qc, "tetra_arm_ase_qc_v2", "ASE", library),
        ])
        reused[library] = {
            "split_ledger": absolute(split_ledger),
            "source_cell_manifest": absolute(cell_manifest),
            "ase": absolute(ase), "ase_qc": absolute(ase_qc),
        }

    payload: dict[str, object] = {
        "schema_version": "tetra_arm_hybrid_source_reuse_v1",
        "orchestrator_release": RELEASE,
        "source_run_root": source,
        "destination_run_root": destination,
        "libraries": list(libraries),
        "reused_stages": ["REFERENCE", "LEDGER", "ASE"],
        "arms_bed": absolute(arms_bed),
        "files": sorted(files, key=lambda record: (
            HYBRID_STAGES.index(str(record["stage"]).replace(
                "PREPARE_CONTEXT", "PREPARE")),
            int(record.get("library", 0)), str(record["path"]))),
        "hash_policy": "NO_HASH_REUSE_GATE",
        "access_policy": "READ_ONLY_SOURCE_NO_WRITES_NO_ARCHIVES_NO_SYMLINKS",
        "status": "PASS",
    }
    return payload, reused


def hybrid_model_bundle_paths(
        run: HybridRunPaths, chromosome: str) -> tuple[str, ...]:
    prefix = os.path.join(run.model, f"tetra_arm_hybrid_{chromosome}")
    return tuple(prefix + suffix for suffix in (
        ".scores.tsv.gz", ".expression_components.tsv.gz",
        ".calibration.tsv.gz", ".qc.tsv", ".contract.json"))


def hybrid_stage_outputs(
        stage: str, run: HybridRunPaths, paths: Sequence[LibraryPaths],
        args: argparse.Namespace,
        chromosome_arms: Mapping[str, Sequence[str]]) -> list[str]:
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
            path.cell_manifest, path.hybrid_cell_manifest,
            path.ambient_sources, path.prepare_qc, path.prepare_contract)]
    if stage == "ASE":
        return [value for path in paths for value in (path.ase, path.ase_qc)]
    if stage == "EXPRESSION":
        return [value for path in paths for value in (
            path.expression, path.expression_model,
            path.expression_qc, path.expression_contract)]
    if stage == "HYBRID_SHARD":
        return [value for path in paths for value in (
            os.path.join(run.shard, f"lib{path.library}.hybrid_shard_qc.tsv"),
            os.path.join(run.shard,
                         f"lib{path.library}.hybrid_shard_contract.json"),
            *(os.path.join(
                run.shard,
                f"lib{path.library}.chr{chromosome}.hybrid_shard.tsv.gz")
              for chromosome in chromosome_arms),
        )]
    if stage == "MODEL":
        return [value for chromosome in chromosome_arms
                for value in hybrid_model_bundle_paths(run, chromosome)]
    if stage == "CALL":
        return list(run.call_outputs)
    if stage == "REPORT":
        return list(run.report_outputs)
    raise ValueError(f"unknown hybrid stage: {stage}")


def _replace_generated_fragment(
        payload: str, old: str, new: str, label: str) -> str:
    if payload.count(old) != 1:
        raise RuntimeError(
            f"cannot add hybrid {label}; legacy script template changed")
    return payload.replace(old, new, 1)


def hybrid_prepare_script(
        args: argparse.Namespace, run: HybridRunPaths,
        manifest: str, tasks: int) -> str:
    """Extend the existing PREPARE job without changing its legacy path."""
    payload = prepare_script(args, run, manifest, tasks)
    payload = _replace_generated_fragment(
        payload,
        'cell_manifest="$output_dir/lib${library}.cell_manifest.tsv.gz"\n'
        'ambient_sources="$output_dir/lib${library}.ambient_sources.tsv.gz"',
        'cell_manifest="$output_dir/lib${library}.cell_manifest.tsv.gz"\n'
        'hybrid_cell_manifest="$output_dir/lib${library}.hybrid_cell_manifest.tsv.gz"\n'
        'ambient_sources="$output_dir/lib${library}.ambient_sources.tsv.gz"',
        "manifest output declaration")
    payload = _replace_generated_fragment(
        payload,
        'for target in "$cell_manifest" "$ambient_sources" "$prepare_qc" "$prepare_contract"; do',
        'for target in "$cell_manifest" "$hybrid_cell_manifest" "$ambient_sources" "$prepare_qc" "$prepare_contract"; do',
        "output conflict guard")
    payload = _replace_generated_fragment(
        payload,
        '    --output-dir "$output_dir"\n)',
        '    --output-dir "$output_dir"\n'
        '    --emit-hybrid-manifest\n)',
        "producer option")
    hybrid_validation = r'''
hybrid_manifest_valid() {
    python3 - "$hybrid_cell_manifest" "$prepare_contract" "$library" <<'PY'
import csv
import gzip
import json
import os
import sys

manifest, contract_path, library = sys.argv[1:]
required = {
    "library", "barcode", "donor_a", "donor_b", "donor_pair",
    "hybrid_target_eligible", "hybrid_target_reasons",
    "expression_reference_eligible", "expression_reference_reasons",
    "cell_group_source", "cell_group_target_chromosome_excluded",
    "schema_version",
}
try:
    if not os.path.isfile(manifest) or os.path.getsize(manifest) <= 0:
        raise ValueError(manifest)
    with gzip.open(manifest, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = list(reader.fieldnames or [])
        if len(fields) != len(set(fields)) or not required <= set(fields):
            raise ValueError(manifest)
        seen = set()
        for row in reader:
            key = row.get("barcode", "")
            if (None in row or any(value is None for value in row.values())
                    or row.get("library") != str(library) or not key
                    or key in seen or row.get("schema_version")
                    != "tetra_arm_hybrid_cell_manifest_v1"):
                raise ValueError(manifest)
            seen.add(key)
    with open(contract_path, "r", encoding="utf-8") as handle:
        contract = json.load(handle)
    hybrid = contract.get("hybrid", {})
    if (hybrid.get("schema_version")
            != "tetra_arm_hybrid_cell_manifest_v1"
            or os.path.abspath(hybrid.get("hybrid_cell_manifest", ""))
               != os.path.abspath(manifest)):
        raise ValueError(contract_path)
except (OSError, EOFError, ValueError, TypeError, csv.Error):
    raise SystemExit(1)
PY
}

hybrid_manifest_valid || {
    echo "ERROR: PREPARE did not produce a valid hybrid cell manifest" >&2
    exit 1
}

'''
    payload = _replace_generated_fragment(
        payload, "prepare_bundle_valid || {", hybrid_validation +
        "prepare_bundle_valid || {", "hybrid validation")
    return payload


def hybrid_expression_script(
        args: argparse.Namespace, run: HybridRunPaths,
        manifest: str, tasks: int) -> str:
    """Extend legacy expression aggregation with independent model evidence."""
    payload = expression_script(args, run, manifest, tasks)
    payload = _replace_generated_fragment(
        payload,
        'expression="$output_dir/lib${library}.arm_expression.tsv.gz"\n'
        'expression_qc="$output_dir/lib${library}.expression_qc.tsv"',
        'expression="$output_dir/lib${library}.arm_expression.tsv.gz"\n'
        'expression_model="$output_dir/lib${library}.arm_expression_model.tsv.gz"\n'
        f'hybrid_cell_manifest={shlex.quote(run.prepare)}/lib${{library}}.hybrid_cell_manifest.tsv.gz\n'
        'expression_qc="$output_dir/lib${library}.expression_qc.tsv"',
        "model output declaration")
    payload = _replace_generated_fragment(
        payload,
        'for target in "$expression" "$expression_qc" "$expression_contract"; do',
        'for target in "$expression" "$expression_model" "$expression_qc" "$expression_contract"; do',
        "output conflict guard")
    payload = _replace_generated_fragment(
        payload,
        f'    --pseudocount {args.expression_pseudocount:.17g}\n)',
        f'    --pseudocount {args.expression_pseudocount:.17g}\n'
        '    --emit-model-evidence\n'
        '    --hybrid-cell-manifest "$hybrid_cell_manifest"\n'
        f'    --expression-scale-factor {args.hybrid_expression_scale_factor:.17g}\n'
        f'    --gene-folds {args.hybrid_gene_folds}\n)',
        "producer options")
    model_validation = r'''
expression_model_valid() {
    python3 - "$expression_model" "$expression_contract" "$library" <<'PY'
import csv
import gzip
import json
import os
import sys

evidence, contract_path, library = sys.argv[1:]
required = {
    "library", "barcode", "chromosome", "p_arm", "q_arm", "p_score",
    "q_score", "hybrid_target_eligible", "expression_reference_eligible",
    "schema_version",
}
try:
    if not os.path.isfile(evidence) or os.path.getsize(evidence) <= 0:
        raise ValueError(evidence)
    with gzip.open(evidence, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = list(reader.fieldnames or [])
        if len(fields) != len(set(fields)) or not required <= set(fields):
            raise ValueError(evidence)
        seen = set()
        for row in reader:
            key = (row.get("barcode", ""), row.get("chromosome", ""))
            if (None in row or any(value is None for value in row.values())
                    or row.get("library") != str(library) or not all(key)
                    or key in seen or row.get("schema_version")
                    != "tetra_arm_expression_model_evidence_v1"):
                raise ValueError(evidence)
            seen.add(key)
    with open(contract_path, "r", encoding="utf-8") as handle:
        contract = json.load(handle)
    model = contract.get("model_evidence", {})
    if (model.get("schema_version")
            != "tetra_arm_expression_model_evidence_v1"
            or os.path.abspath(model.get("path", ""))
               != os.path.abspath(evidence)):
        raise ValueError(contract_path)
except (OSError, EOFError, ValueError, TypeError, csv.Error):
    raise SystemExit(1)
PY
}

expression_model_valid || {
    echo "ERROR: EXPRESSION did not produce valid hybrid model evidence" >&2
    exit 1
}

'''
    payload = _replace_generated_fragment(
        payload, "expression_bundle_valid || {", model_validation +
        "expression_bundle_valid || {", "model validation")
    return payload


def hybrid_shard_script(
        args: argparse.Namespace, run: HybridRunPaths,
        manifest: str, tasks: int,
        chromosome_arms: Mapping[str, Sequence[str]]) -> str:
    payload = sbatch_header(
        "HYBRID_SHARD", run, args, args.hybrid_shard_cpus,
        args.hybrid_shard_memory, SCIENTIFIC_PYTHON_MODULES,
        ("numpy", "scipy"), tasks)
    payload += "\ncommand -v python3 >/dev/null 2>&1\ncommand -v awk >/dev/null 2>&1\n"
    payload += task_row_block(manifest, HYBRID_SHARD_TASK_HEADER)
    chromosomes = " ".join(shlex.quote(value) for value in chromosome_arms)
    checks = "\n".join(
        f'gzip -t "$output_directory/lib${{library}}.chr{chromosome}.hybrid_shard.tsv.gz"'
        for chromosome in chromosome_arms)
    return payload + f"""
HYBRID_SHARD_SCRIPT={shlex.quote(args.hybrid_shard_script)}
test -s "$HYBRID_SHARD_SCRIPT"
python3 "$HYBRID_SHARD_SCRIPT" --version
python3 "$HYBRID_SHARD_SCRIPT" \
    --input-manifest "$input_manifest" \
    --library "$library" \
    --output-dir "$output_directory" \
    --chromosomes {chromosomes} \
    --crossfit-folds {args.hybrid_crossfit_folds}
test -s "$qc"
test -s "$contract"
gzip -t "$output_directory/lib${{library}}.ase_calibration_payload.tsv.gz"
{checks}
echo "COMPLETE: HYBRID_SHARD lib${{library}}"
date
"""


def hybrid_model_script(
        args: argparse.Namespace, run: HybridRunPaths,
        manifest: str, tasks: int) -> str:
    payload = sbatch_header(
        "MODEL", run, args, args.hybrid_model_cpus,
        args.hybrid_model_memory, SCIENTIFIC_PYTHON_MODULES,
        ("numpy", "scipy"), tasks)
    payload += "\ncommand -v python3 >/dev/null 2>&1\ncommand -v awk >/dev/null 2>&1\n"
    payload += task_row_block(manifest, HYBRID_MODEL_TASK_HEADER)
    return payload + f"""
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
HYBRID_MODEL_SCRIPT={shlex.quote(args.hybrid_model_script)}
test -s "$HYBRID_MODEL_SCRIPT"
python3 "$HYBRID_MODEL_SCRIPT" --version
python3 "$HYBRID_MODEL_SCRIPT" \
    --chromosome "$chromosome" \
    --shard-manifest "$shard_manifest" \
    --output-dir "$output_directory" \
    --target-libraries {shlex.quote(','.join(str(value) for value in args.hybrid_output_libraries_resolved))} \
    --min-expression-reference-cells {args.hybrid_min_expression_reference_cells} \
    --min-expression-nonzero-genes {args.hybrid_min_expression_nonzero_genes} \
    --min-ase-calibration-cells {args.hybrid_min_ase_calibration_cells} \
    --min-ase-null-cells {args.hybrid_min_ase_null_cells}
test -s "$qc"
test -s "$contract"
prefix="$output_directory/tetra_arm_hybrid_${{chromosome}}"
gzip -t "${{prefix}}.scores.tsv.gz"
gzip -t "${{prefix}}.expression_components.tsv.gz"
gzip -t "${{prefix}}.calibration.tsv.gz"
echo "COMPLETE: MODEL chromosome $chromosome ($logical_arms)"
date
"""


def hybrid_call_script(
        args: argparse.Namespace, run: HybridRunPaths) -> str:
    payload = sbatch_header(
        "CALL", run, args, args.hybrid_call_cpus,
        args.hybrid_call_memory, SCIENTIFIC_PYTHON_MODULES,
        ("numpy", "scipy"))
    return payload + f"""
command -v python3 >/dev/null 2>&1
HYBRID_AGGREGATE_SCRIPT={shlex.quote(args.hybrid_aggregate_script)}
MODEL_TASK_MANIFEST={shlex.quote(run.model_tasks)}
OUTPUT_PREFIX={shlex.quote(run.call_prefix)}
test -s "$HYBRID_AGGREGATE_SCRIPT"
test -s "$MODEL_TASK_MANIFEST"
python3 "$HYBRID_AGGREGATE_SCRIPT" --version
python3 "$HYBRID_AGGREGATE_SCRIPT" \
    --model-task-manifest "$MODEL_TASK_MANIFEST" \
    --output-prefix "$OUTPUT_PREFIX" \
    --max-q {args.hybrid_max_q:.17g}
gzip -t "${{OUTPUT_PREFIX}}.arm_calls.tsv.gz"
gzip -t "${{OUTPUT_PREFIX}}.calibration.tsv.gz"
gzip -t "${{OUTPUT_PREFIX}}.expression_components.tsv.gz"
gzip -t "${{OUTPUT_PREFIX}}.uid_chromosome_flags.tsv.gz"
gzip -t "${{OUTPUT_PREFIX}}.donor_pair_arm_summary.tsv.gz"
test -s "${{OUTPUT_PREFIX}}.qc.tsv"
test -s "${{OUTPUT_PREFIX}}.contract.json"
echo "COMPLETE: CALL $OUTPUT_PREFIX"
date
"""


def hybrid_report_script(
        args: argparse.Namespace, run: HybridRunPaths,
        qc_files: Sequence[str]) -> str:
    payload = sbatch_header(
        "REPORT", run, args, args.hybrid_report_cpus,
        args.hybrid_report_memory, SCIENTIFIC_PYTHON_MODULES,
        ("csv", "gzip", "json", "numpy"))
    qc_words = "\n".join(f"    {shlex.quote(path)}" for path in qc_files)
    return payload + f"""
command -v python3 >/dev/null 2>&1
REPORT_SCRIPT={shlex.quote(args.report_script)}
OUTPUT_PREFIX={shlex.quote(run.call_prefix)}
REPORT_DIR={shlex.quote(run.report)}
QC_FILES=(
{qc_words}
)
test -s "$REPORT_SCRIPT"
python3 "$REPORT_SCRIPT" --version
for input in \
    "${{OUTPUT_PREFIX}}.arm_calls.tsv.gz" \
    "${{OUTPUT_PREFIX}}.calibration.tsv.gz" \
    "${{OUTPUT_PREFIX}}.uid_chromosome_flags.tsv.gz" \
    "${{OUTPUT_PREFIX}}.donor_pair_arm_summary.tsv.gz" \
    "${{OUTPUT_PREFIX}}.qc.tsv" \
    "${{OUTPUT_PREFIX}}.contract.json"; do
    test -s "$input"
done
for input in "${{QC_FILES[@]}}"; do test -s "$input"; done
python3 "$REPORT_SCRIPT" \
    --calls "${{OUTPUT_PREFIX}}.arm_calls.tsv.gz" \
    --calibration "${{OUTPUT_PREFIX}}.calibration.tsv.gz" \
    --uid-summary "${{OUTPUT_PREFIX}}.uid_chromosome_flags.tsv.gz" \
    --donor-pair-summary "${{OUTPUT_PREFIX}}.donor_pair_arm_summary.tsv.gz" \
    --call-contract "${{OUTPUT_PREFIX}}.contract.json" \
    --output-dir "$REPORT_DIR" \
    --title {shlex.quote(args.report_title)} \
    --top-calls {args.report_top_calls} \
    --qc-files "${{QC_FILES[@]}}"
test -s "$REPORT_DIR/tetra_arm_cnv_hybrid_summary.tsv"
test -s "$REPORT_DIR/tetra_arm_cnv_hybrid_summary.json"
test -s "$REPORT_DIR/tetra_arm_cnv_hybrid_report.html"
echo "COMPLETE: REPORT $REPORT_DIR/tetra_arm_cnv_hybrid_report.html"
date
"""


def validate_hybrid_args(
        args: argparse.Namespace) -> tuple[list[int], tuple[str, ...]]:
    libraries = parse_libraries(args.libraries)
    upstream_cohort = (
        parse_libraries(args.upstream_cohort_libraries)
        if args.upstream_cohort_libraries else list(libraries))
    if not set(libraries) <= set(upstream_cohort):
        raise ValueError(
            "--upstream-cohort-libraries must contain every output library")
    args.upstream_cohort_libraries_resolved = upstream_cohort
    args.hybrid_output_libraries_resolved = list(libraries)

    source_supplied = bool(args.hybrid_source_run_root)
    selected = parse_hybrid_stages(args.stages)
    if source_supplied:
        if args.stages and set(selected) & {"REFERENCE", "LEDGER", "ASE"}:
            raise ValueError(
                "a hybrid source run already supplies REFERENCE, LEDGER, and "
                "ASE; do not select those stages")
        if not args.stages:
            selected = tuple(stage for stage in HYBRID_STAGES
                             if stage not in {"REFERENCE", "LEDGER", "ASE"})
    if args.resume and not args.stages:
        raise ValueError(
            "--resume requires an explicit --stage selection for hybrid runs")
    if args.array_throttle is not None and args.array_throttle < 1:
        raise ValueError("--array-throttle must be positive")
    for name in (
            "hybrid_shard_cpus", "hybrid_model_cpus", "hybrid_call_cpus",
            "hybrid_report_cpus", "ledger_cpus", "ase_threads",
            "expression_cpus"):
        value = int(getattr(args, name))
        if not 1 <= value <= 256:
            raise ValueError(
                f"--{name.replace('_', '-')} must be between 1 and 256")
    for name in (
            "hybrid_shard_memory", "hybrid_model_memory",
            "hybrid_call_memory", "hybrid_report_memory",
            "reference_memory", "ledger_memory", "prepare_memory",
            "ase_memory", "expression_memory"):
        if not re.fullmatch(r"[1-9][0-9]*(?:[KMGTP])?", getattr(args, name)):
            raise ValueError(
                f"--{name.replace('_', '-')} must be a positive SLURM memory token")
    if not re.fullmatch(
            r"[0-9]+-[0-9]{2}:[0-9]{2}:[0-9]{2}|"
            r"[0-9]{1,3}:[0-9]{2}:[0-9]{2}", args.time):
        raise ValueError("--time must use D-HH:MM:SS or HH:MM:SS")
    if not re.fullmatch(r"[A-Za-z0-9_.-]+", args.partition):
        raise ValueError("--partition contains unsupported characters")
    if args.hybrid_crossfit_folds < 2:
        raise ValueError("--hybrid-crossfit-folds must be at least two")
    if args.hybrid_gene_folds != 2:
        raise ValueError("hybrid v1 requires exactly --hybrid-gene-folds 2")
    if (not math.isfinite(args.hybrid_expression_scale_factor) or
            args.hybrid_expression_scale_factor <= 0.0):
        raise ValueError(
            "--hybrid-expression-scale-factor must be finite and positive")
    for name in (
            "hybrid_min_expression_reference_cells",
            "hybrid_min_expression_nonzero_genes",
            "hybrid_min_ase_calibration_cells", "hybrid_min_ase_null_cells"):
        if getattr(args, name) < 1:
            raise ValueError(f"--{name.replace('_', '-')} must be positive")
    if not 0.0 < args.hybrid_max_q <= 1.0:
        raise ValueError("--hybrid-max-q must be in (0,1]")
    if not 1 <= args.report_top_calls <= 10000:
        raise ValueError("--report-top-calls must be between 1 and 10000")

    run_root_was_explicit = bool(args.run_root)
    configure_input_roots(args)
    if not run_root_was_explicit:
        args.run_root = absolute(os.path.join(
            args.upstream_analysis_root, "aggregate_library_analysis",
            "tetra_arm_cnv_hybrid"))
    default_templates(args)
    run = hybrid_run_paths(args.run_root)
    if source_supplied:
        if args.arms_bed or args.gene_annotation or args.hal_file:
            raise ValueError(
                "--hybrid-source-run-root supplies the reference; do not also "
                "set --arms-bed, --gene-annotation, or --hal-file")
        args.reference_mode = "SOURCE_REUSE"
        args.hybrid_source_run_root = absolute(args.hybrid_source_run_root)
        args.arms_bed = os.path.join(
            args.hybrid_source_run_root, "reference", "ancestral_arms.bed")
        args.gene_annotation = ""
        args.reference_fai = ""
        args.hal_file = ""
    else:
        if not (args.arms_bed or args.gene_annotation or args.hal_file):
            raise ValueError(
                "the hybrid workflow requires one reference input or "
                "--hybrid-source-run-root")
        if args.gene_annotation:
            args.reference_mode = "GENE_SYNTENY"
            args.gene_annotation = absolute(args.gene_annotation)
            args.hal_file = ""
            args.reference_fai = (
                absolute(args.reference_fai) if args.reference_fai else "")
            if (args.gene_projection_mode == "anchored-contigs" and
                    not args.reference_fai):
                raise ValueError(
                    "--reference-fai is required with anchored-contigs")
            args.arms_bed = run.generated_arms_bed
        elif args.hal_file:
            args.reference_mode = "HAL_LIFTOVER"
            args.gene_annotation = ""
            args.reference_fai = ""
            args.hal_file = absolute(args.hal_file)
            args.arms_bed = run.generated_arms_bed
        else:
            args.reference_mode = "EXPLICIT_BED"
            args.gene_annotation = ""
            args.reference_fai = ""
            args.hal_file = ""
            args.arms_bed = absolute(args.arms_bed)

    args.run_root = absolute(args.run_root)
    args.gene_arms = absolute(args.gene_arms)
    args.panel_metadata = absolute(args.panel_metadata) if args.panel_metadata else ""
    for name in (
            "prepare_script", "expression_script", "ase_binary",
            "report_script", "arm_builder_script", "hal_reference_script",
            "source_arms_bed", "interindividual_panel",
            "hybrid_shard_script", "hybrid_model_script",
            "hybrid_aggregate_script"):
        setattr(args, name, absolute(getattr(args, name)))
    for field in (
            "run_root", "hybrid_source_run_root", "arms_bed", "gene_arms",
            "gene_annotation", "reference_fai", "hal_file", "prepare_script",
            "expression_script", "ase_binary", "report_script",
            "hybrid_shard_script", "hybrid_model_script",
            "hybrid_aggregate_script"):
        validate_text(str(getattr(args, field)),
                      f"--{field.replace('_', '-')}")
    if any(character.isspace() for character in args.run_root):
        raise ValueError("--run-root cannot contain whitespace in SLURM log paths")
    return libraries, selected


def make_hybrid_library_paths(
        args: argparse.Namespace, run: HybridRunPaths,
        libraries: Sequence[int], selected: Sequence[str],
        reused: Mapping[int, Mapping[str, str]]) -> list[LibraryPaths]:
    paths = make_library_paths(args, run, libraries, selected)
    if not reused:
        return paths
    return [replace(
        path,
        split_ledger=str(reused[path.library]["split_ledger"]),
        ase=str(reused[path.library]["ase"]),
        ase_qc=str(reused[path.library]["ase_qc"]),
    ) for path in paths]


def write_hybrid_manifests(
        args: argparse.Namespace, run: HybridRunPaths,
        paths: Sequence[LibraryPaths], selected: Sequence[str],
        chromosome_arms: Mapping[str, Sequence[str]]) -> dict[str, str]:
    manifests = prepare_manifests(args, run, paths, selected)
    write_task_manifest(run.hybrid_inputs, HYBRID_INPUT_HEADER, ({
        "library": path.library,
        "hybrid_cell_manifest": path.hybrid_cell_manifest,
        "ase": path.ase,
        "expression": path.expression,
        "expression_model": path.expression_model,
    } for path in paths))
    write_task_manifest(run.shard_tasks, HYBRID_SHARD_TASK_HEADER, ({
        "task_index": index,
        "library": path.library,
        "input_manifest": run.hybrid_inputs,
        "output_directory": run.shard,
        "qc": os.path.join(
            run.shard, f"lib{path.library}.hybrid_shard_qc.tsv"),
        "contract": os.path.join(
            run.shard, f"lib{path.library}.hybrid_shard_contract.json"),
    } for index, path in enumerate(paths)))
    model_input_root = os.path.join(run.manifests, "hybrid_model_inputs")
    model_rows = []
    for index, (chromosome, arms) in enumerate(chromosome_arms.items()):
        shard_manifest = os.path.join(
            model_input_root, f"chr{chromosome}.tsv")
        write_task_manifest(
            shard_manifest, ("library", "shard", "qc", "contract"), ({
                "library": path.library,
                "shard": os.path.join(
                    run.shard,
                    f"lib{path.library}.chr{chromosome}.hybrid_shard.tsv.gz"),
                "qc": os.path.join(
                    run.shard, f"lib{path.library}.hybrid_shard_qc.tsv"),
                "contract": os.path.join(
                    run.shard,
                    f"lib{path.library}.hybrid_shard_contract.json"),
            } for path in paths))
        model_bundle = hybrid_model_bundle_paths(run, chromosome)
        model_rows.append({
            "task_index": index, "chromosome": chromosome,
            "logical_arms": ",".join(arms),
            "shard_manifest": shard_manifest,
            "output_directory": run.model,
            "qc": model_bundle[3], "contract": model_bundle[4],
        })
    write_task_manifest(run.model_tasks, HYBRID_MODEL_TASK_HEADER, model_rows)
    manifests.update({
        "HYBRID_INPUTS": run.hybrid_inputs,
        "HYBRID_SHARD": run.shard_tasks,
        "MODEL": run.model_tasks,
    })
    return manifests


def hybrid_omitted_ancestors(selected: Sequence[str]) -> tuple[str, ...]:
    selected_set = set(selected)
    ancestors: set[str] = set()

    def visit(stage: str) -> None:
        for parent in HYBRID_STAGE_PARENTS[stage]:
            if parent not in selected_set:
                ancestors.add(parent)
            visit(parent)

    for stage in selected:
        visit(stage)
    return tuple(stage for stage in HYBRID_STAGES if stage in ancestors)


def validate_completed_hybrid_stage(
        stage: str, run: HybridRunPaths, paths: Sequence[LibraryPaths],
        args: argparse.Namespace,
        chromosome_arms: Mapping[str, Sequence[str]]) -> None:
    outputs = hybrid_stage_outputs(stage, run, paths, args, chromosome_arms)
    missing = [path for path in outputs if not regular_nonempty(path)]
    if missing:
        raise ValueError(
            f"completed hybrid {stage} bundle is missing or empty: "
            + ", ".join(missing[:3]))
    if stage in {"REFERENCE", "LEDGER", "ASE"}:
        validate_completed_stage(stage, run, paths, args)
        return
    if stage == "PREPARE":
        for path in paths:
            fields = set(read_tsv_header(
                path.hybrid_cell_manifest,
                f"hybrid PREPARE lib{path.library}", compressed=True))
            if not {"library", "barcode", "hybrid_target_eligible",
                    "expression_reference_eligible", "schema_version"} <= fields:
                raise ValueError(
                    f"hybrid manifest is incomplete: {path.hybrid_cell_manifest}")
            contract = require_contract(
                path.prepare_contract, f"hybrid PREPARE lib{path.library}",
                "tetra_arm_prepare_contract_v1")
            hybrid = contract.get("hybrid", {})
            if (not isinstance(hybrid, dict) or
                    hybrid.get("schema_version")
                    != "tetra_arm_hybrid_cell_manifest_v1"):
                raise ValueError(
                    f"PREPARE contract lacks hybrid-v1 section: "
                    f"{path.prepare_contract}")
        return
    if stage == "EXPRESSION":
        for path in paths:
            fields = set(read_tsv_header(
                path.expression_model,
                f"hybrid EXPRESSION lib{path.library}", compressed=True))
            if not {"library", "barcode", "chromosome", "p_score", "q_score",
                    "schema_version"} <= fields:
                raise ValueError(
                    f"expression model header is incomplete: "
                    f"{path.expression_model}")
            contract = require_contract(
                path.expression_contract,
                f"hybrid EXPRESSION lib{path.library}",
                "tetra_arm_expression_contract_v1")
            model = contract.get("model_evidence", {})
            if (not isinstance(model, dict) or model.get("schema_version")
                    != "tetra_arm_expression_model_evidence_v1"):
                raise ValueError(
                    f"EXPRESSION contract lacks model evidence: "
                    f"{path.expression_contract}")
        return
    if stage == "HYBRID_SHARD":
        for path in paths:
            contract_path = os.path.join(
                run.shard, f"lib{path.library}.hybrid_shard_contract.json")
            contract = require_contract(
                contract_path, f"HYBRID_SHARD lib{path.library}",
                "tetra_arm_hybrid_shard_contract_v1")
            if int(contract.get("library", -1)) != path.library:
                raise ValueError(f"wrong library in {contract_path}")
            outputs = contract.get("outputs", {})
            calibration_payload = (
                str(outputs.get("ase_calibration") or "").strip()
                if isinstance(outputs, dict) else "")
            if calibration_payload:
                calibration_fields = set(read_tsv_header(
                    calibration_payload,
                    f"HYBRID_SHARD ASE calibration lib{path.library}",
                    compressed=True))
                if not {"library", "barcode",
                        "ase_calibration_observations_ref",
                        "ase_calibration_observations_alt",
                        "ase_calibration_observations_mixed",
                        "schema_version"} <= calibration_fields:
                    raise ValueError(
                        f"ASE calibration payload is incomplete: "
                        f"{calibration_payload}")
            for chromosome in chromosome_arms:
                shard_path = os.path.join(
                    run.shard,
                    f"lib{path.library}.chr{chromosome}.hybrid_shard.tsv.gz")
                fields = set(read_tsv_header(
                    shard_path, f"HYBRID_SHARD {chromosome}", compressed=True))
                if not {"library", "barcode", "chromosome", "arm",
                        "schema_version"} <= fields:
                    raise ValueError(f"hybrid shard header is incomplete: {shard_path}")
        return
    if stage == "MODEL":
        for chromosome in chromosome_arms:
            bundle = hybrid_model_bundle_paths(run, chromosome)
            contract = require_contract(
                bundle[4], f"MODEL chromosome {chromosome}",
                "tetra_arm_hybrid_model_contract_v1")
            if _hybrid_chromosome(contract.get("chromosome")) != chromosome:
                raise ValueError(f"wrong chromosome in {bundle[4]}")
            for data_path in bundle[:3]:
                read_tsv_header(
                    data_path, f"MODEL chromosome {chromosome}", compressed=True)
        return
    if stage == "CALL":
        contract = require_contract(
            run.call_prefix + ".contract.json", "hybrid CALL",
            "tetra_arm_call_contract_v3")
        schemas = contract.get("output_schemas", {})
        if (not isinstance(schemas, dict) or
                schemas.get("calls") != "tetra_arm_cnv_calls_v3"):
            raise ValueError("hybrid CALL output schemas are incompatible")
        return
    if stage == "REPORT":
        report_json = read_json_object(
            os.path.join(run.report, "tetra_arm_cnv_hybrid_summary.json"),
            "hybrid REPORT")
        if report_json.get("schema_version") != "tetra_arm_cnv_report_v4":
            raise ValueError("hybrid REPORT schema is incompatible")
        return
    raise ValueError(f"unsupported completed hybrid stage: {stage}")


def hybrid_resolved_dependencies(
        stage: str, selected: Sequence[str],
        job_ids: Mapping[str, str]) -> list[str]:
    selected_set = set(selected)
    frontier: set[str] = set()

    def visit(node: str) -> None:
        for parent in HYBRID_STAGE_PARENTS[node]:
            if parent in selected_set and parent in job_ids:
                frontier.add(parent)
            else:
                visit(parent)

    visit(stage)
    return [job_ids[parent] for parent in HYBRID_STAGES if parent in frontier]


def hybrid_attempt_token(selected: Sequence[str]) -> str:
    return "_".join(stage.lower() for stage in selected)


def hybrid_generated_script_path(
        run: HybridRunPaths, filename: str, resume: bool) -> str:
    if not resume:
        return os.path.join(run.scripts, filename)
    stem, extension = os.path.splitext(filename)
    return os.path.join(run.scripts, stem + ".resume" + extension)


def hybrid_plan_path(
        run: HybridRunPaths, selected: Sequence[str], resume: bool) -> str:
    if not resume:
        return os.path.join(run.manifests, "hybrid_execution_plan.json")
    return os.path.join(
        run.manifests,
        f"hybrid_execution_plan.resume.{hybrid_attempt_token(selected)}.json")


def hybrid_submission_path(
        run: HybridRunPaths, selected: Sequence[str], resume: bool) -> str:
    if not resume:
        return os.path.join(run.manifests, "hybrid_submission.json")
    return os.path.join(
        run.manifests,
        f"hybrid_submission.resume.{hybrid_attempt_token(selected)}.json")


def replace_mutable_json_atomic(path: str, payload: Mapping[str, object]) -> str:
    """Atomically checkpoint mutable scheduler state after every accepted job."""
    destination = absolute(path)
    os.makedirs(os.path.dirname(destination), exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(
        prefix=os.path.basename(destination) + ".tmp.",
        dir=os.path.dirname(destination), text=True)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8", newline="") as handle:
            json.dump(payload, handle, sort_keys=True, indent=2)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, destination)
        directory_descriptor = os.open(os.path.dirname(destination), os.O_RDONLY)
        try:
            os.fsync(directory_descriptor)
        finally:
            os.close(directory_descriptor)
        return destination
    finally:
        try:
            os.unlink(temporary)
        except FileNotFoundError:
            pass


def build_parser() -> argparse.ArgumentParser:
    """Build the dedicated hybrid-only command line."""
    parser = argparse.ArgumentParser(
        description=(
            "Generate and optionally submit the dedicated CellBouncer "
            "hybrid-v1 expression+ASE chromosome-arm CNV DAG."),
        epilog=(
            "Planning is the default. Inspect the execution plan and rendered "
            "sbatch files, then repeat the identical command with --submit."))
    parser.add_argument(
        "--version", action="version", version=f"%(prog)s {RELEASE}")
    parser.add_argument(
        "--submit", action="store_true",
        help="Submit the reviewed plan; otherwise render it without submission")
    parser.add_argument(
        "--resume", action="store_true",
        help=(
            "Plan an explicit recovery stage selection after validating all "
            "omitted ancestor output bundles"))
    parser.add_argument(
        "--stage", "--stages", dest="stages", action="append",
        help=(
            "Selected stage(s), comma-separated or repeated; default ALL "
            f"({','.join(HYBRID_STAGES)})"))
    parser.add_argument(
        "--libraries", nargs="+", default=["1-40"],
        help="Output libraries as N, libN, N-M, or comma-separated values")
    parser.add_argument(
        "--upstream-cohort-libraries", nargs="+", default=None,
        help=(
            "Libraries used for upstream preparation, ASE, expression, "
            "sharding, calibration, and empirical nulls; must contain every "
            "output library"))
    parser.add_argument(
        "--array-throttle", "--max-concurrent", dest="array_throttle",
        type=int, default=None,
        help="Optional maximum concurrent tasks per array; omitted by default")
    parser.add_argument(
        "--nodelist", default="",
        help="Optional SLURM node list; scheduler placement is used by default")
    parser.add_argument(
        "--run-root", default=None,
        help=(
            "Hybrid output root; defaults to <upstream-analysis-root>/"
            "aggregate_library_analysis/tetra_arm_cnv_hybrid"))
    source = parser.add_argument_group("Hybrid source and upstream inputs")
    source.add_argument(
        "--hybrid-source-run-root", default="",
        help=(
            "Read-only completed native arm-CNV run supplying validated "
            "REFERENCE, LEDGER, and ASE bundles for the complete upstream "
            "cohort"))
    source.add_argument(
        "--mapping-input-root", "--mapping-root",
        dest="mapping_input_root", default=MAPPING_ROOT,
        help="Mapping-output root containing per-library filtered MEX files")
    source.add_argument(
        "--upstream-analysis-root", default=PRODUCTION_ANALYSIS_ROOT,
        help=(
            "Analysis root containing per-library DEMUX/ambient products and "
            "aggregate identity, ploidy, and GEX-calibration products"))
    source.add_argument("--mapping-run-root", default="")
    source.add_argument(
        "--identity-root", default=None,
        help=(
            "Identity-reconciliation root; defaults below the upstream "
            "analysis root"))
    source.add_argument("--identity-validation", default=None)
    source.add_argument("--identity-metadata-manifest", default=None)
    source.add_argument("--ledger-input", default=None)
    source.add_argument("--final-assignments-template", default=None)
    source.add_argument("--mapping-bam-template", default=None)
    source.add_argument("--demux-prefix-template", default=None)
    source.add_argument("--expression-barcodes-template", default=None)
    source.add_argument("--expression-features-template", default=None)
    source.add_argument("--expression-matrix-template", default=None)
    source.add_argument("--ambient-standard-template", default=None)
    source.add_argument("--ambient-arm-a-template", default=None)
    source.add_argument("--ambient-arm-c-template", default=None)
    source.add_argument(
        "--interindividual-panel", default=DEFAULT_INTERINDIVIDUAL_PANEL,
        help="Donor-genotype BCF represented by the ASE pileup sidecars")
    source.add_argument(
        "--panel-distinguishability-binary",
        default=os.path.join(
            DEPLOYED_BIN, "nuclear_panel_distinguishability"))
    source.add_argument(
        "--skip-identity-validation", action="store_true",
        help="Audit override for the upstream identity boundary")
    source.add_argument("--ploidy-nn-template", default=None)
    source.add_argument("--ploidy-input-h5ad", default="")
    source.add_argument(
        "--ploidy-nn-weights", default=DEFAULT_PLOIDY_NN_WEIGHTS)
    source.add_argument(
        "--ploidy-nn-helper", default=DEFAULT_PLOIDY_NN_HELPER)
    source.add_argument(
        "--ploidy-nn-module", default="ploidy-inference/latest")
    source.add_argument("--ambient-condition", default=DEFAULT_CONDITION)
    source.add_argument(
        "--ambient-candidate-set", default="applied",
        choices=("applied", "exploratory"))
    source.add_argument("--cell-groups-template", default="")
    source.add_argument(
        "--gex-ambient-analysis", default=DEFAULT_GEX_AMBIENT_ANALYSIS)
    source.add_argument("--panel-metadata", default=DEFAULT_PANEL_METADATA)

    reference = parser.add_argument_group("Reference construction")
    reference_input = reference.add_mutually_exclusive_group(required=False)
    reference_input.add_argument(
        "--arms-bed", default="",
        help="Existing BED4 in the same coordinates as demux pileup sites")
    reference_input.add_argument(
        "--gene-annotation", default="",
        help="Ancestral-coordinate GTF/GFF used to build the arm BED")
    reference_input.add_argument(
        "--hal-file", default="",
        help="HAL used to lift the source arm BED to ancestral coordinates")
    reference.add_argument("--gene-arms", default=DEFAULT_GENE_ARMS)
    reference.add_argument("--reference-fai", default="")
    reference.add_argument(
        "--gene-projection-mode",
        choices=("anchored-contigs", "gene-spans"),
        default="anchored-contigs")
    reference.add_argument(
        "--arm-builder-script", default=DEFAULT_ARM_BUILDER_SCRIPT)
    reference.add_argument("--source-arms-bed", default=DEFAULT_SOURCE_ARMS)
    reference.add_argument("--hal-source-genome", default="Human")
    reference.add_argument(
        "--hal-target-genome", default="human_chimp_bonobo")
    reference.add_argument(
        "--hal-reference-script", default=DEFAULT_HAL_REFERENCE_SCRIPT)

    workers = parser.add_argument_group("Scientific workers")
    workers.add_argument("--prepare-script", default=DEFAULT_PREPARE_SCRIPT)
    workers.add_argument(
        "--expression-script", default=DEFAULT_EXPRESSION_SCRIPT)
    workers.add_argument("--ase-binary", default=DEFAULT_ASE_BINARY)
    workers.add_argument("--report-script", default=DEFAULT_REPORT_SCRIPT)
    workers.add_argument(
        "--hybrid-shard-script", default=DEFAULT_HYBRID_SHARD_SCRIPT)
    workers.add_argument(
        "--hybrid-model-script", default=DEFAULT_HYBRID_MODEL_SCRIPT)
    workers.add_argument(
        "--hybrid-aggregate-script", default=DEFAULT_HYBRID_AGGREGATE_SCRIPT)

    model = parser.add_argument_group("Hybrid model")
    model.add_argument("--hybrid-crossfit-folds", type=int, default=5)
    model.add_argument(
        "--hybrid-expression-scale-factor", type=float, default=10000.0)
    model.add_argument("--hybrid-gene-folds", type=int, default=2)
    model.add_argument(
        "--hybrid-min-expression-reference-cells", type=int, default=20)
    model.add_argument(
        "--hybrid-min-expression-nonzero-genes", type=int, default=10)
    model.add_argument(
        "--hybrid-min-ase-calibration-cells", type=int, default=12)
    model.add_argument("--hybrid-min-ase-null-cells", type=int, default=20)
    model.add_argument("--hybrid-max-q", type=float, default=0.05)
    model.add_argument(
        "--min-calibration-tet-probability", type=float, default=0.90)
    hard_threshold = model.add_mutually_exclusive_group()
    hard_threshold.add_argument(
        "--hard-allele-threshold", dest="hard_allele_threshold",
        type=float, default=argparse.SUPPRESS)
    hard_threshold.add_argument(
        "--hard-posterior", dest="hard_allele_threshold", type=float,
        default=argparse.SUPPRESS,
        help="Deprecated alias for --hard-allele-threshold")
    parser.set_defaults(hard_allele_threshold=0.80)
    model.add_argument("--expression-pseudocount", type=float, default=0.5)
    model.add_argument("--include-sex-chromosomes", action="store_true")
    model.add_argument(
        "--allow-missing-expression-barcodes", action="store_true")
    model.add_argument(
        "--include-nonheterotypic-expression", action="store_true")
    model.add_argument(
        "--expression-ambient-corrected", action="store_true",
        help="Mark the selected expression MEX as already ambient-corrected")

    resources = parser.add_argument_group("SLURM resources")
    resources.add_argument("--partition", default=DEFAULT_PARTITION)
    resources.add_argument("--time", default=DEFAULT_TIME)
    resources.add_argument("--reference-memory", default="16G")
    resources.add_argument("--ledger-memory", default="256G")
    resources.add_argument("--prepare-memory", default="16G")
    resources.add_argument("--ase-memory", default="256G")
    resources.add_argument("--expression-memory", default="128G")
    resources.add_argument("--hybrid-shard-memory", default="16G")
    resources.add_argument("--hybrid-model-memory", default="32G")
    resources.add_argument("--hybrid-call-memory", default="32G")
    resources.add_argument("--hybrid-report-memory", default="8G")
    resources.add_argument("--ase-threads", type=int, default=16)
    resources.add_argument("--ledger-cpus", type=int, default=16)
    resources.add_argument("--expression-cpus", type=int, default=8)
    resources.add_argument("--hybrid-shard-cpus", type=int, default=4)
    resources.add_argument("--hybrid-model-cpus", type=int, default=4)
    resources.add_argument("--hybrid-call-cpus", type=int, default=8)
    resources.add_argument("--hybrid-report-cpus", type=int, default=2)

    report = parser.add_argument_group("Hybrid report")
    report.add_argument(
        "--report-title", default="Tetraploid chromosome-arm ASE/CNV report")
    report.add_argument("--report-top-calls", type=int, default=100)
    return parser


def validate_scheduler_overrides(args: argparse.Namespace) -> None:
    if not args.nodelist:
        return
    if (not re.fullmatch(r"[A-Za-z0-9_.\[\],-]+", args.nodelist)
            or args.nodelist.startswith((',', '-'))
            or args.nodelist.endswith(',')):
        raise ValueError(
            "--nodelist must be a comma-separated SLURM node list without "
            "whitespace or shell metacharacters")


def render_script(path: str, payload: str, args: argparse.Namespace) -> str:
    if not payload.startswith("#!/bin/bash\n"):
        raise RuntimeError("generated sbatch lacks the required bash shebang")
    required = (
        "set -eo pipefail", "module purge", "module list 2>&1", "command -v",
    )
    missing = [token for token in required if token not in payload]
    if missing:
        raise RuntimeError(
            "generated sbatch lacks required invariant(s): "
            + ", ".join(missing))
    module_loads = re.findall(r"(?m)^module load ([^\s]+)$", payload)
    if not module_loads or any("/" not in module for module in module_loads):
        raise RuntimeError(
            "generated sbatch must load at least one concrete versioned module")
    forbidden = ("conda activate", "source conda.sh", "export MODULEPATH=")
    present = [token for token in forbidden if token in payload]
    if present:
        raise RuntimeError(
            "generated sbatch contains forbidden environment setup: "
            + ", ".join(present))

    result = publish_new_or_identical(
        path, payload, "generated hybrid sbatch")
    os.chmod(result, os.stat(result).st_mode | stat.S_IXUSR)
    checked = subprocess.run(
        ["bash", "-n", result], capture_output=True, text=True, check=False)
    if checked.returncode != 0:
        detail = checked.stderr.strip() or checked.stdout.strip()
        raise RuntimeError(
            f"generated sbatch failed bash -n: {result}: {detail}")
    return result


def submission_path(
        run: HybridRunPaths, selected: Sequence[str], resume: bool
        ) -> str:
    return hybrid_submission_path(run, selected, resume)


def load_submission(
        path: str, run: HybridRunPaths, plan_path: str,
        selected: Sequence[str], scripts: Mapping[str, str],
        ) -> dict[str, object] | None:
    if not os.path.lexists(path):
        return None
    record = read_json_object(path, "hybrid submission record")
    expected = {
        "schema_version": "tetra_arm_hybrid_submission_v1",
        "orchestrator": ORCHESTRATOR_NAME,
        "orchestrator_release": RELEASE,
        "workflow": WORKFLOW,
        "run_root": run.root,
        "plan": plan_path,
        "selected_stages": list(selected),
        "scripts": dict(scripts),
    }
    for field, value in expected.items():
        if record.get(field) != value:
            raise ValueError(
                f"hybrid submission record disagrees on {field}: {path}")
    jobs = record.get("jobs")
    if not isinstance(jobs, dict):
        raise ValueError(f"hybrid submission record has invalid jobs: {path}")
    for stage, entry in jobs.items():
        if (stage not in scripts or not isinstance(entry, dict)
                or not re.fullmatch(r"\d+", str(entry.get("job_id", "")))
                or entry.get("script") != scripts[stage]
                or not isinstance(entry.get("dependencies"), list)):
            raise ValueError(
                f"hybrid submission record has invalid {stage} job: {path}")
    accepted: dict[str, str] = {}
    for stage in STAGES:
        if stage not in jobs:
            continue
        expected_dependencies = hybrid_resolved_dependencies(
            stage, selected, accepted)
        if jobs[stage]["dependencies"] != expected_dependencies:
            raise ValueError(
                f"hybrid submission record has invalid {stage} dependencies: "
                f"{path}")
        accepted[stage] = str(jobs[stage]["job_id"])
    return record


def stage_resources(args: argparse.Namespace) -> dict[str, tuple[int, str]]:
    return {
        "REFERENCE": (2, args.reference_memory),
        "LEDGER": (args.ledger_cpus, args.ledger_memory),
        "PREPARE": (2, args.prepare_memory),
        "ASE": (args.ase_threads, args.ase_memory),
        "EXPRESSION": (args.expression_cpus, args.expression_memory),
        "HYBRID_SHARD": (args.hybrid_shard_cpus, args.hybrid_shard_memory),
        "MODEL": (args.hybrid_model_cpus, args.hybrid_model_memory),
        "CALL": (args.hybrid_call_cpus, args.hybrid_call_memory),
        "REPORT": (args.hybrid_report_cpus, args.hybrid_report_memory),
    }


def run_workflow(args: argparse.Namespace) -> int:
    validate_scheduler_overrides(args)
    libraries, selected = validate_hybrid_args(args)
    calibration_libraries = list(args.upstream_cohort_libraries_resolved)
    run = hybrid_run_paths(args.run_root)
    source_payload: dict[str, object] | None = None
    reused: dict[int, dict[str, str]] = {}
    if args.hybrid_source_run_root:
        source_payload, reused = validate_hybrid_source_run(
            args.hybrid_source_run_root, run.root, calibration_libraries)
        source_payload["orchestrator"] = ORCHESTRATOR_NAME
        source_payload["orchestrator_release"] = RELEASE
        args.arms_bed = str(source_payload["arms_bed"])

    chromosome_arms = hybrid_chromosome_arms(
        args.gene_arms, args.include_sex_chromosomes)
    if len(calibration_libraries) > 999 or len(chromosome_arms) > 999:
        raise ValueError("hybrid array size exceeds the 999-task safety limit")

    stage_directories = {
        "REFERENCE": run.reference, "LEDGER": run.ledger,
        "PREPARE": run.prepare, "ASE": run.ase,
        "EXPRESSION": run.expression, "HYBRID_SHARD": run.shard,
        "MODEL": run.model, "CALL": run.call, "REPORT": run.report,
    }
    directories = [
        run.root, run.logs, run.scripts, run.manifests,
        os.path.join(run.manifests, "hybrid_model_inputs"),
    ]
    directories.extend(stage_directories[stage] for stage in selected)
    for directory in dict.fromkeys(directories):
        os.makedirs(directory, exist_ok=True)
    if source_payload is not None:
        publish_new_or_identical(
            run.source_reuse_manifest,
            json.dumps(source_payload, sort_keys=True, indent=2) + "\n",
            "hybrid source-reuse manifest")

    paths = make_hybrid_library_paths(
        args, run, calibration_libraries, selected, reused)
    supplied_by_source = {"REFERENCE", "LEDGER", "ASE"} if reused else set()
    ancestors = hybrid_omitted_ancestors(selected)
    for stage in ancestors:
        if stage in supplied_by_source:
            continue
        validate_completed_hybrid_stage(
            stage, run, paths, args, chromosome_arms)

    checkpoint_path = submission_path(run, selected, args.resume)
    previously_submitted: set[str] = set()
    if os.path.lexists(checkpoint_path):
        preliminary = read_json_object(
            checkpoint_path, "hybrid submission record")
        if (preliminary.get("run_root") != run.root
                or preliminary.get("selected_stages") != list(selected)
                or preliminary.get("orchestrator") != ORCHESTRATOR_NAME
                or not isinstance(preliminary.get("jobs"), dict)):
            raise ValueError(
                "hybrid submission record does not match this dedicated run: "
                f"{checkpoint_path}")
        previously_submitted = set(preliminary["jobs"])
    conflicts_by_stage = {
        stage: [path for path in hybrid_stage_outputs(
            stage, run, paths, args, chromosome_arms)
                if os.path.lexists(path)]
        for stage in selected if stage not in previously_submitted
    }
    conflicts_by_stage = {
        stage: values for stage, values in conflicts_by_stage.items() if values}
    if conflicts_by_stage:
        stage = next(stage for stage in STAGES if stage in conflicts_by_stage)
        raise ValueError(
            f"selected hybrid stage {stage} already has output path(s): "
            + ", ".join(conflicts_by_stage[stage][:3])
            + "; use a new --run-root")

    manifests = write_hybrid_manifests(
        args, run, paths, selected, chromosome_arms)
    scripts: dict[str, str] = {}
    if "REFERENCE" in selected:
        scripts["REFERENCE"] = render_script(
            hybrid_generated_script_path(
                run, "00_reference.sbatch", args.resume),
            reference_script(args, run), args)
    if "LEDGER" in selected:
        scripts["LEDGER"] = render_script(
            hybrid_generated_script_path(
                run, "01_ledger.sbatch", args.resume),
            ledger_script(args, run, calibration_libraries), args)
    if "PREPARE" in selected:
        scripts["PREPARE"] = render_script(
            hybrid_generated_script_path(
                run, "02_prepare_array.sbatch", args.resume),
            hybrid_prepare_script(
                args, run, manifests["PREPARE"], len(paths)), args)
    if "ASE" in selected:
        scripts["ASE"] = render_script(
            hybrid_generated_script_path(
                run, "03_ase_array.sbatch", args.resume),
            ase_script(args, run, manifests["ASE"], len(paths)), args)
    if "EXPRESSION" in selected:
        scripts["EXPRESSION"] = render_script(
            hybrid_generated_script_path(
                run, "04_expression_array.sbatch", args.resume),
            hybrid_expression_script(
                args, run, manifests["EXPRESSION"], len(paths)), args)
    if "HYBRID_SHARD" in selected:
        scripts["HYBRID_SHARD"] = render_script(
            hybrid_generated_script_path(
                run, "05_hybrid_shard_array.sbatch", args.resume),
            hybrid_shard_script(
                args, run, run.shard_tasks, len(paths), chromosome_arms), args)
    if "MODEL" in selected:
        scripts["MODEL"] = render_script(
            hybrid_generated_script_path(
                run, "06_model_array.sbatch", args.resume),
            hybrid_model_script(
                args, run, run.model_tasks, len(chromosome_arms)), args)
    if "CALL" in selected:
        scripts["CALL"] = render_script(
            hybrid_generated_script_path(
                run, "07_call.sbatch", args.resume),
            hybrid_call_script(args, run), args)
    if "REPORT" in selected:
        report_qc = [
            run.call_prefix + ".qc.tsv",
            *(os.path.join(
                run.shard, f"lib{path.library}.hybrid_shard_qc.tsv")
              for path in paths),
            *(hybrid_model_bundle_paths(run, chromosome)[3]
              for chromosome in chromosome_arms),
        ]
        scripts["REPORT"] = render_script(
            hybrid_generated_script_path(
                run, "08_report.sbatch", args.resume),
            hybrid_report_script(args, run, report_qc), args)

    plan = {
        "schema_version": "tetra_arm_hybrid_execution_plan_v1",
        "orchestrator": ORCHESTRATOR_NAME,
        "orchestrator_release": RELEASE,
        "workflow": WORKFLOW,
        "run_root": run.root,
        "libraries": libraries,
        "upstream_cohort_libraries": args.upstream_cohort_libraries_resolved,
        "selected_stages": list(selected),
        "source_reuse_manifest": (
            run.source_reuse_manifest if source_payload is not None else None),
        "chromosome_arms": {
            chromosome: list(arms)
            for chromosome, arms in chromosome_arms.items()},
        "task_manifests": dict(manifests),
        "scripts": scripts,
        "array_throttle": args.array_throttle,
        "nodelist": args.nodelist or None,
        "resources": {
            stage: {"cpus": cpus, "memory": memory}
            for stage, (cpus, memory) in stage_resources(args).items()
            if stage in selected},
        "partition": args.partition,
        "time": args.time,
        "status": "READY_TO_SUBMIT",
    }
    plan_path = publish_new_or_identical(
        hybrid_plan_path(run, selected, args.resume),
        json.dumps(plan, sort_keys=True, indent=2) + "\n",
        "hybrid execution plan")

    print(f"Tetraploid Arm CNV dedicated hybrid orchestrator {RELEASE}")
    print(f"Run root: {run.root}")
    print(f"Mapping input root: {args.mapping_input_root}")
    print(f"Upstream analysis root: {args.upstream_analysis_root}")
    print("Output libraries: " + ",".join(map(str, libraries)))
    print("Upstream calibration cohort: " + ",".join(
        map(str, args.upstream_cohort_libraries_resolved)))
    print("Selected stages: " + ",".join(selected))
    print("Chromosomes: " + ",".join(chromosome_arms))
    print(f"Arm BED: {args.arms_bed}")
    print(f"Gene-arm map: {args.gene_arms}")
    print(f"Partition/time: {args.partition} / {args.time}")
    print("Node placement: " + (args.nodelist or "scheduler default"))
    print("Array concurrency: " + (
        str(args.array_throttle)
        if args.array_throttle is not None else "scheduler default"))
    if source_payload is not None:
        print(f"Read-only source reuse: {run.source_reuse_manifest}")
    print(f"Execution plan: {plan_path}")
    print("Resolved resources:")
    resources = stage_resources(args)
    for stage in STAGES:
        if stage in selected:
            cpus, memory = resources[stage]
            print(f"  {stage:<14} cpus={cpus} mem={memory}")
    print("Task manifests:")
    for key in ("HYBRID_INPUTS", "PREPARE", "ASE", "EXPRESSION",
                "HYBRID_SHARD", "MODEL"):
        if key in manifests:
            print(f"  {key:<14} {manifests[key]}")
    for stage in STAGES:
        if stage in supplied_by_source:
            state = "validated source reuse"
        elif stage in selected:
            state = "planned"
        elif stage in ancestors:
            state = "validated complete"
        else:
            state = "not selected"
        print(f"  {stage:<14} {state}")
        if stage in scripts:
            print(f"                  script: {scripts[stage]}")

    existing_submission = load_submission(
        checkpoint_path, run, plan_path, selected, scripts)
    if not args.submit:
        if existing_submission is not None:
            print(
                "Submission checkpoint: "
                f"{existing_submission.get('status', 'PARTIAL')} "
                f"({checkpoint_path})")
        print("Planning only: no jobs submitted. Re-run with --submit to launch.")
        return 0

    lock_path = checkpoint_path + ".lock"
    with open(lock_path, "a+", encoding="utf-8") as lock_handle:
        fcntl.flock(lock_handle.fileno(), fcntl.LOCK_EX)
        submission = load_submission(
            checkpoint_path, run, plan_path, selected, scripts)
        if submission is None:
            submission = {
                "schema_version": "tetra_arm_hybrid_submission_v1",
                "orchestrator": ORCHESTRATOR_NAME,
                "orchestrator_release": RELEASE,
                "workflow": WORKFLOW,
                "run_root": run.root,
                "plan": plan_path,
                "selected_stages": list(selected),
                "scripts": dict(scripts),
                "jobs": {},
                "status": "PARTIAL",
            }
        jobs = submission["jobs"]
        job_ids = {
            stage: str(jobs[stage]["job_id"])
            for stage in STAGES if stage in jobs}
        for stage in STAGES:
            if stage not in scripts:
                continue
            if stage in job_ids:
                print(
                    f"Already submitted {stage:<14} job {job_ids[stage]}; "
                    f"checkpoint={checkpoint_path}")
                continue
            dependencies = hybrid_resolved_dependencies(
                stage, selected, job_ids)
            job_id = submit_job(
                scripts[stage], dependencies, run.root)
            job_ids[stage] = job_id
            jobs[stage] = {
                "job_id": job_id,
                "script": scripts[stage],
                "dependencies": dependencies,
            }
            submission["status"] = "PARTIAL"
            replace_mutable_json_atomic(checkpoint_path, submission)
            print(
                f"Submitted {stage:<14} job {job_id}; afterok="
                f"{':'.join(dependencies) if dependencies else 'none'}; "
                f"logs={run.logs}")
        submission["status"] = "COMPLETE"
        replace_mutable_json_atomic(checkpoint_path, submission)
    return 0


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        return run_workflow(args)
    except (OSError, ValueError, RuntimeError,
            subprocess.SubprocessError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
