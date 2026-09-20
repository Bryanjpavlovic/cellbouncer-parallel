#!/usr/bin/env python3
"""Render and submit the Tet 2025 RNA-anchored SCENIC+ job graph."""

from __future__ import annotations

import argparse
import csv
import os
import re
import shlex
import stat
import subprocess
import sys
from pathlib import Path


PROJECT_ROOT = Path("/mnt/beegfs/tetraploid_multiome_cis_trans")
ATAC_MAPPING_ROOT = PROJECT_ROOT / "ATAC" / "mapping_output"
ATAC_ANALYSIS_ROOT = PROJECT_ROOT / "ATAC" / "analysis"
DEFAULT_RUN_LABEL = "scenicplus_batch2_batch3_pools_v1"
DEFAULT_LIBRARIES = (15, 16, 23, 24, 31, 32, 39, 40)
DEFAULT_WORKER = Path(
    "/nvme/software/packages/cellbouncer/dev/bin/tet2025_scenicplus.py"
)
DEFAULT_RESOURCE_ROOT = (
    ATAC_ANALYSIS_ROOT
    / "shared_resources"
    / "cistarget"
    / "hg38_screen_v10_clust"
)
HAL_PATH = Path(
    "/mnt/beegfs/genomes_annotations/ancestral_genomes/Anc1_geno_graph/"
    "SharedData/AncestralGenomeV1/AncestralGenomeV1_primates_eichler-V2.hal"
)
MAO_FASTA = Path(
    "/mnt/beegfs/genomes_annotations/ancestral_genomes/litterbox/"
    "human_chimp_bonobo/human_chimp_bonobo_filt_numtmask.fa.gz"
)
MAO_GTF = Path(
    "/mnt/beegfs/genomes_annotations/ancestral_genomes/litterbox/"
    "human_chimp_bonobo/human_chimp_bonobo.gtf.gz"
)
STAGES = (
    "REFERENCE",
    "RESOURCES",
    "PEAKS",
    "CONSENSUS",
    "CISTOPIC",
    "MERGE",
    "TOPICS",
    "REGIONS",
    "SCENICPLUS",
)
PEAKS_ONWARD_STAGES = STAGES[2:]
CISTOPIC_ONWARD_STAGES = STAGES[4:]
MERGE_ONWARD_STAGES = STAGES[5:]
# All eight 2026-09-17 library jobs exhausted a 96G cgroup while PyRanges was
# materializing the parallel fragment/region join.  pycisTopic's partition
# fallback happens only after that join, so this stage needs node-scale memory.
DEFAULT_CISTOPIC_MEMORY = "500G"


def quote(value: str | Path) -> str:
    return shlex.quote(str(value))


def log(message: str) -> None:
    print(message, flush=True)


def require_file(path: Path, description: str) -> Path:
    path = path.expanduser().resolve()
    if not path.is_file() or path.stat().st_size == 0:
        raise FileNotFoundError(f"{description} is missing or empty: {path}")
    return path


def write_text(path: Path, text: str, executable: bool = False) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + f".tmp.{os.getpid()}")
    try:
        with open(temporary, "w", encoding="utf-8", newline="") as handle:
            handle.write(text)
        os.replace(temporary, path)
        if executable:
            path.chmod(path.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP)
    finally:
        temporary.unlink(missing_ok=True)


def base_header(
    *,
    name: str,
    cpus: int,
    memory: str,
    time: str,
    partition: str,
    stdout: Path,
    stderr: Path,
    array: str | None = None,
) -> str:
    lines = [
        "#!/bin/bash",
        f"#SBATCH --job-name={name}",
        f"#SBATCH --partition={partition}",
        "#SBATCH --nodes=1",
        "#SBATCH --ntasks=1",
        f"#SBATCH --cpus-per-task={cpus}",
        f"#SBATCH --mem={memory}",
        f"#SBATCH --time={time}",
        "#SBATCH --chdir=/tmp",
        f"#SBATCH --output={stdout}",
        f"#SBATCH --error={stderr}",
    ]
    if array is not None:
        lines.append(f"#SBATCH --array={array}")
    lines.extend(
        [
            "",
            "set -euo pipefail",
            "",
        ]
    )
    return "\n".join(lines)


def scenic_modules(with_hal: bool = False) -> str:
    lines = ["module purge"]
    if with_hal:
        lines.extend(
            [
                "module load hal/latest",
                "command -v halStats >/dev/null",
                "command -v halLiftover >/dev/null",
            ]
        )
    lines.extend(
        [
            "module load miniforge/3 scenicplus/latest",
            "module list 2>&1",
            "command -v python >/dev/null",
            "python - <<'PY'",
            "import pycisTopic",
            "import scenicplus",
            "PY",
            "",
            "export OMP_NUM_THREADS=1",
            "export OPENBLAS_NUM_THREADS=1",
            "export MKL_NUM_THREADS=1",
            "export NUMEXPR_NUM_THREADS=1",
            "export BLIS_NUM_THREADS=1",
            "export POLARS_MAX_THREADS=1",
            "",
        ]
    )
    return "\n".join(lines)


def render_jobs(args, root: Path, task_table: Path) -> dict[str, Path]:
    scripts = root / "control" / "slurm_scripts"
    logs = root / "control" / "logs"
    temporary = root / "tmp"
    reference = root / "reference"
    peaks = root / "peaks"
    cistopic = root / "cistopic"
    region_sets = root / "region_sets"
    pipeline = root / "scenicplus_pipeline"
    input_root = root / "input"
    for directory in (
        scripts,
        logs,
        temporary,
        reference,
        peaks,
        cistopic / "per_library",
        region_sets,
        input_root,
    ):
        directory.mkdir(parents=True, exist_ok=True)

    worker = args.worker.expanduser().resolve()
    cell_metadata = input_root / "cell_metadata.tsv.gz"
    rna_h5ad = input_root / "rna_for_scenicplus.h5ad"
    jobs: dict[str, Path] = {}

    path = scripts / "01_reference.sbatch"
    text = base_header(
        name="tet_sp_ref",
        cpus=4,
        memory="32G",
        time="12:00:00",
        partition=args.partition,
        stdout=logs / "tet_sp_ref_%j.out",
        stderr=logs / "tet_sp_ref_%j.err",
    )
    text += scenic_modules(with_hal=True)
    text += (
        f"python {quote(worker)} reference \\\n"
        f"  --hal {quote(args.hal)} \\\n"
        f"  --mao-fasta {quote(args.mao_fasta)} \\\n"
        f"  --mao-gtf {quote(args.mao_gtf)} \\\n"
        f"  --output-dir {quote(reference)}\n"
    )
    write_text(path, text, executable=True)
    jobs["REFERENCE"] = path

    path = scripts / "02_resources.sbatch"
    text = base_header(
        name="tet_sp_db",
        cpus=2,
        memory="16G",
        time="7-00:00:00",
        partition=args.partition,
        stdout=logs / "tet_sp_db_%j.out",
        stderr=logs / "tet_sp_db_%j.err",
    )
    text += scenic_modules()
    text += (
        f"python {quote(worker)} resources \\\n"
        f"  --output-dir {quote(args.resources_root)}\n"
    )
    write_text(path, text, executable=True)
    jobs["RESOURCES"] = path

    path = scripts / "03_pseudobulk_peaks.sbatch"
    text = base_header(
        name="tet_sp_peak",
        cpus=32,
        memory="512G",
        time="2-00:00:00",
        partition=args.partition,
        stdout=logs / "tet_sp_peak_%j.out",
        stderr=logs / "tet_sp_peak_%j.err",
    )
    text += scenic_modules()
    text += (
        "if command -v macs2 >/dev/null; then\n"
        "  :\n"
        "elif command -v macs3 >/dev/null; then\n"
        "  :\n"
        "else\n"
        "  echo 'ERROR: neither macs2 nor macs3 is available' >&2\n"
        "  exit 1\n"
        "fi\n\n"
    )
    text += (
        "PEAK_TMP=$(mktemp -d /tmp/r.XXXXXX)\n"
        "trap 'rm -rf -- \"$PEAK_TMP\"' EXIT\n\n"
    )
    text += (
        f"python {quote(worker)} pseudobulk-peaks \\\n"
        f"  --cell-metadata {quote(cell_metadata)} \\\n"
        f"  --task-table {quote(task_table)} \\\n"
        f"  --mao-chromsizes {quote(reference / 'mao_chromsizes.tsv')} \\\n"
        f"  --output-dir {quote(peaks)} \\\n"
        "  --temp-dir \"$PEAK_TMP\" \\\n"
        "  --cpus ${SLURM_CPUS_PER_TASK}\n"
    )
    write_text(path, text, executable=True)
    jobs["PEAKS"] = path

    path = scripts / "04_consensus_liftover.sbatch"
    text = base_header(
        name="tet_sp_lift",
        cpus=16,
        memory="128G",
        time="1-00:00:00",
        partition=args.partition,
        stdout=logs / "tet_sp_lift_%j.out",
        stderr=logs / "tet_sp_lift_%j.err",
    )
    text += scenic_modules(with_hal=True)
    text += (
        f"python {quote(worker)} consensus-liftover \\\n"
        f"  --hal {quote(args.hal)} \\\n"
        f"  --narrow-peaks {quote(peaks / 'narrow_peaks.pkl')} \\\n"
        f"  --mao-chromsizes {quote(reference / 'mao_chromsizes.tsv')} \\\n"
        f"  --output-dir {quote(peaks / 'consensus')}\n"
    )
    write_text(path, text, executable=True)
    jobs["CONSENSUS"] = path

    path = scripts / "05_cistopic_libraries.sbatch"
    text = base_header(
        name="tet_sp_ct",
        cpus=8,
        memory=args.cistopic_memory,
        time="2-00:00:00",
        partition=args.partition,
        stdout=logs / "tet_sp_ct_%A_%a.out",
        stderr=logs / "tet_sp_ct_%A_%a.err",
        array=f"0-{len(args.libraries) - 1}%{len(args.libraries)}",
    )
    text += scenic_modules()
    text += (
        f"TASK_TABLE={quote(task_table)}\n"
        "LINE_NUMBER=$((SLURM_ARRAY_TASK_ID + 2))\n"
        "LINE=$(sed -n \"${LINE_NUMBER}p\" \"${TASK_TABLE}\")\n"
        "if [[ -z \"${LINE}\" ]]; then\n"
        "  echo \"ERROR: no task row for array index ${SLURM_ARRAY_TASK_ID}\" >&2\n"
        "  exit 1\n"
        "fi\n"
        "IFS=$'\\t' read -r task_index library library_number fragments cistopic_output <<< \"${LINE}\"\n"
        "cistopic_output=${cistopic_output%$'\\r'}\n"
        f"python {quote(worker)} build-cistopic-library \\\n"
        "  --library \"${library}\" \\\n"
        "  --fragments \"${fragments}\" \\\n"
        f"  --regions {quote(peaks / 'consensus' / 'consensus_peaks_mao_retained.bed')} \\\n"
        f"  --cell-metadata {quote(cell_metadata)} \\\n"
        "  --output \"${cistopic_output}\" \\\n"
        "  --cpus ${SLURM_CPUS_PER_TASK}\n"
    )
    write_text(path, text, executable=True)
    jobs["CISTOPIC"] = path

    path = scripts / "06_merge_cistopic.sbatch"
    text = base_header(
        name="tet_sp_merge",
        cpus=16,
        memory="256G",
        time="1-00:00:00",
        partition=args.partition,
        stdout=logs / "tet_sp_merge_%j.out",
        stderr=logs / "tet_sp_merge_%j.err",
    )
    text += scenic_modules()
    text += (
        f"python {quote(worker)} merge-cistopic \\\n"
        f"  --task-table {quote(task_table)} \\\n"
        f"  --peak-crosswalk {quote(peaks / 'consensus' / 'peak_liftover_crosswalk.tsv.gz')} \\\n"
        f"  --cell-metadata {quote(cell_metadata)} \\\n"
        f"  --output {quote(cistopic / 'merged_cistopic_hg38.pkl')}\n"
    )
    write_text(path, text, executable=True)
    jobs["MERGE"] = path

    path = scripts / "07_topic_model.sbatch"
    text = base_header(
        name="tet_sp_topic",
        cpus=4,
        memory="384G",
        time="4-00:00:00",
        partition=args.partition,
        stdout=logs / "tet_sp_topic_%j.out",
        stderr=logs / "tet_sp_topic_%j.err",
    )
    text += scenic_modules()
    text += (
        "TOPIC_TMP=$(mktemp -d /tmp/r.XXXXXX)\n"
        "trap 'rm -rf -- \"$TOPIC_TMP\"' EXIT\n\n"
    )
    text += (
        f"python {quote(worker)} topic-model \\\n"
        f"  --input {quote(cistopic / 'merged_cistopic_hg38.pkl')} \\\n"
        f"  --output {quote(cistopic / 'cistopic_with_topics.pkl')} \\\n"
        "  --temp-dir \"$TOPIC_TMP\" \\\n"
        "  --cpus ${SLURM_CPUS_PER_TASK} \\\n"
        "  --topics 20 30 40 50\n"
    )
    write_text(path, text, executable=True)
    jobs["TOPICS"] = path

    path = scripts / "08_region_sets.sbatch"
    text = base_header(
        name="tet_sp_regions",
        cpus=8,
        memory="384G",
        time="3-00:00:00",
        partition=args.partition,
        stdout=logs / "tet_sp_regions_%j.out",
        stderr=logs / "tet_sp_regions_%j.err",
    )
    text += scenic_modules()
    text += (
        "REGION_TMP=$(mktemp -d /tmp/r.XXXXXX)\n"
        "trap 'rm -rf -- \"$REGION_TMP\"' EXIT\n\n"
    )
    text += (
        f"python {quote(worker)} region-sets \\\n"
        f"  --input {quote(cistopic / 'cistopic_with_topics.pkl')} \\\n"
        f"  --output-dir {quote(region_sets)} \\\n"
        f"  --output-cistopic {quote(cistopic / 'cistopic_scenicplus.pkl')} \\\n"
        "  --temp-dir \"$REGION_TMP\" \\\n"
        "  --cpus ${SLURM_CPUS_PER_TASK} \\\n"
        "  --grouping rna_leiden_reference\n"
    )
    write_text(path, text, executable=True)
    jobs["REGIONS"] = path

    path = scripts / "09_scenicplus.sbatch"
    text = base_header(
        name="tet_scenicplus",
        cpus=48,
        memory="768G",
        time="7-00:00:00",
        partition=args.partition,
        stdout=logs / "tet_scenicplus_%j.out",
        stderr=logs / "tet_scenicplus_%j.err",
    )
    text += scenic_modules()
    text += "command -v scenicplus >/dev/null\ncommand -v snakemake >/dev/null\n\n"
    text += (
        "SCENICPLUS_TMP=$(mktemp -d /tmp/r.XXXXXX)\n"
        "trap 'rm -rf -- \"$SCENICPLUS_TMP\"' EXIT\n\n"
    )
    text += (
        f"python {quote(worker)} run-scenicplus \\\n"
        f"  --cistopic {quote(cistopic / 'cistopic_scenicplus.pkl')} \\\n"
        f"  --rna-h5ad {quote(rna_h5ad)} \\\n"
        f"  --region-sets {quote(region_sets)} \\\n"
        f"  --resources {quote(args.resources_root)} \\\n"
        f"  --genome-annotation {quote(reference / 'genome_annotation.tsv')} \\\n"
        f"  --hg38-chromsizes {quote(reference / 'hg38_chromsizes.tsv')} \\\n"
        f"  --pipeline-root {quote(pipeline)} \\\n"
        "  --temp-dir \"$SCENICPLUS_TMP\" \\\n"
        "  --cpus ${SLURM_CPUS_PER_TASK}\n"
    )
    write_text(path, text, executable=True)
    jobs["SCENICPLUS"] = path
    return jobs


def write_task_table(root: Path, libraries: list[int]) -> Path:
    output = root / "control" / "library_tasks.tsv"
    output.parent.mkdir(parents=True, exist_ok=True)
    rows = []
    for task_index, number in enumerate(libraries):
        library = f"lib{number}"
        fragments = (
            ATAC_MAPPING_ROOT
            / f"Tet_2025_Multiome-ATAC_{number}"
            / "atac_fragments.tsv.gz"
        )
        require_file(fragments, f"{library} ATAC fragments")
        rows.append(
            {
                "task_index": task_index,
                "library": library,
                "library_number": number,
                "fragments": str(fragments),
                "cistopic_output": str(
                    root / "cistopic" / "per_library" / f"{library}.pkl"
                ),
            }
        )
    temporary = output.with_name(output.name + f".tmp.{os.getpid()}")
    try:
        with open(temporary, "w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=list(rows[0]),
                delimiter="\t",
                lineterminator="\n",
            )
            writer.writeheader()
            writer.writerows(rows)
        os.replace(temporary, output)
    finally:
        temporary.unlink(missing_ok=True)
    return output


def submit(script: Path, dependencies: list[str] | None = None) -> str:
    command = ["sbatch", "--parsable"]
    if dependencies:
        command.append("--dependency=afterok:" + ":".join(dependencies))
    command.append(str(script))
    result = subprocess.run(
        command,
        check=True,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    job_id = result.stdout.strip().split(";", 1)[0]
    if not re.fullmatch(r"\d+", job_id):
        raise RuntimeError(f"could not parse sbatch job id from: {result.stdout!r}")
    return job_id


def submit_selected(stage: str, jobs: dict[str, Path]) -> dict[str, str]:
    submitted: dict[str, str] = {}
    if stage not in {
        "ALL",
        "PEAKS_ONWARD",
        "CISTOPIC_ONWARD",
        "MERGE_ONWARD",
    }:
        submitted[stage] = submit(jobs[stage])
        return submitted
    if stage == "MERGE_ONWARD":
        submitted["MERGE"] = submit(jobs["MERGE"])
        submitted["TOPICS"] = submit(jobs["TOPICS"], [submitted["MERGE"]])
        submitted["REGIONS"] = submit(jobs["REGIONS"], [submitted["TOPICS"]])
        submitted["SCENICPLUS"] = submit(
            jobs["SCENICPLUS"], [submitted["REGIONS"]]
        )
        return submitted
    if stage == "CISTOPIC_ONWARD":
        submitted["CISTOPIC"] = submit(jobs["CISTOPIC"])
        submitted["MERGE"] = submit(jobs["MERGE"], [submitted["CISTOPIC"]])
        submitted["TOPICS"] = submit(jobs["TOPICS"], [submitted["MERGE"]])
        submitted["REGIONS"] = submit(jobs["REGIONS"], [submitted["TOPICS"]])
        submitted["SCENICPLUS"] = submit(
            jobs["SCENICPLUS"], [submitted["REGIONS"]]
        )
        return submitted
    if stage == "PEAKS_ONWARD":
        submitted["PEAKS"] = submit(jobs["PEAKS"])
        submitted["CONSENSUS"] = submit(
            jobs["CONSENSUS"], [submitted["PEAKS"]]
        )
        submitted["CISTOPIC"] = submit(
            jobs["CISTOPIC"], [submitted["CONSENSUS"]]
        )
        submitted["MERGE"] = submit(jobs["MERGE"], [submitted["CISTOPIC"]])
        submitted["TOPICS"] = submit(jobs["TOPICS"], [submitted["MERGE"]])
        submitted["REGIONS"] = submit(jobs["REGIONS"], [submitted["TOPICS"]])
        submitted["SCENICPLUS"] = submit(
            jobs["SCENICPLUS"], [submitted["REGIONS"]]
        )
        return submitted
    submitted["REFERENCE"] = submit(jobs["REFERENCE"])
    submitted["RESOURCES"] = submit(jobs["RESOURCES"])
    submitted["PEAKS"] = submit(
        jobs["PEAKS"], [submitted["REFERENCE"]]
    )
    submitted["CONSENSUS"] = submit(
        jobs["CONSENSUS"], [submitted["PEAKS"]]
    )
    submitted["CISTOPIC"] = submit(
        jobs["CISTOPIC"], [submitted["CONSENSUS"]]
    )
    submitted["MERGE"] = submit(jobs["MERGE"], [submitted["CISTOPIC"]])
    submitted["TOPICS"] = submit(jobs["TOPICS"], [submitted["MERGE"]])
    submitted["REGIONS"] = submit(jobs["REGIONS"], [submitted["TOPICS"]])
    submitted["SCENICPLUS"] = submit(
        jobs["SCENICPLUS"], [submitted["REGIONS"], submitted["RESOURCES"]]
    )
    return submitted


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Submit the RNA-anchored Tet 2025 SCENIC+ analysis"
    )
    parser.add_argument("--run-label", default=DEFAULT_RUN_LABEL)
    parser.add_argument(
        "--libraries", type=int, nargs="+", default=list(DEFAULT_LIBRARIES)
    )
    parser.add_argument("--worker", type=Path, default=DEFAULT_WORKER)
    parser.add_argument("--resources-root", type=Path, default=DEFAULT_RESOURCE_ROOT)
    parser.add_argument("--hal", type=Path, default=HAL_PATH)
    parser.add_argument("--mao-fasta", type=Path, default=MAO_FASTA)
    parser.add_argument("--mao-gtf", type=Path, default=MAO_GTF)
    parser.add_argument("--partition", default="compute")
    parser.add_argument(
        "--cistopic-memory",
        default=DEFAULT_CISTOPIC_MEMORY,
        help="SLURM memory requested by each per-library cisTopic job (default: 500G)",
    )
    parser.add_argument(
        "--stage",
        choices=("ALL", "PEAKS_ONWARD", "CISTOPIC_ONWARD", "MERGE_ONWARD")
        + STAGES,
        default="ALL",
    )
    parser.add_argument("--submit", action="store_true")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    try:
        if not re.fullmatch(r"[A-Za-z0-9_.-]+", args.run_label):
            raise ValueError("run label may contain only letters, numbers, dot, dash, underscore")
        if not re.fullmatch(r"[1-9][0-9]*[KMGT]", args.cistopic_memory, re.IGNORECASE):
            raise ValueError(
                "cisTopic memory must be a positive SLURM value such as 500G"
            )
        args.libraries = sorted(set(args.libraries))
        if not args.libraries:
            raise ValueError("at least one library is required")
        root = (ATAC_ANALYSIS_ROOT / args.run_label).resolve()
        args.worker = require_file(args.worker, "deployed SCENIC+ worker")
        args.hal = require_file(args.hal, "authoritative HAL alignment")
        args.mao_fasta = require_file(args.mao_fasta, "Mao FASTA")
        args.mao_gtf = require_file(args.mao_gtf, "Mao GTF")
        args.resources_root = args.resources_root.expanduser().resolve()

        # These are produced by prepare_scenicplus_rna.py on the workstation and
        # copied directly into this run through the existing BeeGFS mount.
        require_file(root / "input" / "rna_for_scenicplus.h5ad", "SCENIC+ RNA export")
        require_file(root / "input" / "cell_metadata.tsv.gz", "paired-cell metadata")
        task_table = write_task_table(root, args.libraries)
        jobs = render_jobs(args, root, task_table)
        if args.stage == "ALL":
            selected = STAGES
        elif args.stage == "PEAKS_ONWARD":
            selected = PEAKS_ONWARD_STAGES
        elif args.stage == "CISTOPIC_ONWARD":
            selected = CISTOPIC_ONWARD_STAGES
        elif args.stage == "MERGE_ONWARD":
            selected = MERGE_ONWARD_STAGES
        else:
            selected = (args.stage,)
        if args.stage in {
            "PEAKS_ONWARD",
            "CISTOPIC_ONWARD",
            "MERGE_ONWARD",
        }:
            for path, description in (
                (root / "reference" / "mao_chromsizes.tsv", "Mao chromosome sizes"),
                (root / "reference" / "hg38_chromsizes.tsv", "hg38 chromosome sizes"),
                (root / "reference" / "genome_annotation.tsv", "projected gene annotation"),
                (
                    args.resources_root
                    / "hg38_screen_v10_clust.regions_vs_motifs.rankings.feather",
                    "cisTarget rankings database",
                ),
                (
                    args.resources_root
                    / "hg38_screen_v10_clust.regions_vs_motifs.scores.feather",
                    "cisTarget scores database",
                ),
                (
                    args.resources_root
                    / "motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl",
                    "motif-to-TF annotation",
                ),
            ):
                require_file(path, description)
        if args.stage in {"CISTOPIC_ONWARD", "MERGE_ONWARD"}:
            require_file(
                root
                / "peaks"
                / "consensus"
                / "peak_liftover_crosswalk.tsv.gz",
                "consensus-peak liftover crosswalk",
            )
        if args.stage == "CISTOPIC_ONWARD":
            require_file(
                root
                / "peaks"
                / "consensus"
                / "consensus_peaks_mao_retained.bed",
                "retained Mao consensus peaks",
            )
        if args.stage == "MERGE_ONWARD":
            for number in args.libraries:
                require_file(
                    root
                    / "cistopic"
                    / "per_library"
                    / f"lib{number}.pkl",
                    f"lib{number} cisTopic object",
                )
        log(f"Libraries: {', '.join(map(str, args.libraries))}")
        log(f"Analysis root: {root}")
        if "CISTOPIC" in selected:
            log(f"cisTopic memory per library: {args.cistopic_memory}")
        log(f"Authoritative HAL: {args.hal}")
        log("HAL direction: human_chimp_bonobo -> Human (reciprocal check enabled)")
        log(f"Task table: {task_table}")
        for stage in selected:
            log(f"Rendered {stage}: {jobs[stage]}")
        if not args.submit:
            log("Jobs were rendered but not submitted; add --submit to send them to SLURM")
            return 0
        submitted = submit_selected(args.stage, jobs)
        for stage, job_id in submitted.items():
            log(f"Submitted {stage}: job {job_id}")
        if args.stage == "ALL":
            log(
                "Dependency chain: REFERENCE and RESOURCES start together; "
                "PEAKS -> CONSENSUS -> CISTOPIC[] -> MERGE -> TOPICS -> "
                "REGIONS -> SCENICPLUS"
            )
        elif args.stage == "PEAKS_ONWARD":
            log(
                "Dependency chain: PEAKS -> CONSENSUS -> CISTOPIC[] -> "
                "MERGE -> TOPICS -> REGIONS -> SCENICPLUS"
            )
        elif args.stage == "CISTOPIC_ONWARD":
            log(
                "Dependency chain: CISTOPIC[] -> MERGE -> TOPICS -> "
                "REGIONS -> SCENICPLUS"
            )
        elif args.stage == "MERGE_ONWARD":
            log(
                "Dependency chain: MERGE -> TOPICS -> REGIONS -> SCENICPLUS"
            )
        return 0
    except subprocess.CalledProcessError as exc:
        print(
            f"ERROR: command failed with exit {exc.returncode}: "
            + " ".join(map(str, exc.cmd)),
            file=sys.stderr,
        )
        if exc.stderr:
            print(exc.stderr, file=sys.stderr)
        return int(exc.returncode or 1)
    except Exception as exc:
        print(f"ERROR: {type(exc).__name__}: {exc}", file=sys.stderr, flush=True)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
