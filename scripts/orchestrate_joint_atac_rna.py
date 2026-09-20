#!/usr/bin/env python3
"""Render and optionally submit the Tet 2025 ATAC/RNA clustering analysis.

The orchestrator discovers paired libraries from the canonical mapping roots,
renders one per-library array plus two cohort jobs, and submits a single SLURM
dependency chain when --submit is supplied.
"""

from __future__ import annotations

import argparse
import csv
import os
from pathlib import Path
import re
import shlex
import stat
import subprocess
import sys
from typing import Iterable, Sequence


VERSION = "1.0.3"

PROJECT_ROOT = Path("/mnt/beegfs/tetraploid_multiome_cis_trans")
ATAC_MAPPING_ROOT = PROJECT_ROOT / "ATAC" / "mapping_output"
RNA_MAPPING_ROOT = PROJECT_ROOT / "3P" / "mapping_output"
ATAC_ANALYSIS_BASE = PROJECT_ROOT / "ATAC" / "analysis"
ATAC_FIGURE_BASE = PROJECT_ROOT / "ATAC" / "figures" / "all40"
RUN_LABEL = "joint_atac_rna_clustering_v1"

IDENTITY_TABLE = (
    PROJECT_ROOT
    / "3P"
    / "analysis"
    / "aggregate_library_analysis"
    / "identity_reconciliation"
    / "aggregate"
    / "identity_assignments.tsv.gz"
)
WORKBOOK = PROJECT_ROOT / "Misc_Metadata" / "Library_conversions.xlsx"
REFERENCE_ROOT = Path(
    "/mnt/beegfs/genomes_annotations/ancestral_genomes/litterbox/"
    "human_chimp_bonobo"
)
REFERENCE_FASTA = REFERENCE_ROOT / "human_chimp_bonobo_filt_numtmask.fa.gz"
REFERENCE_GTF = REFERENCE_ROOT / "human_chimp_bonobo.gtf.gz"

SNAP_MODULES = ("miniforge/3", "snapatac2/2.10.0")


def quote(value: object) -> str:
    return shlex.quote(str(value))


def write_if_changed(path: Path, content: str, executable: bool = False) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not path.exists() or path.read_text(encoding="utf-8") != content:
        temporary = path.with_name(path.name + ".tmp")
        temporary.write_text(content, encoding="utf-8")
        os.replace(temporary, path)
    if executable:
        path.chmod(
            stat.S_IRWXU
            | stat.S_IRGRP
            | stat.S_IXGRP
            | stat.S_IROTH
            | stat.S_IXOTH
        )
    return path


def natural_key(value: str) -> list[object]:
    return [int(x) if x.isdigit() else x.lower() for x in re.split(r"(\d+)", value)]


def discover_numbered_directories(root: Path, prefix: str) -> dict[int, Path]:
    result: dict[int, Path] = {}
    if not root.is_dir():
        raise FileNotFoundError(f"mapping root is missing: {root}")
    pattern = re.compile(rf"^{re.escape(prefix)}_(\d+)$")
    for candidate in root.iterdir():
        if not candidate.is_dir():
            continue
        match = pattern.fullmatch(candidate.name)
        if match:
            result[int(match.group(1))] = candidate.resolve()
    return result


def discover_libraries(
    atac_root: Path,
    rna_root: Path,
    requested: Sequence[int] | None = None,
) -> list[dict[str, object]]:
    """Discover complete ATAC/RNA pairs without a hard-coded library list."""
    atac_dirs = discover_numbered_directories(
        atac_root, "Tet_2025_Multiome-ATAC"
    )
    rna_dirs = discover_numbered_directories(
        rna_root, "Tet_2025_Multiome-RNA"
    )
    numbers = sorted(set(atac_dirs) | set(rna_dirs))
    if requested:
        wanted = set(requested)
        absent = sorted(wanted - set(numbers))
        if absent:
            raise FileNotFoundError(
                "requested libraries were not discovered: "
                + ", ".join(map(str, absent))
            )
        numbers = [number for number in numbers if number in wanted]

    rows: list[dict[str, object]] = []
    failures: list[str] = []
    for number in numbers:
        atac_dir = atac_dirs.get(number)
        rna_dir = rna_dirs.get(number)
        if atac_dir is None:
            failures.append(f"lib{number}: ATAC directory missing")
            continue
        if rna_dir is None:
            failures.append(f"lib{number}: RNA directory missing")
            continue
        fragments = atac_dir / "atac_fragments.tsv.gz"
        filtered = rna_dir / "filtered"
        required = (
            fragments,
            filtered / "barcodes.tsv.gz",
            filtered / "features.tsv.gz",
            filtered / "matrix.mtx.gz",
        )
        missing = [str(path) for path in required if not path.is_file()]
        empty = [str(path) for path in required if path.is_file() and path.stat().st_size == 0]
        if missing or empty:
            detail = []
            if missing:
                detail.append("missing=" + ",".join(missing))
            if empty:
                detail.append("empty=" + ",".join(empty))
            failures.append(f"lib{number}: " + "; ".join(detail))
            continue
        rows.append(
            {
                "library": f"lib{number}",
                "library_number": number,
                "fragments": fragments.resolve(),
                "rna_matrix_dir": filtered.resolve(),
            }
        )

    if failures:
        raise FileNotFoundError(
            "incomplete paired libraries:\n  " + "\n  ".join(failures)
        )
    if not rows:
        raise FileNotFoundError(
            f"no complete paired libraries found under {atac_root} and {rna_root}"
        )
    return rows


def resolve_chrom_sizes(explicit: Path | None, fasta: Path) -> Path:
    candidates: list[Path] = []
    if explicit is not None:
        candidates.append(explicit)
    else:
        candidates.append(Path(str(fasta) + ".fai"))
        if fasta.suffix == ".gz":
            candidates.append(fasta.with_suffix("").with_suffix(".fa.fai"))
            candidates.append(fasta.with_suffix("").with_suffix(".fai"))
        candidates.extend(
            [
                fasta.parent / "chrom.sizes",
                fasta.parent / "chrom_sizes.tsv",
            ]
        )
    for candidate in candidates:
        if candidate.is_file() and candidate.stat().st_size > 0:
            return candidate.resolve()
    raise FileNotFoundError(
        "a chromosome-size table is required; checked: "
        + ", ".join(str(path) for path in candidates)
    )


def write_task_table(path: Path, rows: Sequence[dict[str, object]], analysis_root: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = (
        "library",
        "library_number",
        "fragments",
        "rna_matrix_dir",
        "library_output_dir",
    )
    temporary = path.with_name(path.name + ".tmp")
    with temporary.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            output = analysis_root / "per_library" / str(row["library"])
            writer.writerow(
                {
                    **row,
                    "library_output_dir": output.resolve(),
                }
            )
    os.replace(temporary, path)
    return path


def module_preflight() -> str:
    loads = "\n".join(f"module load {module}" for module in SNAP_MODULES)
    return f"""module purge
{loads}
module list 2>&1
command -v python >/dev/null
command -v macs3 >/dev/null
python - <<'PY'
import anndata
import scanpy
import snapatac2
print("SnapATAC2", snapatac2.__version__)
PY
"""


def sbatch_header(
    *,
    job_name: str,
    log_pattern: Path,
    partition: str,
    cpus: int,
    memory: str,
    duration: str,
    working_directory: Path,
    array: str | None = None,
) -> str:
    array_line = f"#SBATCH --array={array}\n" if array else ""
    return f"""#!/bin/bash
#SBATCH --job-name={job_name}
#SBATCH --output={log_pattern}.out
#SBATCH --error={log_pattern}.err
{array_line}#SBATCH --time={duration}
#SBATCH --cpus-per-task={cpus}
#SBATCH --mem={memory}
#SBATCH --partition={partition}
#SBATCH --nodes=1
#SBATCH --chdir={working_directory}

set -euo pipefail

export OMP_NUM_THREADS="$SLURM_CPUS_PER_TASK"
export OPENBLAS_NUM_THREADS="$SLURM_CPUS_PER_TASK"
export MKL_NUM_THREADS="$SLURM_CPUS_PER_TASK"
export NUMEXPR_NUM_THREADS="$SLURM_CPUS_PER_TASK"
export MPLBACKEND=Agg

""" + module_preflight()


def render_import_script(args: argparse.Namespace, paths: dict[str, Path], task_count: int) -> str:
    concurrency = min(args.array_concurrency, task_count)
    force = " --force" if args.force else ""
    exclusions = ",".join(args.exclude_chrom)
    return sbatch_header(
        job_name="tet_atac_import",
        log_pattern=paths["log_dir"] / "tet_atac_import_%A_%a",
        partition=args.partition,
        cpus=args.import_cpus,
        memory=args.import_memory,
        duration=args.time,
        working_directory=paths["analysis_root"],
        array=f"0-{task_count - 1}%{concurrency}",
    ) + f"""
TASK_TABLE={quote(paths['task_table'])}
WORKER={quote(paths['worker'])}
CHROM_SIZES={quote(paths['chrom_sizes'])}
GTF={quote(args.gtf)}

for required in "$TASK_TABLE" "$WORKER" "$CHROM_SIZES"; do
    if [[ ! -s "$required" ]]; then
        echo "ERROR: required input missing or empty: $required" >&2
        exit 1
    fi
done

TASK_LINE="$(awk -F '\\t' -v task="$SLURM_ARRAY_TASK_ID" 'NR == task + 2 {{print; exit}}' "$TASK_TABLE")"
if [[ -z "$TASK_LINE" ]]; then
    echo "ERROR: no task row for SLURM_ARRAY_TASK_ID=$SLURM_ARRAY_TASK_ID" >&2
    exit 1
fi
IFS=$'\\t' read -r library library_number fragments rna_matrix_dir library_output_dir <<< "$TASK_LINE"

for required in "$fragments" "$rna_matrix_dir/barcodes.tsv.gz" \\
    "$rna_matrix_dir/features.tsv.gz" "$rna_matrix_dir/matrix.mtx.gz"; do
    if [[ ! -s "$required" ]]; then
        echo "ERROR: required library input missing or empty: $required" >&2
        exit 1
    fi
done

task_tmp="${{SLURM_TMPDIR:-$library_output_dir/tmp}}"
mkdir -p "$task_tmp" "$library_output_dir"

python "$WORKER" import-library \\
    --library "$library" \\
    --fragments "$fragments" \\
    --rna-matrix-dir "$rna_matrix_dir" \\
    --chrom-sizes "$CHROM_SIZES" \\
    --gtf "$GTF" \\
    --output-dir "$library_output_dir" \\
    --temp-dir "$task_tmp" \\
    --cpus "$SLURM_CPUS_PER_TASK" \\
    --min-fragments {args.min_fragments} \\
    --min-tsse {args.min_tsse} \\
    --tile-size {args.tile_size} \\
    --exclude-chroms {quote(exclusions)}{force}
"""


def render_atac_script(args: argparse.Namespace, paths: dict[str, Path]) -> str:
    force = " --force" if args.force else ""
    return sbatch_header(
        job_name="tet_atac_cluster",
        log_pattern=paths["log_dir"] / "tet_atac_cluster_%j",
        partition=args.partition,
        cpus=args.atac_cpus,
        memory=args.atac_memory,
        duration=args.time,
        working_directory=paths["analysis_root"],
    ) + f"""
for required in {quote(paths['task_table'])} {quote(paths['worker'])}; do
    if [[ ! -s "$required" ]]; then
        echo "ERROR: required input missing or empty: $required" >&2
        exit 1
    fi
done
mkdir -p {quote(paths['analysis_root'])} {quote(paths['figure_root'])}

python {quote(paths['worker'])} cluster-atac \\
    --task-table {quote(paths['task_table'])} \\
    --analysis-root {quote(paths['analysis_root'])} \\
    --figure-root {quote(paths['figure_root'])} \\
    --identity-table {quote(args.identity_table)} \\
    --workbook {quote(args.workbook)} \\
    --batch-key {quote(args.batch_key)} \\
    --cpus "$SLURM_CPUS_PER_TASK" \\
    --n-features {args.atac_features} \\
    --n-components {args.components} \\
    --neighbors {args.neighbors} \\
    --leiden-resolution {args.leiden_resolution} \\
    --max-plot-cells {args.max_plot_cells}{force}
"""


def render_joint_script(args: argparse.Namespace, paths: dict[str, Path]) -> str:
    force = " --force" if args.force else ""
    return sbatch_header(
        job_name="tet_joint_embed",
        log_pattern=paths["log_dir"] / "tet_joint_embed_%j",
        partition=args.partition,
        cpus=args.joint_cpus,
        memory=args.joint_memory,
        duration=args.time,
        working_directory=paths["analysis_root"],
    ) + f"""
for required in {quote(paths['task_table'])} {quote(paths['worker'])} \\
    {quote(paths['analysis_root'] / 'atac' / 'atac_dataset.h5ads')}; do
    if [[ ! -s "$required" ]]; then
        echo "ERROR: required input missing or empty: $required" >&2
        exit 1
    fi
done
mkdir -p {quote(paths['analysis_root'])} {quote(paths['figure_root'])}

python {quote(paths['worker'])} joint-embed \\
    --task-table {quote(paths['task_table'])} \\
    --analysis-root {quote(paths['analysis_root'])} \\
    --figure-root {quote(paths['figure_root'])} \\
    --batch-key {quote(args.batch_key)} \\
    --cpus "$SLURM_CPUS_PER_TASK" \\
    --rna-features {args.rna_features} \\
    --n-components {args.components} \\
    --neighbors {args.neighbors} \\
    --leiden-resolution {args.leiden_resolution} \\
    --max-plot-cells {args.max_plot_cells}{force}
"""


def validate_bash(path: Path) -> None:
    result = subprocess.run(
        ["bash", "-n", str(path)],
        text=True,
        capture_output=True,
        check=False,
    )
    if result.returncode != 0:
        raise RuntimeError(f"bash -n failed for {path}:\n{result.stderr}")


def submit_job(path: Path, dependency: str | None = None) -> str:
    command = ["sbatch", "--parsable"]
    if dependency:
        command.append(f"--dependency=afterok:{dependency}")
    command.append(str(path))
    result = subprocess.run(command, text=True, capture_output=True, check=False)
    if result.returncode != 0:
        raise RuntimeError(
            f"submission failed for {path}: {result.stderr.strip() or result.stdout.strip()}"
        )
    job_id = result.stdout.strip().split(";", 1)[0]
    if not re.fullmatch(r"\d+(?:_\d+)?", job_id):
        raise RuntimeError(f"could not parse sbatch job id from: {result.stdout!r}")
    return job_id


def positive_int(value: str) -> int:
    parsed = int(value)
    if parsed <= 0:
        raise argparse.ArgumentTypeError("must be positive")
    return parsed


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Discover matched Tet 2025 libraries, render a SnapATAC2 import array, "
            "ATAC clustering job, and RNA-ATAC joint embedding job."
        )
    )
    parser.add_argument("--version", action="version", version=VERSION)
    parser.add_argument("--atac-mapping-root", type=Path, default=ATAC_MAPPING_ROOT)
    parser.add_argument("--rna-mapping-root", type=Path, default=RNA_MAPPING_ROOT)
    parser.add_argument("--run-label", default=RUN_LABEL)
    parser.add_argument("--analysis-root", type=Path, default=None)
    parser.add_argument("--figure-root", type=Path, default=None)
    parser.add_argument("--identity-table", type=Path, default=IDENTITY_TABLE)
    parser.add_argument("--workbook", type=Path, default=WORKBOOK)
    parser.add_argument("--fasta", type=Path, default=REFERENCE_FASTA)
    parser.add_argument("--chrom-sizes", type=Path, default=None)
    parser.add_argument("--gtf", type=Path, default=REFERENCE_GTF)
    parser.add_argument("--worker", type=Path, default=Path(__file__).with_name("joint_atac_rna.py"))
    parser.add_argument("--libraries", type=positive_int, nargs="+", default=None)
    parser.add_argument(
        "--stage",
        choices=(
            "ALL",
            "IMPORT",
            "ATAC_CLUSTER",
            "ATAC_CLUSTER_AND_JOINT",
            "JOINT",
        ),
        default="ALL",
        help=(
            "Workflow segment to submit. ATAC_CLUSTER_AND_JOINT resumes from "
            "completed per-library imports and chains clustering into the "
            "joint embedding."
        ),
    )
    parser.add_argument("--submit", action="store_true")
    parser.add_argument("--force", action="store_true")

    parser.add_argument("--partition", default="compute")
    parser.add_argument("--time", default="7-00:00:00")
    parser.add_argument("--array-concurrency", type=positive_int, default=12)
    parser.add_argument("--import-cpus", type=positive_int, default=8)
    parser.add_argument("--import-memory", default="96G")
    parser.add_argument("--atac-cpus", type=positive_int, default=32)
    parser.add_argument("--atac-memory", default="384G")
    parser.add_argument("--joint-cpus", type=positive_int, default=48)
    parser.add_argument("--joint-memory", default="768G")

    parser.add_argument("--min-fragments", type=positive_int, default=1000)
    parser.add_argument("--min-tsse", type=float, default=5.0)
    parser.add_argument("--tile-size", type=positive_int, default=5000)
    parser.add_argument("--atac-features", type=positive_int, default=50000)
    parser.add_argument("--rna-features", type=positive_int, default=3000)
    parser.add_argument("--components", type=positive_int, default=30)
    parser.add_argument("--neighbors", type=positive_int, default=50)
    parser.add_argument("--leiden-resolution", type=float, default=1.0)
    parser.add_argument("--max-plot-cells", type=positive_int, default=200000)
    parser.add_argument(
        "--batch-key",
        choices=("library", "diff_batch", "none"),
        default="library",
        help="Harmony covariate; raw embeddings are always retained.",
    )
    parser.add_argument(
        "--exclude-chrom",
        action="append",
        default=["chrM", "M", "MT", "chrY", "Y", "human_chimp_bonoborefChr1218"],
        help="Chromosome/scaffold excluded from tile matrices; repeatable.",
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.min_tsse < 0:
        raise SystemExit("--min-tsse must be non-negative")
    if args.leiden_resolution <= 0:
        raise SystemExit("--leiden-resolution must be positive")

    # Resolve user-supplied paths while still on the login node so rendered
    # batch jobs never depend on SLURM's eventual working directory.
    args.identity_table = args.identity_table.expanduser().resolve()
    args.workbook = args.workbook.expanduser().resolve()
    args.fasta = args.fasta.expanduser().resolve()
    args.gtf = args.gtf.expanduser().resolve()

    atac_root = args.atac_mapping_root.expanduser().resolve()
    rna_root = args.rna_mapping_root.expanduser().resolve()
    analysis_root = (
        args.analysis_root.expanduser().resolve()
        if args.analysis_root
        else (ATAC_ANALYSIS_BASE / args.run_label).resolve()
    )
    figure_root = (
        args.figure_root.expanduser().resolve()
        if args.figure_root
        else (ATAC_FIGURE_BASE / args.run_label).resolve()
    )
    worker = args.worker.expanduser().resolve()
    if not worker.is_file():
        raise SystemExit(f"worker script is missing: {worker}")

    rows = discover_libraries(atac_root, rna_root, args.libraries)
    if args.stage in {"ALL", "IMPORT"}:
        chrom_sizes = resolve_chrom_sizes(
            args.chrom_sizes.expanduser().resolve() if args.chrom_sizes else None,
            args.fasta,
        )
    else:
        chrom_sizes = (
            args.chrom_sizes.expanduser().resolve()
            if args.chrom_sizes is not None
            else Path(str(args.fasta) + ".fai")
        )

    control_root = analysis_root / "control"
    script_dir = control_root / "slurm_scripts"
    log_dir = control_root / "logs"
    script_dir.mkdir(parents=True, exist_ok=True)
    log_dir.mkdir(parents=True, exist_ok=True)
    figure_root.mkdir(parents=True, exist_ok=True)
    task_table = write_task_table(control_root / "library_tasks.tsv", rows, analysis_root)

    paths = {
        "analysis_root": analysis_root,
        "figure_root": figure_root,
        "script_dir": script_dir,
        "log_dir": log_dir,
        "task_table": task_table.resolve(),
        "worker": worker,
        "chrom_sizes": chrom_sizes,
    }

    scripts: dict[str, Path] = {}
    if args.stage in {"ALL", "IMPORT"}:
        scripts["IMPORT"] = write_if_changed(
            script_dir / "01_import_libraries.sbatch",
            render_import_script(args, paths, len(rows)),
            executable=True,
        )
    if args.stage in {"ALL", "ATAC_CLUSTER", "ATAC_CLUSTER_AND_JOINT"}:
        scripts["ATAC_CLUSTER"] = write_if_changed(
            script_dir / "02_cluster_atac.sbatch",
            render_atac_script(args, paths),
            executable=True,
        )
    if args.stage in {"ALL", "ATAC_CLUSTER_AND_JOINT", "JOINT"}:
        scripts["JOINT"] = write_if_changed(
            script_dir / "03_joint_embedding.sbatch",
            render_joint_script(args, paths),
            executable=True,
        )
    for script in scripts.values():
        validate_bash(script)

    print(f"Discovered paired libraries: {len(rows)}")
    print("Libraries: " + ", ".join(str(row["library_number"]) for row in rows))
    print(f"Analysis root: {analysis_root}")
    print(f"Figure root: {figure_root}")
    print(f"Task table: {task_table}")
    for stage, script in scripts.items():
        print(f"Rendered {stage}: {script}")

    if not args.submit:
        print("Jobs were rendered and syntax-checked; add --submit to submit them.")
        return 0

    job_ids: dict[str, str] = {}
    dependency: str | None = None
    for stage in ("IMPORT", "ATAC_CLUSTER", "JOINT"):
        script = scripts.get(stage)
        if script is None:
            continue
        job_id = submit_job(script, dependency=dependency)
        job_ids[stage] = job_id
        dependency = job_id
        print(f"Submitted {stage}: job {job_id}")
    if job_ids:
        print("Dependency chain: " + " -> ".join(f"{stage}={job}" for stage, job in job_ids.items()))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (FileNotFoundError, RuntimeError, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1)
