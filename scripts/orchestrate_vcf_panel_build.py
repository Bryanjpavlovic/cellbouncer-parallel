#!/usr/bin/env python3
"""Render and optionally submit RNA/ATAC coverage and VCF-panel jobs.

This is deliberately a small, standard-library-only orchestrator.  The heavy
work is performed by ``bam_window_coverage`` and ``downsample_vcf_parallel`` in
SLURM jobs; this file only resolves the fixed task graph, writes concrete task
lists and sbatch files, syntax-checks them, and optionally submits them.
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import gzip
import hashlib
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys
from typing import Iterable, Sequence


VERSION = "1.9.0-central-figures"

STAGES = (
    "ALL",
    "INCREMENTAL_ALL",
    "COVERAGE_ONLY",
    "AGGREGATE_ONLY",
    "ATAC_CANDIDATES_ONLY",
    "DOWNSAMPLE_ONLY",
    "DOWNSAMPLE_AND_PLOTS",
    "ATAC_RESELECT_AND_PLOTS",
    "PLOTS_ONLY",
)
MODALITIES = ("rna", "atac")
POPULATIONS = ("cells", "all_barcoded", "noncell", "empty", "bam_all")
BASE_MERGE_POPULATIONS = ("bam_all", "all_barcoded", "cells", "noncell")

CURRENT_RNA_RUN_ROOT = (
    "/mnt/beegfs/tetraploid_multiome_cis_trans/3P"
)
CURRENT_RNA_ANALYSIS_ROOT = CURRENT_RNA_RUN_ROOT + "/analysis"
DEFAULT_RNA_ROOT = CURRENT_RNA_RUN_ROOT + "/mapping_output"
DEFAULT_ATAC_ROOT = "/mnt/beegfs/tetraploid_multiome_cis_trans/ATAC/mapping_output"
DEFAULT_OUTPUT_ROOT = (
    CURRENT_RNA_ANALYSIS_ROOT + "/aggregate_library_analysis/vcf_panel_build"
)
DEFAULT_FIGURE_ROOT = CURRENT_RNA_RUN_ROOT + "/figures/all40/vcf_panels"
DEFAULT_SOURCE_BCF = "/nvme/software/shared_data/downsamples/tet.vars.all.use.bcf"
DEFAULT_SOURCE_CSI = DEFAULT_SOURCE_BCF + ".csi"
DEFAULT_GTF = (
    "/mnt/beegfs/genomes_annotations/ancestral_genomes/litterbox/"
    "human_chimp_bonobo/human_chimp_bonobo.gtf.gz"
)
DEFAULT_NUMT_BED = (
    "/mnt/beegfs/genomes_annotations/ancestral_genomes/litterbox/"
    "human_chimp_bonobo/numts.bed"
)
DEFAULT_PANEL_METADATA = (
    "/mnt/beegfs/tetraploid_multiome_cis_trans/Misc_Metadata/panel_metadata.tsv"
)
DEFAULT_POOL_COMBINATIONS = (
    "/mnt/beegfs/tetraploid_multiome_cis_trans/Misc_Metadata/pool_combinations.tsv"
)
DEFAULT_PRODUCTION_RNA_PANEL_ROOT = (
    "/mnt/beegfs/home/b/vcfdownsample/Downsample_ATAC_Species_poolInformer/NoMito"
)
DEFAULT_RNA_DEMUX_PANEL = (
    DEFAULT_PRODUCTION_RNA_PANEL_ROOT + "/tet.vars.downsampled_20M.bcf"
)
DEFAULT_RNA_HET_PANEL = (
    DEFAULT_PRODUCTION_RNA_PANEL_ROOT + "/tet.vars.het_10M.bcf"
)
DEFAULT_RNA_SPECIES_PANEL = (
    DEFAULT_PRODUCTION_RNA_PANEL_ROOT + "/tet.vars.species_20M.bcf"
)

DEFAULT_RNA_BAM_TEMPLATE = (
    "{rna_root}/Tet_2025_Multiome-RNA_{lib}/gex.bam"
)
DEFAULT_ATAC_BAM_TEMPLATE = (
    "{atac_root}/Tet_2025_Multiome-ATAC_{lib}/atac.bam"
)
DEFAULT_CELL_BARCODES_TEMPLATE = (
    "{rna_root}/Tet_2025_Multiome-RNA_{lib}/filtered/barcodes.tsv.gz"
)

DEFAULT_COVERAGE_BINARY = "/nvme/software/packages/cellbouncer/dev/bin/bam_window_coverage"
DEFAULT_DOWNSAMPLE_BINARY = (
    "/nvme/software/packages/cellbouncer/dev/bin/downsample_vcf_parallel"
)
DEFAULT_PANEL_QC_SCRIPT = (
    "/nvme/software/packages/cellbouncer/dev/plot/vcf_panel_qc.py"
)

TASK_COLUMNS = (
    "task_id",
    "library_id",
    "modality",
    "bam",
    "cell_barcodes",
    "empty_barcodes",
    "coverage_output",
    "metadata_output",
)
EXPECTED_COVERAGE_HEADER = (
    "#contig\tstart\tend\tbam_all\tall_barcoded\tcells\tnoncell\tempty"
)


class PipelineError(RuntimeError):
    """An operator-facing configuration or submission error."""


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    try:
        with path.open("rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(chunk)
    except OSError as exc:
        raise PipelineError(f"cannot hash blacklist {path}: {exc}") from exc
    return digest.hexdigest()


def validate_blacklist_syntax(path: Path) -> None:
    try:
        with path.open("r", encoding="utf-8") as handle:
            for line_number, raw_line in enumerate(handle, start=1):
                line = raw_line.strip()
                if not line or line.startswith("#"):
                    continue
                fields = line.split()
                if len(fields) == 1:
                    continue
                if len(fields) == 2:
                    raise PipelineError(
                        f"malformed blacklist line {line_number} in {path}: "
                        "two-column records are not allowed"
                    )
                try:
                    start = int(fields[1], 10)
                    end = int(fields[2], 10)
                except ValueError as exc:
                    raise PipelineError(
                        f"invalid blacklist coordinates at line {line_number} "
                        f"in {path}: {fields[1]!r}, {fields[2]!r}"
                    ) from exc
                if start < 0 or end <= start:
                    raise PipelineError(
                        f"invalid 0-based half-open blacklist interval at line "
                        f"{line_number} in {path}: {start}, {end}"
                    )
    except (OSError, UnicodeError) as exc:
        raise PipelineError(f"cannot read blacklist {path}: {exc}") from exc


def validate_and_stage_blacklist(args: argparse.Namespace, run_dir: Path) -> None:
    args.blacklist_source_path = None
    args.blacklist_staged_path = None
    args.blacklist_size_bytes = None
    args.blacklist_mtime_ns = None
    args.blacklist_sha256 = None
    if not args.blacklist:
        return

    source = Path(args.blacklist)
    require_nonempty(source, "blacklist")
    validate_blacklist_syntax(source)
    source_stat = source.stat()
    source_sha256 = sha256_file(source)

    staged = run_dir / "control" / "inputs" / "blacklist.tsv"
    staged.parent.mkdir(parents=True, exist_ok=True)
    try:
        if source.resolve() != staged.resolve():
            temporary = staged.with_name(
                f"{staged.name}.partial.{os.getpid()}"
            )
            shutil.copy2(source, temporary)
            os.replace(temporary, staged)
    except OSError as exc:
        raise PipelineError(
            f"cannot stage blacklist {source} as {staged}: {exc}"
        ) from exc
    require_nonempty(staged, "staged blacklist")
    if sha256_file(staged) != source_sha256:
        raise PipelineError(
            f"staged blacklist checksum differs from its source: {staged}"
        )

    args.blacklist_source_path = str(source)
    args.blacklist_staged_path = str(staged)
    args.blacklist_size_bytes = source_stat.st_size
    args.blacklist_mtime_ns = source_stat.st_mtime_ns
    args.blacklist_sha256 = source_sha256


def append_blacklist_contract(
    rows: list[tuple[object, object]], args: argparse.Namespace
) -> None:
    if not args.blacklist_staged_path:
        return
    rows.extend(
        (
            ("blacklist_source_path", args.blacklist_source_path),
            ("blacklist_staged_path", args.blacklist_staged_path),
            ("blacklist_size_bytes", args.blacklist_size_bytes),
            ("blacklist_mtime_ns", args.blacklist_mtime_ns),
            ("blacklist_sha256", args.blacklist_sha256),
        )
    )


def shell_join(parts: Iterable[object]) -> str:
    return " ".join(shlex.quote(str(part)) for part in parts)


def positive_int(text: str) -> int:
    try:
        value = int(text)
    except ValueError as exc:
        raise argparse.ArgumentTypeError(f"expected an integer, got {text!r}") from exc
    if value <= 0:
        raise argparse.ArgumentTypeError("value must be greater than zero")
    return value


def nonnegative_int(text: str) -> int:
    try:
        value = int(text)
    except ValueError as exc:
        raise argparse.ArgumentTypeError(f"expected an integer, got {text!r}") from exc
    if value < 0:
        raise argparse.ArgumentTypeError("value must be zero or greater")
    return value


def auto_base_int(text: str) -> int:
    try:
        value = int(text, 16 if text.lower().startswith("0x") else 10)
    except ValueError as exc:
        raise argparse.ArgumentTypeError(
            f"expected an integer (decimal or 0x-prefixed), got {text!r}"
        ) from exc
    if value < 0:
        raise argparse.ArgumentTypeError("value must be zero or greater")
    return value


def nonnegative_float(text: str) -> float:
    try:
        value = float(text)
    except ValueError as exc:
        raise argparse.ArgumentTypeError(f"expected a number, got {text!r}") from exc
    if value < 0:
        raise argparse.ArgumentTypeError("value must be zero or greater")
    return value


def unit_interval_float(text: str) -> float:
    value = nonnegative_float(text)
    if value > 1:
        raise argparse.ArgumentTypeError("value must be between zero and one")
    return value


def parse_modalities(text: str) -> tuple[str, ...]:
    values = [part.strip().lower() for part in text.split(",") if part.strip()]
    if not values:
        raise PipelineError("--modalities must select rna, atac, or both")
    unknown = sorted(set(values) - set(MODALITIES))
    if unknown:
        raise PipelineError(
            "unsupported modality value(s): " + ", ".join(unknown)
        )
    return tuple(modality for modality in MODALITIES if modality in set(values))


def parse_libraries(tokens: Sequence[str]) -> tuple[int, ...]:
    values: set[int] = set()
    for raw_token in tokens:
        for token in raw_token.split(","):
            token = token.strip()
            if not token:
                continue
            match = re.fullmatch(r"(\d+)-(\d+)", token)
            if match:
                start, end = (int(match.group(1)), int(match.group(2)))
                if start <= 0 or end <= 0 or end < start:
                    raise PipelineError(f"invalid library range: {token!r}")
                values.update(range(start, end + 1))
                continue
            if not token.isdigit() or int(token) <= 0:
                raise PipelineError(f"invalid library identifier: {token!r}")
            values.add(int(token))
    if not values:
        raise PipelineError("--libraries resolved to an empty set")
    return tuple(sorted(values))


def ensure_absolute(path_text: str, label: str) -> Path:
    path = Path(path_text)
    if not path.is_absolute():
        raise PipelineError(f"{label} must be an absolute path: {path_text}")
    if "\n" in path_text or "\t" in path_text:
        raise PipelineError(f"{label} contains a newline or tab")
    return path


def require_nonempty(path: Path, label: str) -> None:
    try:
        size = path.stat().st_size
    except FileNotFoundError as exc:
        raise PipelineError(f"missing {label}: {path}") from exc
    if not path.is_file() or size <= 0:
        raise PipelineError(f"{label} is not a nonempty regular file: {path}")


def require_tabix_readable(path: Path, label: str) -> None:
    tabix = shutil.which("tabix")
    if not tabix:
        raise PipelineError(
            f"cannot validate {label} because tabix is not available"
        )
    try:
        completed = subprocess.run(
            [tabix, "-l", str(path)],
            check=False,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
    except OSError as exc:
        raise PipelineError(f"cannot inspect {label} {path}: {exc}") from exc
    if completed.returncode != 0 or not any(
        line.strip() for line in completed.stdout.splitlines()
    ):
        detail = completed.stderr.strip()
        suffix = f": {detail}" if detail else ""
        raise PipelineError(f"{label} is not Tabix-readable: {path}{suffix}")


def require_executable(path: Path, label: str) -> None:
    require_nonempty(path, label)
    if not os.access(path, os.X_OK):
        raise PipelineError(f"{label} is not executable: {path}")


def render_template(template: str, lib: int, args: argparse.Namespace, label: str) -> Path:
    try:
        rendered = template.format(
            lib=lib,
            rna_root=args.rna_root.rstrip("/"),
            atac_root=args.atac_root.rstrip("/"),
        )
    except (KeyError, IndexError, ValueError) as exc:
        raise PipelineError(
            f"invalid {label}; only {{lib}}, {{rna_root}}, and {{atac_root}} are allowed: {exc}"
        ) from exc
    return ensure_absolute(rendered, label)


def write_text(path: Path, text: str, executable: bool = False) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not text.endswith("\n"):
        text += "\n"
    path.write_text(text, encoding="utf-8")
    if executable:
        path.chmod(0o755)


def write_tsv(path: Path, columns: Sequence[str], rows: Iterable[Sequence[object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(columns)
        for row in rows:
            writer.writerow(row)


def bash_syntax_check(path: Path) -> None:
    result = subprocess.run(
        ["/bin/bash", "-n", str(path)],
        text=True,
        capture_output=True,
        check=False,
    )
    if result.returncode:
        detail = result.stderr.strip() or result.stdout.strip() or "unknown syntax error"
        raise PipelineError(f"bash -n failed for {path}: {detail}")


def sbatch_header(
    *,
    job_name: str,
    partition: str,
    threads: int,
    memory: str,
    walltime: str,
    stdout_path: Path,
    stderr_path: Path,
    array_spec: str | None = None,
) -> str:
    lines = [
        "#!/bin/bash",
        f"#SBATCH --job-name={job_name[:120]}",
        f"#SBATCH --partition={partition}",
        "#SBATCH --nodes=1",
        "#SBATCH --ntasks=1",
        f"#SBATCH --cpus-per-task={threads}",
        f"#SBATCH --mem={memory}",
        f"#SBATCH --time={walltime}",
        f"#SBATCH --output={stdout_path}",
        f"#SBATCH --error={stderr_path}",
        "#SBATCH --chdir=/tmp",
    ]
    if array_spec is not None:
        lines.append(f"#SBATCH --array={array_spec}")
    return "\n".join(lines)


def module_block(
    *,
    include_bcftools: bool = False,
    include_genomics_python: bool = False,
) -> str:
    loads = [
        (
            "module load miniforge/3 genomics-base/latest"
            if include_genomics_python
            else "module load miniforge/3"
        ),
        "module load htslib/1.20",
        "module load samtools/1.20",
    ]
    if include_bcftools:
        loads.append("module load bcftools/1.20")
    return "\n".join(
        [
            'command -v module >/dev/null 2>&1 || { echo "ERROR: module command unavailable" >&2; exit 1; }',
            "module purge",
            *loads,
            "module list",
        ]
    )


def make_task_rows(
    args: argparse.Namespace,
    run_dir: Path,
    libraries: Sequence[int],
    modalities: Sequence[str],
    empty_barcode_paths: dict[int, Path],
) -> list[tuple[object, ...]]:
    rows: list[tuple[object, ...]] = []
    task_id = 0
    for lib in libraries:
        cell_barcodes = render_template(
            args.cell_barcodes_template, lib, args, "--cell-barcodes-template"
        )
        empty_barcodes = "NA"
        if empty_barcode_paths:
            empty_barcodes = str(empty_barcode_paths[lib])
        elif args.empty_barcodes_template:
            empty_barcodes = str(
                render_template(
                    args.empty_barcodes_template,
                    lib,
                    args,
                    "--empty-barcodes-template",
                )
            )
        for modality in modalities:
            template = args.rna_bam_template if modality == "rna" else args.atac_bam_template
            bam = render_template(template, lib, args, f"--{modality}-bam-template")
            per_library = run_dir / "coverage" / "per_library" / modality
            output = per_library / f"library_{lib}.{modality}.coverage.tsv.gz"
            metadata = per_library / f"library_{lib}.{modality}.coverage.metadata.tsv"
            rows.append(
                (
                    task_id,
                    lib,
                    modality,
                    str(bam),
                    str(cell_barcodes),
                    empty_barcodes,
                    str(output),
                    str(metadata),
                )
            )
            task_id += 1
    return rows


def has_empty_barcode_source(args: argparse.Namespace) -> bool:
    return bool(args.empty_barcodes_template or args.empty_drop_roster)


def parse_roster_library(value: str) -> int | None:
    value = value.strip()
    for pattern in (
        r"(?:lib)?(\d+)",
        r"Tet_2025_Multiome-(?:RNA|ATAC)_(\d+)",
    ):
        match = re.fullmatch(pattern, value, flags=re.IGNORECASE)
        if match:
            return int(match.group(1))
    return None


def normalize_10x_barcode(value: str) -> str:
    value = value.strip().split()[0] if value.strip() else ""
    if value in {"", "-", "*", "."}:
        return ""
    return re.sub(r"-\d+$", "", value)


def split_empty_drop_roster(
    args: argparse.Namespace, run_dir: Path, libraries: Sequence[int]
) -> dict[int, Path]:
    if not args.empty_drop_roster:
        return {}
    source = Path(args.empty_drop_roster)
    require_nonempty(source, "empty-drop roster")
    opener = gzip.open if source.name.endswith(".gz") else open
    selected = set(libraries)
    barcodes: dict[int, set[str]] = {lib: set() for lib in libraries}
    try:
        with opener(source, "rt", encoding="utf-8") as handle:
            header = handle.readline().rstrip("\r\n").split("\t")
            normalized = {name.strip().lower(): index for index, name in enumerate(header)}
            library_index = normalized.get("library")
            barcode_index = next(
                (
                    normalized[name]
                    for name in ("cell_barcode", "barcode", "cb")
                    if name in normalized
                ),
                None,
            )
            if library_index is None or barcode_index is None:
                raise PipelineError(
                    "--empty-drop-roster must be tab-separated with a library column and "
                    "one of cell_barcode, barcode, or CB"
                )
            for line_number, line in enumerate(handle, start=2):
                if not line.strip() or line.startswith("#"):
                    continue
                fields = line.rstrip("\r\n").split("\t")
                if max(library_index, barcode_index) >= len(fields):
                    raise PipelineError(
                        f"empty-drop roster row {line_number} is incomplete"
                    )
                library = parse_roster_library(fields[library_index])
                barcode = normalize_10x_barcode(fields[barcode_index])
                if library is None:
                    raise PipelineError(
                        f"empty-drop roster row {line_number} has an unsupported library "
                        f"identifier: {fields[library_index]!r}"
                    )
                if not barcode:
                    raise PipelineError(
                        f"empty-drop roster row {line_number} has an empty barcode"
                    )
                if library in selected:
                    barcodes[library].add(barcode)
    except (OSError, UnicodeError) as exc:
        raise PipelineError(f"cannot read empty-drop roster {source}: {exc}") from exc

    output_dir = run_dir / "control" / "tasks" / "empty_barcodes"
    output_dir.mkdir(parents=True, exist_ok=True)
    result: dict[int, Path] = {}
    for lib in libraries:
        if not barcodes[lib]:
            raise PipelineError(
                f"empty-drop roster contains no barcodes for selected library {lib}"
            )
        path = output_dir / f"library_{lib}.empty_barcodes.tsv.gz"
        with gzip.open(path, "wt", encoding="utf-8", newline="") as handle:
            for barcode in sorted(barcodes[lib]):
                handle.write(barcode + "\n")
        result[lib] = path
    return result


def validate_coverage_inputs(rows: Sequence[Sequence[object]]) -> None:
    checked: set[tuple[str, str]] = set()
    for row in rows:
        values = dict(zip(TASK_COLUMNS, row))
        for key, label in (
            ("bam", f"library {values['library_id']} {values['modality']} BAM"),
            ("cell_barcodes", f"library {values['library_id']} cell-barcode list"),
        ):
            path = Path(str(values[key]))
            marker = (key, str(path))
            if marker not in checked:
                require_nonempty(path, label)
                checked.add(marker)
        if values["empty_barcodes"] != "NA":
            path = Path(str(values["empty_barcodes"]))
            marker = ("empty_barcodes", str(path))
            if marker not in checked:
                require_nonempty(path, f"library {values['library_id']} empty-barcode list")
                checked.add(marker)


def render_coverage_script(
    args: argparse.Namespace,
    run_dir: Path,
    task_file: Path,
    task_count: int,
) -> Path:
    if task_count <= 0:
        raise PipelineError("cannot render an empty coverage array")
    concurrency = f"%{args.coverage_concurrency}" if args.coverage_concurrency else ""
    array_spec = f"0-{task_count - 1}{concurrency}"
    script = run_dir / "control" / "slurm" / "coverage_array.sbatch"
    logs = run_dir / "logs"
    force_line = 'BUILD_ARGS+=(--force)' if args.force else ":"
    normalize_line = 'BUILD_ARGS+=(--no-normalize-10x)' if args.no_normalize_10x else ":"
    header = sbatch_header(
        job_name=f"vcov_{args.run_label}",
        partition=args.partition,
        threads=args.coverage_threads,
        memory=args.coverage_memory,
        walltime=args.coverage_time,
        stdout_path=logs / "coverage_%A_%a.out",
        stderr_path=logs / "coverage_%A_%a.err",
        array_spec=array_spec,
    )
    body = f"""{header}

set -euo pipefail

{module_block()}

for executable in samtools python3; do
    command -v "$executable" >/dev/null 2>&1 || {{ echo "ERROR: missing command: $executable" >&2; exit 1; }}
done
[[ -x {shlex.quote(args.bam_window_coverage)} ]] || {{ echo "ERROR: coverage binary is not executable: {args.bam_window_coverage}" >&2; exit 1; }}

TASK_FILE={shlex.quote(str(task_file))}
[[ -s "$TASK_FILE" ]] || {{ echo "ERROR: missing task list: $TASK_FILE" >&2; exit 1; }}
ROW="$(awk -F '\t' -v task="$SLURM_ARRAY_TASK_ID" 'NR == task + 2 {{print; found=1}} END {{if (!found) exit 1}}' "$TASK_FILE")"
IFS=$'\t' read -r TASK_ID LIBRARY_ID MODALITY BAM CELL_BARCODES EMPTY_BARCODES OUTPUT METADATA <<< "$ROW"
[[ "$TASK_ID" == "$SLURM_ARRAY_TASK_ID" ]] || {{ echo "ERROR: task-list index mismatch" >&2; exit 1; }}
for value in "$LIBRARY_ID" "$MODALITY" "$BAM" "$CELL_BARCODES" "$OUTPUT" "$METADATA"; do
    [[ -n "$value" ]] || {{ echo "ERROR: malformed task-list row: $ROW" >&2; exit 1; }}
done

echo "Coverage task: $TASK_ID"
echo "Library: $LIBRARY_ID"
echo "Modality: $MODALITY"
echo "BAM: $BAM"
echo "Cell barcodes (RNA namespace): $CELL_BARCODES"
echo "Empty barcodes: $EMPTY_BARCODES"
echo "Output: $OUTPUT"
if [[ "$MODALITY" == "rna" ]]; then
    MIN_MAPQ={args.rna_min_mapq}
elif [[ "$MODALITY" == "atac" ]]; then
    MIN_MAPQ={args.atac_min_mapq}
else
    echo "ERROR: unsupported modality in task list: $MODALITY" >&2
    exit 1
fi
echo "Minimum MAPQ: $MIN_MAPQ"

[[ -s "$BAM" ]] || {{ echo "ERROR: BAM is missing or empty: $BAM" >&2; exit 1; }}
[[ -s "$CELL_BARCODES" ]] || {{ echo "ERROR: cell-barcode list is missing or empty: $CELL_BARCODES" >&2; exit 1; }}
if [[ "$EMPTY_BARCODES" != "NA" ]]; then
    [[ -s "$EMPTY_BARCODES" ]] || {{ echo "ERROR: empty-barcode list is missing or empty: $EMPTY_BARCODES" >&2; exit 1; }}
fi

samtools quickcheck -v "$BAM"
samtools view -H "$BAM" | awk -F '\t' '
    $1 == "@HD" {{for (i=2; i<=NF; ++i) if ($i == "SO:coordinate") found=1}}
    END {{if (!found) {{print "ERROR: BAM header is not SO:coordinate" > "/dev/stderr"; exit 1}}}}
'
samtools idxstats "$BAM" >/dev/null

set +o pipefail
CB_SAMPLE="$(samtools view -@ 1 "$BAM" | head -n 50000 | awk '
    {{for (i=12; i<=NF; ++i) if ($i ~ /^CB:Z:/) {{print substr($i,6); exit}}}}
')"
set -o pipefail
if [[ -n "$CB_SAMPLE" ]]; then
    echo "Sample CB tag: $CB_SAMPLE"
else
    echo "WARNING: no CB:Z tag occurred in the first 50,000 alignments; the full worker scan will enforce barcode matches" >&2
fi

python3 - "$CELL_BARCODES" "$EMPTY_BARCODES" <<'PY'
import gzip
from pathlib import Path
import sys

for raw in sys.argv[1:]:
    if raw == "NA":
        continue
    path = Path(raw)
    opener = gzip.open if path.suffix == ".gz" else open
    found = False
    with opener(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.strip():
                found = True
                break
    if not found:
        raise SystemExit(f"ERROR: barcode list has no nonempty records: {{path}}")
PY

EXPECTED_EMPTY="$EMPTY_BARCODES"
[[ "$EXPECTED_EMPTY" == "NA" ]] && EXPECTED_EMPTY="NONE"
validate_task_metadata() {{
    python3 - "$METADATA" "$BAM" "$CELL_BARCODES" "$EXPECTED_EMPTY" "$OUTPUT" "$MIN_MAPQ" <<'PY'
import csv
import sys

path, bam, cells, empty, output, min_mapq = sys.argv[1:]
expected = {{
    "schema": "bam_window_coverage_v1",
    "command": "build",
    "status": "complete",
    "bam": bam,
    "cell_barcodes": cells,
    "empty_barcodes": empty,
    "coverage_output": output,
    "window_size": {str(args.coverage_window_size)!r},
    "min_mapq": min_mapq,
    "exclude_flags": {f'0x{args.coverage_exclude_flags:X}'!r},
    "barcode_tag": {args.barcode_tag!r},
    "normalize_10x": {('false' if args.no_normalize_10x else 'true')!r},
}}
with open(path, "r", encoding="utf-8", newline="") as handle:
    reader = csv.reader(handle, delimiter="\t")
    header = next(reader, None)
    if header != ["key", "value"]:
        raise SystemExit(f"ERROR: unexpected metadata header in {{path}}: {{header!r}}")
    observed = {{}}
    for line_number, row in enumerate(reader, 2):
        if len(row) != 2:
            raise SystemExit(f"ERROR: malformed metadata row {{path}}:{{line_number}}")
        observed[row[0]] = row[1]
mismatches = [
    f"{{key}} expected={{value!r}} observed={{observed.get(key)!r}}"
    for key, value in expected.items()
    if observed.get(key) != value
]
if mismatches:
    raise SystemExit(
        "ERROR: existing coverage metadata does not match this task; use a new run label "
        "or rerun with --force: " + "; ".join(mismatches)
    )
PY
}}

mkdir -p "$(dirname "$OUTPUT")" "$(dirname "$METADATA")"
if [[ {str(not args.force).lower()} == true ]] && \
   {{ [[ -e "$OUTPUT" ]] || [[ -e "$METADATA" ]]; }} && \
   ! {{ [[ -s "$OUTPUT" ]] && [[ -s "$METADATA" ]]; }}; then
    echo "ERROR: partial coverage outputs already exist; use a new run label or rerun with --force" >&2
    exit 1
fi
if [[ -s "$OUTPUT" && -s "$METADATA" && {str(not args.force).lower()} == true ]]; then
    validate_task_metadata
    echo "SKIP BUILD: nonempty coverage and metadata outputs already exist"
else
    BUILD_ARGS=(
        {shlex.quote(args.bam_window_coverage)} build
        --bam "$BAM"
        --cell-barcodes "$CELL_BARCODES"
        --output "$OUTPUT"
        --metadata "$METADATA"
        --window-size {args.coverage_window_size}
        --min-mapq "$MIN_MAPQ"
        --exclude-flags {args.coverage_exclude_flags}
        --tag {shlex.quote(args.barcode_tag)}
        --threads "$SLURM_CPUS_PER_TASK"
    )
    if [[ "$EMPTY_BARCODES" != "NA" ]]; then
        BUILD_ARGS+=(--empty-barcodes "$EMPTY_BARCODES")
    fi
    {normalize_line}
    {force_line}

    printf 'Command:'
    printf ' %q' "${{BUILD_ARGS[@]}}"
    printf '\\n'
    "${{BUILD_ARGS[@]}}"
fi

[[ -s "$OUTPUT" ]] || {{ echo "ERROR: coverage output is missing or empty: $OUTPUT" >&2; exit 1; }}
[[ -s "$METADATA" ]] || {{ echo "ERROR: coverage metadata is missing or empty: $METADATA" >&2; exit 1; }}
validate_task_metadata
python3 - "$OUTPUT" <<'PY'
import gzip
import math
import sys

expected = {EXPECTED_COVERAGE_HEADER!r}
with gzip.open(sys.argv[1], "rt", encoding="utf-8") as handle:
    header = None
    first_data = None
    for line in handle:
        stripped = line.rstrip("\\r\\n")
        if stripped.startswith("##") or not stripped:
            continue
        if header is None:
            header = stripped
            continue
        first_data = line
        break
if header != expected:
    raise SystemExit(f"ERROR: unexpected coverage header: {{header!r}}")
if not first_data:
    raise SystemExit("ERROR: coverage output contains no data rows")
PY

echo "PASS: coverage task $TASK_ID completed"
"""
    write_text(script, body, executable=True)
    bash_syntax_check(script)
    return script


def expected_per_library_paths(
    run_dir: Path, libraries: Sequence[int], modality: str
) -> list[tuple[int, Path, Path]]:
    base = run_dir / "coverage" / "per_library" / modality
    return [
        (
            lib,
            base / f"library_{lib}.{modality}.coverage.tsv.gz",
            base / f"library_{lib}.{modality}.coverage.metadata.tsv",
        )
        for lib in libraries
    ]


def metadata_contract_for_task(
    args: argparse.Namespace, row: Sequence[object]
) -> dict[str, str]:
    values = dict(zip(TASK_COLUMNS, (str(value) for value in row)))
    modality = values["modality"]
    min_mapq = args.rna_min_mapq if modality == "rna" else args.atac_min_mapq
    empty = values["empty_barcodes"]
    return {
        "schema": "bam_window_coverage_v1",
        "command": "build",
        "status": "complete",
        "bam": values["bam"],
        "cell_barcodes": values["cell_barcodes"],
        "empty_barcodes": "NONE" if empty == "NA" else empty,
        "coverage_output": values["coverage_output"],
        "window_size": str(args.coverage_window_size),
        "min_mapq": str(min_mapq),
        "exclude_flags": f"0x{args.coverage_exclude_flags:X}",
        "barcode_tag": args.barcode_tag,
        "normalize_10x": "false" if args.no_normalize_10x else "true",
    }


def read_key_value_metadata(path: Path) -> dict[str, str]:
    try:
        with path.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.reader(handle, delimiter="\t")
            header = next(reader, None)
            if header != ["key", "value"]:
                raise PipelineError(
                    f"coverage metadata has an unexpected header in {path}: {header!r}"
                )
            result: dict[str, str] = {}
            for line_number, row in enumerate(reader, start=2):
                if len(row) != 2:
                    raise PipelineError(
                        f"coverage metadata {path}:{line_number} is not key<TAB>value"
                    )
                result[row[0]] = row[1]
            return result
    except (OSError, UnicodeError) as exc:
        raise PipelineError(f"cannot read coverage metadata {path}: {exc}") from exc


def validate_metadata_contract(
    path: Path, expected: dict[str, str], label: str
) -> None:
    observed = read_key_value_metadata(path)
    mismatches = [
        f"{key}: expected {value!r}, observed {observed.get(key)!r}"
        for key, value in expected.items()
        if observed.get(key) != value
    ]
    if mismatches:
        raise PipelineError(
            f"{label} metadata does not match this run contract; use a new run label or "
            "rerun coverage with --force:\n  " + "\n  ".join(mismatches)
        )


def validate_per_library_products(
    args: argparse.Namespace, task_rows: Sequence[Sequence[object]]
) -> None:
    for row in task_rows:
        values = dict(zip(TASK_COLUMNS, row))
        lib = values["library_id"]
        modality = values["modality"]
        coverage = Path(str(values["coverage_output"]))
        metadata = Path(str(values["metadata_output"]))
        require_nonempty(coverage, f"library {lib} {modality} coverage output")
        require_nonempty(metadata, f"library {lib} {modality} coverage metadata")
        validate_metadata_contract(
            metadata,
            metadata_contract_for_task(args, row),
            f"library {lib} {modality}",
        )
        try:
            with gzip.open(coverage, "rt", encoding="utf-8") as handle:
                header = next(
                    (
                        line.rstrip("\r\n")
                        for line in handle
                        if line.strip() and not line.startswith("##")
                    ),
                    None,
                )
        except (OSError, UnicodeError) as exc:
            raise PipelineError(f"cannot read coverage output {coverage}: {exc}") from exc
        if header != EXPECTED_COVERAGE_HEADER:
            raise PipelineError(
                f"unexpected coverage header for {coverage}: {header!r}"
            )


def read_coverage_task_rows(path: Path) -> list[tuple[object, ...]]:
    require_nonempty(path, "base-run coverage task list")
    try:
        with path.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if reader.fieldnames != list(TASK_COLUMNS):
                raise PipelineError(
                    f"base-run coverage task list has an unexpected header in "
                    f"{path}: {reader.fieldnames!r}"
                )
            rows: list[tuple[object, ...]] = []
            for line_number, record in enumerate(reader, start=2):
                try:
                    task_id = int(record["task_id"])
                    library_id = int(record["library_id"])
                except (TypeError, ValueError) as exc:
                    raise PipelineError(
                        f"base-run coverage task list has invalid numeric fields at "
                        f"{path}:{line_number}"
                    ) from exc
                if task_id < 0 or library_id <= 0:
                    raise PipelineError(
                        f"base-run coverage task list has invalid identifiers at "
                        f"{path}:{line_number}"
                    )
                if record["modality"] not in MODALITIES:
                    raise PipelineError(
                        f"base-run coverage task list has unsupported modality at "
                        f"{path}:{line_number}: {record['modality']!r}"
                    )
                if any(record[column] in (None, "") for column in TASK_COLUMNS[3:]):
                    raise PipelineError(
                        f"base-run coverage task list has an incomplete row at "
                        f"{path}:{line_number}"
                    )
                rows.append(
                    (
                        task_id,
                        library_id,
                        record["modality"],
                        record["bam"],
                        record["cell_barcodes"],
                        record["empty_barcodes"],
                        record["coverage_output"],
                        record["metadata_output"],
                    )
                )
    except (OSError, UnicodeError) as exc:
        raise PipelineError(f"cannot read base-run coverage task list {path}: {exc}") from exc
    if not rows:
        raise PipelineError(f"base-run coverage task list is empty: {path}")
    return rows


def reindex_task_rows(
    rows: Sequence[Sequence[object]],
) -> list[tuple[object, ...]]:
    return [
        (task_id, *tuple(row)[1:])
        for task_id, row in enumerate(rows)
    ]


def require_same_underlying_file(
    requested: Path, recorded: Path, label: str
) -> None:
    require_nonempty(requested, f"requested {label}")
    require_nonempty(recorded, f"base-run recorded {label}")
    try:
        same = os.path.samefile(requested, recorded)
    except OSError as exc:
        raise PipelineError(
            f"cannot compare requested and base-run {label}: "
            f"{requested} vs {recorded}: {exc}"
        ) from exc
    if not same:
        raise PipelineError(
            f"cannot reuse coverage because the requested and base-run {label} "
            f"are different files: {requested} vs {recorded}"
        )


def resolve_incremental_coverage_rows(
    args: argparse.Namespace,
    base_run_dir: Path,
    requested_rows: Sequence[Sequence[object]],
    rebuild_libraries: Sequence[int],
) -> tuple[list[tuple[object, ...]], list[tuple[object, ...]], tuple[int, ...]]:
    source_task_file = base_run_dir / "control" / "tasks" / "coverage_tasks.tsv"
    source_rows = read_coverage_task_rows(source_task_file)
    source_by_key: dict[tuple[int, str], tuple[object, ...]] = {}
    for row in source_rows:
        values = dict(zip(TASK_COLUMNS, row))
        key = (int(values["library_id"]), str(values["modality"]))
        if key in source_by_key:
            raise PipelineError(
                "base-run coverage task list contains duplicate library/modality "
                f"rows: library {key[0]} {key[1]}"
            )
        source_by_key[key] = tuple(row)

    rebuild = set(rebuild_libraries)
    coverage_rows: list[tuple[object, ...]] = []
    aggregate_rows: list[tuple[object, ...]] = []
    reused_rows: list[tuple[object, ...]] = []
    reused_libraries: set[int] = set()
    for requested_row in requested_rows:
        requested = dict(zip(TASK_COLUMNS, requested_row))
        library = int(requested["library_id"])
        modality = str(requested["modality"])
        if library in rebuild:
            current_row = tuple(requested_row)
            coverage_rows.append(current_row)
            aggregate_rows.append(current_row)
            continue

        key = (library, modality)
        source_row = source_by_key.get(key)
        if source_row is None:
            raise PipelineError(
                f"base coverage run has no reusable task for library {library} "
                f"{modality}: {source_task_file}"
            )
        source = dict(zip(TASK_COLUMNS, source_row))
        require_same_underlying_file(
            Path(str(requested["bam"])),
            Path(str(source["bam"])),
            f"library {library} {modality} BAM",
        )
        require_same_underlying_file(
            Path(str(requested["cell_barcodes"])),
            Path(str(source["cell_barcodes"])),
            f"library {library} cell-barcode list",
        )
        requested_empty = str(requested["empty_barcodes"])
        source_empty = str(source["empty_barcodes"])
        if requested_empty == "NA" or source_empty == "NA":
            if requested_empty != source_empty:
                raise PipelineError(
                    f"cannot reuse coverage because library {library} {modality} "
                    "empty-barcode semantics changed"
                )
        else:
            require_same_underlying_file(
                Path(requested_empty),
                Path(source_empty),
                f"library {library} empty-barcode list",
            )
        aggregate_rows.append(source_row)
        reused_rows.append(source_row)
        reused_libraries.add(library)

    coverage_rows = reindex_task_rows(coverage_rows)
    aggregate_rows = reindex_task_rows(aggregate_rows)
    validate_per_library_products(args, reused_rows)
    return coverage_rows, aggregate_rows, tuple(sorted(reused_libraries))


def merge_populations(args: argparse.Namespace) -> tuple[str, ...]:
    populations = list(BASE_MERGE_POPULATIONS)
    if has_empty_barcode_source(args):
        populations.append("empty")
    return tuple(populations)


def aggregate_prefix(run_dir: Path, modality: str) -> Path:
    return run_dir / "coverage" / "aggregate" / modality


def aggregate_track(run_dir: Path, modality: str, population: str) -> Path:
    return Path(f"{aggregate_prefix(run_dir, modality)}.{population}.bedGraph.gz")


def aggregate_coverage_override(
    args: argparse.Namespace, modality: str
) -> str | None:
    return getattr(args, f"{modality}_aggregate_coverage")


def resolved_aggregate_track(
    args: argparse.Namespace, coverage_run_dir: Path, modality: str
) -> Path:
    override = aggregate_coverage_override(args, modality)
    if override:
        return Path(override)
    return aggregate_track(
        coverage_run_dir, modality, selected_population(args, modality)
    )


def aggregate_coverage_source(args: argparse.Namespace, modality: str) -> str:
    return (
        "explicit_override"
        if aggregate_coverage_override(args, modality)
        else "managed_run"
    )


def aggregate_coverage_contract_rows(
    args: argparse.Namespace,
    coverage_run_dir: Path,
    modalities: Sequence[str],
) -> list[tuple[object, object]]:
    rows: list[tuple[object, object]] = []
    requested = set(modalities)
    for modality in MODALITIES:
        if modality not in requested:
            rows.extend(
                (
                    (f"{modality}_coverage", "NONE"),
                    (f"{modality}_coverage_source", "NOT_APPLICABLE"),
                    (f"{modality}_coverage_size", "NONE"),
                    (f"{modality}_coverage_mtime_ns", "NONE"),
                )
            )
            continue
        track = resolved_aggregate_track(args, coverage_run_dir, modality)
        track_stat = track.stat()
        rows.extend(
            (
                (f"{modality}_coverage", track),
                (f"{modality}_coverage_source", aggregate_coverage_source(args, modality)),
                (f"{modality}_coverage_size", track_stat.st_size),
                (f"{modality}_coverage_mtime_ns", track_stat.st_mtime_ns),
            )
        )
    return rows


def render_aggregate_script(
    args: argparse.Namespace,
    run_dir: Path,
    modality: str,
    input_list: Path,
    contract_list: Path,
    populations: Sequence[str],
    required_populations: Sequence[str],
) -> Path:
    script = run_dir / "control" / "slurm" / f"aggregate_{modality}.sbatch"
    logs = run_dir / "logs"
    prefix = aggregate_prefix(run_dir, modality)
    metadata = run_dir / "coverage" / "aggregate" / f"{modality}.metadata.tsv"
    output_paths = [aggregate_track(run_dir, modality, population) for population in populations]
    all_data_checks = " && ".join(
        f'[[ -s {shlex.quote(str(path))} ]]' for path in output_paths
    )
    any_existing_checks = " || ".join(
        [f'[[ -e {shlex.quote(str(metadata))} ]]']
        + [
            check
            for path in output_paths
            for check in (
                f'[[ -e {shlex.quote(str(path))} ]]',
                f'[[ -e {shlex.quote(str(path) + ".tbi")} ]]',
            )
        ]
    )
    force_line = 'MERGE_ARGS+=(--force)' if args.force else ":"
    header = sbatch_header(
        job_name=f"vagg_{modality}_{args.run_label}",
        partition=args.partition,
        threads=args.aggregate_threads,
        memory=args.aggregate_memory,
        walltime=args.aggregate_time,
        stdout_path=logs / f"aggregate_{modality}_%j.out",
        stderr_path=logs / f"aggregate_{modality}_%j.err",
    )
    quoted_outputs = " ".join(shlex.quote(str(path)) for path in output_paths)
    quoted_populations = " ".join(shlex.quote(population) for population in populations)
    required_csv = ",".join(required_populations)
    modality_min_mapq = args.rna_min_mapq if modality == "rna" else args.atac_min_mapq
    body = f"""{header}

set -euo pipefail

{module_block()}

for executable in bgzip tabix python3; do
    command -v "$executable" >/dev/null 2>&1 || {{ echo "ERROR: missing command: $executable" >&2; exit 1; }}
done
[[ -x {shlex.quote(args.bam_window_coverage)} ]] || {{ echo "ERROR: coverage binary is not executable: {args.bam_window_coverage}" >&2; exit 1; }}
[[ -s {shlex.quote(str(input_list))} ]] || {{ echo "ERROR: missing merge input list: {input_list}" >&2; exit 1; }}
[[ -s {shlex.quote(str(contract_list))} ]] || {{ echo "ERROR: missing merge contract list: {contract_list}" >&2; exit 1; }}

python3 - {shlex.quote(str(contract_list))} <<'PY'
import csv
from pathlib import Path
import sys

contract_path = Path(sys.argv[1])
with contract_path.open("r", encoding="utf-8", newline="") as handle:
    rows = list(csv.DictReader(handle, delimiter="\t"))
if not rows:
    raise SystemExit(f"ERROR: merge contract list is empty: {{contract_path}}")
for row in rows:
    label = f"library {{row['library_id']}}"
    coverage = Path(row["coverage_path"])
    metadata = Path(row["metadata_path"])
    for path in (coverage, metadata):
        if not path.is_file() or path.stat().st_size == 0:
            raise SystemExit(f"ERROR: {{label}} input is missing or empty: {{path}}")
    expected = {{
        "schema": "bam_window_coverage_v1",
        "command": "build",
        "status": "complete",
        "bam": row["bam"],
        "cell_barcodes": row["cell_barcodes"],
        "empty_barcodes": row["empty_barcodes"],
        "coverage_output": row["coverage_path"],
        "window_size": row["window_size"],
        "min_mapq": row["min_mapq"],
        "exclude_flags": row["exclude_flags"],
        "barcode_tag": row["barcode_tag"],
        "normalize_10x": row["normalize_10x"],
    }}
    with metadata.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader, None)
        if header != ["key", "value"]:
            raise SystemExit(f"ERROR: unexpected metadata header for {{label}}: {{header!r}}")
        observed = {{}}
        for line_number, fields in enumerate(reader, 2):
            if len(fields) != 2:
                raise SystemExit(f"ERROR: malformed metadata row {{metadata}}:{{line_number}}")
            observed[fields[0]] = fields[1]
    mismatches = [
        f"{{key}} expected={{value!r}} observed={{observed.get(key)!r}}"
        for key, value in expected.items()
        if observed.get(key) != value
    ]
    if mismatches:
        raise SystemExit(
            f"ERROR: {{label}} coverage metadata is stale; use a new run label or "
            "rerun coverage with --force: " + "; ".join(mismatches)
        )
PY

mkdir -p {shlex.quote(str(prefix.parent))}
if {all_data_checks} && [[ -s {shlex.quote(str(metadata))} ]] && [[ {str(not args.force).lower()} == true ]]; then
    echo "SKIP MERGE: all {modality} aggregate coverage data products already exist; validating and rebuilding indexes as needed"
elif [[ {str(not args.force).lower()} == true ]] && {{ {any_existing_checks}; }}; then
    echo "ERROR: partial {modality} aggregate outputs exist; use a new run label or rerun with --force" >&2
    exit 1
else
    MERGE_ARGS=(
        {shlex.quote(args.bam_window_coverage)} merge
        --input-list {shlex.quote(str(input_list))}
        --output-prefix {shlex.quote(str(prefix))}
        --metadata {shlex.quote(str(metadata))}
        --populations {shlex.quote(','.join(populations))}
        --threads "$SLURM_CPUS_PER_TASK"
    )
    {force_line}
    printf 'Command:'
    printf ' %q' "${{MERGE_ARGS[@]}}"
    printf '\\n'
    "${{MERGE_ARGS[@]}}"
fi

[[ -s {shlex.quote(str(metadata))} ]] || {{ echo "ERROR: aggregate metadata is missing or empty: {metadata}" >&2; exit 1; }}
python3 - {shlex.quote(str(metadata))} {shlex.quote(str(input_list))} {shlex.quote(str(contract_list))} {shlex.quote(str(prefix))} <<'PY'
import csv
from pathlib import Path
import sys

metadata_path, input_list_path, contract_path, output_prefix = map(Path, sys.argv[1:])
with input_list_path.open("r", encoding="utf-8", newline="") as handle:
    input_rows = list(csv.DictReader(handle, delimiter="\t"))
expected = {{
    "schema": "bam_window_coverage_v1",
    "command": "merge",
    "status": "complete",
    "input_list": str(input_list_path),
    "output_prefix": str(output_prefix),
    "window_size": {str(args.coverage_window_size)!r},
    "libraries": str(len(input_rows)),
    "populations": {','.join(populations)!r},
    "min_mapq": {str(modality_min_mapq)!r},
    "exclude_flags": {str(args.coverage_exclude_flags)!r},
    "barcode_tag": {args.barcode_tag!r},
    "normalize_10x": {('false' if args.no_normalize_10x else 'true')!r},
    "empty_roster": {('present' if has_empty_barcode_source(args) else 'absent')!r},
}}
for index, row in enumerate(input_rows, 1):
    expected[f"library_{{index}}"] = row["library_id"] + r"\\t" + row["coverage_path"]
with metadata_path.open("r", encoding="utf-8", newline="") as handle:
    reader = csv.reader(handle, delimiter="\t")
    header = next(reader, None)
    if header != ["key", "value"]:
        raise SystemExit(f"ERROR: unexpected aggregate metadata header: {{header!r}}")
    observed = {{}}
    for line_number, fields in enumerate(reader, 2):
        if len(fields) != 2:
            raise SystemExit(f"ERROR: malformed aggregate metadata row {{metadata_path}}:{{line_number}}")
        observed[fields[0]] = fields[1]
mismatches = [
    f"{{key}} expected={{value!r}} observed={{observed.get(key)!r}}"
    for key, value in expected.items()
    if observed.get(key) != value
]
if mismatches:
    raise SystemExit(
        "ERROR: aggregate metadata does not match this request; use a new run label or "
        "rerun with --force: " + "; ".join(mismatches)
    )
PY

OUTPUTS=({quoted_outputs})
POPULATION_NAMES=({quoted_populations})
REQUIRED_POPULATIONS=,{required_csv},
for ((OUTPUT_INDEX=0; OUTPUT_INDEX<${{#OUTPUTS[@]}}; ++OUTPUT_INDEX)); do
    COVERAGE="${{OUTPUTS[$OUTPUT_INDEX]}}"
    POPULATION="${{POPULATION_NAMES[$OUTPUT_INDEX]}}"
    REQUIRED=false
    if [[ "$REQUIRED_POPULATIONS" == *",$POPULATION,"* ]]; then REQUIRED=true; fi
    [[ -s "$COVERAGE" ]] || {{ echo "ERROR: aggregate coverage is missing or empty: $COVERAGE" >&2; exit 1; }}
    # The merge binary validates every input row, writes sorted non-overlapping
    # BGZF, and atomically publishes exact row counts in this metadata file.
    # Reuse that count instead of decompressing the entire new aggregate again.
    ROWS="$(awk -F '\t' -v key="output_rows_$POPULATION" '$1 == key {{print $2; exit}}' {shlex.quote(str(metadata))})"
    [[ "$ROWS" =~ ^[0-9]+$ ]] || {{ echo "ERROR: aggregate metadata lacks output_rows_$POPULATION" >&2; exit 1; }}
    if [[ "$REQUIRED" == true && "$ROWS" -eq 0 ]]; then
        echo "ERROR: aggregate coverage contains no data rows: $COVERAGE" >&2
        exit 1
    fi
    if (( ROWS > 0 )); then
        # Merge emits BGZF directly. Preserve a valid, current index so a
        # no-op resume does not make unchanged coverage look newer than panels.
        REBUILD_INDEX=false
        if [[ ! -s "$COVERAGE.tbi" || "$COVERAGE" -nt "$COVERAGE.tbi" ]]; then
            REBUILD_INDEX=true
        elif ! tabix -l "$COVERAGE" | awk 'NF {{found=1}} END {{exit !found}}'; then
            REBUILD_INDEX=true
        fi
        if [[ "$REBUILD_INDEX" == true ]]; then
            tabix -f -p bed "$COVERAGE"
        fi
        [[ -s "$COVERAGE.tbi" ]] || {{ echo "ERROR: tabix index is missing or empty: $COVERAGE.tbi" >&2; exit 1; }}
        tabix -l "$COVERAGE" | awk 'NF {{found=1}} END {{exit !found}}'
    else
        rm -f -- "$COVERAGE.tbi"
        echo "Diagnostic population has zero rows and is not indexed: $POPULATION"
    fi
done

echo "PASS: {modality} aggregate coverage completed"
"""
    write_text(script, body, executable=True)
    bash_syntax_check(script)
    return script


def write_merge_input_lists(
    args: argparse.Namespace,
    run_dir: Path,
    task_rows: Sequence[Sequence[object]],
    modalities: Sequence[str],
) -> dict[str, tuple[Path, Path]]:
    result: dict[str, tuple[Path, Path]] = {}
    for modality in modalities:
        selected_rows = [
            row for row in task_rows if str(row[TASK_COLUMNS.index("modality")]) == modality
        ]
        path = run_dir / "control" / "tasks" / f"aggregate_{modality}_inputs.tsv"
        rows = [
            (
                row[TASK_COLUMNS.index("library_id")],
                row[TASK_COLUMNS.index("coverage_output")],
            )
            for row in selected_rows
        ]
        write_tsv(path, ("library_id", "coverage_path"), rows)
        contract_path = (
            run_dir / "control" / "tasks" / f"aggregate_{modality}_contracts.tsv"
        )
        contract_columns = (
            "library_id",
            "coverage_path",
            "metadata_path",
            "bam",
            "cell_barcodes",
            "empty_barcodes",
            "window_size",
            "min_mapq",
            "exclude_flags",
            "barcode_tag",
            "normalize_10x",
        )
        contract_rows = []
        for row in selected_rows:
            values = dict(zip(TASK_COLUMNS, row))
            contract = metadata_contract_for_task(args, row)
            contract_rows.append(
                (
                    values["library_id"],
                    values["coverage_output"],
                    values["metadata_output"],
                    contract["bam"],
                    contract["cell_barcodes"],
                    contract["empty_barcodes"],
                    contract["window_size"],
                    contract["min_mapq"],
                    contract["exclude_flags"],
                    contract["barcode_tag"],
                    contract["normalize_10x"],
                )
            )
        write_tsv(contract_path, contract_columns, contract_rows)
        result[modality] = (path, contract_path)
    return result


def selected_population(args: argparse.Namespace, modality: str) -> str:
    override = (
        args.rna_coverage_population
        if modality == "rna"
        else args.atac_coverage_population
    )
    return override or args.coverage_population


def panel_label(modality: str, population: str) -> str:
    return f"{modality}_{population}"


def target_label(value: int) -> str:
    if value % 1_000_000 == 0:
        return f"{value // 1_000_000}M"
    if value % 1_000 == 0:
        return f"{value // 1_000}K"
    return str(value)


def atac_candidate_pool_path(args: argparse.Namespace, run_dir: Path) -> Path:
    atac_label = panel_label("atac", selected_population(args, "atac"))
    return run_dir / "panels" / f"tet.vars.{atac_label}.candidates_all.bcf"


def atac_audit_prefix(args: argparse.Namespace, run_dir: Path) -> Path:
    atac_label = panel_label("atac", selected_population(args, "atac"))
    return run_dir / "panels" / f"tet.vars.{atac_label}.capacity_audit"


def atac_audit_outputs(
    args: argparse.Namespace, run_dir: Path, modalities: Sequence[str]
) -> dict[str, Path]:
    if "atac" not in modalities:
        return {}
    prefix = atac_audit_prefix(args, run_dir)
    return {
        "atac_candidate_sfs_audit": Path(f"{prefix}.candidate_sfs.tsv"),
        "atac_pair_capacity_audit": Path(f"{prefix}.pair_capacity.tsv"),
        "atac_species_tiers_audit": Path(f"{prefix}.species_tiers.tsv"),
    }


def atac_diagnostics_enabled(
    args: argparse.Namespace, modalities: Sequence[str]
) -> bool:
    return "atac" in modalities and args.stage != "ATAC_RESELECT_AND_PLOTS"


def mito_panel_enabled(
    args: argparse.Namespace, modalities: Sequence[str]
) -> bool:
    # Full-source ATAC builds always retain a separate mitochondrial panel.
    # A saved nuclear ATAC candidate pool cannot recreate chrM/NUMT outputs.
    return args.mito_panel or atac_diagnostics_enabled(args, modalities)


def panel_outputs(
    args: argparse.Namespace, run_dir: Path, modalities: Sequence[str]
) -> dict[str, Path]:
    panels = run_dir / "panels"
    primary = "rna" if "rna" in modalities else "atac"
    outputs: dict[str, Path] = {}
    if primary == "rna":
        rna_label = panel_label("rna", selected_population(args, "rna"))
        outputs["rna_demux"] = panels / (
            f"tet.vars.{rna_label}.demux_{target_label(args.rna_demux_target)}.bcf"
        )
        outputs["rna_het"] = panels / (
            f"tet.vars.{rna_label}.het_{target_label(args.rna_het_target)}.bcf"
        )
    species_modality = args.species_coverage
    if species_modality == "auto":
        species_modality = primary
    species_population = selected_population(args, species_modality)
    species_label = panel_label(species_modality, species_population)
    outputs["species"] = panels / (
        f"tet.vars.{species_label}.species_{target_label(args.species_target)}.bcf"
    )
    if "atac" in modalities:
        atac_population = selected_population(args, "atac")
        atac_label = panel_label("atac", atac_population)
        outputs["atac_demux"] = panels / (
            f"tet.vars.{atac_label}.demux_{target_label(args.atac_demux_target)}.bcf"
        )
        outputs["atac_het"] = panels / (
            f"tet.vars.{atac_label}.het_{target_label(args.atac_het_target)}.bcf"
        )
        if atac_diagnostics_enabled(args, modalities):
            outputs["atac_numt"] = panels / (
                f"tet.vars.{atac_label}.numt_diagnostic.bcf"
            )
    if mito_panel_enabled(args, modalities):
        mito_label = (
            panel_label("atac", selected_population(args, "atac"))
            if "atac" in modalities and "rna" not in modalities
            else "rna_cells"
        )
        outputs["mt_fusion_ratio"] = panels / (
            f"tet.vars.{mito_label}.mt_fusion_ratio.bcf"
        )
    return outputs


def mito_sidecar_outputs(
    args: argparse.Namespace, run_dir: Path, modalities: Sequence[str]
) -> dict[str, Path]:
    if not mito_panel_enabled(args, modalities):
        return {}
    mito_label = (
        panel_label("atac", selected_population(args, "atac"))
        if "atac" in modalities and "rna" not in modalities
        else "rna_cells"
    )
    prefix = run_dir / "panels" / f"tet.vars.{mito_label}.mt_fusion_ratio"
    return {
        "mt_site_manifest": Path(f"{prefix}.site_manifest.tsv"),
        "mt_pair_audit": Path(f"{prefix}.pair_audit.tsv"),
        "mt_haplotype_groups": Path(f"{prefix}.haplotype_groups.tsv"),
        "mt_haplotype_pairwise": Path(f"{prefix}.haplotype_pairwise.tsv"),
        "mt_sample_audit": Path(f"{prefix}.sample_audit.tsv"),
        "mt_sites_bed": Path(f"{prefix}.sites.bed"),
    }


def validate_downsample_inputs(
    args: argparse.Namespace, run_dir: Path, modalities: Sequence[str]
) -> None:
    require_executable(
        Path(args.downsample_vcf_parallel), "VCF selector binary"
    )
    for path_text, label in (
        (args.source_bcf, "source BCF"),
        (args.source_bcf_index, "source BCF CSI"),
        (args.gtf, "GTF"),
        (args.numt_bed, "NUMT BED"),
        (args.panel_metadata, "panel metadata TSV"),
        (args.pool_combinations, "pool combinations TSV"),
    ):
        require_nonempty(Path(path_text), label)
    if args.atac_regulatory_bed:
        require_nonempty(Path(args.atac_regulatory_bed), "ATAC regulatory BED")
    if args.blacklist:
        require_nonempty(Path(args.blacklist), "blacklist source")
        if not args.blacklist_staged_path:
            raise PipelineError("blacklist was not staged for panel selection")
        require_nonempty(Path(args.blacklist_staged_path), "staged blacklist")
    if "atac" in modalities:
        panel = Path(args.rna_demux_panel)
        require_nonempty(panel, "production RNA demux soft-overlap reference")
        require_nonempty(
            Path(args.rna_demux_panel + ".csi"),
            "production RNA demux soft-overlap reference CSI",
        )
    expected_index = args.source_bcf + ".csi"
    if os.path.normpath(args.source_bcf_index) != os.path.normpath(expected_index):
        raise PipelineError(
            "--source-bcf-index must be the sidecar path <source-bcf>.csi because "
            "downsample_vcf_parallel discovers the index beside the BCF"
        )
    species_modality = args.species_coverage
    if species_modality != "auto" and species_modality not in modalities:
        raise PipelineError(
            f"--species-coverage {species_modality} requires that modality to be selected"
        )


def validate_aggregate_products(
    args: argparse.Namespace,
    coverage_run_dir: Path,
    modalities: Sequence[str],
    libraries: Sequence[int],
) -> None:
    required: dict[str, set[str]] = {
        modality: {selected_population(args, modality)} for modality in modalities
    }
    if args.species_coverage != "auto":
        required.setdefault(args.species_coverage, set()).add(
            selected_population(args, args.species_coverage)
        )
    for modality in sorted(required):
        override = aggregate_coverage_override(args, modality)
        if override:
            track = Path(override)
            require_nonempty(track, f"explicit {modality} aggregate coverage")
            require_nonempty(
                Path(str(track) + ".tbi"),
                f"explicit {modality} aggregate coverage tabix index",
            )
            require_tabix_readable(
                track, f"explicit {modality} aggregate coverage"
            )
            continue

        for population in sorted(required[modality]):
            track = aggregate_track(coverage_run_dir, modality, population)
            require_nonempty(track, f"{modality} {population} aggregate coverage")
            require_nonempty(
                Path(str(track) + ".tbi"),
                f"{modality} {population} tabix index",
            )

        metadata_path = (
            coverage_run_dir
            / "coverage"
            / "aggregate"
            / f"{modality}.metadata.tsv"
        )
        require_nonempty(metadata_path, f"{modality} aggregate metadata")
        observed = read_key_value_metadata(metadata_path)
        input_list_path = (
            coverage_run_dir
            / "control"
            / "tasks"
            / f"aggregate_{modality}_inputs.tsv"
        )
        require_nonempty(input_list_path, f"{modality} aggregate input list")
        try:
            with input_list_path.open("r", encoding="utf-8", newline="") as handle:
                reader = csv.DictReader(handle, delimiter="\t")
                if reader.fieldnames != ["library_id", "coverage_path"]:
                    raise PipelineError(
                        f"{modality} aggregate input list has an unexpected "
                        f"header: {reader.fieldnames!r}"
                    )
                input_rows = list(reader)
        except (OSError, UnicodeError) as exc:
            raise PipelineError(
                f"cannot read {modality} aggregate input list {input_list_path}: {exc}"
            ) from exc
        min_mapq = args.rna_min_mapq if modality == "rna" else args.atac_min_mapq
        expected = {
            "schema": "bam_window_coverage_v1",
            "command": "merge",
            "status": "complete",
            "input_list": str(input_list_path),
            "output_prefix": str(
                coverage_run_dir / "coverage" / "aggregate" / modality
            ),
            "window_size": str(args.coverage_window_size),
            "libraries": str(len(libraries)),
            "min_mapq": str(min_mapq),
            "exclude_flags": str(args.coverage_exclude_flags),
            "barcode_tag": args.barcode_tag,
            "normalize_10x": "false" if args.no_normalize_10x else "true",
        }
        mismatches = [
            f"{key}: expected {value!r}, observed {observed.get(key)!r}"
            for key, value in expected.items()
            if observed.get(key) != value
        ]
        observed_populations = set(observed.get("populations", "").split(","))
        for population in sorted(required[modality]):
            if population not in observed_populations:
                mismatches.append(
                    f"populations: required {population!r}, observed "
                    f"{observed.get('populations')!r}"
                )
        if (
            "empty" in required[modality]
            and observed.get("empty_roster") != "present"
        ):
            mismatches.append(
                "empty_roster: selected empty population requires 'present', "
                f"observed {observed.get('empty_roster')!r}"
            )
        observed_libraries: list[int] = []
        for index, row in enumerate(input_rows, start=1):
            raw_library = row.get("library_id", "")
            coverage_path = row.get("coverage_path", "")
            try:
                library = int(raw_library)
            except (TypeError, ValueError):
                mismatches.append(
                    f"aggregate input row {index}: invalid library_id {raw_library!r}"
                )
                continue
            observed_libraries.append(library)
            if not coverage_path or not Path(coverage_path).is_absolute():
                mismatches.append(
                    f"aggregate input row {index}: coverage_path must be absolute, "
                    f"observed {coverage_path!r}"
                )
                continue
            value = observed.get(f"library_{index}")
            expected_value = f"{library}\\t{coverage_path}"
            if value != expected_value:
                mismatches.append(
                    f"library_{index}: expected {expected_value!r}, observed {value!r}"
                )
        if tuple(observed_libraries) != tuple(libraries):
            mismatches.append(
                f"aggregate input libraries: expected {tuple(libraries)!r}, "
                f"observed {tuple(observed_libraries)!r}"
            )
        if len(input_rows) != len(libraries):
            mismatches.append(
                f"aggregate input row count: expected {len(libraries)}, "
                f"observed {len(input_rows)}"
            )
        if mismatches:
            raise PipelineError(
                f"{modality} aggregate metadata does not match the requested "
                "downsampling coverage contract; use the correct "
                "--coverage-run-label or rebuild coverage:\n  "
                + "\n  ".join(mismatches)
            )


def render_atac_candidate_script(
    args: argparse.Namespace,
    run_dir: Path,
    coverage_run_dir: Path,
    libraries: Sequence[int],
) -> tuple[Path, dict[str, Path]]:
    script = run_dir / "control" / "slurm" / "atac_candidates.sbatch"
    logs = run_dir / "logs"
    output = atac_candidate_pool_path(args, run_dir)
    coverage = aggregate_track(
        coverage_run_dir, "atac", selected_population(args, "atac")
    )
    soft_overlap_panel = args.rna_demux_panel
    audit_prefix = atac_audit_prefix(args, run_dir)
    audit_outputs = atac_audit_outputs(args, run_dir, ("atac",))
    final_audit_array = " ".join(
        shlex.quote(str(path)) for path in audit_outputs.values()
    )
    ram_audit_array = " ".join(
        f'"$RAM_AUDIT_PREFIX.{suffix}"'
        for suffix in ("candidate_sfs.tsv", "pair_capacity.tsv", "species_tiers.tsv")
    )
    audit_complete_test = " && ".join(
        f"-s {shlex.quote(str(path))}" for path in audit_outputs.values()
    )
    audit_partial_test = " || ".join(
        f"-e {shlex.quote(str(path))}" for path in audit_outputs.values()
    )
    regulatory_stage = ""
    regulatory_check = ""
    regulatory_arg = ""
    if args.atac_regulatory_bed:
        regulatory_stage = f"""
STAGED_ATAC_REGULATORY_BED="$RAMDIR/atac_regulatory.bed"
cp -f {shlex.quote(args.atac_regulatory_bed)} "$STAGED_ATAC_REGULATORY_BED"
"""
        regulatory_check = ' "$STAGED_ATAC_REGULATORY_BED"'
        regulatory_arg = '\n    --atac_regulatory_bed "$STAGED_ATAC_REGULATORY_BED"'
    if args.atac_regulatory_required:
        regulatory_arg += "\n    --atac_regulatory_required"
    blacklist_arg = (
        f"\n    --blacklist {shlex.quote(args.blacklist_staged_path)}"
        if args.blacklist_staged_path
        else ""
    )
    request_file = run_dir / "control" / "tasks" / "atac_candidates.request.tsv"
    final_request = run_dir / "panels" / "atac_candidates.request.tsv"
    request_rows: list[tuple[object, object]] = [
        ("schema", "cellbouncer_atac_candidate_request_v2"),
        ("libraries", ",".join(str(lib) for lib in libraries)),
        ("downsample_vcf_parallel", args.downsample_vcf_parallel),
        ("source_bcf", args.source_bcf),
        ("source_bcf_index", args.source_bcf_index),
        ("atac_coverage", coverage),
        ("atac_population", selected_population(args, "atac")),
        ("atac_min_coverage", args.atac_min_coverage),
        ("bin_size", args.bin_size),
        ("gtf", args.gtf),
        ("numt_bed", args.numt_bed),
        (
            "atac_regulatory_bed",
            args.atac_regulatory_bed or "GTF_2KB_UPSTREAM_PROMOTERS",
        ),
        ("atac_regulatory_required", str(args.atac_regulatory_required).lower()),
        ("panel_metadata", args.panel_metadata),
        ("pool_combinations", args.pool_combinations),
        ("pool_objective_policy", "rna_style_scalar_pool_discrimination_no_identity_pair_floors"),
        ("min_pair_score", args.min_pair_score),
        ("min_pairwise", args.min_pairwise),
        ("species_target", args.species_target),
        ("species_outgroup", args.species_outgroup),
        ("species_outgroup_max_fraction", args.species_outgroup_max_fraction),
        ("rna_demux_panel_selection_role", "exact_selected_20M_soft_position_overlap_score_multiplier"),
        ("rna_demux_soft_overlap_panel", soft_overlap_panel),
        ("rna_demux_soft_overlap_multiplier", args.atac_soft_overlap_multiplier),
        ("rna_het_panel_selection_role", "not_used_plot_reference_only"),
        ("rna_species_panel_selection_role", "not_used_plot_reference_only"),
        ("atac_rna_overlap_contract", "only_selected_rna_demux_20M_is_softly_penalized"),
        ("output_atac_candidates", output),
    ]
    for key, path in audit_outputs.items():
        request_rows.append((f"output_{key}", path))
    append_blacklist_contract(request_rows, args)
    for key, path in (
        ("downsample_vcf_parallel", Path(args.downsample_vcf_parallel)),
        ("source_bcf", Path(args.source_bcf)),
        ("source_bcf_index", Path(args.source_bcf_index)),
        ("atac_coverage", coverage),
        ("atac_coverage_index", Path(str(coverage) + ".tbi")),
        ("gtf", Path(args.gtf)),
        ("numt_bed", Path(args.numt_bed)),
        *((
            ("atac_regulatory_bed", Path(args.atac_regulatory_bed)),
        ) if args.atac_regulatory_bed else ()),
        ("panel_metadata", Path(args.panel_metadata)),
        ("pool_combinations", Path(args.pool_combinations)),
        ("rna_demux_soft_overlap_panel", Path(soft_overlap_panel)),
        ("rna_demux_soft_overlap_panel_index", Path(soft_overlap_panel + ".csi")),
    ):
        file_stat = path.stat()
        request_rows.extend(
            ((f"{key}_size", file_stat.st_size),
             (f"{key}_mtime_ns", file_stat.st_mtime_ns))
        )
    write_tsv(request_file, ("key", "value"), request_rows)

    header = sbatch_header(
        job_name=f"vcandidates_{args.run_label}",
        partition=args.partition,
        threads=args.downsample_threads,
        memory=args.downsample_memory,
        walltime=args.downsample_time,
        stdout_path=logs / "atac_candidates_%j.out",
        stderr_path=logs / "atac_candidates_%j.err",
    )
    source_inputs = (
        args.source_bcf,
        args.source_bcf_index,
        str(coverage),
        str(coverage) + ".tbi",
        args.gtf,
        args.numt_bed,
        *((args.blacklist_staged_path,) if args.blacklist_staged_path else ()),
        *((args.atac_regulatory_bed,) if args.atac_regulatory_bed else ()),
        args.panel_metadata,
        args.pool_combinations,
        soft_overlap_panel,
        soft_overlap_panel + ".csi",
    )
    body = f"""{header}

set -euo pipefail

{module_block(include_bcftools=True)}

for executable in bcftools cmp cp df grep gzip mv rm stat tabix; do
    command -v "$executable" >/dev/null 2>&1 || {{ echo "ERROR: missing command: $executable" >&2; exit 1; }}
done
[[ -x {shlex.quote(args.downsample_vcf_parallel)} ]] || {{ echo "ERROR: candidate builder is not executable: {args.downsample_vcf_parallel}" >&2; exit 1; }}
for INPUT in {shell_join(source_inputs)}; do
    [[ -s "$INPUT" ]] || {{ echo "ERROR: required candidate input is missing or empty: $INPUT" >&2; exit 1; }}
done

if [[ {str(not args.force).lower()} == true && -s {shlex.quote(str(output))} && -s {shlex.quote(str(output) + '.csi')} && {audit_complete_test} ]]; then
    [[ -s {shlex.quote(str(final_request))} ]] || {{ echo "ERROR: candidate pool lacks its request contract; use a new run label or --force" >&2; exit 1; }}
    cmp -s {shlex.quote(str(request_file))} {shlex.quote(str(final_request))} || {{ echo "ERROR: candidate pool was built for a different request; use a new run label or --force" >&2; exit 1; }}
    COUNT="$(bcftools index -n {shlex.quote(str(output))})"
    [[ "$COUNT" =~ ^[0-9]+$ && "$COUNT" -gt 0 ]] || {{ echo "ERROR: existing candidate pool is empty or unreadable" >&2; exit 1; }}
    echo "SKIP: reusable ATAC candidate pool is complete ($COUNT variants)"
    exit 0
fi
if [[ {str(not args.force).lower()} == true && ( -e {shlex.quote(str(output))} || -e {shlex.quote(str(output) + '.csi')} || -e {shlex.quote(str(final_request))} || {audit_partial_test} ) ]]; then
    echo "ERROR: incomplete ATAC candidate pool exists; use a new run label or --force" >&2
    exit 1
fi

RAMDIR="/dev/shm/vcf_candidates_${{SLURM_JOB_ID}}"
case "$RAMDIR" in
    /dev/shm/vcf_candidates_[0-9]*) ;;
    *) echo "ERROR: refusing unsafe RAM-work path: $RAMDIR" >&2; exit 1 ;;
esac
cleanup_ramdir() {{ rm -rf -- "$RAMDIR"; }}
trap cleanup_ramdir EXIT
mkdir -p "$RAMDIR"

NEEDED_BYTES=0
for INPUT in {shell_join(source_inputs)}; do
    SIZE="$(stat -c %s "$INPUT")"
    NEEDED_BYTES=$((NEEDED_BYTES + SIZE))
done
SOURCE_BYTES="$(stat -c %s {shlex.quote(args.source_bcf)})"
OUTPUT_RESERVE_BYTES=$((SOURCE_BYTES + 34359738368))
AVAILABLE_BYTES="$(df --output=avail -B1 /dev/shm | awk 'NR == 2 {{print $1}}')"
if (( AVAILABLE_BYTES < NEEDED_BYTES + OUTPUT_RESERVE_BYTES )); then
    echo "ERROR: insufficient /dev/shm capacity for ATAC candidate build" >&2
    exit 1
fi

STAGED_BCF="$RAMDIR/input.bcf"
STAGED_COVERAGE="$RAMDIR/atac.coverage.bedGraph.gz"
STAGED_GTF="$RAMDIR/annotations.gtf"
STAGED_NUMTS="$RAMDIR/numts.bed"
STAGED_PANEL_METADATA="$RAMDIR/panel_metadata.tsv"
STAGED_POOL_COMBINATIONS="$RAMDIR/pool_combinations.tsv"
RNA_DEMUX_SOFT_OVERLAP="$RAMDIR/rna_demux_soft_overlap.bcf"
RAM_CANDIDATES="$RAMDIR/{output.name}"
RAM_AUDIT_PREFIX="$RAMDIR/{audit_prefix.name}"
cp -f {shlex.quote(args.source_bcf)} "$STAGED_BCF"
cp -f {shlex.quote(args.source_bcf_index)} "$STAGED_BCF.csi"
cp -f {shlex.quote(str(coverage))} "$STAGED_COVERAGE"
cp -f {shlex.quote(str(coverage) + '.tbi')} "$STAGED_COVERAGE.tbi"
if [[ {shlex.quote(args.gtf)} == *.gz ]]; then
    gzip -cd {shlex.quote(args.gtf)} > "$STAGED_GTF"
else
    cp -f {shlex.quote(args.gtf)} "$STAGED_GTF"
fi
cp -f {shlex.quote(args.numt_bed)} "$STAGED_NUMTS"
cp -f {shlex.quote(args.panel_metadata)} "$STAGED_PANEL_METADATA"
cp -f {shlex.quote(args.pool_combinations)} "$STAGED_POOL_COMBINATIONS"
cp -f {shlex.quote(soft_overlap_panel)} "$RNA_DEMUX_SOFT_OVERLAP"
cp -f {shlex.quote(soft_overlap_panel + '.csi')} "$RNA_DEMUX_SOFT_OVERLAP.csi"
{regulatory_stage}
for STAGED in "$STAGED_BCF" "$STAGED_BCF.csi" "$STAGED_COVERAGE" "$STAGED_COVERAGE.tbi" "$STAGED_GTF" "$STAGED_NUMTS" "$STAGED_PANEL_METADATA" "$STAGED_POOL_COMBINATIONS" "$RNA_DEMUX_SOFT_OVERLAP" "$RNA_DEMUX_SOFT_OVERLAP.csi"{regulatory_check}; do
    [[ -s "$STAGED" ]] || {{ echo "ERROR: staged candidate input is missing or empty: $STAGED" >&2; exit 1; }}
done

CANDIDATE_ARGS=(
    {shlex.quote(args.downsample_vcf_parallel)}
    --vcf "$STAGED_BCF"
    --gtf "$STAGED_GTF"
    --threads "$SLURM_CPUS_PER_TASK"
    --bin_size {args.bin_size}
    --numts_bed "$STAGED_NUMTS"
    --enable_pool_scoring
    --pool_combinations "$STAGED_POOL_COMBINATIONS"
    --panel_metadata "$STAGED_PANEL_METADATA"
    --min_pair_score {args.min_pair_score}
    --min_pairwise {args.min_pairwise}
    --species_num {args.species_target}
    --species_coverage atac
    --species_outgroup {shlex.quote(args.species_outgroup)}
    --species_outgroup_max_fraction {args.species_outgroup_max_fraction}
    --atac_cov "$STAGED_COVERAGE"
    --atac_min_cov {args.atac_min_coverage}
    --atac_candidates_output "$RAM_CANDIDATES"
    --atac_candidates_only
    --atac_soft_overlap_vcf "$RNA_DEMUX_SOFT_OVERLAP"
    --atac_soft_overlap_multiplier {args.atac_soft_overlap_multiplier}
    --atac_audit_prefix "$RAM_AUDIT_PREFIX"{regulatory_arg}{blacklist_arg}
)
printf 'Command:'
printf ' %q' "${{CANDIDATE_ARGS[@]}}"
printf '\n'
"${{CANDIDATE_ARGS[@]}}"

[[ -s "$RAM_CANDIDATES" && -s "$RAM_CANDIDATES.csi" ]] || {{ echo "ERROR: candidate BCF/CSI was not created" >&2; exit 1; }}
bcftools view -h "$RAM_CANDIDATES" | grep -Fqx '##cellbouncer_atac_candidate_pool_schema=atac_candidate_union_v2' || {{ echo "ERROR: candidate BCF lacks the v2 schema header" >&2; exit 1; }}
bcftools view -h "$RAM_CANDIDATES" | grep -Fqx '##cellbouncer_atac_rna_overlap_policy=soft_penalty' || {{ echo "ERROR: candidate BCF lacks the soft-overlap policy header" >&2; exit 1; }}
bcftools view -h "$RAM_CANDIDATES" | grep -Fqx '##cellbouncer_atac_pool_objective_policy=rna_style_scalar_pool_discrimination_no_identity_pair_floors' || {{ echo "ERROR: candidate BCF lacks the RNA-style scalar pool-ranking policy" >&2; exit 1; }}
bcftools view -h "$RAM_CANDIDATES" | grep -Fqx '##cellbouncer_atac_soft_overlap_multiplier={args.atac_soft_overlap_multiplier:.6f}' || {{ echo "ERROR: candidate BCF soft-overlap multiplier does not match the request" >&2; exit 1; }}
COUNT="$(bcftools index -n "$RAM_CANDIDATES")"
[[ "$COUNT" =~ ^[0-9]+$ && "$COUNT" -gt 0 ]] || {{ echo "ERROR: candidate BCF has no indexed variants" >&2; exit 1; }}
RAM_AUDITS=({ram_audit_array})
FINAL_AUDITS=({final_audit_array})
for AUDIT in "${{RAM_AUDITS[@]}}"; do
    [[ -s "$AUDIT" ]] || {{ echo "ERROR: ATAC capacity audit is missing or empty: $AUDIT" >&2; exit 1; }}
done

mkdir -p {shlex.quote(str(output.parent))}
TEMP_BCF="{str(output)[:-4]}.partial.${{SLURM_JOB_ID}}.bcf"
cp -f "$RAM_CANDIDATES" "$TEMP_BCF"
cp -f "$RAM_CANDIDATES.csi" "$TEMP_BCF.csi"
[[ "$(bcftools index -n "$TEMP_BCF")" == "$COUNT" ]] || {{ echo "ERROR: temporary candidate count mismatch" >&2; exit 1; }}
mv -f "$TEMP_BCF" {shlex.quote(str(output))}
mv -f "$TEMP_BCF.csi" {shlex.quote(str(output) + '.csi')}
for ((INDEX=0; INDEX<${{#RAM_AUDITS[@]}}; ++INDEX)); do
    RAM_AUDIT="${{RAM_AUDITS[$INDEX]}}"
    FINAL_AUDIT="${{FINAL_AUDITS[$INDEX]}}"
    TEMP_AUDIT="$FINAL_AUDIT.partial.${{SLURM_JOB_ID}}"
    cp -f "$RAM_AUDIT" "$TEMP_AUDIT"
    [[ -s "$TEMP_AUDIT" ]] || {{ echo "ERROR: temporary ATAC capacity audit is empty: $TEMP_AUDIT" >&2; exit 1; }}
    mv -f "$TEMP_AUDIT" "$FINAL_AUDIT"
done
TEMP_REQUEST={shlex.quote(str(final_request))}.partial.${{SLURM_JOB_ID}}
cp -f {shlex.quote(str(request_file))} "$TEMP_REQUEST"
mv -f "$TEMP_REQUEST" {shlex.quote(str(final_request))}

echo "PASS: reusable ATAC candidate pool completed"
echo "Candidate variants: $COUNT  {output}"
for AUDIT in "${{FINAL_AUDITS[@]}}"; do echo "ATAC capacity audit: $AUDIT"; done
"""
    write_text(script, body, executable=True)
    bash_syntax_check(script)
    return script, {"atac_candidates": output, **audit_outputs}


def render_downsample_script(
    args: argparse.Namespace,
    run_dir: Path,
    modalities: Sequence[str],
    coverage_run_dir: Path,
    libraries: Sequence[int],
) -> tuple[Path, dict[str, Path]]:
    script = run_dir / "control" / "slurm" / "downsample.sbatch"
    logs = run_dir / "logs"
    bcf_outputs = panel_outputs(args, run_dir, modalities)
    if "atac" in modalities and args.stage != "ATAC_RESELECT_AND_PLOTS":
        bcf_outputs["atac_candidates"] = atac_candidate_pool_path(args, run_dir)
    sidecar_outputs = {
        **mito_sidecar_outputs(args, run_dir, modalities),
        **atac_audit_outputs(args, run_dir, modalities),
    }
    outputs = {**bcf_outputs, **sidecar_outputs}
    primary = "rna" if "rna" in modalities else "atac"
    primary_population = selected_population(args, primary)
    primary_track = resolved_aggregate_track(args, coverage_run_dir, primary)
    primary_demux_target = (
        args.rna_demux_target if primary == "rna" else args.atac_demux_target
    )
    primary_het_target = (
        args.rna_het_target if primary == "rna" else args.atac_het_target
    )
    primary_demux_key = "rna_demux" if primary == "rna" else None
    species_modality = args.species_coverage
    if species_modality == "auto":
        species_modality = primary
    species_balance_args = ""
    if species_modality == "atac":
        species_balance_args = (
            f"\n    --species_outgroup {shlex.quote(args.species_outgroup)}"
            f"\n    --species_outgroup_max_fraction "
            f"{args.species_outgroup_max_fraction}"
        )
    atac_soft_overlap_panels = (
        (args.rna_demux_panel,)
        if "atac" in modalities
        else ()
    )
    header = sbatch_header(
        job_name=f"vpanel_{args.run_label}",
        partition=args.partition,
        threads=args.downsample_threads,
        memory=args.downsample_memory,
        walltime=args.downsample_time,
        stdout_path=logs / "downsample_%j.out",
        stderr_path=logs / "downsample_%j.err",
    )
    selected_tracks = [
        resolved_aggregate_track(args, coverage_run_dir, modality)
        for modality in modalities
    ]
    selected_aggregate_metadata = [
        coverage_run_dir
        / "coverage"
        / "aggregate"
        / f"{modality}.metadata.tsv"
        for modality in modalities
        if not aggregate_coverage_override(args, modality)
    ]
    selected_track_checks = "\n".join(
        f'[[ -s {shlex.quote(str(track))} && -s {shlex.quote(str(track) + ".tbi")} ]] || '
        f'{{ echo "ERROR: aggregate coverage/index is incomplete: {track}" >&2; exit 1; }}'
        for track in selected_tracks
    )
    selected_metadata_checks = "\n".join(
        f'[[ -s {shlex.quote(str(path))} ]] || '
        f'{{ echo "ERROR: aggregate metadata is missing or empty: {path}" >&2; exit 1; }}'
        for path in selected_aggregate_metadata
    )
    upstream_freshness_array = " ".join(
        shlex.quote(str(path))
        for track in selected_tracks
        for path in (track, Path(str(track) + ".tbi"))
    )
    upstream_freshness_array += " " + " ".join(
        shlex.quote(str(path)) for path in selected_aggregate_metadata
    )
    upstream_freshness_array += " " + shell_join(
        (
            args.source_bcf,
            args.source_bcf_index,
            args.gtf,
            args.numt_bed,
            *((args.blacklist_staged_path,) if args.blacklist_staged_path else ()),
            args.panel_metadata,
            args.pool_combinations,
            *((args.atac_regulatory_bed,) if args.atac_regulatory_bed else ()),
            *atac_soft_overlap_panels,
            *(panel + ".csi" for panel in atac_soft_overlap_panels),
        )
    )
    track_contig_checks = "\n".join(
        f"require_contig_overlap {shlex.quote(str(track))}"
        for track in selected_tracks
    )
    track_tabix_checks = "\n".join(
        f"require_tabix_readable {shlex.quote(str(track))}"
        for track in selected_tracks
    )
    mito_contig_check = ""
    if mito_panel_enabled(args, modalities):
        mito_track = primary_track
        mito_contig_check = f"""
[[ -n "${{SOURCE_CONTIGS[chrM]+present}}" ]] || {{ echo "ERROR: source BCF does not contain chrM required by --mito-panel" >&2; exit 1; }}
tabix -l {shlex.quote(str(mito_track))} | awk '$0 == "chrM" {{found=1}} END {{exit !found}}' || {{ echo "ERROR: selected coverage does not contain chrM required by the mitochondrial panel" >&2; exit 1; }}"""
    stage_track_lines: list[str] = []
    staged_track_variables: dict[str, str] = {}
    for modality in modalities:
        population = selected_population(args, modality)
        source_track = resolved_aggregate_track(args, coverage_run_dir, modality)
        variable = f"{modality.upper()}_COVERAGE"
        staged_track_variables[modality] = f"${variable}"
        stage_track_lines.extend(
            [
                f'{variable}="$RAMDIR/{modality}.{population}.bedGraph.gz"',
                f"cp -f {shlex.quote(str(source_track))} \"${variable}\"",
                f"cp -f {shlex.quote(str(source_track) + '.tbi')} \"${variable}.tbi\"",
            ]
        )
    staged_tracks_block = "\n".join(stage_track_lines)
    staged_regulatory_block = ""
    staged_regulatory_check = ""
    if args.atac_regulatory_bed:
        staged_regulatory_block = f"""
STAGED_ATAC_REGULATORY_BED="$RAMDIR/atac_regulatory.bed"
cp -f {shlex.quote(args.atac_regulatory_bed)} "$STAGED_ATAC_REGULATORY_BED"
"""
        staged_regulatory_check = ' "$STAGED_ATAC_REGULATORY_BED"'

    soft_overlap_stage_lines: list[str] = []
    soft_overlap_check_variables: list[str] = []
    soft_overlap_argument_lines: list[str] = []
    for index, source_panel in enumerate(atac_soft_overlap_panels, start=1):
        variable = f"RNA_SOFT_OVERLAP_{index}"
        soft_overlap_check_variables.extend((f'"${variable}"', f'"${variable}.csi"'))
        soft_overlap_stage_lines.extend(
            (
                f'{variable}="$RAMDIR/rna_soft_overlap_{index}.bcf"',
                f"cp -f {shlex.quote(source_panel)} \"${variable}\"",
                f"cp -f {shlex.quote(source_panel + '.csi')} \"${variable}.csi\"",
            )
        )
        soft_overlap_argument_lines.append(
            f'    --atac_soft_overlap_vcf "${variable}"'
        )
    staged_soft_overlaps_block = "\n".join(soft_overlap_stage_lines)
    staged_soft_overlap_checks = " ".join(soft_overlap_check_variables)
    atac_soft_overlap_args = "\n".join(soft_overlap_argument_lines)
    if "atac" in modalities:
        atac_soft_overlap_args += (
            f"\n    --atac_soft_overlap_multiplier "
            f"{args.atac_soft_overlap_multiplier}"
        )
    blacklist_arg = (
        f"\n    --blacklist {shlex.quote(args.blacklist_staged_path)}"
        if args.blacklist_staged_path
        else ""
    )

    output_variable_names: dict[str, str] = {}
    output_assignment_lines: list[str] = []
    for key, final_path in bcf_outputs.items():
        variable = "RAM_" + re.sub(r"[^A-Za-z0-9]", "_", key).upper()
        output_variable_names[key] = variable
        output_assignment_lines.append(
            f'{variable}="$RAMDIR/{shlex.quote(final_path.name)}"'
        )
    output_assignments = "\n".join(output_assignment_lines)
    ram_output_array = " ".join(
        f'"${output_variable_names[key]}"' for key in bcf_outputs
    )
    final_output_array = " ".join(
        shlex.quote(str(path)) for path in bcf_outputs.values()
    )
    # A NUMT diagnostic panel is a valid, useful negative result when no
    # callable ATAC-covered source variants overlap the supplied NUMT BED.
    # All biological selection panels and the reusable candidate pool must
    # still contain at least one indexed record.
    allow_empty_output_array = " ".join(
        "true" if key == "atac_numt" else "false"
        for key in bcf_outputs
    )
    candidate_output_header_check = ""
    if "atac_candidates" in output_variable_names:
        candidate_variable = output_variable_names["atac_candidates"]
        candidate_output_header_check = f"""
bcftools view -h "${candidate_variable}" | grep -Fqx '##cellbouncer_atac_candidate_pool_schema=atac_candidate_union_v2' || {{ echo "ERROR: emitted ATAC candidate pool lacks the v2 schema header" >&2; exit 1; }}
bcftools view -h "${candidate_variable}" | grep -Fqx '##cellbouncer_atac_rna_overlap_policy=soft_penalty' || {{ echo "ERROR: emitted ATAC candidate pool lacks the soft-overlap policy header" >&2; exit 1; }}
bcftools view -h "${candidate_variable}" | grep -Fqx '##cellbouncer_atac_soft_overlap_multiplier={args.atac_soft_overlap_multiplier:.6f}' || {{ echo "ERROR: emitted ATAC candidate-pool multiplier differs from this request" >&2; exit 1; }}
"""

    if primary == "rna":
        primary_output_variable = output_variable_names[primary_demux_key]
        generic_primary_args = f'''\
    --output "${{{primary_output_variable}}}"
    --num {primary_demux_target}
    --cov "{staged_track_variables['rna']}"
    --min_cov {args.min_coverage}'''
    else:
        # Native ATAC selection is self-contained.  Omitting generic --output
        # prevents the selector from constructing and writing a discarded RNA/
        # generic allocation-support panel.
        primary_output_variable = ""
        generic_primary_args = ""

    primary_het_args = ""
    if primary == "rna":
        primary_het_args = f"""
    --het_output \"${output_variable_names['rna_het']}\"
    --het_num {primary_het_target}"""
    native_atac_args = ""
    if "atac" in modalities:
        candidate_output_arg = ""
        if "atac_candidates" in output_variable_names:
            candidate_output_arg = (
                f'\n    --atac_candidates_output '
                f'"${output_variable_names["atac_candidates"]}"'
            )
        regulatory_arg = (
            '\n    --atac_regulatory_bed "$STAGED_ATAC_REGULATORY_BED"'
            if args.atac_regulatory_bed
            else ""
        )
        regulatory_required_arg = (
            "\n    --atac_regulatory_required"
            if args.atac_regulatory_required else ""
        )
        numt_arg = (
            f'\n    --atac_numt_output "${output_variable_names["atac_numt"]}"'
            if "atac_numt" in output_variable_names
            else ""
        )
        audit_arg = '\n    --atac_audit_prefix "$RAM_ATAC_AUDIT_PREFIX"'
        native_atac_args = f"""
    --atac_cov \"{staged_track_variables['atac']}\"
    --atac_output \"${output_variable_names['atac_demux']}\"
    --atac_num {args.atac_demux_target}
    --atac_het_output \"${output_variable_names['atac_het']}\"
    --atac_het_num {args.atac_het_target}
    --atac_min_cov {args.atac_min_coverage}{regulatory_arg}{regulatory_required_arg}{numt_arg}{candidate_output_arg}{audit_arg}"""
    allocation_note = (
        'echo "Selector mode: native ATAC only (generic allocation pass disabled)"'
        if primary == "atac" else ""
    )

    sidecar_variable_names: dict[str, str] = {}
    sidecar_assignment_lines: list[str] = []
    for key, final_path in sidecar_outputs.items():
        variable = "RAM_" + re.sub(r"[^A-Za-z0-9]", "_", key).upper()
        sidecar_variable_names[key] = variable
        sidecar_assignment_lines.append(
            f'{variable}="$RAMDIR/{shlex.quote(final_path.name)}"'
        )
    sidecar_assignments = "\n".join(sidecar_assignment_lines)
    audit_prefix_assignment = ""
    if "atac" in modalities:
        audit_prefix_assignment = (
            f'RAM_ATAC_AUDIT_PREFIX="$RAMDIR/'
            f'{atac_audit_prefix(args, run_dir).name}"'
        )
    mito_args = ""
    if mito_panel_enabled(args, modalities):
        strict_pair_flag = (
            "\n    --mt_require_pair_targets" if args.mt_require_pair_targets else ""
        )
        mito_args = f"""
    --mt_output \"${output_variable_names['mt_fusion_ratio']}\"
    --mt_site_manifest \"${sidecar_variable_names['mt_site_manifest']}\"
    --mt_pair_audit \"${sidecar_variable_names['mt_pair_audit']}\"
    --mt_haplotype_groups \"${sidecar_variable_names['mt_haplotype_groups']}\"
    --mt_haplotype_pairwise \"${sidecar_variable_names['mt_haplotype_pairwise']}\"
    --mt_sample_audit \"${sidecar_variable_names['mt_sample_audit']}\"
    --mt_sites_bed \"${sidecar_variable_names['mt_sites_bed']}\"
    --mt_min_depth {args.mt_min_depth}
    --mt_homoplasmy_af {args.mt_homoplasmy_af}
    --mt_min_coverage {args.mt_min_coverage}
    --mt_min_pair_sites {args.mt_min_pair_sites}
    --mt_min_ambient_sites {args.mt_min_ambient_sites}{strict_pair_flag}"""
    sidecar_publish_block = ""
    final_sidecar_array = ""
    if sidecar_outputs:
        ram_sidecar_array = " ".join(
            f'"${sidecar_variable_names[key]}"' for key in sidecar_outputs
        )
        final_sidecar_array = " ".join(
            shlex.quote(str(path)) for path in sidecar_outputs.values()
        )
        sidecar_publish_block = f"""
RAM_SIDECARS=({ram_sidecar_array})
FINAL_SIDECARS=({final_sidecar_array})
for ((INDEX=0; INDEX<${{#RAM_SIDECARS[@]}}; ++INDEX)); do
    RAM_SIDECAR="${{RAM_SIDECARS[$INDEX]}}"
    FINAL_SIDECAR="${{FINAL_SIDECARS[$INDEX]}}"
    [[ -s "$RAM_SIDECAR" ]] || {{ echo "ERROR: panel sidecar is missing or empty: $RAM_SIDECAR" >&2; exit 1; }}
    TEMP_SIDECAR="$FINAL_SIDECAR.partial.${{SLURM_JOB_ID}}"
    cp -f "$RAM_SIDECAR" "$TEMP_SIDECAR"
    [[ -s "$TEMP_SIDECAR" ]] || {{ echo "ERROR: temporary panel sidecar is empty: $TEMP_SIDECAR" >&2; exit 1; }}
    mv -f "$TEMP_SIDECAR" "$FINAL_SIDECAR"
    [[ -s "$FINAL_SIDECAR" ]] || {{ echo "ERROR: failed to publish panel sidecar: $FINAL_SIDECAR" >&2; exit 1; }}
    echo "Panel sidecar: $FINAL_SIDECAR"
done
"""

    request_file = run_dir / "control" / "tasks" / "downsample_request.tsv"
    final_request_file = run_dir / "panels" / "downsample.request.tsv"
    request_rows = [
        ("schema", "cellbouncer_vcf_panel_request_v2"),
        ("libraries", ",".join(str(lib) for lib in libraries)),
        ("modalities", ",".join(modalities)),
        ("downsample_vcf_parallel", args.downsample_vcf_parallel),
        (
            "source_mode",
            "saved_atac_candidate_pool"
            if args.stage == "ATAC_RESELECT_AND_PLOTS"
            else "full_source_bcf",
        ),
        ("source_bcf", args.source_bcf),
        ("source_bcf_index", args.source_bcf_index),
        ("coverage_run_dir", coverage_run_dir),
        ("gtf", args.gtf),
        ("numt_bed", args.numt_bed),
        ("panel_metadata", args.panel_metadata),
        ("pool_combinations", args.pool_combinations),
        ("pool_objective_policy", "rna_style_scalar_pool_discrimination_no_identity_pair_floors" if "atac" in modalities else "NOT_APPLICABLE"),
        ("rna_demux_panel_selection_role", "exact_selected_20M_soft_position_overlap_score_multiplier" if "atac" in modalities else "NOT_APPLICABLE"),
        ("rna_demux_soft_overlap_panel", args.rna_demux_panel if "atac" in modalities else "NONE"),
        ("rna_demux_soft_overlap_multiplier", args.atac_soft_overlap_multiplier if "atac" in modalities else "NONE"),
        ("rna_het_panel_selection_role", "not_used_plot_reference_only" if "atac" in modalities else "NOT_APPLICABLE"),
        ("rna_het_panel_plot_reference", args.rna_het_panel if "atac" in modalities else "NONE"),
        ("rna_species_panel_selection_role", "not_used_plot_reference_only" if "atac" in modalities else "NOT_APPLICABLE"),
        ("rna_species_panel_plot_reference", args.rna_species_panel if "atac" in modalities else "NONE"),
        ("atac_rna_overlap_contract", "only_selected_rna_demux_20M_is_softly_penalized" if "atac" in modalities else "NOT_APPLICABLE"),
        ("primary_modality", primary),
        ("primary_population", primary_population),
        ("species_modality", species_modality),
        ("rna_population", selected_population(args, "rna") if "rna" in modalities else "NONE"),
        ("atac_population", selected_population(args, "atac") if "atac" in modalities else "NONE"),
        *aggregate_coverage_contract_rows(
            args, coverage_run_dir, modalities
        ),
        ("coverage_window_size", args.coverage_window_size),
        ("rna_min_mapq", args.rna_min_mapq),
        ("atac_min_mapq", args.atac_min_mapq),
        ("coverage_exclude_flags", f"0x{args.coverage_exclude_flags:X}"),
        ("barcode_tag", args.barcode_tag),
        ("normalize_10x", str(not args.no_normalize_10x).lower()),
        ("empty_barcode_source", args.empty_barcodes_template or args.empty_drop_roster or "NONE"),
        ("rna_demux_target", args.rna_demux_target),
        ("rna_het_target", args.rna_het_target),
        ("atac_demux_target", args.atac_demux_target),
        ("atac_het_target", args.atac_het_target),
        ("species_target", args.species_target),
        ("bin_size", args.bin_size),
        ("min_coverage", args.min_coverage),
        ("atac_min_coverage", args.atac_min_coverage),
        ("min_pair_score", args.min_pair_score),
        ("min_pairwise", args.min_pairwise),
        ("max_het_per_bin", args.max_het_per_bin),
        ("species_max_per_bin", args.species_max_per_bin),
        ("seed", args.seed),
        ("atac_regulatory_bed", args.atac_regulatory_bed or "GTF_2KB_UPSTREAM_PROMOTERS"),
        ("atac_regulatory_required", str(args.atac_regulatory_required).lower()),
        ("species_outgroup", args.species_outgroup if species_modality == "atac" else "NONE"),
        ("species_outgroup_max_fraction", args.species_outgroup_max_fraction if species_modality == "atac" else "NONE"),
        ("mito_panel", str(mito_panel_enabled(args, modalities)).lower()),
        ("mt_min_depth", args.mt_min_depth),
        ("mt_homoplasmy_af", args.mt_homoplasmy_af),
        ("mt_min_coverage", args.mt_min_coverage),
        ("mt_min_pair_sites", args.mt_min_pair_sites),
        ("mt_min_ambient_sites", args.mt_min_ambient_sites),
        ("mt_require_pair_targets", str(args.mt_require_pair_targets).lower()),
    ]
    for key, path in outputs.items():
        request_rows.append((f"output_{key}", path))
    append_blacklist_contract(request_rows, args)
    for key, path_text in (
        ("downsample_vcf_parallel", args.downsample_vcf_parallel),
        ("source_bcf", args.source_bcf),
        ("source_bcf_index", args.source_bcf_index),
        ("gtf", args.gtf),
        ("numt_bed", args.numt_bed),
        ("panel_metadata", args.panel_metadata),
        ("pool_combinations", args.pool_combinations),
        *((
            ("atac_regulatory_bed", args.atac_regulatory_bed),
        ) if args.atac_regulatory_bed else ()),
        *((
            ("rna_demux_soft_overlap_panel", args.rna_demux_panel),
            ("rna_demux_soft_overlap_panel_index", args.rna_demux_panel + ".csi"),
        ) if "atac" in modalities else ()),
    ):
        file_stat = Path(path_text).stat()
        request_rows.extend(
            (
                (f"{key}_size", file_stat.st_size),
                (f"{key}_mtime_ns", file_stat.st_mtime_ns),
            )
        )
    write_tsv(request_file, ("key", "value"), request_rows)

    pre_skip_validation = f"""
REQUESTED_BCFS=({final_output_array})
REQUESTED_ALLOW_EMPTY=({allow_empty_output_array})
REQUESTED_SIDECARS=({final_sidecar_array})
if [[ {str(not args.force).lower()} == true ]]; then
    ALL_COMPLETE=true
    COMPLETE_COUNT=0
    MISSING_COUNT=0
    for ((INDEX=0; INDEX<${{#REQUESTED_BCFS[@]}}; ++INDEX)); do
        BCF="${{REQUESTED_BCFS[$INDEX]}}"
        ALLOW_EMPTY="${{REQUESTED_ALLOW_EMPTY[$INDEX]}}"
        if [[ -s "$BCF" && -s "$BCF.csi" ]]; then
            COUNT="$(bcftools index -n "$BCF" 2>/dev/null)" || {{ echo "ERROR: existing BCF/CSI pair is invalid; rerun with --force: $BCF" >&2; exit 1; }}
            [[ "$COUNT" =~ ^[0-9]+$ ]] || {{ echo "ERROR: existing BCF has an invalid indexed count; rerun with --force: $BCF" >&2; exit 1; }}
            if [[ "$ALLOW_EMPTY" != true && "$COUNT" -le 0 ]]; then
                echo "ERROR: existing biological panel has no indexed variants; rerun with --force: $BCF" >&2
                exit 1
            fi
            ((COMPLETE_COUNT+=1))
        elif [[ -e "$BCF" || -e "$BCF.csi" ]]; then
            echo "ERROR: partial BCF/CSI pair exists; rerun with --force: $BCF" >&2
            exit 1
        else
            ALL_COMPLETE=false
            ((MISSING_COUNT+=1))
        fi
    done
    for SIDECAR in "${{REQUESTED_SIDECARS[@]}}"; do
        if [[ -s "$SIDECAR" ]]; then
            ((COMPLETE_COUNT+=1))
        elif [[ -e "$SIDECAR" ]]; then
            echo "ERROR: empty requested sidecar exists; rerun with --force: $SIDECAR" >&2
            exit 1
        else
            ALL_COMPLETE=false
            ((MISSING_COUNT+=1))
        fi
    done
    if [[ "$ALL_COMPLETE" == true ]]; then
        [[ -s {shlex.quote(str(final_request_file))} ]] || {{ echo "ERROR: completed panels lack their request contract; use a new run label or --force" >&2; exit 1; }}
        cmp -s {shlex.quote(str(request_file))} {shlex.quote(str(final_request_file))} || {{ echo "ERROR: completed panels were built for a different selector request; use a new run label or --force" >&2; exit 1; }}
        UPSTREAM_COVERAGE_ARTIFACTS=({upstream_freshness_array})
        for UPSTREAM in "${{UPSTREAM_COVERAGE_ARTIFACTS[@]}}"; do
            for BCF in "${{REQUESTED_BCFS[@]}}"; do
                [[ ! "$UPSTREAM" -nt "$BCF" ]] || {{ echo "ERROR: an upstream input is newer than completed panels; use a new run label or --force: $UPSTREAM" >&2; exit 1; }}
            done
            [[ ! "$UPSTREAM" -nt {shlex.quote(str(final_request_file))} ]] || {{ echo "ERROR: an upstream input is newer than the panel request completion record; use a new run label or --force: $UPSTREAM" >&2; exit 1; }}
        done
        echo "SKIP: all requested panels, indexes, and sidecars are complete and readable"
        exit 0
    fi
    if (( COMPLETE_COUNT > 0 || MISSING_COUNT == 0 )) || [[ -e {shlex.quote(str(final_request_file))} ]]; then
        echo "ERROR: incomplete panel family already exists; use a new run label or --force" >&2
        exit 1
    fi
fi
rm -f {shlex.quote(str(final_request_file))}
"""

    output_reserve_multiplier = 3 if len(bcf_outputs) >= 5 else 2
    candidate_source_header_check = ""
    if args.stage == "ATAC_RESELECT_AND_PLOTS":
        candidate_source_header_check = f"""
bcftools view -h "$STAGED_BCF" | grep -Fqx '##cellbouncer_atac_candidate_pool_schema=atac_candidate_union_v2' || {{
    echo "ERROR: source BCF is not a reusable v2 ATAC candidate pool; rebuild it so RNA auxiliary-panel loci and softly penalized RNA-demux loci are recoverable" >&2
    exit 1
}}
bcftools view -h "$STAGED_BCF" | grep -Fqx '##cellbouncer_atac_rna_overlap_policy=soft_penalty' || {{
    echo "ERROR: ATAC candidate pool lacks the required RNA-demux soft-overlap policy" >&2
    exit 1
}}
bcftools view -h "$STAGED_BCF" | grep -Fqx '##cellbouncer_atac_pool_objective_policy=rna_style_scalar_pool_discrimination_no_identity_pair_floors' || {{
    echo "ERROR: ATAC candidate pool was built with the obsolete derived-pair objective policy; rebuild it from the full source BCF" >&2
    exit 1
}}
bcftools view -h "$STAGED_BCF" | grep -Fqx '##cellbouncer_atac_soft_overlap_multiplier={args.atac_soft_overlap_multiplier:.6f}' || {{
    echo "ERROR: ATAC candidate-pool soft-overlap multiplier differs from this request" >&2
    exit 1
}}
"""
    body = f"""{header}

set -euo pipefail

{module_block(include_bcftools=True)}

for executable in bcftools bgzip tabix cmp cp df grep gzip mv python3 rm stat; do
    command -v "$executable" >/dev/null 2>&1 || {{ echo "ERROR: missing command: $executable" >&2; exit 1; }}
done
[[ -x {shlex.quote(args.downsample_vcf_parallel)} ]] || {{ echo "ERROR: downsampling binary is not executable: {args.downsample_vcf_parallel}" >&2; exit 1; }}

for INPUT in {shell_join((args.source_bcf, args.source_bcf_index, args.gtf, args.numt_bed, *((args.blacklist_staged_path,) if args.blacklist_staged_path else ()), args.panel_metadata, args.pool_combinations, *((args.atac_regulatory_bed,) if args.atac_regulatory_bed else ()), *atac_soft_overlap_panels, *(panel + '.csi' for panel in atac_soft_overlap_panels)))}; do
    [[ -s "$INPUT" ]] || {{ echo "ERROR: required input is missing or empty: $INPUT" >&2; exit 1; }}
done
{selected_track_checks}
{selected_metadata_checks}
bcftools index -n {shlex.quote(args.source_bcf)} >/dev/null
declare -A SOURCE_CONTIGS=()
while IFS=$'\t' read -r CONTIG _REST; do
    [[ -n "$CONTIG" ]] && SOURCE_CONTIGS["$CONTIG"]=1
done < <(bcftools index -s {shlex.quote(args.source_bcf)})
(( ${{#SOURCE_CONTIGS[@]}} > 0 )) || {{ echo "ERROR: source BCF index reports no contigs" >&2; exit 1; }}
require_tabix_readable() {{
    local track="$1"
    tabix -l "$track" | awk 'NF {{found=1}} END {{exit !found}}' || {{
        echo "ERROR: aggregate coverage is not Tabix-readable: $track" >&2
        exit 1
    }}
}}
require_contig_overlap() {{
    local track="$1" contig overlap=""
    while IFS= read -r contig; do
        if [[ -n "${{SOURCE_CONTIGS[$contig]+present}}" ]]; then overlap="$contig"; break; fi
    done < <(tabix -l "$track")
    [[ -n "$overlap" ]] || {{ echo "ERROR: coverage/source BCF contig names have no overlap: $track" >&2; exit 1; }}
    echo "Coverage/source contig overlap: $overlap  $track"
}}
{track_tabix_checks}
{track_contig_checks}
{mito_contig_check}

{pre_skip_validation}

mkdir -p {shlex.quote(str(run_dir / 'panels'))}
echo "Primary selector channel: {primary} ({primary_population})"
echo "Species coverage: {species_modality} ({selected_population(args, species_modality)})"
{allocation_note}

RAMDIR="/dev/shm/vcf_panel_${{SLURM_JOB_ID}}"
case "$RAMDIR" in
    /dev/shm/vcf_panel_[0-9]*) ;;
    *) echo "ERROR: refusing unsafe RAM-work path: $RAMDIR" >&2; exit 1 ;;
esac
cleanup_ramdir() {{ rm -rf -- "$RAMDIR"; }}
trap cleanup_ramdir EXIT
mkdir -p "$RAMDIR"

STAGE_INPUTS=(
    {shell_join((args.source_bcf, args.source_bcf_index, args.gtf, args.numt_bed, *((args.blacklist_staged_path,) if args.blacklist_staged_path else ()), args.panel_metadata, args.pool_combinations, *((args.atac_regulatory_bed,) if args.atac_regulatory_bed else ()), *atac_soft_overlap_panels, *(panel + '.csi' for panel in atac_soft_overlap_panels)))}
    {' '.join(shlex.quote(str(track)) + ' ' + shlex.quote(str(track) + '.tbi') for track in selected_tracks)}
)
NEEDED_BYTES=0
for INPUT in "${{STAGE_INPUTS[@]}}"; do
    SIZE="$(stat -c %s "$INPUT")"
    NEEDED_BYTES=$((NEEDED_BYTES + SIZE))
done
SOURCE_BYTES="$(stat -c %s {shlex.quote(args.source_bcf)})"
OUTPUT_RESERVE_BYTES=$((SOURCE_BYTES * {output_reserve_multiplier} + 34359738368))
AVAILABLE_BYTES="$(df --output=avail -B1 /dev/shm | awk 'NR == 2 {{print $1}}')"
if (( AVAILABLE_BYTES < NEEDED_BYTES + OUTPUT_RESERVE_BYTES )); then
    echo "ERROR: insufficient /dev/shm capacity for staged selector run" >&2
    echo "  available bytes: $AVAILABLE_BYTES" >&2
    echo "  staged input bytes: $NEEDED_BYTES" >&2
    echo "  reserved output bytes: $OUTPUT_RESERVE_BYTES" >&2
    exit 1
fi

STAGED_BCF="$RAMDIR/input.bcf"
STAGED_GTF="$RAMDIR/annotations.gtf"
STAGED_NUMTS="$RAMDIR/numts.bed"
STAGED_PANEL_METADATA="$RAMDIR/panel_metadata.tsv"
STAGED_POOL_COMBINATIONS="$RAMDIR/pool_combinations.tsv"
cp -f {shlex.quote(args.source_bcf)} "$STAGED_BCF"
cp -f {shlex.quote(args.source_bcf_index)} "$STAGED_BCF.csi"
if [[ {shlex.quote(args.gtf)} == *.gz ]]; then
    gzip -cd {shlex.quote(args.gtf)} > "$STAGED_GTF"
else
    cp -f {shlex.quote(args.gtf)} "$STAGED_GTF"
fi
cp -f {shlex.quote(args.numt_bed)} "$STAGED_NUMTS"
cp -f {shlex.quote(args.panel_metadata)} "$STAGED_PANEL_METADATA"
cp -f {shlex.quote(args.pool_combinations)} "$STAGED_POOL_COMBINATIONS"
{staged_tracks_block}
{staged_soft_overlaps_block}
{staged_regulatory_block}
{output_assignments}
{sidecar_assignments}
{audit_prefix_assignment}

for STAGED in "$STAGED_BCF" "$STAGED_BCF.csi" "$STAGED_GTF" "$STAGED_NUMTS" "$STAGED_PANEL_METADATA" "$STAGED_POOL_COMBINATIONS" {staged_soft_overlap_checks}{staged_regulatory_check}; do
    [[ -s "$STAGED" ]] || {{ echo "ERROR: staged input is missing or empty: $STAGED" >&2; exit 1; }}
done
bcftools index -n "$STAGED_BCF" >/dev/null
{candidate_source_header_check}

DOWNSAMPLE_ARGS=(
    {shlex.quote(args.downsample_vcf_parallel)}
    --vcf "$STAGED_BCF"
{generic_primary_args}
    {primary_het_args.lstrip()}
    --gtf "$STAGED_GTF"
    --bin_size {args.bin_size}
    --threads "$SLURM_CPUS_PER_TASK"
    --numts_bed "$STAGED_NUMTS"
    --enable_pool_scoring
    --pool_combinations "$STAGED_POOL_COMBINATIONS"
    --panel_metadata "$STAGED_PANEL_METADATA"
    --species_output "${output_variable_names['species']}"
    --species_num {args.species_target}
    --species_coverage {species_modality}
    --min_pair_score {args.min_pair_score}
    --min_pairwise {args.min_pairwise}
    --max_het_per_bin {args.max_het_per_bin}
    --species_max_per_bin {args.species_max_per_bin}
    --seed {args.seed}{species_balance_args}{native_atac_args}{mito_args}{blacklist_arg}
{atac_soft_overlap_args}
)
printf 'Command:'
printf ' %q' "${{DOWNSAMPLE_ARGS[@]}}"
printf '\\n'
"${{DOWNSAMPLE_ARGS[@]}}"

{candidate_output_header_check}
RAM_OUTPUTS=({ram_output_array})
FINAL_OUTPUTS=({final_output_array})
ALLOW_EMPTY_OUTPUTS=({allow_empty_output_array})
for ((INDEX=0; INDEX<${{#RAM_OUTPUTS[@]}}; ++INDEX)); do
    RAM_BCF="${{RAM_OUTPUTS[$INDEX]}}"
    FINAL_BCF="${{FINAL_OUTPUTS[$INDEX]}}"
    ALLOW_EMPTY="${{ALLOW_EMPTY_OUTPUTS[$INDEX]}}"
    [[ -s "$RAM_BCF" ]] || {{ echo "ERROR: panel BCF is missing or empty: $RAM_BCF" >&2; exit 1; }}
    COUNT=""
    if [[ -s "$RAM_BCF.csi" ]]; then
        COUNT="$(bcftools index -n "$RAM_BCF" 2>/dev/null || true)"
    fi
    if [[ ! "$COUNT" =~ ^[0-9]+$ ]]; then
        bcftools index -f -c "$RAM_BCF"
        COUNT="$(bcftools index -n "$RAM_BCF")"
    fi
    [[ -s "$RAM_BCF.csi" ]] || {{ echo "ERROR: panel CSI is missing or empty: $RAM_BCF.csi" >&2; exit 1; }}
    [[ "$COUNT" =~ ^[0-9]+$ ]] || {{ echo "ERROR: panel has an invalid indexed variant count: $RAM_BCF" >&2; exit 1; }}
    if [[ "$ALLOW_EMPTY" != true && "$COUNT" -le 0 ]]; then
        echo "ERROR: biological panel has no indexed variants: $RAM_BCF" >&2
        exit 1
    fi
    [[ "$FINAL_BCF" == *.bcf ]] || {{ echo "ERROR: requested panel output does not end in .bcf: $FINAL_BCF" >&2; exit 1; }}
    TEMP_BCF="${{FINAL_BCF%.bcf}}.partial.${{SLURM_JOB_ID}}.bcf"
    cp -f "$RAM_BCF" "$TEMP_BCF"
    cp -f "$RAM_BCF.csi" "$TEMP_BCF.csi"
    TEMP_COUNT="$(bcftools index -n "$TEMP_BCF")"
    [[ "$TEMP_COUNT" == "$COUNT" ]] || {{ echo "ERROR: temporary published panel count mismatch: $TEMP_BCF" >&2; exit 1; }}
    mv -f "$TEMP_BCF" "$FINAL_BCF"
    mv -f "$TEMP_BCF.csi" "$FINAL_BCF.csi"
    if [[ -s "$RAM_BCF.clade_audit.tsv" ]]; then
        cp -f "$RAM_BCF.clade_audit.tsv" "$FINAL_BCF.clade_audit.tsv"
    fi
    [[ -s "$FINAL_BCF" && -s "$FINAL_BCF.csi" ]] || {{ echo "ERROR: failed to publish panel/index: $FINAL_BCF" >&2; exit 1; }}
    FINAL_COUNT="$(bcftools index -n "$FINAL_BCF")"
    [[ "$FINAL_COUNT" == "$COUNT" ]] || {{ echo "ERROR: published panel count mismatch: $FINAL_BCF" >&2; exit 1; }}
    echo "Panel variants: $COUNT  $FINAL_BCF"
done
{sidecar_publish_block}

TEMP_REQUEST={shlex.quote(str(final_request_file))}.partial.${{SLURM_JOB_ID}}
cp -f {shlex.quote(str(request_file))} "$TEMP_REQUEST"
cmp -s {shlex.quote(str(request_file))} "$TEMP_REQUEST" || {{ echo "ERROR: temporary panel request contract differs after copy" >&2; exit 1; }}
mv -f "$TEMP_REQUEST" {shlex.quote(str(final_request_file))}

echo "PASS: VCF downsampling completed"
"""
    write_text(script, body, executable=True)
    bash_syntax_check(script)
    return script, outputs


def plot_inputs(
    args: argparse.Namespace,
    run_dir: Path,
    modalities: Sequence[str],
) -> tuple[tuple[str, Path], ...]:
    if "atac" not in modalities:
        raise PipelineError("ATAC/RNA panel plots require --modalities to include atac")
    outputs = panel_outputs(args, run_dir, modalities)
    inputs: list[tuple[str, Path]] = [
        ("rna_demux", Path(args.rna_demux_panel)),
        ("rna_het", Path(args.rna_het_panel)),
        ("rna_species", Path(args.rna_species_panel)),
        ("atac_demux", outputs["atac_demux"]),
        ("atac_het", outputs["atac_het"]),
    ]
    species_modality = args.species_coverage
    if species_modality == "auto":
        species_modality = "rna" if "rna" in modalities else "atac"
    if species_modality == "atac":
        inputs.append(("atac_species", outputs["species"]))
    diagnostics_expected = True
    if args.stage == "PLOTS_ONLY":
        selector_request = run_dir / "panels" / "downsample.request.tsv"
        if selector_request.is_file() and selector_request.stat().st_size > 0:
            source_mode = read_key_value_metadata(selector_request).get("source_mode")
            if source_mode == "saved_atac_candidate_pool":
                diagnostics_expected = False
            elif source_mode == "full_source_bcf":
                diagnostics_expected = True
    if diagnostics_expected and "atac_numt" in outputs:
        inputs.append(("atac_numt", outputs["atac_numt"]))
    if diagnostics_expected and "mt_fusion_ratio" in outputs:
        inputs.append(("atac_mito", outputs["mt_fusion_ratio"]))
    return tuple(inputs)


def validate_plot_inputs(
    args: argparse.Namespace,
    run_dir: Path,
    modalities: Sequence[str],
) -> None:
    validate_static_plot_inputs(args)
    for label, path in plot_inputs(args, run_dir, modalities):
        require_nonempty(path, f"{label} panel")
        require_nonempty(Path(str(path) + ".csi"), f"{label} panel CSI")


def validate_static_plot_inputs(args: argparse.Namespace) -> None:
    """Fail before submission on plot resources that already must exist."""
    require_nonempty(Path(args.panel_qc_script), "panel QC plotting script")
    require_nonempty(Path(args.panel_metadata), "panel metadata TSV")
    require_nonempty(Path(args.gtf), "GTF")
    require_nonempty(Path(args.numt_bed), "NUMT BED")
    if args.atac_regulatory_bed:
        require_nonempty(Path(args.atac_regulatory_bed), "ATAC regulatory BED")
    for label, panel_text in (
        ("RNA demux plot reference", args.rna_demux_panel),
        ("RNA het plot reference", args.rna_het_panel),
        ("RNA species plot reference", args.rna_species_panel),
    ):
        require_nonempty(Path(panel_text), label)
        require_nonempty(Path(panel_text + ".csi"), label + " CSI")


def render_plot_script(
    args: argparse.Namespace,
    run_dir: Path,
    modalities: Sequence[str],
) -> tuple[Path, dict[str, Path]]:
    script = run_dir / "control" / "slurm" / "panel_qc_plots.sbatch"
    logs = run_dir / "logs"
    plot_dir = Path(args.figures_root) / args.run_label
    prefix = plot_dir / "atac_vs_rna_panel_qc"
    cache = Path(str(prefix) + "_cache.npz")
    marker = plot_dir / "plots.complete.tsv"
    inputs = plot_inputs(args, run_dir, modalities)
    header = sbatch_header(
        job_name=f"vplot_{args.run_label}",
        partition=args.partition,
        threads=args.plot_threads,
        memory=args.plot_memory,
        walltime=args.plot_time,
        stdout_path=logs / "panel_qc_plots_%j.out",
        stderr_path=logs / "panel_qc_plots_%j.err",
    )
    input_checks = "\n".join(
        f'[[ -s {shlex.quote(str(path))} && -s {shlex.quote(str(path) + ".csi")} ]] || '
        f'{{ echo "ERROR: panel/index missing for {label}: {path}" >&2; exit 1; }}'
        for label, path in inputs
    )
    input_args = "\n".join(
        f"    --input {shlex.quote(label)} {shlex.quote(str(path))}"
        for label, path in inputs
    )
    regulatory_plot_arg = (
        f"\n    --atac-regulatory-bed {shlex.quote(args.atac_regulatory_bed)}"
        if args.atac_regulatory_bed else ""
    )
    if args.atac_regulatory_required:
        regulatory_plot_arg += "\n    --atac-regulatory-required"
    expected = {
        "plot_cache": cache,
        "plot_summary": Path(str(prefix) + "_summary.tsv"),
        "plot_report": Path(str(prefix) + "_report.txt"),
        "plot_region_policy": Path(str(prefix) + "_region_policy.tsv"),
        "plot_overlap_table": Path(str(prefix) + "_rna_atac_overlap.tsv"),
        "plot_disjointness": Path(str(prefix) + "_rna_atac_disjointness.tsv"),
        "plot_overlap_figure": Path(str(prefix) + "_rna_atac_overlap.png"),
        "plot_structure": Path(str(prefix) + "_structure.png"),
        "plot_sfs": Path(str(prefix) + "_sfs.png"),
        "plot_ibs": Path(str(prefix) + "_ibs.png"),
        "plot_species_enrichment": Path(str(prefix) + "_species_enrichment.png"),
        "plot_selection_scores": Path(str(prefix) + "_selection_scores.png"),
        "plot_annotation": Path(str(prefix) + "_annotation.png"),
        "plot_structure_with_orang": plot_dir / "with_orang" / "atac_vs_rna_panel_qc_structure.png",
        "plot_ibs_with_orang": plot_dir / "with_orang" / "atac_vs_rna_panel_qc_ibs.png",
        "plot_individuals_with_orang": plot_dir / "with_orang" / "atac_vs_rna_panel_qc_individuals_atac_demux.png",
        "plot_completion": marker,
    }
    with_orang_keys = {
        "plot_structure_with_orang",
        "plot_ibs_with_orang",
        "plot_individuals_with_orang",
    }
    expected_checks = "\n".join(
        f'[[ -s {shlex.quote(str(path))} ]] || '
        f'{{ echo "ERROR: expected panel-QC output is missing: {path}" >&2; exit 1; }}'
        for key, path in expected.items()
        if key != "plot_completion" and key not in with_orang_keys
    )
    with_orang_checks = "\n".join(
        f'[[ -s {shlex.quote(str(path))} ]] || '
        f'{{ echo "ERROR: expected with-orangutan panel-QC output is missing: {path}" >&2; exit 1; }}'
        for key, path in expected.items()
        if key in with_orang_keys
    )
    request_file = run_dir / "control" / "tasks" / "panel_qc_plots.request.tsv"
    final_request_file = plot_dir / "panel_qc_plots.request.tsv"
    request_rows: list[tuple[object, object]] = [
        ("schema", "cellbouncer_vcf_panel_plot_request_v3"),
        ("panel_qc_script", args.panel_qc_script),
        ("panel_metadata", args.panel_metadata),
        ("gtf", args.gtf),
        ("numt_bed", args.numt_bed),
        ("atac_regulatory_bed", args.atac_regulatory_bed or "GTF_2KB_UPSTREAM_PROMOTERS"),
        ("atac_regulatory_required", str(args.atac_regulatory_required).lower()),
        ("plot_outgroup_sample", args.plot_outgroup_sample),
        ("rna_atac_overlap_policy", "overlap_allowed_rna_demux_soft_penalty_only"),
        ("rna_demux_soft_overlap_multiplier", args.atac_soft_overlap_multiplier),
        ("rna_het_panel_selection_role", "not_used_plot_reference_only"),
        ("rna_species_panel_selection_role", "not_used_plot_reference_only"),
    ]
    for label, path in inputs:
        request_rows.append((f"input_{label}", path))
    for label, path in expected.items():
        request_rows.append((f"output_{label}", path))
    plot_static_inputs = [
        ("panel_qc_script", args.panel_qc_script),
        ("panel_metadata", args.panel_metadata),
        ("gtf", args.gtf),
        ("numt_bed", args.numt_bed),
        ("rna_demux_panel", args.rna_demux_panel),
        ("rna_demux_panel_index", args.rna_demux_panel + ".csi"),
        ("rna_het_panel", args.rna_het_panel),
        ("rna_het_panel_index", args.rna_het_panel + ".csi"),
        ("rna_species_panel", args.rna_species_panel),
        ("rna_species_panel_index", args.rna_species_panel + ".csi"),
    ]
    if args.atac_regulatory_bed:
        plot_static_inputs.append(("atac_regulatory_bed", args.atac_regulatory_bed))
    for key, path_text in plot_static_inputs:
        file_stat = Path(path_text).stat()
        request_rows.extend(
            (
                (f"{key}_size", file_stat.st_size),
                (f"{key}_mtime_ns", file_stat.st_mtime_ns),
            )
        )
    write_tsv(request_file, ("key", "value"), request_rows)
    plot_upstream_inputs = [
        *(str(path) for _, path in inputs),
        *(str(path) + ".csi" for _, path in inputs),
        *(path for _, path in plot_static_inputs),
    ]
    plot_upstream_array = shell_join(plot_upstream_inputs)
    panel_input_array = shell_join(str(path) for _, path in inputs)
    plot_artifact_array = shell_join(
        (
            *(str(path) for key, path in expected.items() if key != "plot_completion"),
            str(final_request_file),
        )
    )
    obsolete_plot_artifact_array = shell_join(
        (
            str(prefix) + "_pool_identity_qc.tsv",
            str(prefix) + "_pool_donor_pair_qc.tsv",
            str(prefix) + "_selector_capacity_audit.tsv",
            str(prefix) + "_pool_operational_qc.png",
        )
    )
    body = f"""{header}

set -euo pipefail

{module_block(include_bcftools=True, include_genomics_python=True)}

for executable in bcftools cmp cp grep mv python3 rm; do
    command -v "$executable" >/dev/null 2>&1 || {{ echo "ERROR: missing command: $executable" >&2; exit 1; }}
done
[[ -s {shlex.quote(args.panel_qc_script)} ]] || {{ echo "ERROR: panel QC script is missing: {args.panel_qc_script}" >&2; exit 1; }}
[[ -s {shlex.quote(args.panel_metadata)} ]] || {{ echo "ERROR: panel metadata is missing: {args.panel_metadata}" >&2; exit 1; }}
[[ -s {shlex.quote(args.gtf)} ]] || {{ echo "ERROR: GTF is missing: {args.gtf}" >&2; exit 1; }}
[[ -s {shlex.quote(args.numt_bed)} ]] || {{ echo "ERROR: NUMT BED is missing: {args.numt_bed}" >&2; exit 1; }}
{input_checks}

PANEL_INPUTS=({panel_input_array})
OUTGROUP_PRESENT=false
for PANEL in "${{PANEL_INPUTS[@]}}"; do
    if bcftools query -l "$PANEL" | grep -Fqx -- {shlex.quote(args.plot_outgroup_sample)}; then
        OUTGROUP_PRESENT=true
        break
    fi
done

PLOT_REBUILD={str(args.force).lower()}
if [[ ! -s {shlex.quote(str(marker))} && {str(not args.force).lower()} == true ]]; then
    PLOT_EXISTING_ARTIFACTS=({plot_artifact_array})
    for ARTIFACT in "${{PLOT_EXISTING_ARTIFACTS[@]}}"; do
        if [[ -e "$ARTIFACT" ]]; then
            PLOT_REBUILD=true
            echo "REBUILD: plot completion marker is absent but prior artifacts exist"
            break
        fi
    done
fi
if [[ {str(not args.force).lower()} == true && -s {shlex.quote(str(marker))} ]]; then
    if [[ ! -s {shlex.quote(str(final_request_file))} ]] || ! cmp -s {shlex.quote(str(request_file))} {shlex.quote(str(final_request_file))}; then
        PLOT_REBUILD=true
        echo "REBUILD: plot request differs from the completed plot contract"
    else
        PLOT_UPSTREAM_INPUTS=({plot_upstream_array})
        for UPSTREAM in "${{PLOT_UPSTREAM_INPUTS[@]}}"; do
            if [[ "$UPSTREAM" -nt {shlex.quote(str(marker))} ]]; then
                PLOT_REBUILD=true
                echo "REBUILD: plot input is newer than completion marker: $UPSTREAM"
                break
            fi
        done
    fi
    if [[ "$PLOT_REBUILD" == false ]]; then
        echo "SKIP: completed ATAC/RNA panel plots match their request and inputs"
        exit 0
    fi
fi

mkdir -p {shlex.quote(str(plot_dir))}
if [[ "$PLOT_REBUILD" == true ]]; then
    rm -f {shlex.quote(str(marker))}
    OBSOLETE_PLOT_ARTIFACTS=({obsolete_plot_artifact_array})
    rm -f "${{OBSOLETE_PLOT_ARTIFACTS[@]}}"
fi
PLOT_ARGS=(
    {shlex.quote(args.panel_qc_script)}
{input_args}
    --panel-metadata {shlex.quote(args.panel_metadata)}
    --gtf {shlex.quote(args.gtf)}
    --numts-bed {shlex.quote(args.numt_bed)}
    --prefix {shlex.quote(str(prefix))}
    --cache {shlex.quote(str(cache))}
    --threads "$SLURM_CPUS_PER_TASK"
    --outgroup-sample {shlex.quote(args.plot_outgroup_sample)}
    --numt-panel-label atac_numt
    --mito-panel-label atac_mito{regulatory_plot_arg}
)
if [[ "$PLOT_REBUILD" == true ]]; then
    PLOT_ARGS+=(--force)
fi
printf 'Command:'
printf ' %q' "${{PLOT_ARGS[@]}}"
printf '\\n'
python3 "${{PLOT_ARGS[@]}}"

{expected_checks}
if [[ "$OUTGROUP_PRESENT" == true ]]; then
{with_orang_checks}
fi

TEMP_PLOT_REQUEST={shlex.quote(str(final_request_file))}.partial.${{SLURM_JOB_ID}}
cp -f {shlex.quote(str(request_file))} "$TEMP_PLOT_REQUEST"
cmp -s {shlex.quote(str(request_file))} "$TEMP_PLOT_REQUEST" || {{ echo "ERROR: temporary plot request contract differs after copy" >&2; exit 1; }}
mv -f "$TEMP_PLOT_REQUEST" {shlex.quote(str(final_request_file))}
{{
    printf 'schema\\tvcf_panel_plots_v3\\n'
    printf 'status\\tcomplete\\n'
    printf 'completed_utc\\t%s\\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
    printf 'rna_atac_overlap_policy\\toverlap_allowed_rna_demux_soft_penalty_only\\n'
    printf 'rna_demux_soft_overlap_multiplier\\t{args.atac_soft_overlap_multiplier}\\n'
    printf 'rna_demux_soft_overlap_panel\\t{args.rna_demux_panel}\\n'
    printf 'rna_het_panel_selection_role\\tnot_used_plot_reference_only\\n'
    printf 'rna_species_panel_selection_role\\tnot_used_plot_reference_only\\n'
    printf 'outgroup_sample\\t{args.plot_outgroup_sample}\\n'
    printf 'outgroup_present\\t%s\\n' "$OUTGROUP_PRESENT"
    printf 'request_contract\\t{final_request_file}\\n'
}} > {shlex.quote(str(marker))}.partial.${{SLURM_JOB_ID}}
mv -f {shlex.quote(str(marker))}.partial.${{SLURM_JOB_ID}} {shlex.quote(str(marker))}

echo "PASS: ATAC/RNA panel QC plots completed; RNA/ATAC overlap was measured, not forbidden"
echo "Plots: {plot_dir}"
"""
    write_text(script, body, executable=True)
    bash_syntax_check(script)
    return script, expected


def submit_script(sbatch: str, script: Path, dependencies: Sequence[str] = ()) -> str:
    command = [sbatch, "--parsable"]
    if dependencies:
        command.append("--dependency=afterok:" + ":".join(dependencies))
    command.append(str(script))
    result = subprocess.run(command, text=True, capture_output=True, check=False)
    if result.returncode:
        detail = result.stderr.strip() or result.stdout.strip() or "unknown sbatch failure"
        raise PipelineError(f"sbatch failed for {script}: {detail}")
    raw = result.stdout.strip()
    job_id = raw.split(";", 1)[0]
    if not re.fullmatch(r"\d+(?:_\[[^]]+\])?", job_id):
        raise PipelineError(f"unexpected sbatch --parsable response for {script}: {raw!r}")
    return job_id


def record_submissions(
    path: Path,
    rows: Sequence[tuple[str, str, Path, str]],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    needs_header = not path.exists() or path.stat().st_size == 0
    timestamp = dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    with path.open("a", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        if needs_header:
            writer.writerow(("submitted_utc", "job_kind", "job_id", "script", "dependency"))
        for job_kind, job_id, script, dependency in rows:
            writer.writerow((timestamp, job_kind, job_id, str(script), dependency or "NONE"))


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=(
            "Build cell/CB/noncell/BAM-wide coverage from each RNA and/or ATAC BAM, "
            "aggregate selected libraries, and create VCF panels with one selector pass."
        ),
        epilog=(
            "Population meanings: cells = CB values in the RNA filtered-barcode list; "
            "all_barcoded = every alignment with the configured tag; noncell = tagged "
            "alignments not in the cell list; bam_all = all accepted alignments. empty is "
            "created only from --empty-barcodes-template or --empty-drop-roster with an "
            "explicit true-empty list."
        ),
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {VERSION}")
    parser.add_argument("--stage", choices=STAGES, default="ALL")
    parser.add_argument("--run-label", required=True)
    parser.add_argument("--output-root", default=DEFAULT_OUTPUT_ROOT)
    parser.add_argument(
        "--figures-root", default=DEFAULT_FIGURE_ROOT,
        help=("Physical root for centralized VCF-panel figures; each run is "
              "written beneath <figures-root>/<run-label>."),
    )
    parser.add_argument(
        "--coverage-run-label",
        default=None,
        help=(
            "For ATAC_CANDIDATES_ONLY, DOWNSAMPLE_ONLY, DOWNSAMPLE_AND_PLOTS, "
            "or ATAC_RESELECT_AND_PLOTS, reuse aggregate tracks from this run label "
            "while writing panels under --run-label; useful for alternate "
            "cell/noncell/empty panels"
        ),
    )
    parser.add_argument(
        "--coverage-base-run-label",
        default=None,
        help=(
            "With INCREMENTAL_ALL, reuse completed per-library coverage from "
            "this run for every requested library except "
            "--coverage-rebuild-libraries. The old aggregate is never reused."
        ),
    )
    parser.add_argument(
        "--coverage-rebuild-libraries",
        nargs="+",
        default=None,
        metavar="LIB_OR_RANGE",
        help=(
            "With INCREMENTAL_ALL, rebuild coverage only for these libraries; "
            "the fresh aggregate still contains the complete --libraries set"
        ),
    )
    parser.add_argument(
        "--atac-candidate-run-label",
        default=None,
        help=(
            "With ATAC_RESELECT_AND_PLOTS, read the saved reusable ATAC candidate "
            "BCF from this run label; defaults to --run-label"
        ),
    )
    parser.add_argument(
        "--libraries",
        nargs="+",
        default=["1-40"],
        metavar="LIB_OR_RANGE",
        help="Library IDs/ranges; accepts '1-40', '1,3,7', or space-separated values",
    )
    parser.add_argument(
        "--modalities",
        default="rna,atac",
        help="Comma-separated selection: rna, atac, or rna,atac",
    )
    parser.add_argument("--submit", action="store_true", help="Submit the rendered job chain")
    parser.add_argument(
        "--force",
        action="store_true",
        help="Pass --force to coverage build/merge workers and rebuild requested panels",
    )

    paths = parser.add_argument_group("input paths and templates")
    paths.add_argument("--rna-root", default=DEFAULT_RNA_ROOT)
    paths.add_argument("--atac-root", default=DEFAULT_ATAC_ROOT)
    paths.add_argument("--rna-bam-template", default=DEFAULT_RNA_BAM_TEMPLATE)
    paths.add_argument("--atac-bam-template", default=DEFAULT_ATAC_BAM_TEMPLATE)
    paths.add_argument("--cell-barcodes-template", default=DEFAULT_CELL_BARCODES_TEMPLATE)
    empty_source = paths.add_mutually_exclusive_group()
    empty_source.add_argument(
        "--empty-barcodes-template",
        default=None,
        help=(
            "Explicit per-library true-empty barcode list. This is the only way to build/merge "
            "the empty population; never use Cell Ranger's generic non-filtered set here."
        ),
    )
    empty_source.add_argument(
        "--empty-drop-roster",
        default=None,
        help=(
            "Project-native TSV(.gz) with header library and one of barcode, "
            "cell_barcode, or CB; split into selected per-library true-empty lists"
        ),
    )
    paths.add_argument("--bam-window-coverage", default=DEFAULT_COVERAGE_BINARY)
    paths.add_argument("--downsample-vcf-parallel", default=DEFAULT_DOWNSAMPLE_BINARY)
    paths.add_argument("--source-bcf", default=DEFAULT_SOURCE_BCF)
    paths.add_argument("--source-bcf-index", default=DEFAULT_SOURCE_CSI)
    paths.add_argument("--gtf", default=DEFAULT_GTF)
    paths.add_argument("--numt-bed", default=DEFAULT_NUMT_BED)
    paths.add_argument(
        "--rna-aggregate-coverage",
        default=None,
        help=(
            "Explicit indexed RNA aggregate bedGraph.gz for DOWNSAMPLE_ONLY "
            "or DOWNSAMPLE_AND_PLOTS"
        ),
    )
    paths.add_argument(
        "--atac-aggregate-coverage",
        default=None,
        help=(
            "Explicit indexed ATAC aggregate bedGraph.gz for DOWNSAMPLE_ONLY "
            "or DOWNSAMPLE_AND_PLOTS"
        ),
    )
    paths.add_argument(
        "--blacklist",
        default=None,
        help=(
            "Optional mixed exclusion file: one field excludes a complete "
            "scaffold; BED3/BED3+ records exclude 0-based half-open intervals"
        ),
    )
    paths.add_argument("--panel-metadata", default=DEFAULT_PANEL_METADATA)
    paths.add_argument(
        "--pool-combinations",
        default=DEFAULT_POOL_COMBINATIONS,
        help=(
            "Explicit per-library biological identities derived from LibUID: "
            "diploid A, homotypic A+A, and heterotypic A+B are all retained."
        ),
    )
    paths.add_argument(
        "--rna-demux-panel",
        default=DEFAULT_RNA_DEMUX_PANEL,
        help=(
            "Exact selected production RNA demux 20M BCF. Only these selected "
            "loci receive the soft ATAC overlap ranking penalty."
        ),
    )
    paths.add_argument(
        "--rna-het-panel",
        default=DEFAULT_RNA_HET_PANEL,
        help="RNA het panel used for comparison plots only; never a selector penalty.",
    )
    paths.add_argument(
        "--rna-species-panel",
        default=DEFAULT_RNA_SPECIES_PANEL,
        help="RNA species panel used for comparison plots only; never a selector penalty.",
    )
    paths.add_argument("--panel-qc-script", default=DEFAULT_PANEL_QC_SCRIPT)
    paths.add_argument(
        "--atac-regulatory-bed",
        default=None,
        help=(
            "Optional regulatory BED used as a small ATAC ranking prior, never a "
            "hard filter. Without it, strand-aware 2-kb upstream promoters from "
            "--gtf are used."
        ),
    )

    coverage = parser.add_argument_group("coverage definition and resources")
    coverage.add_argument("--coverage-window-size", type=positive_int, default=100)
    coverage.add_argument(
        "--rna-min-mapq",
        type=nonnegative_int,
        default=20,
        help="RNA alignment MAPQ floor; set both modality floors to 0 for a legacy comparison",
    )
    coverage.add_argument(
        "--atac-min-mapq",
        type=nonnegative_int,
        default=30,
        help="ATAC alignment MAPQ floor; set both modality floors to 0 for a legacy comparison",
    )
    coverage.add_argument(
        "--coverage-exclude-flags",
        type=auto_base_int,
        default=0xF04,
        help="SAM flag mask, decimal or 0x-prefixed",
    )
    coverage.add_argument("--barcode-tag", default="CB")
    coverage.add_argument("--no-normalize-10x", action="store_true")
    coverage.add_argument(
        "--coverage-population", choices=POPULATIONS, default="cells",
        help="Default aggregate population used for panel ranking",
    )
    coverage.add_argument("--rna-coverage-population", choices=POPULATIONS, default=None)
    coverage.add_argument("--atac-coverage-population", choices=POPULATIONS, default=None)
    coverage.add_argument(
        "--coverage-concurrency",
        type=nonnegative_int,
        default=0,
        help=(
            "Optional maximum concurrent full-BAM scans; 0 leaves scheduling to SLURM. "
            "Set a cap only after measuring BeeGFS throughput in the pilot"
        ),
    )
    coverage.add_argument("--coverage-threads", type=positive_int, default=8)
    coverage.add_argument("--coverage-memory", default="32G")
    coverage.add_argument("--coverage-time", default="7-00:00:00")
    coverage.add_argument("--aggregate-threads", type=positive_int, default=8)
    coverage.add_argument("--aggregate-memory", default="64G")
    coverage.add_argument("--aggregate-time", default="7-00:00:00")

    selector = parser.add_argument_group("VCF selection and resources")
    selector.add_argument("--rna-demux-target", type=positive_int, default=20_000_000)
    selector.add_argument("--rna-het-target", type=positive_int, default=10_000_000)
    selector.add_argument("--atac-demux-target", type=positive_int, default=20_000_000)
    selector.add_argument("--atac-het-target", type=positive_int, default=10_000_000)
    selector.add_argument("--species-target", type=positive_int, default=20_000_000)
    selector.add_argument("--bin-size", type=positive_int, default=100_000)
    selector.add_argument("--min-coverage", type=nonnegative_float, default=1.0)
    selector.add_argument("--atac-min-coverage", type=nonnegative_float, default=1.0)
    selector.add_argument(
        "--atac-soft-overlap-multiplier",
        type=unit_interval_float,
        default=0.95,
        help=(
            "Score multiplier for ATAC loci also present in --rna-demux-panel. "
            "The default is a 5%% ranking penalty: overlap remains allowed. RNA het "
            "and species panels are plot references only and are never penalized."
        ),
    )
    selector.add_argument(
        "--atac-regulatory-required",
        action="store_true",
        help=(
            "Hard-restrict ATAC nuclear panels to the supplied regulatory BED or "
            "default 2-kb promoters. Off by default so coverage/open chromatin remains "
            "the primary eligibility rule."
        ),
    )
    selector.add_argument("--min-pair-score", type=nonnegative_float, default=0.3)
    selector.add_argument("--min-pairwise", type=nonnegative_int, default=1_500_000)
    selector.add_argument(
        "--species-outgroup",
        default="O",
        help="Species code treated as the singleton outgroup in ATAC species backfill",
    )
    selector.add_argument(
        "--species-outgroup-max-fraction",
        type=unit_interval_float,
        default=0.10,
        help=(
            "Maximum ATAC species-panel fraction assigned only to outgroup-involving "
            "contrasts; sites also helping non-outgroup pairs do not consume the cap"
        ),
    )
    selector.add_argument(
        "--max-het-per-bin",
        type=nonnegative_int,
        default=0,
        help="Explicit integer cap; 0 means unlimited (literal 'auto' is not accepted)",
    )
    selector.add_argument(
        "--species-max-per-bin",
        type=nonnegative_int,
        default=0,
        help="Explicit integer cap; 0 means unlimited (literal 'auto' is not accepted)",
    )
    selector.add_argument(
        "--species-coverage",
        choices=("auto", "rna", "atac"),
        default="auto",
        help="auto uses RNA when selected, otherwise ATAC",
    )
    selector.add_argument("--seed", type=int, default=42)
    selector.add_argument("--downsample-threads", type=positive_int, default=96)
    selector.add_argument("--downsample-memory", default="900G")
    selector.add_argument("--downsample-time", default="7-00:00:00")

    plotting = parser.add_argument_group("ATAC/RNA panel QC plots and resources")
    plotting.add_argument("--plot-threads", type=positive_int, default=48)
    plotting.add_argument("--plot-memory", default="400G")
    plotting.add_argument("--plot-time", default="7-00:00:00")
    plotting.add_argument(
        "--plot-outgroup-sample",
        default="JOS3C1",
        help=(
            "Sample excluded from the normal-scale sample-dependent figures; full "
            "versions are written under plots/with_orang"
        ),
    )

    mito = parser.add_argument_group("mitochondrial panel")
    mito.add_argument(
        "--mito-panel",
        action="store_true",
        help=(
            "Also create the mitochondrial fusion-ratio BCF and audit sidecars in the same "
            "selector pass. Full-source ATAC builds enable this automatically."
        ),
    )
    mito.add_argument("--mt-min-depth", type=positive_int, default=10)
    mito.add_argument("--mt-homoplasmy-af", type=unit_interval_float, default=0.95)
    mito.add_argument("--mt-min-coverage", type=nonnegative_float, default=1.0)
    mito.add_argument("--mt-min-pair-sites", type=positive_int, default=20)
    mito.add_argument("--mt-min-ambient-sites", type=nonnegative_int, default=1)
    mito.add_argument("--mt-require-pair-targets", action="store_true")

    slurm = parser.add_argument_group("SLURM")
    slurm.add_argument("--partition", default="compute")
    return parser


def validate_general_args(args: argparse.Namespace, modalities: Sequence[str]) -> Path:
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9._-]*", args.run_label):
        raise PipelineError(
            "--run-label must start with an alphanumeric character and contain only "
            "letters, digits, '.', '_', or '-'"
        )
    if args.coverage_run_label:
        if not re.fullmatch(
            r"[A-Za-z0-9][A-Za-z0-9._-]*", args.coverage_run_label
        ):
            raise PipelineError(
                "--coverage-run-label must start with an alphanumeric character and "
                "contain only letters, digits, '.', '_', or '-'"
            )
        if args.stage not in (
            "ATAC_CANDIDATES_ONLY",
            "DOWNSAMPLE_ONLY",
            "DOWNSAMPLE_AND_PLOTS",
            "ATAC_RESELECT_AND_PLOTS",
        ):
            raise PipelineError(
                "--coverage-run-label is only valid with --stage "
                "ATAC_CANDIDATES_ONLY, DOWNSAMPLE_ONLY, DOWNSAMPLE_AND_PLOTS, or "
                "ATAC_RESELECT_AND_PLOTS"
            )
    incremental_options_present = bool(
        args.coverage_base_run_label or args.coverage_rebuild_libraries
    )
    if args.stage == "INCREMENTAL_ALL":
        if not args.coverage_base_run_label or not args.coverage_rebuild_libraries:
            raise PipelineError(
                "INCREMENTAL_ALL requires both --coverage-base-run-label and "
                "--coverage-rebuild-libraries"
            )
        if not re.fullmatch(
            r"[A-Za-z0-9][A-Za-z0-9._-]*", args.coverage_base_run_label
        ):
            raise PipelineError(
                "--coverage-base-run-label must start with an alphanumeric "
                "character and contain only letters, digits, '.', '_', or '-'"
            )
        if args.coverage_base_run_label == args.run_label:
            raise PipelineError(
                "--coverage-base-run-label must differ from --run-label"
            )
        if args.coverage_run_label:
            raise PipelineError(
                "INCREMENTAL_ALL builds a fresh aggregate and cannot be combined "
                "with --coverage-run-label"
            )
        if args.empty_drop_roster:
            raise PipelineError(
                "INCREMENTAL_ALL does not reuse per-run split empty-drop rosters; "
                "use stable --empty-barcodes-template inputs instead"
            )
    elif incremental_options_present:
        raise PipelineError(
            "--coverage-base-run-label and --coverage-rebuild-libraries are only "
            "valid with --stage INCREMENTAL_ALL"
        )
    if args.atac_candidate_run_label:
        if not re.fullmatch(
            r"[A-Za-z0-9][A-Za-z0-9._-]*", args.atac_candidate_run_label
        ):
            raise PipelineError(
                "--atac-candidate-run-label must start with an alphanumeric "
                "character and contain only letters, digits, '.', '_', or '-'"
            )
        if args.stage != "ATAC_RESELECT_AND_PLOTS":
            raise PipelineError(
                "--atac-candidate-run-label is only valid with --stage "
                "ATAC_RESELECT_AND_PLOTS"
            )
    if args.stage in ("ATAC_CANDIDATES_ONLY", "ATAC_RESELECT_AND_PLOTS") and tuple(modalities) != ("atac",):
        raise PipelineError(
            f"{args.stage} requires --modalities atac"
        )
    if args.stage == "ATAC_RESELECT_AND_PLOTS" and args.mito_panel:
        raise PipelineError(
            "ATAC_RESELECT_AND_PLOTS reads a nuclear candidate pool and cannot "
            "recreate mitochondrial/NUMT diagnostics; use DOWNSAMPLE_AND_PLOTS "
            "with the full source BCF"
        )
    output_root = ensure_absolute(args.output_root, "--output-root")
    ensure_absolute(args.figures_root, "--figures-root")
    if any(character.isspace() for character in str(output_root)):
        raise PipelineError("--output-root cannot contain whitespace (SBATCH log-path constraint)")
    aggregate_overrides = (
        ("rna", "--rna-aggregate-coverage", args.rna_aggregate_coverage),
        ("atac", "--atac-aggregate-coverage", args.atac_aggregate_coverage),
    )
    supplied_overrides = [item for item in aggregate_overrides if item[2]]
    if supplied_overrides and args.stage not in (
        "DOWNSAMPLE_ONLY",
        "DOWNSAMPLE_AND_PLOTS",
    ):
        raise PipelineError(
            "--rna-aggregate-coverage and --atac-aggregate-coverage are only "
            "valid with --stage DOWNSAMPLE_ONLY or DOWNSAMPLE_AND_PLOTS"
        )
    for modality, option, value in supplied_overrides:
        ensure_absolute(value, option)
        if modality not in modalities:
            raise PipelineError(
                f"{option} requires --modalities to include {modality}"
            )
    for option, value in (
        ("--rna-root", args.rna_root),
        ("--atac-root", args.atac_root),
        ("--bam-window-coverage", args.bam_window_coverage),
        ("--downsample-vcf-parallel", args.downsample_vcf_parallel),
        ("--source-bcf", args.source_bcf),
        ("--source-bcf-index", args.source_bcf_index),
        ("--gtf", args.gtf),
        ("--numt-bed", args.numt_bed),
        ("--panel-metadata", args.panel_metadata),
        ("--pool-combinations", args.pool_combinations),
        ("--rna-demux-panel", args.rna_demux_panel),
        ("--rna-het-panel", args.rna_het_panel),
        ("--rna-species-panel", args.rna_species_panel),
        ("--panel-qc-script", args.panel_qc_script),
    ):
        ensure_absolute(value, option)
    if args.empty_drop_roster:
        ensure_absolute(args.empty_drop_roster, "--empty-drop-roster")
    if args.atac_regulatory_bed:
        ensure_absolute(args.atac_regulatory_bed, "--atac-regulatory-bed")
    if args.blacklist:
        ensure_absolute(args.blacklist, "--blacklist")
    if not re.fullmatch(r"[A-Za-z][A-Za-z0-9]", args.barcode_tag):
        raise PipelineError("--barcode-tag must be exactly two alphanumeric SAM-tag characters")
    if args.stage in (
        "ALL",
        "INCREMENTAL_ALL",
        "ATAC_CANDIDATES_ONLY",
        "DOWNSAMPLE_ONLY",
        "DOWNSAMPLE_AND_PLOTS",
        "ATAC_RESELECT_AND_PLOTS",
    ):
        species = args.species_coverage
        if species != "auto" and species not in modalities:
            raise PipelineError(
                f"--species-coverage {species} requires --modalities to include {species}"
            )
    selected = {selected_population(args, modality) for modality in modalities}
    if (
        "empty" in selected
        and args.stage not in (
            "ATAC_CANDIDATES_ONLY",
            "DOWNSAMPLE_ONLY",
            "DOWNSAMPLE_AND_PLOTS",
            "ATAC_RESELECT_AND_PLOTS",
        )
        and not has_empty_barcode_source(args)
    ):
        raise PipelineError(
            "the empty coverage population requires --empty-barcodes-template or "
            "--empty-drop-roster with an explicit true-empty barcode list"
        )
    return output_root / args.run_label


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        modalities = parse_modalities(args.modalities)
        libraries = parse_libraries(args.libraries)
        run_dir = validate_general_args(args, modalities)
        rebuild_libraries: tuple[int, ...] = ()
        base_coverage_run_dir: Path | None = None
        if args.stage == "INCREMENTAL_ALL":
            rebuild_libraries = parse_libraries(args.coverage_rebuild_libraries)
            unknown_rebuild = sorted(set(rebuild_libraries) - set(libraries))
            if unknown_rebuild:
                raise PipelineError(
                    "--coverage-rebuild-libraries must be a subset of --libraries; "
                    "outside the requested set: "
                    + ",".join(str(library) for library in unknown_rebuild)
                )
            if set(rebuild_libraries) == set(libraries):
                raise PipelineError(
                    "INCREMENTAL_ALL requires at least one reused library; use ALL "
                    "when rebuilding coverage for the complete --libraries set"
                )
            base_coverage_run_dir = run_dir.parent / args.coverage_base_run_label
        candidate_source_run_dir: Path | None = None
        if args.stage == "ATAC_RESELECT_AND_PLOTS":
            candidate_source_run_dir = run_dir.parent / (
                args.atac_candidate_run_label or args.run_label
            )
            candidate_source = atac_candidate_pool_path(
                args, candidate_source_run_dir
            )
            args.source_bcf = str(candidate_source)
            args.source_bcf_index = str(candidate_source) + ".csi"
        if args.coverage_run_label:
            coverage_run_dir = run_dir.parent / args.coverage_run_label
        elif candidate_source_run_dir is not None:
            coverage_run_dir = candidate_source_run_dir
        else:
            coverage_run_dir = run_dir
        for directory in (
            run_dir / "control" / "slurm",
            run_dir / "control" / "tasks",
            run_dir / "logs",
            run_dir / "coverage" / "per_library",
            run_dir / "coverage" / "aggregate",
            run_dir / "panels",
            Path(args.figures_root) / args.run_label,
        ):
            directory.mkdir(parents=True, exist_ok=True)

        blacklist_stages = {
            "ALL",
            "INCREMENTAL_ALL",
            "ATAC_CANDIDATES_ONLY",
            "DOWNSAMPLE_ONLY",
            "DOWNSAMPLE_AND_PLOTS",
            "ATAC_RESELECT_AND_PLOTS",
        }
        if args.stage in blacklist_stages:
            validate_and_stage_blacklist(args, run_dir)
        else:
            args.blacklist_source_path = None
            args.blacklist_staged_path = None
            args.blacklist_size_bytes = None
            args.blacklist_mtime_ns = None
            args.blacklist_sha256 = None

        rendered: list[Path] = []
        coverage_script: Path | None = None
        aggregate_scripts: dict[str, Path] = {}
        candidate_script: Path | None = None
        downsample_script: Path | None = None
        plot_script: Path | None = None
        outputs: dict[str, Path] = {}

        empty_barcode_paths = split_empty_drop_roster(
            args, run_dir, libraries
        )
        requested_task_rows = make_task_rows(
            args, run_dir, libraries, modalities, empty_barcode_paths
        )
        coverage_task_rows = requested_task_rows
        aggregate_task_rows = requested_task_rows
        reused_libraries: tuple[int, ...] = ()
        if base_coverage_run_dir is not None:
            (
                coverage_task_rows,
                aggregate_task_rows,
                reused_libraries,
            ) = resolve_incremental_coverage_rows(
                args,
                base_coverage_run_dir,
                requested_task_rows,
                rebuild_libraries,
            )

        if args.stage in ("ALL", "INCREMENTAL_ALL", "COVERAGE_ONLY"):
            validate_coverage_inputs(coverage_task_rows)
            coverage_catalog = (
                run_dir / "control" / "tasks" / "coverage_tasks.tsv"
            )
            if args.stage == "INCREMENTAL_ALL":
                write_tsv(
                    coverage_catalog,
                    TASK_COLUMNS,
                    aggregate_task_rows,
                )
                task_file = (
                    run_dir
                    / "control"
                    / "tasks"
                    / "coverage_rebuild_tasks.tsv"
                )
                write_tsv(task_file, TASK_COLUMNS, coverage_task_rows)
            else:
                task_file = coverage_catalog
                write_tsv(task_file, TASK_COLUMNS, coverage_task_rows)
            coverage_script = render_coverage_script(
                args, run_dir, task_file, len(coverage_task_rows)
            )
            rendered.append(coverage_script)

        if args.stage in ("ALL", "INCREMENTAL_ALL", "AGGREGATE_ONLY"):
            if args.stage == "AGGREGATE_ONLY":
                validate_per_library_products(args, aggregate_task_rows)
            input_lists = write_merge_input_lists(
                args, run_dir, aggregate_task_rows, modalities
            )
            populations = merge_populations(args)
            for modality in modalities:
                required_populations = {selected_population(args, modality)}
                species_modality = args.species_coverage
                if species_modality == "auto":
                    species_modality = "rna" if "rna" in modalities else "atac"
                if species_modality == modality:
                    required_populations.add(selected_population(args, modality))
                if args.mito_panel and modality == "rna":
                    required_populations.add("cells")
                aggregate_scripts[modality] = render_aggregate_script(
                    args,
                    run_dir,
                    modality,
                    input_lists[modality][0],
                    input_lists[modality][1],
                    populations,
                    tuple(sorted(required_populations)),
                )
                rendered.append(aggregate_scripts[modality])

        if args.stage == "ATAC_CANDIDATES_ONLY":
            validate_downsample_inputs(args, run_dir, modalities)
            validate_aggregate_products(
                args, coverage_run_dir, modalities, libraries
            )
            candidate_script, outputs = render_atac_candidate_script(
                args, run_dir, coverage_run_dir, libraries
            )
            rendered.append(candidate_script)

        if args.stage in (
            "ALL",
            "INCREMENTAL_ALL",
            "DOWNSAMPLE_ONLY",
            "DOWNSAMPLE_AND_PLOTS",
            "ATAC_RESELECT_AND_PLOTS",
        ):
            validate_downsample_inputs(args, run_dir, modalities)
            if args.stage in (
                "DOWNSAMPLE_ONLY",
                "DOWNSAMPLE_AND_PLOTS",
                "ATAC_RESELECT_AND_PLOTS",
            ):
                validate_aggregate_products(
                    args, coverage_run_dir, modalities, libraries
                )
            downsample_script, outputs = render_downsample_script(
                args, run_dir, modalities, coverage_run_dir, libraries
            )
            rendered.append(downsample_script)

        if (
            (args.stage in ("ALL", "INCREMENTAL_ALL") and "atac" in modalities)
            or args.stage in (
                "DOWNSAMPLE_AND_PLOTS",
                "ATAC_RESELECT_AND_PLOTS",
                "PLOTS_ONLY",
            )
        ):
            if args.stage == "PLOTS_ONLY":
                validate_plot_inputs(args, run_dir, modalities)
            else:
                validate_static_plot_inputs(args)
            plot_script, plot_outputs = render_plot_script(
                args, run_dir, modalities
            )
            outputs.update(plot_outputs)
            rendered.append(plot_script)

        print(f"orchestrate_vcf_panel_build.py {VERSION}")
        print(f"Stage: {args.stage}")
        print(f"Libraries: {','.join(str(lib) for lib in libraries)}")
        print(f"Modalities: {','.join(modalities)}")
        if base_coverage_run_dir is not None:
            print(f"Per-library coverage base directory: {base_coverage_run_dir}")
            print(
                "Coverage rebuild libraries: "
                + ",".join(str(library) for library in rebuild_libraries)
            )
            print(
                "Coverage reused libraries: "
                + ",".join(str(library) for library in reused_libraries)
            )
            print(f"Coverage tasks: {len(coverage_task_rows)}")
            print(f"Aggregate libraries: {len(libraries)}")
        if "atac" in modalities:
            print(
                "ATAC/RNA overlap policy: overlap allowed; RNA demux loci receive "
                f"a {args.atac_soft_overlap_multiplier:g} score multiplier"
            )
            print("RNA het/species panels: QC plot references only")
        print(f"Run directory: {run_dir}")
        if candidate_source_run_dir is not None:
            print(f"ATAC candidate source directory: {candidate_source_run_dir}")
        if coverage_run_dir != run_dir:
            print(f"Coverage source directory: {coverage_run_dir}")
        if args.stage in ("DOWNSAMPLE_ONLY", "DOWNSAMPLE_AND_PLOTS"):
            for modality in modalities:
                track = resolved_aggregate_track(
                    args, coverage_run_dir, modality
                )
                track_stat = track.stat()
                modified = dt.datetime.fromtimestamp(
                    track_stat.st_mtime, tz=dt.timezone.utc
                ).isoformat()
                print(
                    f"{modality.upper()} aggregate coverage: {track} "
                    f"(source={aggregate_coverage_source(args, modality)}, "
                    f"size={track_stat.st_size}, modified={modified})"
                )
        print("Rendered and bash -n validated:")
        for path in rendered:
            print(f"  {path}")

        if not args.submit:
            print("No jobs submitted. Re-run this command with --submit to submit the rendered stage.")
            if outputs:
                print("Requested panel outputs:")
                for name, path in outputs.items():
                    print(f"  {name}: {path}")
            return 0

        sbatch = shutil.which("sbatch")
        if not sbatch:
            raise PipelineError("--submit requested but sbatch is not available")
        sbatch = str(ensure_absolute(sbatch, "resolved sbatch executable"))
        submission_rows: list[tuple[str, str, Path, str]] = []
        coverage_job: str | None = None
        aggregate_jobs: dict[str, str] = {}
        candidate_job: str | None = None
        downsample_job: str | None = None
        plot_job: str | None = None

        if coverage_script is not None:
            coverage_job = submit_script(sbatch, coverage_script)
            submission_rows.append(("coverage_array", coverage_job, coverage_script, ""))
            print(f"Coverage array job ID: {coverage_job}")

        for modality in modalities:
            script = aggregate_scripts.get(modality)
            if script is None:
                continue
            dependencies = (coverage_job,) if coverage_job else ()
            job_id = submit_script(sbatch, script, dependencies)
            aggregate_jobs[modality] = job_id
            dependency = ":".join(dependencies)
            submission_rows.append((f"aggregate_{modality}", job_id, script, dependency))
            print(f"{modality.upper()} aggregate job ID: {job_id}")

        if candidate_script is not None:
            dependencies = tuple(
                aggregate_jobs[modality]
                for modality in modalities
                if modality in aggregate_jobs
            )
            candidate_job = submit_script(sbatch, candidate_script, dependencies)
            dependency = ":".join(dependencies)
            submission_rows.append(
                ("atac_candidates", candidate_job, candidate_script, dependency)
            )
            print(f"ATAC candidate-pool job ID: {candidate_job}")

        if downsample_script is not None:
            dependencies = tuple(aggregate_jobs[modality] for modality in modalities if modality in aggregate_jobs)
            downsample_job = submit_script(sbatch, downsample_script, dependencies)
            dependency = ":".join(dependencies)
            submission_rows.append(("downsample", downsample_job, downsample_script, dependency))
            print(f"Downsample job ID: {downsample_job}")

        if plot_script is not None:
            dependencies = (downsample_job,) if downsample_job else ()
            plot_job = submit_script(sbatch, plot_script, dependencies)
            dependency = ":".join(dependencies)
            submission_rows.append(("panel_qc_plots", plot_job, plot_script, dependency))
            print(f"Panel QC plot job ID: {plot_job}")

        record_submissions(run_dir / "control" / "submissions.tsv", submission_rows)
        print(f"Submission record: {run_dir / 'control' / 'submissions.tsv'}")
        print(f"Coverage outputs: {coverage_run_dir / 'coverage'}")
        if outputs:
            print("Panel outputs:")
            for name, path in outputs.items():
                print(f"  {name}: {path}")
        return 0
    except PipelineError as exc:
        parser.error(str(exc))
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
