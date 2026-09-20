#!/usr/bin/env python3
"""Cluster-side workers for the Tet 2025 RNA-anchored SCENIC+ analysis."""

from __future__ import annotations

import argparse
import collections
import csv
import gzip
import hashlib
import os
import pickle
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Iterable


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
MAO_GENOME = "human_chimp_bonobo"
HUMAN_GENOME = "Human"
PROBLEM_CONTIG = "human_chimp_bonoborefChr1218"
PRIMARY_HUMAN_RE = re.compile(r"^chr(?:[1-9]|1[0-9]|2[0-2]|X|Y)$")
HG38_SENTINEL_LENGTHS = {
    "chr1": 248_956_422,
    "chr2": 242_193_529,
    "chrX": 156_040_895,
}
RESOURCE_URLS = {
    "hg38_screen_v10_clust.regions_vs_motifs.rankings.feather": (
        "https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/"
        "screen/mc_v10_clust/region_based/"
        "hg38_screen_v10_clust.regions_vs_motifs.rankings.feather"
    ),
    "hg38_screen_v10_clust.regions_vs_motifs.scores.feather": (
        "https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/"
        "screen/mc_v10_clust/region_based/"
        "hg38_screen_v10_clust.regions_vs_motifs.scores.feather"
    ),
    "motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl": (
        "https://resources.aertslab.org/cistarget/motif2tf/"
        "motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl"
    ),
}
RESOURCE_SHA1_URLS = {
    "hg38_screen_v10_clust.regions_vs_motifs.rankings.feather": (
        RESOURCE_URLS["hg38_screen_v10_clust.regions_vs_motifs.rankings.feather"]
        + ".sha1sum.txt"
    ),
    "hg38_screen_v10_clust.regions_vs_motifs.scores.feather": (
        RESOURCE_URLS["hg38_screen_v10_clust.regions_vs_motifs.scores.feather"]
        + ".sha1sum.txt"
    ),
}


def log(message: str) -> None:
    print(message, flush=True)


def require_file(value: str | Path, description: str) -> Path:
    path = Path(value).expanduser().resolve()
    if not path.is_file() or path.stat().st_size == 0:
        raise FileNotFoundError(f"{description} is missing or empty: {path}")
    return path


def mkdir(path: str | Path) -> Path:
    result = Path(path).expanduser().resolve()
    result.mkdir(parents=True, exist_ok=True)
    return result


def parallel_temp_dir(root: str | Path) -> Path:
    """Return a deliberately short directory for Ray/joblib internals."""
    result = mkdir(Path(root) / "r")
    length = len(os.fsencode(str(result)))
    if length > 32:
        raise ValueError(
            "parallel temporary directory is too long for Ray's Unix sockets "
            f"({length} bytes; maximum supported here is 32): {result}"
        )
    return result


def run(command: list[str], *, cwd: Path | None = None, capture: bool = False):
    log("RUN: " + " ".join(map(str, command)))
    return subprocess.run(
        list(map(str, command)),
        cwd=None if cwd is None else str(cwd),
        check=True,
        text=True,
        stdout=subprocess.PIPE if capture else None,
        stderr=subprocess.PIPE if capture else None,
    )


def atomic_pickle(value, output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_name(output.name + f".tmp.{os.getpid()}")
    try:
        with open(temporary, "wb") as handle:
            pickle.dump(value, handle, protocol=pickle.HIGHEST_PROTOCOL)
        os.replace(temporary, output)
    finally:
        temporary.unlink(missing_ok=True)


def atomic_table(frame, output: Path, *, compression: str | None = None) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_name(output.name + f".tmp.{os.getpid()}")
    try:
        frame.to_csv(
            temporary,
            sep="\t",
            index=False,
            compression=compression,
        )
        os.replace(temporary, output)
    finally:
        temporary.unlink(missing_ok=True)


def atomic_bed(rows: Iterable[tuple], output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_name(output.name + f".tmp.{os.getpid()}")
    try:
        with open(temporary, "w", encoding="utf-8", newline="") as handle:
            for row in rows:
                handle.write("\t".join(map(str, row)) + "\n")
        os.replace(temporary, output)
    finally:
        temporary.unlink(missing_ok=True)


def read_hal_genomes(hal: Path) -> set[str]:
    result = run(["halStats", str(hal), "--genomes"], capture=True)
    genomes = set(result.stdout.split())
    if not genomes:
        raise RuntimeError(f"halStats returned no genomes for {hal}")
    return genomes


def read_hal_lengths(hal: Path, genome: str) -> dict[str, int]:
    result = run(
        ["halStats", str(hal), "--sequenceStats", genome], capture=True
    )
    lengths: dict[str, int] = {}
    for row in csv.reader(result.stdout.splitlines()):
        if not row or row[0] == "SequenceName":
            continue
        if len(row) < 2:
            raise ValueError(f"malformed halStats sequence row for {genome}")
        name, length = row[0], int(row[1])
        if not name or length <= 0 or name in lengths:
            raise ValueError(f"invalid HAL sequence for {genome}: {name}")
        lengths[name] = length
    if not lengths:
        raise RuntimeError(f"halStats returned no sequences for {genome}")
    return lengths


def read_fasta_lengths(path: Path) -> dict[str, int]:
    fai_candidates = [Path(str(path) + ".fai"), path.with_suffix(".fai")]
    for fai in fai_candidates:
        if fai.is_file() and fai.stat().st_size:
            result: dict[str, int] = {}
            with open(fai, "r", encoding="utf-8") as handle:
                for line in handle:
                    fields = line.rstrip("\n").split("\t")
                    if len(fields) >= 2:
                        result[fields[0]] = int(fields[1])
            if result:
                log(f"Read FASTA lengths from {fai}")
                return result
    log(f"No FASTA index found; reading sequence lengths from {path}")
    opener = gzip.open if path.suffix == ".gz" else open
    result: dict[str, int] = {}
    name: str | None = None
    length = 0
    with opener(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(">"):
                if name is not None:
                    result[name] = length
                name = line[1:].split()[0]
                length = 0
            else:
                length += len(line.strip())
    if name is not None:
        result[name] = length
    if not result:
        raise ValueError(f"FASTA contains no sequences: {path}")
    return result


def validate_mao_fasta(
    hal_lengths: dict[str, int], fasta_lengths: dict[str, int]
) -> dict[str, int]:
    shared = sorted(set(fasta_lengths).intersection(hal_lengths))
    if not shared:
        raise ValueError(
            f"Mao FASTA and HAL genome {MAO_GENOME} have no shared contigs"
        )
    mismatched = sorted(
        name
        for name in shared
        if fasta_lengths[name] != hal_lengths[name]
    )
    if mismatched:
        raise ValueError(
            "Mao FASTA is not coordinate-compatible with the HAL genome "
            f"{MAO_GENOME}; length mismatch: "
            + ", ".join(
                f"{name}({fasta_lengths[name]}!={hal_lengths[name]})"
                for name in mismatched[:5]
            )
        )

    fasta_only = sorted(set(fasta_lengths).difference(hal_lengths))
    if fasta_only:
        log(
            "Mapping-reference contigs absent from the HAL will be excluded "
            "from HAL-dependent analysis: " + ", ".join(fasta_only)
        )
    hal_only = sorted(set(hal_lengths).difference(fasta_lengths))
    if hal_only:
        log(
            "HAL contigs absent from the mapping FASTA will be excluded from "
            "fragment and peak analysis: " + ", ".join(hal_only)
        )
    return {name: hal_lengths[name] for name in shared}


def parse_gtf_attributes(text: str) -> dict[str, str]:
    result: dict[str, str] = {}
    for key, quoted, bare in re.findall(
        r"([A-Za-z0-9_.:-]+)\s+(?:\"([^\"]*)\"|([^;\s]+))\s*;?", text
    ):
        result[key] = quoted or bare
    return result


def read_protein_coding_genes(gtf: Path) -> list[dict[str, object]]:
    opener = gzip.open if gtf.suffix == ".gz" else open
    records: list[dict[str, object]] = []
    with opener(gtf, "rt", encoding="utf-8") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "gene":
                continue
            attrs = parse_gtf_attributes(fields[8])
            gene = attrs.get("gene_name") or attrs.get("gene_id")
            gene_type = (
                attrs.get("gene_type")
                or attrs.get("gene_biotype")
                or attrs.get("transcript_type")
                or ""
            )
            if not gene or gene_type != "protein_coding":
                continue
            start = int(fields[3]) - 1
            end = int(fields[4])
            row = {
                "source_chrom": fields[0],
                "source_start": start,
                "source_end": end,
                "source_strand": fields[6],
                "gene": gene,
                "gene_type": gene_type,
            }
            records.append(row)
    if not records:
        raise ValueError(f"no protein-coding gene records found in {gtf}")
    return records


def read_lifted_bed(path: Path) -> dict[str, list[dict[str, object]]]:
    grouped: dict[str, list[dict[str, object]]] = collections.defaultdict(list)
    if not path.exists():
        return grouped
    with open(path, "r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 4:
                raise ValueError(f"{path}:{line_number}: expected BED4 or BED6")
            start, end = int(fields[1]), int(fields[2])
            if start < 0 or end <= start:
                raise ValueError(f"{path}:{line_number}: invalid interval")
            grouped[fields[3]].append(
                {
                    "chrom": fields[0],
                    "start": start,
                    "end": end,
                    "name": fields[3],
                    "strand": fields[5] if len(fields) >= 6 else ".",
                }
            )
    return grouped


def hal_liftover(
    hal: Path,
    source_genome: str,
    source_bed: Path,
    target_genome: str,
    target_bed: Path,
) -> None:
    target_bed.unlink(missing_ok=True)
    run(
        [
            "halLiftover",
            str(hal),
            source_genome,
            str(source_bed),
            target_genome,
            str(target_bed),
        ]
    )
    if not target_bed.is_file() or target_bed.stat().st_size == 0:
        raise RuntimeError(
            f"halLiftover returned no intervals for {source_genome} -> {target_genome}"
        )


def flip_strand(strand: str) -> str:
    if strand == "+":
        return "-"
    if strand == "-":
        return "+"
    raise ValueError(f"cannot flip invalid strand: {strand!r}")


def reference(args) -> int:
    import pandas as pd

    hal = require_file(args.hal, "authoritative HAL alignment")
    fasta = require_file(args.mao_fasta, "Mao reference FASTA")
    gtf = require_file(args.mao_gtf, "Mao gene annotation")
    outdir = mkdir(args.output_dir)

    genomes = read_hal_genomes(hal)
    required_genomes = {args.mao_genome, args.human_genome}
    if not required_genomes.issubset(genomes):
        raise ValueError(
            "HAL lacks required genomes: "
            + ", ".join(sorted(required_genomes.difference(genomes)))
        )
    mao_hal_lengths = read_hal_lengths(hal, args.mao_genome)
    human_lengths = read_hal_lengths(hal, args.human_genome)
    mao_lengths = validate_mao_fasta(
        mao_hal_lengths, read_fasta_lengths(fasta)
    )
    for chrom, expected in HG38_SENTINEL_LENGTHS.items():
        observed = human_lengths.get(chrom)
        if observed != expected:
            raise ValueError(
                f"HAL target {args.human_genome} is not hg38 at {chrom}: "
                f"observed {observed}, expected {expected}"
            )

    mao_chromsizes = pd.DataFrame(
        [
            (name, 0, length)
            for name, length in mao_lengths.items()
        ],
        columns=["Chromosome", "Start", "End"],
    ).sort_values(["Chromosome"], kind="stable")
    hg38_chromsizes = pd.DataFrame(
        [
            (name, 0, length)
            for name, length in human_lengths.items()
            if PRIMARY_HUMAN_RE.fullmatch(name)
        ],
        columns=["Chromosome", "Start", "End"],
    ).sort_values(
        "Chromosome",
        key=lambda series: series.map(
            lambda x: int(x[3:]) if x[3:].isdigit() else 23 + "XY".index(x[3:])
        ),
    )
    if len(hg38_chromsizes) != 24:
        raise ValueError(
            f"expected 24 hg38 primary chromosomes, found {len(hg38_chromsizes)}"
        )
    atomic_table(mao_chromsizes, outdir / "mao_chromsizes.tsv")
    atomic_table(hg38_chromsizes, outdir / "hg38_chromsizes.tsv")

    valid_genes = [
        row
        for row in read_protein_coding_genes(gtf)
        if row["source_chrom"] != args.exclude_contig
        and row["source_chrom"] in mao_lengths
        and int(row["source_end"]) <= mao_lengths[str(row["source_chrom"])]
    ]
    if not valid_genes:
        raise ValueError("no coordinate-valid protein-coding genes remain")
    selected_genes: dict[str, dict[str, object]] = {}
    for row in valid_genes:
        key = str(row["gene"])
        prior = selected_genes.get(key)
        width = int(row["source_end"]) - int(row["source_start"])
        prior_width = (
            -1
            if prior is None
            else int(prior["source_end"]) - int(prior["source_start"])
        )
        if width > prior_width:
            selected_genes[key] = row
    genes = sorted(selected_genes.values(), key=lambda row: str(row["gene"]))
    with tempfile.TemporaryDirectory(prefix="scenicplus_gene_liftover.") as temp:
        tempdir = Path(temp)
        endpoints = tempdir / "mao_gene_endpoints.bed"
        lifted = tempdir / "hg38_gene_endpoints.bed"
        rows = []
        for index, gene in enumerate(genes):
            base = f"gene_{index:07d}"
            left, right = base + "_L", base + "_R"
            rows.append(
                (
                    gene["source_chrom"],
                    gene["source_start"],
                    int(gene["source_start"]) + 1,
                    left,
                    0,
                    "+",
                )
            )
            rows.append(
                (
                    gene["source_chrom"],
                    int(gene["source_end"]) - 1,
                    gene["source_end"],
                    right,
                    0,
                    "+",
                )
            )
        atomic_bed(rows, endpoints)
        hal_liftover(
            hal, args.mao_genome, endpoints, args.human_genome, lifted
        )
        mappings = read_lifted_bed(lifted)

    annotations = []
    crosswalk = []
    for index, gene in enumerate(genes):
        left_name = f"gene_{index:07d}_L"
        right_name = f"gene_{index:07d}_R"
        lefts, rights = mappings.get(left_name, []), mappings.get(right_name, [])
        reason = "accepted"
        if len(lefts) != 1 or len(rights) != 1:
            reason = "endpoint_not_one_to_one"
        elif lefts[0]["chrom"] != rights[0]["chrom"]:
            reason = "endpoints_on_different_chromosomes"
        elif not PRIMARY_HUMAN_RE.fullmatch(str(lefts[0]["chrom"])):
            reason = "nonprimary_human_chromosome"
        elif lefts[0]["strand"] not in {"+", "-"} or (
            lefts[0]["strand"] != rights[0]["strand"]
        ):
            reason = "inconsistent_alignment_orientation"
        if reason == "accepted":
            target_start = min(int(lefts[0]["start"]), int(rights[0]["start"]))
            target_end = max(int(lefts[0]["end"]), int(rights[0]["end"]))
            if target_end <= target_start:
                reason = "invalid_target_span"
            else:
                strand = str(gene["source_strand"])
                if lefts[0]["strand"] == "-":
                    strand = flip_strand(strand)
                # SCENIC+ expects the TSS to equal Start on '+' genes and End
                # on '-' genes, matching the BioMart annotation it normally builds.
                tss = target_start if strand == "+" else target_end
                annotations.append(
                    {
                        "Chromosome": lefts[0]["chrom"],
                        "Start": target_start,
                        "End": target_end,
                        "Strand": strand,
                        "Gene": gene["gene"],
                        "Transcription_Start_Site": tss,
                        "Transcript_type": "protein_coding",
                    }
                )
        crosswalk.append(
            {
                **gene,
                "target_chrom": (
                    lefts[0]["chrom"] if len(lefts) == 1 else ""
                ),
                "target_start": (
                    min(int(lefts[0]["start"]), int(rights[0]["start"]))
                    if len(lefts) == 1 and len(rights) == 1
                    else ""
                ),
                "target_end": (
                    max(int(lefts[0]["end"]), int(rights[0]["end"]))
                    if len(lefts) == 1 and len(rights) == 1
                    else ""
                ),
                "status": reason,
            }
        )
    annotation = pd.DataFrame(annotations).drop_duplicates("Gene", keep="first")
    if len(annotation) < 10_000:
        raise ValueError(
            f"only {len(annotation):,} protein-coding genes mapped to hg38 primary "
            "chromosomes; refusing an incomplete annotation"
        )
    atomic_table(annotation, outdir / "genome_annotation.tsv")
    atomic_table(
        pd.DataFrame(crosswalk),
        outdir / "gene_liftover_crosswalk.tsv.gz",
        compression="gzip",
    )
    log(
        f"Reference preparation complete: {len(mao_chromsizes):,} Mao contigs; "
        f"{len(annotation):,} hg38 protein-coding genes"
    )
    return 0


def download_one(url: str, output: Path, minimum_bytes: int) -> None:
    import requests

    output.parent.mkdir(parents=True, exist_ok=True)
    if output.is_file() and output.stat().st_size >= minimum_bytes:
        log(f"Resource already present: {output}")
        return
    partial = output.with_name(output.name + ".part")
    existing = partial.stat().st_size if partial.exists() else 0
    headers = {"Range": f"bytes={existing}-"} if existing else {}
    log(f"Downloading {url} -> {output} (resume at {existing:,} bytes)")
    with requests.get(url, headers=headers, stream=True, timeout=(30, 600)) as response:
        response.raise_for_status()
        append = existing > 0 and response.status_code == 206
        mode = "ab" if append else "wb"
        with open(partial, mode) as handle:
            for chunk in response.iter_content(chunk_size=16 * 1024 * 1024):
                if chunk:
                    handle.write(chunk)
    if partial.stat().st_size < minimum_bytes:
        raise ValueError(
            f"download is unexpectedly small ({partial.stat().st_size:,} bytes): {url}"
        )
    os.replace(partial, output)


def resources(args) -> int:
    root = mkdir(args.output_dir)
    for name, url in RESOURCE_URLS.items():
        minimum = 1_000_000 if name.endswith(".tbl") else 1_000_000_000
        download_one(url, root / name, minimum)
    for name, checksum_url in RESOURCE_SHA1_URLS.items():
        checksum_path = root / f"{name}.sha1sum.txt"
        download_one(checksum_url, checksum_path, 20)
        expected = checksum_path.read_text(encoding="utf-8").split()[0].lower()
        if not re.fullmatch(r"[0-9a-f]{40}", expected):
            raise ValueError(f"invalid SHA-1 sidecar: {checksum_path}")
        digest = hashlib.sha1()
        with open(root / name, "rb") as handle:
            for chunk in iter(lambda: handle.read(16 * 1024 * 1024), b""):
                digest.update(chunk)
        observed = digest.hexdigest()
        if observed != expected:
            raise ValueError(
                f"SHA-1 mismatch for {root / name}: {observed} != {expected}"
            )
        log(f"Checksum verified: {name}")
    log(f"cisTarget resources ready under {root}")
    return 0


def read_cell_metadata(path: Path):
    import pandas as pd

    frame = pd.read_csv(path, sep="\t", low_memory=False)
    required = {"cell_id", "barcode", "sample_id", "peak_group", "rna_leiden_reference"}
    missing = required.difference(frame.columns)
    if missing:
        raise KeyError("cell metadata is missing: " + ", ".join(sorted(missing)))
    frame["cell_id"] = frame["cell_id"].astype(str)
    frame["barcode"] = frame["barcode"].astype(str)
    frame["sample_id"] = frame["sample_id"].astype(str)
    frame["peak_group"] = frame["peak_group"].astype(str)
    frame["rna_leiden_reference"] = frame["rna_leiden_reference"].astype(str)
    if frame["cell_id"].duplicated().any():
        raise ValueError("cell metadata has duplicate cell_id values")
    return frame


def read_fragment_manifest(path: Path) -> dict[str, str]:
    import pandas as pd

    frame = pd.read_csv(path, sep="\t", dtype=str)
    required = {"library", "fragments"}
    missing = required.difference(frame.columns)
    if missing:
        raise KeyError("library task table is missing: " + ", ".join(sorted(missing)))
    result = {}
    for row in frame.itertuples(index=False):
        library = str(getattr(row, "library"))
        fragments = require_file(getattr(row, "fragments"), f"{library} fragments")
        result[library] = str(fragments)
    return result


def pseudobulk_peaks(args) -> int:
    import pandas as pd
    from pycisTopic.pseudobulk_peak_calling import export_pseudobulk, peak_calling
    from scatac_fragment_tools.library.split.split_fragments_by_cell_type import (
        _santize_string_for_filename,
    )

    metadata = read_cell_metadata(require_file(args.cell_metadata, "cell metadata"))
    fragments = read_fragment_manifest(require_file(args.task_table, "library task table"))
    chromsizes = pd.read_csv(
        require_file(args.mao_chromsizes, "Mao chromosome sizes"), sep="\t"
    )
    counts = metadata.groupby("peak_group", observed=True).size().rename("n_cells")
    kept_groups = sorted(
        map(str, counts[counts >= args.min_cells_per_group].index.tolist())
    )
    filtered = metadata.loc[metadata["peak_group"].isin(kept_groups)].copy()
    if filtered.empty:
        raise ValueError("no RNA-cluster-by-library groups pass the cell threshold")
    missing_fragments = set(filtered["sample_id"]).difference(fragments)
    if missing_fragments:
        raise KeyError(
            "fragment paths are missing for: " + ", ".join(sorted(missing_fragments))
        )
    outdir = mkdir(args.output_dir)
    atomic_table(
        counts.reset_index().sort_values(["n_cells", "peak_group"], ascending=[False, True]),
        outdir / "peak_group_counts.tsv",
    )
    filtered = filtered.set_index("cell_id", drop=False)
    bed_root = mkdir(outdir / "pseudobulk_bed")
    bigwig_root = mkdir(outdir / "pseudobulk_bigwig")
    sanitized = {
        group: _santize_string_for_filename(group) for group in kept_groups
    }
    if len(set(sanitized.values())) != len(sanitized):
        raise ValueError("peak-group names collide after filename sanitization")

    expected_beds = {
        group: bed_root / f"{name}.fragments.tsv.gz"
        for group, name in sanitized.items()
    }
    expected_bigwigs = {
        group: bigwig_root / f"{name}.bw" for group, name in sanitized.items()
    }
    completed = all(
        path.is_file() and path.stat().st_size > 0
        for path in (*expected_beds.values(), *expected_bigwigs.values())
    )
    if completed:
        beds = {group: str(path) for group, path in expected_beds.items()}
        bigwigs = {
            group: str(path) for group, path in expected_bigwigs.items()
        }
        log(
            f"Reusing {len(beds):,} completed pseudobulk BED/bigWig pairs; "
            "fragment splitting and bigWig generation will not be repeated"
        )
    else:
        missing_beds = sum(
            not path.is_file() or path.stat().st_size == 0
            for path in expected_beds.values()
        )
        missing_bigwigs = sum(
            not path.is_file() or path.stat().st_size == 0
            for path in expected_bigwigs.values()
        )
        log(
            "Pseudobulk output set is incomplete; regenerating it "
            f"({missing_beds} BED and {missing_bigwigs} bigWig files missing/empty)"
        )
        bigwigs, beds = export_pseudobulk(
            input_data=filtered,
            variable="peak_group",
            chromsizes=chromsizes,
            bed_path=str(bed_root),
            bigwig_path=str(bigwig_root),
            path_to_fragments=fragments,
            sample_id_col="sample_id",
            n_cpu=args.cpus,
            normalize_bigwig=True,
            split_pattern="___",
            temp_dir=str(mkdir(Path(args.temp_dir) / "p")),
        )
        missing_after_export = [
            str(path)
            for path in (*expected_beds.values(), *expected_bigwigs.values())
            if not path.is_file() or path.stat().st_size == 0
        ]
        if missing_after_export:
            raise RuntimeError(
                "pseudobulk export did not create every expected nonempty file; "
                f"first missing: {missing_after_export[:3]}"
            )
        beds = {group: str(path) for group, path in expected_beds.items()}
        bigwigs = {
            group: str(path) for group, path in expected_bigwigs.items()
        }
    macs = shutil.which("macs2") or shutil.which("macs3")
    if not macs:
        raise FileNotFoundError("neither macs2 nor macs3 is available on PATH")
    genome_size = str(int(chromsizes["End"].sum()))
    narrow_peaks = peak_calling(
        macs_path=macs,
        bed_paths=beds,
        outdir=str(mkdir(outdir / "macs2")),
        genome_size=genome_size,
        n_cpu=args.cpus,
        input_format="BEDPE",
        shift=73,
        ext_size=146,
        keep_dup="all",
        q_value=args.q_value,
        nolambda=True,
        skip_empty_peaks=True,
        _temp_dir=str(parallel_temp_dir(args.temp_dir)),
    )
    if not narrow_peaks:
        raise RuntimeError("MACS produced no nonempty peak sets")
    atomic_pickle(narrow_peaks, outdir / "narrow_peaks.pkl")
    log(
        f"Pseudobulk peak calling complete: {len(filtered):,} cells, "
        f"{len(beds):,} RNA-cluster-by-library pseudobulks, "
        f"{len(narrow_peaks):,} nonempty peak sets, {len(bigwigs):,} bigWigs"
    )
    return 0


def interval_overlap(start_a: int, end_a: int, start_b: int, end_b: int) -> int:
    return max(0, min(end_a, end_b) - max(start_a, start_b))


def consensus_liftover(args) -> int:
    import pandas as pd
    from pycisTopic.iterative_peak_calling import get_consensus_peaks

    hal = require_file(args.hal, "authoritative HAL alignment")
    narrow_path = require_file(args.narrow_peaks, "pseudobulk narrow-peak pickle")
    chromsizes = pd.read_csv(
        require_file(args.mao_chromsizes, "Mao chromosome sizes"), sep="\t"
    )
    with open(narrow_path, "rb") as handle:
        narrow_peaks = pickle.load(handle)
    consensus = get_consensus_peaks(
        narrow_peaks_dict=narrow_peaks,
        peak_half_width=args.peak_half_width,
        chromsizes=chromsizes,
        path_to_blacklist=None,
    )
    source = consensus.df.loc[:, ["Chromosome", "Start", "End"]].copy()
    source["Chromosome"] = source["Chromosome"].astype(str)
    source = source.loc[source["Chromosome"] != args.exclude_contig].copy()
    source["Start"] = source["Start"].astype(int)
    source["End"] = source["End"].astype(int)
    source = source.loc[source["End"] > source["Start"]].copy()
    source = source.drop_duplicates(["Chromosome", "Start", "End"])
    source = source.sort_values(["Chromosome", "Start", "End"], kind="stable")
    source["peak_id"] = [f"peak_{index:08d}" for index in range(len(source))]
    if source.empty:
        raise RuntimeError("consensus peak creation produced no valid Mao intervals")

    outdir = mkdir(args.output_dir)
    all_source_bed = outdir / "consensus_peaks_mao_all.bed"
    forward_bed = outdir / "consensus_peaks_hg38_raw.bed"
    reciprocal_input = outdir / "consensus_peaks_hg38_reciprocal_input.bed"
    reciprocal_bed = outdir / "consensus_peaks_mao_reciprocal_raw.bed"
    atomic_bed(
        (
            (row.Chromosome, row.Start, row.End, row.peak_id, 0, "+")
            for row in source.itertuples(index=False)
        ),
        all_source_bed,
    )
    hal_liftover(
        hal,
        args.mao_genome,
        all_source_bed,
        args.human_genome,
        forward_bed,
    )
    forward = read_lifted_bed(forward_bed)
    source_by_id = {
        row.peak_id: row for row in source.itertuples(index=False)
    }

    candidates: dict[str, dict[str, object]] = {}
    decisions: dict[str, str] = {}
    for peak_id, row in source_by_id.items():
        mappings = forward.get(peak_id, [])
        if len(mappings) != 1:
            decisions[peak_id] = "forward_not_one_to_one"
            continue
        target = mappings[0]
        if not PRIMARY_HUMAN_RE.fullmatch(str(target["chrom"])):
            decisions[peak_id] = "nonprimary_human_chromosome"
            continue
        source_width = int(row.End) - int(row.Start)
        target_width = int(target["end"]) - int(target["start"])
        ratio = target_width / source_width
        if not args.min_length_ratio <= ratio <= args.max_length_ratio:
            decisions[peak_id] = "forward_length_ratio"
            continue
        candidates[peak_id] = target
        decisions[peak_id] = "awaiting_reciprocal_check"

    if not candidates:
        raise RuntimeError("no consensus peaks passed the forward Mao-to-hg38 lift")
    atomic_bed(
        (
            (
                row["chrom"],
                row["start"],
                row["end"],
                peak_id,
                0,
                "+",
            )
            for peak_id, row in candidates.items()
        ),
        reciprocal_input,
    )
    hal_liftover(
        hal,
        args.human_genome,
        reciprocal_input,
        args.mao_genome,
        reciprocal_bed,
    )
    reciprocal = read_lifted_bed(reciprocal_bed)
    accepted: dict[str, dict[str, object]] = {}
    for peak_id, target in candidates.items():
        mappings = reciprocal.get(peak_id, [])
        source_row = source_by_id[peak_id]
        if len(mappings) != 1:
            decisions[peak_id] = "reciprocal_not_one_to_one"
            continue
        returned = mappings[0]
        if returned["chrom"] != source_row.Chromosome:
            decisions[peak_id] = "reciprocal_different_mao_contig"
            continue
        source_width = int(source_row.End) - int(source_row.Start)
        overlap = interval_overlap(
            int(source_row.Start),
            int(source_row.End),
            int(returned["start"]),
            int(returned["end"]),
        )
        reciprocal_fraction = overlap / source_width
        if reciprocal_fraction < args.min_reciprocal_overlap:
            decisions[peak_id] = "insufficient_reciprocal_overlap"
            continue
        accepted[peak_id] = {
            **target,
            "reciprocal_overlap_fraction": reciprocal_fraction,
        }
        decisions[peak_id] = "accepted"

    target_coordinates: dict[tuple[str, int, int], list[str]] = collections.defaultdict(list)
    for peak_id, target in accepted.items():
        target_coordinates[
            (str(target["chrom"]), int(target["start"]), int(target["end"]))
        ].append(peak_id)
    for peak_ids in target_coordinates.values():
        if len(peak_ids) > 1:
            for peak_id in peak_ids:
                accepted.pop(peak_id, None)
                decisions[peak_id] = "duplicate_hg38_interval"

    records = []
    for peak_id, row in source_by_id.items():
        target = candidates.get(peak_id)
        accepted_target = accepted.get(peak_id)
        source_region = f"{row.Chromosome}:{int(row.Start)}-{int(row.End)}"
        target_region = ""
        if target is not None:
            target_region = (
                f"{target['chrom']}:{int(target['start'])}-{int(target['end'])}"
            )
        records.append(
            {
                "peak_id": peak_id,
                "mao_chrom": row.Chromosome,
                "mao_start": int(row.Start),
                "mao_end": int(row.End),
                "mao_region": source_region,
                "hg38_chrom": "" if target is None else target["chrom"],
                "hg38_start": "" if target is None else int(target["start"]),
                "hg38_end": "" if target is None else int(target["end"]),
                "hg38_region": target_region,
                "reciprocal_overlap_fraction": (
                    ""
                    if accepted_target is None
                    else accepted_target["reciprocal_overlap_fraction"]
                ),
                "status": decisions[peak_id],
            }
        )
    crosswalk = pd.DataFrame(records)
    retained = crosswalk.loc[crosswalk["status"] == "accepted"].copy()
    if len(retained) < args.min_retained_peaks:
        raise ValueError(
            f"only {len(retained):,} peaks passed reciprocal HAL projection; "
            f"minimum is {args.min_retained_peaks:,}"
        )
    atomic_table(
        crosswalk,
        outdir / "peak_liftover_crosswalk.tsv.gz",
        compression="gzip",
    )
    atomic_bed(
        retained[["mao_chrom", "mao_start", "mao_end"]].itertuples(
            index=False, name=None
        ),
        outdir / "consensus_peaks_mao_retained.bed",
    )
    target_sorted = retained.sort_values(
        ["hg38_chrom", "hg38_start", "hg38_end"], kind="stable"
    )
    atomic_bed(
        target_sorted[["hg38_chrom", "hg38_start", "hg38_end"]].itertuples(
            index=False, name=None
        ),
        outdir / "consensus_peaks_hg38_retained.bed",
    )
    log(
        f"Consensus/liftover complete: {len(source):,} Mao peaks; "
        f"{len(retained):,} unique reciprocal hg38 mappings "
        f"({100 * len(retained) / len(source):.2f}%)"
    )
    return 0


def build_cistopic_library(args) -> int:
    from pycisTopic.cistopic_class import create_cistopic_object_from_fragments

    metadata = read_cell_metadata(require_file(args.cell_metadata, "cell metadata"))
    selected = metadata.loc[metadata["sample_id"] == args.library].copy()
    if selected.empty:
        raise ValueError(f"cell metadata has no cells for {args.library}")
    fragments = require_file(args.fragments, f"{args.library} fragments")
    regions = require_file(args.regions, "retained Mao consensus peaks")
    output = Path(args.output).expanduser().resolve()
    obj = create_cistopic_object_from_fragments(
        path_to_fragments=str(fragments),
        path_to_regions=str(regions),
        path_to_blacklist=None,
        metrics=None,
        valid_bc=selected["barcode"].astype(str).tolist(),
        n_cpu=args.cpus,
        min_frag=1,
        min_cell=1,
        is_acc=1,
        check_for_duplicates=False,
        project=args.library,
        partition=args.partitions,
        split_pattern="___",
        use_polars=True,
    )
    if not obj.cell_names or not obj.region_names:
        raise RuntimeError(f"cisTopic object is empty for {args.library}")
    atomic_pickle(obj, output)
    log(
        f"{args.library}: cisTopic object contains {len(obj.cell_names):,} cells "
        f"and {len(obj.region_names):,} accessible regions"
    )
    return 0


def merge_cistopic(args) -> int:
    import pandas as pd
    from pycisTopic.cistopic_class import merge
    from pycisTopic.utils import region_names_to_coordinates

    task_table = pd.read_csv(
        require_file(args.task_table, "library task table"), sep="\t", dtype=str
    )
    if "cistopic_output" not in task_table.columns:
        raise KeyError("library task table lacks cistopic_output")
    objects = []
    for row in task_table.itertuples(index=False):
        path = require_file(row.cistopic_output, f"{row.library} cisTopic object")
        with open(path, "rb") as handle:
            objects.append(pickle.load(handle))
    if len(objects) < 2:
        raise ValueError("at least two library cisTopic objects are required")
    obj = merge(objects, project=args.project, split_pattern="___")

    crosswalk = pd.read_csv(
        require_file(args.peak_crosswalk, "peak liftover crosswalk"), sep="\t"
    )
    crosswalk = crosswalk.loc[crosswalk["status"] == "accepted"].copy()
    coordinate_map = dict(zip(crosswalk["mao_region"], crosswalk["hg38_region"]))
    missing = sorted(set(obj.region_names).difference(coordinate_map))
    if missing:
        raise ValueError(
            f"{len(missing):,} cisTopic regions are absent from the accepted "
            f"peak crosswalk; first: {missing[:3]}"
        )
    original_names = list(obj.region_names)
    target_names = [coordinate_map[name] for name in original_names]
    if len(set(target_names)) != len(target_names):
        raise ValueError("hg38 region names are not unique after coordinate translation")
    region_data = obj.region_data.loc[original_names].copy()
    region_data["Mao_region"] = original_names
    coordinates = region_names_to_coordinates(target_names)
    region_data.index = target_names
    for column in ("Chromosome", "Start", "End"):
        region_data[column] = coordinates.loc[target_names, column]
    region_data["Width"] = region_data["End"] - region_data["Start"]
    obj.region_names = target_names
    obj.region_data = region_data.loc[target_names]

    metadata = read_cell_metadata(require_file(args.cell_metadata, "cell metadata"))
    metadata = metadata.set_index("cell_id", drop=False)
    absent_cells = sorted(set(obj.cell_names).difference(metadata.index))
    if absent_cells:
        raise ValueError(
            f"{len(absent_cells):,} cisTopic cells lack RNA metadata; first: "
            f"{absent_cells[:3]}"
        )
    obj.add_cell_data(metadata, split_pattern="___")
    missing_groups = obj.cell_data["rna_leiden_reference"].isna().sum()
    if missing_groups:
        raise ValueError(f"{missing_groups:,} merged cells lack RNA Leiden labels")
    atomic_pickle(obj, Path(args.output).expanduser().resolve())
    log(
        f"Merged cisTopic object: {len(obj.cell_names):,} cells, "
        f"{len(obj.region_names):,} reciprocal hg38 regions"
    )
    return 0


def topic_model(args) -> int:
    from pycisTopic.lda_models import evaluate_models, run_cgs_models

    input_path = require_file(args.input, "merged cisTopic object")
    with open(input_path, "rb") as handle:
        obj = pickle.load(handle)
    output = Path(args.output).expanduser().resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    model_dir = mkdir(output.parent / "candidate_topic_models")
    models = run_cgs_models(
        obj,
        n_topics=list(args.topics),
        n_cpu=min(args.cpus, len(args.topics)),
        n_iter=args.iterations,
        random_state=args.seed,
        save_path=str(model_dir),
        _temp_dir=str(parallel_temp_dir(args.temp_dir)),
    )
    selected = evaluate_models(
        models,
        select_model=args.select_topics,
        return_model=True,
        plot=False,
        save=str(output.parent / "topic_model_selection.pdf"),
    )
    obj.add_LDA_model(selected)
    atomic_pickle(obj, output)
    log(f"Selected cisTopic model with {selected.n_topic} topics")
    return 0


def safe_filename(value: object) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value)).strip("_") or "set"


def write_region_sets(region_sets: dict, directory: Path) -> int:
    directory.mkdir(parents=True, exist_ok=True)
    for stale in directory.glob("*.bed"):
        stale.unlink()
    written = 0
    for name, frame in region_sets.items():
        if frame is None or frame.empty:
            continue
        rows = []
        for region in map(str, frame.index):
            match = re.fullmatch(r"([^:]+):(\d+)-(\d+)", region)
            if match is None:
                raise ValueError(f"invalid hg38 region name: {region}")
            rows.append((match.group(1), int(match.group(2)), int(match.group(3))))
        atomic_bed(rows, directory / f"{safe_filename(name)}.bed")
        written += 1
    return written


def region_sets(args) -> int:
    import numpy as np
    import ray
    from pycisTopic.diff_features import (
        find_highly_variable_features,
        impute_accessibility,
        markers,
        normalize_scores,
    )
    from pycisTopic.topic_binarization import binarize_topics

    input_path = require_file(args.input, "topic-modelled cisTopic object")
    with open(input_path, "rb") as handle:
        obj = pickle.load(handle)
    if not getattr(obj, "selected_model", None):
        raise ValueError("cisTopic object has no selected topic model")
    root = mkdir(args.output_dir)
    otsu = binarize_topics(obj, method="otsu", plot=False)
    top_n = min(args.top_regions_per_topic, max(1, len(obj.region_names) - 1))
    top = binarize_topics(obj, method="ntop", ntop=top_n, plot=False)
    n_otsu = write_region_sets(otsu, root / "Topics_otsu")
    n_top = write_region_sets(top, root / "Topics_top_3k")

    imputed = impute_accessibility(
        obj,
        selected_cells=None,
        selected_regions=None,
        scale_factor=10**6,
        chunk_size=args.imputation_chunk_size,
    )
    normalized = normalize_scores(imputed, scale_factor=10**4)
    n_variable = min(args.variable_regions, max(1, len(normalized.feature_names) - 1))
    variable = find_highly_variable_features(
        normalized,
        n_top_features=n_variable,
        plot=False,
    )
    group_var = obj.cell_data[args.grouping].dropna().astype(str)
    group_counts = group_var.value_counts().sort_index()
    levels = sorted(
        group_counts.loc[group_counts >= args.min_cells_per_dar].index.tolist()
    )
    if len(levels) < 2:
        raise ValueError(
            f"{args.grouping} must contain at least two groups with at least "
            f"{args.min_cells_per_dar} cells; found {levels}"
        )
    obj.cell_data[args.grouping] = obj.cell_data[args.grouping].astype(str)
    atomic_table(
        group_counts.rename("n_cells").reset_index().rename(
            columns={args.grouping: "group"}
        ),
        root / f"DAR_cell_counts_{safe_filename(args.grouping)}.tsv",
    )

    # pycisTopic's current find_diff_features implementation constructs an
    # incorrect foreground-versus-foreground comparison.  Call its tested
    # marker engine directly with the intended one-group-versus-rest cells.
    selected_imputed = imputed.subset(
        cells=None,
        features=variable,
        copy=True,
        split_pattern="___",
    )
    available = set(selected_imputed.cell_names)
    group_var = group_var.loc[group_var.index.map(lambda cell: cell in available)]
    marker_sets = {}
    ray_temp = parallel_temp_dir(args.temp_dir)
    ray.init(num_cpus=args.cpus, _temp_dir=str(ray_temp))
    try:
        for level in levels:
            foreground = group_var.index[group_var == level].tolist()
            background = group_var.index[group_var != level].tolist()
            if len(foreground) < args.min_cells_per_dar:
                continue
            if len(background) < args.min_cells_per_dar:
                raise ValueError(
                    f"RNA group {level!r} has only {len(background):,} background cells"
                )
            log(
                f"RNA-defined DARs for {level}: {len(foreground):,} cells versus "
                f"{len(background):,} remaining cells"
            )
            marker_sets[level] = markers(
                selected_imputed,
                [foreground, background],
                level,
                adjpval_thr=args.adjusted_pvalue,
                log2fc_thr=np.log2(args.fold_change),
                n_cpu=args.cpus,
            )
    finally:
        ray.shutdown()
    n_dar = write_region_sets(
        marker_sets, root / f"DARs_{safe_filename(args.grouping)}"
    )
    if n_otsu == 0 or n_top == 0 or n_dar == 0:
        raise RuntimeError(
            f"incomplete region sets: Otsu={n_otsu}, top={n_top}, DAR={n_dar}"
        )
    atomic_pickle(obj, Path(args.output_cistopic).expanduser().resolve())
    log(
        f"Region sets complete: {n_otsu} Otsu topics, {n_top} top-region "
        f"topics, {n_dar} RNA-defined DAR sets"
    )
    return 0


def configure_and_run_scenicplus(args) -> int:
    import gc

    import anndata as ad
    import pandas as pd
    import yaml

    cistopic = require_file(args.cistopic, "final cisTopic object")
    rna = require_file(args.rna_h5ad, "SCENIC+ RNA H5AD")
    region_root = Path(args.region_sets).expanduser().resolve()
    if not region_root.is_dir() or not any(region_root.glob("*/*.bed")):
        raise FileNotFoundError(f"region-set folder has no BED files: {region_root}")
    resources_root = Path(args.resources).expanduser().resolve()
    ranking = require_file(
        resources_root / "hg38_screen_v10_clust.regions_vs_motifs.rankings.feather",
        "hg38 cisTarget rankings database",
    )
    scores = require_file(
        resources_root / "hg38_screen_v10_clust.regions_vs_motifs.scores.feather",
        "hg38 cisTarget scores database",
    )
    motifs = require_file(
        resources_root / "motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl",
        "hg38 motif-to-TF annotation",
    )
    annotation = require_file(args.genome_annotation, "HAL-projected hg38 annotation")
    chromsizes = require_file(args.hg38_chromsizes, "HAL hg38 chromosome sizes")

    with open(cistopic, "rb") as handle:
        cistopic_object = pickle.load(handle)
    rna_object = ad.read_h5ad(rna, backed="r")
    try:
        if rna_object.raw is None:
            raise ValueError("SCENIC+ RNA H5AD has no .raw count matrix")
        rna_cells = set(map(str, rna_object.obs_names))
        rna_genes = set(map(str, rna_object.raw.var_names))
    finally:
        rna_object.file.close()
    common_cells = rna_cells.intersection(map(str, cistopic_object.cell_names))
    cell_fraction = len(common_cells) / max(1, len(cistopic_object.cell_names))
    if len(common_cells) < 100 or cell_fraction < args.min_cell_overlap:
        raise ValueError(
            f"only {len(common_cells):,}/{len(cistopic_object.cell_names):,} "
            f"cisTopic cells ({100 * cell_fraction:.2f}%) occur in the RNA H5AD"
        )
    annotation_frame = pd.read_csv(annotation, sep="\t", usecols=["Gene"])
    annotation_genes = set(annotation_frame["Gene"].dropna().astype(str))
    common_genes = rna_genes.intersection(annotation_genes)
    if len(common_genes) < args.min_gene_overlap:
        raise ValueError(
            f"only {len(common_genes):,} RNA genes match the HAL-projected gene "
            f"annotation; minimum is {args.min_gene_overlap:,}"
        )
    log(
        f"SCENIC+ input overlap: {len(common_cells):,} paired cells "
        f"({100 * cell_fraction:.2f}% of cisTopic) and {len(common_genes):,} "
        "RNA genes with projected hg38 annotations"
    )
    del cistopic_object, annotation_frame, annotation_genes, common_cells
    del common_genes, rna_cells, rna_genes
    gc.collect()

    pipeline_root = Path(args.pipeline_root).expanduser().resolve()
    snake_root = pipeline_root / "Snakemake"
    snakefile = snake_root / "workflow" / "Snakefile"
    config_path = snake_root / "config" / "config.yaml"
    if not snakefile.is_file():
        if snake_root.exists():
            raise FileExistsError(
                f"incomplete SCENIC+ Snakemake directory exists: {snake_root}"
            )
        pipeline_root.parent.mkdir(parents=True, exist_ok=True)
        run(["scenicplus", "init_snakemake", "--out_dir", str(pipeline_root)])
    require_file(snakefile, "SCENIC+ Snakefile")
    require_file(config_path, "SCENIC+ configuration template")
    with open(config_path, "r", encoding="utf-8") as handle:
        config = yaml.safe_load(handle)

    results = mkdir(pipeline_root / "results")
    temp_dir = parallel_temp_dir(args.temp_dir)
    config["input_data"].update(
        {
            "cisTopic_obj_fname": str(cistopic),
            "GEX_anndata_fname": str(rna),
            "region_set_folder": str(region_root),
            "ctx_db_fname": str(ranking),
            "dem_db_fname": str(scores),
            "path_to_motif_annotations": str(motifs),
        }
    )
    output_names = {
        "combined_GEX_ACC_mudata": "ACC_GEX.h5mu",
        "dem_result_fname": "dem_results.hdf5",
        "ctx_result_fname": "ctx_results.hdf5",
        "output_fname_dem_html": "dem_results.html",
        "output_fname_ctx_html": "ctx_results.html",
        "cistromes_direct": "cistromes_direct.h5ad",
        "cistromes_extended": "cistromes_extended.h5ad",
        "tf_names": "tf_names.txt",
        "search_space": "search_space.tsv",
        "tf_to_gene_adjacencies": "tf_to_gene_adj.tsv",
        "region_to_gene_adjacencies": "region_to_gene_adj.tsv",
        "eRegulons_direct": "eRegulon_direct.tsv",
        "eRegulons_extended": "eRegulons_extended.tsv",
        "AUCell_direct": "AUCell_direct.h5mu",
        "AUCell_extended": "AUCell_extended.h5mu",
        "scplus_mdata": "scplusmdata.h5mu",
    }
    config["output_data"].update(
        {key: str(results / name) for key, name in output_names.items()}
    )
    # These two files already exist and were built from the exact HAL/reference.
    # Snakemake therefore skips its network-dependent BioMart rule.
    config["output_data"]["genome_annotation"] = str(annotation)
    config["output_data"]["chromsizes"] = str(chromsizes)
    config["params_general"].update(
        {"temp_dir": str(temp_dir), "n_cpu": args.cpus, "seed": args.seed}
    )
    config["params_data_preparation"].update(
        {
            "bc_transform_func": '"lambda x: x"',
            "is_multiome": True,
            "species": "hsapiens",
            "search_space_upstream": "1000 150000",
            "search_space_downstream": "1000 150000",
            "search_space_extend_tss": "10 10",
        }
    )
    config["params_motif_enrichment"].update(
        {
            "species": "homo_sapiens",
            "annotation_version": "v10nr_clust",
            "annotations_to_use": "Direct_annot Orthology_annot",
        }
    )
    temporary_config = config_path.with_name(config_path.name + f".tmp.{os.getpid()}")
    with open(temporary_config, "w", encoding="utf-8") as handle:
        yaml.safe_dump(config, handle, sort_keys=False)
    os.replace(temporary_config, config_path)

    run(
        [
            "snakemake",
            "--snakefile",
            "workflow/Snakefile",
            "--configfile",
            "config/config.yaml",
            "--cores",
            str(args.cpus),
            "--rerun-incomplete",
            "--printshellcmds",
            "--latency-wait",
            "120",
        ],
        cwd=snake_root,
    )
    final = require_file(results / "scplusmdata.h5mu", "final SCENIC+ MuData")
    log(f"SCENIC+ complete: {final}")
    return 0


def positive_int(value: str) -> int:
    result = int(value)
    if result <= 0:
        raise argparse.ArgumentTypeError("must be positive")
    return result


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Tet 2025 RNA-anchored SCENIC+ worker"
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    command = subparsers.add_parser("reference")
    command.add_argument("--hal", type=Path, default=HAL_PATH)
    command.add_argument("--mao-fasta", type=Path, default=MAO_FASTA)
    command.add_argument("--mao-gtf", type=Path, default=MAO_GTF)
    command.add_argument("--mao-genome", default=MAO_GENOME)
    command.add_argument("--human-genome", default=HUMAN_GENOME)
    command.add_argument("--exclude-contig", default=PROBLEM_CONTIG)
    command.add_argument("--output-dir", type=Path, required=True)
    command.set_defaults(function=reference)

    command = subparsers.add_parser("resources")
    command.add_argument("--output-dir", type=Path, required=True)
    command.set_defaults(function=resources)

    command = subparsers.add_parser("pseudobulk-peaks")
    command.add_argument("--cell-metadata", type=Path, required=True)
    command.add_argument("--task-table", type=Path, required=True)
    command.add_argument("--mao-chromsizes", type=Path, required=True)
    command.add_argument("--output-dir", type=Path, required=True)
    command.add_argument("--temp-dir", type=Path, required=True)
    command.add_argument("--cpus", type=positive_int, required=True)
    command.add_argument("--min-cells-per-group", type=positive_int, default=50)
    command.add_argument("--q-value", type=float, default=0.05)
    command.set_defaults(function=pseudobulk_peaks)

    command = subparsers.add_parser("consensus-liftover")
    command.add_argument("--hal", type=Path, default=HAL_PATH)
    command.add_argument("--mao-genome", default=MAO_GENOME)
    command.add_argument("--human-genome", default=HUMAN_GENOME)
    command.add_argument("--exclude-contig", default=PROBLEM_CONTIG)
    command.add_argument("--narrow-peaks", type=Path, required=True)
    command.add_argument("--mao-chromsizes", type=Path, required=True)
    command.add_argument("--output-dir", type=Path, required=True)
    command.add_argument("--peak-half-width", type=positive_int, default=250)
    command.add_argument("--min-length-ratio", type=float, default=0.80)
    command.add_argument("--max-length-ratio", type=float, default=1.25)
    command.add_argument("--min-reciprocal-overlap", type=float, default=0.80)
    command.add_argument("--min-retained-peaks", type=positive_int, default=10_000)
    command.set_defaults(function=consensus_liftover)

    command = subparsers.add_parser("build-cistopic-library")
    command.add_argument("--library", required=True)
    command.add_argument("--fragments", type=Path, required=True)
    command.add_argument("--regions", type=Path, required=True)
    command.add_argument("--cell-metadata", type=Path, required=True)
    command.add_argument("--output", type=Path, required=True)
    command.add_argument("--cpus", type=positive_int, required=True)
    command.add_argument("--partitions", type=positive_int, default=10)
    command.set_defaults(function=build_cistopic_library)

    command = subparsers.add_parser("merge-cistopic")
    command.add_argument("--task-table", type=Path, required=True)
    command.add_argument("--peak-crosswalk", type=Path, required=True)
    command.add_argument("--cell-metadata", type=Path, required=True)
    command.add_argument("--output", type=Path, required=True)
    command.add_argument("--project", default="Tet2025_batch2_batch3_pools")
    command.set_defaults(function=merge_cistopic)

    command = subparsers.add_parser("topic-model")
    command.add_argument("--input", type=Path, required=True)
    command.add_argument("--output", type=Path, required=True)
    command.add_argument("--temp-dir", type=Path, required=True)
    command.add_argument("--cpus", type=positive_int, required=True)
    command.add_argument("--topics", type=positive_int, nargs="+", default=[20, 30, 40, 50])
    command.add_argument("--iterations", type=positive_int, default=150)
    command.add_argument("--seed", type=int, default=2025)
    command.add_argument("--select-topics", type=positive_int)
    command.set_defaults(function=topic_model)

    command = subparsers.add_parser("region-sets")
    command.add_argument("--input", type=Path, required=True)
    command.add_argument("--output-dir", type=Path, required=True)
    command.add_argument("--output-cistopic", type=Path, required=True)
    command.add_argument("--temp-dir", type=Path, required=True)
    command.add_argument("--cpus", type=positive_int, required=True)
    command.add_argument("--grouping", default="rna_leiden_reference")
    command.add_argument("--top-regions-per-topic", type=positive_int, default=3000)
    command.add_argument("--variable-regions", type=positive_int, default=20_000)
    command.add_argument("--imputation-chunk-size", type=positive_int, default=20_000)
    command.add_argument("--adjusted-pvalue", type=float, default=0.05)
    command.add_argument("--fold-change", type=float, default=1.5)
    command.add_argument("--min-cells-per-dar", type=positive_int, default=50)
    command.set_defaults(function=region_sets)

    command = subparsers.add_parser("run-scenicplus")
    command.add_argument("--cistopic", type=Path, required=True)
    command.add_argument("--rna-h5ad", type=Path, required=True)
    command.add_argument("--region-sets", type=Path, required=True)
    command.add_argument("--resources", type=Path, required=True)
    command.add_argument("--genome-annotation", type=Path, required=True)
    command.add_argument("--hg38-chromsizes", type=Path, required=True)
    command.add_argument("--pipeline-root", type=Path, required=True)
    command.add_argument("--temp-dir", type=Path, required=True)
    command.add_argument("--cpus", type=positive_int, required=True)
    command.add_argument("--seed", type=int, default=2025)
    command.add_argument("--min-cell-overlap", type=float, default=0.90)
    command.add_argument("--min-gene-overlap", type=positive_int, default=10_000)
    command.set_defaults(function=configure_and_run_scenicplus)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    try:
        return int(args.function(args))
    except subprocess.CalledProcessError as exc:
        print(
            f"ERROR: command failed with exit {exc.returncode}: "
            + " ".join(map(str, exc.cmd)),
            file=sys.stderr,
            flush=True,
        )
        if exc.stdout:
            print(exc.stdout, file=sys.stderr)
        if exc.stderr:
            print(exc.stderr, file=sys.stderr)
        return int(exc.returncode or 1)
    except Exception as exc:
        print(f"ERROR: {type(exc).__name__}: {exc}", file=sys.stderr, flush=True)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
