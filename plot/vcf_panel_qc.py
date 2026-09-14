#!/usr/bin/env python3
"""
Unified SNP-panel QC runner for CellBouncer downsampling outputs.

This single script replaces the separate:
  - alternative-BCF preparation shell script,
  - validation cache builder,
  - validation plot renderer,
  - VCF quality reporter, and
  - GTF annotation analyzer.

The default full run performs one indexed BCF stream per input and derives:
  - genotype/sample QC and pairwise distinguishing power,
  - SFS, IBS, PCA, and NJ-tree validation,
  - INFO score distributions,
  - species-pair summaries,
  - GTF location summaries, and
  - chrM/NUMT panel-policy checks.

The compact NPZ cache supports fast plot-only reruns without reopening BCFs.
Optional sample harmonization replaces prep_bcfs_for_alt_comparison_V1_R1.sh.
When RNA and ATAC labels are supplied, the runner also computes exact
position-level overlap, renders an overlap heatmap, and can require complete
RNA-versus-ATAC disjointness.

Typical RNA-versus-ATAC run:
  python3 vcf_panel_qc.py \
    --input rna_demux /path/rna.demux.bcf \
    --input rna_het /path/rna.het.bcf \
    --input rna_species /path/rna.species.bcf \
    --input atac_demux /path/atac.demux.bcf \
    --input atac_het /path/atac.het.bcf \
    --input atac_species /path/atac.species.bcf \
    --panel-metadata /path/panel_metadata.tsv \
    --gtf /path/reference.gtf.gz \
    --numts-bed /path/numts.bed \
    --prefix /path/qc/snp_qc \
    --require-rna-atac-disjoint \
    --threads 16

Alternative comparison with automatic reheader + sample intersection:
  python3 snp_qc_unified_V1_R1.py \
    --input input /path/input.bcf \
    --input species /path/species.bcf \
    --alt-input species species_alt /path/nathan_species.bcf \
    --rename-sample C6007B=C6007 \
    --rename-sample Chinobo=Chinobo-mCherry \
    --harmonize-samples \
    --panel-metadata /path/panel_metadata.tsv \
    --prefix /path/qc/snp_qc \
    --threads 16

Revision history is at the end of the file.
"""

import argparse
import collections
from array import array
import datetime as _datetime
import gzip
import json
import math
import multiprocessing
import os
import re
import shutil
import subprocess
import sys
import tempfile
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

try:
    import numpy as np
except ImportError:
    np = None

try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.gridspec import GridSpec
    from matplotlib.patches import Patch
except ImportError:
    matplotlib = None
    plt = None
    GridSpec = None
    Patch = None

VERSION = '1.2.1'
CACHE_SCHEMA_VERSION = 3

DEFAULT_GTF = (
    '/mnt/beegfs/genomes_annotations/ancestral_genomes/litterbox/'
    'human_chimp_bonobo/human_chimp_bonobo.gtf.gz'
)
DEFAULT_NUMTS_BED = (
    '/mnt/beegfs/genomes_annotations/ancestral_genomes/litterbox/'
    'human_chimp_bonobo/numts.bed'
)

NUMERIC_SCORE_FIELDS = [
    'DEMUX_SCORE', 'HET_SCORE', 'HET_FREQ', 'ANNOT_SCORE',
    'COV_SCORE', 'HET_SEL_SCORE',
    'ATAC_REGULATORY',
    'MT_COV_SCORE', 'MT_LIBRARY_COUNT', 'MT_PAIR_USE_COUNT',
    'MT_RATIO_USE_COUNT', 'MT_AMBIENT_USE_COUNT',
    'MT_DEPTH_FILTER', 'MT_HOMOPLASMY_AF_FILTER',
]

LABEL_COLORS = {
    'rna_demux': '#F18F01',
    'rna_het': '#2E86AB',
    'atac_demux': '#2A9D8F',
    'atac_het': '#17BEBB',
    'rna_species': '#8E44AD',
    'atac_species': '#5B8C5A',
    'atac_numt': '#6C757D',
    'atac_mito': '#D1495B',
    'species': '#7B2D8E',
    'mito': '#D1495B',
    'mt': '#D1495B',
    'input': '#888888',
    'all': '#888888',
}
ALT_LABEL_COLORS = {
    'rna_demux': '#D62828',
    'rna_het': '#06AED5',
    'atac_demux': '#5A189A',
    'atac_het': '#FFB703',
    'species': '#FFB400',
    'mito': '#5F0F40',
    'mt': '#5F0F40',
    'input': '#3D3D3D',
    'all': '#3D3D3D',
}
DEFAULT_COLOR = '#888888'
C_ACCENT = '#C73E1D'
C_GOOD = '#2A9D8F'
C_HOM_REF = '#3288bd'
C_HET_GT = '#99d594'
C_HOM_ALT = '#d53e4f'
C_EXONIC = '#2A9D8F'
C_INTRONIC = '#E9C46A'
C_INTERGENIC = '#264653'
C_5PRIME = '#E76F51'
C_3PRIME = '#457B9D'

SPECIES_COLORS = {
    'human': '#E63946', 'h': '#E63946',
    'chimp': '#457B9D', 'chimpanzee': '#457B9D', 'c': '#457B9D',
    'bonobo': '#2A9D8F', 'b': '#2A9D8F',
    'orangutan': '#F18F01', 'orang': '#F18F01', 'o': '#F18F01',
    'chinobo': '#7B2D8E', 'hybrid': '#7B2D8E', 'hy': '#7B2D8E',
}
SPECIES_PRETTY = {
    'h': 'Human', 'human': 'Human',
    'c': 'Chimp', 'chimp': 'Chimp', 'chimpanzee': 'Chimp',
    'b': 'Bonobo', 'bonobo': 'Bonobo',
    'o': 'Orangutan', 'orangutan': 'Orangutan', 'orang': 'Orangutan',
    'hy': 'Chinobo (F1)', 'chinobo': 'Chinobo (F1)', 'hybrid': 'Chinobo (F1)',
}
LABEL_PRETTY = {
    'rna_demux': 'RNA demux',
    'rna_het': 'RNA het',
    'atac_demux': 'ATAC demux',
    'atac_het': 'ATAC het',
    'rna_species': 'RNA species',
    'atac_species': 'ATAC species',
    'atac_numt': 'ATAC NUMT diagnostic',
    'atac_mito': 'ATAC mitochondrial',
    'species': 'Species',
    'mito': 'Mitochondrial',
    'mt': 'Mitochondrial',
    'input': 'Input',
    'all': 'Input',
}
EXPECTED_GROUPS = {
    'human': 0, 'h': 0,
    'chimp': 1, 'chimpanzee': 1, 'c': 1,
    'bonobo': 1, 'b': 1,
    'chinobo': 1, 'hybrid': 1, 'hy': 1,
    'orangutan': 2, 'orang': 2, 'o': 2,
}

GENE_REL_BINS = np.linspace(0.0, 1.0, 51) if np is not None else None
INTERGENIC_LOG_BINS = np.linspace(0.0, 8.5, 86) if np is not None else None


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


def require_numpy():
    if np is None:
        raise RuntimeError('numpy is required. Load the project Python/miniforge module first.')


def require_matplotlib():
    if plt is None:
        raise RuntimeError('matplotlib is required for plotting. Load the project Python/miniforge module first.')


def require_executable(name):
    path = shutil.which(name)
    if not path:
        raise RuntimeError(f'Required executable not found on PATH: {name}')
    return path


def run_checked(cmd, stdout=None, text=True, description=None):
    if description:
        print(description, flush=True)
    result = subprocess.run(cmd, stdout=stdout, stderr=subprocess.PIPE,
                            text=text, check=False)
    if result.returncode != 0:
        stderr = result.stderr.strip() if isinstance(result.stderr, str) else ''
        raise RuntimeError(
            f"Command failed ({result.returncode}): {' '.join(map(str, cmd))}"
            + (f"\n{stderr}" if stderr else '')
        )
    return result


def parse_flexible_input(values, expected, flag_name):
    if len(values) == expected:
        return tuple(values)
    if len(values) == 1:
        parts = values[0].split(None, expected - 1)
        if len(parts) == expected:
            return tuple(parts)
    raise ValueError(
        f"{flag_name} requires {expected} values; received: {values!r}"
    )


def normalize_input_specs(raw_inputs, raw_alts):
    specs = []
    labels = set()
    for raw in raw_inputs or []:
        label, path = parse_flexible_input(raw, 2, '--input')
        if label in labels:
            raise ValueError(f'Duplicate dataset label: {label}')
        labels.add(label)
        specs.append({'label': label, 'path': os.path.abspath(path), 'alt_of': None})
    primary_labels = set(labels)
    for raw in raw_alts or []:
        orig, label, path = parse_flexible_input(raw, 3, '--alt-input')
        if orig not in primary_labels:
            raise ValueError(
                f"--alt-input original label '{orig}' is not a primary --input label"
            )
        if label in labels:
            raise ValueError(f'Duplicate dataset label: {label}')
        labels.add(label)
        specs.append({'label': label, 'path': os.path.abspath(path), 'alt_of': orig})
    primary = [s for s in specs if s['alt_of'] is None]
    alts_by_primary = defaultdict(list)
    for s in specs:
        if s['alt_of'] is not None:
            alts_by_primary[s['alt_of']].append(s)
    ordered = []
    for s in primary:
        ordered.append(s)
        ordered.extend(alts_by_primary.get(s['label'], []))
    return ordered


def parse_rename_map(values):
    mapping = {}
    for item in values or []:
        if '=' not in item:
            raise ValueError(f"--rename-sample requires OLD=NEW, received: {item}")
        old, new = item.split('=', 1)
        old = old.strip()
        new = new.strip()
        if not old or not new:
            raise ValueError(f"Invalid --rename-sample value: {item}")
        if old in mapping and mapping[old] != new:
            raise ValueError(f"Conflicting rename for {old}: {mapping[old]} vs {new}")
        mapping[old] = new
    return mapping


def get_samples(vcf_file):
    result = run_checked(['bcftools', 'query', '-l', vcf_file], stdout=subprocess.PIPE)
    return [x for x in result.stdout.splitlines() if x]


def get_header_text(vcf_file):
    result = run_checked(['bcftools', 'view', '-h', vcf_file], stdout=subprocess.PIPE)
    return result.stdout


def parse_info_header(vcf_file):
    info = {}
    pattern = re.compile(r'^##INFO=<ID=([^,>]+),Number=([^,>]+),Type=([^,>]+)')
    for line in get_header_text(vcf_file).splitlines():
        m = pattern.match(line)
        if m:
            info[m.group(1)] = {'number': m.group(2), 'type': m.group(3)}
    return info


def get_available_info_fields(vcf_file):
    return set(parse_info_header(vcf_file))


def get_chromosomes(vcf_file):
    result = subprocess.run(['bcftools', 'index', '-s', vcf_file],
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                            text=True, check=False)
    if result.returncode == 0 and result.stdout.strip():
        return [line.split('\t', 1)[0] for line in result.stdout.splitlines() if line]
    chroms = []
    for line in get_header_text(vcf_file).splitlines():
        if line.startswith('##contig=<ID='):
            value = line[len('##contig=<ID='):].split(',', 1)[0].split('>', 1)[0]
            if value:
                chroms.append(value)
    return chroms


def load_panel_metadata(panel_file):
    species_map = {}
    with open(panel_file, 'rt') as handle:
        header = handle.readline()
        if not header:
            raise RuntimeError(f'Empty panel metadata: {panel_file}')
        for line_no, line in enumerate(handle, start=2):
            if not line.strip() or line.startswith('#'):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 2:
                eprint(f'WARNING: panel metadata line {line_no} has fewer than 2 columns')
                continue
            sample = parts[0].strip()
            species = parts[1].strip().lower()
            if sample:
                species_map[sample] = species
    if not species_map:
        raise RuntimeError(f'No sample/species mappings loaded from {panel_file}')
    return species_map


def check_panel_coverage(samples, sample_species, label):
    matched = [s for s in samples if s in sample_species]
    unmatched = [s for s in samples if s not in sample_species]
    if not matched:
        raise RuntimeError(
            f"[{label}] none of the {len(samples)} BCF samples are present in panel metadata"
        )
    if unmatched:
        eprint(
            f"WARNING [{label}]: {len(matched)}/{len(samples)} samples matched panel metadata; "
            f"unmatched: {', '.join(unmatched[:12])}"
            + (' ...' if len(unmatched) > 12 else '')
        )
    else:
        print(f'  [{label}] panel metadata coverage: {len(samples)}/{len(samples)}', flush=True)


def ensure_bcf_index(vcf_file, threads, create_missing=False):
    result = subprocess.run(['bcftools', 'index', '-s', vcf_file],
                            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
                            check=False)
    if result.returncode == 0:
        return
    if not create_missing:
        raise RuntimeError(
            f'BCF/VCF is not indexed or index is unreadable: {vcf_file}. '
            f'Create an index or pass --index-missing.'
        )
    print(f'  Indexing missing input: {vcf_file}', flush=True)
    run_checked(['bcftools', 'index', '--csi', '--force', '--threads', str(threads), vcf_file])


def harmonize_inputs(specs, rename_map, out_dir, threads, min_intersection, force):
    """Reheader alternative inputs, intersect samples across all BCFs, and subset.

    Original inputs are never modified. Renames are applied only to --alt-input files.
    """
    require_executable('bcftools')
    out_dir = os.path.abspath(out_dir)
    os.makedirs(out_dir, exist_ok=True)
    print('=' * 70)
    print('Preparing sample-harmonized comparison BCFs')
    print('=' * 70)

    working = []
    for spec in specs:
        src = spec['path']
        if not os.path.isfile(src):
            raise RuntimeError(f'Missing input BCF: {src}')
        work_path = src
        if spec['alt_of'] is not None and rename_map:
            present = set(get_samples(src))
            active = {old: new for old, new in rename_map.items() if old in present}
            if active:
                map_path = os.path.join(out_dir, f"{spec['label']}.rename_map.tsv")
                renamed = os.path.join(out_dir, f"{spec['label']}.renamed.bcf")
                if (os.path.exists(renamed) or os.path.exists(renamed + '.csi')) and not force:
                    raise RuntimeError(f'Prepared file exists; pass --force to replace: {renamed}')
                with open(map_path, 'wt') as handle:
                    for old, new in active.items():
                        handle.write(f'{old}\t{new}\n')
                run_checked([
                    'bcftools', 'reheader', '--samples', map_path,
                    '--threads', str(threads), '-o', renamed, src,
                ], description=f"  Reheadering {spec['label']}")
                run_checked(['bcftools', 'index', '--csi', '--force',
                             '--threads', str(threads), renamed])
                work_path = renamed
            else:
                print(f"  {spec['label']}: none of the requested rename keys were present", flush=True)
        item = dict(spec)
        item['_working_path'] = work_path
        item['_samples'] = get_samples(work_path)
        working.append(item)

    if not working:
        raise RuntimeError('No inputs supplied for harmonization')
    common = set(working[0]['_samples'])
    for item in working[1:]:
        common &= set(item['_samples'])
    canonical_order = [s for s in working[0]['_samples'] if s in common]
    if len(canonical_order) < min_intersection:
        raise RuntimeError(
            f'Sample intersection has only {len(canonical_order)} samples; '
            f'minimum is {min_intersection}'
        )
    sample_file = os.path.join(out_dir, 'intersect_samples.txt')
    with open(sample_file, 'wt') as handle:
        handle.write('\n'.join(canonical_order) + '\n')
    print(f'  Shared samples: {len(canonical_order)}', flush=True)

    manifest_rows = []
    result_specs = []
    for item in working:
        label = item['label']
        out_bcf = os.path.join(out_dir, f'{label}.intersect.bcf')
        if (os.path.exists(out_bcf) or os.path.exists(out_bcf + '.csi')) and not force:
            raise RuntimeError(f'Prepared file exists; pass --force to replace: {out_bcf}')
        for stale in (out_bcf, out_bcf + '.csi'):
            if os.path.exists(stale):
                os.remove(stale)
        run_checked([
            'bcftools', 'view', '--samples-file', sample_file,
            '--threads', str(threads), '-O', 'b', '-o', out_bcf,
            item['_working_path'],
        ], description=f'  Subsetting {label}')
        run_checked(['bcftools', 'index', '--csi', '--force',
                     '--threads', str(threads), out_bcf])
        observed = get_samples(out_bcf)
        if set(observed) != set(canonical_order):
            raise RuntimeError(f'Post-subset sample mismatch for {out_bcf}')
        out_spec = {'label': label, 'path': out_bcf, 'alt_of': item['alt_of']}
        result_specs.append(out_spec)
        manifest_rows.append((label, item['path'], out_bcf, len(observed)))

    manifest = os.path.join(out_dir, 'prepared_inputs.tsv')
    with open(manifest, 'wt') as handle:
        handle.write('label\toriginal_path\tprepared_path\tn_samples\n')
        for row in manifest_rows:
            handle.write('\t'.join(map(str, row)) + '\n')
    print(f'Prepared inputs manifest: {manifest}')
    return result_specs


def load_numt_bed(path):
    """Load BED intervals and convert 0-based half-open BED to 1-based inclusive."""
    if not path:
        return {}
    if not os.path.isfile(path):
        raise RuntimeError(f'NUMT BED not found: {path}')
    intervals = defaultdict(list)
    opener = gzip.open if path.endswith('.gz') else open
    mode = 'rt'
    with opener(path, mode) as handle:
        for line_no, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith('#'):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 3:
                raise RuntimeError(f'Invalid BED line {line_no} in {path}')
            chrom = parts[0]
            try:
                start0 = int(parts[1])
                end0 = int(parts[2])
            except ValueError:
                raise RuntimeError(f'Invalid BED coordinates at line {line_no} in {path}')
            if end0 <= start0:
                continue
            intervals[chrom].append((start0 + 1, end0))
    merged = {}
    for chrom, vals in intervals.items():
        vals.sort()
        out = []
        for start, end in vals:
            if not out or start > out[-1][1] + 1:
                out.append([start, end])
            else:
                out[-1][1] = max(out[-1][1], end)
        merged[chrom] = {
            'starts': np.array([x[0] for x in out], dtype=np.int64),
            'ends': np.array([x[1] for x in out], dtype=np.int64),
        }
    return merged


class GeneIndex:
    """Per-chromosome gene/exon arrays for batched 1-based position annotation."""

    def __init__(self):
        self.genes = defaultdict(list)
        self.exons = defaultdict(list)
        self.chrom_payloads = {}

    def add_gene(self, chrom, start, end, strand, gene_type):
        self.genes[chrom].append((start, end, strand, gene_type))

    def add_exon(self, chrom, start, end):
        self.exons[chrom].append((start, end))

    @staticmethod
    def _prefix_max(ends):
        if len(ends) == 0:
            return np.zeros(0, dtype=np.int64), np.zeros(0, dtype=np.int32)
        values = np.empty(len(ends), dtype=np.int64)
        indices = np.empty(len(ends), dtype=np.int32)
        best_end = -1
        best_idx = -1
        for i, end in enumerate(ends):
            if int(end) > best_end:
                best_end = int(end)
                best_idx = i
            values[i] = best_end
            indices[i] = best_idx
        return values, indices

    def finalize(self):
        chroms = set(self.genes) | set(self.exons)
        for chrom in chroms:
            genes = sorted(self.genes.get(chrom, []), key=lambda x: (x[0], x[1]))
            exons = sorted(self.exons.get(chrom, []), key=lambda x: (x[0], x[1]))
            g_starts = np.array([x[0] for x in genes], dtype=np.int64)
            g_ends = np.array([x[1] for x in genes], dtype=np.int64)
            g_strands = np.array([1 if x[2] == '+' else -1 for x in genes], dtype=np.int8)
            g_types = [x[3] for x in genes]
            g_prefix_end, g_prefix_idx = self._prefix_max(g_ends)
            e_starts = np.array([x[0] for x in exons], dtype=np.int64)
            e_ends = np.array([x[1] for x in exons], dtype=np.int64)
            e_prefix_end, _ = self._prefix_max(e_ends)
            self.chrom_payloads[chrom] = {
                'gene_starts': g_starts,
                'gene_ends': g_ends,
                'gene_strands': g_strands,
                'gene_types': g_types,
                'gene_prefix_max_end': g_prefix_end,
                'gene_prefix_max_idx': g_prefix_idx,
                'exon_starts': e_starts,
                'exon_ends': e_ends,
                'exon_prefix_max_end': e_prefix_end,
            }

    def payload(self, chrom):
        return self.chrom_payloads.get(chrom)


def load_gtf(gtf_file):
    print(f'Loading GTF: {gtf_file}', flush=True)
    if not os.path.isfile(gtf_file):
        raise RuntimeError(f'GTF not found: {gtf_file}')
    index = GeneIndex()
    opener = gzip.open if gtf_file.endswith('.gz') else open
    gene_count = 0
    exon_count = 0
    with opener(gtf_file, 'rt') as handle:
        for line in handle:
            if not line.strip() or line.startswith('#'):
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 9:
                continue
            chrom, _, feature, start, end, _, strand, _, attributes = fields
            try:
                start_i = int(start)
                end_i = int(end)
            except ValueError:
                continue
            attrs = {}
            for attr in attributes.split(';'):
                attr = attr.strip()
                if not attr or ' ' not in attr:
                    continue
                key, value = attr.split(' ', 1)
                attrs[key] = value.strip().strip('"')
            gene_type = attrs.get('gene_type', attrs.get('gene_biotype', 'unknown'))
            if feature == 'gene':
                index.add_gene(chrom, start_i, end_i, strand, gene_type)
                gene_count += 1
            elif feature == 'exon':
                index.add_exon(chrom, start_i, end_i)
                exon_count += 1
    index.finalize()
    print(f'  Loaded {gene_count:,} genes and {exon_count:,} exons '
          f'across {len(index.chrom_payloads)} contigs', flush=True)
    return index


def _empty_annotation_summary():
    return {
        'exonic': 0,
        'intronic': 0,
        'intergenic': 0,
        'gene_rel_hist': np.zeros(len(GENE_REL_BINS) - 1, dtype=np.int64),
        'intergenic_log_hist': np.zeros(len(INTERGENIC_LOG_BINS) - 1, dtype=np.int64),
        'gene_5prime': 0,
        'gene_middle': 0,
        'gene_3prime': 0,
        'intergenic_distance_count': 0,
        'gene_type_counts': Counter(),
    }


def _merge_annotation_summary(target, source):
    for key in ('exonic', 'intronic', 'intergenic', 'gene_5prime',
                'gene_middle', 'gene_3prime', 'intergenic_distance_count'):
        target[key] += int(source.get(key, 0))
    target['gene_rel_hist'] += source['gene_rel_hist']
    target['intergenic_log_hist'] += source['intergenic_log_hist']
    target['gene_type_counts'].update(source.get('gene_type_counts', {}))


def annotate_positions_batch(positions, payload):
    """Annotate one chromosome batch without the nested-gene early-break bug."""
    out = _empty_annotation_summary()
    n = len(positions)
    if n == 0:
        return out
    if payload is None or len(payload['gene_starts']) == 0:
        out['intergenic'] = n
        return out

    pos_arr = np.asarray(positions, dtype=np.int64)
    starts = payload['gene_starts']
    ends = payload['gene_ends']
    strands = payload['gene_strands']
    types = payload['gene_types']
    prefix_end = payload['gene_prefix_max_end']
    prefix_idx = payload['gene_prefix_max_idx']
    e_starts = payload['exon_starts']
    e_ends = payload['exon_ends']
    e_prefix_end = payload['exon_prefix_max_end']

    rel_values = []
    intergenic_log = []
    type_counts = Counter()

    for pos in pos_arr:
        idx = int(np.searchsorted(starts, pos, side='right'))
        chosen = -1
        j = idx - 1
        while j >= 0 and prefix_end[j] >= pos:
            if ends[j] >= pos:
                chosen = j
                break
            j -= 1

        if chosen >= 0:
            is_exon = False
            if len(e_starts):
                ei = int(np.searchsorted(e_starts, pos, side='right')) - 1
                while ei >= 0 and e_prefix_end[ei] >= pos:
                    if e_ends[ei] >= pos:
                        is_exon = True
                        break
                    ei -= 1
            if is_exon:
                out['exonic'] += 1
            else:
                out['intronic'] += 1
            gene_len = int(ends[chosen] - starts[chosen])
            if gene_len > 0:
                rel = float(pos - starts[chosen]) / float(gene_len)
                if strands[chosen] < 0:
                    rel = 1.0 - rel
                rel = min(1.0, max(0.0, rel))
                rel_values.append(rel)
                if rel < 0.25:
                    out['gene_5prime'] += 1
                elif rel > 0.75:
                    out['gene_3prime'] += 1
                else:
                    out['gene_middle'] += 1
            if chosen < len(types):
                type_counts[types[chosen]] += 1
            continue

        out['intergenic'] += 1
        best_dist = None
        best_gene = -1
        if idx > 0:
            up_idx = int(prefix_idx[idx - 1])
            up_dist = int(pos - ends[up_idx])
            if up_dist > 0:
                best_dist = up_dist
                best_gene = up_idx
        if idx < len(starts):
            down_dist = int(starts[idx] - pos)
            if down_dist > 0 and (best_dist is None or down_dist < best_dist):
                best_dist = down_dist
                best_gene = idx
        if best_dist is not None:
            out['intergenic_distance_count'] += 1
            intergenic_log.append(math.log10(max(1, best_dist)))
            if best_dist <= 10000 and 0 <= best_gene < len(types):
                type_counts[types[best_gene]] += 1

    if rel_values:
        out['gene_rel_hist'] += np.histogram(rel_values, bins=GENE_REL_BINS)[0]
    if intergenic_log:
        vals = np.clip(np.asarray(intergenic_log), INTERGENIC_LOG_BINS[0],
                       np.nextafter(INTERGENIC_LOG_BINS[-1], -np.inf))
        out['intergenic_log_hist'] += np.histogram(vals, bins=INTERGENIC_LOG_BINS)[0]
    out['gene_type_counts'].update(type_counts)
    return out


def count_interval_overlaps(positions, interval_payload):
    if interval_payload is None or len(positions) == 0:
        return 0
    starts = interval_payload['starts']
    ends = interval_payload['ends']
    if len(starts) == 0:
        return 0
    pos = np.asarray(positions, dtype=np.int64)
    idx = np.searchsorted(starts, pos, side='right') - 1
    valid = idx >= 0
    result = np.zeros(len(pos), dtype=bool)
    result[valid] = ends[idx[valid]] >= pos[valid]
    return int(np.sum(result))


def process_region_worker(payload):
    (vcf_file, region, n_samples, score_fields, aux_numeric_fields,
     string_fields, subsample_stride, score_sample_stride,
     annotation_payload, numt_payload) = payload

    gt_map = {
        '.': -1, './.': -1, '.|.': -1,
        '0/0': 0, '0|0': 0, '0': 0,
        '0/1': 1, '1/0': 1, '0|1': 1, '1|0': 1,
        '1/1': 2, '1|1': 2, '1': 2,
    }
    all_numeric = list(score_fields) + list(aux_numeric_fields)
    all_info = all_numeric + list(string_fields)
    fmt = '%CHROM\t%POS'
    for field in all_info:
        fmt += f'\t%INFO/{field}'
    fmt += '\t[%GT\t]\n'
    cmd = ['bcftools', 'query', '-f', fmt, '-r', region, vcf_file]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                            text=True, bufsize=1024 * 1024)

    sample_ac = np.zeros(n_samples, dtype=np.int64)
    sample_called = np.zeros(n_samples, dtype=np.int64)
    ibs_match = np.zeros((n_samples, n_samples), dtype=np.int64)
    ibs_total = np.zeros((n_samples, n_samples), dtype=np.int64)
    diff_count = np.zeros((n_samples, n_samples), dtype=np.int64)
    pair_matrix = np.zeros((n_samples, n_samples), dtype=np.int64)
    het_counts = np.zeros(n_samples, dtype=np.int64)
    hom_ref_counts = np.zeros(n_samples, dtype=np.int64)
    hom_alt_counts = np.zeros(n_samples, dtype=np.int64)
    missing_counts = np.zeros(n_samples, dtype=np.int64)
    sfs_counts = np.zeros(2 * n_samples + 1, dtype=np.int64)
    score_samples = {f: [] for f in score_fields}
    score_stats = {f: {'count': 0, 'sum': 0.0, 'min': None, 'max': None}
                   for f in score_fields}
    aux_stats = {f: {'count': 0, 'sum': 0.0} for f in aux_numeric_fields}
    string_counts = {f: Counter() for f in string_fields}
    string_uniques = {f: set() for f in string_fields}
    gt_subsample = []
    annotation = _empty_annotation_summary()
    numt_overlap = 0
    total_sites = 0
    locus_positions = array('I')
    batch_gt = []
    batch_pos = []
    batch_size = 10000
    n_info = len(all_info)

    def flush_batch():
        nonlocal batch_gt, batch_pos, numt_overlap
        if not batch_gt:
            return
        arr = np.asarray(batch_gt, dtype=np.int8)
        called = arr >= 0
        sample_ac[:] += np.sum(np.where(called, arr, 0), axis=0)
        sample_called[:] += 2 * np.sum(called, axis=0)
        het_counts[:] += np.sum(arr == 1, axis=0)
        hom_ref_counts[:] += np.sum(arr == 0, axis=0)
        hom_alt_counts[:] += np.sum(arr == 2, axis=0)
        missing_counts[:] += np.sum(arr == -1, axis=0)

        site_ac = np.sum(np.where(called, arr, 0), axis=1)
        site_called = np.sum(called, axis=1)
        full_mask = site_called == arr.shape[1]
        if np.any(full_mask):
            full_ac = np.clip(site_ac[full_mask], 0, 2 * n_samples)
            np.add.at(sfs_counts, full_ac, 1)

        a = arr[:, :, np.newaxis]
        b = arr[:, np.newaxis, :]
        both_called = (a >= 0) & (b >= 0)
        non_overlap = ((a == 0) & (b == 2)) | ((a == 2) & (b == 0))
        ibs_total[:] += np.sum(both_called, axis=0)
        ibs_match[:] += np.sum(both_called & ~non_overlap, axis=0)
        diff_count[:] += np.sum(non_overlap, axis=0)
        pair_matrix[:] += np.sum(both_called & (a != b), axis=0)

        if annotation_payload is not None:
            _merge_annotation_summary(
                annotation,
                annotate_positions_batch(batch_pos, annotation_payload),
            )
        numt_overlap += count_interval_overlaps(batch_pos, numt_payload)
        batch_gt = []
        batch_pos = []

    for line in proc.stdout:
        fields = line.rstrip('\n').rstrip('\t').split('\t')
        expected_min = 2 + n_info + n_samples
        if len(fields) < expected_min:
            continue
        try:
            pos = int(fields[1])
        except ValueError:
            continue
        info_values = fields[2:2 + n_info]
        for idx, field in enumerate(all_numeric):
            raw = info_values[idx]
            try:
                value = float(raw) if raw not in ('', '.') and ',' not in raw else math.nan
            except ValueError:
                value = math.nan
            if math.isfinite(value):
                if field in score_stats:
                    stat = score_stats[field]
                    stat['count'] += 1
                    stat['sum'] += value
                    stat['min'] = value if stat['min'] is None else min(stat['min'], value)
                    stat['max'] = value if stat['max'] is None else max(stat['max'], value)
                    if total_sites % score_sample_stride == 0:
                        score_samples[field].append(value)
                else:
                    aux_stats[field]['count'] += 1
                    aux_stats[field]['sum'] += value
        for sidx, field in enumerate(string_fields):
            raw = info_values[len(all_numeric) + sidx]
            if raw not in ('', '.'):
                if field == 'BIN_ID':
                    string_uniques[field].add(raw)
                else:
                    string_counts[field][raw] += 1

        gt_start = 2 + n_info
        gts = [gt_map.get(fields[gt_start + k].strip(), -1)
               for k in range(n_samples)]
        batch_gt.append(gts)
        batch_pos.append(pos)
        locus_positions.append(pos)
        if total_sites % subsample_stride == 0:
            gt_subsample.append(gts)
        total_sites += 1
        if len(batch_gt) >= batch_size:
            flush_batch()

    flush_batch()
    stderr = proc.stderr.read()
    return_code = proc.wait()
    if return_code != 0:
        raise RuntimeError(
            f"bcftools query failed for {vcf_file} region {region}: "
            f"{stderr.strip()}"
        )

    return {
        'region': region,
        'total_sites': total_sites,
        'sample_ac': sample_ac,
        'sample_called': sample_called,
        'ibs_match': ibs_match,
        'ibs_total': ibs_total,
        'diff_count': diff_count,
        'pair_matrix': pair_matrix,
        'het_counts': het_counts,
        'hom_ref_counts': hom_ref_counts,
        'hom_alt_counts': hom_alt_counts,
        'missing_counts': missing_counts,
        'sfs_counts': sfs_counts,
        'score_samples': {k: np.asarray(v, dtype=np.float32)
                          for k, v in score_samples.items()},
        'score_stats': score_stats,
        'aux_stats': aux_stats,
        'string_counts': {k: dict(v) for k, v in string_counts.items()},
        'string_unique_values': {k: sorted(v) for k, v in string_uniques.items()},
        'gt_subsample': (np.asarray(gt_subsample, dtype=np.int8)
                         if gt_subsample else np.zeros((0, n_samples), dtype=np.int8)),
        'annotation': annotation,
        'numt_overlap': numt_overlap,
        'locus_positions': np.unique(np.frombuffer(locus_positions, dtype=np.uint32)),
    }


def analyze_vcf(vcf_file, label, alt_of, sample_species, gene_index,
                numt_intervals, subsample_stride, score_sample_stride,
                n_threads):
    print('\n' + '=' * 70)
    print(f'Analyzing {label}: {vcf_file}')
    print('=' * 70, flush=True)
    samples = get_samples(vcf_file)
    if not samples:
        raise RuntimeError(f'No samples found in {vcf_file}')
    check_panel_coverage(samples, sample_species, label)
    chroms = get_chromosomes(vcf_file)
    if not chroms:
        raise RuntimeError(f'No indexed contigs found in {vcf_file}')
    info_header = parse_info_header(vcf_file)
    score_fields = [f for f in NUMERIC_SCORE_FIELDS if f in info_header]
    aux_numeric_fields = sorted(f for f in info_header if f.startswith('PAIR_DISCRIM_'))
    string_fields = [f for f in ('PAIR_ASSIGNED', 'BIN_ID') if f in info_header]
    print(f'  Samples: {len(samples)}')
    print(f'  Contigs: {len(chroms)}')
    print(f"  Numeric score fields: {', '.join(score_fields) if score_fields else 'none'}")
    print(f'  Threads: {n_threads}', flush=True)

    args = []
    for chrom in chroms:
        args.append((
            vcf_file, chrom, len(samples), score_fields, aux_numeric_fields,
            string_fields, subsample_stride, score_sample_stride,
            gene_index.payload(chrom) if gene_index is not None else None,
            numt_intervals.get(chrom),
        ))

    total_sites = 0
    chrom_counts = {}
    loci_by_chrom = {}
    n = len(samples)
    combined = {
        'sample_ac': np.zeros(n, dtype=np.int64),
        'sample_called': np.zeros(n, dtype=np.int64),
        'ibs_match': np.zeros((n, n), dtype=np.int64),
        'ibs_total': np.zeros((n, n), dtype=np.int64),
        'diff_count': np.zeros((n, n), dtype=np.int64),
        'pair_matrix': np.zeros((n, n), dtype=np.int64),
        'het_counts': np.zeros(n, dtype=np.int64),
        'hom_ref_counts': np.zeros(n, dtype=np.int64),
        'hom_alt_counts': np.zeros(n, dtype=np.int64),
        'missing_counts': np.zeros(n, dtype=np.int64),
        'sfs_counts': np.zeros(2 * n + 1, dtype=np.int64),
    }
    score_chunks = {f: [] for f in score_fields}
    score_stats = {f: {'count': 0, 'sum': 0.0, 'min': None, 'max': None}
                   for f in score_fields}
    aux_stats = {f: {'count': 0, 'sum': 0.0} for f in aux_numeric_fields}
    string_counts = {f: Counter() for f in string_fields}
    string_unique_values = {f: set() for f in string_fields}
    subsample_chunks = []
    annotation = _empty_annotation_summary()
    numt_overlap = 0

    with ProcessPoolExecutor(max_workers=n_threads) as executor:
        futures = {executor.submit(process_region_worker, arg): arg[1] for arg in args}
        completed = 0
        try:
            for future in as_completed(futures):
                chrom = futures[future]
                result = future.result()
                chrom_counts[chrom] = int(result['total_sites'])
                loci_by_chrom[chrom] = result['locus_positions']
                total_sites += int(result['total_sites'])
                for key in combined:
                    combined[key] += result[key]
                for field in score_fields:
                    if len(result['score_samples'][field]):
                        score_chunks[field].append(result['score_samples'][field])
                    src = result['score_stats'][field]
                    dst = score_stats[field]
                    dst['count'] += src['count']
                    dst['sum'] += src['sum']
                    if src['min'] is not None:
                        dst['min'] = src['min'] if dst['min'] is None else min(dst['min'], src['min'])
                    if src['max'] is not None:
                        dst['max'] = src['max'] if dst['max'] is None else max(dst['max'], src['max'])
                for field in aux_numeric_fields:
                    aux_stats[field]['count'] += result['aux_stats'][field]['count']
                    aux_stats[field]['sum'] += result['aux_stats'][field]['sum']
                for field in string_fields:
                    string_counts[field].update(result['string_counts'][field])
                    string_unique_values[field].update(result['string_unique_values'][field])
                if result['gt_subsample'].shape[0]:
                    subsample_chunks.append(result['gt_subsample'])
                _merge_annotation_summary(annotation, result['annotation'])
                numt_overlap += int(result['numt_overlap'])
                completed += 1
                if completed % 100 == 0 or completed == len(chroms):
                    print(f'    {completed}/{len(chroms)} contigs; '
                          f'{total_sites:,} variants', flush=True)
        except Exception:
            for future in futures:
                future.cancel()
            raise

    scores = {
        f: (np.concatenate(score_chunks[f]) if score_chunks[f]
            else np.zeros(0, dtype=np.float32))
        for f in score_fields
    }
    gt_subsample = (
        np.vstack(subsample_chunks)
        if subsample_chunks else np.zeros((0, n), dtype=np.int8)
    )
    pair_assigned_counts = dict(string_counts.get('PAIR_ASSIGNED', {}))
    bin_count = len(string_unique_values.get('BIN_ID', set()))
    pair_discrim_means = {
        field: (stat['sum'] / stat['count'] if stat['count'] else math.nan)
        for field, stat in aux_stats.items()
    }
    print(f'  Total variants: {total_sites:,}')
    print(f'  GT subsample: {gt_subsample.shape[0]:,} x {n}')
    if numt_intervals:
        print(f'  NUMT-overlap records: {numt_overlap:,}')

    data = {
        'label': label,
        'alt_of': alt_of,
        '_vcf_file': vcf_file,
        'samples': samples,
        'total_sites': total_sites,
        'chrom_counts': chrom_counts,
        'scores': scores,
        'score_stats': score_stats,
        'pair_assigned_counts': pair_assigned_counts,
        'pair_discrim_means': pair_discrim_means,
        'bin_count': bin_count,
        'gt_subsample': gt_subsample,
        'annotation': annotation,
        'numt_overlap': numt_overlap,
        '_loci_by_chrom': loci_by_chrom,
    }
    data.update(combined)
    return data


def save_cache(cache_file, all_data, sample_species, run_meta):
    cache_file = os.path.abspath(cache_file)
    os.makedirs(os.path.dirname(cache_file) or '.', exist_ok=True)
    arrays = {
        'species_map_keys': np.asarray(list(sample_species.keys()), dtype='U'),
        'species_map_values': np.asarray(list(sample_species.values()), dtype='U'),
    }
    meta = {
        'schema_version': CACHE_SCHEMA_VERSION,
        'script_version': VERSION,
        'created_at': _datetime.datetime.now().isoformat(timespec='seconds'),
        'run_meta': run_meta,
        'datasets': [],
    }
    for i, data in enumerate(all_data):
        score_field_order = sorted(data['scores'])
        gene_type_counts = dict(data['annotation'].get('gene_type_counts', {}))
        dmeta = {
            'label': data['label'],
            'alt_of': data.get('alt_of'),
            'path': data.get('_vcf_file', ''),
            'total_sites': int(data['total_sites']),
            'chrom_counts': {k: int(v) for k, v in data.get('chrom_counts', {}).items()},
            'score_fields': score_field_order,
            'score_stats': data.get('score_stats', {}),
            'pair_assigned_counts': data.get('pair_assigned_counts', {}),
            'pair_discrim_means': data.get('pair_discrim_means', {}),
            'bin_count': int(data.get('bin_count', 0)),
            'numt_overlap': int(data.get('numt_overlap', 0)),
            'annotation': {
                key: int(data['annotation'].get(key, 0))
                for key in ('exonic', 'intronic', 'intergenic', 'gene_5prime',
                            'gene_middle', 'gene_3prime', 'intergenic_distance_count')
            },
            'gene_type_counts': gene_type_counts,
        }
        meta['datasets'].append(dmeta)
        arrays[f'samples_{i}'] = np.asarray(data['samples'], dtype='U')
        for key in ('sample_ac', 'sample_called', 'ibs_match', 'ibs_total',
                    'diff_count', 'pair_matrix', 'het_counts', 'hom_ref_counts',
                    'hom_alt_counts', 'missing_counts', 'sfs_counts',
                    'gt_subsample'):
            arrays[f'{key}_{i}'] = data[key]
        arrays[f'annotation_gene_rel_hist_{i}'] = data['annotation']['gene_rel_hist']
        arrays[f'annotation_intergenic_log_hist_{i}'] = data['annotation']['intergenic_log_hist']
        for j, field in enumerate(score_field_order):
            arrays[f'score_{i}_{j}'] = data['scores'][field]
    arrays['metadata_json'] = np.asarray(json.dumps(meta, sort_keys=True), dtype='U')

    tmp_file = cache_file + '.tmp.npz'
    if os.path.exists(tmp_file):
        os.remove(tmp_file)
    np.savez_compressed(tmp_file, **arrays)
    os.replace(tmp_file, cache_file)
    print(f'Cache written: {cache_file} ({os.path.getsize(cache_file) / 1e6:.1f} MB)')


def load_cache(cache_file):
    print(f'Loading cache: {cache_file}', flush=True)
    z = np.load(cache_file, allow_pickle=False)
    meta = json.loads(str(z['metadata_json']))
    if int(meta.get('schema_version', -1)) != CACHE_SCHEMA_VERSION:
        raise RuntimeError(
            f"Unsupported cache schema {meta.get('schema_version')}; "
            f'expected {CACHE_SCHEMA_VERSION}'
        )
    species_map = dict(zip([str(x) for x in z['species_map_keys']],
                           [str(x) for x in z['species_map_values']]))
    all_data = []
    for i, dmeta in enumerate(meta['datasets']):
        annotation = dict(dmeta.get('annotation', {}))
        annotation['gene_rel_hist'] = z[f'annotation_gene_rel_hist_{i}']
        annotation['intergenic_log_hist'] = z[f'annotation_intergenic_log_hist_{i}']
        annotation['gene_type_counts'] = Counter(dmeta.get('gene_type_counts', {}))
        scores = {}
        for j, field in enumerate(dmeta.get('score_fields', [])):
            scores[field] = z[f'score_{i}_{j}']
        data = {
            'label': dmeta['label'],
            'alt_of': dmeta.get('alt_of'),
            '_vcf_file': dmeta.get('path', ''),
            'samples': [str(x) for x in z[f'samples_{i}']],
            'total_sites': int(dmeta['total_sites']),
            'chrom_counts': dmeta.get('chrom_counts', {}),
            'scores': scores,
            'score_stats': dmeta.get('score_stats', {}),
            'pair_assigned_counts': dmeta.get('pair_assigned_counts', {}),
            'pair_discrim_means': dmeta.get('pair_discrim_means', {}),
            'bin_count': int(dmeta.get('bin_count', 0)),
            'numt_overlap': int(dmeta.get('numt_overlap', 0)),
            'annotation': annotation,
        }
        for key in ('sample_ac', 'sample_called', 'ibs_match', 'ibs_total',
                    'diff_count', 'pair_matrix', 'het_counts', 'hom_ref_counts',
                    'hom_alt_counts', 'missing_counts', 'sfs_counts',
                    'gt_subsample'):
            data[key] = z[f'{key}_{i}']
        all_data.append(data)
    return all_data, species_map, meta.get('run_meta', {})


def label_modality(label):
    label = str(label).lower()
    if label in ('atac_numt', 'atac_mito'):
        return None
    if label.startswith('rna_'):
        return 'rna'
    if label.startswith('atac_'):
        return 'atac'
    return None


def subset_dataset_samples(data, excluded_samples):
    """Return a shallow dataset copy with sample-indexed arrays subsetted."""
    excluded = set(excluded_samples)
    keep = [i for i, sample in enumerate(data['samples']) if sample not in excluded]
    if len(keep) == len(data['samples']):
        return data
    out = dict(data)
    out['samples'] = [data['samples'][i] for i in keep]
    for key in ('sample_ac', 'sample_called', 'het_counts', 'hom_ref_counts',
                'hom_alt_counts', 'missing_counts'):
        out[key] = np.asarray(data[key])[keep]
    for key in ('ibs_match', 'ibs_total', 'diff_count', 'pair_matrix'):
        matrix = np.asarray(data[key])
        out[key] = matrix[np.ix_(keep, keep)]
    out['gt_subsample'] = np.asarray(data['gt_subsample'])[:, keep]
    return out


def _locus_count(loci_by_chrom):
    return sum(len(values) for values in loci_by_chrom.values())


def _intersection_count(left, right):
    total = 0
    for chrom in set(left) & set(right):
        total += int(np.intersect1d(
            left[chrom], right[chrom], assume_unique=True,
        ).size)
    return total


def _union_loci(datasets):
    chunks = defaultdict(list)
    for data in datasets:
        for chrom, values in data.get('_loci_by_chrom', {}).items():
            if len(values):
                chunks[chrom].append(values)
    union = {}
    for chrom, arrays in chunks.items():
        union[chrom] = (
            arrays[0] if len(arrays) == 1
            else np.unique(np.concatenate(arrays))
        )
    return union


def compute_locus_overlap(all_data):
    """Compute exact position-level overlap and RNA/ATAC union disjointness."""
    rows = []
    for i, left in enumerate(all_data):
        left_loci = left.get('_loci_by_chrom', {})
        n_left = _locus_count(left_loci)
        for j in range(i, len(all_data)):
            right = all_data[j]
            right_loci = right.get('_loci_by_chrom', {})
            n_right = _locus_count(right_loci)
            shared = n_left if i == j else _intersection_count(left_loci, right_loci)
            union = n_left + n_right - shared
            rows.append({
                'label_a': left['label'],
                'label_b': right['label'],
                'modality_a': label_modality(left['label']) or 'other',
                'modality_b': label_modality(right['label']) or 'other',
                'n_a': n_left,
                'n_b': n_right,
                'shared_loci': shared,
                'a_unique': n_left - shared,
                'b_unique': n_right - shared,
                'overlap_fraction_smaller': (
                    shared / min(n_left, n_right) if min(n_left, n_right) else 0.0
                ),
                'jaccard': shared / union if union else 0.0,
            })

    rna_data = [d for d in all_data if label_modality(d['label']) == 'rna']
    atac_data = [d for d in all_data if label_modality(d['label']) == 'atac']
    rna_union = _union_loci(rna_data)
    atac_union = _union_loci(atac_data)
    n_rna = _locus_count(rna_union)
    n_atac = _locus_count(atac_union)
    shared = _intersection_count(rna_union, atac_union)
    union_total = n_rna + n_atac - shared
    summary = {
        'rna_labels': [d['label'] for d in rna_data],
        'atac_labels': [d['label'] for d in atac_data],
        'rna_union_loci': n_rna,
        'atac_union_loci': n_atac,
        'shared_loci': shared,
        'rna_unique_loci': n_rna - shared,
        'atac_unique_loci': n_atac - shared,
        'overlap_fraction_smaller': (
            shared / min(n_rna, n_atac) if min(n_rna, n_atac) else 0.0
        ),
        'jaccard': shared / union_total if union_total else 0.0,
        'status': 'PASS' if rna_data and atac_data and shared == 0 else (
            'FAIL' if rna_data and atac_data else 'NOT_TESTED'
        ),
    }
    return rows, summary


def write_overlap_tsv(rows, output_file):
    columns = [
        'label_a', 'label_b', 'modality_a', 'modality_b', 'n_a', 'n_b',
        'shared_loci', 'a_unique', 'b_unique', 'overlap_fraction_smaller',
        'jaccard',
    ]
    with open(output_file, 'wt') as handle:
        handle.write('\t'.join(columns) + '\n')
        for row in rows:
            handle.write('\t'.join(
                f"{row[column]:.10g}" if isinstance(row[column], float)
                else str(row[column])
                for column in columns
            ) + '\n')


def write_disjointness_tsv(summary, output_file):
    columns = [
        'status', 'rna_labels', 'atac_labels', 'rna_union_loci',
        'atac_union_loci', 'shared_loci', 'rna_unique_loci',
        'atac_unique_loci', 'overlap_fraction_smaller', 'jaccard',
    ]
    row = dict(summary)
    row['rna_labels'] = ','.join(row.get('rna_labels', []))
    row['atac_labels'] = ','.join(row.get('atac_labels', []))
    with open(output_file, 'wt') as handle:
        handle.write('\t'.join(columns) + '\n')
        handle.write('\t'.join(
            f"{row[column]:.10g}" if isinstance(row.get(column), float)
            else str(row.get(column, ''))
            for column in columns
        ) + '\n')


def create_overlap_figure(rows, summary, output_file):
    require_matplotlib()
    labels = []
    for row in rows:
        for key in ('label_a', 'label_b'):
            if row[key] not in labels:
                labels.append(row[key])
    if not labels:
        return False
    index = {label: i for i, label in enumerate(labels)}
    matrix = np.zeros((len(labels), len(labels)), dtype=float)
    counts = np.zeros((len(labels), len(labels)), dtype=np.int64)
    for row in rows:
        i = index[row['label_a']]
        j = index[row['label_b']]
        matrix[i, j] = matrix[j, i] = row['overlap_fraction_smaller']
        counts[i, j] = counts[j, i] = row['shared_loci']

    width = max(8.5, 1.05 * len(labels) + 3.0)
    fig, ax = plt.subplots(figsize=(width, width * 0.82))
    image = ax.imshow(matrix * 100.0, vmin=0, vmax=100, cmap='viridis')
    for i in range(len(labels)):
        for j in range(len(labels)):
            color = 'white' if matrix[i, j] < 0.55 else 'black'
            ax.text(j, i, f"{100 * matrix[i, j]:.2f}%\n{counts[i, j]:,}",
                    ha='center', va='center', fontsize=8, color=color)
    pretty = [get_label_pretty(label) for label in labels]
    ax.set_xticks(range(len(labels)), pretty, rotation=35, ha='right')
    ax.set_yticks(range(len(labels)), pretty)
    ax.set_title('SNP-locus overlap across RNA and ATAC panels',
                 fontsize=15, fontweight='bold', pad=14)
    subtitle = (
        f"RNA union: {summary.get('rna_union_loci', 0):,} loci   |   "
        f"ATAC union: {summary.get('atac_union_loci', 0):,} loci   |   "
        f"shared: {summary.get('shared_loci', 0):,}   |   "
        f"disjointness: {summary.get('status', 'NOT_TESTED')}"
    )
    ax.text(0.5, 1.01, subtitle, transform=ax.transAxes,
            ha='center', va='bottom', fontsize=10)
    colorbar = fig.colorbar(image, ax=ax, fraction=0.046, pad=0.04)
    colorbar.set_label('Shared loci as % of smaller panel')
    fig.tight_layout()
    fig.savefig(output_file, dpi=180, bbox_inches='tight')
    plt.close(fig)
    return True


def get_label_pretty(label, alt_of=None):
    """Pretty-print a dataset label. If alt_of is given, the label is
    treated as a paired alternative and rendered as 'orig_pretty (alt)'
    when the user-supplied alt_label is just '<orig>_alt' or similar."""
    pretty = LABEL_PRETTY.get(label, label)
    return pretty


def get_color(label, alt_of=None):
    """Return color for a dataset label. If alt_of is given, the label is
    treated as an alternative of alt_of and gets a fully distinct color
    from ALT_LABEL_COLORS so paired comparisons read clearly."""
    if alt_of is not None:
        return ALT_LABEL_COLORS.get(alt_of, DEFAULT_COLOR)
    return LABEL_COLORS.get(label, DEFAULT_COLOR)


def get_species_color(species):
    if species is None:
        return DEFAULT_COLOR
    return SPECIES_COLORS.get(species.lower(), DEFAULT_COLOR)


def get_species_pretty(species):
    if species is None:
        return 'unknown'
    return SPECIES_PRETTY.get(species.lower(), species)


def get_clade(species):
    """Return clade index (0/1/2) or None if unknown."""
    if species is None:
        return None
    return EXPECTED_GROUPS.get(species.lower())


def compute_genotype_distance(gt_matrix):
    """Pairwise allelic distance from genotype matrix (sites x samples)."""
    n = gt_matrix.shape[1]
    dist = np.zeros((n, n), dtype=np.float64)
    called = gt_matrix >= 0
    gt_f = gt_matrix.astype(np.float32)
    gt_f[~called] = 0.0
    for i in range(n):
        for j in range(i + 1, n):
            joint = called[:, i] & called[:, j]
            njoint = np.sum(joint)
            if njoint == 0:
                d = 0.0
            else:
                diff = np.abs(gt_f[joint, i] - gt_f[joint, j])
                d = np.mean(diff) / 2.0
            dist[i, j] = d
            dist[j, i] = d
    return dist


def neighbor_join(dist, labels):
    """Saitou-Nei neighbor joining."""
    n = dist.shape[0]
    D = dist.astype(np.float64).copy()
    active = list(range(n))
    edges = []
    next_node = n

    while len(active) > 2:
        m = len(active)
        sub = np.array([[D[a, b] for b in active] for a in active])
        r = np.sum(sub, axis=1) / (m - 2)
        Q = sub - r[:, None] - r[None, :]
        np.fill_diagonal(Q, np.inf)
        ij = np.unravel_index(np.argmin(Q), Q.shape)
        i, j = ij
        a, b = active[i], active[j]
        d_ab = D[a, b]
        d_au = 0.5 * d_ab + 0.5 * (r[i] - r[j])
        d_bu = d_ab - d_au
        u = next_node
        next_node += 1
        edges.append((u, a, max(d_au, 0.0)))
        edges.append((u, b, max(d_bu, 0.0)))
        new_size = D.shape[0] + 1
        new_D = np.zeros((new_size, new_size), dtype=np.float64)
        new_D[:-1, :-1] = D
        for k in range(D.shape[0]):
            if k == a or k == b:
                continue
            d_ku = 0.5 * (D[a, k] + D[b, k] - d_ab)
            new_D[u, k] = d_ku
            new_D[k, u] = d_ku
        D = new_D
        active = [x for x in active if x != a and x != b] + [u]

    a, b = active
    d_ab = D[a, b]
    edges.append((a, b, max(d_ab, 0.0)))
    return edges, a


def topology_score(gt_matrix, samples, sample_species):
    """Fraction of within-clade pairs closer than nearest cross-clade pair."""
    if gt_matrix.shape[0] < 100:
        return float('nan')
    dist = compute_genotype_distance(gt_matrix)
    n = len(samples)
    within_correct = 0
    within_total = 0
    for i in range(n):
        clade_i = get_clade(sample_species.get(samples[i]))
        if clade_i is None:
            continue
        min_cross = np.inf
        for j in range(n):
            if i == j:
                continue
            clade_j = get_clade(sample_species.get(samples[j]))
            if clade_j is None:
                continue
            if clade_i != clade_j:
                if dist[i, j] < min_cross:
                    min_cross = dist[i, j]
        for j in range(n):
            if i == j:
                continue
            clade_j = get_clade(sample_species.get(samples[j]))
            if clade_j is None:
                continue
            if clade_i == clade_j:
                within_total += 1
                if dist[i, j] < min_cross:
                    within_correct += 1
    if within_total == 0:
        return float('nan')
    return within_correct / within_total


def midpoint_root(edges, n_leaves):
    """Find the midpoint of the longest path in the tree and return a node
    closest to it as the root. The midpoint root puts the deepest split
    at the top of the layout, fixing the issue where a single clade gets
    split across the y-axis just because the NJ termination node happened
    to sit inside it.

    Algorithm:
      1. Pick any leaf, find the leaf farthest from it (call it A).
      2. From A, find the leaf farthest from it (call it B). Path A-B is
         the diameter of the tree.
      3. Walk along path A->B, accumulating distance, until we cross half
         the path's length. Return the node we're at (or just past).
    """
    adj = defaultdict(list)
    for u, v, l in edges:
        adj[u].append((v, l))
        adj[v].append((u, l))

    if not adj:
        return 0

    def farthest_from(start):
        """BFS-like search returning (farthest_node, distance, parent_dict)."""
        dist_to = {start: 0.0}
        parent_of = {start: None}
        stack = [start]
        visited = {start}
        max_dist = 0.0
        max_node = start
        # DFS over tree (no cycles in trees, so this is fine)
        while stack:
            cur = stack.pop()
            for v, l in adj[cur]:
                if v not in visited:
                    visited.add(v)
                    dist_to[v] = dist_to[cur] + l
                    parent_of[v] = cur
                    if dist_to[v] > max_dist:
                        max_dist = dist_to[v]
                        max_node = v
                    stack.append(v)
        return max_node, max_dist, parent_of

    # Step 1: from any leaf (use 0), find farthest
    a, _, _ = farthest_from(0)
    # Step 2: from a, find farthest = b. parent_of gives path a -> b
    b, total_len, parent_of = farthest_from(a)

    # Step 3: walk from b back toward a, accumulating distance
    target = total_len / 2.0
    path = []
    cur = b
    while cur is not None:
        path.append(cur)
        cur = parent_of[cur]
    # path is b -> ... -> a; we want to walk from a toward b
    path.reverse()  # now a -> ... -> b
    # Compute cumulative distance along path
    cum = 0.0
    for i in range(len(path) - 1):
        u = path[i]
        v = path[i + 1]
        # Find edge length between u and v in adj
        edge_len = 0.0
        for nbr, l in adj[u]:
            if nbr == v:
                edge_len = l
                break
        if cum + edge_len >= target:
            # The midpoint falls on edge u-v. Return the closer endpoint.
            half_remaining = target - cum
            return v if half_remaining > edge_len / 2 else u
        cum += edge_len

    return path[len(path) // 2]


def _draw_phylogeny_panel(ax, data, sample_species, shared_max_x=None):
    """Draw a single NJ phylogeny panel."""
    gt = data['gt_subsample']
    samples = data['samples']
    n = len(samples)
    if gt.shape[0] < 50:
        ax.text(0.5, 0.5, 'Insufficient sites for tree',
                ha='center', va='center', transform=ax.transAxes,
                fontsize=12, color='gray')
        ax.axis('off')
        return None

    dist = compute_genotype_distance(gt)
    edges, _nj_root = neighbor_join(dist, samples)
    root = midpoint_root(edges, n)
    coords, leaf_order, parent = build_rectangular_layout(
        edges, root, n, sample_species, samples)

    # Light banded background grouping by clade
    clade_colors = {0: '#fef3f3', 1: '#f3f7fa', 2: '#fef8f0'}
    cur_clade = None
    band_start = None
    for i, leaf in enumerate(leaf_order):
        sp = sample_species.get(samples[leaf])
        clade = get_clade(sp)
        if clade != cur_clade:
            if cur_clade is not None and band_start is not None:
                ax.axhspan(band_start - 0.5, i - 0.5,
                           color=clade_colors.get(cur_clade, '#ffffff'),
                           alpha=0.5, zorder=0)
            cur_clade = clade
            band_start = i
    if band_start is not None and cur_clade is not None:
        ax.axhspan(band_start - 0.5, len(leaf_order) - 0.5,
                   color=clade_colors.get(cur_clade, '#ffffff'),
                   alpha=0.5, zorder=0)

    # Draw edges as right-angle elbows
    for u, v, l in edges:
        x0, y0 = coords[u]
        x1, y1 = coords[v]
        ax.plot([x0, x0], [y0, y1], color='#444444', lw=0.9, alpha=0.8, zorder=1)
        ax.plot([x0, x1], [y1, y1], color='#444444', lw=0.9, alpha=0.8, zorder=1)

    panel_max_x = max(c[0] for c in coords.values())
    max_x = shared_max_x if shared_max_x is not None else panel_max_x
    # Place leaf markers at panel-local max so trees show their actual extent
    leaf_x = panel_max_x * 1.02

    for leaf in leaf_order:
        x, y = coords[leaf]
        sp = sample_species.get(samples[leaf])
        c = get_species_color(sp)
        ax.plot([x, leaf_x], [y, y], color='#cccccc', lw=0.5,
                linestyle=':', alpha=0.6, zorder=1)
        ax.scatter([leaf_x], [y], s=70, color=c, edgecolor='black',
                   lw=0.6, zorder=3)
        ax.text(leaf_x + max_x * 0.012, y, samples[leaf],
                fontsize=8, ha='left', va='center',
                color=c, fontweight='bold')

    score = topology_score(gt, samples, sample_species)

    ax.set_xlabel('Genotype distance from root', fontsize=10)
    ax.set_yticks([])
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    # Use shared_max_x for x-axis bound if provided, else panel-local
    label_room = max_x * 0.30
    ax.set_xlim(-max_x * 0.02, max(leaf_x, max_x) + label_room)
    ax.set_ylim(-1, n)
    ax.invert_yaxis()

    return {'panel_max_x': panel_max_x, 'topology_score': score, 'n_sites': gt.shape[0]}


def _draw_pca_panel(ax, data, sample_species):
    """Draw a single PCA panel. Returns dict with metrics, or None if skipped."""
    gt = data['gt_subsample']
    samples = data['samples']
    if gt.shape[0] < 50:
        ax.text(0.5, 0.5, 'Insufficient sites for PCA',
                ha='center', va='center', transform=ax.transAxes,
                fontsize=12, color='gray')
        return None

    X = gt.T.astype(np.float32).copy()
    for j in range(X.shape[1]):
        col = X[:, j]
        mask = col >= 0
        if np.any(mask):
            col[~mask] = np.mean(col[mask])
        else:
            col[:] = 0
        X[:, j] = col
    X = X - np.mean(X, axis=0, keepdims=True)
    try:
        U, S, Vt = np.linalg.svd(X, full_matrices=False)
        pc_scores = U * S
        var_explained = (S ** 2) / np.sum(S ** 2)
    except np.linalg.LinAlgError:
        ax.text(0.5, 0.5, 'SVD failed', ha='center', va='center',
                transform=ax.transAxes, fontsize=12, color='gray')
        return None

    pc1 = pc_scores[:, 0] if pc_scores.shape[1] > 0 else np.zeros(len(samples))
    pc2 = pc_scores[:, 1] if pc_scores.shape[1] > 1 else np.zeros(len(samples))

    species_groups = defaultdict(list)
    for i, s in enumerate(samples):
        sp = sample_species.get(s)
        if sp is not None:
            species_groups[sp].append(pc1[i])
    if len(species_groups) >= 2:
        grand_mean = np.mean(pc1)
        ss_between = sum(len(v) * (np.mean(v) - grand_mean) ** 2
                         for v in species_groups.values())
        ss_total = np.sum((pc1 - grand_mean) ** 2)
        r2_pc1 = ss_between / ss_total if ss_total > 0 else 0.0
    else:
        r2_pc1 = float('nan')

    used_species = sorted(set(sample_species.get(s) for s in samples
                              if sample_species.get(s) is not None))

    # Plot per-species with convex-hull-like background fill if >=3 individuals
    for sp in used_species:
        xs = []
        ys = []
        for s, x, y in zip(samples, pc1, pc2):
            if sample_species.get(s) == sp:
                xs.append(x)
                ys.append(y)
        if len(xs) >= 3:
            # Background ellipse approximation: 2-sigma cluster halo
            mx, my = np.mean(xs), np.mean(ys)
            sx, sy = np.std(xs), np.std(ys)
            if sx > 0 and sy > 0:
                circle = plt.matplotlib.patches.Ellipse(
                    (mx, my), 4 * sx + 0.01 * abs(mx),
                    4 * sy + 0.01 * abs(my),
                    facecolor=get_species_color(sp), alpha=0.12,
                    edgecolor=get_species_color(sp),
                    linewidth=1.0, linestyle='--', zorder=1)
                ax.add_patch(circle)
        ax.scatter(xs, ys, s=70, color=get_species_color(sp),
                   edgecolor='black', lw=0.5, alpha=0.9,
                   label=get_species_pretty(sp), zorder=3)

    # Per-panel legend removed - shared figure legend at top
    v1 = var_explained[0] * 100 if len(var_explained) > 0 else 0
    v2 = var_explained[1] * 100 if len(var_explained) > 1 else 0
    ax.set_xlabel(f'PC1 ({v1:.1f}%)', fontsize=11)
    ax.set_ylabel(f'PC2 ({v2:.1f}%)', fontsize=11)
    ax.grid(True, alpha=0.3)
    ax.axhline(0, color='#cccccc', lw=0.5, zorder=0)
    ax.axvline(0, color='#cccccc', lw=0.5, zorder=0)

    return {'pc1_var': v1, 'pc2_var': v2, 'r2_pc1': r2_pc1, 'n_sites': gt.shape[0]}


def create_structure_figure(all_data, sample_species, output_file):
    """Combined phylogeny + PCA figure: one row per input set, two columns."""
    print(f"Creating: {output_file}", flush=True)
    n_inputs = len(all_data)
    # Estimate row height from sample count: ~0.18 inch per leaf, min 4.5
    n_samples_max = max(len(d['samples']) for d in all_data)
    row_height = max(4.5, n_samples_max * 0.20)
    fig = plt.figure(figsize=(20, row_height * n_inputs))
    gs = GridSpec(n_inputs, 2, figure=fig, hspace=0.4, wspace=0.18,
                  width_ratios=[1.4, 1.0])

    # First pass: compute max tree x across all panels for shared scale,
    # and PC1/PC2 axis ranges for cross-row PCA comparability
    shared_max_x = 0.0
    pca_ranges = []  # (pc1_range, pc2_range) per dataset
    for data in all_data:
        gt = data['gt_subsample']
        samples = data['samples']
        if gt.shape[0] < 50:
            continue
        dist = compute_genotype_distance(gt)
        edges, _nj_root = neighbor_join(dist, samples)
        root = midpoint_root(edges, len(samples))
        coords, _, _ = build_rectangular_layout(
            edges, root, len(samples), sample_species, samples)
        m = max(c[0] for c in coords.values())
        if m > shared_max_x:
            shared_max_x = m

    # Single shared species legend at top to avoid five identical legends
    used_species_all = set()
    for data in all_data:
        for s in data['samples']:
            sp = sample_species.get(s)
            if sp is not None:
                used_species_all.add(sp)
    clade_order_keys = ['h', 'human', 'c', 'chimp', 'chimpanzee',
                        'b', 'bonobo', 'hy', 'chinobo', 'hybrid',
                        'o', 'orangutan', 'orang']
    order_idx = {sp: i for i, sp in enumerate(clade_order_keys)}
    used_species_sorted = sorted(used_species_all,
                                 key=lambda s: order_idx.get(s.lower(), 99))
    species_legend_handles = [
        plt.Line2D([0], [0], marker='o', color='w',
                   markerfacecolor=get_species_color(sp),
                   markeredgecolor='black', markeredgewidth=0.6,
                   markersize=10, label=get_species_pretty(sp))
        for sp in used_species_sorted
    ]

    for idx, data in enumerate(all_data):
        ax_tree = fig.add_subplot(gs[idx, 0])
        ax_pca = fig.add_subplot(gs[idx, 1])

        tree_info = _draw_phylogeny_panel(ax_tree, data, sample_species,
                                          shared_max_x=shared_max_x)
        pca_info = _draw_pca_panel(ax_pca, data, sample_species)

        # Row title spans both panels via a centered text (left of tree)
        row_label = get_label_pretty(data['label'])
        alt_suffix = ''
        if data.get('alt_of') is not None:
            alt_suffix = f" [alt of {get_label_pretty(data['alt_of'])}]"
        if tree_info is not None:
            score_str = (f"{tree_info['topology_score']:.0%}"
                         if not np.isnan(tree_info['topology_score']) else 'n/a')
            tree_title = (f"Phylogeny - {row_label}{alt_suffix}\n"
                          f"({tree_info['n_sites']:,} sites, topology accuracy: {score_str})")
        else:
            tree_title = f"Phylogeny - {row_label}{alt_suffix}"
        ax_tree.set_title(tree_title, fontsize=12, fontweight='bold')

        if pca_info is not None:
            r2_str = (f"{pca_info['r2_pc1']:.0%}"
                      if not np.isnan(pca_info['r2_pc1']) else 'n/a')
            pca_title = (f"PCA - {row_label}{alt_suffix}\n"
                         f"({pca_info['n_sites']:,} sites, "
                         f"PC1 species R^2: {r2_str})")
        else:
            pca_title = f"PCA - {row_label}{alt_suffix}"
        ax_pca.set_title(pca_title, fontsize=12, fontweight='bold')

    # Figure-level species legend at top (replaces 4 per-panel legends)
    fig.legend(handles=species_legend_handles, loc='upper center',
               bbox_to_anchor=(0.5, 0.98), ncol=len(species_legend_handles),
               fontsize=11, framealpha=0.95, title='Species',
               title_fontsize=11)

    fig.suptitle('Population Structure: NJ Phylogeny + PCA',
                 fontsize=15, fontweight='bold', y=1.0)
    fig.text(0.5, -0.005,
             'Each row represents one SNP selection set. Left: NJ phylogeny from genotype '
             'distance (shared x-axis scale across rows for cross-set comparison). '
             'Right: PCA of mean-imputed genotype matrix; ellipses are 2-sigma cluster halos. '
             'Topology accuracy = fraction of within-clade pairs closer than nearest cross-clade pair. '
             'Higher PC1 R^2 = species more cleanly separated on the leading axis.',
             ha='center', va='top', fontsize=9, style='italic',
             color='#444444', wrap=True)
    plt.savefig(output_file, dpi=150, bbox_inches='tight', facecolor='white')
    plt.close()


def create_sfs_figure(all_data, output_file):
    """Two-panel SFS figure:
       Left: per-MAC grouped bar chart (log scale, side-by-side bars per MAC)
       Right: stacked MAF bin proportions per dataset (most actionable view)
    """
    print(f"Creating: {output_file}", flush=True)
    fig = plt.figure(figsize=(17, 6.5))
    gs = GridSpec(1, 2, figure=fig, wspace=0.25, width_ratios=[1.5, 1.0])

    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])

    # ---- Panel 1: per-MAC fractions, grouped bars ----
    n_inputs = len(all_data)
    if n_inputs == 0:
        return

    # Compute folded SFS per dataset
    folded_per_dataset = []
    for data in all_data:
        sfs = data['sfs_counts']
        n_samples = len(data['samples'])
        sfs_max = 2 * n_samples
        folded = np.zeros(sfs_max // 2 + 1, dtype=np.int64)
        for i in range(sfs_max + 1):
            j = min(i, sfs_max - i)
            folded[j] += sfs[i]
        folded_per_dataset.append((data, folded, sfs_max))

    # Common x-axis: max MAC across datasets
    max_mac = max(len(f[1]) - 1 for f in folded_per_dataset)
    bar_width = 0.8 / n_inputs
    x_centers = np.arange(1, max_mac + 1)

    for offset, (data, folded, sfs_max) in enumerate(folded_per_dataset):
        y = folded[1:].astype(np.float64)
        if len(y) < max_mac:
            y = np.concatenate([y, np.zeros(max_mac - len(y))])
        y_total = y.sum()
        y_frac = y / y_total if y_total > 0 else y
        positions = x_centers + (offset - n_inputs / 2 + 0.5) * bar_width
        color = get_color(data['label'], alt_of=data.get('alt_of'))
        legend_text = get_label_pretty(data['label'])
        if data.get('alt_of') is not None:
            legend_text = f"{legend_text} (alt of {get_label_pretty(data['alt_of'])})"
        ax1.bar(positions, y_frac, width=bar_width, color=color,
                alpha=0.85, edgecolor='white', lw=0.3,
                label=f"{legend_text} "
                      f"({int(y_total):,} fully-called of {data['total_sites']:,})")

    # Get final y-axis limits to position labels properly
    ax1.set_xlabel('Minor Allele Count (across all individuals)', fontsize=11)
    ax1.set_ylabel('Fraction of Variant Sites', fontsize=11)
    ax1.set_title('Folded Site Frequency Spectrum (per-MAC fractions)',
                  fontsize=12, fontweight='bold')
    ax1.set_xlim(0.5, max_mac + 0.5)
    ax1.grid(True, alpha=0.3, axis='y')

    # Shaded MAF region bands
    sfs_max_ref = folded_per_dataset[0][2]
    ymax = ax1.get_ylim()[1]
    # Reserve top 12% for legend
    ax1.set_ylim(0, ymax * 1.15)
    ymax_new = ymax * 1.15
    label_y = ymax * 1.07  # below the legend
    for xlo_frac, xhi_frac, fill, lbl in [
        (0.0, 0.05, '#fff5f0', 'rare\n<5%'),
        (0.05, 0.20, '#fee8c8', 'low\n5-20%'),
        (0.20, 0.40, '#fdbb84', 'mid\n20-40%'),
        (0.40, 0.50, '#e34a33', 'common\n>40%'),
    ]:
        xlo = xlo_frac * sfs_max_ref
        xhi = min(xhi_frac * sfs_max_ref, max_mac + 0.5)
        ax1.axvspan(xlo, xhi, color=fill, alpha=0.25, zorder=0)
        if xhi - xlo > 0.5:  # only label visible regions
            ax1.text((xlo + xhi) / 2, label_y, lbl, ha='center', va='center',
                     fontsize=9, color='#555555', fontweight='bold', zorder=2)

    ax1.legend(fontsize=9, loc='upper right', framealpha=0.95,
               bbox_to_anchor=(0.99, 0.92))

    # ---- Panel 2: stacked MAF bin summary (this is the actionable view) ----
    bin_names = ['rare\n(<5%)', 'low\n(5-20%)', 'mid\n(20-40%)', 'common\n(>40%)']
    bin_colors = ['#fdd0a2', '#fdae6b', '#e6550d', '#a63603']
    labels = []
    rare_pcts = []
    low_pcts = []
    mid_pcts = []
    common_pcts = []
    for data, folded, sfs_max in folded_per_dataset:
        var_sites = folded[1:].sum()
        if var_sites == 0:
            continue
        maf = np.arange(1, len(folded)) / sfs_max
        rare = folded[1:][maf < 0.05].sum() / var_sites
        low = folded[1:][(maf >= 0.05) & (maf < 0.20)].sum() / var_sites
        mid = folded[1:][(maf >= 0.20) & (maf < 0.40)].sum() / var_sites
        common = folded[1:][maf >= 0.40].sum() / var_sites
        label_text = get_label_pretty(data['label'])
        if data.get('alt_of') is not None:
            label_text = f"{label_text}\n(alt)"
        labels.append(label_text)
        rare_pcts.append(rare * 100)
        low_pcts.append(low * 100)
        mid_pcts.append(mid * 100)
        common_pcts.append(common * 100)

    x = np.arange(len(labels))
    bottoms = np.zeros(len(labels))
    bin_data = [rare_pcts, low_pcts, mid_pcts, common_pcts]
    bin_lbls = ['rare (<5%)', 'low (5-20%)', 'mid (20-40%)', 'common (>40%)']
    for vals, name, color in zip(bin_data, bin_lbls, bin_colors):
        vals_arr = np.array(vals)
        ax2.bar(x, vals_arr, bottom=bottoms, width=0.65, color=color,
                edgecolor='white', lw=1.0, label=name)
        # Annotate cells with %
        for i, v in enumerate(vals_arr):
            if v >= 4:  # only show if visible
                ax2.text(x[i], bottoms[i] + v / 2, f"{v:.1f}%",
                         ha='center', va='center', fontsize=9,
                         color='white' if v > 12 else 'black',
                         fontweight='bold')
        bottoms += vals_arr

    ax2.set_xticks(x)
    ax2.set_xticklabels(labels, fontsize=10, rotation=0)
    ax2.set_ylabel('Percent of Variant Sites', fontsize=11)
    ax2.set_title('MAF Composition by Selection Set',
                  fontsize=12, fontweight='bold')
    ax2.set_ylim(0, 100)
    ax2.legend(fontsize=9, loc='center left',
               bbox_to_anchor=(1.02, 0.5), framealpha=0.95,
               title='MAF bin', title_fontsize=10)
    ax2.grid(True, alpha=0.3, axis='y')

    fig.text(0.5, -0.04,
             'Left: bars are the fraction of sites at each minor allele count, side-by-side per dataset. '
             'Right: same data binned into MAF categories. Demux sets should be enriched for mid/common '
             'sites (informative for demultiplexing). Het sets should be dominated by low MAF (within-species '
             'heterozygosity). Species sets should be bimodal: rare singletons (within-species) plus '
             'common cross-species fixed differences.',
             ha='center', va='top', fontsize=9, style='italic',
             color='#444444', wrap=True)

    fig.suptitle('Site Frequency Spectrum Comparison', fontsize=14,
                 fontweight='bold', y=1.02)
    plt.savefig(output_file, dpi=150, bbox_inches='tight', facecolor='white')
    plt.close()


def create_ibs_figure(all_data, sample_species, output_file):
    """Two-panel IBS figure:
       Left: per-dataset side-by-side violin plots of within and between clade pairs
       Right: bar chart of separation magnitude (within - between mean)
    """
    print(f"Creating: {output_file}", flush=True)
    fig = plt.figure(figsize=(20, 7.5))
    gs = GridSpec(1, 2, figure=fig, wspace=0.34, width_ratios=[1.35, 1.15])

    ax_main = fig.add_subplot(gs[0, 0])
    ax_summary = fig.add_subplot(gs[0, 1])

    # Compute IBS per dataset
    summary_rows = []
    dataset_distributions = []
    for data in all_data:
        samples = data['samples']
        n = len(samples)
        ibs_match = data['ibs_match']
        ibs_total = data['ibs_total']

        within_ibs = []
        between_ibs = []
        for i in range(n):
            clade_i = get_clade(sample_species.get(samples[i]))
            for j in range(i + 1, n):
                clade_j = get_clade(sample_species.get(samples[j]))
                if clade_i is None or clade_j is None:
                    continue
                if ibs_total[i, j] == 0:
                    continue
                ratio = ibs_match[i, j] / ibs_total[i, j]
                if clade_i == clade_j:
                    within_ibs.append(ratio)
                else:
                    between_ibs.append(ratio)

        dataset_distributions.append({
            'label': data['label'],
            'alt_of': data.get('alt_of'),
            'within': within_ibs,
            'between': between_ibs,
        })
        if within_ibs and between_ibs:
            wm, bm = np.mean(within_ibs), np.mean(between_ibs)
            summary_rows.append((data['label'], wm, bm, wm - bm,
                                 len(within_ibs), len(between_ibs),
                                 data.get('alt_of')))

    # ---- Left panel: side-by-side violin plots, colored by within/between ----
    n_inputs = len(dataset_distributions)
    if n_inputs > 0 and any(d['within'] for d in dataset_distributions):
        positions = []
        violin_data = []
        violin_colors = []
        x_labels = []
        for i, d in enumerate(dataset_distributions):
            base_x = i * 3
            if d['within']:
                positions.append(base_x)
                violin_data.append(d['within'])
                violin_colors.append(C_GOOD)
            if d['between']:
                positions.append(base_x + 1)
                violin_data.append(d['between'])
                violin_colors.append(C_ACCENT)
            label_text = get_label_pretty(d['label'])
            if d.get('alt_of') is not None:
                label_text = f"{label_text}\n(alt of {get_label_pretty(d['alt_of'])})"
            x_labels.append((base_x + 0.5, label_text))

        if violin_data:
            parts = ax_main.violinplot(violin_data, positions=positions,
                                       widths=0.85, showmeans=True,
                                       showmedians=False, showextrema=False)
            for pc, color in zip(parts['bodies'], violin_colors):
                pc.set_facecolor(color)
                pc.set_edgecolor('black')
                pc.set_linewidth(0.6)
                pc.set_alpha(0.8)
            if 'cmeans' in parts:
                parts['cmeans'].set_color('black')
                parts['cmeans'].set_linewidth(1.5)

        # X labels at midpoints
        ax_main.set_xticks([t[0] for t in x_labels])
        ax_main.set_xticklabels([t[1] for t in x_labels], fontsize=9,
                                rotation=25, ha='right')

        # Custom legend
        legend_handles = [
            Patch(color=C_GOOD, alpha=0.8, label='Within clade'),
            Patch(color=C_ACCENT, alpha=0.8, label='Between clade'),
        ]
        ax_main.legend(handles=legend_handles, fontsize=10, loc='lower left',
                       framealpha=0.95)
        ax_main.set_ylabel('IBS Sharing Fraction', fontsize=11)
        ax_main.set_title('IBS Distribution per Dataset (within vs between clades)',
                          fontsize=12, fontweight='bold')
        # Auto-zoom y-axis to relevant range based on actual data
        all_vals = []
        for d in dataset_distributions:
            all_vals.extend(d['within'])
            all_vals.extend(d['between'])
        if all_vals:
            data_min = min(all_vals)
            y_lo = max(0, data_min - 0.05)
        else:
            y_lo = 0
        ax_main.set_ylim(y_lo, 1.02)
        ax_main.grid(True, alpha=0.3, axis='y')
    else:
        ax_main.text(0.5, 0.5, 'No clade pair data\n(check species mapping)',
                     ha='center', va='center', transform=ax_main.transAxes,
                     fontsize=11, color='gray')
        ax_main.axis('off')

    # ---- Right panel: separation magnitude as colored bars ----
    if summary_rows:
        labels = []
        for r in summary_rows:
            lbl = get_label_pretty(r[0])
            if r[6] is not None:
                lbl = f"{lbl}\n(alt of {get_label_pretty(r[6])})"
            labels.append(lbl)
        seps = [r[3] for r in summary_rows]
        # Color bars by dataset's label color for consistency; alts get
        # a darker shade of the original's color
        bar_colors = [get_color(r[0], alt_of=r[6]) for r in summary_rows]

        x = np.arange(len(labels))
        bars = ax_summary.bar(x, seps, width=0.6, color=bar_colors,
                              alpha=0.85, edgecolor='black', lw=0.6)
        ax_summary.bar_label(
            bars, labels=[f'{sep:+.3f}' for sep in seps], padding=4,
            fontsize=9, fontweight='bold', rotation=0,
        )

        ax_summary.set_xticks(x)
        ax_summary.set_xticklabels(labels, fontsize=9, rotation=30, ha='right')
        ax_summary.set_ylabel('Within - Between mean IBS', fontsize=11)
        ax_summary.set_title('Clade Separation Magnitude',
                             fontsize=12, fontweight='bold')
        ax_summary.axhline(0, color='black', lw=0.8, zorder=2)
        ax_summary.grid(True, alpha=0.3, axis='y')
        # Round up to nearest 0.05 with headroom for labels
        max_abs = max(abs(s) for s in seps)
        ymax = max(0.05, np.ceil(max_abs * 1.65 / 0.05) * 0.05)
        ax_summary.set_ylim(-ymax * 0.2, ymax)
    else:
        ax_summary.text(0.5, 0.5, 'No data', ha='center', va='center',
                        transform=ax_summary.transAxes, fontsize=11, color='gray')
        ax_summary.axis('off')

    fig.text(0.5, -0.05,
             'IBS sharing fraction = fraction of comparison sites where two individuals share at '
             'least one allele. Within-clade pairs are within human, within Pan (chimp/bonobo/chinobo), '
             'or within orangutan. Larger within - between separation = the SNP set captures '
             'species-level differentiation more sharply.',
             ha='center', va='top', fontsize=9, style='italic',
             color='#444444', wrap=True)

    fig.suptitle('Identity-by-State Within vs Between Species Clades',
                 fontsize=14, fontweight='bold', y=1.02)
    plt.savefig(output_file, dpi=150, bbox_inches='tight', facecolor='white')
    plt.close()


def compute_species_pair_diff_rate(data, sample_species):
    """Per species pair: fraction of comparison sites where the pair has
    non-overlapping alleles."""
    samples = data['samples']
    n = len(samples)
    diff_count = data['diff_count']
    ibs_total = data['ibs_total']

    species_indices = defaultdict(list)
    for i, s in enumerate(samples):
        sp = sample_species.get(s)
        if sp:
            species_indices[sp].append(i)

    species_list = sorted(species_indices.keys())
    pair_diff = {}
    for ia, sp_a in enumerate(species_list):
        for sp_b in species_list[ia:]:
            if sp_a == sp_b:
                idx = species_indices[sp_a]
                if len(idx) < 2:
                    continue
                diff_sum = 0
                tot_sum = 0
                for i in idx:
                    for j in idx:
                        if i < j:
                            diff_sum += diff_count[i, j]
                            tot_sum += ibs_total[i, j]
                if tot_sum > 0:
                    pair_diff[(sp_a, sp_b)] = diff_sum / tot_sum
            else:
                idx_a = species_indices[sp_a]
                idx_b = species_indices[sp_b]
                diff_sum = 0
                tot_sum = 0
                for i in idx_a:
                    for j in idx_b:
                        diff_sum += diff_count[i, j]
                        tot_sum += ibs_total[i, j]
                if tot_sum > 0:
                    pair_diff[(sp_a, sp_b)] = diff_sum / tot_sum
    return pair_diff


def create_species_enrichment_figure(all_data, sample_species, output_file):
    """Two-section figure:
       Top: between-species pair-difference enrichment heatmaps (diagonal masked),
            shared color scale.
       Bottom: within-species variation as fold-change relative to baseline,
               grouped bar chart so all 4 datasets are comparable on one panel.
    """
    print(f"Creating: {output_file}", flush=True)
    baseline = None
    for d in all_data:
        if d['label'].lower() in ('input', 'all'):
            baseline = d
            break
    if baseline is None:
        baseline = all_data[0]
    print(f"  Using baseline: {baseline['label']}", flush=True)

    baseline_pairs = compute_species_pair_diff_rate(baseline, sample_species)
    if not baseline_pairs:
        print("  WARNING: Baseline has no species pair data", file=sys.stderr)
        return

    species_set = set()
    for sp_a, sp_b in baseline_pairs:
        species_set.add(sp_a)
        species_set.add(sp_b)
    clade_order = ['h', 'human', 'c', 'chimp', 'chimpanzee', 'b', 'bonobo',
                   'hy', 'chinobo', 'hybrid', 'o', 'orangutan', 'orang']
    order_idx = {sp: i for i, sp in enumerate(clade_order)}
    species_list = sorted(species_set, key=lambda s: order_idx.get(s.lower(), 99))
    species_pretty = [get_species_pretty(s) for s in species_list]
    nsp = len(species_list)
    sp_to_idx = {sp: i for i, sp in enumerate(species_list)}

    # Compute off-diagonal (between-species) enrichments per dataset
    n_inputs = len(all_data)
    enrichment_mats = []  # list of (label, matrix, alt_of) - off-diagonal only
    within_fold = []      # list of (label, sp -> fold_change, alt_of)

    for data in all_data:
        pair_diff = compute_species_pair_diff_rate(data, sample_species)
        mat = np.full((nsp, nsp), np.nan, dtype=np.float64)
        within_dict = {}
        for (sp_a, sp_b), val in pair_diff.items():
            base_val = baseline_pairs.get((sp_a, sp_b))
            if base_val is None or base_val == 0:
                continue
            ratio = val / base_val
            ia = sp_to_idx[sp_a]
            ib = sp_to_idx[sp_b]
            if sp_a == sp_b:
                within_dict[sp_a] = ratio
            else:
                mat[ia, ib] = ratio
                mat[ib, ia] = ratio
        enrichment_mats.append((data['label'], mat, data.get('alt_of')))
        within_fold.append((data['label'], within_dict, data.get('alt_of')))

    # Compute shared vmax across off-diagonal cells of non-baseline datasets
    all_off_diag = []
    for lbl, mat, _ in enrichment_mats:
        if lbl == baseline['label']:
            continue
        for i in range(nsp):
            for j in range(nsp):
                if i != j and not np.isnan(mat[i, j]):
                    all_off_diag.append(mat[i, j])
    vmax = max(max(all_off_diag) if all_off_diag else 2.0, 2.0)
    vmin = 0.0

    # Layout: 2 rows of heatmap panels (n_inputs cols + 1 for colorbar) + 1 wide bottom row
    cols = n_inputs
    fig = plt.figure(figsize=(5 * cols + 1, 5 + 5))
    gs = GridSpec(2, cols + 1, figure=fig, hspace=0.5, wspace=0.4,
                  height_ratios=[1.0, 0.8],
                  width_ratios=[1.0] * cols + [0.08])

    # ---- Top row: heatmaps ----
    last_im = None
    for idx, (label, mat, alt_of) in enumerate(enrichment_mats):
        ax = fig.add_subplot(gs[0, idx])
        masked = np.ma.masked_invalid(mat)
        im = ax.imshow(masked, cmap='viridis', aspect='equal',
                       vmin=vmin, vmax=vmax)
        last_im = im

        ax.set_xticks(np.arange(nsp))
        ax.set_yticks(np.arange(nsp))
        ax.set_xticklabels(species_pretty, rotation=45, ha='right', fontsize=9)
        ax.set_yticklabels(species_pretty, fontsize=9)

        for i in range(nsp):
            for j in range(nsp):
                if i == j:
                    ax.add_patch(plt.matplotlib.patches.Rectangle(
                        (j - 0.5, i - 0.5), 1, 1, facecolor='#dddddd',
                        edgecolor='white', zorder=2))
                    ax.text(j, i, '—', ha='center', va='center',
                            color='#888888', fontsize=10, zorder=3)
                elif not np.isnan(mat[i, j]):
                    val = mat[i, j]
                    text = f"{val:.2f}x"
                    color = "white" if val < (vmax * 0.5) else "black"
                    ax.text(j, i, text, ha="center", va="center",
                            color=color, fontsize=9, fontweight='bold')

        title = get_label_pretty(label)
        if alt_of is not None:
            title = f"{title}\n(alt of {get_label_pretty(alt_of)})"
        ax.set_title(title, fontsize=12, fontweight='bold')

    # Single shared colorbar in dedicated last column
    if last_im is not None:
        cax = fig.add_subplot(gs[0, cols])
        cbar = fig.colorbar(last_im, cax=cax)
        cbar.ax.set_ylabel(f'Fold enrichment vs {get_label_pretty(baseline["label"])}',
                           rotation=-90, va='bottom', fontsize=10)

    fig.text(0.5, 0.96,
             'Between-species pair difference enrichment (off-diagonal only)',
             ha='center', va='top', fontsize=12, fontweight='bold',
             color='#222222')

    # ---- Bottom row: within-species fold change as grouped bars ----
    ax_within = fig.add_subplot(gs[1, :])

    species_with_data = [sp for sp in species_list
                         if any(sp in w for _, w, _ in within_fold)]
    species_with_data_pretty = [get_species_pretty(s) for s in species_with_data]
    n_sp_with_data = len(species_with_data)
    bar_width = 0.85 / n_inputs
    x = np.arange(n_sp_with_data)

    # Track max height for legend placement and build legend handles manually
    max_height = 0.0
    legend_handles = []
    for offset, (label, w_dict, alt_of) in enumerate(within_fold):
        vals = [w_dict.get(sp, np.nan) for sp in species_with_data]
        positions = x + (offset - n_inputs / 2 + 0.5) * bar_width
        color = get_color(label, alt_of=alt_of)
        for pos, val in zip(positions, vals):
            if not np.isnan(val):
                bar_height = max(val, 0)
                ax_within.bar(pos, bar_height, width=bar_width * 0.9,
                              color=color, alpha=0.9,
                              edgecolor='black', lw=0.5)
                label_y = max(bar_height + 0.04, 0.06)
                ax_within.text(pos, label_y, f"{val:.2f}x",
                               ha='center', va='bottom', fontsize=8,
                               color='#222222', fontweight='bold')
                if bar_height > max_height:
                    max_height = bar_height
        legend_label = get_label_pretty(label)
        if alt_of is not None:
            legend_label = f"{legend_label} (alt of {get_label_pretty(alt_of)})"
        legend_handles.append(Patch(color=color, alpha=0.9,
                                     label=legend_label))

    ax_within.axhline(1.0, color='black', lw=0.8, linestyle='--',
                      alpha=0.6, zorder=1)
    ax_within.text(n_sp_with_data - 0.4, 1.02, '1.0x = baseline',
                   ha='right', va='bottom', fontsize=9, color='#666666',
                   style='italic')

    ax_within.set_xticks(x)
    ax_within.set_xticklabels(species_with_data_pretty, fontsize=10)
    ax_within.set_ylabel(f'Within-species fold change\nvs {get_label_pretty(baseline["label"])}',
                         fontsize=11)
    ax_within.set_title('Within-species variation per dataset',
                        fontsize=12, fontweight='bold')
    ax_within.legend(handles=legend_handles, fontsize=10, loc='upper center',
                     bbox_to_anchor=(0.5, -0.15), framealpha=0.95,
                     ncol=n_inputs, title='Dataset', title_fontsize=10)
    ax_within.grid(True, alpha=0.3, axis='y')
    # Set y-limit with headroom for labels
    ax_within.set_ylim(0, max(1.3, max_height * 1.2))

    fig.text(0.5, -0.02,
             'Top heatmaps: each off-diagonal cell is the fraction of sites where individuals from species A '
             'and B carry non-overlapping alleles, as fold change vs the baseline. Higher = better species '
             'discrimination. Diagonal cells (within-species) are scaled differently (the diagonal is '
             'within-species variation, not between-species), so are shown separately below. '
             'Bottom panel: within-species fold change. < 1.0 = selection has less within-species '
             'variation than the input (expected for species set); 1.0 = same; > 1.0 = enriched for '
             'within-species variants (expected for het set).',
             ha='center', va='top', fontsize=9, style='italic',
             color='#444444', wrap=True)

    fig.suptitle(f'Cross-Species Differentiation Enrichment '
                 f'(vs {get_label_pretty(baseline["label"])})',
                 fontsize=14, fontweight='bold', y=1.02)
    plt.savefig(output_file, dpi=150, bbox_inches='tight', facecolor='white')
    plt.close()




def build_rectangular_layout(edges, root, n_leaves, sample_species, samples):
    """Rooted rectangular layout with correct descendant-only internal y values."""
    adj = defaultdict(list)
    for u, v, length in edges:
        adj[u].append((v, length))
        adj[v].append((u, length))

    parent = {root: None}
    node_depth = {root: 0.0}
    order = [root]
    for node in order:
        for child, length in adj[node]:
            if child in parent:
                continue
            parent[child] = node
            node_depth[child] = node_depth[node] + length
            order.append(child)

    children = defaultdict(list)
    edge_length = {}
    for child, par in parent.items():
        if par is not None:
            children[par].append(child)
            for nbr, length in adj[par]:
                if nbr == child:
                    edge_length[(par, child)] = length
                    break

    descendant_cache = {}

    def descendant_leaves(node):
        if node in descendant_cache:
            return descendant_cache[node]
        if node < n_leaves:
            result = [node]
        else:
            result = []
            for child in children[node]:
                result.extend(descendant_leaves(child))
        descendant_cache[node] = result
        return result

    def child_sort_key(child):
        leaves = descendant_leaves(child)
        keys = []
        for leaf in leaves:
            sp = sample_species.get(samples[leaf], 'zzz')
            clade = get_clade(sp)
            keys.append((clade if clade is not None else 9, sp, samples[leaf]))
        return min(keys) if keys else (9, 'zzz', 'zzz')

    for node in list(children):
        children[node].sort(key=child_sort_key)

    leaf_order = []

    def visit(node):
        if node < n_leaves:
            leaf_order.append(node)
            return
        for child in children[node]:
            visit(child)

    visit(root)
    leaf_y = {leaf: float(i) for i, leaf in enumerate(leaf_order)}
    node_y = dict(leaf_y)
    for node in reversed(order):
        if node in node_y:
            continue
        leaves = descendant_leaves(node)
        node_y[node] = float(np.mean([leaf_y[x] for x in leaves])) if leaves else 0.0
    coords = {node: (node_depth[node], node_y[node]) for node in order}
    return coords, leaf_order, parent


def get_title_label(data):
    """Pretty label for use in plot titles, including alt-of suffix if any."""
    s = get_label_pretty(data['label'])
    if data.get('alt_of') is not None:
        s = f"{s} [alt of {get_label_pretty(data['alt_of'])}]"
    return s


def create_pairwise_heatmap(data, output_file):
    print(f"Creating: {output_file}", flush=True)
    samples = data['samples']
    n = len(samples)
    pair_matrix = data['pair_matrix']
    total_snps = data['total_sites']
    matrix_millions = pair_matrix / 1_000_000
    fig_size = max(12, n * 0.5)
    fig, ax = plt.subplots(figsize=(fig_size, fig_size))
    im = ax.imshow(matrix_millions, cmap='viridis', aspect='equal')
    cbar = fig.colorbar(im, ax=ax, shrink=0.8)
    cbar.ax.set_ylabel('Distinguishing SNPs (millions)', rotation=-90, va="bottom", fontsize=12)
    ax.set_xticks(np.arange(n))
    ax.set_yticks(np.arange(n))
    ax.set_xticklabels(samples, fontsize=9)
    ax.set_yticklabels(samples, fontsize=9)
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")
    for i in range(n):
        for j in range(n):
            if i == j:
                text, color = "-", "white"
            else:
                val = matrix_millions[i, j]
                text = f"{val:.1f}M" if val >= 1 else f"{int(val * 1000)}K"
                color = "white" if val < (matrix_millions.max() * 0.5) else "black"
            ax.text(j, i, text, ha="center", va="center", color=color, fontsize=7)
    ax.set_title(f'Pairwise Distinguishing SNPs - {get_title_label(data)}\n({total_snps:,} total)', fontsize=14)
    plt.tight_layout()
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    plt.close()


def plot_score_comparison(ax, datasets, field, title=''):
    all_valid = []
    for label, data, color in datasets:
        scores = data['scores'].get(field, np.array([]))
        valid = scores[~np.isnan(scores)]
        valid = valid[valid > 0]
        if len(valid) > 100:
            all_valid.append((label, np.log10(valid), color, len(valid)))
    if all_valid:
        all_concat = np.concatenate([v[1] for v in all_valid])
        bin_min, bin_max = np.percentile(all_concat, [1, 99])
        bins = np.linspace(bin_min, bin_max, 50)
        for label, log_vals, color, n in all_valid:
            ax.hist(log_vals, bins=bins, alpha=0.5, color=color,
                    label=label, density=True, edgecolor='none')
        ax.legend(fontsize=8, loc='upper right', ncol=2, framealpha=0.9)
        ax.set_xlabel(f'log10({field})', fontsize=9)
        ax.set_ylabel('Density', fontsize=9)
    else:
        ax.text(0.5, 0.5, f'No {field} data', ha='center', va='center',
                transform=ax.transAxes, fontsize=12, color='gray')
    ax.set_title(title, fontsize=10, fontweight='bold')


def plot_cdf_comparison(ax, datasets, field, title=''):
    has_data = False
    for label, data, color in datasets:
        scores = data['scores'].get(field, np.array([]))
        valid = scores[~np.isnan(scores)]
        valid = valid[valid > 0]
        if len(valid) > 100:
            has_data = True
            log_vals = np.log10(valid)
            sorted_vals = np.sort(log_vals)
            cdf = np.arange(1, len(sorted_vals) + 1) / len(sorted_vals)
            ax.plot(sorted_vals, cdf, color=color, lw=2, label=label)
    if has_data:
        ax.legend(fontsize=8, loc='lower right', framealpha=0.9)
        ax.set_xlabel(f'log10({field})', fontsize=9)
        ax.set_ylabel('Cumulative Fraction', fontsize=9)
        ax.grid(True, alpha=0.3)
    else:
        ax.text(0.5, 0.5, f'No {field} data', ha='center', va='center',
                transform=ax.transAxes, fontsize=12, color='gray')
    ax.set_title(title, fontsize=10, fontweight='bold')


def create_selection_figure(all_data, output_file):
    print(f"Creating: {output_file}", flush=True)

    def _legend_label(d):
        s = get_label_pretty(d['label'])
        if d.get('alt_of') is not None:
            s = f"{s} (alt of {get_label_pretty(d['alt_of'])})"
        return s

    datasets = [(_legend_label(d), d, get_color(d['label'], alt_of=d.get('alt_of')))
                for d in all_data]
    fig = plt.figure(figsize=(16, 10))
    gs = GridSpec(2, 3, figure=fig, hspace=0.4, wspace=0.35)

    ax1 = fig.add_subplot(gs[0, 0])
    plot_score_comparison(ax1, datasets, 'DEMUX_SCORE', 'DEMUX_SCORE (lower = rarer clade)')
    ax2 = fig.add_subplot(gs[0, 1])
    plot_score_comparison(ax2, datasets, 'HET_SCORE', 'HET_SCORE (higher = more hets)')
    ax3 = fig.add_subplot(gs[0, 2])
    plot_score_comparison(ax3, datasets, 'COV_SCORE', 'COV_SCORE (higher = more coverage)')
    ax4 = fig.add_subplot(gs[1, 0])
    plot_cdf_comparison(ax4, datasets, 'DEMUX_SCORE', 'DEMUX_SCORE CDF')
    ax5 = fig.add_subplot(gs[1, 1])
    plot_cdf_comparison(ax5, datasets, 'HET_SCORE', 'HET_SCORE CDF')
    ax6 = fig.add_subplot(gs[1, 2])
    plot_cdf_comparison(ax6, datasets, 'COV_SCORE', 'COV_SCORE CDF')

    fig.suptitle('Selection Quality Comparison', fontsize=14, fontweight='bold', y=0.99)
    plt.savefig(output_file, dpi=150, bbox_inches='tight', facecolor='white')
    plt.close()


def create_individuals_figure(data, output_file):
    print(f"Creating: {output_file}", flush=True)
    samples = data['samples']
    n_samples = len(samples)
    color = get_color(data['label'], alt_of=data.get('alt_of'))

    def short_name(s, max_len=12):
        return s[:max_len] if len(s) > max_len else s

    fig = plt.figure(figsize=(18, 12))
    gs = GridSpec(2, 2, figure=fig, hspace=0.4, wspace=0.25)
    x = np.arange(n_samples)

    # Panel 1: Het counts
    ax1 = fig.add_subplot(gs[0, 0])
    sort_idx = np.argsort(data['het_counts'])[::-1]
    sorted_samples = [short_name(samples[i]) for i in sort_idx]
    het_h = data['het_counts'][sort_idx]
    ax1.bar(x, het_h, color=color, alpha=0.85)
    ax1.axhline(np.median(het_h), color=C_ACCENT, linestyle='--', alpha=0.6, lw=1.5)
    ax1.set_xticks(x)
    ax1.set_xticklabels(sorted_samples, rotation=60, ha='right', fontsize=8)
    ax1.set_ylabel('Het Sites', fontsize=11)
    ax1.set_title(f'Per-Individual Het Counts - {get_title_label(data)}', fontsize=12, fontweight='bold')
    ax1.set_xlim(-0.5, n_samples - 0.5)

    # Panel 2: Genotype composition
    ax2 = fig.add_subplot(gs[0, 1])
    sort_idx2 = np.argsort(data['hom_alt_counts'])[::-1]
    sorted_samples2 = [short_name(samples[i]) for i in sort_idx2]
    called = data['hom_ref_counts'] + data['het_counts'] + data['hom_alt_counts']
    called = np.maximum(called, 1)
    hom_ref_frac = data['hom_ref_counts'][sort_idx2] / called[sort_idx2]
    het_frac = data['het_counts'][sort_idx2] / called[sort_idx2]
    hom_alt_frac = data['hom_alt_counts'][sort_idx2] / called[sort_idx2]
    ax2.bar(x, hom_ref_frac, label='Hom Ref', color=C_HOM_REF, alpha=0.85)
    ax2.bar(x, het_frac, bottom=hom_ref_frac, label='Het', color=C_HET_GT, alpha=0.85)
    ax2.bar(x, hom_alt_frac, bottom=hom_ref_frac + het_frac, label='Hom Alt', color=C_HOM_ALT, alpha=0.85)
    ax2.set_xticks(x)
    ax2.set_xticklabels(sorted_samples2, rotation=60, ha='right', fontsize=8)
    ax2.set_ylabel('Fraction', fontsize=11)
    ax2.set_title(f'Genotype Composition - {get_title_label(data)}', fontsize=12, fontweight='bold')
    ax2.legend(fontsize=8, loc='lower right', ncol=3, framealpha=0.9)
    ax2.set_ylim(0, 1)
    ax2.set_xlim(-0.5, n_samples - 0.5)

    # Panel 3: Missing rate
    ax3 = fig.add_subplot(gs[1, 0])
    total_per = data['hom_ref_counts'] + data['het_counts'] + data['hom_alt_counts'] + data['missing_counts']
    total_per = np.maximum(total_per, 1)
    miss_rate = 100 * data['missing_counts'] / total_per
    sort_idx3 = np.argsort(miss_rate)[::-1]
    sorted_samples3 = [short_name(samples[i]) for i in sort_idx3]
    miss_sorted = miss_rate[sort_idx3]
    colors_miss = [C_ACCENT if r > 5 else color for r in miss_sorted]
    ax3.barh(x, miss_sorted, color=colors_miss, alpha=0.85)
    ax3.axvline(5, color=C_ACCENT, linestyle='--', lw=1.5)
    ax3.set_yticks(x)
    ax3.set_yticklabels(sorted_samples3, fontsize=8)
    ax3.set_xlabel('Missing Rate (%)', fontsize=11)
    ax3.set_title(f'Missing Data - {get_title_label(data)}', fontsize=12, fontweight='bold')
    ax3.invert_yaxis()
    ax3.set_xlim(0, max(miss_sorted.max() * 1.1, 5.5))

    # Panel 4: Min pairwise distinguishing power
    ax4 = fig.add_subplot(gs[1, 1])
    pair_matrix = data['pair_matrix']
    min_pairs = []
    for i in range(n_samples):
        row = [(j, pair_matrix[i, j]) for j in range(n_samples) if i != j]
        if row:
            min_pairs.append(min(row, key=lambda value: value[1])[1])
        else:
            min_pairs.append(0)
    min_pairs = np.asarray(min_pairs)
    sort_idx4 = np.argsort(min_pairs)[::-1]
    sorted_samples4 = [short_name(samples[i]) for i in sort_idx4]
    min_sorted = min_pairs[sort_idx4]
    median_min = np.median(min_sorted)
    colors_pair = [C_ACCENT if value < median_min * 0.5 else color
                   for value in min_sorted]
    ax4.bar(x, min_sorted / 1e6, color=colors_pair, alpha=0.85)
    ax4.axhline(median_min / 1e6, color=C_GOOD, linestyle='--', lw=1.5)
    ax4.set_xticks(x)
    ax4.set_xticklabels(sorted_samples4, rotation=60, ha='right', fontsize=8)
    ax4.set_ylabel('Min Distinguishing SNPs (M)', fontsize=11)
    ax4.set_title(f'Min Pairwise Power - {get_title_label(data)}', fontsize=12, fontweight='bold')
    ax4.set_xlim(-0.5, n_samples - 0.5)

    fig.suptitle(f'Per-Individual Analysis - {get_title_label(data)}', fontsize=14, fontweight='bold', y=0.99)
    plt.savefig(output_file, dpi=150, bbox_inches='tight', facecolor='white')
    plt.close()




def create_species_pairs_figure_cached(data, output_file):
    pair_counts = data.get('pair_assigned_counts', {})
    discrim = data.get('pair_discrim_means', {})
    if not pair_counts and not discrim and not data.get('bin_count'):
        return False
    print(f'Creating: {output_file}', flush=True)
    fig = plt.figure(figsize=(16, 6))
    gs = GridSpec(1, 3, figure=fig, wspace=0.35)

    ax1 = fig.add_subplot(gs[0, 0])
    if pair_counts:
        pairs = sorted(pair_counts)
        counts = [pair_counts[x] for x in pairs]
        assigned_total = max(1, sum(counts))
        colors = plt.cm.Set3(np.linspace(0, 1, len(pairs)))
        ax1.barh(np.arange(len(pairs)), counts, color=colors)
        ax1.set_yticks(np.arange(len(pairs)))
        ax1.set_yticklabels(pairs, fontsize=9)
        ax1.invert_yaxis()
        for i, count in enumerate(counts):
            ax1.text(count, i,
                     f' {count:,} ({100 * count / assigned_total:.1f}%)',
                     va='center', fontsize=8)
    else:
        ax1.text(0.5, 0.5, 'No PAIR_ASSIGNED data', ha='center', va='center',
                 transform=ax1.transAxes, color='gray')
    ax1.set_xlabel('SNPs assigned')
    ax1.set_title('Assigned selection budget by species pair', fontweight='bold')

    ax2 = fig.add_subplot(gs[0, 1])
    finite = [(key.replace('PAIR_DISCRIM_', ''), value)
              for key, value in sorted(discrim.items()) if math.isfinite(value)]
    if finite:
        names = [x[0] for x in finite]
        values = [x[1] for x in finite]
        ax2.barh(np.arange(len(names)), values, color='#7B2D8E', alpha=0.8)
        ax2.set_yticks(np.arange(len(names)))
        ax2.set_yticklabels(names, fontsize=8)
        ax2.invert_yaxis()
        for i, value in enumerate(values):
            ax2.text(value, i, f' {value:.3f}', va='center', fontsize=8)
    else:
        ax2.text(0.5, 0.5, 'No PAIR_DISCRIM data', ha='center', va='center',
                 transform=ax2.transAxes, color='gray')
    ax2.set_xlabel('Mean discrimination score')
    ax2.set_title('Mean pair discrimination', fontweight='bold')

    ax3 = fig.add_subplot(gs[0, 2])
    total_sites = int(data['total_sites'])
    total_bins = int(data.get('bin_count', 0))
    mean_per_bin = total_sites / total_bins if total_bins else math.nan
    text = f'Total sites: {total_sites:,}\nUnique bins: {total_bins:,}'
    if math.isfinite(mean_per_bin):
        text += f'\nMean SNPs/bin: {mean_per_bin:.1f}'
    ax3.text(0.5, 0.5, text, ha='center', va='center',
             transform=ax3.transAxes, fontsize=14,
             bbox=dict(boxstyle='round', facecolor='#f0f0f0', alpha=0.8))
    ax3.set_title('Bin coverage', fontweight='bold')
    ax3.axis('off')

    fig.suptitle(f'Species-pair analysis - {get_title_label(data)}',
                 fontsize=14, fontweight='bold', y=1.02)
    plt.savefig(output_file, dpi=150, bbox_inches='tight', facecolor='white')
    plt.close()
    return True


def _hist_median_from_log_bins(hist, bins):
    total = int(np.sum(hist))
    if total <= 0:
        return math.nan
    idx = int(np.searchsorted(np.cumsum(hist), (total + 1) / 2.0))
    idx = min(idx, len(bins) - 2)
    center = 0.5 * (bins[idx] + bins[idx + 1])
    return 10.0 ** center


def create_annotation_figure(all_data, output_file):
    if not any(d.get('annotation') for d in all_data):
        return False
    print(f'Creating: {output_file}', flush=True)
    n_inputs = len(all_data)
    fig_width = max(16, 4.8 * n_inputs)
    fig = plt.figure(figsize=(fig_width, 14))
    gs = GridSpec(3, n_inputs, figure=fig, hspace=0.48, wspace=0.35)

    ax_stack = fig.add_subplot(gs[0, :])
    x = np.arange(n_inputs)
    labels = []
    exonic = []
    intronic = []
    intergenic = []
    for data in all_data:
        label = get_label_pretty(data['label'])
        if data.get('alt_of') is not None:
            label += f"\n(alt of {get_label_pretty(data['alt_of'])})"
        labels.append(label)
        ann = data['annotation']
        exonic.append(int(ann.get('exonic', 0)))
        intronic.append(int(ann.get('intronic', 0)))
        intergenic.append(int(ann.get('intergenic', 0)))
    ax_stack.bar(x, exonic, 0.65, label='Exonic', color=C_EXONIC)
    ax_stack.bar(x, intronic, 0.65, bottom=exonic, label='Intronic', color=C_INTRONIC)
    bottoms = np.asarray(exonic) + np.asarray(intronic)
    ax_stack.bar(x, intergenic, 0.65, bottom=bottoms,
                 label='Intergenic', color=C_INTERGENIC)
    for i, values in enumerate(zip(exonic, intronic, intergenic)):
        total = sum(values)
        running = 0
        colors = ('white', 'black', 'white')
        for value, color in zip(values, colors):
            if value and total:
                ax_stack.text(i, running + value / 2.0,
                              f'{100.0 * value / total:.1f}%',
                              ha='center', va='center', fontsize=8,
                              color=color, fontweight='bold')
            running += value
        regulatory = data.get('score_stats', {}).get('ATAC_REGULATORY', {})
        regulatory_n = int(regulatory.get('count', 0) or 0)
        regulatory_fraction = (
            float(regulatory.get('sum', 0.0)) / regulatory_n
            if regulatory_n else math.nan
        )
        top_label = f'{total:,}'
        if math.isfinite(regulatory_fraction):
            top_label += f'\nRegulatory prior: {100 * regulatory_fraction:.1f}%'
        ax_stack.text(i, total, top_label, ha='center', va='bottom', fontsize=8)
    ax_stack.set_xticks(x)
    ax_stack.set_xticklabels(labels)
    ax_stack.set_ylabel('Variants')
    ax_stack.set_title(
        'Genomic location classification (descriptive only; not an ATAC eligibility filter)',
        fontweight='bold')
    ax_stack.legend(loc='upper left')

    rel_centers = 0.5 * (GENE_REL_BINS[:-1] + GENE_REL_BINS[1:])
    dist_centers = 0.5 * (INTERGENIC_LOG_BINS[:-1] + INTERGENIC_LOG_BINS[1:])
    for idx, data in enumerate(all_data):
        ann = data['annotation']
        color = get_color(data['label'], alt_of=data.get('alt_of'))
        title = get_title_label(data)

        ax = fig.add_subplot(gs[1, idx])
        hist = ann['gene_rel_hist']
        ax.bar(rel_centers, hist, width=np.diff(GENE_REL_BINS),
               color=color, alpha=0.85, align='center')
        ax.axvline(0.25, color=C_5PRIME, linestyle='--', lw=1.2)
        ax.axvline(0.75, color=C_3PRIME, linestyle='--', lw=1.2)
        n_rel = int(np.sum(hist))
        ax.text(0.03, 0.96,
                f"5': {int(ann.get('gene_5prime', 0)):,}\n"
                f"Middle: {int(ann.get('gene_middle', 0)):,}\n"
                f"3': {int(ann.get('gene_3prime', 0)):,}",
                transform=ax.transAxes, va='top', fontsize=8,
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.85))
        ax.set_xlabel("Relative gene position (5' to 3')")
        ax.set_ylabel('Variants')
        ax.set_title(f'Position in genes - {title}\n(n={n_rel:,})',
                     fontsize=10, fontweight='bold')

        ax2 = fig.add_subplot(gs[2, idx])
        dhist = ann['intergenic_log_hist']
        ax2.bar(dist_centers, dhist, width=np.diff(INTERGENIC_LOG_BINS),
                color=color, alpha=0.85, align='center')
        median_bp = _hist_median_from_log_bins(dhist, INTERGENIC_LOG_BINS)
        if math.isfinite(median_bp):
            ax2.axvline(math.log10(median_bp), color='black', linestyle='--', lw=1.5)
            ax2.text(0.03, 0.95, f'Approx. median: {median_bp / 1000.0:.1f} kb',
                     transform=ax2.transAxes, va='top', fontsize=8,
                     bbox=dict(boxstyle='round', facecolor='white', alpha=0.85))
        ax2.set_xlabel('log10(distance to nearest gene, bp)')
        ax2.set_ylabel('Variants')
        ax2.set_title(f'Intergenic distance - {title}\n'
                      f"(n={int(ann.get('intergenic_distance_count', 0)):,})",
                      fontsize=10, fontweight='bold')

    fig.suptitle('GTF-based annotation analysis', fontsize=14,
                 fontweight='bold', y=0.995)
    plt.savefig(output_file, dpi=150, bbox_inches='tight', facecolor='white')
    plt.close()
    return True


def pca_summary(data, sample_species):
    gt = data['gt_subsample']
    samples = data['samples']
    if gt.shape[0] < 50:
        return {'pc1_var': math.nan, 'pc2_var': math.nan, 'pc1_species_r2': math.nan}
    x = gt.T.astype(np.float32).copy()
    for j in range(x.shape[1]):
        col = x[:, j]
        mask = col >= 0
        col[~mask] = np.mean(col[mask]) if np.any(mask) else 0.0
        x[:, j] = col
    x -= np.mean(x, axis=0, keepdims=True)
    try:
        u, s, _ = np.linalg.svd(x, full_matrices=False)
    except np.linalg.LinAlgError:
        return {'pc1_var': math.nan, 'pc2_var': math.nan, 'pc1_species_r2': math.nan}
    var = (s ** 2) / np.sum(s ** 2) if np.sum(s ** 2) > 0 else np.zeros_like(s)
    pc1 = (u * s)[:, 0] if len(s) else np.zeros(len(samples))
    groups = defaultdict(list)
    for i, sample in enumerate(samples):
        species = sample_species.get(sample)
        if species is not None:
            groups[species].append(pc1[i])
    r2 = math.nan
    if len(groups) >= 2:
        gm = np.mean(pc1)
        ss_between = sum(len(v) * (np.mean(v) - gm) ** 2 for v in groups.values())
        ss_total = np.sum((pc1 - gm) ** 2)
        r2 = float(ss_between / ss_total) if ss_total > 0 else 0.0
    return {
        'pc1_var': float(var[0]) if len(var) > 0 else math.nan,
        'pc2_var': float(var[1]) if len(var) > 1 else math.nan,
        'pc1_species_r2': r2,
    }


def ibs_summary(data, sample_species):
    within = []
    between = []
    samples = data['samples']
    for i in range(len(samples)):
        ci = get_clade(sample_species.get(samples[i]))
        for j in range(i + 1, len(samples)):
            cj = get_clade(sample_species.get(samples[j]))
            total = data['ibs_total'][i, j]
            if ci is None or cj is None or total <= 0:
                continue
            value = data['ibs_match'][i, j] / total
            (within if ci == cj else between).append(value)
    wm = float(np.mean(within)) if within else math.nan
    bm = float(np.mean(between)) if between else math.nan
    return {
        'within_ibs': wm,
        'between_ibs': bm,
        'ibs_separation': wm - bm if math.isfinite(wm) and math.isfinite(bm) else math.nan,
    }


def sfs_summary(data):
    sfs = data['sfs_counts']
    n_samples = len(data['samples'])
    maximum = 2 * n_samples
    folded = np.zeros(maximum // 2 + 1, dtype=np.int64)
    for i in range(maximum + 1):
        folded[min(i, maximum - i)] += sfs[i]
    total = int(np.sum(folded[1:]))
    if total == 0:
        return {k: math.nan for k in ('rare', 'low', 'mid', 'common')}
    maf = np.arange(1, len(folded)) / maximum
    values = folded[1:]
    return {
        'rare': float(np.sum(values[maf < 0.05]) / total),
        'low': float(np.sum(values[(maf >= 0.05) & (maf < 0.20)]) / total),
        'mid': float(np.sum(values[(maf >= 0.20) & (maf < 0.40)]) / total),
        'common': float(np.sum(values[maf >= 0.40]) / total),
    }


def region_policy_rows(all_data, mito_contig, mito_labels, numt_labels,
                       audit_labels, numt_checked=True):
    rows = []
    mito_labels = set(mito_labels)
    numt_labels = set(numt_labels)
    audit_labels = set(audit_labels)
    for data in all_data:
        label = data['label']
        mito_count = int(data.get('chrom_counts', {}).get(mito_contig, 0))
        non_mito = int(data['total_sites']) - mito_count
        numt = int(data.get('numt_overlap', 0))
        if label in audit_labels:
            role = 'audit_source'
            status = 'INFO_ONLY'
            reason = 'Audit/source VCF is allowed to retain all records'
        elif label in mito_labels:
            role = 'mitochondrial_panel'
            violations = []
            if non_mito != 0:
                violations.append(f'{non_mito} non-{mito_contig} records')
            status = 'PASS' if not violations else 'FAIL'
            reason = '; '.join(violations) if violations else f'All records are on {mito_contig}'
        elif label in numt_labels:
            role = 'numt_diagnostic_panel'
            violations = []
            if mito_count != 0:
                violations.append(f'{mito_count} {mito_contig} records')
            if numt_checked and numt != int(data['total_sites']):
                violations.append(
                    f'{int(data["total_sites"]) - numt} records do not overlap the NUMT BED'
                )
            status = 'PASS' if not violations else 'FAIL'
            reason = ('; '.join(violations) if violations else
                      'All records overlap the NUMT BED and none are mitochondrial')
        else:
            role = 'nuclear_panel'
            violations = []
            if mito_count != 0:
                violations.append(f'{mito_count} {mito_contig} records')
            if numt_checked and numt != 0:
                violations.append(f'{numt} NUMT-overlap records')
            if violations:
                status = 'FAIL'
                reason = '; '.join(violations)
            elif numt_checked:
                status = 'PASS'
                reason = 'No mitochondrial or NUMT-overlap records'
            else:
                status = 'PARTIAL'
                reason = f'No {mito_contig} records; NUMT overlap was not checked'
        rows.append({
            'label': label, 'role': role, 'total_sites': int(data['total_sites']),
            'mito_records': mito_count, 'non_mito_records': non_mito,
            'numt_overlap_records': numt, 'status': status, 'reason': reason,
        })
    return rows


def write_region_policy(rows, output_file):
    with open(output_file, 'wt') as handle:
        columns = ['label', 'role', 'total_sites', 'mito_records',
                   'non_mito_records', 'numt_overlap_records', 'status', 'reason']
        handle.write('\t'.join(columns) + '\n')
        for row in rows:
            handle.write('\t'.join(str(row[c]) for c in columns) + '\n')


def format_float(value, digits=4):
    return 'n/a' if value is None or not math.isfinite(float(value)) else f'{float(value):.{digits}f}'


def summarize_dataset(data, sample_species, mito_contig):
    pca = pca_summary(data, sample_species)
    ibs = ibs_summary(data, sample_species)
    sfs = sfs_summary(data)
    topology = topology_score(data['gt_subsample'], data['samples'], sample_species)
    ann = data.get('annotation', {})
    return {
        'label': data['label'],
        'alt_of': data.get('alt_of') or '',
        'total_sites': int(data['total_sites']),
        'n_samples': len(data['samples']),
        'mito_records': int(data.get('chrom_counts', {}).get(mito_contig, 0)),
        'numt_overlap_records': int(data.get('numt_overlap', 0)),
        'topology_accuracy': float(topology) if not np.isnan(topology) else math.nan,
        **pca,
        **ibs,
        **{f'sfs_{k}': v for k, v in sfs.items()},
        'exonic': int(ann.get('exonic', 0)),
        'intronic': int(ann.get('intronic', 0)),
        'intergenic': int(ann.get('intergenic', 0)),
        'regulatory_prior_fraction': (
            float(data.get('score_stats', {}).get('ATAC_REGULATORY', {}).get('sum', 0.0)) /
            int(data.get('score_stats', {}).get('ATAC_REGULATORY', {}).get('count', 0))
            if int(data.get('score_stats', {}).get('ATAC_REGULATORY', {}).get('count', 0))
            else math.nan
        ),
    }


def write_summary_tsv(summaries, output_file):
    columns = [
        'label', 'alt_of', 'total_sites', 'n_samples', 'mito_records',
        'numt_overlap_records', 'topology_accuracy', 'pc1_var', 'pc2_var',
        'pc1_species_r2', 'within_ibs', 'between_ibs', 'ibs_separation',
        'sfs_rare', 'sfs_low', 'sfs_mid', 'sfs_common',
        'exonic', 'intronic', 'intergenic', 'regulatory_prior_fraction',
    ]
    with open(output_file, 'wt') as handle:
        handle.write('\t'.join(columns) + '\n')
        for row in summaries:
            values = []
            for column in columns:
                value = row[column]
                if isinstance(value, float):
                    values.append('NA' if not math.isfinite(value) else f'{value:.8g}')
                else:
                    values.append(str(value))
            handle.write('\t'.join(values) + '\n')


def write_combined_report(all_data, sample_species, summaries, policy_rows,
                          output_file, run_meta):
    by_label = {row['label']: row for row in summaries}
    lines = []
    lines.append('=' * 78)
    lines.append('UNIFIED SNP-PANEL QC REPORT')
    lines.append('=' * 78)
    lines.append(f"Generated: {_datetime.datetime.now().isoformat(timespec='seconds')}")
    lines.append(f"Script: {VERSION}")
    lines.append(f"GTF: {run_meta.get('gtf') or 'not used'}")
    lines.append(f"NUMT BED: {run_meta.get('numts_bed') or 'not used'}")
    lines.append(
        f"ATAC regulatory prior: {run_meta.get('atac_regulatory_source') or 'not recorded'}"
    )
    if run_meta.get('atac_regulatory_required'):
        lines.append(
            'ATAC policy: coverage plus the explicitly requested regulatory-region '
            'restriction are hard eligibility gates; gene bodies are not required.'
        )
    else:
        lines.append(
            'ATAC policy: coverage is the hard eligibility gate; gene/regulatory '
            'annotations are ranking/display metadata, not hard filters.'
        )
    lines.append(
        f"Normal-scale sample plots exclude: {run_meta.get('outgroup_sample') or 'none'}; "
        "full sample plots are written under with_orang/."
    )
    lines.append('')

    overlap_summary = run_meta.get('rna_atac_disjointness', {})
    lines.append('RNA / ATAC LOCUS DISJOINTNESS')
    lines.append('-' * 78)
    lines.append(
        f"  Status: {overlap_summary.get('status', 'NOT_TESTED')}; "
        f"RNA union={int(overlap_summary.get('rna_union_loci', 0)):,}; "
        f"ATAC union={int(overlap_summary.get('atac_union_loci', 0)):,}; "
        f"shared loci={int(overlap_summary.get('shared_loci', 0)):,}; "
        f"Jaccard={format_float(overlap_summary.get('jaccard'), 8)}"
    )
    lines.append('  RNA labels: ' + ', '.join(overlap_summary.get('rna_labels', [])))
    lines.append('  ATAC labels: ' + ', '.join(overlap_summary.get('atac_labels', [])))
    lines.append('')

    lines.append('DATASET SUMMARY')
    lines.append('-' * 78)
    for data in all_data:
        row = by_label[data['label']]
        alt = f" [alt of {data['alt_of']}]" if data.get('alt_of') else ''
        lines.append(f"{data['label']}{alt}")
        lines.append(f"  Path: {data.get('_vcf_file', '')}")
        lines.append(f"  Variants: {data['total_sites']:,}; samples: {len(data['samples'])}")
        lines.append(
            '  Structure: topology=' + format_float(row['topology_accuracy'])
            + ', PC1 species R^2=' + format_float(row['pc1_species_r2'])
            + ', IBS separation=' + format_float(row['ibs_separation'])
        )
        lines.append(
            '  SFS fractions: rare=' + format_float(row['sfs_rare'])
            + ', low=' + format_float(row['sfs_low'])
            + ', mid=' + format_float(row['sfs_mid'])
            + ', common=' + format_float(row['sfs_common'])
        )
        if data.get('annotation'):
            total_ann = row['exonic'] + row['intronic'] + row['intergenic']
            if total_ann:
                lines.append(
                    f"  Annotation: exonic={row['exonic']:,} "
                    f"({100 * row['exonic'] / total_ann:.2f}%), "
                    f"intronic={row['intronic']:,} "
                    f"({100 * row['intronic'] / total_ann:.2f}%), "
                    f"intergenic={row['intergenic']:,} "
                    f"({100 * row['intergenic'] / total_ann:.2f}%)"
                )
            top_types = data['annotation'].get('gene_type_counts', Counter()).most_common(8)
            if top_types:
                lines.append('  Top gene types: ' + ', '.join(f'{k}={v:,}' for k, v in top_types))
        if math.isfinite(row.get('regulatory_prior_fraction', math.nan)):
            lines.append(
                '  Regulatory-prior overlap: '
                f"{100 * row['regulatory_prior_fraction']:.2f}%"
            )
        if data.get('score_stats'):
            lines.append('  INFO score statistics (exact mean; median from stored stride sample):')
            for field in sorted(data['score_stats']):
                stat = data['score_stats'][field]
                mean = stat['sum'] / stat['count'] if stat['count'] else math.nan
                sample = data['scores'].get(field, np.zeros(0))
                median = float(np.median(sample)) if len(sample) else math.nan
                lines.append(
                    f"    {field}: n={int(stat['count']):,}, mean={format_float(mean, 6)}, "
                    f"sample_median={format_float(median, 6)}, sample_n={len(sample):,}"
                )
        samples = data['samples']
        pairs = []
        matrix = data['pair_matrix']
        for i in range(len(samples)):
            for j in range(i + 1, len(samples)):
                pairs.append((int(matrix[i, j]), samples[i], samples[j]))
        pairs.sort()
        if pairs:
            lines.append('  Hardest genotype pairs:')
            for count, a, b in pairs[:10]:
                lines.append(f'    {a} vs {b}: {count:,} distinguishing genotypes')
        lines.append('')

    lines.append('MITOCHONDRIAL / NUMT REGION POLICY')
    lines.append('-' * 78)
    for row in policy_rows:
        lines.append(
            f"  {row['label']}: {row['status']} ({row['role']}) - {row['reason']}"
        )
    lines.append('')
    lines.append('INTERPRETATION NOTES')
    lines.append('-' * 78)
    lines.append('  - Input/all audit VCFs may retain chrM and NUMT-overlap records.')
    lines.append('  - Nuclear panels should contain zero chrM and zero NUMT-overlap records.')
    lines.append('  - The dedicated mitochondrial panel should contain chrM records only.')
    lines.append('  - Score distribution plots use a deterministic per-contig stride sample;')
    lines.append('    exact score counts and means are retained in this report.')
    lines.append('  - The intergenic median shown in the annotation plot is histogram-based.')
    lines.append('')
    with open(output_file, 'wt') as handle:
        handle.write('\n'.join(lines) + '\n')


def render_outputs(all_data, sample_species, prefix, compact,
                   pair_comparison_plots, run_meta, outgroup_sample=None):
    require_matplotlib()
    os.makedirs(os.path.dirname(prefix) or '.', exist_ok=True)
    outputs = []

    diagnostic_labels = set(run_meta.get('mito_labels', [])) | set(
        run_meta.get('numt_labels', []))
    analytical_data = [d for d in all_data if d['label'] not in diagnostic_labels]
    if not analytical_data:
        analytical_data = all_data
    has_outgroup = bool(outgroup_sample) and any(
        outgroup_sample in d['samples'] for d in analytical_data)
    display_data = (
        [subset_dataset_samples(d, [outgroup_sample]) for d in analytical_data]
        if has_outgroup else analytical_data
    )

    combined = [
        (create_structure_figure, (display_data, sample_species,
                                   f'{prefix}_structure.png')),
        # SFS and selection/annotation are site-level summaries and therefore
        # stay exact; only sample-dependent figures are rescaled without JOS3C1.
        (create_sfs_figure, (analytical_data, f'{prefix}_sfs.png')),
        (create_ibs_figure, (display_data, sample_species, f'{prefix}_ibs.png')),
        (create_species_enrichment_figure,
         (display_data, sample_species, f'{prefix}_species_enrichment.png')),
        (create_selection_figure, (analytical_data, f'{prefix}_selection_scores.png')),
    ]
    for func, args in combined:
        func(*args)
        outputs.append(args[-1])
    if run_meta.get('gtf'):
        ann_path = f'{prefix}_annotation.png'
        if create_annotation_figure(analytical_data, ann_path):
            outputs.append(ann_path)

    overlap_path = f'{prefix}_rna_atac_overlap.png'
    if create_overlap_figure(
            run_meta.get('locus_overlap_rows', []),
            run_meta.get('rna_atac_disjointness', {}), overlap_path):
        outputs.append(overlap_path)

    if not compact:
        by_label_full = {d['label']: d for d in analytical_data}
        for data in display_data:
            label = data['label']
            path = f'{prefix}_individuals_{label}.png'
            create_individuals_figure(data, path)
            outputs.append(path)
            path = f'{prefix}_pairwise_{label}.png'
            create_pairwise_heatmap(data, path)
            outputs.append(path)
            full_data = by_label_full[label]
            if 'species' in label or full_data.get('pair_assigned_counts'):
                path = f'{prefix}_species_pairs_{label}.png'
                if create_species_pairs_figure_cached(full_data, path):
                    outputs.append(path)

    if has_outgroup:
        prefix_path = Path(prefix)
        with_orang_dir = prefix_path.parent / 'with_orang'
        with_orang_dir.mkdir(parents=True, exist_ok=True)
        full_prefix = str(with_orang_dir / prefix_path.name)
        for func, call_args in [
            (create_structure_figure,
             (analytical_data, sample_species, f'{full_prefix}_structure.png')),
            (create_ibs_figure,
             (analytical_data, sample_species, f'{full_prefix}_ibs.png')),
            (create_species_enrichment_figure,
             (analytical_data, sample_species,
              f'{full_prefix}_species_enrichment.png')),
        ]:
            func(*call_args)
            outputs.append(call_args[-1])
        if not compact:
            for data in analytical_data:
                label = data['label']
                path = f'{full_prefix}_individuals_{label}.png'
                create_individuals_figure(data, path)
                outputs.append(path)
                path = f'{full_prefix}_pairwise_{label}.png'
                create_pairwise_heatmap(data, path)
                outputs.append(path)

    if pair_comparison_plots:
        primary = {d['label']: d for d in display_data if d.get('alt_of') is None}
        for alt in [d for d in display_data if d.get('alt_of') is not None]:
            if alt['alt_of'] not in primary:
                continue
            pair = [primary[alt['alt_of']], alt]
            base = f"{prefix}_paircompare_{alt['alt_of']}_vs_{alt['label']}"
            for func, args in [
                (create_structure_figure, (pair, sample_species, f'{base}_structure.png')),
                (create_sfs_figure, (pair, f'{base}_sfs.png')),
                (create_ibs_figure, (pair, sample_species, f'{base}_ibs.png')),
                (create_species_enrichment_figure,
                 (pair, sample_species, f'{base}_species_enrichment.png')),
                (create_selection_figure, (pair, f'{base}_selection_scores.png')),
            ]:
                func(*args)
                outputs.append(args[-1])
            if run_meta.get('gtf'):
                path = f'{base}_annotation.png'
                if create_annotation_figure(pair, path):
                    outputs.append(path)
    return outputs


def build_parser():
    parser = argparse.ArgumentParser(
        description='Unified SNP-panel QC: one BCF pass, one cache, one runner.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument('--input', action='append', nargs='+', default=[],
                        metavar='VALUE',
                        help='Primary input. Repeatable. Accepts LABEL PATH or one quoted "LABEL PATH".')
    parser.add_argument('--alt-input', action='append', nargs='+', default=[],
                        metavar='VALUE',
                        help='Paired alternative input. Repeatable.')
    parser.add_argument('--panel-metadata',
                        help='TSV whose first two columns are sample ID and species.')
    parser.add_argument('--gtf', default=DEFAULT_GTF,
                        help='GTF/GTF.GZ for annotation analysis.')
    parser.add_argument('--atac-regulatory-bed', default=None,
                        help='Optional ATAC regulatory BED; otherwise report the default strand-aware 2-kb promoter prior.')
    parser.add_argument('--atac-regulatory-required', action='store_true',
                        help='Record that the regulatory prior was used as a hard ATAC filter.')
    parser.add_argument('--skip-annotation', action='store_true',
                        help='Do not load the GTF or compute annotation QC.')
    parser.add_argument('--numts-bed', default=DEFAULT_NUMTS_BED,
                        help='Reference-specific nuclear NUMT BED for panel-policy QC.')
    parser.add_argument('--skip-numt-check', action='store_true',
                        help='Do not count variants overlapping the NUMT BED.')
    parser.add_argument('--prefix', default='snp_qc', help='Output path prefix.')
    parser.add_argument('--cache', help='NPZ cache path; default is PREFIX_cache.npz.')
    parser.add_argument('--threads', type=int, default=max(1, multiprocessing.cpu_count()),
                        help='Maximum parallel contig workers per BCF.')
    parser.add_argument('--subsample-stride', type=int, default=1000,
                        help='Keep every Nth site per contig for PCA and tree.')
    parser.add_argument('--score-sample-stride', type=int, default=100,
                        help='Keep every Nth site per contig for score distributions.')
    parser.add_argument('--cache-only', action='store_true',
                        help='Build the shared cache but do not render plots.')
    parser.add_argument('--plot-only', action='store_true',
                        help='Read --cache and render without reopening BCFs.')
    parser.add_argument('--compact', action='store_true',
                        help='Create combined figures only; skip per-input detail figures.')
    parser.add_argument('--pair-comparison-plots', action='store_true',
                        help='Also create isolated primary-vs-alt figure sets.')
    parser.add_argument('--harmonize-samples', action='store_true',
                        help='Reheader alt inputs, intersect samples, and subset all inputs first.')
    parser.add_argument('--rename-sample', action='append', default=[], metavar='OLD=NEW',
                        help='Sample rename applied to alternative BCFs during harmonization.')
    parser.add_argument('--prepared-dir',
                        help='Harmonized BCF output directory; default PREFIX_prepared.')
    parser.add_argument('--min-intersection', type=int, default=5,
                        help='Minimum shared sample count for harmonization.')
    parser.add_argument('--index-missing', action='store_true',
                        help='Create missing CSI indexes in place.')
    parser.add_argument('--mito-contig', default='chrM')
    parser.add_argument('--mito-panel-label', action='append', default=[],
                        help='Dataset label treated as dedicated mitochondrial panel. Repeatable.')
    parser.add_argument('--numt-panel-label', action='append', default=[],
                        help='Dataset label treated as dedicated NUMT diagnostic panel. Repeatable.')
    parser.add_argument('--audit-label', action='append', default=[],
                        help='Dataset label allowed to retain all regions. Repeatable.')
    parser.add_argument('--allow-region-policy-failures', action='store_true',
                        help='Report chrM/NUMT policy failures but exit zero.')
    parser.add_argument('--require-rna-atac-disjoint', action='store_true',
                        help='Fail unless RNA and ATAC panel unions share exactly zero loci.')
    parser.add_argument('--outgroup-sample', default='JOS3C1',
                        help='Sample excluded from normal-scale sample-dependent plots; full plots go under with_orang/.')
    parser.add_argument('--force', action='store_true',
                        help='Allow replacement of prepared BCFs and output cache.')
    return parser


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.cache_only and args.plot_only:
        parser.error('--cache-only and --plot-only are mutually exclusive')
    if args.threads < 1:
        parser.error('--threads must be >= 1')
    if args.subsample_stride < 1 or args.score_sample_stride < 1:
        parser.error('stride arguments must be >= 1')

    cache_file = os.path.abspath(args.cache or f'{args.prefix}_cache.npz')
    prefix = os.path.abspath(args.prefix)
    mito_labels = args.mito_panel_label or ['mito', 'mt']
    numt_labels = args.numt_panel_label or []
    audit_labels = args.audit_label or ['input', 'all']

    if not args.plot_only and os.path.exists(cache_file) and not args.force:
        parser.error(f'Cache exists; pass --force to replace: {cache_file}')

    if args.plot_only:
        require_numpy()
        if not os.path.isfile(cache_file):
            parser.error(f'Cache not found: {cache_file}')
        all_data, sample_species, run_meta = load_cache(cache_file)
    else:
        require_numpy()
        require_executable('bcftools')
        if not args.input:
            parser.error('At least one --input is required unless --plot-only is used')
        if not args.panel_metadata:
            parser.error('--panel-metadata is required unless --plot-only is used')
        specs = normalize_input_specs(args.input, args.alt_input)
        for spec in specs:
            if not os.path.isfile(spec['path']):
                parser.error(f"Input file not found: {spec['path']}")
        rename_map = parse_rename_map(args.rename_sample)
        if args.harmonize_samples:
            prepared_dir = args.prepared_dir or f'{prefix}_prepared'
            specs = harmonize_inputs(
                specs, rename_map, prepared_dir, args.threads,
                args.min_intersection, args.force,
            )
        elif rename_map:
            eprint('WARNING: --rename-sample has no effect without --harmonize-samples')

        for spec in specs:
            ensure_bcf_index(spec['path'], args.threads, args.index_missing)
        sample_species = load_panel_metadata(args.panel_metadata)
        gene_index = None
        gtf_path = None
        if not args.skip_annotation:
            gtf_path = os.path.abspath(args.gtf)
            gene_index = load_gtf(gtf_path)
        numts_path = None
        numt_intervals = {}
        if not args.skip_numt_check:
            numts_path = os.path.abspath(args.numts_bed)
            numt_intervals = load_numt_bed(numts_path)
            print(f'Loaded NUMT intervals on {len(numt_intervals)} contigs', flush=True)

        all_data = []
        for spec in specs:
            all_data.append(analyze_vcf(
                spec['path'], spec['label'], spec['alt_of'], sample_species,
                gene_index, numt_intervals, args.subsample_stride,
                args.score_sample_stride, args.threads,
            ))
        overlap_rows, disjointness = compute_locus_overlap(all_data)
        run_meta = {
            'panel_metadata': os.path.abspath(args.panel_metadata),
            'gtf': gtf_path,
            'numts_bed': numts_path,
            'mito_contig': args.mito_contig,
            'mito_labels': mito_labels,
            'numt_labels': numt_labels,
            'outgroup_sample': args.outgroup_sample,
            'atac_regulatory_source': (
                os.path.abspath(args.atac_regulatory_bed)
                if args.atac_regulatory_bed
                else 'GTF_2KB_UPSTREAM_PROMOTERS'
            ),
            'atac_regulatory_required': bool(args.atac_regulatory_required),
            'audit_labels': audit_labels,
            'subsample_stride': args.subsample_stride,
            'score_sample_stride': args.score_sample_stride,
            'locus_overlap_rows': overlap_rows,
            'rna_atac_disjointness': disjointness,
        }
        for data in all_data:
            data.pop('_loci_by_chrom', None)
        save_cache(cache_file, all_data, sample_species, run_meta)

    # Command-line policy labels override cached defaults when explicitly supplied.
    mito_contig = args.mito_contig or run_meta.get('mito_contig', 'chrM')
    if not args.mito_panel_label:
        mito_labels = run_meta.get('mito_labels', mito_labels)
    if not args.numt_panel_label:
        numt_labels = run_meta.get('numt_labels', numt_labels)
    if not args.audit_label:
        audit_labels = run_meta.get('audit_labels', audit_labels)

    policy_rows = region_policy_rows(
        all_data, mito_contig, mito_labels, numt_labels, audit_labels,
        numt_checked=bool(run_meta.get('numts_bed')),
    )
    summaries = [summarize_dataset(d, sample_species, mito_contig) for d in all_data]
    policy_file = f'{prefix}_region_policy.tsv'
    summary_file = f'{prefix}_summary.tsv'
    report_file = f'{prefix}_report.txt'
    overlap_file = f'{prefix}_rna_atac_overlap.tsv'
    disjointness_file = f'{prefix}_rna_atac_disjointness.tsv'
    write_region_policy(policy_rows, policy_file)
    write_summary_tsv(summaries, summary_file)
    write_overlap_tsv(run_meta.get('locus_overlap_rows', []), overlap_file)
    write_disjointness_tsv(
        run_meta.get('rna_atac_disjointness', {}), disjointness_file)
    write_combined_report(all_data, sample_species, summaries, policy_rows,
                          report_file, run_meta)

    rendered = []
    if not args.cache_only:
        rendered = render_outputs(
            all_data, sample_species, prefix, args.compact,
            args.pair_comparison_plots, run_meta,
            outgroup_sample=args.outgroup_sample or run_meta.get('outgroup_sample'),
        )

    failed_policy = [row for row in policy_rows if row['status'] == 'FAIL']
    print('\n' + '=' * 70)
    print('Unified QC complete')
    print('=' * 70)
    print(f'Cache:         {cache_file}')
    print(f'Summary:       {summary_file}')
    print(f'Report:        {report_file}')
    print(f'Region policy: {policy_file}')
    print(f'Locus overlap: {overlap_file}')
    print(f'Disjointness:  {disjointness_file}')
    if rendered:
        print(f'Figures:       {len(rendered)}')
        for path in rendered:
            print(f'  {path}')
    if failed_policy:
        print('Region-policy failures:')
        for row in failed_policy:
            print(f"  {row['label']}: {row['reason']}")
        if not args.allow_region_policy_failures:
            return 2
    disjointness = run_meta.get('rna_atac_disjointness', {})
    if args.require_rna_atac_disjoint and disjointness.get('status') != 'PASS':
        print(
            'RNA/ATAC disjointness failure: '
            f"status={disjointness.get('status', 'NOT_TESTED')}, "
            f"shared_loci={int(disjointness.get('shared_loci', 0))}"
        )
        return 3
    return 0


if __name__ == '__main__':
    try:
        raise SystemExit(main())
    except KeyboardInterrupt:
        eprint('\nInterrupted')
        raise SystemExit(130)
    except Exception as exc:
        eprint(f'ERROR: {exc}')
        if os.environ.get('SNP_QC_TRACEBACK') == '1':
            import traceback
            traceback.print_exc()
        raise SystemExit(1)


# Revision history
# V1_R1 (2026-08-02): Consolidated five QC/preparation scripts into one
# standalone runner. Shared one-pass genotype/score/annotation accumulation;
# optional alt-BCF reheader and sample harmonization; cache-only and plot-only
# modes; cached species-pair metrics; chrM/NUMT policy enforcement; compact
# output mode; optional isolated pair-comparison plots. Fixed fail-open worker
# handling, nested-gene/exon overlap lookup, and rooted tree internal-node
# layout while preserving the existing QC metrics and figure families.
# 1.1.0 (2026-08-28): Adapted the consolidated RNA panel QC for joint
# RNA/ATAC validation; added modality-aware labels, exact position-level
# pairwise overlap tables, RNA/ATAC union disjointness, overlap heatmap, and a
# zero-shared-locus failure gate.
