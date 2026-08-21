#!/usr/bin/env python3
"""CellBouncer identity-reconciliation command line.

Subcommands keep the SLURM stages independent without deploying a separate
Python executable for every small operation. Ambient RNA is intentionally not
an evidence channel. Decision thresholds are versioned in
identity_reconciliation_common.py rather than a separate YAML file.
"""
from __future__ import annotations

import sys


# -----------------------------------------------------------------------------
# metadata
# -----------------------------------------------------------------------------

import argparse
import json
import os
import re
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Tuple

from identity_reconciliation_common import (
    SCHEMA_VERSION, canonical_genotype, canonical_uid_set, clean, json_dump_atomic,
    natural_key, sha256_file, write_tsv,
)

EXPECTED_FIELDS = [
    "library", "canonical_genotype", "donor_components", "expected_ploidy_class",
    "n_metadata_rows", "uid_candidate_count", "uid_candidates", "reconciled_uid",
    "uid_resolution_status", "uid_resolution_basis", "uid_resolution_scope",
    "uid_member_details", "fz_batches", "corrected_fzgrps", "line_labels", "source_sheets",
]
UID_MEMBER_FIELDS = [
    "library", "uid_resolution_scope", "canonical_genotype", "uid", "ordered_wgs_key",
    "ordered_line_label", "denotion_r1", "fz_batch", "corrected_fzgrp", "source_sheet", "source_row",
]
WARNING_FIELDS = ["library", "canonical_genotype", "warning", "detail", "source_sheet", "source_row"]
GLOBAL_LINE_FIELDS = [
    "canonical_genotype", "donor_components", "biological_ploidy_class",
    "n_line_meta_rows", "uid_candidates", "line_labels", "fz_batches",
    "corrected_fzgrps", "source_sheet", "source_workbook_sha256",
]
GLOBAL_DONOR_FIELDS = [
    "donor_id", "n_global_biological_lines", "global_biological_lines",
    "source_sheet", "source_workbook_sha256",
]


def metadata_parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--xlsx", required=True)
    p.add_argument("--panel-metadata", required=True)
    p.add_argument("--output-root", required=True)
    return p.parse_args()


def norm_col(df, *names):
    by = {re.sub(r"[^a-z0-9]", "", str(c).lower()): c for c in df.columns}
    for name in names:
        k = re.sub(r"[^a-z0-9]", "", name.lower())
        if k in by:
            return by[k]
    return None


def text(v):
    try:
        import pandas as pd
        if pd.isna(v):
            return ""
    except Exception:
        pass
    return clean(v)


def row_value(row, *names):
    for name in names:
        if name in row.index:
            v = text(row[name])
            if v:
                return v
    return ""


def build_aliases(convert, convert_extra):
    label_col = norm_col(convert, "Library", "library_label", "line")
    vcf_col = norm_col(convert, "VCF_ID", "vcf_id", "canonical_vcf_id")
    if not label_col or not vcf_col:
        raise SystemExit("convert sheet requires Library and VCF_ID columns")
    aliases: Dict[str, str] = {}
    sources: Dict[str, str] = {}
    supplemental = defaultdict(set)
    spsx_col = norm_col(convert, "SPSX", "spsx")
    spsx: Dict[str, str] = {}
    for _, row in convert.iterrows():
        label = text(row[label_col]); vcf = text(row[vcf_col])
        if label and vcf:
            aliases[label] = vcf; aliases.setdefault(vcf, vcf); sources[label] = "convert"
            if spsx_col:
                spsx[label] = text(row[spsx_col])
    if convert_extra is not None:
        ex_label = norm_col(convert_extra, "Library", "library_label", "line")
        ex_vcf = norm_col(convert_extra, "VCF_ID", "vcf_id", "canonical_vcf_id")
        if ex_label and ex_vcf:
            for _, row in convert_extra.iterrows():
                label = text(row[ex_label]); vcf = text(row[ex_vcf])
                if not label or not vcf:
                    continue
                if label not in aliases:
                    aliases[label] = vcf; sources[label] = "convert_extra"
                elif aliases[label] != vcf:
                    supplemental[label].add(vcf)
                aliases.setdefault(vcf, aliases.get(label, vcf))
    return aliases, sources, supplemental, spsx


def load_panel_ids(path: str):
    import csv
    ids = set()
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if reader.fieldnames:
            id_col = None
            for c in reader.fieldnames:
                if re.sub(r"[^a-z0-9]", "", c.lower()) in {"individ", "vcfid", "sample", "id"}:
                    id_col = c; break
            id_col = id_col or reader.fieldnames[0]
            for row in reader:
                v = clean(row.get(id_col))
                if v: ids.add(v)
    return ids


def infer_libnum(value: str):
    m = re.search(r"(?:lib|RNA[_-]?)(\d+)$", clean(value), re.I)
    if m: return int(m.group(1))
    if clean(value).isdigit(): return int(clean(value))
    return None


def canonical_from_metadata(raw: str, aliases: Dict[str, str]):
    parts = [clean(x) for x in re.split(r"[+x]", clean(raw)) if clean(x)]
    mapped = []
    missing = []
    for part in parts:
        if part in aliases:
            mapped.append(aliases[part])
        elif part in aliases.values():
            mapped.append(part)
        else:
            missing.append(part); mapped.append(part)
    return canonical_genotype("+".join(mapped)), missing


def metadata_main():
    args = metadata_parse_args()
    try:
        import pandas as pd
    except ImportError as exc:
        raise SystemExit("pandas/openpyxl are required") from exc

    xls = pd.ExcelFile(args.xlsx)
    required = {"convert", "libs", "Multiome3Prime10X", "2025_LineMeta", "LibUID"}
    missing_sheets = sorted(required - set(xls.sheet_names))
    if missing_sheets:
        raise SystemExit("missing required workbook sheet(s): " + ", ".join(missing_sheets))
    convert = pd.read_excel(args.xlsx, sheet_name="convert", dtype=str)
    convert_extra = pd.read_excel(args.xlsx, sheet_name="convert_extra", dtype=str) if "convert_extra" in xls.sheet_names else None
    libs = pd.read_excel(args.xlsx, sheet_name="libs", dtype=str)
    multi = pd.read_excel(args.xlsx, sheet_name="Multiome3Prime10X", dtype=str)
    line_meta = pd.read_excel(args.xlsx, sheet_name="2025_LineMeta", dtype=str)
    libuid = pd.read_excel(args.xlsx, sheet_name="LibUID", dtype=str)

    aliases, alias_sources, supplemental, spsx = build_aliases(convert, convert_extra)
    panel_ids = load_panel_ids(args.panel_metadata)
    out = Path(args.output_root); out.mkdir(parents=True, exist_ok=True)
    workbook_hash = sha256_file(args.xlsx)
    warnings: List[dict] = []

    alias_rows = []
    for label in sorted(alias_sources):
        canonical = aliases[label]
        eq = sorted(supplemental.get(label, set()))
        alias_rows.append({
            "library_label": label, "canonical_vcf_id": canonical, "spsx": spsx.get(label, ""),
            "canonical_source": alias_sources[label],
            "equivalent_or_reporter_vcf_ids": ",".join([canonical] + eq),
            "noncanonical_equivalent_ids": ",".join(eq),
            "resolution_note": "SUPPLEMENTAL_EQUIVALENT_PRESENT" if eq else "UNAMBIGUOUS",
            "source_workbook_sha256": workbook_hash,
        })
    write_tsv(str(out / "donor_aliases.tsv"), alias_rows, [
        "library_label", "canonical_vcf_id", "spsx", "canonical_source",
        "equivalent_or_reporter_vcf_ids", "noncanonical_equivalent_ids", "resolution_note",
        "source_workbook_sha256",
    ])

    # Global biological-state universe.  2025_LineMeta answers a different
    # question from LibUID/libs: did this exact biological line ever exist
    # anywhere in the project?  Unconstrained demux hypotheses must pass this
    # whitelist before they can become biological single-cell candidates.
    lm_wgs = norm_col(line_meta, "WGS_Key", "wgskey", "genotype")
    if not lm_wgs:
        raise SystemExit("2025_LineMeta requires WGS_Key")
    lm_uid = norm_col(line_meta, "UID", "uid")
    lm_den = norm_col(line_meta, "DenotionR1", "denotion", "line")
    lm_fz = norm_col(line_meta, "FzBatch", "fz_batch")
    lm_grp = norm_col(line_meta, "CorrectedFZGRP", "corrected_fzgrp")
    global_members = defaultdict(list)
    for idx, row in line_meta.iterrows():
        raw = text(row[lm_wgs])
        if not raw:
            continue
        genotype, missing = canonical_from_metadata(raw, aliases)
        if not genotype:
            continue
        global_members[genotype].append({
            "uid": text(row[lm_uid]) if lm_uid else "",
            "line_label": text(row[lm_den]) if lm_den else "",
            "fz_batch": text(row[lm_fz]) if lm_fz else "",
            "corrected_fzgrp": text(row[lm_grp]) if lm_grp else "",
            "source_row": idx + 2,
        })
        for token in missing:
            warnings.append({
                "library": "GLOBAL", "canonical_genotype": genotype,
                "warning": "MISSING_DONOR_ALIAS", "detail": token,
                "source_sheet": "2025_LineMeta", "source_row": idx + 2,
            })

    global_line_rows = []
    donor_to_lines = defaultdict(set)
    for genotype in sorted(global_members):
        members_here = global_members[genotype]
        comps = [x for x in genotype.split("+") if x]
        if len(comps) == 1:
            ploidy = "DIPLOID"
        elif len(comps) == 2 and len(set(comps)) == 1:
            ploidy = "HOMOTYPIC_TETRAPLOID"
        elif len(comps) == 2:
            ploidy = "HETEROTYPIC_TETRAPLOID"
        else:
            ploidy = "OTHER"
        global_line_rows.append({
            "canonical_genotype": genotype,
            "donor_components": ",".join(comps),
            "biological_ploidy_class": ploidy,
            "n_line_meta_rows": len(members_here),
            "uid_candidates": canonical_uid_set(m["uid"] for m in members_here),
            "line_labels": "|".join(sorted({m["line_label"] for m in members_here if m["line_label"]})),
            "fz_batches": canonical_uid_set(m["fz_batch"] for m in members_here),
            "corrected_fzgrps": canonical_uid_set(m["corrected_fzgrp"] for m in members_here),
            "source_sheet": "2025_LineMeta",
            "source_workbook_sha256": workbook_hash,
        })
        for donor in comps:
            donor_to_lines[donor].add(genotype)

    global_donor_rows = [
        {
            "donor_id": donor,
            "n_global_biological_lines": len(lines),
            "global_biological_lines": "|".join(sorted(lines)),
            "source_sheet": "2025_LineMeta",
            "source_workbook_sha256": workbook_hash,
        }
        for donor, lines in sorted(donor_to_lines.items())
    ]
    write_tsv(str(out / "global_biological_lines.tsv"), global_line_rows, GLOBAL_LINE_FIELDS)
    write_tsv(str(out / "global_donors.tsv"), global_donor_rows, GLOBAL_DONOR_FIELDS)

    # Expected genotypes from libs. This establishes library presence even when no UID row exists.
    lib_col = norm_col(libs, "lib", "library", "LibNum")
    line_col = norm_col(libs, "line", "Line", "WGS_Key")
    if not lib_col or not line_col:
        raise SystemExit("libs sheet requires lib and line columns")
    expected = defaultdict(set)
    for idx, row in libs.iterrows():
        n = infer_libnum(text(row[lib_col])); raw = text(row[line_col])
        if n is None or not raw: continue
        genotype, missing = canonical_from_metadata(raw, aliases)
        if genotype: expected[n].add(genotype)
        for token in missing:
            warnings.append({"library": f"lib{n}", "canonical_genotype": genotype, "warning": "MISSING_DONOR_ALIAS", "detail": token, "source_sheet": "libs", "source_row": idx + 2})

    # Normalize LibUID rows; prefer its library+genotype+UID relation as the physical roster.
    c_libnum = norm_col(libuid, "LibNum", "lib", "library")
    c_libid = norm_col(libuid, "LibID", "CatCoreSampleName", "sample")
    c_wgs = norm_col(libuid, "WGS_Key", "wgskey", "genotype")
    c_uid = norm_col(libuid, "UID", "uid")
    if not c_wgs or not c_uid or (not c_libnum and not c_libid):
        raise SystemExit("LibUID requires LibNum/LibID, WGS_Key, and UID columns")
    c_den = norm_col(libuid, "DenotionR1", "denotion", "line")
    c_fz = norm_col(libuid, "FzBatch", "fz_batch")
    c_grp = norm_col(libuid, "CorrectedFZGRP", "corrected_fzgrp")

    members: List[dict] = []
    grouped = defaultdict(list)
    seen_member = set()
    for idx, row in libuid.iterrows():
        n = infer_libnum(text(row[c_libnum])) if c_libnum else infer_libnum(text(row[c_libid]))
        if n is None: continue
        raw_wgs = text(row[c_wgs]); uid = text(row[c_uid])
        genotype, missing = canonical_from_metadata(raw_wgs, aliases)
        scope = f"library:{n}"
        den = text(row[c_den]) if c_den else ""
        fz = text(row[c_fz]) if c_fz else ""
        grp = text(row[c_grp]) if c_grp else ""
        key = (n, genotype, uid, raw_wgs, den, fz, grp)
        if key in seen_member:
            warnings.append({"library": f"lib{n}", "canonical_genotype": genotype, "warning": "DUPLICATE_METADATA_ROW", "detail": f"UID={uid}", "source_sheet": "LibUID", "source_row": idx + 2})
            continue
        seen_member.add(key)
        member = {
            "library": f"lib{n}", "uid_resolution_scope": scope, "canonical_genotype": genotype,
            "uid": uid, "ordered_wgs_key": raw_wgs, "ordered_line_label": den,
            "denotion_r1": den, "fz_batch": fz, "corrected_fzgrp": grp,
            "source_sheet": "LibUID", "source_row": idx + 2,
        }
        members.append(member); grouped[(n, genotype)].append(member); expected[n].add(genotype)
        for token in missing:
            warnings.append({"library": f"lib{n}", "canonical_genotype": genotype, "warning": "MISSING_DONOR_ALIAS", "detail": token, "source_sheet": "LibUID", "source_row": idx + 2})
        for donor in [x for x in genotype.split("+") if x]:
            if donor not in panel_ids:
                warnings.append({"library": f"lib{n}", "canonical_genotype": genotype, "warning": "DONOR_NOT_IN_NUCLEAR_PANEL", "detail": donor, "source_sheet": "LibUID", "source_row": idx + 2})
        if not uid:
            warnings.append({"library": f"lib{n}", "canonical_genotype": genotype, "warning": "MISSING_UID", "detail": raw_wgs, "source_sheet": "LibUID", "source_row": idx + 2})

    write_tsv(str(out / "library_uid_members.tsv"), members, UID_MEMBER_FIELDS)

    expected_rows = []
    resolution_audit = []
    genotype_to_physical = []
    for n in sorted(expected):
        for genotype in sorted(expected[n]):
            rows = grouped.get((n, genotype), [])
            uids = canonical_uid_set(r["uid"] for r in rows)
            uid_list = [x for x in uids.split("|") if x]
            if not rows:
                status = "NO_LIBRARY_METADATA_MATCH"; basis = "library_metadata"; reconciled = ""
            elif not uid_list:
                status = "MISSING_UID_MAPPING"; basis = "library_metadata"; reconciled = ""
            elif len(uid_list) == 1:
                status = "EXACT_LIBRARY_METADATA_MATCH"; basis = "library_metadata"; reconciled = uids
            else:
                status = "MULTIPLE_EXPECTED_UIDS_SAME_GENOTYPE"; basis = "physical_pool_metadata"; reconciled = uids
            detail = []
            for r in sorted(rows, key=lambda x: natural_key(x["uid"])):
                detail.append({"uid": r["uid"], "ordered_wgs_key": r["ordered_wgs_key"], "line_label": r["ordered_line_label"], "fz_batch": r["fz_batch"], "corrected_fzgrp": r["corrected_fzgrp"]})
            comps = [x for x in genotype.split("+") if x]
            ploidy = "DIPLOID" if len(comps) == 1 else ("HOMOTYPIC_TETRAPLOID" if len(set(comps)) == 1 else "HETEROTYPIC_TETRAPLOID")
            aggregate = {
                "library": f"lib{n}", "canonical_genotype": genotype,
                "donor_components": ",".join(comps), "expected_ploidy_class": ploidy,
                "n_metadata_rows": len(rows), "uid_candidate_count": len(uid_list),
                "uid_candidates": uids, "reconciled_uid": reconciled,
                "uid_resolution_status": status, "uid_resolution_basis": basis,
                "uid_resolution_scope": f"library:{n}",
                "uid_member_details": json.dumps(detail, separators=(",", ":")),
                "fz_batches": canonical_uid_set(r["fz_batch"] for r in rows),
                "corrected_fzgrps": canonical_uid_set(r["corrected_fzgrp"] for r in rows),
                "line_labels": "|".join(sorted({r["ordered_line_label"] for r in rows if r["ordered_line_label"]})),
                "source_sheets": "LibUID" if rows else "libs",
            }
            expected_rows.append(aggregate)
            resolution_audit.append(dict(aggregate))
            if rows:
                genotype_to_physical.append({"canonical_genotype": genotype, "library": f"lib{n}", "uid_resolution_scope": aggregate["uid_resolution_scope"], "uid_candidates": uids, "uid_resolution_status": status})
            else:
                warnings.append({"library": f"lib{n}", "canonical_genotype": genotype, "warning": "NO_LIBRARY_METADATA_MATCH", "detail": "genotype appears in libs but has no LibUID member", "source_sheet": "libs", "source_row": ""})

    write_tsv(str(out / "library_expected_genotypes.tsv"), expected_rows, EXPECTED_FIELDS)
    write_tsv(str(out / "library_resolution_audit.tsv"), resolution_audit, EXPECTED_FIELDS)
    write_tsv(str(out / "genotype_to_physical_lines.tsv"), genotype_to_physical, ["canonical_genotype", "library", "uid_resolution_scope", "uid_candidates", "uid_resolution_status"])
    write_tsv(str(out / "metadata_warnings.tsv"), warnings, WARNING_FIELDS)

    # Preparation relationships are descriptive context only; never broaden UID scope here.
    prep_rows = []
    multi_sample = norm_col(multi, "CatCoreSampleName", "sample", "LibID")
    for idx, row in multi.iterrows():
        lib = infer_libnum(text(row[multi_sample])) if multi_sample else None
        if lib is None: continue
        prep_rows.append({
            "library": f"lib{lib}", "relationship_scope": f"library:{lib}",
            "diff_batch": row_value(row, "DiffBatch"), "seed_batch": row_value(row, "SeedBatch"),
            "treatment_patterning": row_value(row, "Treatment_Patterning"),
            "source_sheet": "Multiome3Prime10X", "source_row": idx + 2,
        })
    write_tsv(str(out / "preparation_relationships.tsv"), prep_rows, ["library", "relationship_scope", "diff_batch", "seed_batch", "treatment_patterning", "source_sheet", "source_row"])

    json_dump_atomic(str(out / "metadata_manifest.json"), {
        "schema_version": SCHEMA_VERSION, "workbook": os.path.abspath(args.xlsx),
        "workbook_sha256": workbook_hash, "panel_metadata": os.path.abspath(args.panel_metadata),
        "panel_metadata_sha256": sha256_file(args.panel_metadata), "libraries": sorted(expected),
        "global_biological_line_source": "2025_LineMeta",
        "n_global_biological_lines": len(global_line_rows),
        "n_global_donors": len(global_donor_rows),
        "ambient_rna_evaluated": False,
    })
    print(f"Wrote identity metadata contract to {out}")
    return 0



# -----------------------------------------------------------------------------
# candidates
# -----------------------------------------------------------------------------

import argparse
import glob
import math
import os
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Set, Tuple

from identity_reconciliation_common import (
    SCHEMA_VERSION, biological_state, canonical_genotype, clean, donor_components,
    ffloat, parse_library_spec, read_assignments, read_refined_assignments, read_tsv,
    write_tsv,
)

FIELDS = [
    "library", "barcode", "hypothesis_id", "state_notation", "donor_genotype",
    "donor_components", "biological_ploidy", "droplet_state", "droplet_constituents",
    "candidate_origin", "current_state_notation", "current_donor_genotype",
    "expected_genotype_status", "project_genotype_status", "biological_admissibility",
    "physical_resolution_status", "source_identity",
    "replaced_component", "replacement_component", "preserved_partner",
    "candidate_priority", "candidate_generation_notes", "schema_version",
]
CANDIDATE_EVENT_FIELDS = ["library", "unexpected_component", "source_identity", "replaced_component", "preserved_partner", "n_cells_nominated", "candidate_origin"]

TECHNICAL_MULTIPLET_FIELDS = [
    "library", "barcode", "technical_hypothesis", "technical_class",
    "donor_composition", "current_donor_genotype", "additional_donor",
    "candidate_origin", "global_biological_status", "demux_type",
    "audit_delta_ll_best_vs_expected", "current_droplet_flag", "quad_pattern_score",
    "identification_status", "schema_version",
]

DOUBLET_SUMMARY_FIELDS = [
    "library", "status", "eligible_singlet_cells", "usable_doublet_cells",
    "excluded_real_biological_pair_cells", "excluded_unresolvable_doublet_cells",
    "doublet_rate", "input_assignments", "dd_all", "dd_indv", "detail",
]
DOUBLET_PAIR_FIELDS = [
    "library", "donor_pair", "observed_count", "expected_count", "binomial_p",
    "observed_expected_ratio", "doublet_rate",
]


def args_parser():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--libraries", nargs="+", required=True)
    p.add_argument("--demux-root", required=True, help="mapping_output root containing Tet_2025_Multiome-RNA_N")
    p.add_argument("--audit-root", required=True)
    p.add_argument("--refined-root", default="")
    p.add_argument("--nn-root", required=True)
    p.add_argument("--metadata-root", required=True)
    p.add_argument("--output-root", required=True)
    p.add_argument("--top-k", type=int, default=5)
    p.add_argument("--max-candidates", type=int, default=12)
    p.add_argument("--nn-homotypic-threshold", type=float, default=0.90)
    return p.parse_args()


def demux_prefix(root: str, n: int) -> str:
    return os.path.join(root, f"Tet_2025_Multiome-RNA_{n}", "demux_nomito", f"lib{n}_demuxed")


def audit_path(root: str, n: int, suffix: str) -> str:
    return os.path.join(root, f"lib{n}", f"lib{n}.{suffix}")


def read_optional_tsv(path: str) -> List[dict]:
    try:
        return read_tsv(path) if path and os.path.isfile(path) else []
    except Exception:
        return []


def read_runnerups(path: str) -> Dict[str, List[str]]:
    out = defaultdict(list)
    for row in read_optional_tsv(path):
        bc = clean(row.get("barcode")); ident = canonical_genotype(row.get("identity", ""))
        if bc and ident and ident not in out[bc]:
            out[bc].append(ident)
    return out


def discover_audit_assignment(libdir: str, lib: str, mode: str) -> str:
    patterns = [
        os.path.join(libdir, f"{lib}_audit_{mode}.assignments"),
        os.path.join(libdir, f"{lib}.audit_{mode}.assignments"),
        os.path.join(libdir, f"*audit_{mode}*.assignments"),
    ]
    for pat in patterns:
        hits = sorted(glob.glob(pat))
        if hits:
            return hits[0]
    return ""


def discover_audit_runnerups(libdir: str, lib: str, mode: str) -> str:
    patterns = [
        os.path.join(libdir, f"{lib}_audit_{mode}.runner_ups.gz"),
        os.path.join(libdir, f"{lib}.audit_{mode}.runner_ups.gz"),
        os.path.join(libdir, f"*audit_{mode}*.runner_ups.gz"),
    ]
    for pat in patterns:
        hits = sorted(glob.glob(pat))
        if hits:
            return hits[0]
    return ""


def expected_metadata(metadata_root: str):
    rows = read_tsv(os.path.join(metadata_root, "library_expected_genotypes.tsv"))
    by_lib = defaultdict(dict)
    for row in rows:
        lib = clean(row.get("library")); g = canonical_genotype(row.get("canonical_genotype", ""))
        if lib and g:
            by_lib[lib][g] = row
    return by_lib


def project_biological_universe(metadata_root: str) -> Tuple[Set[str], Set[str]]:
    """Return exact biological lines and donor identities from 2025_LineMeta."""
    line_path = os.path.join(metadata_root, "global_biological_lines.tsv")
    donor_path = os.path.join(metadata_root, "global_donors.tsv")
    if not os.path.isfile(line_path) or not os.path.isfile(donor_path):
        raise RuntimeError(
            "identity metadata is missing global_biological_lines.tsv/global_donors.tsv; "
            "rerun IDENTITY_METADATA with the current code"
        )
    lines = {
        canonical_genotype(r.get("canonical_genotype", ""))
        for r in read_tsv(line_path)
        if canonical_genotype(r.get("canonical_genotype", ""))
    }
    donors = {
        clean(r.get("donor_id"))
        for r in read_tsv(donor_path)
        if clean(r.get("donor_id"))
    }
    return lines, donors


def candidate_biological_status(
    genotype: str, expected_here: Mapping[str, dict], global_lines: Set[str], global_donors: Set[str]
) -> Tuple[str, str]:
    """Classify an unconstrained hypothesis against global and library metadata."""
    genotype = canonical_genotype(genotype)
    comps = donor_components(genotype)
    if not comps:
        return "UNKNOWN", "NOT_BIOLOGICALLY_ADMISSIBLE"
    if len(comps) == 1:
        donor = comps[0]
        if donor not in global_donors:
            return "NOT_REAL_GLOBAL_DONOR", "NOT_BIOLOGICALLY_ADMISSIBLE"
        if genotype in expected_here:
            return "GLOBAL_REAL_DONOR_LIBRARY_EXPECTED", "SINGLET_IDENTITY_CANDIDATE"
        return "GLOBAL_REAL_DONOR_LIBRARY_UNEXPECTED", "SINGLET_IDENTITY_CANDIDATE"
    if len(comps) == 2:
        if genotype not in global_lines:
            return "NOT_REAL_GLOBAL_BIOLOGICAL_LINE", "TECHNICAL_MULTIPLET_ONLY"
        if genotype in expected_here:
            return "GLOBAL_REAL_LINE_LIBRARY_EXPECTED", "BIOLOGICAL_SINGLE_CELL_ALLOWED"
        return "GLOBAL_REAL_LINE_LIBRARY_UNEXPECTED", "BIOLOGICAL_SINGLE_CELL_ALLOWED"
    return "NOT_REAL_GLOBAL_BIOLOGICAL_LINE", "TECHNICAL_MULTIPLET_ONLY"


def nn_rows(path: str) -> Dict[str, dict]:
    out = {}
    for row in read_optional_tsv(path):
        bc = clean(row.get("barcode") or row.get("cell") or row.get("Barcode"))
        if bc: out[bc] = row
    return out


def swap_rows(path: str) -> Dict[str, dict]:
    out = {}
    for row in read_optional_tsv(path):
        bc = clean(row.get("barcode"))
        if bc: out[bc] = row
    return out


def swap_report(path: str) -> dict:
    rows = read_optional_tsv(path)
    return rows[0] if rows else {}


def replace_component(current: str, replacement: str) -> List[Tuple[str, str, str]]:
    """Return (candidate,replaced,preserved) partner-preserving alternatives."""
    comps = donor_components(current)
    repl = clean(replacement)
    if not repl:
        return []
    if len(comps) == 1:
        if comps[0] == repl: return []
        return [(canonical_genotype(repl), comps[0], "")]
    if len(comps) == 2:
        out = []
        for i, old in enumerate(comps):
            if old == repl: continue
            other = comps[1 - i]
            cand = canonical_genotype(f"{repl}+{other}")
            if cand != canonical_genotype(current):
                out.append((cand, old, other))
        return out
    return []


def ploidy_for(genotype: str, homotypic_requested: bool = False) -> str:
    comps = donor_components(genotype)
    if homotypic_requested or (len(comps) == 2 and len(set(comps)) == 1):
        return "TETRAPLOID"
    if len(comps) == 2:
        return "TETRAPLOID"
    return "DIPLOID"


def state_for(genotype: str, homotypic_requested: bool = False) -> str:
    comps = donor_components(genotype)
    if homotypic_requested and len(comps) == 1:
        return f"T[{comps[0]}+{comps[0]}]"
    return biological_state(genotype, "TETRAPLOID" if len(comps) == 2 else "DIPLOID")


def candidates_main():
    args = args_parser(); libs = parse_library_spec(args.libraries)
    meta = expected_metadata(args.metadata_root)
    global_lines, global_donors = project_biological_universe(args.metadata_root)
    out_root = Path(args.output_root); out_root.mkdir(parents=True, exist_ok=True)
    event_counts = defaultdict(int)

    for n in libs:
        lib = f"lib{n}"; prefix = demux_prefix(args.demux_root, n); libdir = os.path.join(args.audit_root, lib)
        current = read_assignments(prefix + ".assignments")
        refined_path = os.path.join(args.refined_root, lib, f"{lib}.refined_assignments") if args.refined_root else ""
        refined = read_refined_assignments(refined_path)
        swaps = swap_rows(audit_path(args.audit_root, n, "swap_scores.tsv"))
        report = swap_report(audit_path(args.audit_root, n, "swap_report.tsv"))
        nn = nn_rows(os.path.join(args.nn_root, f"{lib}.ploidy_calls_nn.tsv"))
        uncon_path = discover_audit_assignment(libdir, lib, "unconstrained")
        con_path = discover_audit_assignment(libdir, lib, "constrained")
        uncon = read_assignments(uncon_path) if uncon_path else {}
        constrained = read_assignments(con_path) if con_path else {}
        uncon_ru = read_runnerups(discover_audit_runnerups(libdir, lib, "unconstrained"))
        con_ru = read_runnerups(discover_audit_runnerups(libdir, lib, "constrained"))
        expected = meta.get(lib, {})
        event_donor = canonical_genotype(clean(report.get("top_unexpected_component_supported")))
        if "+" in event_donor:
            event_donor = ""  # event component must be one donor token

        rows = []
        technical_by_key: Dict[Tuple[str, str, str], dict] = {}
        universe: Set[str] = set()
        for bc, cur in current.items():
            demux_g = canonical_genotype(cur["assignment"])
            ref = refined.get(bc, {})
            refined_g = canonical_genotype(ref.get("refined_assignment", ""))
            # The refined call is the current production interpretation when present.
            # Original demux remains a candidate/audit field, not the state we silently revert to.
            cur_g = refined_g or demux_g
            sw = swaps.get(bc, {})
            u_g = canonical_genotype(uncon.get(bc, {}).get("assignment", "") or sw.get("audit_unconstrained_best", ""))
            c_g = canonical_genotype(constrained.get(bc, {}).get("assignment", "") or sw.get("audit_constrained_best", ""))
            origins: Dict[str, Set[str]] = defaultdict(set)
            notes: Dict[str, List[str]] = defaultdict(list)
            partner_info: Dict[str, Tuple[str, str, str]] = {}

            filtered_global = 0

            def record_technical(g: str, origin: str, technical_class: str = "", additional_donor: str = "") -> None:
                g = canonical_genotype(g)
                comps = donor_components(g)
                if not g or not comps or any(x not in global_donors for x in comps):
                    return
                if not technical_class:
                    technical_class = "UNMADE_DONOR_PAIR" if len(comps) == 2 else "MULTI_DONOR_COMPOSITION"
                if technical_class == "KNOWN_TETRAPLOID_PLUS_EXTRA_DONOR":
                    hypothesis = f"M{{T[{cur_g}]|D[{additional_donor}]}}"
                else:
                    hypothesis = "M{" + "|".join(f"D[{x}]" for x in comps) + "}"
                key = (bc, g, technical_class)
                identification = "TECHNICAL_MULTIPLET_CANDIDATE"
                if (g == cur_g and clean(cur.get("type")).upper() == "D" and len(comps) == 2
                        and g not in global_lines):
                    identification = "GENOTYPE_IDENTIFIABLE_CURRENT_DOUBLET"
                row = technical_by_key.get(key)
                if row is None:
                    row = {
                        "library": lib, "barcode": bc, "technical_hypothesis": hypothesis,
                        "technical_class": technical_class, "donor_composition": g,
                        "current_donor_genotype": cur_g, "additional_donor": additional_donor,
                        "candidate_origin": origin, "global_biological_status": "NOT_REAL_GLOBAL_BIOLOGICAL_LINE",
                        "demux_type": clean(cur.get("type")),
                        "audit_delta_ll_best_vs_expected": clean(sw.get("delta_ll_best_vs_expected")),
                        "current_droplet_flag": clean(ref.get("droplet_flag")),
                        "quad_pattern_score": clean(ref.get("quad_pattern_score")),
                        "identification_status": identification, "schema_version": SCHEMA_VERSION,
                    }
                    technical_by_key[key] = row
                else:
                    have = set(x for x in clean(row.get("candidate_origin")).split(",") if x)
                    have.add(origin)
                    row["candidate_origin"] = ",".join(sorted(have))
                    if identification == "GENOTYPE_IDENTIFIABLE_CURRENT_DOUBLET":
                        row["identification_status"] = identification

            def add(g: str, origin: str, note: str = "") -> bool:
                nonlocal filtered_global
                g = canonical_genotype(g)
                if not g:
                    return False
                if len(donor_components(g)) > 4:
                    return False
                project_status, admissibility = candidate_biological_status(
                    g, expected, global_lines, global_donors
                )
                # Keep globally impossible combinations out of the biological score
                # manifest, but preserve real-donor compositions as a separate
                # technical-multiplet hypothesis.  The frozen H0 remains scoreable
                # even when it is an old demux doublet rather than a real cell line.
                if admissibility not in {
                    "SINGLET_IDENTITY_CANDIDATE", "BIOLOGICAL_SINGLE_CELL_ALLOWED"
                }:
                    record_technical(g, origin)
                    if g != cur_g:
                        filtered_global += 1
                        return False
                origins[g].add(origin)
                if note:
                    notes[g].append(note)
                return True

            add(demux_g, "CURRENT_DEMUX")
            if refined_g: add(refined_g, "CURRENT_REFINED")
            if cur_g not in origins: add(cur_g, "CURRENT_REFINED" if refined_g else "CURRENT_DEMUX")
            if u_g: add(u_g, "AUDIT_UNCONSTRAINED_WINNER")
            if c_g: add(c_g, "AUDIT_CONSTRAINED_WINNER")
            for g in uncon_ru.get(bc, [])[:args.top_k]: add(g, "AUDIT_UNCONSTRAINED_TOP_RUNNER")
            for g in con_ru.get(bc, [])[:args.top_k]: add(g, "AUDIT_CONSTRAINED_TOP_RUNNER")

            # Expected library identities are kept as a bounded context set. Prefer those sharing a current component.
            expected_sorted = sorted(expected)
            shared_expected = [g for g in expected_sorted if set(donor_components(g)) & set(donor_components(cur_g))]
            for g in (shared_expected + expected_sorted)[:args.top_k]: add(g, "EXPECTED_LIBRARY_IDENTITY")

            unexpected = clean(sw.get("unexpected_components"))
            unexpected_tokens = [canonical_genotype(x) for x in unexpected.replace(";", ",").split(",") if clean(x)]
            nominated = unexpected_tokens[:]
            if event_donor and event_donor not in nominated: nominated.append(event_donor)
            for donor in nominated[:args.top_k]:
                for cand, replaced, preserved in replace_component(cur_g, donor):
                    if add(cand, "PARTNER_PRESERVING_REPLACEMENT"):
                        partner_info[cand] = (cur_g, replaced, donor)
                        if preserved:
                            notes[cand].append(f"preserved_partner={preserved}")
                        event_counts[(lib, donor, cur_g, replaced, preserved, "PARTNER_PRESERVING_REPLACEMENT")] += 1
                add(donor, "EVENT_NOMINATED_DONOR")
                # If the current refined state is an established biological
                # tetraploid and audit evidence nominates an additional real donor,
                # retain T[current]+D[donor] as a technical-multiplet hypothesis.
                # We do not invent a three-way nuclear likelihood: the current
                # two-donor scorer cannot score that state.
                cur_comps = donor_components(cur_g)
                if len(cur_comps) == 2 and cur_g in global_lines and donor in global_donors and donor not in cur_comps:
                    composition = canonical_genotype(cur_g + "+" + donor)
                    record_technical(composition, "KNOWN_TETRAPLOID_PLUS_EXTRA_DONOR",
                                     "KNOWN_TETRAPLOID_PLUS_EXTRA_DONOR", donor)

            nnrow = nn.get(bc, {})
            ptet = ffloat(nnrow.get("prob_tetraploid"))
            qc = clean(nnrow.get("qc_pass") or nnrow.get("QC_pass") or nnrow.get("status"))
            nn_good = math.isfinite(ptet) and ptet >= args.nn_homotypic_threshold and qc.lower() not in {"0", "false", "fail", "failed"}
            if nn_good and len(donor_components(cur_g)) == 1:
                homotypic_g = canonical_genotype(f"{cur_g}+{cur_g}")
                # NN can nominate A+A only when that exact homotet line exists in
                # the global 2025_LineMeta universe. It need not be expected in
                # this library because unexpected real lines are valid swap hypotheses.
                if homotypic_g in global_lines:
                    add(homotypic_g, "HOMOTYPIC_FROM_NN", "global_line_eligible_homotypic_state")

            # Sparse multiplet residual nomination: only when the existing refinement/audit already carries a trigger.
            droplet_flag = clean(ref.get("droplet_flag"))
            quad = ffloat(ref.get("quad_pattern_score"))
            trigger = droplet_flag.lower() not in {"", "none", "single", "single_cell", "0", "no"} or (math.isfinite(quad) and quad > 0)
            if trigger:
                component_pool = []
                for g in [cur_g, demux_g, refined_g, u_g, c_g] + nominated:
                    for d in donor_components(g):
                        if d not in component_pool: component_pool.append(d)
                if 3 <= len(component_pool) <= 4:
                    add("+".join(sorted(component_pool)), "MULTIPLET_RESIDUAL_CANDIDATE", "existing_multiplet_trigger")

            # Determine candidate priority and cap while never dropping current/refined.
            priority_order = {
                "CURRENT_DEMUX": 0, "CURRENT_REFINED": 1, "PARTNER_PRESERVING_REPLACEMENT": 2,
                "AUDIT_UNCONSTRAINED_WINNER": 3, "AUDIT_CONSTRAINED_WINNER": 4,
                "EVENT_NOMINATED_DONOR": 5, "HOMOTYPIC_FROM_NN": 6,
                "AUDIT_UNCONSTRAINED_TOP_RUNNER": 7, "EXPECTED_LIBRARY_IDENTITY": 8,
                "MULTIPLET_RESIDUAL_CANDIDATE": 9,
            }
            ordered = sorted(origins, key=lambda g: (min(priority_order.get(o, 99) for o in origins[g]), g))
            # Invalid global hypotheses were already removed by add().  Cap only
            # the biologically admissible set, and never drop the frozen H0.
            kept = [cur_g] if cur_g in origins else []
            for g in ordered:
                if g == cur_g:
                    continue
                if len(kept) >= args.max_candidates:
                    break
                kept.append(g)

            current_state = state_for(cur_g)
            for rank, g in enumerate(kept, 1):
                state = state_for(g, False)
                comps = donor_components(g)
                meta_row = expected.get(g, {})
                expected_status = "EXPECTED" if g in expected else "UNEXPECTED_OR_EVENT_NOMINATED"
                project_status, admissibility = candidate_biological_status(g, expected, global_lines, global_donors)
                source_identity, replaced, replacement = partner_info.get(g, ("", "", ""))
                preserved = ""
                if source_identity and len(donor_components(source_identity)) == 2 and replaced:
                    preserved = next((x for x in donor_components(source_identity) if x != replaced), "")
                if admissibility == "TECHNICAL_MULTIPLET_ONLY":
                    droplet_state = "TECHNICAL_MULTIPLET_CANDIDATE"
                    droplet_constituents = "|".join(f"D[{x}]" for x in comps)
                    biological_ploidy = "UNRESOLVED_MULTIPLET"
                    state = "M{" + "|".join(f"D[{x}]" for x in comps) + "}"
                else:
                    droplet_state = "SINGLE_CELL_CANDIDATE" if len(comps) <= 2 else "TECHNICAL_MULTIPLET_CANDIDATE"
                    droplet_constituents = "|".join(f"D[{x}]" for x in comps) if len(comps) >= 3 else ""
                    biological_ploidy = ploidy_for(g, False)
                rows.append({
                    "library": lib, "barcode": bc, "hypothesis_id": f"{lib}:{bc}:H{rank:02d}",
                    "state_notation": state, "donor_genotype": g, "donor_components": ",".join(comps),
                    "biological_ploidy": biological_ploidy, "droplet_state": droplet_state,
                    "droplet_constituents": droplet_constituents,
                    "candidate_origin": ",".join(sorted(origins[g], key=lambda x: priority_order.get(x, 99))),
                    "current_state_notation": current_state, "current_donor_genotype": cur_g,
                    "expected_genotype_status": expected_status,
                    "project_genotype_status": project_status,
                    "biological_admissibility": admissibility,
                    "physical_resolution_status": clean(meta_row.get("uid_resolution_status")) or "NO_LIBRARY_METADATA_MATCH",
                    "source_identity": source_identity, "replaced_component": replaced,
                    "replacement_component": replacement, "preserved_partner": preserved,
                    "candidate_priority": rank,
                    "candidate_generation_notes": ";".join(sorted(set(notes[g])))
                    + (";CANDIDATE_CAP_APPLIED" if len(origins) > len(kept) else "")
                    + (f";GLOBAL_UNIVERSE_FILTERED={filtered_global}" if filtered_global else ""),
                    "schema_version": SCHEMA_VERSION,
                })
                if len(comps) <= 2: universe.add(g)

        out_path = out_root / f"{lib}.identity_candidates.tsv.gz"
        write_tsv(str(out_path), rows, FIELDS)
        technical_rows = list(technical_by_key.values())
        write_tsv(str(out_root / f"{lib}.technical_multiplet_candidates.tsv.gz"),
                  technical_rows, TECHNICAL_MULTIPLET_FIELDS)
        (out_root / f"{lib}.target_identity_universe.txt").write_text("".join(g + "\n" for g in sorted(universe)))
        print(f"{lib}: {len(rows)} biological candidate rows, {len(technical_rows)} technical-multiplet candidates, {len(universe)} scoreable identities")

    event_rows = []
    for (lib, donor, source, replaced, preserved, origin), count in sorted(event_counts.items()):
        event_rows.append({"library": lib, "unexpected_component": donor, "source_identity": source, "replaced_component": replaced, "preserved_partner": preserved, "n_cells_nominated": count, "candidate_origin": origin})
    write_tsv(str(out_root / "all_libraries.identity_events_candidates.tsv"), event_rows, CANDIDATE_EVENT_FIELDS)
    return 0


# -----------------------------------------------------------------------------
# doublet-context
# -----------------------------------------------------------------------------

def doublet_context_main():
    """Run Doublet Dragon on the diploid-resolvable subset only.

    Real biological two-donor lines from 2025_LineMeta are deliberately excluded
    from Doublet Dragon's D category.  This prevents true tetraploid/fusion lines
    from inflating the estimated technical-doublet rate.  The resulting rate is a
    population context diagnostic, never a per-cell identity decision.
    """
    import subprocess

    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--libraries", nargs="+", required=True)
    p.add_argument("--metadata-root", required=True)
    p.add_argument("--demux-root", required=True)
    p.add_argument("--doublet-dragon", required=True)
    p.add_argument("--output-root", required=True)
    a = p.parse_args()

    libs = parse_library_spec(a.libraries)
    global_lines, global_donors = project_biological_universe(a.metadata_root)
    out_root = Path(a.output_root); out_root.mkdir(parents=True, exist_ok=True)

    for n in libs:
        lib = f"lib{n}"
        prefix = demux_prefix(a.demux_root, n)
        original = read_assignments(prefix + ".assignments")
        singlet_rows = []
        raw_doublets = []
        excluded_real_pair = 0
        excluded_unresolvable = 0

        for bc, row in original.items():
            g = canonical_genotype(row.get("assignment", ""))
            comps = donor_components(g)
            score = abs(ffloat(row.get("score"), 1.0))
            if not math.isfinite(score) or score <= 0:
                score = 1.0
            if len(comps) == 1 and comps[0] in global_donors:
                singlet_rows.append((bc, comps[0], "S", score))
            elif len(comps) == 2 and all(x in global_donors for x in comps):
                if g in global_lines:
                    excluded_real_pair += 1
                else:
                    raw_doublets.append((bc, canonical_genotype(g), "D", score))
            elif len(comps) >= 2:
                excluded_unresolvable += 1

        singlet_ids = {row[1] for row in singlet_rows}
        usable_doublets = []
        for row in raw_doublets:
            comps = donor_components(row[1])
            if len(comps) == 2 and all(x in singlet_ids for x in comps):
                usable_doublets.append(row)
            else:
                excluded_unresolvable += 1

        inp = out_root / f"{lib}.doublet_dragon_input.assignments"
        write_headerless_tsv(str(inp), singlet_rows + usable_doublets)
        dd_prefix = out_root / f"{lib}.doublet_dragon"
        dd_all = Path(str(dd_prefix) + ".dd.all")
        dd_indv = Path(str(dd_prefix) + ".dd.indv")
        for stale in (dd_all, dd_indv):
            if stale.exists():
                stale.unlink()

        summary = {
            "library": lib, "status": "INSUFFICIENT_DIPLOID_RESOLVABLE_INPUT",
            "eligible_singlet_cells": len(singlet_rows),
            "usable_doublet_cells": len(usable_doublets),
            "excluded_real_biological_pair_cells": excluded_real_pair,
            "excluded_unresolvable_doublet_cells": excluded_unresolvable,
            "doublet_rate": math.nan, "input_assignments": str(inp),
            "dd_all": str(dd_all), "dd_indv": str(dd_indv), "detail": "",
        }
        pair_rows = []

        if usable_doublets and len(singlet_ids) >= 1:
            if not (os.path.isfile(a.doublet_dragon) and os.access(a.doublet_dragon, os.X_OK)):
                summary["status"] = "DOUBLET_DRAGON_UNAVAILABLE"
                summary["detail"] = a.doublet_dragon
            else:
                result = subprocess.run(
                    [a.doublet_dragon, str(dd_prefix), str(inp)],
                    stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=False,
                )
                if result.returncode != 0:
                    summary["status"] = "DOUBLET_DRAGON_FAILED"
                    summary["detail"] = (result.stderr.strip() or result.stdout.strip())[-2000:]
                elif not dd_all.is_file() or not dd_indv.is_file():
                    summary["status"] = "DOUBLET_DRAGON_MISSING_OUTPUT"
                    summary["detail"] = "doublet_dragon returned success without .dd.all/.dd.indv"
                else:
                    rate = math.nan
                    with open(dd_all) as fh:
                        for raw in fh:
                            f = raw.rstrip("\n").split("\t")
                            if len(f) >= 3 and f[0] == "all" and f[1] == "doublet_rate":
                                rate = ffloat(f[2])
                                break
                    with open(dd_indv) as fh:
                        for raw in fh:
                            f = raw.rstrip("\n").split("\t")
                            if len(f) < 6 or f[1] != "D":
                                continue
                            pair = canonical_genotype(f[2])
                            obs = fint(f[3]); exp = ffloat(f[4]); pv = ffloat(f[5])
                            ratio = obs / exp if math.isfinite(exp) and exp > 0 else math.nan
                            pair_rows.append({
                                "library": lib, "donor_pair": pair, "observed_count": obs,
                                "expected_count": exp, "binomial_p": pv,
                                "observed_expected_ratio": ratio, "doublet_rate": rate,
                            })
                    summary["status"] = "PASS"
                    summary["doublet_rate"] = rate

        write_tsv(str(out_root / f"{lib}.doublet_dragon_summary.tsv"), [summary], DOUBLET_SUMMARY_FIELDS)
        write_tsv(str(out_root / f"{lib}.doublet_dragon_pairs.tsv"), pair_rows, DOUBLET_PAIR_FIELDS)
        print(
            f"{lib}: Doublet Dragon context {summary['status']}; "
            f"singlets={len(singlet_rows)} usable_D={len(usable_doublets)} "
            f"real_fusion_D_excluded={excluded_real_pair}"
        )
    return 0


# -----------------------------------------------------------------------------
# reconcile
# -----------------------------------------------------------------------------

import argparse
import json
import math
import os
from collections import Counter, defaultdict
from pathlib import Path
from statistics import median
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

from identity_reconciliation_common import (
    DEFAULT_IDENTITY_POLICY, DEFAULT_IDENTITY_POLICY_SHA256, POLICY_VERSION, SCHEMA_VERSION, biological_state, canonical_genotype, clean,
    classify_singlet_library_relationship, donor_components, expected_library_context, ffloat, fint, json_dump_atomic,
    natural_key, pairwise_library_roster,
    parse_library_spec, read_assignments, read_refined_assignments, read_tsv,
    sha256_file, write_headerless_tsv, write_tsv,
)

CELL_FIELDS = [
    "library", "barcode", "original_demux_assignment", "original_demux_type",
    "current_refined_assignment", "current_ploidy", "current_droplet_flag",
    "current_quad_pattern_score", "explicit_multiplet_evidence",
    "proposed_state", "proposed_donor_genotype", "proposed_biological_ploidy",
    "proposed_droplet_state", "proposed_droplet_constituents",
    "proposed_global_biological_status", "proposed_library_expected_status",
    "proposed_biological_admissibility", "singlet_library_relationship",
    "expected_composite_context", "library_exchange_evidence_eligible",
    "occupancy_resolution_status",
    "competing_technical_state",
    "proposed_uid", "proposed_uid_candidate_count", "proposed_uid_candidates",
    "proposed_uid_resolution_status", "proposed_uid_resolution_basis", "proposed_uid_resolution_scope",
    "reconciled_state", "reconciled_donor_genotype", "reconciled_biological_ploidy",
    "reconciled_droplet_state", "original_uid", "reconciled_uid",
    "uid_candidate_count", "uid_candidates", "uid_resolution_status",
    "uid_resolution_basis", "uid_resolution_scope", "uid_member_details",
    "fz_batch_candidates", "corrected_fzgrp_candidates",
    "nuclear_current_ll", "nuclear_alternative_ll", "nuclear_delta_ll", "nuclear_rank",
    "nuclear_informative_depth", "nuclear_dosage_concordance_current",
    "nuclear_dosage_concordance_alternative", "nuclear_residual_current",
    "nuclear_residual_alternative", "nuclear_fold_support_fraction", "nuclear_status",
    "nn_prob_tetraploid", "nn_ploidy_call", "nn_qc_pass",
    "species_current_status", "species_alternative_status", "species_relation",
    "mt_current_status", "mt_alternative_status", "mt_verification_mode", "mt_proposed_component", "mt_current_component",
    "mt_best_identity", "mt_second_identity", "mt_delta_ll", "mt_sites_used", "mt_molecules_used", "mt_haplotype_resolution", "mt_fit_status",
    "atac_current_status", "atac_alternative_status", "atac_delta_ll", "atac_informative_depth", "atac_status",
    "source_identity", "replaced_component", "replacement_component", "preserved_partner",
    "event_id", "event_class", "event_confidence", "alternative_line_event_mass", "final_action", "decision_confidence",
    "reassignment_applied", "decision_reason_codes", "decision_reason", "policy_version", "schema_version",
]
RECONCILE_EVENT_FIELDS = [
    "event_id", "library", "event_scope", "unexpected_component", "n_implicated_cells", "fraction_library_implicated",
    "primary_source_identity", "primary_source_fraction", "source_identity_entropy", "primary_replaced_component",
    "fraction_primary_source_displaced", "n_primary_source_remaining", "fraction_event_from_primary_source",
    "partner_preservation_fraction", "transition_type_summary", "fraction_nn_diploid", "fraction_nn_high_tet",
    "median_nn_prob_tet", "fraction_mt_supports_alternative", "fraction_mt_supports_current", "fraction_mt_unresolved",
    "fraction_atac_supports_alternative", "fraction_atac_supports_current", "fraction_atac_unavailable",
    "fraction_site_fold_replicated", "median_depth_implicated", "median_depth_background",
    "fraction_low_information_implicated", "fraction_low_information_background", "species_consistency",
    "related_library_recurrence", "reciprocal_event_id", "nuclear_confusability_status",
    "reconciled_uid", "uid_candidate_count", "uid_candidates", "uid_resolution_status", "uid_resolution_basis",
    "uid_resolution_scope", "uid_member_details", "singlet_library_relationship", "expected_composite_context",
    "contributes_to_library_exchange_evidence", "event_class", "event_confidence", "event_reason",
]
LIBRARY_EXCHANGE_FIELDS = [
    "library_a", "library_b", "roster_relation", "pair_discriminability", "shared_donors",
    "a_specific_donors", "b_specific_donors", "n_shared_donors", "n_a_specific_donors", "n_b_specific_donors",
    "a_specific_detected_in_a", "a_specific_detected_in_b", "b_specific_detected_in_b", "b_specific_detected_in_a",
    "a_specific_cells_in_a", "a_specific_cells_in_b", "b_specific_cells_in_b", "b_specific_cells_in_a",
    "a_signature_coverage_in_b", "b_signature_coverage_in_a", "a_native_retention_fraction", "b_native_retention_fraction",
    "a_native_displacement_fraction", "b_native_displacement_fraction", "a_to_b_exchange_strength", "b_to_a_exchange_strength",
    "reciprocal_exchange_status", "exchange_interpretation", "exchange_confidence", "supporting_event_ids", "notes",
]
LIBRARY_EXCHANGE_DONOR_FIELDS = [
    "source_library", "target_library", "diagnostic_donor", "diagnostic_for_library",
    "present_in_source_expected_components", "present_in_target_expected_components", "shared_or_specific",
    "source_observed_robust_singlet_cells", "target_observed_robust_singlet_cells",
    "source_singlet_relationship", "target_singlet_relationship", "source_event_id", "target_event_id",
]
META_AMEND_FIELDS = [
    "library", "proposed_donor_genotype", "uid_resolution_scope", "reconciled_uid",
    "uid_candidate_count", "uid_candidates", "uid_resolution_status", "proposed_metadata_event",
    "supporting_libraries", "confidence", "reason",
]


def reconcile_parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--libraries", nargs="+", required=True)
    p.add_argument("--metadata-root", required=True)
    p.add_argument("--candidate-root", required=True)
    p.add_argument("--doublet-context-root", default="")
    p.add_argument("--nuclear-root", required=True)
    p.add_argument("--mt-root", required=True)
    p.add_argument("--atac-root", required=True)
    p.add_argument("--nn-root", required=True)
    p.add_argument("--refined-root", default="")
    p.add_argument("--refined-optional", action="store_true")
    p.add_argument("--audit-root", required=True)
    p.add_argument("--demux-root", required=True)
    p.add_argument("--panel-metadata", default="")
    p.add_argument("--panel-distinguishability", default="")
    p.add_argument("--output-root", required=True)
    p.add_argument("--reports-root", required=True)
    p.add_argument("--evidence-mode", choices=("rna", "rna-atac"), default="rna")
    p.add_argument("--write-reconciled-assignments", action="store_true")
    p.add_argument("--auto-apply", action="store_true", default=None)
    return p.parse_args()


def score_map(path: str) -> Dict[str, Dict[str, dict]]:
    out = defaultdict(dict)
    if not path or not os.path.isfile(path): return out
    for row in read_tsv(path):
        bc = clean(row.get("barcode")); hid = clean(row.get("hypothesis_id"))
        if bc and hid: out[bc][hid] = row
    return out


def mt_donor_score_map(path: str) -> Dict[str, Dict[str, dict]]:
    """Best donor-level MT probe row per barcode/donor.

    The MT scorer emits only single-donor probes.  Some rows correspond to
    ordinary singlet hypotheses and others are synthetic MT_COMPONENT probes
    added solely so a nuclear-proposed real line can be orthogonally checked at
    its constituent-donor level.  Compound mitochondrial genotypes are never
    created or scored here.
    """
    out = defaultdict(dict)
    if not path or not os.path.isfile(path):
        return out
    for row in read_tsv(path):
        bc = clean(row.get("barcode")); donor = canonical_genotype(row.get("donor_genotype", ""))
        if not bc or len(donor_components(donor)) != 1:
            continue
        prior = out[bc].get(donor)
        if prior is None:
            out[bc][donor] = row
            continue
        prior_ll = ffloat(prior.get("mt_log_likelihood"))
        row_ll = ffloat(row.get("mt_log_likelihood"))
        prior_pass = clean(prior.get("mt_fit_status")) in {"PASS", ""} and math.isfinite(prior_ll)
        row_pass = clean(row.get("mt_fit_status")) in {"PASS", ""} and math.isfinite(row_ll)
        if row_pass and not prior_pass:
            out[bc][donor] = row
    return out


def fold_map(path: str) -> Dict[Tuple[str, str], dict]:
    by = defaultdict(list)
    if not path or not os.path.isfile(path): return {}
    for row in read_tsv(path):
        bc = clean(row.get("barcode")); hid = clean(row.get("hypothesis_id"))
        if bc and hid: by[(bc, hid)].append(row)
    out = {}
    for key, rows in by.items():
        vals = [ffloat(r.get("fold_delta_ll")) for r in rows]
        vals = [x for x in vals if math.isfinite(x)]
        supporting = sum(x > 0 for x in vals)
        current = sum(x < 0 for x in vals)
        out[key] = {
            "folds_evaluable": len(vals), "folds_supporting_alternative": supporting,
            "folds_supporting_current": current,
            "fold_support_fraction": supporting / len(vals) if vals else math.nan,
        }
    return out


def candidate_map(path: str):
    by = defaultdict(list)
    for row in read_tsv(path):
        bc = clean(row.get("barcode"))
        if bc: by[bc].append(row)
    return by


def metadata_map(root: str):
    by = defaultdict(dict)
    for row in read_tsv(os.path.join(root, "library_expected_genotypes.tsv")):
        lib = clean(row.get("library")); g = canonical_genotype(row.get("canonical_genotype", ""))
        if lib and g: by[lib][g] = row
    return by


def global_metadata_map(root: str):
    """Global physical-line metadata keyed by canonical genotype.

    This is deliberately separate from the per-library roster.  It is used to
    describe a proposed unexpected identity without pretending that the line
    was expected in the current library.
    """
    out = {}
    for row in read_tsv(os.path.join(root, "global_biological_lines.tsv")):
        g = canonical_genotype(row.get("canonical_genotype", ""))
        if g:
            out[g] = row
    return out


def proposed_uid_metadata(lib: str, genotype: str, library_meta: Mapping[str, Mapping[str, str]], global_meta: Mapping[str, Mapping[str, str]]):
    """Resolve the proposed identity first in-library, then globally.

    Library-exact resolution remains the strongest scope.  A globally real but
    library-unexpected line keeps its global UID candidates and an explicit
    global scope instead of inheriting the UID of the current assignment.
    Real donor components without a standalone diploid line remain resolved as
    donor identities but have no physical-line UID.
    """
    g = canonical_genotype(genotype)
    local = library_meta.get(lib, {}).get(g, {})
    if local:
        return {
            "uid": clean(local.get("reconciled_uid")),
            "uid_candidate_count": fint(local.get("uid_candidate_count")),
            "uid_candidates": clean(local.get("uid_candidates")),
            "uid_resolution_status": clean(local.get("uid_resolution_status")) or "NO_LIBRARY_METADATA_MATCH",
            "uid_resolution_basis": clean(local.get("uid_resolution_basis")) or "library_metadata",
            "uid_resolution_scope": clean(local.get("uid_resolution_scope")) or lib,
        }
    glob = global_meta.get(g, {})
    if glob:
        uids = clean(glob.get("uid_candidates"))
        nuid = len([x for x in uids.split("|") if x])
        return {
            "uid": uids if nuid == 1 else "",
            "uid_candidate_count": nuid,
            "uid_candidates": uids,
            "uid_resolution_status": "EXACT_GLOBAL_METADATA_MATCH" if nuid == 1 else ("MULTIPLE_GLOBAL_UIDS_SAME_GENOTYPE" if nuid > 1 else "GLOBAL_LINE_MISSING_UID"),
            "uid_resolution_basis": "global_2025_LineMeta",
            "uid_resolution_scope": "global:2025_LineMeta",
        }
    return {
        "uid": "",
        "uid_candidate_count": 0,
        "uid_candidates": "",
        "uid_resolution_status": "GLOBAL_DONOR_NO_STANDALONE_LINE_UID",
        "uid_resolution_basis": "global_donor_universe",
        "uid_resolution_scope": "global:donor",
    }


def nn_map(path: str):
    out = {}
    if not os.path.isfile(path): return out
    for row in read_tsv(path):
        bc = clean(row.get("barcode") or row.get("cell"))
        if bc: out[bc] = row
    return out


def call_qc_map(path: str):
    out = {}
    if not os.path.isfile(path): return out
    for row in read_tsv(path):
        bc = clean(row.get("barcode"))
        if bc: out[bc] = row
    return out


def panel_species(path: str) -> Dict[str, str]:
    if not path or not os.path.isfile(path): return {}
    rows = read_tsv(path)
    out = {}
    for row in rows:
        donor = clean(row.get("indiv_id") or row.get("VCF_ID") or row.get("sample"))
        sp = clean(row.get("species") or row.get("Species"))
        if donor and sp: out[donor] = sp
    return out


def panel_distinguishability(path: str) -> Dict[frozenset, dict]:
    out = {}
    if not path or not os.path.isfile(path):
        return out
    for row in read_tsv(path):
        a = clean(row.get("donor1")); b = clean(row.get("donor2"))
        if a and b:
            out[frozenset((a, b))] = row
    return out


def transition_confusability(current: str, alternative: str, lookup: Mapping[frozenset, dict]) -> str:
    if not current or not alternative or current == alternative:
        return "NOT_APPLICABLE"
    cur = Counter(donor_components(current)); alt = Counter(donor_components(alternative))
    common = cur & alt
    cur -= common; alt -= common
    removed = list(cur.elements()); added = list(alt.elements())
    if not removed or not added:
        return "NOT_APPLICABLE"
    statuses = []
    for a in removed:
        for b in added:
            if a == b:
                continue
            row = lookup.get(frozenset((a, b)))
            statuses.append(clean(row.get("confusability_status")) if row else "NOT_EVALUATED")
    if not statuses:
        return "NOT_APPLICABLE"
    if "INDISTINGUISHABLE_ON_PANEL" in statuses:
        return "INDISTINGUISHABLE_ON_PANEL"
    if all(x == "DISTINGUISHABLE_ON_PANEL" for x in statuses):
        return "DISTINGUISHABLE_ON_PANEL"
    if "NO_JOINT_CALLABLE_SITES" in statuses:
        return "NO_JOINT_CALLABLE_SITES"
    return "MIXED_OR_NOT_EVALUATED"


def species_set(genotype: str, panel: Mapping[str, str]) -> set:
    out = set()
    for donor in donor_components(genotype):
        sp = panel.get(donor, "")
        if sp in {"Hy", "Chinobo", "Chinobo-mCherry"} or "chinobo" in sp.lower() or "chinobo" in donor.lower():
            out.update(["B", "C"])
        elif sp:
            out.add(sp)
    return out


def observed_species_set(qc: Mapping[str, str]) -> set:
    raw = clean(qc.get("species_best_identity") or qc.get("observed_species_best") or qc.get("species_audit_best"))
    return set(x for x in raw.replace(";", "+").split("+") if x and x != "NA")


def species_status(genotype: str, qc: Mapping[str, str], panel: Mapping[str, str]) -> str:
    expected = species_set(genotype, panel)
    observed = observed_species_set(qc)
    if not expected or not observed: return "UNINFORMATIVE"
    if expected == observed: return "SUPPORTS"
    if expected & observed: return "PARTIAL"
    return "CONTRADICTS"


def entropy(counts: Counter) -> float:
    total = sum(counts.values())
    if total <= 0: return math.nan
    e = 0.0
    for n in counts.values():
        if n:
            p = n / total; e -= p * math.log(p)
    return e


def frac(items: Sequence[object], pred) -> float:
    return sum(1 for x in items if pred(x)) / len(items) if items else math.nan


def med(vals: Iterable[float]) -> float:
    xs = [x for x in vals if math.isfinite(x)]
    return median(xs) if xs else math.nan


def robust_nuclear_support(delta: float, norm_delta: float, depth: float, fold_fraction: float, folds_evaluable: int, policy: Mapping[str, object]) -> bool:
    return (
        math.isfinite(delta)
        and delta >= float(policy.get("nuclear_strong_delta_ll", 100.0))
        and math.isfinite(norm_delta)
        and norm_delta >= float(policy.get("nuclear_strong_depth_normalized_delta", 0.01))
        and depth >= float(policy.get("nuclear_min_informative_depth", 5000))
        and math.isfinite(fold_fraction)
        and fold_fraction >= float(policy.get("site_fold_support_fraction", 0.80))
        and folds_evaluable >= int(policy.get("site_fold_min_evaluable", 4))
    )


def classify_event(n_cells: int, n_lib: int, source_frac: float, displaced: float, partner_frac: float, recurrence: bool, policy: Mapping[str, object]) -> Tuple[str, str, str]:
    fraction = n_cells / n_lib if n_lib else 0.0
    min_cells = int(policy.get("event_min_cells", 10))
    if n_cells < min_cells:
        return "BELOW_EVENT_MASS_THRESHOLD", "INSUFFICIENT", f"{n_cells} implicated cells is below event_min_cells={min_cells}"
    full = float(policy.get("full_replacement_fraction", 0.90))
    partial = float(policy.get("partial_replacement_fraction", 0.20))
    carry = float(policy.get("foreign_carryover_max_fraction", 0.10))
    if displaced >= full and source_frac >= float(policy.get("event_primary_source_fraction", 0.60)):
        return "LIKELY_COMPLETE_REPLACEMENT", "DECISIVE" if recurrence else "STRONG", "primary source is nearly completely displaced"
    if displaced >= partial and source_frac >= 0.5:
        return "LIKELY_PARTIAL_REPLACEMENT", "STRONG" if recurrence else "SUGGESTIVE", "substantial fraction of one source population is displaced"
    if fraction <= carry and source_frac < 0.5:
        return "LIKELY_FOREIGN_INTACT_CELL_CROSS_CONTAMINATION", "SUGGESTIVE", "coherent low-frequency foreign-cell population without dominant source displacement"
    if partner_frac >= float(policy.get("event_partner_preservation_fraction", 0.60)):
        return "LIKELY_WRONG_FUSION_COMPONENT_IN_METADATA", "STRONG" if recurrence else "SUGGESTIVE", "repeated partner-preserving component replacement"
    return "UNRESOLVED_EVENT", "SUGGESTIVE", "event structure does not meet a calibrated replacement class"


def directional_exchange_strength(n_detected: int, n_diagnostic: int) -> str:
    if n_diagnostic <= 0 or n_detected <= 0:
        return "NONE"
    coverage = n_detected / n_diagnostic
    if n_detected == 1:
        return "WEAK"
    if n_detected >= 3 and coverage >= 0.75:
        return "STRONG"
    return "MODERATE"


def classify_library_exchange(
    roster_relation: str,
    pair_discriminability: str,
    a_to_b_strength: str,
    b_to_a_strength: str,
    a_native_retention: float,
    b_native_retention: float,
) -> Tuple[str, str, str]:
    if roster_relation == "ROSTER_EQUIVALENT_NONDISCRIMINATING":
        return "NONDISCRIMINATING", "ROSTER_EQUIVALENT_NONDISCRIMINATING", "NONDISCRIMINATING"

    a_to_b = a_to_b_strength != "NONE"
    b_to_a = b_to_a_strength != "NONE"
    if not a_to_b and not b_to_a:
        return "NONE", "NO_LIBRARY_EXCHANGE_SIGNAL", "NONE"

    if a_to_b and b_to_a:
        reciprocal = "RECIPROCAL"
        a_displacement = 1.0 - a_native_retention if math.isfinite(a_native_retention) else 0.0
        b_displacement = 1.0 - b_native_retention if math.isfinite(b_native_retention) else 0.0
        strong_reciprocal_signature = (
            pair_discriminability in {"MODERATE", "STRONG"}
            and "STRONG" in {a_to_b_strength, b_to_a_strength}
        )
        if strong_reciprocal_signature and a_displacement >= 0.75 and b_displacement >= 0.75:
            interpretation = "LIKELY_RECIPROCAL_LIBRARY_SWAP"
        else:
            interpretation = "RECIPROCAL_LIBRARY_MIXING"
        if pair_discriminability == "WEAK":
            confidence = "SUGGESTIVE"
        elif interpretation == "LIKELY_RECIPROCAL_LIBRARY_SWAP" and pair_discriminability == "STRONG":
            confidence = "STRONG"
        elif a_to_b_strength == "STRONG" and b_to_a_strength == "STRONG":
            confidence = "STRONG"
        else:
            confidence = "MODERATE"
        return reciprocal, interpretation, confidence

    if a_to_b:
        reciprocal = "ONE_WAY_A_TO_B"
        strength = a_to_b_strength
        target_retention = b_native_retention
    else:
        reciprocal = "ONE_WAY_B_TO_A"
        strength = b_to_a_strength
        target_retention = a_native_retention

    target_displacement = 1.0 - target_retention if math.isfinite(target_retention) else 0.0
    if strength in {"MODERATE", "STRONG"} and target_displacement >= 0.50:
        interpretation = "LIKELY_PARTIAL_LIBRARY_REPLACEMENT"
    elif strength in {"MODERATE", "STRONG"} and math.isfinite(target_retention) and target_retention >= 0.75:
        interpretation = "LIKELY_CROSS_CONTAMINATION"
    else:
        interpretation = "ONE_WAY_FOREIGN_SIGNATURE"

    if pair_discriminability == "WEAK" or strength == "WEAK":
        confidence = "SUGGESTIVE"
    elif strength == "STRONG":
        confidence = "STRONG"
    else:
        confidence = "MODERATE"
    return reciprocal, interpretation, confidence


def decision_for_cell(cell: dict, policy: Mapping[str, object], auto_apply: bool):
    """Conservative adjudication of one already-nominated alternative.

    Audit/unconstrained results are hypotheses, never sufficient authority for a
    final identity.  Biological admissibility, robust nuclear support and the
    channel-specific rules below are required before a change can even be called
    DECISIVE.  Auto-application is separately controlled and defaults off.
    """
    reasons = []
    current = cell["current_donor"]
    alternative = cell["best_alt"]
    if not alternative or alternative == current:
        return "KEEP", "INSUFFICIENT", False, ["NO_SUPPORTED_ALTERNATIVE"]

    alt_comps = donor_components(alternative)
    cur_comps = donor_components(current)
    cand = cell.get("candidate", {})
    admissibility = clean(cand.get("biological_admissibility"))
    project_status = clean(cand.get("project_genotype_status"))

    # A two-donor combination that was never made anywhere in the project is not
    # a biological tetraploid line.  It may be evidence for a technical D+D
    # doublet, but this first-version reconciler does not pretend to resolve that
    # occupancy question from a pair likelihood alone.
    if admissibility == "TECHNICAL_MULTIPLET_ONLY" or project_status == "NOT_REAL_GLOBAL_BIOLOGICAL_LINE":
        return (
            "KEEP_CURRENT_CONFLICTED", "CONFLICTED", False,
            ["UNMADE_DONOR_PAIR_NOT_BIOLOGICAL_SINGLE_CELL", "TECHNICAL_MULTIPLET_REQUIRES_OCCUPANCY_ADJUDICATION"],
        )
    if len(alt_comps) >= 3:
        return "KEEP_CURRENT_CONFLICTED", "CONFLICTED", False, ["MULTI_DONOR_SIGNAL_REQUIRES_TECHNICAL_MULTIPLET_MODEL"]

    # Homotypic A -> A+A is a ploidy decision, not a nuclear donor substitution.
    # A high NN P(tet) remains decisive only when the local occupancy context does
    # not leave an unresolved same-donor technical-doublet explanation.
    if len(cur_comps) == 1 and len(alt_comps) == 2 and len(set(alt_comps)) == 1 and alt_comps[0] == cur_comps[0]:
        ptet = cell["nn_prob"]
        eligible = (
            admissibility == "BIOLOGICAL_SINGLE_CELL_ALLOWED"
            and project_status in {"GLOBAL_REAL_LINE_LIBRARY_EXPECTED", "GLOBAL_REAL_LINE_LIBRARY_UNEXPECTED"}
        )
        high_nn = (
            eligible
            and math.isfinite(ptet)
            and ptet >= float(policy.get("nn_homotypic_prob_tet", 0.90))
            and cell["nn_qc"]
        )
        if high_nn:
            reasons = ["NN_HIGH_TETRAPLOID_SUPPORT", "GLOBAL_REAL_HOMOTYPIC_LINE"]
            if cell.get("explicit_multiplet_evidence"):
                return (
                    "REVIEW_HOMOTET_OCCUPANCY", "CONFLICTED", False,
                    reasons + ["EXPLICIT_MULTIPLET_EVIDENCE_PRESENT", "HOMOTET_CELLULAR_ORIGIN_UNRESOLVED"],
                )
            if (
                bool(policy.get("unexpected_homotet_occupancy_guard", True))
                and cell.get("occupancy_resolution_status") == "HOMOTET_VS_SAME_DONOR_DOUBLET_UNRESOLVED"
            ):
                return (
                    "REVIEW_HOMOTET_OCCUPANCY", "CONFLICTED", False,
                    reasons + [
                        "HOMOTET_LIBRARY_UNEXPECTED",
                        "SAME_DONOR_DIPLOID_POPULATION_PRESENT",
                        "HOMOTET_VS_SAME_DONOR_DOUBLET_UNRESOLVED",
                    ],
                )
            return "RECLASSIFY_PLOIDY", "DECISIVE", auto_apply, reasons
        return "KEEP", "STRONG_NOT_AUTOAPPLIED", False, ["HOMOTYPIC_PLOIDY_NOT_SUFFICIENTLY_SUPPORTED"]

    nuclear_delta = cell["nuclear_delta"]
    nuclear_norm = cell.get("nuclear_norm_delta", math.nan)
    nuclear_depth = cell["nuclear_depth"]
    folds = cell["fold_fraction"]
    strong = robust_nuclear_support(
        nuclear_delta, nuclear_norm, nuclear_depth, folds, cell["folds_evaluable"], policy
    )
    suggestive = math.isfinite(nuclear_delta) and nuclear_delta >= float(policy.get("nuclear_suggestive_delta_ll", 20.0))
    fold_ok = (
        math.isfinite(folds)
        and folds >= float(policy.get("site_fold_support_fraction", 0.80))
        and cell["folds_evaluable"] >= int(policy.get("site_fold_min_evaluable", 4))
    )
    species_ok = cell["species_alt"] != "CONTRADICTS"
    atac_alt = cell["atac_status"] == "ATAC_SUPPORTS_ALTERNATIVE"
    atac_current = cell["atac_status"] == "ATAC_SUPPORTS_CURRENT"
    event_ok = cell["event_confidence"] in {"DECISIVE", "STRONG"}

    if strong:
        reasons += ["NUCLEAR_ALTERNATIVE_STRONGLY_SUPPORTED", "NUCLEAR_EFFECT_SIZE_SUFFICIENT"]
    if fold_ok:
        reasons.append("NUCLEAR_SUPPORT_REPLICATED_ACROSS_SITE_FOLDS")
    if cell.get("preserved_partner"):
        reasons.append("PARTNER_PRESERVED_IN_REPLACEMENT")
    if not species_ok:
        reasons.append("SPECIES_CONTRADICTS_ALTERNATIVE")
    if atac_alt:
        reasons.append("ATAC_SUPPORTS_ALTERNATIVE")
    if atac_current:
        reasons.append("ATAC_SUPPORTS_CURRENT")

    # MT may decisively verify a proposed singlet.  For an exact real two-donor
    # line with zero donor overlap versus the current call, donor-level MT is
    # strictly corroborative: it can confirm that a proposed-unique component is
    # present, but it cannot nominate, construct, promote, or veto the fusion.
    alt_is_singlet = len(alt_comps) == 1
    mt_line_mode = cell.get("mt_verification_mode") == "DISJOINT_REAL_LINE_COMPONENT"
    mt_alt = (alt_is_singlet or mt_line_mode) and cell["mt_alt"] == "SUPPORTS_ALTERNATIVE"
    mt_current = (alt_is_singlet or mt_line_mode) and cell["mt_alt"] == "SUPPORTS_CURRENT"
    mt_contra = alt_is_singlet and cell["mt_alt"] == "CONTRADICTS_ALTERNATIVE"
    if alt_is_singlet and mt_alt:
        reasons.append("MITOCHONDRIA_VERIFY_PROPOSED_SINGLET")
    if alt_is_singlet and mt_current:
        reasons.append("MITOCHONDRIA_SUPPORT_CURRENT_SINGLET")
    if alt_is_singlet and mt_contra:
        reasons.append("MITOCHONDRIA_CONTRADICT_PROPOSED_SINGLET")
    # Do not label donor-level MT as corroborating a multi-donor line until the
    # proposed line has independently cleared the robust nuclear + site-fold gate.
    # Raw MT fields remain available on weak proposals, but MT cannot make a weak
    # fusion proposal look supported simply because one discriminating haplotype
    # happens to be observed.

    if not species_ok or atac_current or (alt_is_singlet and (mt_current or mt_contra)):
        return "KEEP_CURRENT_CONFLICTED", "CONFLICTED", False, reasons or ["ORTHOGONAL_EVIDENCE_CONFLICT"]

    if not strong or not fold_ok:
        if suggestive:
            return "KEEP", "STRONG_NOT_AUTOAPPLIED", False, reasons + ["NUCLEAR_SUPPORT_NOT_STRONG_AND_REPLICATED_ENOUGH"]
        if nuclear_depth < float(policy.get("nuclear_low_information_depth", 500)):
            return "UNRESOLVED_INSUFFICIENT_EVIDENCE", "INSUFFICIENT", False, reasons + ["LOW_INFORMATION_NUCLEAR_DEPTH"]
        return "KEEP", "SUGGESTIVE", False, reasons + ["UNRESOLVED_WEAK_OR_DIFFUSE_NUCLEAR_SIGNAL"]

    # Disjoint-line MT is corroboration only. Attach its interpretation after
    # strong nuclear support and replicated site-fold support already exist.
    if mt_line_mode and mt_alt:
        reasons.append("MITOCHONDRIA_SUPPORT_DISJOINT_PROPOSED_LINE_COMPONENT")
    if mt_line_mode and mt_current:
        reasons.append("MITOCHONDRIA_SUPPORT_DISJOINT_CURRENT_LINE_COMPONENT_NONDECISIONAL")

    if alt_is_singlet:
        # A library-unexpected donor can be a very strong cell-level observation,
        # but production reassignment requires substantially more exact-line mass
        # than is needed merely to report a coherent event.  Independent MT/ATAC
        # verification does not bypass either population-level gate.
        singlet_relationship = clean(cell.get("singlet_library_relationship"))
        unexpected_here = singlet_relationship == "UNEXPECTED_SINGLET"
        line_mass = int(cell.get("alternative_line_event_mass", 0))
        min_event_cells = int(policy.get("event_min_cells", 10))
        min_autoapply_cells = int(policy.get("unexpected_line_autoapply_min_cells", 100))
        if unexpected_here and line_mass < min_event_cells:
            return "REVIEW_UNEXPECTED_IDENTITY", "STRONG_NOT_AUTOAPPLIED", False, reasons + ["ALTERNATIVE_LINE_BELOW_EVENT_MASS_THRESHOLD"]
        if unexpected_here and line_mass < min_autoapply_cells:
            return "REVIEW_UNEXPECTED_IDENTITY", "STRONG_NOT_AUTOAPPLIED", False, reasons + ["ALTERNATIVE_LINE_CLEARS_EVENT_MASS_THRESHOLD", "ALTERNATIVE_LINE_BELOW_AUTOAPPLY_MASS_THRESHOLD"]
        if mt_alt:
            return "REASSIGN_GENOTYPE", "DECISIVE", auto_apply, reasons + (["ALTERNATIVE_LINE_CLEARS_EVENT_MASS_THRESHOLD", "ALTERNATIVE_LINE_CLEARS_AUTOAPPLY_MASS_THRESHOLD"] if unexpected_here else [])
        if atac_alt:
            return "REASSIGN_GENOTYPE", "DECISIVE", auto_apply, reasons + ["SINGLET_VERIFIED_BY_INDEPENDENT_ATAC"] + (["ALTERNATIVE_LINE_CLEARS_EVENT_MASS_THRESHOLD", "ALTERNATIVE_LINE_CLEARS_AUTOAPPLY_MASS_THRESHOLD"] if unexpected_here else [])
        return "KEEP", "STRONG_NOT_AUTOAPPLIED", False, reasons + ["PROPOSED_SINGLET_NOT_INDEPENDENTLY_VERIFIED"]

    # A heterotypic single-cell reassignment is only admissible if the genotype
    # exists as a physical biological line somewhere in the project.  For a
    # library-unexpected line, >=event_min_cells is sufficient to report a
    # coherent event, but >=unexpected_line_autoapply_min_cells exact-line robust
    # cells are required before production assignments are rewritten.  Homotets
    # are handled above by the independent NN ploidy route and are exempt.
    if len(alt_comps) == 2:
        if admissibility != "BIOLOGICAL_SINGLE_CELL_ALLOWED":
            return "KEEP_CURRENT_CONFLICTED", "CONFLICTED", False, reasons + ["ALTERNATIVE_NOT_BIOLOGICALLY_ADMISSIBLE"]
        unexpected_here = clean(cand.get("expected_genotype_status")) != "EXPECTED"
        line_mass = int(cell.get("alternative_line_event_mass", 0))
        min_event_cells = int(policy.get("event_min_cells", 10))
        min_autoapply_cells = int(policy.get("unexpected_line_autoapply_min_cells", 100))

        # Strong A+B nuclear evidence establishes donor composition, not whether
        # those donors occupied one intact biological cell.  Explicit existing
        # multiplet structure is therefore a veto, and when both A and B are
        # independently present as diploid populations the competing
        # M{D[A]|D[B]} technical-doublet state remains observationally plausible.
        explicit_multiplet = bool(cell.get("explicit_multiplet_evidence"))
        local_pair_ambiguity = (
            bool(policy.get("heterotypic_pair_occupancy_guard", True))
            and len(set(alt_comps)) == 2
            and cell.get("occupancy_resolution_status") == "DONOR_PAIR_CELLULAR_ORIGIN_UNRESOLVED"
        )
        if explicit_multiplet or local_pair_ambiguity:
            occ_reasons = reasons + ["NUCLEAR_DONOR_COMPOSITION_AB_SUPPORTED"]
            if explicit_multiplet:
                occ_reasons += ["EXPLICIT_MULTIPLET_EVIDENCE_PRESENT", "TECHNICAL_MULTIPLET_OCCUPANCY_SUPPORTED"]
            if local_pair_ambiguity:
                occ_reasons += [
                    "BOTH_DONORS_PRESENT_AS_DIPLOID_POPULATIONS",
                    "TECHNICAL_DOUBLET_PLAUSIBLE",
                    "DONOR_PAIR_CELLULAR_ORIGIN_UNRESOLVED",
                ]
            if unexpected_here:
                occ_reasons.append(
                    "ALTERNATIVE_LINE_CLEARS_EVENT_MASS_THRESHOLD"
                    if line_mass >= min_event_cells
                    else "ALTERNATIVE_LINE_BELOW_EVENT_MASS_THRESHOLD"
                )
                if line_mass >= min_autoapply_cells:
                    occ_reasons.append("ALTERNATIVE_LINE_CLEARS_AUTOAPPLY_MASS_THRESHOLD")
                else:
                    occ_reasons.append("ALTERNATIVE_LINE_BELOW_AUTOAPPLY_MASS_THRESHOLD")
            return "REVIEW_CELLULAR_ORIGIN", "CONFLICTED", False, occ_reasons

        if unexpected_here and line_mass < min_event_cells:
            return "REVIEW_UNEXPECTED_IDENTITY", "STRONG_NOT_AUTOAPPLIED", False, reasons + ["ALTERNATIVE_LINE_BELOW_EVENT_MASS_THRESHOLD"]
        if unexpected_here and line_mass < min_autoapply_cells:
            return "REVIEW_UNEXPECTED_IDENTITY", "STRONG_NOT_AUTOAPPLIED", False, reasons + ["ALTERNATIVE_LINE_CLEARS_EVENT_MASS_THRESHOLD", "ALTERNATIVE_LINE_BELOW_AUTOAPPLY_MASS_THRESHOLD"]
        if event_ok:
            return "REASSIGN_GENOTYPE", "DECISIVE", auto_apply, reasons + ["COHERENT_REPEATED_NUCLEAR_EVENT"] + (["ALTERNATIVE_LINE_CLEARS_EVENT_MASS_THRESHOLD", "ALTERNATIVE_LINE_CLEARS_AUTOAPPLY_MASS_THRESHOLD"] if unexpected_here else [])
        if atac_alt:
            return "REASSIGN_GENOTYPE", "DECISIVE", auto_apply, reasons + ["INDEPENDENT_ATAC_SUPPORT_FOR_KNOWN_BIOLOGICAL_LINE"] + (["ALTERNATIVE_LINE_CLEARS_EVENT_MASS_THRESHOLD", "ALTERNATIVE_LINE_CLEARS_AUTOAPPLY_MASS_THRESHOLD"] if unexpected_here else [])
        return "KEEP", "STRONG_NOT_AUTOAPPLIED", False, reasons + ["KNOWN_BIOLOGICAL_PAIR_BUT_EVENT_SUPPORT_NOT_DECISIVE"]

    return "KEEP", "INSUFFICIENT", False, reasons + ["NO_DECISIVE_ROUTE"]


def reconcile_main():
    args = reconcile_parse_args(); libs = parse_library_spec(args.libraries)
    policy = dict(DEFAULT_IDENTITY_POLICY); policy_hash = DEFAULT_IDENTITY_POLICY_SHA256
    auto_apply = bool(policy.get("auto_apply_decisive", True)) if args.auto_apply is None else args.auto_apply
    meta = metadata_map(args.metadata_root); panel = panel_species(args.panel_metadata)
    explicit_expected_singlets, expected_donor_components, expected_composites_by_component = expected_library_context(meta)
    global_meta = global_metadata_map(args.metadata_root)
    global_lines, global_donors = project_biological_universe(args.metadata_root)
    panel_confusability = panel_distinguishability(args.panel_distinguishability)
    out_root = Path(args.output_root); rep_root = Path(args.reports_root)
    out_root.mkdir(parents=True, exist_ok=True); rep_root.mkdir(parents=True, exist_ok=True)

    # Frozen per-cell evidence assembled first; event context is derived without iteration.
    frozen_by_lib: Dict[str, List[dict]] = {}
    for n in libs:
        lib = f"lib{n}"; demux = os.path.join(args.demux_root, f"Tet_2025_Multiome-RNA_{n}", "demux_nomito", f"{lib}_demuxed")
        original = read_assignments(demux + ".assignments")
        refined = read_refined_assignments(os.path.join(args.refined_root, lib, f"{lib}.refined_assignments")) if args.refined_root else {}
        candidates = candidate_map(os.path.join(args.candidate_root, f"{lib}.identity_candidates.tsv.gz"))
        nuclear = score_map(os.path.join(args.nuclear_root, f"{lib}.identity_hypothesis_scores.tsv.gz"))
        folds = fold_map(os.path.join(args.nuclear_root, f"{lib}.identity_site_fold_scores.tsv.gz"))
        mt_path = os.path.join(args.mt_root, f"{lib}.mt_identity_scores.tsv.gz")
        mt = score_map(mt_path)
        mt_donors = mt_donor_score_map(mt_path)
        atac_path = os.path.join(args.atac_root, f"{lib}.atac_identity_scores.tsv.gz")
        atac = score_map(atac_path) if args.evidence_mode == "rna-atac" else defaultdict(dict)
        nn = nn_map(os.path.join(args.nn_root, f"{lib}.ploidy_calls_nn.tsv"))
        qc = call_qc_map(os.path.join(args.audit_root, lib, f"{lib}.call_qc.tsv.gz"))
        dd_summary_rows = read_optional_tsv(os.path.join(args.doublet_context_root, f"{lib}.doublet_dragon_summary.tsv")) if args.doublet_context_root else []
        dd_summary = dd_summary_rows[0] if dd_summary_rows else {}
        dd_pairs = {
            canonical_genotype(r.get("donor_pair", "")): r
            for r in (read_optional_tsv(os.path.join(args.doublet_context_root, f"{lib}.doublet_dragon_pairs.tsv")) if args.doublet_context_root else [])
            if canonical_genotype(r.get("donor_pair", ""))
        }
        cells = []
        for bc, orig in original.items():
            demux_g = canonical_genotype(orig["assignment"])
            ref = refined.get(bc, {})
            ref_g = canonical_genotype(ref.get("refined_assignment", ""))
            cur_g = ref_g or demux_g
            cand_rows = candidates.get(bc, [])
            nuc_rows = nuclear.get(bc, {})
            current_score = None
            best_alt = None
            best_alt_score = None
            best_alt_row = None
            for cand in cand_rows:
                hid = clean(cand.get("hypothesis_id")); g = canonical_genotype(cand.get("donor_genotype", ""))
                sc = nuc_rows.get(hid)
                if not sc:
                    continue
                ll = ffloat(sc.get("log_likelihood") or sc.get("nuclear_log_likelihood"))
                if g == cur_g and (current_score is None or ll > ffloat(current_score.get("log_likelihood"))):
                    current_score = sc
                    continue
                if g == cur_g:
                    continue
                # Defense in depth: the candidate generator already removes globally
                # impossible biological alternatives, but reconciliation refuses them
                # again so stale manifests cannot recreate the previous failure mode.
                admissibility = clean(cand.get("biological_admissibility"))
                if admissibility not in {"SINGLET_IDENTITY_CANDIDATE", "BIOLOGICAL_SINGLE_CELL_ALLOWED"}:
                    continue
                if math.isfinite(ll) and (best_alt_score is None or ll > best_alt_score):
                    best_alt_score = ll; best_alt = g; best_alt_row = (cand, sc)
            cur_ll = ffloat(current_score.get("log_likelihood")) if current_score else math.nan
            alt_ll = best_alt_score if best_alt_score is not None else math.nan
            delta = alt_ll - cur_ll if math.isfinite(alt_ll) and math.isfinite(cur_ll) else (ffloat(best_alt_row[1].get("delta_ll_vs_current")) if best_alt_row else math.nan)
            cand = best_alt_row[0] if best_alt_row else {}
            sc = best_alt_row[1] if best_alt_row else {}
            fold = folds.get((bc, clean(cand.get("hypothesis_id"))), {})
            nnrow = nn.get(bc, {}); ptet = ffloat(nnrow.get("prob_tetraploid")); qctxt = clean(nnrow.get("qc_pass") or nnrow.get("status"))
            nn_qc = qctxt.lower() not in {"0", "false", "fail", "failed", "qc_fail"}

            # If the best donor-substitution hypothesis does not clear the strong,
            # replicated nuclear bar, fall back to the same-donor homotet question.
            # The NN is used only in its strong direction (diploid -> tetraploid),
            # and only when the exact A+A line is globally real.
            alt_robust = robust_nuclear_support(
                delta,
                ffloat(sc.get("nuclear_depth_normalized_delta") or sc.get("depth_normalized_delta")),
                ffloat(sc.get("nuclear_informative_depth") or sc.get("informative_depth"), 0.0),
                ffloat(fold.get("fold_support_fraction")),
                fint(fold.get("folds_evaluable")),
                policy,
            ) if best_alt_row else False
            cur_comps_for_fallback = donor_components(cur_g)
            if (
                len(cur_comps_for_fallback) == 1
                and nn_qc
                and math.isfinite(ptet)
                and ptet >= float(policy.get("nn_homotypic_prob_tet", 0.90))
                and not alt_robust
            ):
                homotet_g = canonical_genotype(cur_comps_for_fallback[0] + "+" + cur_comps_for_fallback[0])
                if homotet_g in global_lines:
                    homotet_choice = None
                    for hcand in cand_rows:
                        if canonical_genotype(hcand.get("donor_genotype", "")) != homotet_g:
                            continue
                        hhid = clean(hcand.get("hypothesis_id"))
                        hsc = nuc_rows.get(hhid)
                        if hsc:
                            homotet_choice = (hcand, hsc)
                            break
                    if homotet_choice is not None:
                        cand, sc = homotet_choice
                        best_alt = homotet_g
                        best_alt_row = homotet_choice
                        alt_ll = ffloat(sc.get("log_likelihood") or sc.get("nuclear_log_likelihood"))
                        delta = alt_ll - cur_ll if math.isfinite(alt_ll) and math.isfinite(cur_ll) else ffloat(sc.get("delta_ll_vs_current"))
                        fold = folds.get((bc, clean(cand.get("hypothesis_id"))), {})

            qcrow = qc.get(bc, {})
            sp_cur = species_status(cur_g, qcrow, panel); sp_alt = species_status(best_alt or cur_g, qcrow, panel)

            hid = clean(cand.get("hypothesis_id"))
            alt_comps = donor_components(best_alt or "")
            cur_comps_mt = donor_components(cur_g)
            mt_verification_mode = "NOT_APPLICABLE"
            mt_proposed_component = ""
            mt_current_component = ""

            # MT has two deliberately narrow roles:
            #   1) verify a nuclear-proposed singlet donor; or
            #   2) after nuclear evidence independently proposes an exact real
            #      two-donor biological line that shares ZERO donors with the
            #      current call, ask whether the best resolved mitochondrial
            #      donor belongs uniquely to the proposed or current identity.
            # MT never constructs, nominates, or scores a compound genotype.
            if len(alt_comps) == 1 and hid:
                mt_verification_mode = "SINGLET_IDENTITY"
                mt_proposed_component = alt_comps[0]
                mt_current_component = cur_comps_mt[0] if len(cur_comps_mt) == 1 else ""
                mtrow = mt.get(bc, {}).get(hid, {})
                mt_delta = ffloat(mtrow.get("mt_delta_vs_best_other_singlet") or mtrow.get("mt_delta_hypothesis_vs_current") or mtrow.get("mt_delta_ll_vs_current"))
                mt_resolution = clean(mtrow.get("mt_haplotype_resolution"))
                mt_best = canonical_genotype(mtrow.get("mt_best_identity", ""))
                mt_fit = clean(mtrow.get("mt_fit_status"))
                if mt_resolution == "MT_HAPLOTYPE_UNRESOLVED":
                    mt_status = "NONDISCRIMINATING"
                elif mt_fit not in {"PASS", ""}:
                    mt_status = "MISSING" if not mtrow else "NONDISCRIMINATING"
                elif math.isfinite(mt_delta) and mt_delta >= float(policy.get("mt_support_delta_ll", 4.0)) and mt_best == best_alt:
                    mt_status = "SUPPORTS_ALTERNATIVE"
                elif math.isfinite(mt_delta) and mt_delta >= float(policy.get("mt_support_delta_ll", 4.0)) and mt_best:
                    mt_status = "SUPPORTS_CURRENT" if len(cur_comps_mt) == 1 and mt_best == cur_g else "CONTRADICTS_ALTERNATIVE"
                elif mtrow:
                    mt_status = "NONDISCRIMINATING"
                else:
                    mt_status = "MISSING"
            elif (
                len(alt_comps) == 2
                and canonical_genotype(best_alt or "") in global_lines
                and cur_comps_mt
                and set(alt_comps).isdisjoint(set(cur_comps_mt))
            ):
                mt_verification_mode = "DISJOINT_REAL_LINE_COMPONENT"
                donor_rows = mt_donors.get(bc, {})

                def mt_usable_component(donor):
                    row = donor_rows.get(donor, {})
                    ll = ffloat(row.get("mt_log_likelihood"))
                    fit = clean(row.get("mt_fit_status"))
                    if row and fit in {"PASS", ""} and math.isfinite(ll):
                        return (donor, row, ll)
                    return None

                proposed_scored = [x for x in (mt_usable_component(d) for d in sorted(set(alt_comps))) if x]
                current_scored = [x for x in (mt_usable_component(d) for d in sorted(set(cur_comps_mt))) if x]
                all_scored = [x for x in (mt_usable_component(d) for d in sorted(donor_rows)) if x]
                if proposed_scored and current_scored and all_scored:
                    best_proposed = max(proposed_scored, key=lambda x: x[2])
                    best_current_mt = max(current_scored, key=lambda x: x[2])
                    best_global_mt = max(all_scored, key=lambda x: x[2])
                    mt_proposed_component = best_proposed[0]
                    mt_current_component = best_current_mt[0]
                    mt_delta = best_proposed[2] - best_current_mt[2]
                    mtrow = best_global_mt[1]
                    mt_resolution = clean(mtrow.get("mt_haplotype_resolution"))
                    threshold = float(policy.get("mt_support_delta_ll", 4.0))
                    if mt_resolution == "MT_HAPLOTYPE_UNRESOLVED":
                        mt_status = "NONDISCRIMINATING"
                    elif best_global_mt[0] in set(alt_comps) and mt_delta >= threshold:
                        mt_status = "SUPPORTS_ALTERNATIVE"
                        mt_proposed_component = best_global_mt[0]
                    elif best_global_mt[0] in set(cur_comps_mt) and mt_delta <= -threshold:
                        mt_status = "SUPPORTS_CURRENT"
                        mt_current_component = best_global_mt[0]
                    else:
                        # A third donor may be the best MT haplotype, or the
                        # proposed/current groups may be too close.  Neither is
                        # allowed to steer a multi-donor identity decision.
                        mt_status = "NONDISCRIMINATING"
                else:
                    mtrow = {}
                    mt_delta = math.nan
                    mt_status = "MISSING"
            else:
                mtrow = {}
                mt_delta = math.nan
                mt_status = "NOT_APPLICABLE_NON_SINGLET"

            if args.evidence_mode == "rna":
                atacrow = {}; atac_status = "ATAC_NOT_REQUESTED"; atac_delta = math.nan
            else:
                atacrow = atac.get(bc, {}).get(hid, {}) if hid else {}
                atac_delta = ffloat(atacrow.get("atac_delta_ll_vs_current") or atacrow.get("delta_ll_vs_current"))
                raw_status = clean(atacrow.get("atac_status"))
                if raw_status: atac_status = raw_status
                elif math.isfinite(atac_delta) and atac_delta >= float(policy.get("atac_support_delta_ll", 8.0)): atac_status = "ATAC_SUPPORTS_ALTERNATIVE"
                elif math.isfinite(atac_delta) and atac_delta <= -float(policy.get("atac_support_delta_ll", 8.0)): atac_status = "ATAC_SUPPORTS_CURRENT"
                elif atacrow: atac_status = "ATAC_NONDISCRIMINATING"
                else: atac_status = "ATAC_UNAVAILABLE"
            comps = donor_components(best_alt or cur_g)
            multiplet_trigger = len(comps) >= 3 or clean(cand.get("candidate_origin")).find("MULTIPLET_RESIDUAL_CANDIDATE") >= 0
            droplet_flag = clean(ref.get("droplet_flag"))
            quad = ffloat(ref.get("quad_pattern_score"))
            # Hard occupancy veto requires explicit droplet-level evidence.
            # quad_pattern_score remains useful genotype-composition evidence
            # (and may nominate residual multiplet candidates upstream), but a
            # positive value alone cannot establish that multiple cells
            # occupied the droplet.
            explicit_multiplet_evidence = (
                droplet_flag.lower() not in {"", "none", "single", "single_cell", "0", "no"}
            )
            singlet_relationship, singlet_composite_context = classify_singlet_library_relationship(
                best_alt or cur_g,
                explicit_expected_singlets.get(lib, set()),
                expected_donor_components.get(lib, set()),
                expected_composites_by_component.get(lib, {}),
                globally_valid=(len(donor_components(best_alt or cur_g)) != 1 or (donor_components(best_alt or cur_g)[0] in global_donors)),
            )
            cells.append({
                "library": lib, "barcode": bc, "orig": orig, "ref": ref, "original_demux": demux_g, "current_donor": cur_g,
                "best_alt": best_alt, "candidate": cand, "score": sc, "current_score": current_score or {},
                "nuclear_current_ll": cur_ll, "nuclear_alt_ll": alt_ll, "nuclear_delta": delta,
                "nuclear_depth": ffloat(sc.get("nuclear_informative_depth") or sc.get("informative_depth"), 0.0),
                "nuclear_norm_delta": ffloat(sc.get("nuclear_depth_normalized_delta") or sc.get("depth_normalized_delta")),
                "fold_fraction": ffloat(fold.get("fold_support_fraction")), "folds_evaluable": fint(fold.get("folds_evaluable")),
                "nn_prob": ptet, "nn_qc": nn_qc, "nnrow": nnrow, "qc": qcrow,
                "species_cur": sp_cur, "species_alt": sp_alt, "mtrow": mtrow, "mt_alt": mt_status,
                "mt_delta": mt_delta, "mt_verification_mode": mt_verification_mode,
                "mt_proposed_component": mt_proposed_component, "mt_current_component": mt_current_component,
                "atacrow": atacrow, "atac_status": atac_status, "atac_delta": atac_delta,
                "multiplet_trigger": multiplet_trigger,
                "explicit_multiplet_evidence": explicit_multiplet_evidence,
                "current_identified_technical_multiplet": (
                    clean(orig.get("type")).upper() == "D" and len(donor_components(cur_g)) == 2 and cur_g not in global_lines
                ),
                "dd_summary": dd_summary, "dd_pair": dd_pairs.get(cur_g, {}),
                "source_identity": clean(cand.get("source_identity")), "replaced_component": clean(cand.get("replaced_component")),
                "replacement_component": clean(cand.get("replacement_component")), "preserved_partner": clean(cand.get("preserved_partner")),
                "nuclear_confusability_status": transition_confusability(cur_g, best_alt or cur_g, panel_confusability),
                "singlet_library_relationship": singlet_relationship,
                "expected_composite_context": singlet_composite_context,
            })
        frozen_by_lib[lib] = cells

    # Pass A: event context from frozen evidence.
    event_context: Dict[Tuple[str, str], dict] = {}
    line_event_context: Dict[Tuple[str, str], dict] = {}
    component_libraries = defaultdict(set)
    line_libraries = defaultdict(set)
    min_event_cells = int(policy.get("event_min_cells", 10))
    # Recurrence is itself event-level evidence: a library only counts toward
    # cross-library recurrence after it clears the same absolute mass threshold.
    robust_component_counts = defaultdict(Counter)
    robust_line_counts = defaultdict(Counter)
    for lib, cells in frozen_by_lib.items():
        for c in cells:
            replacement = c["replacement_component"] or (donor_components(c["best_alt"])[0] if c["best_alt"] else "")
            robust = (
                math.isfinite(c["nuclear_delta"]) and c["nuclear_delta"] >= float(policy.get("nuclear_strong_delta_ll", 100.0))
                and math.isfinite(c["nuclear_norm_delta"]) and c["nuclear_norm_delta"] >= float(policy.get("nuclear_strong_depth_normalized_delta", 0.01))
                and c["nuclear_depth"] >= float(policy.get("nuclear_min_informative_depth", 5000))
                and math.isfinite(c["fold_fraction"]) and c["fold_fraction"] >= float(policy.get("site_fold_support_fraction", 0.80))
                and c["folds_evaluable"] >= int(policy.get("site_fold_min_evaluable", 4))
            )
            if robust and c["best_alt"]:
                robust_line_counts[lib][canonical_genotype(c["best_alt"])] += 1
            if (
                replacement
                and robust
                and c.get("singlet_library_relationship") != "EXPECTED_COMPOSITE_COMPONENT_SINGLET"
            ):
                robust_component_counts[lib][replacement] += 1
    for lib, counts in robust_component_counts.items():
        for component, count in counts.items():
            if count >= min_event_cells:
                component_libraries[component].add(lib)
    for lib, counts in robust_line_counts.items():
        for genotype, count in counts.items():
            if count >= min_event_cells:
                comps = donor_components(genotype)
                if len(comps) == 1:
                    relationship, _ = classify_singlet_library_relationship(
                        genotype,
                        explicit_expected_singlets.get(lib, set()),
                        expected_donor_components.get(lib, set()),
                        expected_composites_by_component.get(lib, {}),
                        globally_valid=(comps[0] in global_donors),
                    )
                    if relationship != "UNEXPECTED_SINGLET":
                        continue
                line_libraries[genotype].add(lib)

    # Library-local diploid populations define the competing technical-doublet
    # occupancy model.  Expected singlet lines are physical pool members by
    # metadata; robust unexpected singlet events are added once they clear the
    # same event_min_cells threshold used for coherent event reporting.
    local_diploid_donors = defaultdict(set)
    for lib in frozen_by_lib:
        for genotype in meta.get(lib, {}):
            comps = donor_components(genotype)
            if len(comps) == 1:
                local_diploid_donors[lib].add(comps[0])
        current_component_singlet_counts = Counter()
        for c in frozen_by_lib[lib]:
            current = canonical_genotype(c.get("current_donor") or "")
            comps = donor_components(current)
            if len(comps) != 1:
                continue
            current_ploidy = clean(c.get("ref", {}).get("ploidy")).upper()
            if c.get("explicit_multiplet_evidence") or current_ploidy in {"T", "TET", "TETRAPLOID", "HOMOTYPIC"}:
                continue
            relationship, _ = classify_singlet_library_relationship(
                current,
                explicit_expected_singlets.get(lib, set()),
                expected_donor_components.get(lib, set()),
                expected_composites_by_component.get(lib, {}),
                globally_valid=(comps[0] in global_donors),
            )
            if relationship == "EXPECTED_COMPOSITE_COMPONENT_SINGLET":
                current_component_singlet_counts[comps[0]] += 1
        for donor, count in current_component_singlet_counts.items():
            if count >= min_event_cells:
                local_diploid_donors[lib].add(donor)
        for genotype, count in robust_line_counts.get(lib, {}).items():
            comps = donor_components(genotype)
            if len(comps) == 1 and count >= min_event_cells:
                local_diploid_donors[lib].add(comps[0])

    for lib, cells in frozen_by_lib.items():
        local_diploids = local_diploid_donors[lib]
        for c in cells:
            alt = canonical_genotype(c.get("best_alt") or "")
            comps = donor_components(alt)
            expected_here = clean(c.get("candidate", {}).get("expected_genotype_status")) == "EXPECTED"
            if len(comps) == 2 and len(set(comps)) == 2 and all(d in local_diploids for d in comps):
                c["occupancy_resolution_status"] = "DONOR_PAIR_CELLULAR_ORIGIN_UNRESOLVED"
                c["competing_technical_state"] = "M{" + "|".join(f"D[{d}]" for d in comps) + "}"
            elif (
                len(comps) == 2
                and len(set(comps)) == 1
                and not expected_here
                and comps[0] in local_diploids
            ):
                c["occupancy_resolution_status"] = "HOMOTET_VS_SAME_DONOR_DOUBLET_UNRESOLVED"
                c["competing_technical_state"] = f"M{{D[{comps[0]}]|D[{comps[0]}]}}"
            elif len(comps) == 2:
                c["occupancy_resolution_status"] = "NO_LOCAL_DIPLOID_DOUBLET_CONFOUND"
                c["competing_technical_state"] = ""
            else:
                c["occupancy_resolution_status"] = "NOT_APPLICABLE"
                c["competing_technical_state"] = ""

    for lib, cells in frozen_by_lib.items():
        nlib = len(cells); by_component = defaultdict(list); source_counts = defaultdict(Counter)
        background_depth = [c["nuclear_depth"] for c in cells]
        for c in cells:
            replacement = c["replacement_component"] or (next((d for d in donor_components(c["best_alt"] or "") if d not in donor_components(c["current_donor"])), ""))
            robust = (
                math.isfinite(c["nuclear_delta"]) and c["nuclear_delta"] >= float(policy.get("nuclear_strong_delta_ll", 100.0))
                and math.isfinite(c["nuclear_norm_delta"]) and c["nuclear_norm_delta"] >= float(policy.get("nuclear_strong_depth_normalized_delta", 0.01))
                and c["nuclear_depth"] >= float(policy.get("nuclear_min_informative_depth", 5000))
                and math.isfinite(c["fold_fraction"]) and c["fold_fraction"] >= float(policy.get("site_fold_support_fraction", 0.80))
                and c["folds_evaluable"] >= int(policy.get("site_fold_min_evaluable", 4))
            )
            if not replacement or not c["best_alt"] or not robust:
                continue
            if c.get("singlet_library_relationship") == "EXPECTED_COMPOSITE_COMPONENT_SINGLET":
                continue
            by_component[replacement].append(c); source_counts[replacement][c["source_identity"] or c["current_donor"]] += 1
        for component, implicated in by_component.items():
            scount = source_counts[component]; primary, pn = scount.most_common(1)[0]; source_frac = pn / len(implicated)
            source_total = sum(1 for c in cells if c["current_donor"] == primary)
            source_implicated = sum(1 for c in implicated if (c["source_identity"] or c["current_donor"]) == primary)
            displaced = source_implicated / source_total if source_total else 0.0
            preserved = frac(implicated, lambda c: bool(c["preserved_partner"]))
            recurrence = len(component_libraries[component]) > 1
            event_class, conf, reason = classify_event(len(implicated), nlib, source_frac, displaced, preserved, recurrence, policy)
            eid = f"{lib}:{component}"
            event_context[(lib, component)] = {
                "event_id": eid, "event_scope": "COMPONENT", "event_class": event_class, "event_confidence": conf, "event_reason": reason,
                "unexpected_component": component, "n_implicated": len(implicated), "nlib": nlib,
                "primary_source": primary, "source_frac": source_frac, "source_entropy": entropy(scount),
                "source_total": source_total, "source_implicated": source_implicated, "displaced": displaced,
                "preserved": preserved, "recurrence": recurrence, "implicated": implicated, "background_depth": background_depth,
            }

    # Exact-identity events take precedence over component-level interpretations.
    # This prevents a coherent real line such as A+B from being reported merely
    # as an A component replacement when the same intact A+B identity dominates.
    for lib, cells in frozen_by_lib.items():
        nlib = len(cells)
        background_depth = [c["nuclear_depth"] for c in cells]
        by_line = defaultdict(list)
        for c in cells:
            g = canonical_genotype(c.get("best_alt") or "")
            if not g or clean(c.get("candidate", {}).get("expected_genotype_status")) == "EXPECTED":
                continue
            if clean(c.get("candidate", {}).get("biological_admissibility")) not in {"SINGLET_IDENTITY_CANDIDATE", "BIOLOGICAL_SINGLE_CELL_ALLOWED"}:
                continue
            if robust_nuclear_support(
                c["nuclear_delta"], c["nuclear_norm_delta"], c["nuclear_depth"],
                c["fold_fraction"], c["folds_evaluable"], policy
            ):
                by_line[g].append(c)
        for genotype, implicated in by_line.items():
            n = len(implicated)
            comps = donor_components(genotype)
            singlet_relationship, singlet_composite_context = classify_singlet_library_relationship(
                genotype,
                explicit_expected_singlets.get(lib, set()),
                expected_donor_components.get(lib, set()),
                expected_composites_by_component.get(lib, {}),
                globally_valid=(len(comps) != 1 or comps[0] in global_donors),
            )
            scount = Counter(c["source_identity"] or c["current_donor"] for c in implicated)
            primary, pn = scount.most_common(1)[0]
            source_frac = pn / n if n else 0.0
            source_total = sum(1 for c in cells if c["current_donor"] == primary)
            source_implicated = sum(1 for c in implicated if (c["source_identity"] or c["current_donor"]) == primary)
            displaced = source_implicated / source_total if source_total else 0.0
            preserved = frac(implicated, lambda c: bool(c["preserved_partner"]))
            recurrence = len(line_libraries[genotype]) > 1
            mt_line_support = sum(
                c.get("mt_verification_mode") == "DISJOINT_REAL_LINE_COMPONENT" and c.get("mt_alt") == "SUPPORTS_ALTERNATIVE"
                for c in implicated
            )
            mt_line_current = sum(
                c.get("mt_verification_mode") == "DISJOINT_REAL_LINE_COMPONENT" and c.get("mt_alt") == "SUPPORTS_CURRENT"
                for c in implicated
            )
            if n < min_event_cells:
                event_class = "BELOW_EVENT_MASS_THRESHOLD"
                conf = "INSUFFICIENT"
                if singlet_relationship == "EXPECTED_COMPOSITE_COMPONENT_SINGLET":
                    reason = (
                        f"{n} robust cells support residual component singlet {genotype}, explained by "
                        f"library-expected composite genotype(s) {singlet_composite_context}; below event_min_cells={min_event_cells}"
                    )
                else:
                    reason = f"{n} robust cells support exact unexpected identity {genotype}; below event_min_cells={min_event_cells}"
            elif len(comps) == 2:
                occupancy_unresolved = any(
                    c.get("occupancy_resolution_status") == "DONOR_PAIR_CELLULAR_ORIGIN_UNRESOLVED"
                    for c in implicated
                )
                if occupancy_unresolved:
                    event_class = "DONOR_PAIR_CELLULAR_ORIGIN_UNRESOLVED"
                    conf = "STRONG" if recurrence else "SUGGESTIVE"
                    reason = (
                        f"{n} robust cells support donor composition {genotype}, but both donor singlets "
                        "are independently present in this library; intact T[A+B] and technical "
                        "M{D[A]|D[B]} occupancy are not distinguished by the available evidence"
                    )
                else:
                    event_class = "LIKELY_UNEXPECTED_INTACT_BIOLOGICAL_LINE"
                    conf = "DECISIVE" if recurrence else "STRONG"
                    reason = f"{n} robust cells converge on the same globally real, library-unexpected biological line without a local two-singlet occupancy confound"
                if mt_line_support or mt_line_current:
                    reason += f"; donor-level MT supports a proposed-unique component in {mt_line_support}/{n} cells and a current-unique component in {mt_line_current}/{n} cells"
            elif singlet_relationship == "EXPECTED_COMPOSITE_COMPONENT_SINGLET":
                event_class = "EXPECTED_COMPOSITE_COMPONENT_SINGLET_POPULATION"
                conf = "STRONG"
                reason = (
                    f"{n} robust cells converge on donor {genotype}, which is not explicitly expected as a singlet "
                    f"but is a component of library-expected composite genotype(s): {singlet_composite_context}"
                )
            else:
                event_class = "LIKELY_UNEXPECTED_SINGLET_POPULATION"
                conf = "DECISIVE" if recurrence else "STRONG"
                reason = f"{n} robust cells converge on the same globally real donor identity absent from the library expected donor-component universe"
            line_event_context[(lib, genotype)] = {
                "event_id": f"{lib}:identity:{genotype}", "event_scope": "EXACT_IDENTITY",
                "event_class": event_class, "event_confidence": conf, "event_reason": reason,
                "unexpected_component": genotype, "n_implicated": n, "nlib": nlib,
                "primary_source": primary, "source_frac": source_frac, "source_entropy": entropy(scount),
                "source_total": source_total, "source_implicated": source_implicated, "displaced": displaced,
                "preserved": preserved, "recurrence": recurrence, "implicated": implicated, "background_depth": background_depth,
                "singlet_library_relationship": singlet_relationship,
                "expected_composite_context": singlet_composite_context,
                "contributes_to_library_exchange_evidence": (
                    singlet_relationship == "UNEXPECTED_SINGLET" and n >= min_event_cells
                ),
            }

    all_rows = []; event_rows = []; amendments = []
    for lib, cells in frozen_by_lib.items():
        lib_rows = []
        for c in cells:
            replacement = c["replacement_component"] or (next((d for d in donor_components(c["best_alt"] or "") if d not in donor_components(c["current_donor"])), ""))
            exact_g = canonical_genotype(c.get("best_alt") or "")
            ev = line_event_context.get((lib, exact_g))
            if ev is None:
                ev = event_context.get((lib, replacement), {"event_id": "", "event_scope": "NONE", "event_class": "NO_METADATA_PROBLEM", "event_confidence": "INSUFFICIENT"})
            c["event_class"] = ev.get("event_class", "NO_METADATA_PROBLEM"); c["event_confidence"] = ev.get("event_confidence", "INSUFFICIENT")
            c["alternative_line_event_mass"] = robust_line_counts[lib].get(exact_g, 0)
            action, confidence, applied, reasons = decision_for_cell(c, policy, auto_apply)
            proposed_g = c["best_alt"] or c["current_donor"]
            cand_ploidy = clean(c["candidate"].get("biological_ploidy"))
            cand_state = clean(c["candidate"].get("state_notation"))
            cand_droplet = clean(c["candidate"].get("droplet_state"))
            current_ploidy = clean(c["ref"].get("ploidy"))
            current_identified_multiplet = bool(c.get("current_identified_technical_multiplet"))
            if current_identified_multiplet:
                comps_now = donor_components(c["current_donor"])
                current_ploidy = "UNRESOLVED_MULTIPLET"
                current_state = "M{" + "|".join(f"D[{x}]" for x in comps_now) + "}"
            else:
                if not current_ploidy:
                    current_ploidy = "TETRAPLOID" if len(donor_components(c["current_donor"])) == 2 else "DIPLOID"
                current_state = biological_state(c["current_donor"], current_ploidy)
            reconciled_g = proposed_g if applied and action in {"REASSIGN_GENOTYPE", "RECLASSIFY_PLOIDY"} else c["current_donor"]
            if applied and action == "REASSIGN_GENOTYPE":
                reconciled_ploidy = cand_ploidy if cand_ploidy in {"DIPLOID", "TETRAPLOID"} else current_ploidy
                reconciled_state = cand_state or biological_state(reconciled_g, reconciled_ploidy)
            elif action == "RECLASSIFY_PLOIDY" and applied:
                reconciled_state = biological_state(reconciled_g, "TETRAPLOID")
                reconciled_ploidy = "TETRAPLOID"
            else:
                reconciled_state = current_state
                reconciled_ploidy = current_ploidy
            droplet_state = "TECHNICAL_MULTIPLET" if current_identified_multiplet and not applied else "SINGLE_CELL"
            proposed_uid_meta = proposed_uid_metadata(lib, proposed_g, meta, global_meta)
            reconciled_uid_meta = proposed_uid_metadata(lib, reconciled_g, meta, global_meta)
            uid_status = clean(reconciled_uid_meta.get("uid_resolution_status")) or "NO_LIBRARY_METADATA_MATCH"
            uid = clean(reconciled_uid_meta.get("uid"))
            proposed_uid_status = clean(proposed_uid_meta.get("uid_resolution_status"))
            if proposed_uid_status == "EXACT_LIBRARY_METADATA_MATCH": reasons.append("UID_EXACT_LIBRARY_METADATA_MATCH")
            elif proposed_uid_status == "MULTIPLE_EXPECTED_UIDS_SAME_GENOTYPE": reasons.append("UID_SET_FROM_COLOADED_EQUIVALENT_LINES")
            elif proposed_uid_status == "EXACT_GLOBAL_METADATA_MATCH": reasons.append("UID_EXACT_GLOBAL_METADATA_MATCH")
            elif proposed_uid_status == "MULTIPLE_GLOBAL_UIDS_SAME_GENOTYPE": reasons.append("UID_GLOBAL_SET_FOR_UNEXPECTED_LINE")
            else: reasons.append("UID_NOT_FOUND_FOR_PROPOSED_IDENTITY")
            singlet_relationship, singlet_composite_context = classify_singlet_library_relationship(
                proposed_g,
                explicit_expected_singlets.get(lib, set()),
                expected_donor_components.get(lib, set()),
                expected_composites_by_component.get(lib, {}),
                globally_valid=(len(donor_components(proposed_g)) != 1 or donor_components(proposed_g)[0] in global_donors),
            )
            if singlet_relationship == "EXPECTED_COMPOSITE_COMPONENT_SINGLET":
                reasons.append("SINGLET_IS_COMPONENT_OF_LIBRARY_EXPECTED_COMPOSITE")
            exchange_eligible = bool(
                ev.get("contributes_to_library_exchange_evidence")
                and exact_g == proposed_g
                and robust_nuclear_support(
                    c["nuclear_delta"], c["nuclear_norm_delta"], c["nuclear_depth"],
                    c["fold_fraction"], c["folds_evaluable"], policy,
                )
            )
            meta_row = meta.get(lib, {}).get(reconciled_g, {})
            global_row = global_meta.get(reconciled_g, {})
            score = c["score"]; current_score = c["current_score"]; q = c["qc"]; mtrow = c["mtrow"]; atacrow = c["atacrow"]
            row = {
                "library": lib, "barcode": c["barcode"], "original_demux_assignment": c["original_demux"],
                "original_demux_type": c["orig"].get("type", ""), "current_refined_assignment": canonical_genotype(c["ref"].get("refined_assignment", "")),
                "current_ploidy": clean(c["ref"].get("ploidy")), "current_droplet_flag": clean(c["ref"].get("droplet_flag")),
                "current_quad_pattern_score": clean(c["ref"].get("quad_pattern_score")),
                "explicit_multiplet_evidence": bool(c.get("explicit_multiplet_evidence")),
                "proposed_state": cand_state or biological_state(proposed_g, cand_ploidy),
                "proposed_donor_genotype": proposed_g, "proposed_biological_ploidy": cand_ploidy or "UNRESOLVED",
                "proposed_droplet_state": "TECHNICAL_MULTIPLET_CANDIDATE" if cand_droplet == "TECHNICAL_MULTIPLET_CANDIDATE" else "SINGLE_CELL",
                "proposed_droplet_constituents": clean(c["candidate"].get("droplet_constituents")),
                "proposed_global_biological_status": clean(c["candidate"].get("project_genotype_status")) or "CURRENT_STATE",
                "proposed_library_expected_status": clean(c["candidate"].get("expected_genotype_status")) or "CURRENT_STATE",
                "proposed_biological_admissibility": clean(c["candidate"].get("biological_admissibility")) or "CURRENT_STATE",
                "singlet_library_relationship": singlet_relationship,
                "expected_composite_context": singlet_composite_context,
                "library_exchange_evidence_eligible": exchange_eligible,
                "occupancy_resolution_status": c.get("occupancy_resolution_status", "NOT_APPLICABLE"),
                "competing_technical_state": c.get("competing_technical_state", ""),
                "proposed_uid": proposed_uid_meta.get("uid", ""),
                "proposed_uid_candidate_count": proposed_uid_meta.get("uid_candidate_count", 0),
                "proposed_uid_candidates": proposed_uid_meta.get("uid_candidates", ""),
                "proposed_uid_resolution_status": proposed_uid_meta.get("uid_resolution_status", ""),
                "proposed_uid_resolution_basis": proposed_uid_meta.get("uid_resolution_basis", ""),
                "proposed_uid_resolution_scope": proposed_uid_meta.get("uid_resolution_scope", ""),
                "reconciled_state": reconciled_state, "reconciled_donor_genotype": reconciled_g,
                "reconciled_biological_ploidy": reconciled_ploidy, "reconciled_droplet_state": droplet_state,
                "original_uid": "", "reconciled_uid": uid, "uid_candidate_count": reconciled_uid_meta.get("uid_candidate_count", 0),
                "uid_candidates": reconciled_uid_meta.get("uid_candidates", ""), "uid_resolution_status": uid_status,
                "uid_resolution_basis": reconciled_uid_meta.get("uid_resolution_basis", ""), "uid_resolution_scope": reconciled_uid_meta.get("uid_resolution_scope", ""),
                "uid_member_details": meta_row.get("uid_member_details", ""),
                "fz_batch_candidates": meta_row.get("fz_batches", global_row.get("fz_batches", "")),
                "corrected_fzgrp_candidates": meta_row.get("corrected_fzgrps", global_row.get("corrected_fzgrps", "")),
                "nuclear_current_ll": c["nuclear_current_ll"], "nuclear_alternative_ll": c["nuclear_alt_ll"], "nuclear_delta_ll": c["nuclear_delta"],
                "nuclear_rank": score.get("rank_within_candidates", ""), "nuclear_informative_depth": c["nuclear_depth"],
                "nuclear_dosage_concordance_current": current_score.get("dosage_concordance", ""), "nuclear_dosage_concordance_alternative": score.get("dosage_concordance", ""),
                "nuclear_residual_current": current_score.get("residual_mismatch", ""), "nuclear_residual_alternative": score.get("residual_mismatch", ""),
                "nuclear_fold_support_fraction": c["fold_fraction"], "nuclear_status": score.get("score_status", "MISSING") if score else "MISSING",
                "nn_prob_tetraploid": c["nn_prob"], "nn_ploidy_call": c["nnrow"].get("ploidy_call", c["nnrow"].get("predicted_ploidy", "")), "nn_qc_pass": c["nn_qc"],
                "species_current_status": c["species_cur"], "species_alternative_status": c["species_alt"], "species_relation": q.get("species_relation", ""),
                "mt_current_status": "SUPPORTS_CURRENT" if c["mt_alt"] == "SUPPORTS_CURRENT" else "", "mt_alternative_status": c["mt_alt"],
                "mt_verification_mode": c.get("mt_verification_mode", "NOT_APPLICABLE"),
                "mt_proposed_component": c.get("mt_proposed_component", ""), "mt_current_component": c.get("mt_current_component", ""),
                "mt_best_identity": mtrow.get("mt_best_identity", ""), "mt_second_identity": mtrow.get("mt_second_identity", ""), "mt_delta_ll": c["mt_delta"],
                "mt_sites_used": mtrow.get("mt_sites_used", ""), "mt_molecules_used": mtrow.get("mt_molecules_used", ""),
                "mt_haplotype_resolution": mtrow.get("mt_haplotype_resolution", ""), "mt_fit_status": mtrow.get("mt_fit_status", ""),
                "atac_current_status": "ATAC_SUPPORTS_CURRENT" if c["atac_status"] == "ATAC_SUPPORTS_CURRENT" else "",
                "atac_alternative_status": c["atac_status"], "atac_delta_ll": c["atac_delta"],
                "atac_informative_depth": atacrow.get("atac_informative_depth", atacrow.get("informative_depth", "")), "atac_status": c["atac_status"],
                "source_identity": c["source_identity"], "replaced_component": c["replaced_component"], "replacement_component": replacement,
                "preserved_partner": c["preserved_partner"], "event_id": ev.get("event_id", ""), "event_class": ev.get("event_class", "NO_METADATA_PROBLEM"),
                "event_confidence": ev.get("event_confidence", "INSUFFICIENT"), "alternative_line_event_mass": c.get("alternative_line_event_mass", 0), "final_action": action, "decision_confidence": confidence,
                "reassignment_applied": applied, "decision_reason_codes": ",".join(sorted(set(reasons))),
                "decision_reason": "; ".join(r.replace("_", " ").lower() for r in reasons), "policy_version": clean(policy.get("policy_version")) or POLICY_VERSION,
                "schema_version": SCHEMA_VERSION,
            }
            lib_rows.append(row); all_rows.append(row)
            if (
                applied
                and reconciled_g != c["current_donor"]
                and singlet_relationship != "EXPECTED_COMPOSITE_COMPONENT_SINGLET"
                and uid_status in {"NO_LIBRARY_METADATA_MATCH", "METADATA_CONFLICT", "MISSING_UID_MAPPING"}
            ):
                amendments.append({
                    "library": lib, "proposed_donor_genotype": reconciled_g, "uid_resolution_scope": row["uid_resolution_scope"],
                    "reconciled_uid": uid, "uid_candidate_count": row["uid_candidate_count"], "uid_candidates": row["uid_candidates"],
                    "uid_resolution_status": uid_status, "proposed_metadata_event": row["event_class"],
                    "supporting_libraries": lib, "confidence": confidence, "reason": row["decision_reason"],
                })

        lib_dir = out_root
        write_tsv(str(lib_dir / f"{lib}.reconciled_cells.tsv.gz"), lib_rows, CELL_FIELDS)
        write_tsv(str(lib_dir / f"{lib}.changed_cells.tsv.gz"), [r for r in lib_rows if r["reassignment_applied"]], CELL_FIELDS)
        write_tsv(str(lib_dir / f"{lib}.technical_multiplets.tsv.gz"), [r for r in lib_rows if r["reconciled_droplet_state"] == "TECHNICAL_MULTIPLET"], CELL_FIELDS)
        write_tsv(str(lib_dir / f"{lib}.review_queue.tsv.gz"), [r for r in lib_rows if (not r["reassignment_applied"]) and r["final_action"] != "KEEP"], CELL_FIELDS)
        if args.write_reconciled_assignments:
            simple = []
            for r in lib_rows:
                if r["reconciled_droplet_state"] != "SINGLE_CELL" or not clean(r["reconciled_donor_genotype"]): continue
                comps = donor_components(r["reconciled_donor_genotype"])
                typ = "D" if clean(r["reconciled_biological_ploidy"]).upper() == "TETRAPLOID" else "S"
                score = r["nuclear_delta_ll"] if r["reassignment_applied"] else "NA"
                simple.append((r["barcode"], r["reconciled_donor_genotype"], typ, score))
            write_headerless_tsv(str(lib_dir / f"{lib}.reconciled_single_cell.assignments"), simple)
        # concise per-library report
        changed = sum(bool(r["reassignment_applied"]) for r in lib_rows)
        decisive_proposed = sum(r["decision_confidence"] == "DECISIVE" and not r["reassignment_applied"] and r["final_action"] != "KEEP" for r in lib_rows)
        below_event_mass_review = sum(
            r["final_action"] == "REVIEW_UNEXPECTED_IDENTITY"
            and r.get("singlet_library_relationship") != "EXPECTED_COMPOSITE_COMPONENT_SINGLET"
            and "ALTERNATIVE_LINE_BELOW_EVENT_MASS_THRESHOLD" in r["decision_reason_codes"]
            for r in lib_rows
        )
        below_autoapply_mass_review = sum(
            r["final_action"] == "REVIEW_UNEXPECTED_IDENTITY"
            and r.get("singlet_library_relationship") != "EXPECTED_COMPOSITE_COMPONENT_SINGLET"
            and "ALTERNATIVE_LINE_BELOW_AUTOAPPLY_MASS_THRESHOLD" in r["decision_reason_codes"]
            for r in lib_rows
        )
        exact_events_clearing_mass = sum(
            1 for (elib, _), e in line_event_context.items()
            if elib == lib and e.get("event_confidence") in {"STRONG", "DECISIVE"}
        )
        residual_component_singlet_events = sum(
            1 for (elib, _), e in line_event_context.items()
            if elib == lib and e.get("event_class") == "EXPECTED_COMPOSITE_COMPONENT_SINGLET_POPULATION"
        )
        unexpected_singlet_events = sum(
            1 for (elib, _), e in line_event_context.items()
            if elib == lib and e.get("event_class") == "LIKELY_UNEXPECTED_SINGLET_POPULATION"
        )
        multip = sum(r["reconciled_droplet_state"] == "TECHNICAL_MULTIPLET" for r in lib_rows)
        pair_origin_reviews = sum(r["final_action"] == "REVIEW_CELLULAR_ORIGIN" for r in lib_rows)
        homotet_origin_reviews = sum(r["final_action"] == "REVIEW_HOMOTET_OCCUPANCY" for r in lib_rows)
        tech_candidates_path = os.path.join(args.candidate_root, f"{lib}.technical_multiplet_candidates.tsv.gz")
        tech_candidates_n = len(read_optional_tsv(tech_candidates_path))
        dd = cells[0].get("dd_summary", {}) if cells else {}
        dd_status = clean(dd.get("status")) or "NOT_RUN"
        dd_rate = clean(dd.get("doublet_rate")) or "NA"
        mt_disjoint_support = sum(
            c.get("mt_verification_mode") == "DISJOINT_REAL_LINE_COMPONENT" and c.get("mt_alt") == "SUPPORTS_ALTERNATIVE"
            for c in cells
        )
        (rep_root / f"{lib}.md").write_text(
            f"# {lib} identity reconciliation\n\nCells: {len(lib_rows)}  \nApplied changes: {changed}  \nDecisive proposals held for review: {decisive_proposed}  \nExpected singlets in library metadata: {len(explicit_expected_singlets.get(lib, set()))}  \nResidual expected-composite component singlet populations: {residual_component_singlet_events}  \nGenuinely unexpected singlet populations: {unexpected_singlet_events}  \nDonor-pair cellular-origin ambiguities held: {pair_origin_reviews}  \nHomotet/same-donor-doublet ambiguities held: {homotet_origin_reviews}  \nUnexpected-identity observations below event mass (<{int(policy.get('event_min_cells', 10))}): {below_event_mass_review}  \nCoherent unexpected identities held below auto-apply mass (<{int(policy.get('unexpected_line_autoapply_min_cells', 100))}): {below_autoapply_mass_review}  \nExact identity events clearing event mass: {exact_events_clearing_mass}  \nDisjoint-line MT corroborations: {mt_disjoint_support}  \nGenotype-identifiable current technical multiplets: {multip}  \nAdditional technical-multiplet candidates: {tech_candidates_n}  \nDoublet Dragon context: {dd_status} (diploid-resolvable rate={dd_rate})  \nEvidence mode: {args.evidence_mode}  \nAuto-apply enabled: {'yes' if auto_apply else 'no'}  \nAmbient RNA evaluated: no\n\nResidual parental/component singlets are donors not explicitly expected as singlets but represented in an expected composite genotype in this library.  \nUnexpected singlets are globally real donors absent from the library expected donor-component universe.\n"
        )

    # Event rows after decisions. Exact-identity events are emitted separately
    # from component events so downstream review can distinguish the two scopes.
    combined_events = list(event_context.items()) + list(line_event_context.items())
    for (lib, component), ev in sorted(combined_events, key=lambda x: (x[0][0], x[1].get("event_scope", ""), x[0][1])):
        imp = ev["implicated"]; n = len(imp)
        event_identity = canonical_genotype(ev.get("unexpected_component") or component)
        uid_meta = proposed_uid_metadata(lib, event_identity, meta, global_meta)
        event_rows.append({
            "event_id": ev["event_id"], "library": lib, "event_scope": ev.get("event_scope", "COMPONENT"), "unexpected_component": event_identity, "n_implicated_cells": n,
            "fraction_library_implicated": n / ev["nlib"] if ev["nlib"] else math.nan,
            "primary_source_identity": ev["primary_source"], "primary_source_fraction": ev["source_frac"], "source_identity_entropy": ev["source_entropy"],
            "primary_replaced_component": Counter(c["replaced_component"] for c in imp if c["replaced_component"]).most_common(1)[0][0] if any(c["replaced_component"] for c in imp) else "",
            "fraction_primary_source_displaced": ev["displaced"], "n_primary_source_remaining": max(0, ev["source_total"] - ev["source_implicated"]),
            "fraction_event_from_primary_source": ev["source_frac"], "partner_preservation_fraction": ev["preserved"],
            "transition_type_summary": ",".join(f"{k}:{v}" for k,v in Counter((c["source_identity"] or c["current_donor"]) + "->" + (c["best_alt"] or "") for c in imp).most_common()),
            "fraction_nn_diploid": frac(imp, lambda c: clean(c["nnrow"].get("ploidy_call")).lower().startswith("dip")),
            "fraction_nn_high_tet": frac(imp, lambda c: math.isfinite(c["nn_prob"]) and c["nn_prob"] >= float(policy.get("nn_homotypic_prob_tet", 0.9))),
            "median_nn_prob_tet": med(c["nn_prob"] for c in imp),
            "fraction_mt_supports_alternative": frac(imp, lambda c: c["mt_alt"] == "SUPPORTS_ALTERNATIVE"),
            "fraction_mt_supports_current": frac(imp, lambda c: c["mt_alt"] == "SUPPORTS_CURRENT"),
            "fraction_mt_unresolved": frac(imp, lambda c: c["mt_alt"] in {"MISSING", "NONDISCRIMINATING"}),
            "fraction_atac_supports_alternative": frac(imp, lambda c: c["atac_status"] == "ATAC_SUPPORTS_ALTERNATIVE"),
            "fraction_atac_supports_current": frac(imp, lambda c: c["atac_status"] == "ATAC_SUPPORTS_CURRENT"),
            "fraction_atac_unavailable": frac(imp, lambda c: c["atac_status"] in {"ATAC_NOT_REQUESTED", "ATAC_UNAVAILABLE", "ATAC_INSUFFICIENT"}),
            "fraction_site_fold_replicated": frac(imp, lambda c: math.isfinite(c["fold_fraction"]) and c["fold_fraction"] >= float(policy.get("site_fold_support_fraction", 0.6))),
            "median_depth_implicated": med(c["nuclear_depth"] for c in imp), "median_depth_background": med(ev["background_depth"]),
            "fraction_low_information_implicated": frac(imp, lambda c: c["nuclear_depth"] < float(policy.get("nuclear_low_information_depth", 10))),
            "fraction_low_information_background": frac(ev["background_depth"], lambda d: d < float(policy.get("nuclear_low_information_depth", 10))),
            "species_consistency": frac(imp, lambda c: c["species_alt"] != "CONTRADICTS"),
            "related_library_recurrence": "YES" if ev.get("recurrence") else "NO", "reciprocal_event_id": "",
            "nuclear_confusability_status": (
                "INDISTINGUISHABLE_ON_PANEL" if any(c.get("nuclear_confusability_status") == "INDISTINGUISHABLE_ON_PANEL" for c in imp)
                else "DISTINGUISHABLE_ON_PANEL" if imp and all(c.get("nuclear_confusability_status") in {"DISTINGUISHABLE_ON_PANEL", "NOT_APPLICABLE"} for c in imp)
                else "MIXED_OR_NOT_EVALUATED"
            ), "reconciled_uid": uid_meta.get("uid", ""),
            "uid_candidate_count": uid_meta.get("uid_candidate_count", 0), "uid_candidates": uid_meta.get("uid_candidates", ""),
            "uid_resolution_status": uid_meta.get("uid_resolution_status", "NO_LIBRARY_METADATA_MATCH"), "uid_resolution_basis": uid_meta.get("uid_resolution_basis", ""),
            "uid_resolution_scope": uid_meta.get("uid_resolution_scope", ""), "uid_member_details": "",
            "singlet_library_relationship": ev.get("singlet_library_relationship", "NOT_SINGLET"),
            "expected_composite_context": ev.get("expected_composite_context", ""),
            "contributes_to_library_exchange_evidence": bool(ev.get("contributes_to_library_exchange_evidence", False)),
            "event_class": ev["event_class"], "event_confidence": ev["event_confidence"], "event_reason": ev["event_reason"],
        })

    # Pairwise library-exchange inference.  Only robust exact-identity singlet
    # events that are genuinely absent from the target library donor-component
    # universe contribute foreign-signature evidence.  Shared expected donors
    # are excluded by construction through the pairwise unique-donor sets.
    component_barcodes = defaultdict(lambda: defaultdict(set))
    singlet_cell_counts = defaultdict(Counter)
    for row in all_rows:
        if clean(row.get("reconciled_droplet_state")) != "SINGLE_CELL":
            continue
        lib = clean(row.get("library")); bc = clean(row.get("barcode"))
        genotype = canonical_genotype(row.get("reconciled_donor_genotype", ""))
        comps = donor_components(genotype)
        for donor in set(comps):
            component_barcodes[lib][donor].add(bc)
        if len(comps) == 1:
            singlet_cell_counts[lib][comps[0]] += 1

    foreign_event_cells = defaultdict(Counter)
    foreign_event_ids = defaultdict(dict)
    for ev in event_rows:
        if clean(ev.get("event_scope")) != "EXACT_IDENTITY":
            continue
        if clean(ev.get("singlet_library_relationship")) != "UNEXPECTED_SINGLET":
            continue
        if clean(ev.get("contributes_to_library_exchange_evidence")).upper() != "TRUE" and ev.get("contributes_to_library_exchange_evidence") is not True:
            continue
        donor = canonical_genotype(ev.get("unexpected_component", ""))
        if len(donor_components(donor)) != 1:
            continue
        lib = clean(ev.get("library"))
        n_cells = fint(ev.get("n_implicated_cells"))
        if n_cells < min_event_cells:
            continue
        foreign_event_cells[lib][donor] = max(foreign_event_cells[lib][donor], n_cells)
        foreign_event_ids[lib][donor] = clean(ev.get("event_id"))

    def joined(values):
        return ";".join(sorted(set(values), key=natural_key))

    def coverage(numerator: int, denominator: int) -> float:
        return numerator / denominator if denominator else math.nan

    def native_support(lib: str, donors: set):
        detected = {d for d in donors if len(component_barcodes[lib].get(d, set())) >= min_event_cells}
        cells = set()
        for donor in donors:
            cells.update(component_barcodes[lib].get(donor, set()))
        return detected, len(cells)

    exchange_rows = []
    exchange_donor_rows = []
    lib_names = [f"lib{n}" for n in libs]
    for i, library_a in enumerate(lib_names):
        for library_b in lib_names[i + 1:]:
            roster = pairwise_library_roster(library_a, library_b, expected_donor_components)
            shared = set(roster["shared"])
            a_specific = set(roster["a_specific"])
            b_specific = set(roster["b_specific"])

            a_native_detected, a_native_cells = native_support(library_a, a_specific)
            b_native_detected, b_native_cells = native_support(library_b, b_specific)
            a_foreign_in_b = {d for d in a_specific if foreign_event_cells[library_b].get(d, 0) >= min_event_cells}
            b_foreign_in_a = {d for d in b_specific if foreign_event_cells[library_a].get(d, 0) >= min_event_cells}
            a_foreign_cells_in_b = sum(foreign_event_cells[library_b].get(d, 0) for d in a_foreign_in_b)
            b_foreign_cells_in_a = sum(foreign_event_cells[library_a].get(d, 0) for d in b_foreign_in_a)

            a_cov_in_b = coverage(len(a_foreign_in_b), len(a_specific))
            b_cov_in_a = coverage(len(b_foreign_in_a), len(b_specific))
            a_native_retention = coverage(len(a_native_detected), len(a_specific))
            b_native_retention = coverage(len(b_native_detected), len(b_specific))
            a_native_displacement = 1.0 - a_native_retention if math.isfinite(a_native_retention) else math.nan
            b_native_displacement = 1.0 - b_native_retention if math.isfinite(b_native_retention) else math.nan
            a_to_b_strength = directional_exchange_strength(len(a_foreign_in_b), len(a_specific))
            b_to_a_strength = directional_exchange_strength(len(b_foreign_in_a), len(b_specific))
            reciprocal, interpretation, confidence = classify_library_exchange(
                str(roster["roster_relation"]), str(roster["pair_discriminability"]),
                a_to_b_strength, b_to_a_strength, a_native_retention, b_native_retention,
            )
            supporting_ids = {
                foreign_event_ids[library_b].get(d, "") for d in a_foreign_in_b
            } | {
                foreign_event_ids[library_a].get(d, "") for d in b_foreign_in_a
            }
            supporting_ids.discard("")
            notes = "Shared donors are excluded from library-exchange evidence."
            if (not a_specific or not b_specific) and roster["roster_relation"] != "ROSTER_EQUIVALENT_NONDISCRIMINATING":
                notes += " Only one directional donor fingerprint is identifiable for this nested roster pair."

            exchange_rows.append({
                "library_a": library_a, "library_b": library_b,
                "roster_relation": roster["roster_relation"], "pair_discriminability": roster["pair_discriminability"],
                "shared_donors": joined(shared), "a_specific_donors": joined(a_specific), "b_specific_donors": joined(b_specific),
                "n_shared_donors": len(shared), "n_a_specific_donors": len(a_specific), "n_b_specific_donors": len(b_specific),
                "a_specific_detected_in_a": len(a_native_detected), "a_specific_detected_in_b": len(a_foreign_in_b),
                "b_specific_detected_in_b": len(b_native_detected), "b_specific_detected_in_a": len(b_foreign_in_a),
                "a_specific_cells_in_a": a_native_cells, "a_specific_cells_in_b": a_foreign_cells_in_b,
                "b_specific_cells_in_b": b_native_cells, "b_specific_cells_in_a": b_foreign_cells_in_a,
                "a_signature_coverage_in_b": a_cov_in_b, "b_signature_coverage_in_a": b_cov_in_a,
                "a_native_retention_fraction": a_native_retention, "b_native_retention_fraction": b_native_retention,
                "a_native_displacement_fraction": a_native_displacement, "b_native_displacement_fraction": b_native_displacement,
                "a_to_b_exchange_strength": a_to_b_strength, "b_to_a_exchange_strength": b_to_a_strength,
                "reciprocal_exchange_status": reciprocal, "exchange_interpretation": interpretation,
                "exchange_confidence": confidence, "supporting_event_ids": joined(supporting_ids), "notes": notes,
            })

            for source_library, target_library, donors in (
                (library_a, library_b, a_specific),
                (library_b, library_a, b_specific),
            ):
                for donor in sorted(donors, key=natural_key):
                    source_relationship, _ = classify_singlet_library_relationship(
                        donor,
                        explicit_expected_singlets.get(source_library, set()),
                        expected_donor_components.get(source_library, set()),
                        expected_composites_by_component.get(source_library, {}),
                        globally_valid=(donor in global_donors),
                    )
                    target_relationship, _ = classify_singlet_library_relationship(
                        donor,
                        explicit_expected_singlets.get(target_library, set()),
                        expected_donor_components.get(target_library, set()),
                        expected_composites_by_component.get(target_library, {}),
                        globally_valid=(donor in global_donors),
                    )
                    source_singlet_cells = singlet_cell_counts[source_library].get(donor, 0)
                    if source_singlet_cells < min_event_cells:
                        source_singlet_cells = 0
                    target_singlet_cells = foreign_event_cells[target_library].get(donor, 0)
                    exchange_donor_rows.append({
                        "source_library": source_library, "target_library": target_library,
                        "diagnostic_donor": donor, "diagnostic_for_library": source_library,
                        "present_in_source_expected_components": True, "present_in_target_expected_components": False,
                        "shared_or_specific": "SPECIFIC",
                        "source_observed_robust_singlet_cells": source_singlet_cells,
                        "target_observed_robust_singlet_cells": target_singlet_cells,
                        "source_singlet_relationship": source_relationship, "target_singlet_relationship": target_relationship,
                        "source_event_id": foreign_event_ids[source_library].get(donor, ""),
                        "target_event_id": foreign_event_ids[target_library].get(donor, ""),
                    })

    exchange_candidates_by_lib = Counter()
    for row in exchange_rows:
        if row["exchange_interpretation"] in {"NO_LIBRARY_EXCHANGE_SIGNAL", "ROSTER_EQUIVALENT_NONDISCRIMINATING"}:
            continue
        exchange_candidates_by_lib[row["library_a"]] += 1
        exchange_candidates_by_lib[row["library_b"]] += 1
    for lib in lib_names:
        report_path = rep_root / f"{lib}.md"
        with open(report_path, "a") as fh:
            fh.write(f"\nLibrary-exchange candidate pairs: {exchange_candidates_by_lib.get(lib, 0)}  \n")

    write_tsv(str(out_root / "all_libraries.reconciled_cells.tsv.gz"), all_rows, CELL_FIELDS)
    write_tsv(str(out_root / "all_libraries.changed_cells.tsv.gz"), [r for r in all_rows if r["reassignment_applied"]], CELL_FIELDS)
    write_tsv(str(out_root / "all_libraries.identity_events.tsv"), event_rows, RECONCILE_EVENT_FIELDS)
    write_tsv(str(out_root / "all_libraries.library_exchange_events.tsv"), exchange_rows, LIBRARY_EXCHANGE_FIELDS)
    write_tsv(str(out_root / "all_libraries.library_exchange_donor_evidence.tsv"), exchange_donor_rows, LIBRARY_EXCHANGE_DONOR_FIELDS)
    write_tsv(str(out_root / "all_libraries.metadata_amendments_proposed.tsv"), amendments, META_AMEND_FIELDS)
    json_dump_atomic(str(out_root / "reconciliation_manifest.json"), {
        "schema_version": SCHEMA_VERSION, "policy_version": clean(policy.get("policy_version")) or POLICY_VERSION,
        "policy_sha256": policy_hash, "evidence_mode": args.evidence_mode, "libraries": libs,
        "ambient_rna_evaluated": False, "contam_dependency": False, "empty_drops_dependency": False,
        "tetra_refine_rerun_dependency": False, "atac_requested": args.evidence_mode == "rna-atac",
        "auto_apply_enabled": bool(auto_apply),
        "global_biological_line_source": "2025_LineMeta",
        "doublet_dragon_context_used": bool(args.doublet_context_root),
        "doublet_dragon_role": "population_context_only_diploid_resolvable_subset",
        "event_min_cells": int(policy.get("event_min_cells", 10)),
        "unexpected_line_autoapply_min_cells": int(policy.get("unexpected_line_autoapply_min_cells", 100)),
        "occupancy_aware_reconciliation": True,
        "diploid_population_definition": "expected_singlets_plus_robust_singlet_events_plus_component_derived_singlet_populations",
        "heterotypic_pair_occupancy_guard": bool(policy.get("heterotypic_pair_occupancy_guard", True)),
        "unexpected_homotet_occupancy_guard": bool(policy.get("unexpected_homotet_occupancy_guard", True)),
        "explicit_multiplet_evidence_source": clean(policy.get("explicit_multiplet_evidence_source")) or "droplet_flag_only",
        "component_singlet_classification": True,
        "library_exchange_inference": True,
        "library_exchange_evidence_source": "robust_unexpected_singlet_populations_unique_to_source_roster",
        "library_exchange_shared_donors_excluded": True,
    })
    n_changed = sum(bool(r["reassignment_applied"]) for r in all_rows)
    n_residual_component_events = sum(r["event_class"] == "EXPECTED_COMPOSITE_COMPONENT_SINGLET_POPULATION" for r in event_rows)
    n_unexpected_singlet_events = sum(r["event_class"] == "LIKELY_UNEXPECTED_SINGLET_POPULATION" for r in event_rows)
    n_exchange_candidates = sum(
        r["exchange_interpretation"] not in {"NO_LIBRARY_EXCHANGE_SIGNAL", "ROSTER_EQUIVALENT_NONDISCRIMINATING"}
        for r in exchange_rows
    )
    summary = (
        "# Identity reconciliation summary\n\n"
        f"Libraries: {len(libs)}  \nCells: {len(all_rows)}  \nApplied changes: {n_changed}  \n"
        f"Events: {len(event_rows)}  \nExpected singlet definitions: {sum(len(v) for v in explicit_expected_singlets.values())}  \n"
        f"Residual expected-composite component singlet populations: {n_residual_component_events}  \n"
        f"Genuinely unexpected singlet populations: {n_unexpected_singlet_events}  \n"
        f"Library-exchange candidate pairs: {n_exchange_candidates}  \nEvidence mode: {args.evidence_mode}  \n"
        "Ambient RNA evaluated: no  \nOriginal assignments modified: no\n\n"
        "Residual parental/component singlets are donors not explicitly expected as singlets but represented in an expected composite genotype in the same library.  \n"
        "Unexpected singlets are globally real donors absent from the library's complete expected donor-component universe.\n"
    )
    (rep_root / "summary.md").write_text(summary)
    print(f"Wrote {len(all_rows)} reconciled cell rows; {n_changed} applied changes")
    return 0



# -----------------------------------------------------------------------------
# atac-barcode-map
# -----------------------------------------------------------------------------
import argparse, gzip, hashlib
from pathlib import Path


def atac_open_text(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def atac_barcode_map_main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--rna-whitelist", required=True); p.add_argument("--atac-whitelist", required=True); p.add_argument("--output", required=True)
    a=p.parse_args()
    with atac_open_text(a.rna_whitelist) as r: rna=[x.strip().split()[0] for x in r if x.strip()]
    with atac_open_text(a.atac_whitelist) as r: atac=[x.strip().split()[0] for x in r if x.strip()]
    if len(rna)!=len(atac): raise SystemExit(f"whitelist lengths differ: RNA={len(rna)} ATAC={len(atac)}")
    if len(set(rna))!=len(rna) or len(set(atac))!=len(atac): raise SystemExit("whitelist contains duplicate barcodes")
    out=Path(a.output); out.parent.mkdir(parents=True,exist_ok=True)
    with open(out,"w") as fh:
        fh.write("atac_barcode\trna_barcode\n")
        for x,y in zip(atac,rna): fh.write(f"{x}\t{y}\n")
    print(f"Wrote {len(rna)} barcode pairs to {out}")
    return 0

# -----------------------------------------------------------------------------
# atac-finalize
# -----------------------------------------------------------------------------
import argparse, math
from pathlib import Path
from identity_reconciliation_common import clean, ffloat, read_tsv, sha256_file, write_tsv


def atac_finalize_main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--scores', required=True); p.add_argument('--qc', required=True); p.add_argument('--output', required=True)
    p.add_argument('--barcode-map', default=''); p.add_argument('--support-delta-ll', type=float, default=8.0); p.add_argument('--min-depth', type=float, default=20.0)
    a=p.parse_args(); rows=read_tsv(a.scores); qrows=read_tsv(a.qc); qc=qrows[0] if qrows else {}
    qc_map=clean(qc.get('atac_barcode_map_path'))
    map_path=qc_map if qc_map not in {'', 'NA', '.'} else ''
    map_hash=sha256_file(map_path) if map_path and Path(map_path).is_file() else 'NA'
    out=[]
    for r in rows:
        delta=ffloat(r.get('atac_delta_ll_vs_current')); depth=ffloat(r.get('atac_informative_depth'),0.0); base=clean(r.get('atac_score_status'))
        if base not in {'PASS',''} or depth < a.min_depth: status='ATAC_INSUFFICIENT'
        elif math.isfinite(delta) and delta >= a.support_delta_ll: status='ATAC_SUPPORTS_ALTERNATIVE'
        elif math.isfinite(delta) and delta <= -a.support_delta_ll: status='ATAC_SUPPORTS_CURRENT'
        else: status='ATAC_NONDISCRIMINATING'
        z=dict(r); z.update({
            'atac_requested':'1','atac_barcode_mode':clean(qc.get('atac_barcode_mode')) or 'direct',
            'atac_barcode_map_path':map_path or 'NA','atac_barcode_map_sha256':map_hash,
            'atac_direct_overlap_fraction':qc.get('atac_direct_overlap_fraction','NA'),'atac_mapped_overlap_fraction':qc.get('atac_mapped_overlap_fraction','NA'),
            'atac_status':status,
        }); out.append(z)
    fields=list(out[0].keys()) if out else ['library','barcode','hypothesis_id','atac_status']
    write_tsv(a.output,out,fields)
    return 0

# -----------------------------------------------------------------------------
# optional-status
# -----------------------------------------------------------------------------
import argparse
from identity_reconciliation_common import read_tsv, write_tsv

def optional_status_main():
    p=argparse.ArgumentParser(description=__doc__); p.add_argument('--candidate-manifest',required=True);p.add_argument('--modality',choices=['mt','atac'],required=True);p.add_argument('--status',required=True);p.add_argument('--output',required=True);a=p.parse_args()
    c=read_tsv(a.candidate_manifest); out=[]
    for r in c:
        base={k:r.get(k,'') for k in ('library','barcode','hypothesis_id','state_notation','donor_genotype')}
        if a.modality=='mt':
            if len(donor_components(canonical_genotype(r.get('donor_genotype','')))) != 1:
                continue
            base.update({'mt_log_likelihood':'NA','mt_delta_hypothesis_vs_current':'NA','mt_delta_vs_best_other_singlet':'NA','mt_rank_within_singlet_candidates':'0','mt_sites_used':'0','mt_molecules_used':'0','mt_best_identity':'','mt_second_identity':'','mt_haplotype_resolution':'MT_HAPLOTYPE_UNRESOLVED' if a.status=='MT_HAPLOTYPE_UNRESOLVED' else 'MT_UNAVAILABLE','mt_fit_status':a.status,'schema_version':'mt_identity_scores_v3_donor_probe'})
        else: base.update({'atac_log_likelihood':'NA','atac_delta_ll_vs_current':'NA','atac_rank_within_candidates':'0','atac_informative_depth':'0','atac_informative_units':'0','atac_dosage_concordance':'NA','atac_residual_mismatch':'NA','atac_depth_normalized_delta':'NA','comparison_state':'UNAVAILABLE','atac_score_status':'UNAVAILABLE','atac_requested':'0' if a.status=='ATAC_NOT_REQUESTED' else '1','atac_barcode_mode':'not_requested' if a.status=='ATAC_NOT_REQUESTED' else 'unavailable','atac_barcode_map_path':'NA','atac_barcode_map_sha256':'NA','atac_direct_overlap_fraction':'NA','atac_mapped_overlap_fraction':'NA','atac_status':a.status,'schema_version':'identity_hypothesis_scores_v1'})
        out.append(base)
    write_tsv(a.output,out,list(out[0].keys()) if out else ['library','barcode','hypothesis_id'])
    return 0


_COMMANDS = {
    "metadata": metadata_main,
    "candidates": candidates_main,
    "doublet-context": doublet_context_main,
    "reconcile": reconcile_main,
    "optional-status": optional_status_main,
    "atac-barcode-map": atac_barcode_map_main,
    "atac-finalize": atac_finalize_main,
}

def main():
    if len(sys.argv) < 2 or sys.argv[1] in {"-h", "--help"}:
        print("usage: identity_reconciliation.py <command> [options]")
        print("commands: " + ", ".join(_COMMANDS))
        return 0
    command = sys.argv[1]
    func = _COMMANDS.get(command)
    if func is None:
        print(f"unknown command: {command}", file=sys.stderr)
        print("commands: " + ", ".join(_COMMANDS), file=sys.stderr)
        return 2
    sys.argv = [f"{sys.argv[0]} {command}"] + sys.argv[2:]
    return func()

if __name__ == "__main__":
    raise SystemExit(main())
