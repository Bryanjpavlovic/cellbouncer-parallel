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
import csv
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

CANDIDATE_AUDIT_REQUIRED_FIELDS = [
    "library", "barcode", "event_id", "candidate_raw",
    "candidate_canonical", "candidate_kind", "candidate_components",
    "candidate_source", "candidate_tier", "candidate_rank_within_source",
    "candidate_eligibility", "candidate_set_aside_reason",
    "selected_as_proposal", "selected_for_candidate_axis",
    "axis_endpoint_role", "lower_rank_considered_after_set_aside",
]
CANDIDATE_AUDIT_FIELDS = list(dict.fromkeys(
    CANDIDATE_AUDIT_REQUIRED_FIELDS + FIELDS + TECHNICAL_MULTIPLET_FIELDS
))

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
        candidate_audit_rows = []
        candidate_audit_by_key = {}
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
            raw_labels: Dict[str, Set[str]] = defaultdict(set)
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
                raw = clean(g)
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
                raw_labels[g].add(raw or g)
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
            output_by_candidate = {}
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
                candidate_output = {
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
                }
                rows.append(candidate_output)
                output_by_candidate[g] = candidate_output
                if len(comps) <= 2: universe.add(g)

            source_rank = defaultdict(int)
            ranks_by_candidate = defaultdict(list)
            for g in ordered:
                for source in sorted(origins[g], key=lambda x: (priority_order.get(x, 99), x)):
                    source_rank[source] += 1
                    ranks_by_candidate[g].append(f"{source}:{source_rank[source]}")
            for material_rank, g in enumerate(ordered, 1):
                project_status, admissibility = candidate_biological_status(
                    g, expected, global_lines, global_donors
                )
                candidate_kind = (
                    "CURRENT_TECHNICAL_COMPOSITION"
                    if admissibility not in {
                        "SINGLET_IDENTITY_CANDIDATE",
                        "BIOLOGICAL_SINGLE_CELL_ALLOWED",
                    } else
                    "EXPECTED_LIBRARY_BIOLOGICAL_IDENTITY"
                    if g in expected else
                    "GLOBALLY_REAL_UNEXPECTED_BIOLOGICAL_IDENTITY"
                )
                evidence = dict(output_by_candidate.get(g, {}))
                if not evidence:
                    evidence.update({
                        "library": lib,
                        "barcode": bc,
                        "state_notation": state_for(g),
                        "donor_genotype": g,
                        "donor_components": ",".join(donor_components(g)),
                        "biological_ploidy": ploidy_for(g),
                        "droplet_state": "SINGLE_CELL_CANDIDATE",
                        "candidate_origin": ",".join(sorted(
                            origins[g], key=lambda x: (priority_order.get(x, 99), x))),
                        "current_state_notation": current_state,
                        "current_donor_genotype": cur_g,
                        "expected_genotype_status": (
                            "EXPECTED" if g in expected else
                            "UNEXPECTED_OR_EVENT_NOMINATED"),
                        "project_genotype_status": project_status,
                        "biological_admissibility": admissibility,
                        "candidate_priority": material_rank,
                        "candidate_generation_notes": ";".join(
                            sorted(set(notes[g]))),
                        "schema_version": SCHEMA_VERSION,
                    })
                evidence.update({
                    "library": lib,
                    "barcode": bc,
                    "event_id": "",
                    "candidate_raw": "|".join(sorted(raw_labels[g], key=natural_key)) or g,
                    "candidate_canonical": g,
                    "candidate_kind": candidate_kind,
                    "candidate_components": ",".join(donor_components(g)),
                    "candidate_source": ",".join(sorted(
                        origins[g], key=lambda x: (priority_order.get(x, 99), x))),
                    "candidate_tier": str(min(
                        priority_order.get(source, 99) for source in origins[g])),
                    "candidate_rank_within_source": ";".join(ranks_by_candidate[g]),
                    "candidate_eligibility": (
                        "SCOREABLE" if g in output_by_candidate else "SET_ASIDE"),
                    "candidate_set_aside_reason": (
                        "" if g in output_by_candidate else "MAX_CANDIDATES_CAP"),
                    "selected_as_proposal": False,
                    "selected_for_candidate_axis": False,
                    "axis_endpoint_role": "",
                    "lower_rank_considered_after_set_aside": False,
                })
                candidate_audit_rows.append(evidence)
                candidate_audit_by_key[(bc, g)] = evidence

        out_path = out_root / f"{lib}.identity_candidates.tsv.gz"
        write_tsv(str(out_path), rows, FIELDS)
        technical_rows = list(technical_by_key.values())
        write_tsv(str(out_root / f"{lib}.technical_multiplet_candidates.tsv.gz"),
                  technical_rows, TECHNICAL_MULTIPLET_FIELDS)
        technical_rank = Counter()
        for technical in sorted(
                technical_rows,
                key=lambda row: (
                    natural_key(clean(row.get("barcode"))),
                    natural_key(clean(row.get("donor_composition"))),
                    clean(row.get("technical_class")))):
            barcode = clean(technical.get("barcode"))
            technical_rank[barcode] += 1
            canonical = canonical_genotype(technical.get("donor_composition", ""))
            audit = dict(technical)
            audit.update({
                "library": lib,
                "barcode": barcode,
                "event_id": "",
                "candidate_raw": clean(technical.get("donor_composition")),
                "candidate_canonical": canonical,
                "candidate_kind": "TECHNICAL_MULTIPLET_HYPOTHESIS",
                "candidate_components": ",".join(donor_components(canonical)),
                "candidate_source": clean(technical.get("candidate_origin")),
                "candidate_tier": "TECHNICAL",
                "candidate_rank_within_source": (
                    f"TECHNICAL:{technical_rank[barcode]}"),
                "candidate_eligibility": "TECHNICAL_ONLY",
                "candidate_set_aside_reason": "NOT_A_BIOLOGICAL_SINGLE_CELL_LINE",
                "selected_as_proposal": False,
                "selected_for_candidate_axis": False,
                "axis_endpoint_role": "",
                "lower_rank_considered_after_set_aside": False,
            })
            existing = candidate_audit_by_key.get((barcode, canonical))
            if existing is None:
                candidate_audit_rows.append(audit)
                candidate_audit_by_key[(barcode, canonical)] = audit
            else:
                existing["candidate_raw"] = "|".join(sorted({
                    value for value in (
                        clean(existing.get("candidate_raw")),
                        clean(audit.get("candidate_raw")))
                    if value
                }, key=natural_key))
                existing["candidate_source"] = ",".join(sorted({
                    value
                    for source in (
                        clean(existing.get("candidate_source")),
                        clean(audit.get("candidate_source")))
                    for value in source.split(",") if value
                }, key=natural_key))
                existing["candidate_rank_within_source"] = ";".join(sorted({
                    value
                    for ranks in (
                        clean(existing.get("candidate_rank_within_source")),
                        clean(audit.get("candidate_rank_within_source")))
                    for value in ranks.split(";") if value
                }, key=natural_key))
                for field, value in technical.items():
                    if not clean(existing.get(field)) and clean(value):
                        existing[field] = value
        audit_keys = [
            (clean(row.get("barcode")),
             canonical_genotype(row.get("candidate_canonical", "")))
            for row in candidate_audit_rows
        ]
        if len(audit_keys) != len(set(audit_keys)):
            raise ValueError(f"{lib} candidate audit contains duplicate canonical candidates")
        write_tsv(
            str(out_root / f"{lib}.identity_candidate_audit.tsv.gz"),
            candidate_audit_rows, CANDIDATE_AUDIT_FIELDS,
        )
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
    p.add_argument("--atac-root", default="")
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
    if args.evidence_mode == "rna-atac" and not args.atac_root:
        raise SystemExit("--atac-root is required only for rna-atac evidence mode")
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
        atac = (
            score_map(os.path.join(
                args.atac_root, f"{lib}.atac_identity_scores.tsv.gz"))
            if args.evidence_mode == "rna-atac" else defaultdict(dict)
        )
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


# -----------------------------------------------------------------------------
# score-pairs
# -----------------------------------------------------------------------------

SCORE_PAIR_SCHEMA = "identity_reconciliation_score_pair_manifest_v2"
SCORE_PAIR_CONTRACT = (
    "ORIGINAL_ALLOWED_DEMUX_VS_RECONCILIATION_NOMINATED_SWAP_ONLY"
)
SCORE_PAIR_ROLES = (
    "ORIGINAL_ALLOWED_DEMUX",
    "RECONCILIATION_NOMINATED_SWAP",
)
SCORE_PAIR_SUPPORTED_SWAP_EVENT_CLASSES = {
    "LIKELY_UNEXPECTED_INTACT_BIOLOGICAL_LINE",
    "LIKELY_UNEXPECTED_SINGLET_POPULATION",
}
SCORE_PAIR_APPLIED_REASSIGNMENT = "APPLIED_REASSIGNMENT"
SCORE_PAIR_RECOMMENDED_NOT_APPLIED = "RECOMMENDED_NOT_APPLIED"
SCORE_PAIR_SUPPORTED_EVENT_HELD_CELL = "SUPPORTED_EVENT_HELD_CELL"
SCORE_PAIR_REVIEW_ONLY = "REVIEW_ONLY_UNEXPECTED_IDENTITY"
SCORE_PAIR_REASSIGNMENT_SCOPES = {
    SCORE_PAIR_APPLIED_REASSIGNMENT,
    SCORE_PAIR_RECOMMENDED_NOT_APPLIED,
}
SCORE_PAIR_EXTRA_FIELDS = [
    "score_pair_id", "score_pair_role", "score_pair_source",
    "score_population_scope", "supported_event_key",
    "reconciliation_event_id", "reconciliation_event_class",
    "reconciliation_event_confidence", "reconciliation_final_action",
    "reconciliation_decision_confidence",
    "reconciliation_reassignment_applied",
    "reconciliation_current_refined_assignment",
    "original_demux_assignment", "proposed_donor_genotype",
    "score_scope_contract",
]
SCORE_PAIR_SUMMARY_FIELDS = [
    "library", "source_decisions", "source_decisions_sha256",
    "validation_summary", "validation_summary_sha256",
    "score_pair_builder", "score_pair_builder_sha256",
    "metadata_expected_genotypes_sha256",
    "metadata_global_biological_lines_sha256", "metadata_global_donors_sha256",
    "n_decision_rows", "n_not_reconciliation_nominated",
    "n_candidate_rows_considered", "n_supported_event_keys",
    "n_supported_reassignment_pairs", "n_supported_event_held_pairs",
    "n_review_only_pairs", "n_score_pairs", "n_manifest_rows",
    "n_candidate_rows_excluded", "exclusion_reason_counts",
    "score_scope_contract", "schema_version",
]
SCORE_PAIR_EXCLUSION_FIELDS = [
    "library", "barcode", "original_demux_assignment",
    "proposed_donor_genotype", "reconciliation_event_id",
    "reconciliation_event_class", "reconciliation_final_action",
    "exclusion_reason", "detail", "score_scope_contract", "schema_version",
]


def score_pairs_parse_args():
    p = argparse.ArgumentParser(
        description=(
            "Build exact post-reconciliation identity-score pairs. Only the "
            "original allowed demux winner and the reconciliation-nominated "
            "biological swap are emitted."
        )
    )
    p.add_argument("--libraries", nargs="+", required=True)
    p.add_argument("--decisions-root", required=True)
    p.add_argument("--metadata-root", required=True)
    p.add_argument("--validation-summary", required=True)
    p.add_argument("--output-root", required=True)
    return p.parse_args()


def _score_pair_ploidy(genotype: str) -> str:
    components = donor_components(genotype)
    if len(components) == 1:
        return "DIPLOID"
    if len(components) == 2:
        return "TETRAPLOID"
    return "UNRESOLVED"


def _score_pair_same_snp_identity(left: str, right: str) -> bool:
    """Return true for A versus A+A, which nuclear SNP dosage cannot resolve."""
    left_components = donor_components(left)
    right_components = donor_components(right)
    return (
        bool(left_components)
        and bool(right_components)
        and len(set(left_components)) == 1
        and len(set(right_components)) == 1
        and left_components[0] == right_components[0]
    )


def _score_pair_event_key(source: Mapping[str, object]) -> Tuple[str, str]:
    """Return the frozen reconciliation event and exact proposed identity."""
    return (
        clean(source.get("event_id")),
        canonical_genotype(source.get("proposed_donor_genotype", "")),
    )


def _score_pair_event_key_text(
    library: str, event_key: Tuple[str, str]
) -> str:
    event_id, proposed = event_key
    return f"{library}|{event_id}|{proposed}"


def _score_pair_nonidentity_event_reason(
    source: Mapping[str, object],
) -> Tuple[str, str]:
    event_class = clean(source.get("event_class")).upper()
    action = clean(source.get("final_action")).upper()
    if "HOMOTET" in event_class or "OCCUPANCY" in event_class:
        return (
            "HOMOTET_OCCUPANCY_REVIEW_NOT_IDENTITY_SWAP",
            "same-donor ploidy/occupancy ambiguity is handled by the NN route",
        )
    if "PLOIDY" in event_class or action == "RECLASSIFY_PLOIDY":
        return (
            "PLOIDY_RECLASSIFICATION_NOT_IDENTITY_SWAP",
            "ploidy reclassification is not an original-versus-swap comparison",
        )
    if event_class in {
        "NO_METADATA_PROBLEM",
        "EXPECTED_COMPOSITE_COMPONENT_SINGLET_POPULATION",
    }:
        return (
            "NOT_A_SAMPLE_SWAP_RECONCILIATION_EVENT",
            f"event_class={event_class or 'MISSING'}",
        )
    return "", ""


def _score_pair_population_scope(
    source: Mapping[str, object],
    supported_event_keys: Set[Tuple[str, str]],
) -> str:
    """Apply the established v10 decision boundary without new thresholds.

    The supported-event report includes REASSIGN_GENOTYPE decisions plus
    REVIEW_CELLULAR_ORIGIN context cells attached to the same event and exact
    proposed identity.  Applied, recommended-not-applied, and attached-held
    scopes remain distinct downstream; held cells never vote in the
    authoritative reassignment verdict.  REVIEW_UNEXPECTED_IDENTITY remains
    scoreable only as a separate review-only population.  KEEP, unresolved,
    conflicted, ploidy, homotet-occupancy, and standalone cellular-origin
    proposals are outside the identity-swap scoring universe.
    """
    action = clean(source.get("final_action")).upper()
    event_class = clean(source.get("event_class")).upper()
    event_key = _score_pair_event_key(source)
    if (
        action == "REASSIGN_GENOTYPE"
        and event_class in SCORE_PAIR_SUPPORTED_SWAP_EVENT_CLASSES
    ):
        applied = clean(source.get("reassignment_applied")).upper()
        if applied in {"TRUE", "1", "YES", "Y"}:
            return SCORE_PAIR_APPLIED_REASSIGNMENT
        if applied in {"FALSE", "0", "NO", "N"}:
            return SCORE_PAIR_RECOMMENDED_NOT_APPLIED
        raise ValueError(
            f"{clean(source.get('library'))} {clean(source.get('barcode'))}: "
            "reassignment_applied must be an explicit TRUE or FALSE value"
        )
    if action == "REVIEW_CELLULAR_ORIGIN" and event_key in supported_event_keys:
        return SCORE_PAIR_SUPPORTED_EVENT_HELD_CELL
    if action == "REVIEW_UNEXPECTED_IDENTITY":
        return SCORE_PAIR_REVIEW_ONLY
    return ""


def _score_pair_scope_exclusion_reason(
    source: Mapping[str, object],
    supported_event_keys: Set[Tuple[str, str]],
) -> Tuple[str, str]:
    action = clean(source.get("final_action")).upper()
    event_class = clean(source.get("event_class")).upper()
    event_key = _score_pair_event_key(source)
    if (
        action == "REASSIGN_GENOTYPE"
        and event_class not in SCORE_PAIR_SUPPORTED_SWAP_EVENT_CLASSES
    ):
        return (
            "REASSIGNMENT_NOT_A_SUPPORTED_SAMPLE_SWAP_EVENT",
            "cell-level reassignment is outside the established v10 "
            f"population-supported swap classes; event_class={event_class or 'MISSING'}",
        )
    if action == "KEEP":
        return (
            "FINAL_ACTION_KEEP_NOT_A_SUPPORTED_SWAP",
            "the final v10 reconciliation decision retained the original identity",
        )
    if action == "REVIEW_CELLULAR_ORIGIN":
        return (
            "STANDALONE_CELLULAR_ORIGIN_REVIEW_NOT_SUPPORTED_SWAP",
            "no REASSIGN_GENOTYPE cell exists for the same event and exact proposed identity"
            if event_key not in supported_event_keys
            else "cellular-origin review was not admitted by the supported-event gate",
        )
    if action == "REVIEW_HOMOTET_OCCUPANCY":
        return (
            "HOMOTET_OCCUPANCY_REVIEW_NOT_IDENTITY_SWAP",
            "same-donor ploidy/occupancy ambiguity is handled by the NN route",
        )
    if action == "RECLASSIFY_PLOIDY":
        return (
            "PLOIDY_RECLASSIFICATION_NOT_IDENTITY_SWAP",
            "ploidy reclassification is not an original-versus-swap comparison",
        )
    if action == "KEEP_CURRENT_CONFLICTED":
        return (
            "CONFLICTED_PROPOSAL_NOT_SUPPORTED_SWAP",
            "the final reconciliation decision found conflicting evidence",
        )
    if action == "UNRESOLVED_INSUFFICIENT_EVIDENCE":
        return (
            "UNRESOLVED_PROPOSAL_NOT_SUPPORTED_SWAP",
            "the final reconciliation decision found insufficient evidence",
        )
    return (
        "FINAL_ACTION_OUTSIDE_SUPPORTED_SWAP_SCOPE",
        f"final_action={action or 'MISSING'}",
    )


def _score_pair_candidate_row(
    source: Mapping[str, object],
    library: str,
    barcode: str,
    pair_id: str,
    role: str,
    genotype: str,
    original: str,
    proposed: str,
    project_status: str,
    admissibility: str,
    expected_status: str,
    score_population_scope: str,
    supported_event_key: str,
) -> dict:
    ploidy = _score_pair_ploidy(genotype)
    origin = (
        "ORIGINAL_ALLOWED_DEMUX"
        if role == "ORIGINAL_ALLOWED_DEMUX"
        else "RECONCILIATION_NOMINATED_SWAP"
    )
    return {
        "library": library,
        "barcode": barcode,
        "hypothesis_id": pair_id + (
            ":ORIGINAL" if role == "ORIGINAL_ALLOWED_DEMUX" else ":SWAP"
        ),
        "state_notation": biological_state(genotype, ploidy),
        "donor_genotype": genotype,
        "donor_components": ",".join(donor_components(genotype)),
        "biological_ploidy": ploidy,
        "droplet_state": "SINGLE_CELL_CANDIDATE",
        "droplet_constituents": "",
        "candidate_origin": origin,
        "current_state_notation": biological_state(original, _score_pair_ploidy(original)),
        "current_donor_genotype": original,
        "expected_genotype_status": expected_status,
        "project_genotype_status": project_status,
        "biological_admissibility": admissibility,
        "physical_resolution_status": (
            "CONSTRAINED_DEMUX_TOP_WINNER"
            if role == "ORIGINAL_ALLOWED_DEMUX"
            else clean(source.get("proposed_uid_resolution_status"))
            or "RECONCILIATION_NOMINATED"
        ),
        "source_identity": clean(source.get("source_identity")),
        "replaced_component": clean(source.get("replaced_component")),
        "replacement_component": clean(source.get("replacement_component")),
        "preserved_partner": clean(source.get("preserved_partner")),
        "candidate_priority": 1 if role == "ORIGINAL_ALLOWED_DEMUX" else 2,
        "candidate_generation_notes": (
            "FROZEN_ORIGINAL_ALLOWED_DEMUX_TOP_WINNER"
            if role == "ORIGINAL_ALLOWED_DEMUX"
            else "EXACT_RECONCILIATION_PROPOSAL;UNCONSTRAINED_DISCOVERY_NOT_DIRECTLY_SCORED"
        ),
        "schema_version": SCORE_PAIR_SCHEMA,
        "score_pair_id": pair_id,
        "score_pair_role": role,
        "score_pair_source": "POST_RECONCILIATION_DECISION_TABLE",
        "score_population_scope": score_population_scope,
        "supported_event_key": supported_event_key,
        "reconciliation_event_id": clean(source.get("event_id")),
        "reconciliation_event_class": clean(source.get("event_class")),
        "reconciliation_event_confidence": clean(source.get("event_confidence")),
        "reconciliation_final_action": clean(source.get("final_action")),
        "reconciliation_decision_confidence": clean(source.get("decision_confidence")),
        "reconciliation_reassignment_applied": clean(source.get("reassignment_applied")),
        "reconciliation_current_refined_assignment": canonical_genotype(
            source.get("current_refined_assignment", "")
        ),
        "original_demux_assignment": original,
        "proposed_donor_genotype": proposed,
        "score_scope_contract": SCORE_PAIR_CONTRACT,
    }


def score_pairs_main():
    args = score_pairs_parse_args()
    libraries = parse_library_spec(args.libraries)
    decisions_root = Path(args.decisions_root)
    output_root = Path(args.output_root)
    output_root.mkdir(parents=True, exist_ok=True)

    validation_path = Path(args.validation_summary)
    if not validation_path.is_file():
        raise FileNotFoundError(
            f"required validation summary is missing: {validation_path}"
        )
    validation_rows = read_tsv(str(validation_path))
    if not validation_rows or any(
        clean(row.get("status")).upper() != "PASS"
        or int(float(clean(row.get("n_failures")) or -1)) != 0
        for row in validation_rows
    ):
        raise ValueError(
            f"validation summary is not an all-PASS ledger: {validation_path}"
        )
    validation_sha256 = sha256_file(str(validation_path))

    expected_by_library = expected_metadata(args.metadata_root)
    global_lines, global_donors = project_biological_universe(args.metadata_root)
    manifest_fields = FIELDS + SCORE_PAIR_EXTRA_FIELDS
    all_summaries = []
    all_exclusions = []

    for library_number in libraries:
        library = f"lib{library_number}"
        decisions_path = decisions_root / f"{library}.reconciled_cells.tsv.gz"
        if not decisions_path.is_file():
            raise FileNotFoundError(
                f"{library}: required reconciliation decision table is missing: "
                f"{decisions_path}"
            )
        if validation_path.stat().st_mtime_ns < decisions_path.stat().st_mtime_ns:
            raise ValueError(
                f"{library}: validation summary predates the decision table; "
                "rerun validation before score-pair construction"
            )
        rows = read_tsv(str(decisions_path))
        required = {
            "library", "barcode", "original_demux_assignment",
            "proposed_donor_genotype", "proposed_biological_admissibility",
            "proposed_global_biological_status",
            "proposed_library_expected_status", "proposed_droplet_state",
            "event_id", "event_class", "event_confidence", "final_action",
            "decision_confidence", "reassignment_applied",
        }
        if rows:
            missing = sorted(required - set(rows[0]))
            if missing:
                raise ValueError(
                    f"{decisions_path}: required reconciliation columns missing: "
                    + ", ".join(missing)
                )

        expected_here = expected_by_library.get(library, {})
        supported_event_keys = {
            _score_pair_event_key(source)
            for source in rows
            if clean(source.get("final_action")).upper() == "REASSIGN_GENOTYPE"
            and not _score_pair_nonidentity_event_reason(source)[0]
            and clean(source.get("event_class")).upper()
            in SCORE_PAIR_SUPPORTED_SWAP_EVENT_CLASSES
            and clean(source.get("event_id"))
            and canonical_genotype(source.get("proposed_donor_genotype", ""))
            and canonical_genotype(source.get("proposed_donor_genotype", ""))
            != canonical_genotype(source.get("original_demux_assignment", ""))
        }
        manifest_rows = []
        exclusions = []
        summary_reasons = Counter()
        n_not_nominated = 0
        n_considered = 0

        def exclude(source, reason, detail):
            summary_reasons[reason] += 1
            exclusions.append({
                "library": library,
                "barcode": clean(source.get("barcode")),
                "original_demux_assignment": canonical_genotype(
                    source.get("original_demux_assignment", "")
                ),
                "proposed_donor_genotype": canonical_genotype(
                    source.get("proposed_donor_genotype", "")
                ),
                "reconciliation_event_id": clean(source.get("event_id")),
                "reconciliation_event_class": clean(source.get("event_class")),
                "reconciliation_final_action": clean(source.get("final_action")),
                "exclusion_reason": reason,
                "detail": detail,
                "score_scope_contract": SCORE_PAIR_CONTRACT,
                "schema_version": SCORE_PAIR_SCHEMA,
            })

        seen_barcodes = set()
        for source in rows:
            barcode = clean(source.get("barcode"))
            if not barcode:
                raise ValueError(f"{decisions_path}: blank barcode")
            if barcode in seen_barcodes:
                raise ValueError(f"{decisions_path}: duplicate barcode {barcode}")
            seen_barcodes.add(barcode)

            original = canonical_genotype(source.get("original_demux_assignment", ""))
            proposed = canonical_genotype(source.get("proposed_donor_genotype", ""))
            event_id = clean(source.get("event_id"))
            event_class = clean(source.get("event_class"))
            if not event_id or not proposed or proposed == original:
                n_not_nominated += 1
                continue
            n_considered += 1

            nonidentity_reason, nonidentity_detail = (
                _score_pair_nonidentity_event_reason(source)
            )
            if nonidentity_reason:
                exclude(source, nonidentity_reason, nonidentity_detail)
                continue

            score_population_scope = _score_pair_population_scope(
                source, supported_event_keys
            )
            if not score_population_scope:
                reason, detail = _score_pair_scope_exclusion_reason(
                    source, supported_event_keys
                )
                exclude(source, reason, detail)
                continue
            supported_event_key = (
                _score_pair_event_key_text(
                    library, _score_pair_event_key(source)
                )
                if score_population_scope in (
                    SCORE_PAIR_REASSIGNMENT_SCOPES
                    | {SCORE_PAIR_SUPPORTED_EVENT_HELD_CELL}
                )
                else "NA"
            )
            if not original:
                exclude(source, "MISSING_ORIGINAL_ALLOWED_DEMUX", "original assignment is blank")
                continue
            original_status, original_admissibility = candidate_biological_status(
                original, expected_here, global_lines, global_donors
            )
            if original not in expected_here:
                exclude(
                    source,
                    "ORIGINAL_NOT_IN_ALLOWED_LIBRARY_ROSTER",
                    f"original={original}; project_status={original_status}",
                )
                continue
            if original_admissibility not in {
                "SINGLET_IDENTITY_CANDIDATE", "BIOLOGICAL_SINGLE_CELL_ALLOWED"
            }:
                exclude(
                    source,
                    "ORIGINAL_NOT_A_BIOLOGICAL_SINGLE_CELL_IDENTITY",
                    f"original={original}; admissibility={original_admissibility}",
                )
                continue

            proposed_status, proposed_admissibility = candidate_biological_status(
                proposed, expected_here, global_lines, global_donors
            )
            declared_admissibility = clean(source.get("proposed_biological_admissibility"))
            if proposed_admissibility not in {
                "SINGLET_IDENTITY_CANDIDATE", "BIOLOGICAL_SINGLE_CELL_ALLOWED"
            }:
                exclude(
                    source,
                    "PROPOSED_SWAP_NOT_A_REAL_BIOLOGICAL_LINE",
                    f"proposed={proposed}; project_status={proposed_status}",
                )
                continue
            if declared_admissibility and declared_admissibility != proposed_admissibility:
                exclude(
                    source,
                    "RECONCILIATION_METADATA_ADMISSIBILITY_MISMATCH",
                    f"declared={declared_admissibility}; recomputed={proposed_admissibility}",
                )
                continue
            if proposed in expected_here or clean(
                    source.get("proposed_library_expected_status")).upper() == "EXPECTED":
                exclude(
                    source,
                    "PROPOSED_IDENTITY_ALREADY_ALLOWED_IN_LIBRARY",
                    "all-expected roster alternatives are not sample-swap evidence",
                )
                continue
            if "TECHNICAL" in clean(source.get("proposed_droplet_state")).upper():
                exclude(
                    source,
                    "PROPOSED_TECHNICAL_MULTIPLET_NOT_AN_IDENTITY",
                    f"proposed_droplet_state={clean(source.get('proposed_droplet_state'))}",
                )
                continue
            if _score_pair_same_snp_identity(original, proposed):
                exclude(
                    source,
                    "SAME_DONOR_PLOIDY_NOT_NUCLEAR_SCOREABLE",
                    "A versus A+A is preserved for the NN/ploidy decision route",
                )
                continue

            pair_id = f"{library}:{barcode}:RECONCILIATION_SWAP"
            manifest_rows.append(_score_pair_candidate_row(
                source, library, barcode, pair_id,
                "ORIGINAL_ALLOWED_DEMUX", original, original, proposed,
                original_status, original_admissibility, "EXPECTED",
                score_population_scope, supported_event_key,
            ))
            manifest_rows.append(_score_pair_candidate_row(
                source, library, barcode, pair_id,
                "RECONCILIATION_NOMINATED_SWAP", proposed, original, proposed,
                proposed_status, proposed_admissibility,
                "UNEXPECTED_RECONCILIATION_NOMINATION",
                score_population_scope, supported_event_key,
            ))

        manifest_path = (
            output_root / f"{library}.reconciliation_score_pairs.tsv.gz"
        )
        exclusion_path = (
            output_root / f"{library}.reconciliation_score_pair_exclusions.tsv.gz"
        )
        write_tsv(str(manifest_path), manifest_rows, manifest_fields)
        write_tsv(str(exclusion_path), exclusions, SCORE_PAIR_EXCLUSION_FIELDS)
        pair_count = len(manifest_rows) // 2
        population_scope_counts = Counter(
            row["score_population_scope"]
            for row in manifest_rows
            if row["score_pair_role"] == "ORIGINAL_ALLOWED_DEMUX"
        )
        summary = {
            "library": library,
            "source_decisions": str(decisions_path.resolve()),
            "source_decisions_sha256": sha256_file(str(decisions_path)),
            "validation_summary": str(validation_path.resolve()),
            "validation_summary_sha256": validation_sha256,
            "score_pair_builder": str(Path(__file__).resolve()),
            "score_pair_builder_sha256": sha256_file(str(Path(__file__).resolve())),
            "metadata_expected_genotypes_sha256": sha256_file(str(
                Path(args.metadata_root) / "library_expected_genotypes.tsv")),
            "metadata_global_biological_lines_sha256": sha256_file(str(
                Path(args.metadata_root) / "global_biological_lines.tsv")),
            "metadata_global_donors_sha256": sha256_file(str(
                Path(args.metadata_root) / "global_donors.tsv")),
            "n_decision_rows": len(rows),
            "n_not_reconciliation_nominated": n_not_nominated,
            "n_candidate_rows_considered": n_considered,
            "n_supported_event_keys": len(supported_event_keys),
            "n_supported_reassignment_pairs": sum(
                population_scope_counts.get(scope, 0)
                for scope in SCORE_PAIR_REASSIGNMENT_SCOPES
            ),
            "n_supported_event_held_pairs": population_scope_counts.get(
                SCORE_PAIR_SUPPORTED_EVENT_HELD_CELL, 0
            ),
            "n_review_only_pairs": population_scope_counts.get(
                SCORE_PAIR_REVIEW_ONLY, 0
            ),
            "n_score_pairs": pair_count,
            "n_manifest_rows": len(manifest_rows),
            "n_candidate_rows_excluded": len(exclusions),
            "exclusion_reason_counts": ",".join(
                f"{reason}:{count}" for reason, count in sorted(summary_reasons.items())
            ) or "NONE",
            "score_scope_contract": SCORE_PAIR_CONTRACT,
            "schema_version": SCORE_PAIR_SCHEMA,
        }
        write_tsv(
            str(output_root / f"{library}.reconciliation_score_pair_summary.tsv"),
            [summary], SCORE_PAIR_SUMMARY_FIELDS,
        )
        all_summaries.append(summary)
        all_exclusions.extend(exclusions)
        print(
            f"{library}: {pair_count} decision-scoped score pairs "
            f"({population_scope_counts.get(SCORE_PAIR_APPLIED_REASSIGNMENT, 0)} applied, "
            f"{population_scope_counts.get(SCORE_PAIR_RECOMMENDED_NOT_APPLIED, 0)} recommended, "
            f"{population_scope_counts.get(SCORE_PAIR_SUPPORTED_EVENT_HELD_CELL, 0)} attached held, "
            f"{population_scope_counts.get(SCORE_PAIR_REVIEW_ONLY, 0)} review-only); "
            f"{len(exclusions)} candidate rows excluded by the biological/scope gate"
        )

    write_tsv(
        str(output_root / "all_libraries.reconciliation_score_pair_summary.tsv"),
        all_summaries, SCORE_PAIR_SUMMARY_FIELDS,
    )
    write_tsv(
        str(output_root / "all_libraries.reconciliation_score_pair_exclusions.tsv.gz"),
        all_exclusions, SCORE_PAIR_EXCLUSION_FIELDS,
    )
    return 0


# -----------------------------------------------------------------------------
# candidate-axis-pairs
# -----------------------------------------------------------------------------

CANDIDATE_AXIS_PAIR_SCHEMA = "identity_candidate_axis_pair_manifest_v1"
CANDIDATE_AXIS_OPERATIONAL_CONTRACT = (
    "ORIGINAL_ALLOWED_DEMUX_VS_RECONCILIATION_NOMINATED_SWAP_ONLY"
)
CANDIDATE_AXIS_RETAINED_CONTRACT = (
    "ORIGINAL_ALLOWED_DEMUX_VS_FROZEN_SUPPORTED_EVENT_PROPOSAL_CONTRAST_ONLY"
)
CANDIDATE_AXIS_RETAINED = "RETAINED_ORIGINAL_CONTRAST_ONLY"
CANDIDATE_AXIS_POPULATIONS = {
    SCORE_PAIR_APPLIED_REASSIGNMENT,
    SCORE_PAIR_RECOMMENDED_NOT_APPLIED,
    SCORE_PAIR_SUPPORTED_EVENT_HELD_CELL,
    SCORE_PAIR_REVIEW_ONLY,
    CANDIDATE_AXIS_RETAINED,
}
CANDIDATE_AXIS_EXTRA_FIELDS = [
    "score_pair_id", "score_pair_role", "score_pair_source",
    "score_population_scope", "population_votes_in_authoritative_event",
    "supported_event_key", "selected_supported_event_id",
    "selected_supported_event_proposal", "reconciliation_event_id",
    "reconciliation_event_class", "reconciliation_event_confidence",
    "reconciliation_final_action", "reconciliation_decision_confidence",
    "reconciliation_reassignment_applied",
    "reconciliation_current_refined_assignment", "reconciled_donor_genotype",
    "reconciled_droplet_state", "original_demux_assignment",
    "proposed_donor_genotype", "reconciliation_nominated_swap",
    "candidate_b_fixed_identity", "source_reconciliation_event_id",
    "source_reconciliation_proposed_identity", "pair_construction_mode",
    "score_scope_contract",
]
CANDIDATE_AXIS_EXCLUSION_FIELDS = [
    "exclusion_stage", "library", "barcode", "original_demux_assignment",
    "proposed_donor_genotype", "source_reconciliation_event_id",
    "source_reconciliation_proposed_identity", "reconciliation_event_id",
    "reconciliation_event_class", "reconciliation_final_action",
    "exclusion_reason", "detail", "selected_supported_event_id",
    "selected_supported_event_proposal", "score_scope_contract",
    "schema_version",
]
CANDIDATE_AXIS_AUDIT_FIELDS = [
    "status", "check", "value", "detail", "path", "size_bytes",
    "mtime_ns", "schema_version",
]
CANDIDATE_AXIS_SUMMARY_FIELDS = [
    "library", "selected_supported_event_id",
    "selected_supported_event_proposal", "supported_event_key",
    "selected_event_class", "selected_event_confidence",
    "authoritative_original_assignment_strata", "n_decision_rows",
    "n_applied_reassignment_pairs", "n_recommended_not_applied_pairs",
    "n_supported_event_held_pairs", "n_review_only_pairs",
    "n_retained_original_contrast_pairs", "n_score_pairs", "n_manifest_rows",
    "n_pair_construction_exclusions", "n_explicitly_non_nominated_rows",
    "population_scope_counts", "exclusion_reason_counts", "error_ref",
    "error_alt", "min_evidence", "min_evidence_source",
    "poor_fit_residual", "poor_fit_residual_source", "source_decisions",
    "validation_summary", "samples_path", "pileup_sites_path",
    "pileup_observations_path", "frozen_v3_probability_path",
    "frozen_v3_provenance_path", "frozen_v6_3_cell_path",
    "frozen_v6_3_review_only_path", "v6_3_aggregator_path",
    "pair_builder_path", "warnings",
    "score_scope_contract_counts", "schema_version",
]


def candidate_axis_pairs_parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Build one library's fixed A/B candidate-axis manifest from "
            "already-finalized supported reconciliation events. "
            "Retained cells are nonvoting, selection-conditioned contrasts."
        )
    )
    parser.add_argument("--libraries", nargs="+", required=True)
    parser.add_argument("--decisions-root", required=True)
    parser.add_argument("--metadata-root", required=True)
    parser.add_argument("--validation-summary", required=True)
    parser.add_argument("--event-id", default="")
    parser.add_argument("--proposed-identity", default="")
    parser.add_argument("--output-root", required=True)
    parser.add_argument("--samples", default="")
    parser.add_argument("--pileup-sites", default="")
    parser.add_argument("--pileup-observations", default="")
    parser.add_argument("--frozen-v3-probability", default="")
    parser.add_argument("--frozen-v3-provenance", default="")
    parser.add_argument("--frozen-v6-3-cell", default="")
    parser.add_argument("--frozen-v6-3-review-only", default="")
    parser.add_argument(
        "--v6-3-source",
        default=str(Path(__file__).resolve().with_name(
            "identity_probability_aggregate.py"
        )),
    )
    parser.add_argument("--min-evidence", type=int, default=10)
    parser.add_argument(
        "--min-evidence-source",
        default="AUDITED_TETRA_SCORE_CALLS_V3_DEFAULT",
    )
    parser.add_argument("--poor-fit-residual", type=float, default=None)
    return parser.parse_args()


def _candidate_axis_atomic_tsv(path: Path, rows, fields):
    path.parent.mkdir(parents=True, exist_ok=True)
    # identity_reconciliation_common.write_tsv already stages and atomically
    # replaces while preserving gzip based on the final filename.
    write_tsv(str(path), rows, fields)


def _candidate_axis_file_record(path):
    if not path:
        return {"path": "NA", "size": "NA", "mtime": "NA"}
    resolved = Path(path).resolve()
    if not resolved.is_file():
        raise FileNotFoundError(f"required candidate-axis input is missing: {resolved}")
    stat_result = resolved.stat()
    return {
        "path": str(resolved),
        "size": str(stat_result.st_size),
        "mtime": str(stat_result.st_mtime_ns),
    }


def _candidate_axis_optional_file_record(path):
    if not path or not Path(path).is_file():
        return {"path": "NA", "size": "NA", "mtime": "NA"}
    return _candidate_axis_file_record(path)


def _candidate_axis_unanimous_finite(rows, field):
    values = set()
    for row in rows:
        raw = clean(row.get(field))
        if not raw:
            raise ValueError(f"frozen V3 row is missing {field}")
        value = float(raw)
        if not math.isfinite(value):
            raise ValueError(f"frozen V3 {field} is nonfinite")
        values.add(value)
    if len(values) != 1:
        raise ValueError(
            f"frozen V3 {field} must be one unanimous run-level value; saw {sorted(values)}"
        )
    return next(iter(values))


def _candidate_axis_resolve_inputs(args, library):
    decisions_root = Path(args.decisions_root).resolve()
    identity_root = decisions_root.parent
    nuclear_root = identity_root / "nuclear"
    v3_probability = Path(args.frozen_v3_probability).resolve() if (
        args.frozen_v3_probability
    ) else nuclear_root / f"{library}.identity_pair_probabilities.tsv.gz"
    v3_provenance = Path(args.frozen_v3_provenance).resolve() if (
        args.frozen_v3_provenance
    ) else nuclear_root / f"{library}.identity_pair_probability_provenance.tsv"
    provenance = {}
    if v3_provenance.is_file():
        provenance_rows = read_tsv(str(v3_provenance))
        if len(provenance_rows) == 1:
            provenance = provenance_rows[0]
    samples = args.samples or clean(provenance.get("samples_path"))
    sites = args.pileup_sites or clean(provenance.get("pileup_sites_path"))
    observations = args.pileup_observations or clean(
        provenance.get("pileup_observations_path")
    )
    if not samples or not sites or not observations:
        raise ValueError(
            "samples, pileup sites, and pileup observations must be supplied "
            "explicitly or recorded by frozen V3 provenance"
        )
    return {
        "v3_probability": str(v3_probability),
        "v3_provenance": str(v3_provenance) if v3_provenance.is_file() else "",
        "v3_provenance_row": provenance,
        "samples": str(Path(samples).resolve()),
        "sites": str(Path(sites).resolve()),
        "observations": str(Path(observations).resolve()),
    }


def _candidate_axis_bool(value):
    normalized = clean(value).upper()
    if normalized in {"TRUE", "1", "YES", "Y"}:
        return True
    if normalized in {"FALSE", "0", "NO", "N"}:
        return False
    raise ValueError(f"expected an explicit TRUE/FALSE value, saw {value!r}")


def _candidate_axis_source_audit_row(status, check, value, detail="", record=None):
    record = record or {"path": "NA", "size": "NA", "mtime": "NA"}
    return {
        "status": status,
        "check": check,
        "value": value if value not in {None, ""} else "NA",
        "detail": detail or "NA",
        "path": record["path"],
        "size_bytes": record["size"],
        "mtime_ns": record["mtime"],
        "schema_version": "identity_candidate_axis_pair_source_audit_v1",
    }


def _candidate_axis_manifest_row(
    source,
    event_metadata,
    library,
    barcode,
    pair_id,
    role,
    genotype,
    original,
    fixed_proposal,
    project_status,
    admissibility,
    expected_status,
    population,
    retained,
    event_key,
):
    ploidy = _score_pair_ploidy(genotype)
    candidate_origin = (
        "ORIGINAL_ALLOWED_DEMUX"
        if role == "ORIGINAL_ALLOWED_DEMUX"
        else (
            "FROZEN_SUPPORTED_EVENT_PROPOSAL_CONTRAST"
            if retained
            else "RECONCILIATION_NOMINATED_SWAP"
        )
    )
    pair_source = (
        "RETAINED_KEEP_CELL_VS_FROZEN_SUPPORTED_EVENT_PROPOSAL"
        if retained
        else "POST_RECONCILIATION_DECISION_TABLE"
    )
    contract = (
        CANDIDATE_AXIS_RETAINED_CONTRACT
        if retained
        else CANDIDATE_AXIS_OPERATIONAL_CONTRACT
    )
    source_event = clean(source.get("event_id"))
    source_proposal = canonical_genotype(source.get("proposed_donor_genotype", ""))
    return {
        "library": library,
        "barcode": barcode,
        "hypothesis_id": pair_id + (
            ":ORIGINAL" if role == "ORIGINAL_ALLOWED_DEMUX" else ":FIXED_B"
        ),
        "state_notation": biological_state(genotype, ploidy),
        "donor_genotype": genotype,
        "donor_components": ",".join(donor_components(genotype)),
        "biological_ploidy": ploidy,
        "droplet_state": "SINGLE_CELL_CANDIDATE",
        "droplet_constituents": "",
        "candidate_origin": candidate_origin,
        "current_state_notation": biological_state(
            original, _score_pair_ploidy(original)
        ),
        "current_donor_genotype": original,
        "expected_genotype_status": expected_status,
        "project_genotype_status": project_status,
        "biological_admissibility": admissibility,
        "physical_resolution_status": (
            "CONSTRAINED_DEMUX_TOP_WINNER"
            if role == "ORIGINAL_ALLOWED_DEMUX"
            else (
                "FROZEN_SUPPORTED_EVENT_CONTRAST"
                if retained
                else clean(source.get("proposed_uid_resolution_status"))
                or "RECONCILIATION_NOMINATED"
            )
        ),
        "source_identity": clean(source.get("source_identity")),
        "replaced_component": clean(source.get("replaced_component")),
        "replacement_component": clean(source.get("replacement_component")),
        "preserved_partner": clean(source.get("preserved_partner")),
        "candidate_priority": 1 if role == "ORIGINAL_ALLOWED_DEMUX" else 2,
        "candidate_generation_notes": (
            "FROZEN_ORIGINAL_ALLOWED_DEMUX_TOP_WINNER"
            if role == "ORIGINAL_ALLOWED_DEMUX"
            else (
                "FROZEN_SUPPORTED_EVENT_PROPOSAL_CONTRAST_ONLY_SELECTION_CONDITIONED_NONVOTING"
                if retained
                else "EXACT_RECONCILIATION_PROPOSAL"
            )
        ),
        "schema_version": CANDIDATE_AXIS_PAIR_SCHEMA,
        "score_pair_id": pair_id,
        "score_pair_role": role,
        "score_pair_source": pair_source,
        "score_population_scope": population,
        "population_votes_in_authoritative_event": (
            "TRUE" if population in SCORE_PAIR_REASSIGNMENT_SCOPES else "FALSE"
        ),
        "supported_event_key": event_key,
        "selected_supported_event_id": clean(event_metadata.get("event_id")),
        "selected_supported_event_proposal": fixed_proposal,
        "reconciliation_event_id": clean(event_metadata.get("event_id")),
        "reconciliation_event_class": clean(event_metadata.get("event_class")),
        "reconciliation_event_confidence": clean(
            event_metadata.get("event_confidence")
        ),
        "reconciliation_final_action": clean(source.get("final_action")),
        "reconciliation_decision_confidence": clean(
            source.get("decision_confidence")
        ),
        "reconciliation_reassignment_applied": clean(
            source.get("reassignment_applied")
        ),
        "reconciliation_current_refined_assignment": canonical_genotype(
            source.get("current_refined_assignment", "")
        ),
        "reconciled_donor_genotype": canonical_genotype(
            source.get("reconciled_donor_genotype", "")
        ),
        "reconciled_droplet_state": clean(source.get("reconciled_droplet_state")),
        "original_demux_assignment": original,
        "proposed_donor_genotype": fixed_proposal,
        "reconciliation_nominated_swap": "" if retained else fixed_proposal,
        "candidate_b_fixed_identity": fixed_proposal,
        "source_reconciliation_event_id": source_event,
        "source_reconciliation_proposed_identity": source_proposal,
        "pair_construction_mode": (
            "SUPPORTED_EVENT_CHALLENGE"
            if retained
            else "RECONCILIATION_NOMINATED_SWAP"
        ),
        "score_scope_contract": contract,
    }


def candidate_axis_pairs_main():
    args = candidate_axis_pairs_parse_args()
    libraries = parse_library_spec(args.libraries)
    if len(libraries) != 1:
        raise ValueError("candidate-axis pair construction requires exactly one library")
    library = f"lib{libraries[0]}"
    targeted_event_id = clean(args.event_id)
    targeted_proposal = canonical_genotype(args.proposed_identity)
    if bool(targeted_event_id) != bool(targeted_proposal):
        raise ValueError(
            "candidate-axis event ID and proposal must be supplied together"
        )
    targeted_mode = bool(targeted_event_id)
    output_root = Path(args.output_root).resolve()
    output_root.mkdir(parents=True, exist_ok=True)
    audit_path = output_root / f"{library}.candidate_axis_pair_source_audit.tsv"
    audit_rows = []
    warnings = []

    def audit(check, value, detail="", record=None, status="PASS"):
        audit_rows.append(_candidate_axis_source_audit_row(
            status, check, value, detail, record
        ))

    try:
        decisions_path = (
            Path(args.decisions_root).resolve()
            / f"{library}.reconciled_cells.tsv.gz"
        )
        validation_path = Path(args.validation_summary).resolve()
        metadata_root = Path(args.metadata_root).resolve()
        metadata_paths = {
            "expected_genotypes": metadata_root / "library_expected_genotypes.tsv",
            "global_biological_lines": metadata_root / "global_biological_lines.tsv",
            "global_donors": metadata_root / "global_donors.tsv",
        }
        input_paths = _candidate_axis_resolve_inputs(args, library)
        records = {
            "decisions": _candidate_axis_file_record(decisions_path),
            "validation": _candidate_axis_file_record(validation_path),
            "samples": _candidate_axis_file_record(input_paths["samples"]),
            "pileup_sites": _candidate_axis_file_record(input_paths["sites"]),
            "pileup_observations": _candidate_axis_file_record(
                input_paths["observations"]
            ),
            "frozen_v3_provenance": _candidate_axis_optional_file_record(
                input_paths["v3_provenance"]
            ),
        }
        for name, path in metadata_paths.items():
            records[name] = _candidate_axis_file_record(path)
        records["frozen_v6_3_cell"] = _candidate_axis_optional_file_record(
            args.frozen_v6_3_cell
        )
        records["frozen_v6_3_review_only"] = _candidate_axis_optional_file_record(
            args.frozen_v6_3_review_only
        )
        for name, record in records.items():
            audit(
                f"INPUT_{name.upper()}",
                "PRESENT" if record["path"] != "NA" else "OPTIONAL_ABSENT",
                record=record,
            )

        validation_rows = read_tsv(str(validation_path))
        if not validation_rows:
            raise ValueError("validation summary is empty")
        for row in validation_rows:
            if clean(row.get("status")).upper() != "PASS" or int(float(
                clean(row.get("n_failures")) or -1
            )) != 0:
                raise ValueError("validation summary is not an all-PASS zero-failure ledger")
        if validation_path.stat().st_mtime_ns < decisions_path.stat().st_mtime_ns:
            raise ValueError("validation summary predates the selected decision table")
        audit("VALIDATION_LEDGER", "ALL_PASS_ZERO_FAILURES")

        decisions = read_tsv(str(decisions_path))
        if not decisions:
            raise ValueError("selected decision table is empty")
        seen_barcodes = set()
        for row in decisions:
            barcode = clean(row.get("barcode"))
            if not barcode or barcode in seen_barcodes:
                raise ValueError(
                    f"decision table has blank or duplicate barcode: {barcode!r}"
                )
            seen_barcodes.add(barcode)
            if clean(row.get("library")) not in {library, str(libraries[0])}:
                raise ValueError(
                    f"decision table contains a cross-library row: {row.get('library')!r}"
                )

        authoritative_by_event = defaultdict(list)
        all_supported_by_stratum = defaultdict(set)
        for row in decisions:
            event_key = _score_pair_event_key(row)
            original = canonical_genotype(
                row.get("original_demux_assignment", "")
            )
            nonidentity_reason, _ = _score_pair_nonidentity_event_reason(row)
            if (
                clean(row.get("final_action")).upper() == "REASSIGN_GENOTYPE"
                and clean(row.get("event_class")).upper()
                in SCORE_PAIR_SUPPORTED_SWAP_EVENT_CLASSES
                and event_key[0] and event_key[1] and original
                and event_key[1] != original
                and not nonidentity_reason
            ):
                authoritative_by_event[event_key].append(row)
                all_supported_by_stratum[original].add(event_key)

        supported_event_keys = set(authoritative_by_event)
        requested_key = (targeted_event_id, targeted_proposal)
        if targeted_mode:
            if requested_key not in supported_event_keys:
                raise ValueError(
                    "targeted candidate-axis event is not an already-finalized "
                    f"supported event: {library}|{targeted_event_id}|{targeted_proposal}"
                )
            selected_event_keys = {requested_key}
        else:
            selected_event_keys = set(supported_event_keys)

        event_metadata = {}
        source_strata = {}
        metadata_keys = (
            "event_id", "event_class", "event_confidence",
        )
        for selected_key in sorted(selected_event_keys):
            rows = authoritative_by_event[selected_key]
            for field in metadata_keys:
                values = {clean(row.get(field)) for row in rows}
                if len(values) != 1 or not next(iter(values)):
                    raise ValueError(
                        "finalized event has inconsistent or blank "
                        f"{field}: {library}|{selected_key[0]}|{selected_key[1]} "
                        f"values={sorted(values)}"
                    )
            event_metadata[selected_key] = rows[0]
            strata = {
                canonical_genotype(row.get("original_demux_assignment", ""))
                for row in rows
            }
            if "" in strata or not strata:
                raise ValueError(
                    "finalized event has blank/no authoritative source strata: "
                    f"{library}|{selected_key[0]}|{selected_key[1]}"
                )
            source_strata[selected_key] = strata
        audit(
            "SELECTED_FINALIZED_EVENT_KEYS",
            ";".join(
                _score_pair_event_key_text(library, key)
                for key in sorted(selected_event_keys)
            ) or "NONE",
            "TARGETED_DIAGNOSTIC" if targeted_mode else "ROUTINE_AUTOMATIC",
        )
        audit(
            "AUTHORITATIVE_SOURCE_STRATA",
            ";".join(
                f"{_score_pair_event_key_text(library, key)}="
                + ",".join(sorted(source_strata[key]))
                for key in sorted(selected_event_keys)
            ) or "NONE",
        )

        records["frozen_v3_probability"] = _candidate_axis_optional_file_record(
            input_paths["v3_probability"]
        )
        records["v6_3_source"] = (
            _candidate_axis_file_record(args.v6_3_source)
            if selected_event_keys else
            _candidate_axis_optional_file_record(args.v6_3_source)
        )
        for name in ("frozen_v3_probability", "v6_3_source"):
            record = records[name]
            if record["path"] != "NA":
                input_status = "PRESENT"
            elif selected_event_keys and name == "frozen_v3_probability":
                input_status = "OPTIONAL_ABSENT"
            else:
                input_status = "NOT_APPLICABLE_ZERO_EVENT"
            audit(
                f"INPUT_{name.upper()}",
                input_status,
                record=record,
            )

        expected_by_library = expected_metadata(str(metadata_root))
        global_lines, global_donors = project_biological_universe(str(metadata_root))
        expected_here = expected_by_library.get(library, {})

        if (selected_event_keys and
                records["frozen_v3_probability"]["path"] != "NA"):
            v3_rows = read_tsv(input_paths["v3_probability"])
            if not v3_rows:
                raise ValueError("frozen V3 probability table is empty")
            error_ref = _candidate_axis_unanimous_finite(v3_rows, "error_ref")
            error_alt = _candidate_axis_unanimous_finite(v3_rows, "error_alt")
            if not (0 <= error_ref <= 1 and 0 <= error_alt <= 1 and
                    error_ref + error_alt < 1):
                raise ValueError("frozen V3 error pair is outside the valid transform domain")
            audit("ERROR_REF", repr(error_ref), "unanimous across finite frozen V3 rows")
            audit("ERROR_ALT", repr(error_alt), "unanimous across finite frozen V3 rows")
        elif selected_event_keys:
            error_ref = error_alt = 0.001
            audit(
                "ERROR_REF", repr(error_ref),
                "frozen tetra_score_calls default; optional V3 table absent",
            )
            audit(
                "ERROR_ALT", repr(error_alt),
                "frozen tetra_score_calls default; optional V3 table absent",
            )
        else:
            error_ref = error_alt = 0.001
            audit("ERROR_REF", "NOT_APPLICABLE_ZERO_EVENT")
            audit("ERROR_ALT", "NOT_APPLICABLE_ZERO_EVENT")
        provenance = input_paths["v3_provenance_row"]
        if "demux_nomito" not in records["pileup_sites"]["path"].lower():
            raise ValueError(
                "pileup-site path does not prove the production demux_nomito nuclear source"
            )
        audit(
            "NUCLEAR_EVIDENCE_IDENTITY", "EXPLICIT_DEMUX_NOMITO_INPUTS",
            "primary mitochondrial exclusion is independently enforced by the scorer",
        )
        poor_fit = args.poor_fit_residual
        poor_fit_source = "EXPLICIT_CLI"
        if poor_fit is None:
            raw = clean(provenance.get("poor_fit_residual"))
            if raw:
                poor_fit = float(raw)
                poor_fit_source = "FROZEN_V3_PROVENANCE"
            else:
                poor_fit = 0.30
                poor_fit_source = "DEFAULT_0.30"
        if not math.isfinite(poor_fit) or not 0 <= poor_fit <= 1:
            raise ValueError("poor_fit_residual must be finite within [0,1]")
        if args.min_evidence < 0:
            raise ValueError("min_evidence must be nonnegative")
        if not selected_event_keys:
            poor_fit_source = "NOT_APPLICABLE_ZERO_EVENT"
        audit("MIN_EVIDENCE", str(args.min_evidence), args.min_evidence_source)
        audit("POOR_FIT_RESIDUAL", repr(poor_fit), poor_fit_source)

        samples = [line.strip() for line in Path(records["samples"]["path"]).read_text().splitlines() if line.strip()]
        if not samples or len(samples) != len(set(samples)):
            raise ValueError("nuclear sample vector is empty or contains duplicates")
        sample_set = set(samples)
        for _, selected_proposal in sorted(selected_event_keys):
            for component in donor_components(selected_proposal):
                if component not in sample_set:
                    raise ValueError(
                        "selected proposal component is absent from nuclear "
                        f"sample vector: {component}"
                    )
        audit("NUCLEAR_SAMPLE_VECTOR", "UNIQUE_NONBLANK")

        molecule_path = clean(provenance.get("pileup_molecules_path"))
        molecule_basis = Counter()
        if molecule_path and Path(molecule_path).is_file():
            import gzip
            with gzip.open(molecule_path, "rt") as handle:
                for index, line in enumerate(handle):
                    fields = line.rstrip("\n\r").split("\t")
                    if len(fields) >= 3:
                        molecule_basis[fields[2].strip() or "BLANK"] += 1
                    if index >= 99999:
                        break
            audit(
                "MOLECULE_SIDECAR_SCHEMA_INSPECTION",
                ",".join(f"{key}:{value}" for key, value in sorted(molecule_basis.items())) or "EMPTY",
                "inspection only; molecule data never changes the fixed site-balanced primary axis",
            )
        else:
            audit("MOLECULE_SIDECAR_SCHEMA_INSPECTION", "ABSENT")

        manifest_rows = []
        exclusions = []
        exclusion_counts = Counter()
        population_counts = Counter()
        explicitly_non_nominated = 0

        def exclude(
            source, reason, detail, contract=CANDIDATE_AXIS_OPERATIONAL_CONTRACT,
            selected_key=None,
        ):
            if selected_key is None:
                selected_key = (
                    next(iter(selected_event_keys))
                    if len(selected_event_keys) == 1
                    else ("MULTIPLE", "MULTIPLE")
                )
            selected_metadata = event_metadata.get(selected_key, {})
            exclusion_counts[reason] += 1
            exclusions.append({
                "exclusion_stage": "PAIR_CONSTRUCTION",
                "library": library,
                "barcode": clean(source.get("barcode")),
                "original_demux_assignment": canonical_genotype(
                    source.get("original_demux_assignment", "")
                ),
                "proposed_donor_genotype": canonical_genotype(
                    source.get("proposed_donor_genotype", "")
                ),
                "source_reconciliation_event_id": clean(source.get("event_id")),
                "source_reconciliation_proposed_identity": canonical_genotype(
                    source.get("proposed_donor_genotype", "")
                ),
                "reconciliation_event_id": selected_key[0],
                "reconciliation_event_class": (
                    clean(selected_metadata.get("event_class")) or "MULTIPLE"
                ),
                "reconciliation_final_action": clean(source.get("final_action")),
                "exclusion_reason": reason,
                "detail": detail,
                "selected_supported_event_id": selected_key[0],
                "selected_supported_event_proposal": selected_key[1],
                "score_scope_contract": contract,
                "schema_version": CANDIDATE_AXIS_PAIR_SCHEMA,
            })

        admitted_barcodes = set()

        def admit(source, population, retained, selected_key):
            barcode = clean(source.get("barcode"))
            if barcode in admitted_barcodes:
                raise ValueError(
                    f"barcode admitted to more than one candidate-axis population: {barcode}"
                )
            original = canonical_genotype(source.get("original_demux_assignment", ""))
            event_id, proposal = selected_key
            metadata = event_metadata[selected_key]
            event_key = _score_pair_event_key_text(library, selected_key)
            original_status, original_admissibility = candidate_biological_status(
                original, expected_here, global_lines, global_donors
            )
            proposed_status, proposed_admissibility = candidate_biological_status(
                proposal, expected_here, global_lines, global_donors
            )
            if not original:
                exclude(source, "MISSING_ORIGINAL_ALLOWED_DEMUX", "original assignment is blank", selected_key=selected_key)
                return
            if original not in expected_here:
                exclude(
                    source, "ORIGINAL_NOT_IN_ALLOWED_LIBRARY_ROSTER",
                    f"original={original}; project_status={original_status}", selected_key=selected_key,
                )
                return
            allowed_admissibility = {
                "SINGLET_IDENTITY_CANDIDATE", "BIOLOGICAL_SINGLE_CELL_ALLOWED"
            }
            if original_admissibility not in allowed_admissibility:
                exclude(
                    source, "ORIGINAL_NOT_A_BIOLOGICAL_SINGLE_CELL_IDENTITY",
                    f"admissibility={original_admissibility}", selected_key=selected_key,
                )
                return
            if proposed_admissibility not in allowed_admissibility:
                exclude(
                    source, "PROPOSED_SWAP_NOT_A_REAL_BIOLOGICAL_LINE",
                    f"project_status={proposed_status}; admissibility={proposed_admissibility}", selected_key=selected_key,
                )
                return
            if proposal in expected_here:
                exclude(
                    source, "PROPOSED_IDENTITY_ALREADY_ALLOWED_IN_LIBRARY",
                    "fixed proposal must remain library-unexpected", selected_key=selected_key,
                )
                return
            if original == proposal:
                exclude(source, "IDENTICAL_FIXED_CANDIDATES", "A and B are identical", selected_key=selected_key)
                return
            if _score_pair_same_snp_identity(original, proposal):
                exclude(
                    source, "SAME_DONOR_PLOIDY_NOT_NUCLEAR_SCOREABLE",
                    "A versus A+A is not admitted to the nuclear candidate axis", selected_key=selected_key,
                )
                return
            missing_components = [
                component
                for genotype in (original, proposal)
                for component in donor_components(genotype)
                if component not in sample_set
            ]
            if missing_components:
                exclude(
                    source, "CANDIDATE_COMPONENT_MISSING_FROM_NUCLEAR_SAMPLE_VECTOR",
                    ",".join(sorted(set(missing_components))), selected_key=selected_key,
                )
                return
            pair_id = f"{library}:{barcode}:CANDIDATE_AXIS_PILOT"
            b_role = (
                "FROZEN_SUPPORTED_EVENT_PROPOSAL_CONTRAST"
                if retained
                else "RECONCILIATION_NOMINATED_SWAP"
            )
            rows = [
                _candidate_axis_manifest_row(
                    source, metadata, library, barcode, pair_id,
                    "ORIGINAL_ALLOWED_DEMUX", original, original, proposal,
                    original_status, original_admissibility, "EXPECTED",
                    population, retained, event_key,
                ),
                _candidate_axis_manifest_row(
                    source, metadata, library, barcode, pair_id,
                    b_role, proposal, original, proposal, proposed_status,
                    proposed_admissibility,
                    "UNEXPECTED_RECONCILIATION_NOMINATION" if not retained
                    else "FROZEN_SUPPORTED_EVENT_PROPOSAL_CONTRAST",
                    population, retained, event_key,
                ),
            ]
            shared_fields = [
                field for field in CANDIDATE_AXIS_EXTRA_FIELDS
                if field not in {"score_pair_role"}
            ]
            for field in shared_fields:
                if rows[0].get(field) != rows[1].get(field):
                    raise ValueError(
                        f"internal candidate-axis pair metadata mismatch: {barcode} {field}"
                    )
            manifest_rows.extend(rows)
            population_counts[population] += 1
            admitted_barcodes.add(barcode)

        for source in decisions:
            row_event = clean(source.get("event_id"))
            row_proposal = canonical_genotype(source.get("proposed_donor_genotype", ""))
            action = clean(source.get("final_action")).upper()
            row_key = (row_event, row_proposal)
            if row_key in selected_event_keys:
                population = _score_pair_population_scope(source, supported_event_keys)
                if population in {
                    SCORE_PAIR_APPLIED_REASSIGNMENT,
                    SCORE_PAIR_RECOMMENDED_NOT_APPLIED,
                    SCORE_PAIR_SUPPORTED_EVENT_HELD_CELL,
                    SCORE_PAIR_REVIEW_ONLY,
                }:
                    admit(source, population, retained=False, selected_key=row_key)
                    continue
                if action != "KEEP":
                    reason, detail = _score_pair_scope_exclusion_reason(
                        source, supported_event_keys
                    )
                    exclude(source, reason, detail, selected_key=row_key)
                    continue

            original = canonical_genotype(source.get("original_demux_assignment", ""))
            if action != "KEEP":
                explicitly_non_nominated += 1
                continue
            retained_contract = CANDIDATE_AXIS_RETAINED_CONTRACT
            selected_key = None
            if row_event or row_proposal:
                if row_key in selected_event_keys and original in source_strata[row_key]:
                    selected_key = row_key
                elif any(
                    original in source_strata[key] for key in selected_event_keys
                ):
                    exclude(
                        source, "RETAINED_CELL_HAS_DIFFERENT_NOMINATED_EVENT",
                        f"source_event={row_event or 'NA'}; source_proposal={row_proposal or 'NA'}",
                        retained_contract, row_key if row_key in selected_event_keys else None,
                    )
                    continue
                else:
                    explicitly_non_nominated += 1
                    continue
            else:
                mappings = all_supported_by_stratum.get(original, set())
                selected_mappings = mappings & selected_event_keys
                if len(mappings) > 1 and selected_mappings:
                    exclude(
                        source, "AMBIGUOUS_RETAINED_CONTRAST_EVENT_MAPPING",
                        ";".join(
                            f"{eid}|{candidate}"
                            for eid, candidate in sorted(mappings)
                        ),
                        retained_contract,
                    )
                    continue
                if len(mappings) == 1 and len(selected_mappings) == 1:
                    selected_key = next(iter(selected_mappings))
                else:
                    explicitly_non_nominated += 1
                    continue
            try:
                applied = _candidate_axis_bool(source.get("reassignment_applied"))
            except ValueError:
                exclude(
                    source, "RETAINED_REASSIGNMENT_APPLIED_NOT_EXPLICIT_FALSE",
                    f"value={clean(source.get('reassignment_applied')) or 'MISSING'}",
                    retained_contract, selected_key,
                )
                continue
            if applied:
                exclude(
                    source, "RETAINED_REASSIGNMENT_APPLIED_NOT_FALSE",
                    "final_action KEEP cannot be a retained contrast when reassignment_applied is true",
                    retained_contract, selected_key,
                )
                continue
            reconciled = canonical_genotype(source.get("reconciled_donor_genotype", ""))
            if reconciled != original:
                exclude(
                    source, "RETAINED_RECONCILED_IDENTITY_DIFFERS_FROM_ORIGINAL",
                    f"reconciled={reconciled or 'NA'}; original={original}",
                    retained_contract, selected_key,
                )
                continue
            droplet = clean(source.get("reconciled_droplet_state")).upper()
            if droplet != "SINGLE_CELL":
                exclude(
                    source, "RETAINED_NOT_RECONCILED_SINGLE_CELL",
                    f"reconciled_droplet_state={droplet or 'MISSING'}",
                    retained_contract, selected_key,
                )
                continue
            admit(
                source, CANDIDATE_AXIS_RETAINED, retained=True,
                selected_key=selected_key,
            )

        manifest_fields = FIELDS + [
            field for field in CANDIDATE_AXIS_EXTRA_FIELDS if field not in FIELDS
        ]
        if len(manifest_rows) % 2:
            raise ValueError("candidate-axis manifest row count is not divisible by two")
        pair_ids = Counter(row["score_pair_id"] for row in manifest_rows)
        if any(count != 2 for count in pair_ids.values()):
            raise ValueError("candidate-axis manifest contains a pair with other than two rows")
        retained_count = population_counts[CANDIDATE_AXIS_RETAINED]
        if retained_count == 0:
            warnings.append("NO_UNAMBIGUOUS_RETAINED_CONTRASTS")
        if not selected_event_keys:
            warnings.append("NO_FINALIZED_SUPPORTED_EVENTS")

        manifest_path = output_root / f"{library}.candidate_axis_pairs.tsv.gz"
        exclusion_path = output_root / f"{library}.candidate_axis_pair_exclusions.tsv.gz"
        summary_path = output_root / f"{library}.candidate_axis_pair_summary.tsv"
        _candidate_axis_atomic_tsv(manifest_path, manifest_rows, manifest_fields)
        _candidate_axis_atomic_tsv(
            exclusion_path, exclusions, CANDIDATE_AXIS_EXCLUSION_FIELDS
        )
        if len(selected_event_keys) == 1:
            summary_key = next(iter(selected_event_keys))
            summary_metadata = event_metadata[summary_key]
            summary_event_id, summary_proposal = summary_key
            summary_event_key = _score_pair_event_key_text(library, summary_key)
            summary_event_class = clean(summary_metadata.get("event_class"))
            summary_event_confidence = clean(summary_metadata.get("event_confidence"))
            summary_strata = ",".join(sorted(source_strata[summary_key]))
        elif selected_event_keys:
            summary_event_id = summary_proposal = summary_event_key = "MULTIPLE"
            summary_event_class = summary_event_confidence = summary_strata = "MULTIPLE"
        else:
            summary_event_id = summary_proposal = summary_event_key = "NONE"
            summary_event_class = summary_event_confidence = summary_strata = "NONE"
        summary = {
            "library": library,
            "selected_supported_event_id": summary_event_id,
            "selected_supported_event_proposal": summary_proposal,
            "supported_event_key": summary_event_key,
            "selected_event_class": summary_event_class,
            "selected_event_confidence": summary_event_confidence,
            "authoritative_original_assignment_strata": summary_strata,
            "n_decision_rows": len(decisions),
            "n_applied_reassignment_pairs": population_counts[SCORE_PAIR_APPLIED_REASSIGNMENT],
            "n_recommended_not_applied_pairs": population_counts[SCORE_PAIR_RECOMMENDED_NOT_APPLIED],
            "n_supported_event_held_pairs": population_counts[SCORE_PAIR_SUPPORTED_EVENT_HELD_CELL],
            "n_review_only_pairs": population_counts[SCORE_PAIR_REVIEW_ONLY],
            "n_retained_original_contrast_pairs": retained_count,
            "n_score_pairs": len(pair_ids),
            "n_manifest_rows": len(manifest_rows),
            "n_pair_construction_exclusions": len(exclusions),
            "n_explicitly_non_nominated_rows": explicitly_non_nominated,
            "population_scope_counts": ",".join(
                f"{key}:{value}" for key, value in sorted(population_counts.items())
            ) or "NONE",
            "exclusion_reason_counts": ",".join(
                f"{key}:{value}" for key, value in sorted(exclusion_counts.items())
            ) or "NONE",
            "error_ref": repr(error_ref) if selected_event_keys else "NA",
            "error_alt": repr(error_alt) if selected_event_keys else "NA",
            "min_evidence": args.min_evidence,
            "min_evidence_source": (
                args.min_evidence_source if selected_event_keys
                else "NOT_APPLICABLE_ZERO_EVENT"
            ),
            "poor_fit_residual": repr(poor_fit),
            "poor_fit_residual_source": poor_fit_source,
            "source_decisions": records["decisions"]["path"],
            "validation_summary": records["validation"]["path"],
            "samples_path": records["samples"]["path"],
            "pileup_sites_path": records["pileup_sites"]["path"],
            "pileup_observations_path": records["pileup_observations"]["path"],
            "frozen_v3_probability_path": records["frozen_v3_probability"]["path"],
            "frozen_v3_provenance_path": records["frozen_v3_provenance"]["path"],
            "frozen_v6_3_cell_path": records["frozen_v6_3_cell"]["path"],
            "frozen_v6_3_review_only_path": records["frozen_v6_3_review_only"]["path"],
            "v6_3_aggregator_path": records["v6_3_source"]["path"],
            "pair_builder_path": str(Path(__file__).resolve()),
            "warnings": ";".join(warnings) or "NONE",
            "score_scope_contract_counts": ",".join(
                f"{key}:{value}" for key, value in sorted(Counter(
                    row["score_scope_contract"] for row in manifest_rows
                    if row["score_pair_role"] == "ORIGINAL_ALLOWED_DEMUX"
                ).items())
            ) or "NONE",
            "schema_version": CANDIDATE_AXIS_PAIR_SCHEMA,
        }
        _candidate_axis_atomic_tsv(
            summary_path, [summary], CANDIDATE_AXIS_SUMMARY_FIELDS
        )
        audit("PAIR_COUNT", str(len(pair_ids)))
        audit("MANIFEST_ROW_COUNT", str(len(manifest_rows)))
        audit("PAIR_CONSTRUCTION_EXCLUSION_COUNT", str(len(exclusions)))
        audit("POPULATION_COUNTS", summary["population_scope_counts"])
        audit("PAIR_SOURCE_GATE", "PASS")
        _candidate_axis_atomic_tsv(
            audit_path, audit_rows, CANDIDATE_AXIS_AUDIT_FIELDS
        )
        print(
            f"{library}: candidate-axis finalized events={len(selected_event_keys)}; "
            f"pairs={len(pair_ids)} retained={retained_count} "
            f"exclusions={len(exclusions)}"
        )
        return 0
    except Exception as exc:
        audit("PAIR_SOURCE_GATE", "STOP", str(exc), status="STOP")
        _candidate_axis_atomic_tsv(
            audit_path, audit_rows, CANDIDATE_AXIS_AUDIT_FIELDS
        )
        raise


# -----------------------------------------------------------------------------
# finalize
# -----------------------------------------------------------------------------

FINAL_SCHEMA_VERSION = "identity_reconciliation_final_v2_phase3_dispositions"
FINAL_REVIEW_DISPOSITIONS = {
    "ACCEPT_PROPOSAL", "KEEP_CURRENT", "LEAVE_UNRESOLVED",
}
FINAL_NUCLEAR_PROPOSAL_PATTERNS = {
    "ALL_PROPOSAL", "PROPOSAL_PLURALITY_WITH_CURRENT_EXCEPTIONS",
}
FINAL_NUCLEAR_CURRENT_PATTERNS = {
    "ALL_CURRENT", "CURRENT_PLURALITY_WITH_PROPOSAL_EXCEPTIONS",
}
FINAL_PHASE3_MODIFIER_ORDER = (
    "AMBIENT_REFIT_SENSITIVE",
    "AMBIENT_REFIT_UNAVAILABLE",
    "NUCLEAR_CELL_HETEROGENEITY",
    "SAMPLING_ADJUSTED_FIT_DISAGREEMENT",
    "NUCLEAR_NONDISCRIMINATING",
    "NUCLEAR_UNSTABLE",
    "PAIR_INADEQUATE",
    "MITO_CONCORDANT",
    "MITO_DISCORDANT",
    "MITO_MIXED",
    "MITO_UNAVAILABLE",
    "PLOIDY_UNRESOLVED",
    "OCCUPANCY_UNRESOLVED",
    "METADATA_OR_LIBRARY_EXCHANGE_UNRESOLVED",
    "CHINOBO_CONGO_RELATED_PAIR",
)
FINAL_EVENT_DISPOSITION_FIELDS = [
    "library", "event_id", "event_class", "event_confidence",
    "n_initial_cells", "n_interpreted_cells",
    "n_nuclear_proposal_direction", "n_nuclear_current_direction",
    "n_nuclear_nondirectional", "nuclear_event_pattern",
    "controlled_ambient_pattern", "refitted_ambient_pattern",
    "identity_evidence_disposition", "evidence_modifiers",
    "mechanism_disposition", "review_scope", "review_reasons",
    "production_assignment_effect", "final_schema_version",
]
FINAL_ACTIONABLE_REVIEW_FIELDS = [
    "review_level", "library", "event_id", "barcode",
    "identity_evidence_disposition", "mechanism_disposition",
    "review_scope", "review_reasons", "review_disposition",
    "review_rationale", "comparison_current_assignment",
    "nominated_proposal", "production_assignment",
    "production_assignment_effect", "final_schema_version",
]
FINAL_AMBIENT_FIELDS = [
    "ambient_frozen_current_c", "ambient_frozen_proposal_c",
    "ambient_frozen_proposal_minus_current_c",
    "ambient_frozen_current_fit", "ambient_frozen_proposal_fit",
    "ambient_frozen_fit_delta",
    "ambient_frozen_current_exact_candidate_burden",
    "ambient_frozen_proposal_exact_candidate_burden",
    "ambient_frozen_profile_id", "ambient_frozen_evaluation_status",
    "ambient_arm_a_c", "ambient_arm_b_c", "ambient_arm_c_c",
    "ambient_arm_d_c", "ambient_roster_effect_b_minus_a",
    "ambient_assignment_effect_c_minus_b",
    "ambient_replacement_effect_d_minus_c",
    "ambient_combined_augmented_c_minus_a", "ambient_production_arm",
    "ambient_production_c", "ambient_production_minus_original_c",
    "ambient_exact_donor_burden_fields", "ambient_background_shift_fields",
    "ambient_evaluation_status",
]
FINAL_CELL_FIELDS = [
    "library", "barcode", "demux_original_assignment",
    "demux_original_assignment_raw", "refined_assignment",
    "comparison_current_assignment", "nominated_proposal", "event_id",
    "preliminary_reconciliation_action", "preliminary_action_applied",
    "preliminary_reconciled_assignment", "preliminary_decision_confidence",
    "preliminary_decision_reason_codes", "interpreted_identity",
    "scientific_recommendation", "recommendation_basis", "review_required",
    "review_reasons", "review_disposition", "review_rationale",
    "review_record_scope", "review_record_target",
    "application_state", "application_reason", "production_assignment",
    "production_assignment_source", "candidate_roster_relationship",
    "proposal_kind", "proposal_components", "axis_candidate_a_assignment",
    "axis_candidate_b_assignment",
    "axis_pair_relationship_to_current_proposal", "candidate_axis_scope_status",
    "candidate_axis_exclusion_reason", "nuclear_reconciliation_status",
    "nuclear_warning_reasons", "current_donor_composition",
    "proposal_donor_composition", "current_ploidy_state",
    "proposal_ploidy_state", "ploidy_evidence_status", "occupancy_state",
    "occupancy_evidence_status", "known_line_relationship", "technical_state",
    "species_evidence_status", "mitochondrial_evidence_status",
    "mitochondrial_resolution_status", "atac_evidence_status",
] + FINAL_AMBIENT_FIELDS + [
    "uid_resolution_status", "uid_or_uid_set", "metadata_event_status",
    "metadata_amendment_proposal", "library_exchange_status",
    "downstream_release_status", "downstream_exclusion_reason",
    "event_identity_evidence_disposition", "event_review_scope",
    "review_scope", "cell_exception_reasons",
    "policy_version", "run_id", "final_schema_version",
]


def finalize_parse_args():
    p = argparse.ArgumentParser(
        description="Finalize validated reconciliation with frozen nuclear and ambient evidence."
    )
    p.add_argument("--libraries", nargs="+", required=True)
    p.add_argument("--demux-root", required=True)
    p.add_argument("--candidate-root", required=True)
    p.add_argument("--decisions-root", required=True)
    p.add_argument("--candidate-axis-root", default="")
    p.add_argument("--frozen-ambient-root", default="")
    p.add_argument("--four-arm-root", default="")
    p.add_argument("--review-input", default="")
    p.add_argument("--evidence-mode", choices=("rna", "rna-atac"), default="rna")
    p.add_argument("--run-id", default="")
    p.add_argument("--output-root", required=True)
    return p.parse_args()


def _final_library(value):
    value = clean(value)
    match = re.fullmatch(r"(?:lib)?(\d+)", value, flags=re.I)
    return f"lib{int(match.group(1))}" if match else value


def _final_bool(value):
    return clean(value).lower() in {"1", "true", "yes", "y"}


def _final_na(value):
    value = clean(value)
    return value if value and value.upper() != "NAN" else "NA"


def _final_float(value):
    value = ffloat(value)
    return value if math.isfinite(value) else math.nan


def _final_counter(values):
    counts = Counter(_final_na(value) for value in values)
    return ";".join(f"{key}:{counts[key]}" for key in sorted(counts, key=natural_key)) or "NONE"


def _final_set(values):
    return ";".join(sorted({clean(value) for value in values if clean(value)}, key=natural_key))


def _final_counter_keys(value):
    keys = set()
    for item in clean(value).split(";"):
        item = clean(item)
        if not item or item == "NONE":
            continue
        keys.add(item.rsplit(":", 1)[0] if ":" in item else item)
    return keys


def _final_nuclear_event_pattern(voting_rows, exclusion_rows, applicable):
    if not applicable:
        return "NOT_APPLICABLE", 0, 0, 0
    statuses = [
        clean(row.get("nuclear_reconciliation_status")).upper()
        for row in voting_rows
    ]
    proposal = statuses.count("NUCLEAR_SUPPORTS_PROPOSAL")
    current = statuses.count("NUCLEAR_SUPPORTS_CURRENT")
    nondirectional = len(statuses) - proposal - current
    if proposal and not current:
        return "ALL_PROPOSAL", proposal, current, nondirectional
    if current and not proposal:
        return "ALL_CURRENT", proposal, current, nondirectional
    if proposal > current:
        return (
            "PROPOSAL_PLURALITY_WITH_CURRENT_EXCEPTIONS",
            proposal, current, nondirectional,
        )
    if current > proposal:
        return (
            "CURRENT_PLURALITY_WITH_PROPOSAL_EXCEPTIONS",
            proposal, current, nondirectional,
        )
    if proposal and proposal == current:
        return (
            "DIRECTION_TIE_OR_NONDISCRIMINATING",
            proposal, current, nondirectional,
        )
    unavailable_only = {
        "NUCLEAR_UNAVAILABLE",
        "NUCLEAR_PAIR_DOES_NOT_MATCH_CURRENT_COMPARISON",
        "NUCLEAR_OUTSIDE_PAIR_CONTRACT",
    }
    nondiscriminating = {
        "NUCLEAR_NONDISCRIMINATING", "NUCLEAR_PAIR_INADEQUATE",
    }
    if any(status in nondiscriminating for status in statuses):
        return (
            "DIRECTION_TIE_OR_NONDISCRIMINATING",
            proposal, current, nondirectional,
        )
    exclusion_text = ";".join(
        clean(row.get("exclusion_reason")).upper() for row in exclusion_rows
    )
    if any(token in exclusion_text for token in (
            "NONDISCRIMINATING", "NON_DISCRIMINATING",
            "INSUFFICIENT_CANDIDATE_SEPARATION", "SAME_SNP_IDENTITY",
            "PAIR_INADEQUATE")):
        return (
            "DIRECTION_TIE_OR_NONDISCRIMINATING",
            proposal, current, nondirectional,
        )
    if statuses and all(status == "NUCLEAR_NOT_APPLICABLE" for status in statuses):
        return "NOT_APPLICABLE", proposal, current, nondirectional
    if not statuses or all(
            status in unavailable_only | {"", "NA"} for status in statuses):
        return "UNAVAILABLE", proposal, current, nondirectional
    return "UNAVAILABLE", proposal, current, nondirectional


def _final_ambient_pattern(value, applicable, prefix):
    if not applicable:
        return f"{prefix}_NOT_APPLICABLE"
    number = _final_float(value)
    if not math.isfinite(number):
        return f"{prefix}_UNAVAILABLE"
    if number < 0:
        return f"{prefix}_SUPPORTS_PROPOSAL"
    if number > 0:
        return f"{prefix}_SUPPORTS_CURRENT"
    return f"{prefix}_NEUTRAL"


def _final_identity_evidence_disposition(
        event_class, applicable, nuclear_pattern, controlled_pattern):
    if clean(event_class).upper() == "BELOW_EVENT_MASS_THRESHOLD":
        return "BELOW_EVENT_MASS_THRESHOLD"
    if not applicable:
        return "NOT_APPLICABLE"
    nuclear_proposal = nuclear_pattern in FINAL_NUCLEAR_PROPOSAL_PATTERNS
    nuclear_current = nuclear_pattern in FINAL_NUCLEAR_CURRENT_PATTERNS
    controlled_proposal = controlled_pattern == "CONTROLLED_SUPPORTS_PROPOSAL"
    controlled_current = controlled_pattern == "CONTROLLED_SUPPORTS_CURRENT"
    if nuclear_proposal and controlled_proposal:
        return "RNA_SUPPORTED_RECONCILIATION"
    if nuclear_current and controlled_current:
        return "RNA_SUPPORTS_CURRENT"
    if ((nuclear_proposal and controlled_current)
            or (nuclear_current and controlled_proposal)):
        return "RNA_DISCORDANT"
    if controlled_proposal and not (nuclear_proposal or nuclear_current):
        return "AMBIENT_SUPPORTED_NUCLEAR_INCONCLUSIVE"
    if nuclear_proposal and controlled_pattern in {
            "CONTROLLED_NEUTRAL", "CONTROLLED_UNAVAILABLE",
            "CONTROLLED_NOT_APPLICABLE"}:
        return "NUCLEAR_SUPPORTED_AMBIENT_INCONCLUSIVE"
    return "RNA_INSUFFICIENT"


def _final_status_has_unresolved(value):
    value = clean(value).upper()
    return any(token in value for token in (
        "UNRESOLVED", "AMBIG", "CONFLICT", "MISSING", "UNAVAILABLE",
        "NO_LIBRARY_METADATA_MATCH", "MULTIPLE_",
    ))


def _final_mechanism_disposition(event, linked, applicable):
    if not applicable or clean(event.get("event_class")).upper() == "BELOW_EVENT_MASS_THRESHOLD":
        return "NOT_APPLICABLE"
    event_class = clean(event.get("event_class")).upper()
    technical = any(
        clean(row.get("downstream_release_status")).upper()
        == "EXCLUDED_TECHNICAL_MULTIPLET"
        for row in linked
    )
    if technical or "TECHNICAL" in event_class or "MULTIPLET" in event_class:
        return "TECHNICAL_MIXTURE_OR_MULTIPLET"
    exchange_unresolved = any(
        any(token in clean(row.get("library_exchange_status")).upper()
            for token in ("PARTIAL_RECIPROCAL", "UNRESOLVED", "ROSTER_EQUIVALENT"))
        for row in linked
    )
    uid_unresolved = _final_status_has_unresolved(event.get("uid_resolution_status"))
    if exchange_unresolved or uid_unresolved:
        return "METADATA_OR_LIBRARY_EXCHANGE_UNRESOLVED"
    occupancy_unresolved = (
        "OCCUPANCY" in event_class
        or "CELLULAR_ORIGIN" in event_class
        or any(_final_status_has_unresolved(row.get("occupancy_evidence_status"))
               for row in linked)
    )
    if occupancy_unresolved:
        return "OCCUPANCY_OR_CELLULAR_ORIGIN_UNRESOLVED"
    if clean(event.get("event_confidence")).upper() in {"STRONG", "DECISIVE"}:
        return "MECHANISM_SUPPORTED"
    return "MECHANISM_INSUFFICIENT"


def _final_mito_modifier(linked):
    statuses = {
        clean(row.get("mitochondrial_evidence_status")).upper()
        for row in linked
    }
    proposal = bool(statuses & {
        "SUPPORTS_ALTERNATIVE", "SUPPORTS_PROPOSAL", "CONCORDANT",
    })
    current = bool(statuses & {"SUPPORTS_CURRENT", "CONTRADICTS"})
    if proposal and current:
        return "MITO_MIXED"
    if proposal:
        return "MITO_CONCORDANT"
    if current:
        return "MITO_DISCORDANT"
    return "MITO_UNAVAILABLE"


def _final_evidence_modifiers(
        event, linked, voting_rows, nuclear_pattern, controlled_pattern,
        refitted_pattern, mechanism_disposition):
    modifiers = set()
    if controlled_pattern == "CONTROLLED_SUPPORTS_PROPOSAL":
        if refitted_pattern in {"REFIT_NEUTRAL", "REFIT_SUPPORTS_CURRENT"}:
            modifiers.add("AMBIENT_REFIT_SENSITIVE")
        elif refitted_pattern in {
                "REFIT_UNAVAILABLE", "REFIT_NOT_APPLICABLE"}:
            modifiers.add("AMBIENT_REFIT_UNAVAILABLE")
    if nuclear_pattern in {
            "PROPOSAL_PLURALITY_WITH_CURRENT_EXCEPTIONS",
            "CURRENT_PLURALITY_WITH_PROPOSAL_EXCEPTIONS"}:
        modifiers.add("NUCLEAR_CELL_HETEROGENEITY")
    candidate_a_excess = [
        _final_float(row.get("candidate_a_excess_brier_mean"))
        for row in voting_rows
    ]
    candidate_b_excess = [
        _final_float(row.get("candidate_b_excess_brier_mean"))
        for row in voting_rows
    ]
    candidate_a_excess = [
        value for value in candidate_a_excess if math.isfinite(value)
    ]
    candidate_b_excess = [
        value for value in candidate_b_excess if math.isfinite(value)
    ]
    sampling_direction = "TIE_OR_UNAVAILABLE"
    if candidate_a_excess and candidate_b_excess:
        median_a = median(candidate_a_excess)
        median_b = median(candidate_b_excess)
        sampling_direction = (
            "CURRENT" if median_a < median_b else
            "PROPOSAL" if median_b < median_a else
            "TIE_OR_UNAVAILABLE"
        )
    if ((nuclear_pattern in FINAL_NUCLEAR_PROPOSAL_PATTERNS
         and sampling_direction == "CURRENT")
            or (nuclear_pattern in FINAL_NUCLEAR_CURRENT_PATTERNS
                and sampling_direction == "PROPOSAL")):
        modifiers.add("SAMPLING_ADJUSTED_FIT_DISAGREEMENT")
    nuclear_statuses = {
        clean(row.get("nuclear_reconciliation_status")).upper()
        for row in voting_rows
    }
    if (nuclear_pattern == "DIRECTION_TIE_OR_NONDISCRIMINATING"
            or "NUCLEAR_NONDISCRIMINATING" in nuclear_statuses):
        modifiers.add("NUCLEAR_NONDISCRIMINATING")
    if "NUCLEAR_UNSTABLE" in nuclear_statuses:
        modifiers.add("NUCLEAR_UNSTABLE")
    if nuclear_statuses & {
            "NUCLEAR_PAIR_INADEQUATE",
            "NUCLEAR_PAIR_DOES_NOT_MATCH_CURRENT_COMPARISON",
            "NUCLEAR_OUTSIDE_PAIR_CONTRACT"}:
        modifiers.add("PAIR_INADEQUATE")
    exclusion_text = ";".join(
        clean(row.get("candidate_axis_exclusion_reason")).upper()
        for row in linked
    )
    if any(token in exclusion_text for token in (
            "INADEQUATE", "OUTSIDE_PAIR_CONTRACT", "PAIR_DOES_NOT_MATCH",
            "NO_COMMON_EVIDENCE")):
        modifiers.add("PAIR_INADEQUATE")
    modifiers.add(_final_mito_modifier(linked))
    ploidy_statuses = {
        clean(row.get("ploidy_evidence_status")).upper() for row in linked
    }
    if ploidy_statuses & {"UNAVAILABLE", "CONFLICTED", "UNRESOLVED"}:
        modifiers.add("PLOIDY_UNRESOLVED")
    if any(_final_status_has_unresolved(row.get("occupancy_evidence_status"))
           for row in linked):
        modifiers.add("OCCUPANCY_UNRESOLVED")
    if mechanism_disposition == "METADATA_OR_LIBRARY_EXCHANGE_UNRESOLVED":
        modifiers.add("METADATA_OR_LIBRARY_EXCHANGE_UNRESOLVED")
    relationship_pairs = [
        (event.get("primary_source_identity"), event.get("unexpected_component")),
        *(
            (row.get("comparison_current_assignment"),
             row.get("nominated_proposal"))
            for row in linked
        ),
    ]
    if any(
            {"CHINOBO-MCHERRY", "CONGOA4B"} <= {
                component.upper()
                for value in pair
                for component in donor_components(
                    canonical_genotype(value or ""))
            }
            for pair in relationship_pairs):
        modifiers.add("CHINOBO_CONGO_RELATED_PAIR")
    return ";".join(
        modifier for modifier in FINAL_PHASE3_MODIFIER_ORDER
        if modifier in modifiers
    ) or "NONE"


def _final_event_review_scope(
        event, identity_disposition, mechanism_disposition, modifiers,
        explicit_event_review):
    if identity_disposition == "BELOW_EVENT_MASS_THRESHOLD":
        if explicit_event_review:
            return "EVENT_REVIEW", "EXPLICIT_EVENT_REVIEW_RECORD"
        return "NO_IMMEDIATE_REVIEW", "NONE"
    reasons = set()
    if identity_disposition in {
            "RNA_SUPPORTS_CURRENT", "RNA_DISCORDANT", "RNA_INSUFFICIENT",
            "AMBIENT_SUPPORTED_NUCLEAR_INCONCLUSIVE",
            "NUCLEAR_SUPPORTED_AMBIENT_INCONCLUSIVE"}:
        reasons.add("IDENTITY_EVIDENCE=" + identity_disposition)
    if mechanism_disposition not in {"MECHANISM_SUPPORTED", "NOT_APPLICABLE"}:
        reasons.add("MECHANISM=" + mechanism_disposition)
    modifier_set = set(clean(modifiers).split(";"))
    for modifier in (
            "AMBIENT_REFIT_SENSITIVE", "NUCLEAR_CELL_HETEROGENEITY",
            "SAMPLING_ADJUSTED_FIT_DISAGREEMENT", "NUCLEAR_UNSTABLE",
            "PAIR_INADEQUATE", "MITO_DISCORDANT", "MITO_MIXED"):
        if modifier in modifier_set:
            reasons.add(modifier)
    raw_reasons = _final_counter_keys(event.get("review_reason_counts"))
    for reason in raw_reasons:
        upper = reason.upper()
        if (upper in {"EVENT_MECHANISM_REQUIRES_REVIEW", "LIBRARY_EXCHANGE_AMBIGUITY"}
                or "METADATA_UID_" in upper
                or upper.startswith("OCCUPANCY_")):
            reasons.add(reason)
    if explicit_event_review:
        reasons.add("EXPLICIT_EVENT_REVIEW_RECORD")
    if reasons:
        return "EVENT_REVIEW", _final_set(reasons)
    return "NO_IMMEDIATE_REVIEW", "NONE"


def _final_cell_exception_reasons(row):
    reasons = set()
    production = canonical_genotype(row.get("production_assignment", ""))
    current = canonical_genotype(row.get("comparison_current_assignment", ""))
    proposal = canonical_genotype(row.get("nominated_proposal", ""))
    nuclear = clean(row.get("nuclear_reconciliation_status")).upper()
    if proposal and proposal != current:
        if nuclear == "NUCLEAR_SUPPORTS_CURRENT" and production == proposal:
            reasons.add("STABLE_CURRENT_SIDE_NUCLEAR_EVIDENCE_OPPOSES_PRODUCTION")
        elif nuclear == "NUCLEAR_SUPPORTS_PROPOSAL" and production == current:
            reasons.add("CELL_LOCAL_NUCLEAR_ASSIGNMENT_CONFLICT")
    ploidy = clean(row.get("ploidy_evidence_status")).upper()
    if ploidy in {"CONFLICTED", "SUPPORTS_CURRENT"} and production == proposal:
        reasons.add("CELL_LOCAL_PLOIDY_CONFLICT")
    occupancy = clean(row.get("occupancy_evidence_status")).upper()
    if "CONFLICT" in occupancy and "UNRESOLVED" not in occupancy:
        reasons.add("CELL_LOCAL_OCCUPANCY_CONFLICT")
    if (clean(row.get("downstream_release_status")).upper()
            == "EXCLUDED_TECHNICAL_MULTIPLET"):
        reasons.add("TECHNICAL_MIXTURE_OR_MULTIPLET")
    warnings = clean(row.get("nuclear_warning_reasons")).upper()
    if "UPSTREAM_" in warnings:
        reasons.add("CELL_LOCAL_UPSTREAM_ASSIGNMENT_WARNING")
    if clean(row.get("review_record_scope")).upper() == "CELL":
        reasons.add("EXPLICIT_CELL_REVIEW_RECORD")
    return _final_set(reasons) or "NONE"


def _final_markdown_counts(values):
    counts = Counter(clean(value) or "NA" for value in values)
    return "\n".join(
        f"- `{key}`: {counts[key]}"
        for key in sorted(counts, key=natural_key)
    ) or "- `NONE`: 0"


def _final_phase3_summary_text(
        libraries, final_cells, final_events, actionable_review, review_input,
        used_review_inputs, evidence_mode, nonfinalized_worsened):
    library_numbers = [int(re.sub(r"\D", "", library)) for library in libraries]
    deferred = [number for number in (38, 40) if number not in library_numbers]
    assignment_changes = sum(
        row["production_assignment"] != row["preliminary_reconciled_assignment"]
        for row in final_cells
    )
    modifier_values = [
        modifier
        for row in final_events
        for modifier in clean(row.get("evidence_modifiers")).split(";")
        if modifier and modifier != "NONE"
    ]

    def event_text(row):
        return (
            f"{_final_library(row.get('library'))}/{clean(row.get('event_id'))}: "
            f"{clean(row.get('identity_evidence_disposition'))}; "
            f"mechanism={clean(row.get('mechanism_disposition'))}; "
            f"modifiers={clean(row.get('evidence_modifiers'))}"
        )

    def exact_identity_event(library, candidate):
        event_id = f"{library}:identity:{candidate}"
        return [
            row for row in final_events
            if _final_library(row.get("library")) == library
            and clean(row.get("event_id")) == event_id
            and canonical_genotype(row.get("unexpected_component", ""))
            == candidate
        ]

    def has_directional_candidate_axis(row):
        return (
            fint(row.get("n_nuclear_proposal_direction"), 0)
            + fint(row.get("n_nuclear_current_direction"), 0)
        ) > 0

    named_specs = [
        ("lib12 C40210+CongoA4B", "lib12", canonical_genotype("C40210+CongoA4B")),
        ("lib4 C8861+H20961", "lib4", canonical_genotype("C8861+H20961")),
        ("lib36 C3651+H20961", "lib36", canonical_genotype("C3651+H20961")),
        ("lib36 Chinobo-mCherry+H20961", "lib36", canonical_genotype("Chinobo-mCherry+H20961")),
        ("lib11 Chinobo-mCherry+CongoA4B", "lib11", canonical_genotype("Chinobo-mCherry+CongoA4B")),
    ]
    special_lines = []
    for label, library, candidate in named_specs:
        matches = exact_identity_event(library, candidate)
        special_lines.append(
            f"- {label}: "
            + (" | ".join(event_text(row) for row in matches)
               if matches else "NOT_FOUND_IN_SELECTED_EVENT_LEDGER")
        )
    for proposal in ("H20157", "H29089", "H21194", "H27322"):
        matches = exact_identity_event("lib7", canonical_genotype(proposal))
        special_lines.append(
            f"- lib7 proposal group {proposal}: "
            + (" | ".join(event_text(row) for row in matches)
               if matches else "NOT_FOUND_IN_SELECTED_EVENT_LEDGER")
        )

    sampling_proposal_reversal_rows = sorted((
        row for row in final_events
        if "SAMPLING_ADJUSTED_FIT_DISAGREEMENT" in
        clean(row.get("evidence_modifiers")).split(";")
        and row.get("nuclear_event_pattern") in
        FINAL_NUCLEAR_PROPOSAL_PATTERNS
    ), key=lambda row: (
        natural_key(_final_library(row.get("library"))),
        natural_key(clean(row.get("event_id"))),
    ))
    sampling_current_reversal_rows = sorted((
        row for row in final_events
        if "SAMPLING_ADJUSTED_FIT_DISAGREEMENT" in
        clean(row.get("evidence_modifiers")).split(";")
        and row.get("nuclear_event_pattern") in
        FINAL_NUCLEAR_CURRENT_PATTERNS
    ), key=lambda row: (
        natural_key(_final_library(row.get("library"))),
        natural_key(clean(row.get("event_id"))),
    ))
    refit_reversal_rows = [
        row for row in final_events
        if has_directional_candidate_axis(row)
        and row.get("controlled_ambient_pattern") == "CONTROLLED_SUPPORTS_PROPOSAL"
        and row.get("refitted_ambient_pattern") == "REFIT_SUPPORTS_CURRENT"
    ]
    supported_rows = sorted(
        (row for row in final_events
         if row.get("identity_evidence_disposition") ==
         "RNA_SUPPORTED_RECONCILIATION"
         and row.get("mechanism_disposition") == "MECHANISM_SUPPORTED"
         and row.get("review_scope") == "NO_IMMEDIATE_REVIEW"
         and set(clean(row.get("evidence_modifiers")).split(";"))
         <= {"NONE", "MITO_CONCORDANT"}),
        key=lambda row: (
            natural_key(_final_library(row.get("library"))),
            natural_key(clean(row.get("event_id"))),
        ),
    )
    event_libraries = {
        _final_library(row.get("library")) for row in final_events
    }
    zero_event_libraries = [
        library for library in libraries if library not in event_libraries
    ]
    not_applicable_rows = sorted(
        (row for row in final_events
         if row.get("identity_evidence_disposition") == "NOT_APPLICABLE"),
        key=lambda row: (
            natural_key(_final_library(row.get("library"))),
            natural_key(clean(row.get("event_id"))),
        ),
    )
    zero_or_not_applicable_example = (
        zero_event_libraries[0] if zero_event_libraries else
        event_text(not_applicable_rows[0]) if not_applicable_rows else
        "NONE_AVAILABLE"
    )
    special_lines.extend([
        "- Proposal-side median sampling-adjusted fit reversals: "
        + (" | ".join(
            event_text(row) for row in sampling_proposal_reversal_rows)
           if sampling_proposal_reversal_rows
           else "NONE_IN_SELECTED_EVENT_LEDGER"),
        "- Separate current-side median sampling-adjusted disagreement "
        "(excluded from the proposal-side reversal set): "
        + (" | ".join(
            event_text(row) for row in sampling_current_reversal_rows)
           if sampling_current_reversal_rows
           else "NONE_IN_SELECTED_EVENT_LEDGER"),
        "- Finalized refitted-ambient reversals: "
        + (" | ".join(event_text(row) for row in refit_reversal_rows)
           if refit_reversal_rows else "NONE_IN_SELECTED_EVENT_LEDGER"),
        "- Non-finalized controlled-ambient worsened candidate identities: "
        + ("; ".join(nonfinalized_worsened)
           if nonfinalized_worsened else "NONE_IN_AVAILABLE_CONTROLLED_EVIDENCE"),
        "- Clean strongly supported example: "
        + (event_text(supported_rows[0]) if supported_rows else "NONE_AVAILABLE"),
        "- Zero-event/not-applicable example: "
        + zero_or_not_applicable_example,
    ])

    return (
        "# Identity Reconciliation Phase 3 Summary\n\n"
        f"Selected libraries: {' '.join(str(number) for number in library_numbers)}  \n"
        f"Deferred libraries: {' '.join(str(number) for number in deferred) or 'NONE'}  \n"
        f"Evidence mode: {evidence_mode.upper()}  \n"
        f"ATAC: {'ATAC_NOT_REQUESTED' if evidence_mode == 'rna' else 'RNA_ATAC'}  \n"
        f"Review records present: {len(review_input)}  \n"
        f"Review records applied: {len(used_review_inputs)}  \n"
        f"Assignments changed relative to accepted Phase 2 production assignments: "
        f"{assignment_changes}\n\n"
        "## Accounting\n\n"
        f"- Canonical ledger barcodes: {len(final_cells)}\n"
        f"- Event/population rows: {len(final_events)}\n"
        f"- Actionable event-review rows: "
        f"{sum(row['review_level'] == 'EVENT_REVIEW' for row in actionable_review)}\n"
        f"- Actionable cell-exception rows: "
        f"{sum(row['review_level'] == 'CELL_EXCEPTION_REVIEW' for row in actionable_review)}\n\n"
        "## Identity evidence dispositions\n\n"
        + _final_markdown_counts(
            row.get("identity_evidence_disposition") for row in final_events
        ) + "\n\n"
        "## Mechanism dispositions\n\n"
        + _final_markdown_counts(
            row.get("mechanism_disposition") for row in final_events
        ) + "\n\n"
        "## Evidence modifiers\n\n"
        + _final_markdown_counts(modifier_values) + "\n\n"
        "## Review scope counts\n\n"
        "Event/population scope:\n\n"
        + _final_markdown_counts(row.get("review_scope") for row in final_events)
        + "\n\nCell scope:\n\n"
        + _final_markdown_counts(row.get("review_scope") for row in final_cells)
        + "\n\nActionable rows:\n\n"
        + _final_markdown_counts(
            row.get("review_level") for row in actionable_review
        ) + "\n\n"
        "## Fixed special-resolution set\n\n"
        + "\n".join(special_lines) + "\n\n"
        "The preliminary all-40 `reports/summary.md` is not authoritative for "
        "this selected-scope Phase 3 result.\n"
    )


def _final_nonfinalized_controlled_worsened(
        controlled, controlled_planned, selected_libraries,
        finalized_candidate_axis_event_keys):
    grouped = defaultdict(list)
    for (library, barcode), rows in controlled.items():
        if library not in selected_libraries:
            continue
        for row in rows:
            event_id = clean(row.get("event_id"))
            candidate = canonical_genotype(row.get("candidate_identity", ""))
            plan_key = (library, barcode, event_id, candidate)
            if (not event_id or event_id.upper() == "NA"
                    or not candidate or candidate.upper() == "NA"
                    or plan_key not in controlled_planned):
                continue
            event_key = (library, event_id, candidate)
            value = _final_float(row.get("fixed_delta_c"))
            if (event_key not in finalized_candidate_axis_event_keys
                    and math.isfinite(value)):
                grouped[event_key].append(value)
    return [
        f"{key[0]}/{key[1]}/{key[2]}"
        for key, values in sorted(
            grouped.items(),
            key=lambda item: tuple(natural_key(value) for value in item[0]),
        )
        if values and median(values) > 0
    ]


def _final_read(path):
    if not path:
        return []
    path = str(path)
    return read_tsv(path) if os.path.isfile(path) else []


def _final_merge_candidate_evidence(candidate, evidence, source):
    for field, value in evidence.items():
        if field not in candidate or not clean(candidate.get(field)):
            candidate[field] = value
            continue
        if clean(candidate.get(field)) == clean(value):
            continue
        candidate[f"{source}_{field}"] = value


def _final_axis_status(axis, exclusion, current, proposal, applicable):
    warnings = set()
    if not applicable:
        return "NUCLEAR_NOT_APPLICABLE", warnings
    if exclusion:
        return "NUCLEAR_OUTSIDE_PAIR_CONTRACT", warnings
    if not axis:
        return "NUCLEAR_UNAVAILABLE", warnings
    status = clean(axis.get("candidate_axis_status")).upper()
    comparison = clean(axis.get("comparison_status_legacy")).upper()
    structural_nonseparation = (
        status == "INSUFFICIENT_CANDIDATE_SEPARATION"
        or comparison == "PANEL_NONDISCRIMINATING")
    if ((status not in {"AVAILABLE", "INSUFFICIENT_CANDIDATE_SEPARATION"}
         and not structural_nonseparation)
            or comparison in {"NO_COMMON_EVIDENCE", ""}):
        return "NUCLEAR_UNAVAILABLE", warnings
    candidate_a = canonical_genotype(axis.get("candidate_a", ""))
    candidate_b = canonical_genotype(axis.get("candidate_b", ""))
    if candidate_a != current or candidate_b != proposal:
        return "NUCLEAR_PAIR_DOES_NOT_MATCH_CURRENT_COMPARISON", warnings
    if structural_nonseparation:
        return "NUCLEAR_NONDISCRIMINATING", warnings
    direction = clean(axis.get("candidate_axis_direction")).upper()
    if comparison == "LOW_EVIDENCE":
        return "NUCLEAR_PAIR_INADEQUATE", warnings
    if comparison == "NO_CANDIDATE_FITS":
        warnings.add("LEGACY_RAW_RESIDUAL_NO_CANDIDATE_FITS")
    raw_flag = clean(axis.get("raw_residual_threshold_flag")).upper()
    if raw_flag not in {"", "NA", "UNAVAILABLE", "NEITHER_ABOVE_LEGACY_THRESHOLD"}:
        warnings.add("RAW_RESIDUAL_THRESHOLD_WARNING=" + raw_flag)
    unstable = (
        clean(axis.get("candidate_axis_fold_direction_stability_status")).upper()
        in {"DIRECTION_CHANGED_OR_TIED", "FOLD_RECONSTRUCTION_MISMATCH"}
        or clean(axis.get("candidate_axis_direction_preserved_without_top_primary_unit")).upper()
        == "FALSE"
        or clean(axis.get("candidate_axis_direction_preserved_without_top_five_primary_units")).upper()
        == "FALSE"
    )
    if unstable:
        return "NUCLEAR_UNSTABLE", warnings
    a_excess = _final_float(axis.get("candidate_a_excess_brier_mean"))
    b_excess = _final_float(axis.get("candidate_b_excess_brier_mean"))
    if math.isfinite(a_excess) and math.isfinite(b_excess) and a_excess != b_excess:
        adjusted = "ORIGINAL_SIDE" if a_excess < b_excess else "PROPOSAL_SIDE"
        if adjusted != direction:
            warnings.add("SAMPLING_ADJUSTED_BRIER_DIRECTION_DISAGREEMENT")
    if direction == "PROPOSAL_SIDE":
        return "NUCLEAR_SUPPORTS_PROPOSAL", warnings
    if direction == "ORIGINAL_SIDE":
        return "NUCLEAR_SUPPORTS_CURRENT", warnings
    if direction == "TIE":
        return "NUCLEAR_NONDISCRIMINATING", warnings
    return "NUCLEAR_UNAVAILABLE", warnings


def _final_inferred_ploidy(genotype):
    n_components = len(donor_components(genotype))
    if n_components == 1:
        return "DIPLOID"
    if n_components == 2:
        return "TETRAPLOID"
    if n_components >= 3:
        return "UNRESOLVED_MULTIPLET"
    return "UNAVAILABLE"


def _final_ploidy_status(preliminary, current, proposal):
    if not proposal or proposal == current:
        return "NOT_APPLICABLE"
    current_state = (
        clean(preliminary.get("current_ploidy"))
        or _final_inferred_ploidy(current))
    proposal_state = (
        clean(preliminary.get("proposed_biological_ploidy"))
        or _final_inferred_ploidy(proposal))
    if current_state.upper() == proposal_state.upper():
        return "NOT_APPLICABLE"
    if not _final_bool(preliminary.get("nn_qc_pass")):
        return "UNAVAILABLE"
    call = clean(preliminary.get("nn_ploidy_call")).upper()
    if not call:
        return "UNAVAILABLE"
    if proposal_state.upper().startswith(call[:3]) or call.startswith(proposal_state.upper()[:3]):
        return "SUPPORTS_PROPOSAL"
    if current_state.upper().startswith(call[:3]) or call.startswith(current_state.upper()[:3]):
        return "SUPPORTS_CURRENT"
    return "CONFLICTED"


def _final_review_input(path):
    if not path:
        return {}
    if not os.path.isfile(path):
        raise FileNotFoundError(f"review input does not exist: {path}")
    rows = read_tsv(path)
    out = {}
    for row in rows:
        library = _final_library(row.get("library"))
        target = clean(row.get("barcode_or_event_id"))
        disposition = clean(row.get("review_disposition")).upper()
        if not library or not target or disposition not in FINAL_REVIEW_DISPOSITIONS:
            raise ValueError("review input has an invalid library, target, or disposition")
        key = (library, target)
        if key in out:
            raise ValueError(f"duplicate review input target: {library}/{target}")
        out[key] = {
            "disposition": disposition,
            "rationale": clean(row.get("rationale")),
        }
    return out


def _final_demux_warning_rows(path):
    warnings = {}
    if not path or not os.path.isfile(path):
        return warnings, "UNAVAILABLE"
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = set(reader.fieldnames or ())
        if "maximin_score" not in fields or not ({"barcode", "cell"} & fields):
            return warnings, "UNAVAILABLE"
        for row in reader:
            barcode = clean(row.get("barcode") or row.get("cell"))
            score = _final_float(row.get("maximin_score"))
            if barcode and math.isfinite(score) and score <= 0:
                winner = canonical_genotype(
                    row.get("assignment", "")
                    or row.get("maximin_winner", ""))
                warnings[barcode] = (
                    "UPSTREAM_B041_NONPOSITIVE_MAXIMIN"
                    if "B041" in donor_components(winner) else
                    "UPSTREAM_ACTIVE_NONPOSITIVE_MAXIMIN"
                )
    return warnings, "AVAILABLE"


def _final_controlled_ambient(root):
    root = Path(root) if root else None
    rows = _final_read(root / "data" / "ambient_swap_cell_contrasts.tsv.gz") if root else []
    planned_rows = _final_read(root / "data" / "ambient_swap_candidate_cells.tsv") if root else []
    applicability = _final_read(root / "data" / "ambient_swap_library_applicability.tsv") if root else []
    profile_by_library = {
        _final_library(row.get("library")): _final_na(row.get("plan_fingerprint"))
        for row in applicability
    }
    indexed = defaultdict(list)
    for row in rows:
        indexed[(_final_library(row.get("library")), clean(row.get("barcode")))].append(row)
    planned = {
        (
            _final_library(row.get("library")), clean(row.get("barcode")),
            clean(row.get("event_id")),
            canonical_genotype(row.get("candidate_identity", "")),
        )
        for row in planned_rows
        if clean(row.get("barcode")) and clean(row.get("event_id"))
        and canonical_genotype(row.get("candidate_identity", ""))
    }
    return indexed, profile_by_library, planned


def _final_four_arm_ambient(root):
    root = Path(root) if root else None
    data = root / "data" if root else None
    contrasts = _final_read(data / "reconciliation_planned_contrast_cells.tsv") if data else []
    summaries = _final_read(data / "reconciliation_planned_contrasts.tsv") if data else []
    burdens = _final_read(data / "reconciliation_exact_donor_burden.tsv") if data else []
    pivot = defaultdict(dict)
    strata = {}
    for row in contrasts:
        key = (_final_library(row.get("library")), clean(row.get("barcode")))
        name = clean(row.get("contrast"))
        pivot[key][name] = row
        strata[key] = clean(row.get("stratum")) or "full_library"
    summary_index = defaultdict(list)
    for row in summaries:
        summary_index[_final_library(row.get("library"))].append(row)
    burden_index = defaultdict(list)
    for row in burdens:
        burden_index[_final_library(row.get("library"))].append(row)
    return pivot, strata, summary_index, burden_index


def _final_four_arm_fields(library, barcode, event_bearing, pivot, strata,
                           summaries, burdens):
    key = (library, barcode)
    comparisons = pivot.get(key, {})
    roster = comparisons.get("roster_effect_demux", {})
    assignment = comparisons.get("assignment_effect_augmented", {})
    replacement = comparisons.get("replacement_sensitivity", {})
    combined = comparisons.get("combined_production_change", {})
    arm_a = roster.get("left_rate") or combined.get("left_rate")
    arm_b = roster.get("right_rate") or assignment.get("left_rate")
    arm_c = assignment.get("right_rate") or replacement.get("left_rate") or combined.get("right_rate")
    arm_d = replacement.get("right_rate")
    production_arm = "C" if event_bearing else "A"
    production_c = arm_c if production_arm == "C" else arm_a
    population = strata.get(key, "full_library")
    matching_burdens = [
        row for row in burdens.get(library, [])
        if clean(row.get("arm")) == production_arm
        and clean(row.get("population")) in {population, "full_library"}
    ]
    preferred_population = population if any(
        clean(row.get("population")) == population for row in matching_burdens
    ) else "full_library"
    burden_text = ";".join(
        f"{clean(row.get('source_label'))}={_final_na(row.get('mean_exact_contam_burden'))}"
        for row in sorted(matching_burdens, key=lambda item: natural_key(clean(item.get("source_label"))))
        if clean(row.get("population")) == preferred_population
    ) or "NA"
    background_text = ";".join(
        f"{clean(row.get('contrast'))}:n={_final_na(row.get('n_common'))},mean={_final_na(row.get('mean_delta'))},median={_final_na(row.get('median_delta'))}"
        for row in sorted(summaries.get(library, []), key=lambda item: natural_key(clean(item.get("contrast"))))
        if clean(row.get("population")) == "background"
    ) or "NA"
    required = {"roster_effect_demux", "assignment_effect_augmented", "combined_production_change"}
    if not event_bearing:
        status = "NOT_APPLICABLE_ZERO_EVENT"
    elif not comparisons:
        status = "PLANNED_CELL_MISSING_FROM_ALL_RECONCILIATION_ARMS"
    elif not required <= set(comparisons):
        status = "PARTIALLY_PAIRED_RECONCILIATION_ARMS"
    elif replacement:
        status = "PAIRED_A_B_C_D"
    else:
        status = "PAIRED_A_B_C_D_NOT_APPLICABLE"
    return {
        "ambient_arm_a_c": _final_na(arm_a),
        "ambient_arm_b_c": _final_na(arm_b),
        "ambient_arm_c_c": _final_na(arm_c),
        "ambient_arm_d_c": _final_na(arm_d),
        "ambient_roster_effect_b_minus_a": _final_na(roster.get("right_minus_left")),
        "ambient_assignment_effect_c_minus_b": _final_na(assignment.get("right_minus_left")),
        "ambient_replacement_effect_d_minus_c": _final_na(replacement.get("right_minus_left")),
        "ambient_combined_augmented_c_minus_a": _final_na(combined.get("right_minus_left")),
        "ambient_production_arm": production_arm,
        "ambient_production_c": _final_na(production_c),
        "ambient_production_minus_original_c": _final_na(
            combined.get("right_minus_left") if production_arm == "C" else 0),
        "ambient_exact_donor_burden_fields": burden_text,
        "ambient_background_shift_fields": background_text,
        "ambient_evaluation_status": status,
    }


def finalize_main():
    args = finalize_parse_args()
    library_numbers = parse_library_spec(args.libraries)
    libraries = [f"lib{number}" for number in library_numbers]
    output_root = Path(args.output_root)
    output_root.mkdir(parents=True, exist_ok=True)

    decisions_root = Path(args.decisions_root)
    candidate_root = Path(args.candidate_root)
    event_rows = _final_read(decisions_root / "all_libraries.identity_events.tsv")
    exchange_rows = _final_read(
        decisions_root / "all_libraries.library_exchange_events.tsv")
    amendment_rows = _final_read(
        decisions_root / "all_libraries.metadata_amendments_proposed.tsv")
    selected_libraries = set(libraries)
    event_rows = [
        row for row in event_rows
        if _final_library(row.get("library")) in selected_libraries
    ]
    final_event_keys = {
        (
            _final_library(row.get("library")),
            clean(row.get("event_id")),
            canonical_genotype(row.get("unexpected_component", "")),
        )
        for row in event_rows
        if clean(row.get("event_id"))
        and canonical_genotype(row.get("unexpected_component", ""))
    }
    event_bearing = {
        _final_library(row.get("library")) for row in event_rows
        if _final_library(row.get("library")) in selected_libraries
    }

    exchange_by_library = defaultdict(list)
    for row in exchange_rows:
        for own, other in (("library_a", "library_b"), ("library_b", "library_a")):
            library = _final_library(row.get(own))
            if library:
                exchange_by_library[library].append(
                    f"{_final_library(row.get(other))}:"
                    f"{_final_na(row.get('exchange_interpretation'))}:"
                    f"{_final_na(row.get('reciprocal_exchange_status'))}"
                )
    amendment_by_key = defaultdict(list)
    for row in amendment_rows:
        key = (
            _final_library(row.get("library")),
            canonical_genotype(row.get("proposed_donor_genotype", "")),
        )
        amendment_by_key[key].append(
            f"{_final_na(row.get('proposed_metadata_event'))}:"
            f"{_final_na(row.get('uid_resolution_status'))}"
        )

    axis_root = Path(args.candidate_axis_root) if args.candidate_axis_root else None
    if axis_root and (axis_root / "aggregate").is_dir():
        axis_root = axis_root / "aggregate"
    axis_rows = _final_read(axis_root / "candidate_axis_cell_scores.tsv.gz") if axis_root else []
    exclusion_rows = _final_read(axis_root / "candidate_axis_pair_exclusions.tsv.gz") if axis_root else []
    axis_rows = [
        row for row in axis_rows
        if (
            _final_library(row.get("library")),
            clean(row.get("selected_supported_event_id")),
            canonical_genotype(
                row.get("selected_supported_event_proposal", "")),
        ) in final_event_keys
    ]
    exclusion_rows = [
        row for row in exclusion_rows
        if (
            _final_library(row.get("library")),
            clean(row.get("selected_supported_event_id")),
            canonical_genotype(
                row.get("selected_supported_event_proposal", "")),
        ) in final_event_keys
    ]
    axis_fields = []
    for row in axis_rows:
        for field in row:
            if field not in axis_fields:
                axis_fields.append(field)
    axis_by_cell = defaultdict(list)
    for row in axis_rows:
        axis_by_cell[(_final_library(row.get("library")), clean(row.get("barcode")))].append(row)
    exclusion_by_cell = defaultdict(list)
    for row in exclusion_rows:
        exclusion_by_cell[(_final_library(row.get("library")), clean(row.get("barcode")))].append(row)

    controlled, frozen_profile_by_library, controlled_planned = (
        _final_controlled_ambient(args.frozen_ambient_root))
    four_pivot, four_strata, four_summaries, four_burdens = (
        _final_four_arm_ambient(args.four_arm_root))
    review_input = _final_review_input(args.review_input)

    preliminary_by_cell = {}
    demux_by_library = {}
    demux_warnings = {}
    demux_warning_provenance = {}
    all_candidate_audit = []
    candidate_audit_fields = []
    for number, library in zip(library_numbers, libraries):
        prefix = demux_prefix(args.demux_root, number)
        demux = read_assignments(prefix + ".assignments")
        demux_by_library[library] = demux
        diagnostic_path = prefix + ".diagnostics.gz"
        (demux_warnings[library], demux_warning_provenance[library]) = (
            _final_demux_warning_rows(diagnostic_path))
        decision_path = decisions_root / f"{library}.reconciled_cells.tsv.gz"
        rows = _final_read(decision_path)
        indexed = {}
        for row in rows:
            barcode = clean(row.get("barcode"))
            if not barcode or barcode in indexed:
                raise ValueError(f"{library} preliminary decision table has a blank or duplicate barcode")
            indexed[barcode] = row
            preliminary_by_cell[(library, barcode)] = row
        if set(indexed) != set(demux):
            missing = sorted(set(demux) - set(indexed), key=natural_key)
            extra = sorted(set(indexed) - set(demux), key=natural_key)
            raise ValueError(
                f"{library} preliminary/demux barcode accounting mismatch; "
                f"missing={missing[:10]} extra={extra[:10]}"
            )
        candidate_rows = _final_read(candidate_root / f"{library}.identity_candidate_audit.tsv.gz")
        identity_root = decisions_root.parent
        evidence_sources = [
            ("nuclear", identity_root / "nuclear" /
             f"{library}.identity_hypothesis_scores.tsv.gz"),
            ("mt", identity_root / "mt" /
             f"{library}.mt_identity_scores.tsv.gz"),
        ]
        if args.evidence_mode == "rna-atac":
            evidence_sources.append((
                "atac", identity_root / "atac" /
                f"{library}.atac_identity_scores.tsv.gz"))
        evidence_by_source = {}
        for source, path in evidence_sources:
            indexed_evidence = {}
            for evidence in _final_read(path):
                evidence_key = (
                    clean(evidence.get("barcode")),
                    clean(evidence.get("hypothesis_id")),
                )
                if not all(evidence_key):
                    continue
                if evidence_key in indexed_evidence:
                    raise ValueError(
                        f"{library} {source} candidate evidence has a duplicate "
                        f"barcode/hypothesis key: {evidence_key[0]}/{evidence_key[1]}")
                indexed_evidence[evidence_key] = evidence
            evidence_by_source[source] = indexed_evidence
        candidate_folds = fold_map(str(
            identity_root / "nuclear" /
            f"{library}.identity_site_fold_scores.tsv.gz"))
        for row in candidate_rows:
            evidence_key = (
                clean(row.get("barcode")), clean(row.get("hypothesis_id")))
            if all(evidence_key):
                for source, indexed_evidence in evidence_by_source.items():
                    evidence = indexed_evidence.get(evidence_key)
                    if evidence:
                        _final_merge_candidate_evidence(row, evidence, source)
                fold = candidate_folds.get(evidence_key)
                if fold:
                    _final_merge_candidate_evidence(row, fold, "nuclear_fold")
            for field in row:
                if field not in candidate_audit_fields:
                    candidate_audit_fields.append(field)
        all_candidate_audit.extend(candidate_rows)

    selected_axis = {}
    selected_exclusion = {}
    used_review_inputs = set()
    final_cells = []
    for library in libraries:
        demux = demux_by_library[library]
        for barcode in sorted(demux, key=natural_key):
            preliminary = preliminary_by_cell[(library, barcode)]
            demux_raw = clean(demux[barcode].get("assignment"))
            original = canonical_genotype(
                preliminary.get("original_demux_assignment", ""))
            if original != canonical_genotype(demux_raw):
                raise ValueError(
                    f"{library}/{barcode} immutable demux assignment mismatch")
            refined = canonical_genotype(
                preliminary.get("current_refined_assignment", ""))
            current = refined or original
            proposal = canonical_genotype(
                preliminary.get("proposed_donor_genotype", ""))
            event_id = clean(preliminary.get("event_id"))
            action = clean(preliminary.get("final_action")).upper()
            prelim_applied = _final_bool(preliminary.get("reassignment_applied"))
            prelim_production = canonical_genotype(
                preliminary.get("reconciled_donor_genotype", "")) or current
            active_proposal = bool(proposal and proposal != current)

            available_axis = axis_by_cell.get((library, barcode), [])
            matching_axis = [
                row for row in available_axis
                if (not event_id or clean(row.get("selected_supported_event_id")) == event_id)
                and (not proposal or canonical_genotype(row.get("candidate_b", "")) == proposal)
            ]
            axis = sorted(
                matching_axis or available_axis,
                key=lambda row: natural_key(clean(row.get("score_pair_id"))))[0] \
                if (matching_axis or available_axis) else {}
            available_exclusions = exclusion_by_cell.get((library, barcode), [])
            matching_exclusions = [
                row for row in available_exclusions
                if (not event_id or clean(row.get("selected_supported_event_id")) == event_id)
                and (not proposal or canonical_genotype(
                    row.get("selected_supported_event_proposal", "")) == proposal)
            ]
            exclusion = sorted(
                matching_exclusions,
                key=lambda row: natural_key(clean(row.get("exclusion_reason"))))[0] \
                if matching_exclusions else {}
            selected_axis[(library, barcode)] = axis
            selected_exclusion[(library, barcode)] = exclusion
            nuclear_status, nuclear_warnings = _final_axis_status(
                axis, exclusion, current, proposal,
                active_proposal and bool(event_id),
            )
            if len(available_axis) > 1:
                nuclear_warnings.add("MULTIPLE_CANDIDATE_AXIS_ROWS")
            if demux_warnings[library].get(barcode):
                nuclear_warnings.add(demux_warnings[library][barcode])

            axis_a = canonical_genotype(axis.get("candidate_a", ""))
            axis_b = canonical_genotype(axis.get("candidate_b", ""))
            if axis_a or axis_b:
                relationship = ";".join((
                    "A=CURRENT" if axis_a == current else
                    "A=ORIGINAL_DEMUX" if axis_a == original else "A=OTHER",
                    "B=PROPOSAL" if axis_b == proposal else "B=OTHER",
                ))
            else:
                relationship = "NOT_APPLICABLE" if not active_proposal else "NO_EMITTED_PAIR"

            frozen_candidates = controlled.get((library, barcode), [])
            frozen_matching = [
                row for row in frozen_candidates
                if (not event_id or clean(row.get("event_id")) == event_id)
                and (not proposal or canonical_genotype(row.get("candidate_identity", "")) == proposal)
            ]
            frozen = sorted(
                frozen_matching,
                key=lambda row: natural_key(clean(row.get("condition"))))[0] \
                if frozen_matching else {}
            frozen_plan_key = (library, barcode, event_id, proposal)
            if not active_proposal:
                frozen_status = "NOT_APPLICABLE"
            elif frozen:
                frozen_current = canonical_genotype(frozen.get("original_identity", ""))
                frozen_proposal = canonical_genotype(frozen.get("candidate_identity", ""))
                frozen_status = (
                    "PAIRED_FIXED_PROFILE_FIT_NOT_EMITTED"
                    if frozen_current == current and frozen_proposal == proposal else
                    "PAIRED_ENDPOINTS_DO_NOT_MATCH_CURRENT_COMPARISON"
                )
            elif frozen_plan_key in controlled_planned:
                frozen_status = "PLANNED_NOT_EMITTED"
            else:
                frozen_status = "NOT_PLANNED_OR_NOT_APPLICABLE"
            frozen_fields = {
                "ambient_frozen_current_c": _final_na(frozen.get("fixed_original_c")),
                "ambient_frozen_proposal_c": _final_na(frozen.get("fixed_proposed_c")),
                "ambient_frozen_proposal_minus_current_c": _final_na(frozen.get("fixed_delta_c")),
                "ambient_frozen_current_fit": "NA",
                "ambient_frozen_proposal_fit": "NA",
                "ambient_frozen_fit_delta": "NA",
                "ambient_frozen_current_exact_candidate_burden": "NA",
                "ambient_frozen_proposal_exact_candidate_burden": "NA",
                "ambient_frozen_profile_id": frozen_profile_by_library.get(library, "NA"),
                "ambient_frozen_evaluation_status": frozen_status,
            }
            four_fields = _final_four_arm_fields(
                library, barcode, library in event_bearing, four_pivot,
                four_strata, four_summaries, four_burdens,
            )

            ploidy_status = _final_ploidy_status(
                preliminary, current, proposal)
            occupancy = clean(preliminary.get("occupancy_resolution_status")) or "NOT_APPLICABLE"
            technical_state = clean(preliminary.get("competing_technical_state"))
            if (not technical_state and
                    clean(preliminary.get("proposed_droplet_state")) ==
                    "TECHNICAL_MULTIPLET_CANDIDATE"):
                technical_state = (
                    clean(preliminary.get("proposed_droplet_constituents"))
                    or "TECHNICAL_MULTIPLET_CANDIDATE")
            technical_state = technical_state or "NOT_APPLICABLE"

            review_reasons = set()
            if active_proposal and nuclear_status in {
                    "NUCLEAR_SUPPORTS_CURRENT", "NUCLEAR_UNSTABLE",
                    "NUCLEAR_NONDISCRIMINATING", "NUCLEAR_PAIR_INADEQUATE",
                    "NUCLEAR_PAIR_DOES_NOT_MATCH_CURRENT_COMPARISON",
                    "NUCLEAR_OUTSIDE_PAIR_CONTRACT", "NUCLEAR_UNAVAILABLE"}:
                review_reasons.add(nuclear_status)
            if active_proposal and ploidy_status in {"UNAVAILABLE", "CONFLICTED", "SUPPORTS_CURRENT"}:
                review_reasons.add("PLOIDY_" + ploidy_status)
            if "UNRESOLVED" in occupancy or "AMBIG" in occupancy:
                review_reasons.add("OCCUPANCY_" + occupancy)
            if technical_state != "NOT_APPLICABLE" or _final_bool(
                    preliminary.get("explicit_multiplet_evidence")):
                review_reasons.add("TECHNICAL_MIXTURE_OR_MULTIPLET")
            mt_status = clean(preliminary.get("mt_alternative_status")) or "MISSING"
            if active_proposal and mt_status in {"SUPPORTS_CURRENT", "CONTRADICTS"}:
                review_reasons.add("MITOCHONDRIAL_" + mt_status)
            uid_status = clean(preliminary.get("uid_resolution_status")) or "UNAVAILABLE"
            if any(token in uid_status for token in ("CONFLICT", "MISSING")):
                review_reasons.add("METADATA_UID_" + uid_status)
            upstream_warning = demux_warnings[library].get(barcode)
            if upstream_warning:
                review_reasons.add(upstream_warning)
            event_class = clean(preliminary.get("event_class"))
            if active_proposal and (
                    "UNRESOLVED" in event_class
                    or clean(preliminary.get("event_confidence")) in {"INSUFFICIENT", "SUGGESTIVE"}):
                review_reasons.add("EVENT_MECHANISM_REQUIRES_REVIEW")
            library_exchange_status = _final_set(exchange_by_library.get(library, [])) or "NONE"
            if active_proposal and any(token in library_exchange_status for token in (
                    "PARTIAL_RECIPROCAL", "UNRESOLVED", "ROSTER_EQUIVALENT")):
                review_reasons.add("LIBRARY_EXCHANGE_AMBIGUITY")
            frozen_delta = _final_float(frozen_fields[
                "ambient_frozen_proposal_minus_current_c"])
            refit_delta = _final_float(four_fields[
                "ambient_assignment_effect_c_minus_b"])
            if (math.isfinite(frozen_delta) and math.isfinite(refit_delta)
                    and frozen_delta != 0 and refit_delta != 0
                    and (frozen_delta < 0) != (refit_delta < 0)):
                review_reasons.add(
                    "AMBIENT_CONTROLLED_REFITTED_DIRECTION_CONFLICT")

            if not active_proposal:
                recommendation = "NOT_APPLICABLE"
                interpreted = current
                recommendation_basis = "NO_ACTIVE_RECONCILIATION_PROPOSAL"
            elif prelim_applied or action in {"REASSIGN_GENOTYPE", "RECLASSIFY_PLOIDY"}:
                recommendation = "USE_PROPOSAL"
                interpreted = proposal
                recommendation_basis = "PRELIMINARY_RECONCILIATION_ACTION"
            elif action == "KEEP":
                recommendation = "USE_CURRENT"
                interpreted = current
                recommendation_basis = "PRELIMINARY_RECONCILIATION_KEEP"
            else:
                recommendation = "UNRESOLVED"
                interpreted = "UNRESOLVED"
                recommendation_basis = "PRELIMINARY_RECONCILIATION_REVIEW"

            cell_review_key = (library, barcode)
            event_review_key = (library, event_id) if event_id else None
            review = review_input.get(cell_review_key, {})
            review_record_scope = "CELL" if review else "NONE"
            review_record_target = barcode if review else "NA"
            review_key = cell_review_key if review else None
            if not review and event_review_key:
                review = review_input.get(event_review_key, {})
                if review:
                    review_record_scope = "EVENT"
                    review_record_target = event_id
                    review_key = event_review_key
            if review:
                used_review_inputs.add(review_key)
            disposition = review.get("disposition", "")
            rationale = review.get("rationale", "")
            production = prelim_production
            production_source = (
                "PRELIMINARY_RECONCILIATION_APPLIED"
                if prelim_applied else
                "REFINED_ASSIGNMENT" if refined else "DEMUX_ORIGINAL_ASSIGNMENT"
            )
            if disposition == "ACCEPT_PROPOSAL":
                if not active_proposal:
                    raise ValueError(
                        f"review ACCEPT_PROPOSAL has no active proposal: {library}/{barcode}")
                recommendation = "USE_PROPOSAL"
                interpreted = proposal
                production = proposal
                production_source = "EXPLICIT_REVIEW_ACCEPT_PROPOSAL"
                application_state = "APPLIED"
                application_reason = rationale or disposition
            elif disposition == "KEEP_CURRENT":
                recommendation = "USE_CURRENT"
                interpreted = current
                production = current
                production_source = "EXPLICIT_REVIEW_KEEP_CURRENT"
                application_state = "NOT_APPLIED"
                application_reason = rationale or disposition
            elif disposition == "LEAVE_UNRESOLVED":
                recommendation = "UNRESOLVED"
                interpreted = "UNRESOLVED"
                application_state = "HELD_FOR_REVIEW"
                application_reason = rationale or disposition
            elif not active_proposal:
                application_state = "NOT_APPLICABLE"
                application_reason = "NO_ACTIVE_RECONCILIATION_PROPOSAL"
            elif prelim_applied:
                application_state = "APPLIED"
                application_reason = "PRESERVED_PRELIMINARY_APPLICATION"
            elif review_reasons or action != "KEEP":
                application_state = "HELD_FOR_REVIEW"
                application_reason = "PRESERVED_PRELIMINARY_NONAPPLICATION"
            else:
                application_state = "NOT_APPLIED"
                application_reason = "PRESERVED_PRELIMINARY_CURRENT_STATE"

            review_required = bool(review_reasons) and not disposition
            current_droplet = clean(preliminary.get("reconciled_droplet_state"))
            if current_droplet == "TECHNICAL_MULTIPLET":
                release = "EXCLUDED_TECHNICAL_MULTIPLET"
                exclusion_reason = "TECHNICAL_MULTIPLET"
            elif disposition == "LEAVE_UNRESOLVED" or review_required:
                release = "HELD_FOR_REVIEW"
                exclusion_reason = _final_set(review_reasons) or "UNRESOLVED"
            elif not production:
                release = "EXCLUDED_NO_PRODUCTION_ASSIGNMENT"
                exclusion_reason = "NO_PRODUCTION_ASSIGNMENT"
            else:
                release = "READY"
                exclusion_reason = "NONE"

            proposal_kind = (
                "NOT_APPLICABLE" if not proposal else
                "TECHNICAL_MULTIPLET_HYPOTHESIS"
                if clean(preliminary.get("proposed_droplet_state")) ==
                "TECHNICAL_MULTIPLET_CANDIDATE" else
                "BIOLOGICAL_SINGLE_CELL_IDENTITY")
            row = {
                "library": library, "barcode": barcode,
                "demux_original_assignment": original,
                "demux_original_assignment_raw": demux_raw,
                "refined_assignment": refined or "NA",
                "comparison_current_assignment": current,
                "nominated_proposal": proposal or "NA", "event_id": event_id or "NA",
                "preliminary_reconciliation_action": action or "NA",
                "preliminary_action_applied": "TRUE" if prelim_applied else "FALSE",
                "preliminary_reconciled_assignment": prelim_production,
                "preliminary_decision_confidence": _final_na(
                    preliminary.get("decision_confidence")),
                "preliminary_decision_reason_codes": _final_na(
                    preliminary.get("decision_reason_codes")),
                "interpreted_identity": interpreted,
                "scientific_recommendation": recommendation,
                "recommendation_basis": recommendation_basis,
                "review_required": "TRUE" if review_required else "FALSE",
                "review_reasons": _final_set(review_reasons) or "NONE",
                "review_disposition": disposition or (
                    "PENDING" if review_required else "NONE"),
                "review_rationale": rationale or "NA",
                "review_record_scope": review_record_scope,
                "review_record_target": review_record_target,
                "application_state": application_state,
                "application_reason": application_reason,
                "production_assignment": production,
                "production_assignment_source": production_source,
                "candidate_roster_relationship": _final_na(
                    preliminary.get("singlet_library_relationship")
                    or preliminary.get("proposed_library_expected_status")),
                "proposal_kind": proposal_kind,
                "proposal_components": ",".join(donor_components(proposal)) or "NA",
                "axis_candidate_a_assignment": axis_a or "NA",
                "axis_candidate_b_assignment": axis_b or "NA",
                "axis_pair_relationship_to_current_proposal": relationship,
                "candidate_axis_scope_status": (
                    _final_na(axis.get("candidate_axis_status")) if axis else
                    "EXCLUDED" if exclusion else
                    "NOT_APPLICABLE" if not active_proposal else "UNAVAILABLE"),
                "candidate_axis_exclusion_reason": _final_na(
                    exclusion.get("exclusion_reason")),
                "nuclear_reconciliation_status": nuclear_status,
                "nuclear_warning_reasons": _final_set(nuclear_warnings) or "NONE",
                "current_donor_composition": ",".join(donor_components(current)),
                "proposal_donor_composition": ",".join(donor_components(proposal)) or "NA",
                "current_ploidy_state": _final_na(
                    preliminary.get("current_ploidy")
                    or _final_inferred_ploidy(current)),
                "proposal_ploidy_state": _final_na(
                    preliminary.get("proposed_biological_ploidy")
                    or (_final_inferred_ploidy(proposal) if proposal else "")),
                "ploidy_evidence_status": ploidy_status,
                "occupancy_state": _final_na(
                    preliminary.get("proposed_droplet_state")
                    or preliminary.get("reconciled_droplet_state")),
                "occupancy_evidence_status": occupancy,
                "known_line_relationship": _final_na("|".join(filter(None, (
                    clean(preliminary.get("proposed_global_biological_status")),
                    clean(preliminary.get("proposed_library_expected_status")))))),
                "technical_state": technical_state,
                "species_evidence_status": "CURRENT=" + _final_na(
                    preliminary.get("species_current_status")) + ";PROPOSAL=" +
                    _final_na(preliminary.get("species_alternative_status")),
                "mitochondrial_evidence_status": _final_na(mt_status),
                "mitochondrial_resolution_status": _final_na(
                    preliminary.get("mt_haplotype_resolution")
                    or preliminary.get("mt_fit_status")),
                "atac_evidence_status": (
                    "ATAC_NOT_REQUESTED" if args.evidence_mode == "rna" else
                    _final_na(preliminary.get("atac_status"))),
                **frozen_fields, **four_fields,
                "uid_resolution_status": uid_status,
                "uid_or_uid_set": _final_na(
                    preliminary.get("reconciled_uid")
                    or preliminary.get("uid_candidates")),
                "metadata_event_status": event_class or "NO_EVENT",
                "metadata_amendment_proposal": _final_set(
                    amendment_by_key.get((library, proposal), [])) or "NONE",
                "library_exchange_status": library_exchange_status,
                "downstream_release_status": release,
                "downstream_exclusion_reason": exclusion_reason,
                "event_identity_evidence_disposition": "NOT_DERIVED",
                "event_review_scope": "NOT_DERIVED",
                "review_scope": "NO_IMMEDIATE_REVIEW",
                "cell_exception_reasons": "NONE",
                "policy_version": _final_na(preliminary.get("policy_version")),
                "run_id": args.run_id or "NA",
                "final_schema_version": FINAL_SCHEMA_VERSION,
            }
            for field in axis_fields:
                if field not in row:
                    row[field] = axis.get(field, "NA") if axis else "NA"
            final_cells.append(row)

    unused_review_inputs = sorted(
        set(review_input) - used_review_inputs,
        key=lambda key: (natural_key(key[0]), natural_key(key[1])))
    if unused_review_inputs:
        raise ValueError(
            "review input contains targets that do not map to an exact selected "
            "library/barcode-or-event-id: "
            + ", ".join(f"{library}/{target}"
                        for library, target in unused_review_inputs[:10]))

    final_cell_index = {
        (row["library"], row["barcode"]): row for row in final_cells
    }
    if len(final_cell_index) != len(final_cells):
        raise ValueError("canonical final-cell ledger contains duplicate keys")

    audit_by_cell = defaultdict(list)
    for row in all_candidate_audit:
        audit_by_cell[(_final_library(row.get("library")), clean(row.get("barcode")))].append(row)
    for key in sorted(
            final_cell_index,
            key=lambda item: (natural_key(item[0]), natural_key(item[1]))):
        rows = audit_by_cell.get(key, [])
        if not rows:
            raise ValueError(
                f"candidate audit has no rows for selected cell: {key[0]}/{key[1]}")
        preliminary = preliminary_by_cell.get(key, {})
        current = canonical_genotype(
            preliminary.get("current_refined_assignment", "")
            or preliminary.get("original_demux_assignment", ""))
        proposal = canonical_genotype(preliminary.get("proposed_donor_genotype", ""))
        active = bool(proposal and proposal != current)
        axis = selected_axis.get(key, {})
        exclusion = selected_exclusion.get(key, {})
        endpoints = {
            "A": canonical_genotype(
                axis.get("candidate_a", "")
                or exclusion.get("original_demux_assignment", "")),
            "B": canonical_genotype(
                axis.get("candidate_b", "")
                or exclusion.get("selected_supported_event_proposal", "")
                or exclusion.get("proposed_donor_genotype", "")),
        }
        for role, endpoint in endpoints.items():
            if not endpoint or any(
                    canonical_genotype(row.get("candidate_canonical", "")) == endpoint
                    for row in rows):
                continue
            source_field = f"candidate_{role.lower()}_origin"
            raw_field = f"candidate_{role.lower()}"
            synthetic = {
                "library": key[0],
                "barcode": key[1],
                "event_id": clean(preliminary.get("event_id")) or "NA",
                "candidate_raw": clean(axis.get(raw_field)) or endpoint,
                "candidate_canonical": endpoint,
                "candidate_kind": (
                    "EXPECTED_LIBRARY_BIOLOGICAL_IDENTITY" if role == "A" else
                    "GLOBALLY_REAL_UNEXPECTED_BIOLOGICAL_IDENTITY"),
                "candidate_components": ",".join(donor_components(endpoint)),
                "candidate_source": (
                    clean(axis.get(source_field))
                    or f"CANDIDATE_AXIS_FIXED_ENDPOINT_{role}"),
                "candidate_tier": "CANDIDATE_AXIS_FIXED_ENDPOINT",
                "candidate_rank_within_source": f"CANDIDATE_AXIS:{role}",
                "candidate_eligibility": "SCORED" if axis else "SET_ASIDE",
                "candidate_set_aside_reason": (
                    "" if axis else
                    clean(exclusion.get("exclusion_reason"))
                    or "CANDIDATE_AXIS_PAIR_EXCLUDED"),
                "selected_as_proposal": "FALSE",
                "selected_for_candidate_axis": "FALSE",
                "axis_endpoint_role": "NA",
                "lower_rank_considered_after_set_aside": "FALSE",
                "candidate_priority": 10**9,
                "score_pair_id": clean(axis.get("score_pair_id")) or "NA",
                "schema_version": FINAL_SCHEMA_VERSION,
            }
            rows.append(synthetic)
            all_candidate_audit.append(synthetic)
        matching_proposal = [
            row for row in rows
            if canonical_genotype(row.get("candidate_canonical", "")) == proposal
        ]
        def audit_selection_key(row):
            return (
                clean(row.get("candidate_eligibility")) != "SCOREABLE",
                "TECHNICAL" in clean(row.get("candidate_kind")),
                fint(row.get("candidate_priority"), 10**9),
                natural_key(clean(row.get("candidate_source"))),
            )

        matching_proposal.sort(key=audit_selection_key)
        chosen_proposal = matching_proposal[0] if active and matching_proposal else None
        chosen_endpoint_roles = defaultdict(list)
        for role, endpoint in endpoints.items():
            matching_endpoint = [
                row for row in rows
                if endpoint and canonical_genotype(
                    row.get("candidate_canonical", "")) == endpoint
            ]
            matching_endpoint.sort(key=audit_selection_key)
            if matching_endpoint:
                chosen_endpoint = (
                    chosen_proposal
                    if endpoint == proposal and chosen_proposal in matching_endpoint
                    else matching_endpoint[0]
                )
                chosen_endpoint_roles[id(chosen_endpoint)].append(role)
        set_aside_priorities = [
            fint(row.get("candidate_priority"), 10**9) for row in rows
            if clean(row.get("candidate_eligibility")) == "SET_ASIDE"
        ]
        for row in rows:
            row["event_id"] = clean(preliminary.get("event_id")) or "NA"
            row["selected_as_proposal"] = "TRUE" if row is chosen_proposal else "FALSE"
            roles = chosen_endpoint_roles.get(id(row), [])
            row["selected_for_candidate_axis"] = "TRUE" if roles else "FALSE"
            row["axis_endpoint_role"] = ",".join(roles) or "NA"
            selected_priority = fint(row.get("candidate_priority"), 10**9)
            row["lower_rank_considered_after_set_aside"] = (
                "TRUE" if row is chosen_proposal and any(
                    priority < selected_priority for priority in set_aside_priorities)
                else "FALSE"
            )
        if active and chosen_proposal is None:
            raise ValueError(
                f"selected proposal has no candidate-audit row: {key[0]}/{key[1]} {proposal}")

    all_candidate_audit.sort(key=lambda row: (
        natural_key(_final_library(row.get("library"))),
        natural_key(clean(row.get("barcode"))),
        fint(row.get("candidate_priority"), 10**9),
        natural_key(clean(row.get("candidate_canonical"))),
        natural_key(clean(row.get("candidate_kind"))),
    ))

    event_extra_fields = [
        "initial_supporting_barcodes", "initial_supporting_count_from_cells",
        "initial_supporting_count_from_event", "interpreted_barcodes",
        "interpreted_barcode_count", "movement_reason_counts",
        "original_assignment_strata", "comparison_current_assignment_strata",
        "proposal_composition", "nuclear_status_distribution",
        "nuclear_stability_distribution", "ploidy_status_distribution",
        "occupancy_status_distribution", "library_exchange_status_final",
        "ambient_frozen_delta_median", "ambient_assignment_effect_c_minus_b_median",
        "review_required_cells", "review_reason_counts", "event_review_state",
        "n_nuclear_proposal_direction", "n_nuclear_current_direction",
        "n_nuclear_nondirectional", "nuclear_event_pattern",
        "controlled_ambient_pattern", "refitted_ambient_pattern",
        "identity_evidence_disposition", "evidence_modifiers",
        "mechanism_disposition", "review_scope", "phase3_review_reasons",
        "production_assignment_effect",
        "final_schema_version",
    ]
    voting_rows_by_event = defaultdict(list)
    for row in final_cells:
        if clean(row.get("population_votes_in_authoritative_event")).upper() != "TRUE":
            continue
        key = (
            row["library"], clean(row.get("selected_supported_event_id")),
            canonical_genotype(row.get("selected_supported_event_proposal", "")),
        )
        if key in final_event_keys:
            voting_rows_by_event[key].append(row)
    exclusions_by_event = defaultdict(list)
    for row in exclusion_rows:
        key = (
            _final_library(row.get("library")),
            clean(row.get("selected_supported_event_id")),
            canonical_genotype(row.get("selected_supported_event_proposal", "")),
        )
        if key in final_event_keys:
            exclusions_by_event[key].append(row)
    finalized_candidate_axis_event_keys = set(voting_rows_by_event)
    final_events = []
    for event in event_rows:
        event_id = clean(event.get("event_id"))
        event_library = _final_library(event.get("library"))
        event_proposal = canonical_genotype(
            event.get("unexpected_component", ""))
        event_key = (event_library, event_id, event_proposal)
        linked = [
            row for row in final_cells
            if row["library"] == event_library and row["event_id"] == event_id
        ]
        interpreted = [
            row for row in linked if row["scientific_recommendation"]
            in {"USE_PROPOSAL", "UNRESOLVED", "NO_IDENTITY_PREFERENCE"}
        ]
        movement = []
        for row in linked:
            movement.append({
                "USE_PROPOSAL": "INTERPRETED_AS_PROPOSAL",
                "USE_CURRENT": "RETAINED_CURRENT",
                "UNRESOLVED": "INTERPRETATION_UNRESOLVED",
                "NO_IDENTITY_PREFERENCE": "NO_IDENTITY_PREFERENCE",
                "NOT_APPLICABLE": "NO_ACTIVE_CELL_PROPOSAL",
            }.get(row["scientific_recommendation"], "OTHER"))
        frozen_values = [
            _final_float(row["ambient_frozen_proposal_minus_current_c"])
            for row in linked
        ]
        frozen_values = [value for value in frozen_values if math.isfinite(value)]
        refit_values = [
            _final_float(row["ambient_assignment_effect_c_minus_b"])
            for row in linked
        ]
        refit_values = [value for value in refit_values if math.isfinite(value)]
        applicable = bool(
            voting_rows_by_event.get(event_key)
            or any(
                canonical_genotype(row.get("nominated_proposal", ""))
                and canonical_genotype(row.get("nominated_proposal", ""))
                != canonical_genotype(row.get("comparison_current_assignment", ""))
                for row in linked
            )
        )
        (nuclear_pattern, nuclear_proposal_count, nuclear_current_count,
         nuclear_nondirectional_count) = _final_nuclear_event_pattern(
            voting_rows_by_event.get(event_key, []),
            exclusions_by_event.get(event_key, []), applicable,
        )
        frozen_median = median(frozen_values) if frozen_values else "NA"
        refit_median = median(refit_values) if refit_values else "NA"
        controlled_pattern = _final_ambient_pattern(
            frozen_median, applicable, "CONTROLLED")
        refitted_pattern = _final_ambient_pattern(
            refit_median, applicable, "REFIT")
        identity_disposition = _final_identity_evidence_disposition(
            event.get("event_class"), applicable, nuclear_pattern,
            controlled_pattern,
        )
        mechanism_disposition = _final_mechanism_disposition(
            event, linked, applicable)
        modifiers = _final_evidence_modifiers(
            event, linked, voting_rows_by_event.get(event_key, []),
            nuclear_pattern, controlled_pattern, refitted_pattern,
            mechanism_disposition,
        )
        production_changes = sum(
            row["production_assignment"]
            != row["preliminary_reconciled_assignment"] for row in linked
        )
        production_effect = (
            f"EXPLICIT_REVIEW_CHANGED_CELLS:{production_changes}"
            if production_changes else "NO_CHANGE_FROM_PHASE2_PRELIMINARY"
        )
        augmented = dict(event)
        augmented.update({
            "initial_supporting_barcodes": _final_set(row["barcode"] for row in linked) or "NOT_EMITTED",
            "initial_supporting_count_from_cells": len(linked),
            "initial_supporting_count_from_event": _final_na(event.get("n_implicated_cells")),
            "interpreted_barcodes": _final_set(row["barcode"] for row in interpreted) or "NONE",
            "interpreted_barcode_count": len(interpreted),
            "movement_reason_counts": _final_counter(movement),
            "original_assignment_strata": _final_counter(
                row["demux_original_assignment"] for row in linked),
            "comparison_current_assignment_strata": _final_counter(
                row["comparison_current_assignment"] for row in linked),
            "proposal_composition": _final_counter(
                row["nominated_proposal"] for row in linked),
            "nuclear_status_distribution": _final_counter(
                row["nuclear_reconciliation_status"] for row in linked),
            "nuclear_stability_distribution": _final_counter(
                row.get("candidate_axis_fold_direction_stability_status", "NA")
                for row in linked),
            "ploidy_status_distribution": _final_counter(
                row["ploidy_evidence_status"] for row in linked),
            "occupancy_status_distribution": _final_counter(
                row["occupancy_evidence_status"] for row in linked),
            "library_exchange_status_final": _final_set(
                row["library_exchange_status"] for row in linked) or "NONE",
            "ambient_frozen_delta_median": frozen_median,
            "ambient_assignment_effect_c_minus_b_median": refit_median,
            "review_required_cells": sum(
                row["review_required"] == "TRUE" for row in linked),
            "review_reason_counts": _final_counter(
                reason for row in linked
                for reason in row["review_reasons"].split(";")
                if reason and reason != "NONE"),
            "event_review_state": (
                "REVIEW_REQUIRED" if any(
                    row["review_required"] == "TRUE" for row in linked)
                else "NO_PENDING_CELL_REVIEW"),
            "n_nuclear_proposal_direction": nuclear_proposal_count,
            "n_nuclear_current_direction": nuclear_current_count,
            "n_nuclear_nondirectional": nuclear_nondirectional_count,
            "nuclear_event_pattern": nuclear_pattern,
            "controlled_ambient_pattern": controlled_pattern,
            "refitted_ambient_pattern": refitted_pattern,
            "identity_evidence_disposition": identity_disposition,
            "evidence_modifiers": modifiers,
            "mechanism_disposition": mechanism_disposition,
            "production_assignment_effect": production_effect,
            "final_schema_version": FINAL_SCHEMA_VERSION,
        })
        event_review_scope, phase3_review_reasons = _final_event_review_scope(
            augmented, identity_disposition, mechanism_disposition, modifiers,
            (event_library, event_id) in review_input,
        )
        augmented["review_scope"] = event_review_scope
        augmented["phase3_review_reasons"] = phase3_review_reasons
        final_events.append(augmented)

    event_disposition_by_id = {
        (_final_library(row.get("library")), clean(row.get("event_id"))): row
        for row in final_events
    }
    for row in final_cells:
        event = event_disposition_by_id.get((row["library"], row["event_id"]))
        row["event_identity_evidence_disposition"] = (
            event["identity_evidence_disposition"] if event else
            "NOT_APPLICABLE"
        )
        row["event_review_scope"] = (
            event["review_scope"] if event else "NO_IMMEDIATE_REVIEW"
        )
        cell_reasons = _final_cell_exception_reasons(row)
        row["cell_exception_reasons"] = cell_reasons
        row["review_scope"] = (
            "CELL_EXCEPTION_REVIEW" if cell_reasons != "NONE" else
            "NO_IMMEDIATE_REVIEW"
        )

    run_summary_fields = [
        "library", "input_barcodes", "output_ledger_rows", "event_count",
        "candidate_axis_planned", "candidate_axis_scored",
        "candidate_axis_excluded", "candidate_axis_unavailable",
        "candidate_axis_reason_counts", "candidate_axis_frozen_match_status",
        "ambient_frozen_planned", "ambient_frozen_emitted",
        "ambient_frozen_paired", "ambient_frozen_missing",
        "ambient_four_arm_comparison_counts", "evidence_status_counts",
        "review_required_cells", "review_reason_counts",
        "identity_disposition_counts", "mechanism_disposition_counts",
        "event_review_scope_counts", "cell_review_scope_counts",
        "actionable_event_review_rows", "actionable_cell_exception_rows",
        "review_records_present", "review_records_applied",
        "changes_demux_to_refined", "changes_current_to_preliminary_production",
        "changes_preliminary_to_final_production", "zero_event_status",
        "atac_evidence_mode", "accounting_status", "warnings",
        "final_schema_version",
    ]
    run_summary = []
    for library in libraries:
        rows = [row for row in final_cells if row["library"] == library]
        library_events = [
            row for row in final_events
            if _final_library(row.get("library")) == library
        ]
        active = [
            row for row in rows
            if row["nominated_proposal"] not in {"NA", row["comparison_current_assignment"]}
        ]
        scored_keys = {
            (_final_library(row.get("library")), clean(row.get("barcode")))
            for row in axis_rows if _final_library(row.get("library")) == library
        }
        excluded_keys = {
            (_final_library(row.get("library")), clean(row.get("barcode")))
            for row in exclusion_rows if _final_library(row.get("library")) == library
        }
        planned_keys = scored_keys | excluded_keys
        frozen_planned_keys = {
            key for key in controlled_planned if key[0] == library}
        frozen_emitted_keys = {
            (
                key[0], key[1], clean(row.get("event_id")),
                canonical_genotype(row.get("candidate_identity", "")),
            )
            for key, values in controlled.items() if key[0] == library
            for row in values
            if clean(row.get("event_id"))
            and canonical_genotype(row.get("candidate_identity", ""))
        }
        four_counts = Counter()
        for (lib, _), comparisons in four_pivot.items():
            if lib == library:
                four_counts.update(comparisons.keys())
        review_reason_values = [
            reason for row in rows for reason in row["review_reasons"].split(";")
            if reason and reason != "NONE"
        ]
        warnings = []
        if demux_warning_provenance[library] != "AVAILABLE":
            warnings.append("DEMUX_MAXIMIN_PROVENANCE_UNAVAILABLE")
        run_summary.append({
            "library": library,
            "input_barcodes": len(demux_by_library[library]),
            "output_ledger_rows": len(rows),
            "event_count": sum(
                _final_library(event.get("library")) == library
                for event in event_rows),
            "candidate_axis_planned": len(planned_keys),
            "candidate_axis_scored": len(scored_keys),
            "candidate_axis_excluded": len(excluded_keys),
            "candidate_axis_unavailable": sum(
                row["nuclear_reconciliation_status"] == "NUCLEAR_UNAVAILABLE"
                for row in active),
            "candidate_axis_reason_counts": _final_counter(
                row.get("exclusion_reason") for row in exclusion_rows
                if _final_library(row.get("library")) == library),
            "candidate_axis_frozen_match_status": (
                "RAW_FIELDS_COPIED_VERBATIM"
                if scored_keys else
                "NOT_APPLICABLE_ZERO_EVENT" if library not in event_bearing else
                "NO_SCORE_ROWS_EMITTED"),
            "ambient_frozen_planned": len(frozen_planned_keys),
            "ambient_frozen_emitted": len(frozen_emitted_keys),
            "ambient_frozen_paired": len(frozen_planned_keys & frozen_emitted_keys),
            "ambient_frozen_missing": len(frozen_planned_keys - frozen_emitted_keys),
            "ambient_four_arm_comparison_counts": ";".join(
                f"{key}:{four_counts[key]}" for key in sorted(four_counts)) or "NONE",
            "evidence_status_counts": _final_counter(
                row["nuclear_reconciliation_status"] for row in rows),
            "review_required_cells": sum(
                row["review_required"] == "TRUE" for row in rows),
            "review_reason_counts": _final_counter(review_reason_values),
            "identity_disposition_counts": _final_counter(
                row["identity_evidence_disposition"] for row in library_events),
            "mechanism_disposition_counts": _final_counter(
                row["mechanism_disposition"] for row in library_events),
            "event_review_scope_counts": _final_counter(
                row["review_scope"] for row in library_events),
            "cell_review_scope_counts": _final_counter(
                row["review_scope"] for row in rows),
            "actionable_event_review_rows": sum(
                row["review_scope"] == "EVENT_REVIEW"
                for row in library_events),
            "actionable_cell_exception_rows": sum(
                row["review_scope"] == "CELL_EXCEPTION_REVIEW"
                for row in rows),
            "review_records_present": sum(
                key[0] == library for key in review_input),
            "review_records_applied": sum(
                key[0] == library for key in used_review_inputs),
            "changes_demux_to_refined": sum(
                row["refined_assignment"] != "NA"
                and row["refined_assignment"] != row["demux_original_assignment"]
                for row in rows),
            "changes_current_to_preliminary_production": sum(
                row["preliminary_reconciled_assignment"] !=
                row["comparison_current_assignment"] for row in rows),
            "changes_preliminary_to_final_production": sum(
                row["production_assignment"] !=
                row["preliminary_reconciled_assignment"] for row in rows),
            "zero_event_status": (
                "NOT_APPLICABLE_ZERO_EVENT_SUCCESS" if library not in event_bearing
                else "EVENT_BEARING"),
            "atac_evidence_mode": (
                "ATAC_NOT_REQUESTED" if args.evidence_mode == "rna" else "RNA_ATAC"),
            "accounting_status": "PASS",
            "warnings": _final_set(warnings) or "NONE",
            "final_schema_version": FINAL_SCHEMA_VERSION,
        })

    overall_reasons = [
        reason for row in final_cells for reason in row["review_reasons"].split(";")
        if reason and reason != "NONE"
    ]
    run_summary.append({
        "library": "ALL",
        "input_barcodes": sum(len(value) for value in demux_by_library.values()),
        "output_ledger_rows": len(final_cells),
        "event_count": len(event_rows),
        "candidate_axis_planned": sum(int(row["candidate_axis_planned"]) for row in run_summary),
        "candidate_axis_scored": sum(int(row["candidate_axis_scored"]) for row in run_summary),
        "candidate_axis_excluded": sum(int(row["candidate_axis_excluded"]) for row in run_summary),
        "candidate_axis_unavailable": sum(int(row["candidate_axis_unavailable"]) for row in run_summary),
        "candidate_axis_reason_counts": _final_counter(
            row.get("exclusion_reason") for row in exclusion_rows),
        "candidate_axis_frozen_match_status": (
            "RAW_FIELDS_COPIED_VERBATIM" if axis_rows else "NO_SCORE_ROWS_EMITTED"),
        "ambient_frozen_planned": sum(int(row["ambient_frozen_planned"]) for row in run_summary),
        "ambient_frozen_emitted": sum(int(row["ambient_frozen_emitted"]) for row in run_summary),
        "ambient_frozen_paired": sum(int(row["ambient_frozen_paired"]) for row in run_summary),
        "ambient_frozen_missing": sum(int(row["ambient_frozen_missing"]) for row in run_summary),
        "ambient_four_arm_comparison_counts": _final_counter(
            contrast for (library, _), comparisons in four_pivot.items()
            if library in selected_libraries for contrast in comparisons),
        "evidence_status_counts": _final_counter(
            row["nuclear_reconciliation_status"] for row in final_cells),
        "review_required_cells": sum(
            row["review_required"] == "TRUE" for row in final_cells),
        "review_reason_counts": _final_counter(overall_reasons),
        "identity_disposition_counts": _final_counter(
            row["identity_evidence_disposition"] for row in final_events),
        "mechanism_disposition_counts": _final_counter(
            row["mechanism_disposition"] for row in final_events),
        "event_review_scope_counts": _final_counter(
            row["review_scope"] for row in final_events),
        "cell_review_scope_counts": _final_counter(
            row["review_scope"] for row in final_cells),
        "actionable_event_review_rows": sum(
            row["review_scope"] == "EVENT_REVIEW" for row in final_events),
        "actionable_cell_exception_rows": sum(
            row["review_scope"] == "CELL_EXCEPTION_REVIEW"
            for row in final_cells),
        "review_records_present": len(review_input),
        "review_records_applied": len(used_review_inputs),
        "changes_demux_to_refined": sum(int(row["changes_demux_to_refined"]) for row in run_summary),
        "changes_current_to_preliminary_production": sum(
            int(row["changes_current_to_preliminary_production"]) for row in run_summary),
        "changes_preliminary_to_final_production": sum(
            int(row["changes_preliminary_to_final_production"]) for row in run_summary),
        "zero_event_status": f"ZERO_EVENT_LIBRARIES={sum(row['zero_event_status'].startswith('NOT_APPLICABLE') for row in run_summary)}",
        "atac_evidence_mode": (
            "ATAC_NOT_REQUESTED" if args.evidence_mode == "rna" else "RNA_ATAC"),
        "accounting_status": "PASS",
        "warnings": _final_set(
            row["warnings"] for row in run_summary if row["warnings"] != "NONE") or "NONE",
        "final_schema_version": FINAL_SCHEMA_VERSION,
    })

    final_cells.sort(key=lambda row: (
        natural_key(row["library"]), natural_key(row["barcode"])))
    final_events.sort(key=lambda row: (
        natural_key(_final_library(row.get("library"))),
        natural_key(clean(row.get("event_id")))))
    phase3_assignment_changes = [
        row for row in final_cells
        if row["production_assignment"] != row["preliminary_reconciled_assignment"]
    ]
    if not review_input and phase3_assignment_changes:
        raise ValueError(
            "production assignments changed without an explicit review record")
    untraceable_changes = [
        row for row in phase3_assignment_changes
        if row["review_record_scope"] not in {"CELL", "EVENT"}
    ]
    if untraceable_changes:
        raise ValueError(
            "production assignment change is not traceable to an explicit "
            "review record: "
            + ", ".join(
                f"{row['library']}/{row['barcode']}"
                for row in untraceable_changes[:10]
            )
        )

    review_queue = [row for row in final_cells if row["review_required"] == "TRUE"]
    event_dispositions = []
    for event in final_events:
        event_dispositions.append({
            "library": _final_library(event.get("library")),
            "event_id": clean(event.get("event_id")) or "NA",
            "event_class": _final_na(event.get("event_class")),
            "event_confidence": _final_na(event.get("event_confidence")),
            "n_initial_cells": _final_na(
                event.get("initial_supporting_count_from_event")),
            "n_interpreted_cells": _final_na(
                event.get("interpreted_barcode_count")),
            "n_nuclear_proposal_direction": event[
                "n_nuclear_proposal_direction"],
            "n_nuclear_current_direction": event[
                "n_nuclear_current_direction"],
            "n_nuclear_nondirectional": event["n_nuclear_nondirectional"],
            "nuclear_event_pattern": event["nuclear_event_pattern"],
            "controlled_ambient_pattern": event["controlled_ambient_pattern"],
            "refitted_ambient_pattern": event["refitted_ambient_pattern"],
            "identity_evidence_disposition": event[
                "identity_evidence_disposition"],
            "evidence_modifiers": event["evidence_modifiers"],
            "mechanism_disposition": event["mechanism_disposition"],
            "review_scope": event["review_scope"],
            "review_reasons": event["phase3_review_reasons"],
            "production_assignment_effect": event[
                "production_assignment_effect"],
            "final_schema_version": FINAL_SCHEMA_VERSION,
        })

    actionable_review = []
    for event in final_events:
        if event["review_scope"] != "EVENT_REVIEW":
            continue
        key = (_final_library(event.get("library")), clean(event.get("event_id")))
        explicit = review_input.get(key, {})
        actionable_review.append({
            "review_level": "EVENT_REVIEW",
            "library": key[0],
            "event_id": key[1] or "NA",
            "barcode": "NA",
            "identity_evidence_disposition": event[
                "identity_evidence_disposition"],
            "mechanism_disposition": event["mechanism_disposition"],
            "review_scope": "EVENT_REVIEW",
            "review_reasons": event["phase3_review_reasons"],
            "review_disposition": explicit.get("disposition", "PENDING"),
            "review_rationale": explicit.get("rationale") or "NA",
            "comparison_current_assignment": "EVENT_LEVEL",
            "nominated_proposal": canonical_genotype(
                event.get("unexpected_component", "")) or "NA",
            "production_assignment": "EVENT_LEVEL",
            "production_assignment_effect": event[
                "production_assignment_effect"],
            "final_schema_version": FINAL_SCHEMA_VERSION,
        })
    for row in final_cells:
        if row["review_scope"] != "CELL_EXCEPTION_REVIEW":
            continue
        actionable_review.append({
            "review_level": "CELL_EXCEPTION_REVIEW",
            "library": row["library"],
            "event_id": row["event_id"],
            "barcode": row["barcode"],
            "identity_evidence_disposition": row[
                "event_identity_evidence_disposition"],
            "mechanism_disposition": (
                event_disposition_by_id.get(
                    (row["library"], row["event_id"]), {}
                ).get("mechanism_disposition", "NOT_APPLICABLE")
            ),
            "review_scope": "CELL_EXCEPTION_REVIEW",
            "review_reasons": row["cell_exception_reasons"],
            "review_disposition": row["review_disposition"],
            "review_rationale": row["review_rationale"],
            "comparison_current_assignment": row[
                "comparison_current_assignment"],
            "nominated_proposal": row["nominated_proposal"],
            "production_assignment": row["production_assignment"],
            "production_assignment_effect": (
                "CHANGED_FROM_PHASE2_PRELIMINARY"
                if row["production_assignment"]
                != row["preliminary_reconciled_assignment"] else
                "NO_CHANGE_FROM_PHASE2_PRELIMINARY"
            ),
            "final_schema_version": FINAL_SCHEMA_VERSION,
        })
    actionable_review.sort(key=lambda row: (
        0 if row["review_level"] == "EVENT_REVIEW" else 1,
        natural_key(row["library"]), natural_key(row["event_id"]),
        natural_key(row["barcode"]),
    ))

    nonfinalized_worsened = _final_nonfinalized_controlled_worsened(
        controlled, controlled_planned, selected_libraries,
        finalized_candidate_axis_event_keys,
    )

    final_cell_fields = FINAL_CELL_FIELDS + [
        field for field in axis_fields if field not in FINAL_CELL_FIELDS
    ]
    event_fields = RECONCILE_EVENT_FIELDS + [
        field for field in event_extra_fields if field not in RECONCILE_EVENT_FIELDS
    ]
    for row in all_candidate_audit:
        for field in row:
            if field not in candidate_audit_fields:
                candidate_audit_fields.append(field)
    candidate_fields = CANDIDATE_AUDIT_REQUIRED_FIELDS + [
        field for field in candidate_audit_fields
        if field not in CANDIDATE_AUDIT_REQUIRED_FIELDS
    ]
    write_tsv(
        str(output_root / "identity_reconciliation_final_cells.tsv.gz"),
        final_cells, final_cell_fields)
    write_tsv(
        str(output_root / "identity_reconciliation_candidate_audit.tsv.gz"),
        all_candidate_audit, candidate_fields)
    write_tsv(
        str(output_root / "identity_reconciliation_final_events.tsv"),
        final_events, event_fields)
    write_tsv(
        str(output_root / "identity_reconciliation_review_queue.tsv.gz"),
        review_queue, final_cell_fields)
    write_tsv(
        str(output_root / "identity_reconciliation_run_summary.tsv"),
        run_summary, run_summary_fields)
    write_tsv(
        str(output_root / "identity_reconciliation_event_dispositions.tsv"),
        event_dispositions, FINAL_EVENT_DISPOSITION_FIELDS)
    write_tsv(
        str(output_root / "identity_reconciliation_actionable_review.tsv.gz"),
        actionable_review, FINAL_ACTIONABLE_REVIEW_FIELDS)
    (output_root / "identity_reconciliation_phase3_summary.md").write_text(
        _final_phase3_summary_text(
            libraries, final_cells, final_events, actionable_review,
            review_input, used_review_inputs, args.evidence_mode,
            nonfinalized_worsened,
        ),
        encoding="utf-8",
    )

    assignments_root = output_root.parent / "final_assignments"
    assignments_root.mkdir(parents=True, exist_ok=True)
    for library in libraries:
        rows = [row for row in final_cells if row["library"] == library]
        demux = demux_by_library[library]
        assignments = []
        for row in rows:
            identity = canonical_genotype(row["production_assignment"])
            if not identity:
                raise ValueError(
                    f"{library}/{row['barcode']} has no production assignment")
            assignment_type = (
                "D" if identity.startswith("M{")
                or len(donor_components(identity)) >= 2 else "S")
            score = clean(demux[row["barcode"]].get("score")) or "NA"
            assignments.append((row["barcode"], identity, assignment_type, score))
        if len(assignments) != len(demux):
            raise ValueError(f"{library} compatibility assignment accounting failed")
        write_headerless_tsv(
            str(assignments_root / f"{library}.reconciled.assignments"),
            assignments)
    print(
        f"Finalized {len(final_cells)} cells, {len(final_events)} events, "
        f"{len(actionable_review)} actionable Phase 3 review rows, and "
        f"{len(review_queue)} legacy pending review rows")
    return 0


_COMMANDS = {
    "metadata": metadata_main,
    "candidates": candidates_main,
    "doublet-context": doublet_context_main,
    "reconcile": reconcile_main,
    "optional-status": optional_status_main,
    "score-pairs": score_pairs_main,
    "candidate-axis-pairs": candidate_axis_pairs_main,
    "finalize": finalize_main,
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
