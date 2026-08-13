#!/usr/bin/env python3
"""Shared helpers for CellBouncer post-hoc QC / swap-audit scripts.

These helpers are intentionally dependency-light. They only use the Python
standard library so the audit layer can run in whatever environment launches the
pipeline.
"""
from __future__ import annotations

import csv
import gzip
import json
import math
import os
import re
import shutil
import sys
import time
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Mapping, Optional, Sequence, Set, Tuple

NA = "NA"


def now_iso() -> str:
    return time.strftime("%Y-%m-%dT%H:%M:%S%z")


def open_text(path: os.PathLike | str, mode: str = "rt"):
    p = str(path)
    if "b" in mode:
        raise ValueError("open_text only supports text modes")
    if p.endswith(".gz"):
        return gzip.open(p, mode)
    return open(p, mode, newline="")


def read_tsv(path: os.PathLike | str) -> List[Dict[str, str]]:
    with open_text(path, "rt") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write_tsv(path: os.PathLike | str, rows: Sequence[Mapping[str, object]], fields: Sequence[str]) -> None:
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open_text(path, "wt") as fh:
        w = csv.DictWriter(fh, fieldnames=list(fields), delimiter="\t", lineterminator="\n", extrasaction="ignore")
        w.writeheader()
        for row in rows:
            out = {f: stringify(row.get(f, "")) for f in fields}
            w.writerow(out)


def stringify(x: object) -> str:
    if x is None:
        return NA
    if isinstance(x, bool):
        return "true" if x else "false"
    if isinstance(x, float):
        if math.isnan(x):
            return NA
        return f"{x:.6g}"
    if isinstance(x, (list, tuple, set)):
        return ",".join(str(v) for v in x)
    return str(x)


def read_samples(path: os.PathLike | str) -> List[str]:
    samples: List[str] = []
    with open_text(path, "rt") as fh:
        for line in fh:
            s = line.strip().split()
            if s:
                samples.append(s[0])
    return samples


def canonical_identity(identity: str, samples: Optional[Set[str]] = None) -> str:
    """Validate and canonicalize one identity.

    Singlets are preserved. Homotypics are preserved as A+A. Heterotypics are
    sorted lexicographically so B+A and A+B become A+B. Unknown components fail
    when a sample set is provided.
    """
    raw = identity.strip()
    if not raw:
        raise ValueError("empty identity")
    raw = re.sub(r"\s+", "", raw)
    if raw.count("+") == 0:
        parts = [raw]
    elif raw.count("+") == 1:
        parts = raw.split("+")
    else:
        raise ValueError(f"malformed identity with more than one '+': {identity!r}")
    if any(not p for p in parts):
        raise ValueError(f"malformed identity with empty component: {identity!r}")
    if samples is not None:
        missing = [p for p in parts if p not in samples]
        if missing:
            raise ValueError(f"unknown identity component(s) {missing} in {identity!r}")
    if len(parts) == 1:
        return parts[0]
    if parts[0] == parts[1]:
        return f"{parts[0]}+{parts[1]}"
    a, b = sorted(parts)
    return f"{a}+{b}"


def snp_resolvable_identity(identity: str, samples: Optional[Set[str]] = None) -> str:
    """Return the identity that SNP evidence can actually distinguish.

    CellBouncer genotype evidence cannot distinguish a diploid/singlet A call
    from a homotypic tetraploid A+A call.  For audit composition comparisons,
    expected-vs-observed matching, and runner-up lookup, collapse A+A to A.

    This deliberately does *not* replace ``canonical_identity`` globally because
    downstream reporting may still want to preserve the biological/raw label
    A+A when it came from expected metadata or tetra_refine.
    """
    cid = canonical_identity(identity, samples)
    if "+" not in cid:
        return cid
    a, b = cid.split("+", 1)
    if a == b:
        return a
    return cid


def identity_components(identity: str) -> Tuple[str, ...]:
    cid = canonical_identity(identity, None)
    if "+" not in cid:
        return (cid,)
    a, b = cid.split("+", 1)
    if a == b:
        return (a,)
    return (a, b)


def identity_type(identity: str) -> str:
    cid = canonical_identity(identity, None)
    if "+" not in cid:
        return "singlet"
    a, b = cid.split("+", 1)
    return "homotypic" if a == b else "heterotypic"


def sort_identity_key(identity: str) -> Tuple[int, str]:
    t = identity_type(identity)
    rank = {"singlet": 0, "homotypic": 1, "heterotypic": 2}.get(t, 9)
    return (rank, identity)




# -----------------------------------------------------------------------------
# Refined assignment sidecar helpers
# -----------------------------------------------------------------------------

def read_refined_assignments(path: os.PathLike | str | None) -> Dict[str, Dict[str, str]]:
    """Read tetra_refine rich sidecar or demux-compatible refined assignment file.

    Preferred input is ``*.refined_assignments`` with a header containing at
    least ``barcode`` and ``refined_assignment``.  For robustness this also
    accepts the older/simple four-column assignment-like output:

        barcode  assignment  S/D  llr

    The return value is keyed by barcode and always contains:

        original_assignment, refined_assignment, refined_assignment_source,
        refined_assignment_confidence, droplet_doublet_flag,
        quad_pattern_score, changed

    Missing optional fields are populated with ``NA``.  The parser never treats
    the header row as a cell.
    """
    if not path:
        return {}
    p = Path(path)
    if not p.exists():
        return {}
    rows: Dict[str, Dict[str, str]] = {}
    with open_text(p, "rt") as fh:
        first = fh.readline()
        if not first:
            return {}
        first_parts = first.rstrip("\n\r").split("\t")
        has_header = "barcode" in first_parts and (
            "refined_assignment" in first_parts or "assignment" in first_parts
        )
        if has_header:
            header = first_parts
            lines = fh
        else:
            header = []
            lines = [first] + list(fh)

        if has_header:
            reader = csv.DictReader(lines, fieldnames=header, delimiter="\t")
            for r in reader:
                bc = (r.get("barcode") or "").strip()
                if not bc:
                    continue
                refined = (r.get("refined_assignment") or r.get("assignment") or "").strip()
                if not refined:
                    continue
                original = (r.get("original_assignment") or refined).strip()
                source = (r.get("ploidy_method") or r.get("refined_assignment_source") or "demux_only").strip() or "demux_only"
                confidence = (r.get("ploidy_confidence") or r.get("overall_confidence") or r.get("refined_assignment_confidence") or NA).strip() or NA
                droplet = (r.get("droplet_flag") or r.get("droplet_doublet_flag") or "NONE").strip() or "NONE"
                rows[bc] = {
                    "barcode": bc,
                    "original_assignment": original,
                    "refined_assignment": refined,
                    "refined_assignment_source": source,
                    "refined_assignment_confidence": confidence,
                    "ploidy": (r.get("ploidy") or NA).strip() or NA,
                    "ploidy_method": source,
                    "droplet_doublet_flag": droplet,
                    "droplet_candidates": (r.get("droplet_candidates") or NA).strip() or NA,
                    "quad_pattern_score": (r.get("quad_pattern_score") or NA).strip() or NA,
                    "changed": (r.get("changed") or str(original != refined)).strip() or str(original != refined),
                    "cells_in_droplet": (r.get("cells_in_droplet") or NA).strip() or NA,
                }
        else:
            for line in lines:
                if not line.strip():
                    continue
                parts = line.rstrip("\n\r").split("\t")
                if len(parts) < 2:
                    continue
                bc, refined = parts[0].strip(), parts[1].strip()
                if not bc or bc.lower() == "barcode" or not refined:
                    continue
                rows[bc] = {
                    "barcode": bc,
                    "original_assignment": refined,
                    "refined_assignment": refined,
                    "refined_assignment_source": "demux_compatible_refined_file",
                    "refined_assignment_confidence": NA,
                    "ploidy": NA,
                    "ploidy_method": "demux_compatible_refined_file",
                    "droplet_doublet_flag": "NONE",
                    "droplet_candidates": NA,
                    "quad_pattern_score": NA,
                    "changed": "NA",
                    "cells_in_droplet": NA,
                }
    return rows


def refined_assignment_counts(refined: Mapping[str, Mapping[str, str]]) -> Counter:
    """Count refined biological assignment classes from read_refined_assignments().

    Important terminology:
      * heterotypic/homotypic/plus_identity describe biological plus identities
        (A+B or A+A) in the refined label space.
      * droplet_flagged describes evidence for a possible multi-cell droplet from
        tetra_refine's dedicated droplet/quad-pattern QC.

    The legacy ``doublet`` counter is retained for backward compatibility, but
    consumers should prefer ``plus_identity`` / ``biological_fusion`` and
    ``droplet_flagged`` to avoid confusing a biological A+A or A+B label with a
    droplet doublet.
    """
    counts: Counter = Counter()
    for row in refined.values():
        ident = row.get("refined_assignment", "") or ""
        if not ident or ident == NA:
            continue
        t = identity_type(ident)
        if t == "singlet":
            counts["singlet"] += 1
        elif t in {"heterotypic", "homotypic"}:
            counts[t] += 1
            counts["plus_identity"] += 1
            counts["biological_fusion"] += 1
            # Backward-compatible legacy name. This does NOT imply droplet doublet.
            counts["doublet"] += 1
        else:
            counts[t] += 1
        if str(row.get("changed", "")).upper() in {"TRUE", "1", "YES"}:
            counts["changed"] += 1
        if (row.get("droplet_doublet_flag") or "NONE") not in {"", ".", "NA", "NONE"}:
            counts["droplet_flagged"] += 1
    counts["total"] = len(refined)
    return counts


def composition_distances_raw_identity(expected_ids: Sequence[str], observed_ids: Iterable[str]) -> Dict[str, float]:
    """Composition distances in biological/raw identity space.

    Unlike snp_resolvable_identity(), this preserves A+A as A+A.  Use it only
    for refined biological reporting, never for SNP-resolvable audit matching.
    """
    obs_counts = Counter(canonical_identity(x, None) for x in observed_ids if x and x != NA)
    total = sum(obs_counts.values()) or 1
    obs = {k: v / total for k, v in obs_counts.items()}
    exp_ids = sorted({canonical_identity(x, None) for x in expected_ids}, key=sort_identity_key)
    exp = {k: 1.0 / len(exp_ids) for k in exp_ids} if exp_ids else {}
    keys = sorted(set(obs) | set(exp), key=sort_identity_key)
    unexpected = sum(v for k, v in obs.items() if k not in exp)
    missing = sum(1 for k in exp if obs.get(k, 0.0) < 0.005) / len(exp) if exp else math.nan
    return {
        "unexpected_identity_fraction": unexpected,
        "missing_expected_identity_fraction": missing,
    }


def read_allowed_identities(path: os.PathLike | str, samples: Sequence[str]) -> Tuple[List[str], Dict[str, object]]:
    sample_set = set(samples)
    seen: Set[str] = set()
    out: List[str] = []
    n_input = 0
    n_comments_blank = 0
    duplicate_inputs: List[str] = []
    errors: List[str] = []
    with open_text(path, "rt") as fh:
        for lineno, line in enumerate(fh, 1):
            raw = line.strip()
            if not raw or raw.startswith("#"):
                n_comments_blank += 1
                continue
            n_input += 1
            # Permit accidental comma-separated lines by failing loudly rather
            # than silently treating them as one malformed sample name.
            if "," in raw:
                errors.append(f"line {lineno}: expected one identity per line, got comma(s): {raw!r}")
                continue
            try:
                cid = canonical_identity(raw, sample_set)
            except ValueError as e:
                errors.append(f"line {lineno}: {e}")
                continue
            if cid in seen:
                duplicate_inputs.append(raw)
                continue
            seen.add(cid)
            out.append(cid)
    if errors:
        raise ValueError("invalid allowed identities file:\n" + "\n".join(errors[:100]))
    out.sort(key=sort_identity_key)
    counts = Counter(identity_type(x) for x in out)
    warnings: List[str] = []
    if counts.get("singlet", 0) == 0:
        warnings.append("ALLOWED_IDENTITIES_NO_SINGLETS")
    if counts.get("homotypic", 0) == 0:
        warnings.append("ALLOWED_IDENTITIES_NO_HOMOTYPICS")
    if counts.get("heterotypic", 0) == 0:
        warnings.append("ALLOWED_IDENTITIES_NO_HETEROTYPICS")
    if duplicate_inputs:
        warnings.append("ALLOWED_IDENTITIES_DEDUPLICATED")
    meta = {
        "n_input_identity_lines": n_input,
        "n_comments_blank_lines": n_comments_blank,
        "n_unique_after_canonicalization": len(out),
        "n_duplicate_or_redundant_lines": len(duplicate_inputs),
        "n_allowed_singlets": counts.get("singlet", 0),
        "n_allowed_homotypics": counts.get("homotypic", 0),
        "n_allowed_heterotypics": counts.get("heterotypic", 0),
        "warnings": warnings,
    }
    return out, meta


def parse_expected_metadata(path: os.PathLike | str, samples: Optional[Sequence[str]] = None) -> Dict[str, List[str]]:
    sample_set = set(samples) if samples is not None else None
    rows = read_tsv(path)
    if not rows:
        raise ValueError(f"expected metadata is empty: {path}")
    header = set(rows[0].keys())
    lib_col = None
    ids_col = None
    for c in ("library_id", "library", "lib", "Lib"):
        if c in header:
            lib_col = c
            break
    for c in ("cell_identities", "expected_identities", "expected_identity", "identity", "identities"):
        if c in header:
            ids_col = c
            break
    if lib_col is None or ids_col is None:
        raise ValueError(f"expected metadata needs library_id/lib and cell_identities columns; saw {sorted(header)}")
    out: Dict[str, List[str]] = {}
    errors: List[str] = []
    for i, row in enumerate(rows, 2):
        lib = normalize_library(row.get(lib_col, ""))
        raw = row.get(ids_col, "") or ""
        vals = [x.strip() for x in raw.split(",") if x.strip()]
        can: List[str] = []
        for v in vals:
            try:
                can.append(canonical_identity(v, sample_set))
            except ValueError as e:
                errors.append(f"row {i} library {lib}: {e}")
        out[lib] = sorted(set(can), key=sort_identity_key)
    if errors:
        raise ValueError("invalid expected metadata:\n" + "\n".join(errors[:100]))
    return out


def build_allowed_from_expected(expected: Mapping[str, Sequence[str]], samples: Sequence[str], include_all_singlets: bool = True) -> Tuple[List[str], Dict[str, object]]:
    seen: Set[str] = set()
    if include_all_singlets:
        seen.update(samples)
    for vals in expected.values():
        for v in vals:
            seen.add(canonical_identity(v, set(samples)))
    out = sorted(seen, key=sort_identity_key)
    counts = Counter(identity_type(x) for x in out)
    meta = {
        "n_input_identity_lines": sum(len(v) for v in expected.values()),
        "n_comments_blank_lines": 0,
        "n_unique_after_canonicalization": len(out),
        "n_duplicate_or_redundant_lines": 0,
        "n_allowed_singlets": counts.get("singlet", 0),
        "n_allowed_homotypics": counts.get("homotypic", 0),
        "n_allowed_heterotypics": counts.get("heterotypic", 0),
        "warnings": [],
    }
    return out, meta


def normalize_library(x: object) -> str:
    s = str(x).strip()
    if not s:
        return s
    m = re.match(r"(?i)^lib(\d+)$", s)
    if m:
        return f"lib{int(m.group(1))}"
    if re.match(r"^\d+$", s):
        return f"lib{int(s)}"
    return s


def library_number(lib: str) -> Optional[int]:
    m = re.match(r"(?i)^lib(\d+)$", normalize_library(lib))
    return int(m.group(1)) if m else None


def parse_library_list(spec: str, expected: Optional[Mapping[str, Sequence[str]]] = None) -> List[str]:
    spec = (spec or "all").strip()
    if spec.lower() == "all":
        if expected:
            return sorted(expected.keys(), key=lambda x: (library_number(x) is None, library_number(x) or 10**9, x))
        return []
    libs: List[str] = []
    for part in spec.split(","):
        part = part.strip()
        if not part:
            continue
        if "-" in part and re.match(r"(?i)^(lib)?\d+-(lib)?\d+$", part):
            a, b = part.split("-", 1)
            a_n = library_number(a) if not a.isdigit() else int(a)
            b_n = library_number(b) if not b.isdigit() else int(b)
            if a_n is None or b_n is None:
                libs.append(normalize_library(part))
            else:
                for n in range(min(a_n, b_n), max(a_n, b_n) + 1):
                    libs.append(f"lib{n}")
        else:
            libs.append(normalize_library(part))
    return libs


def discover_demux_prefixes(demux_root: os.PathLike | str) -> Dict[str, str]:
    """Discover filtered CellBouncer demux prefixes under *demux_root*.

    The demux directories used by this project can contain multiple count
    prefixes for the same library, for example::

        lib8_demuxed.counts
        lib8_raw.counts
        lib8_species_counts.counts

    The swap-audit POSTHOC stage requires the filtered demux prefix with the
    full sidecar set (.assignments, .diagnostics.gz, .runner_ups.gz).  The old
    first-match discovery could select libN_raw or libN_species_counts
    depending on filesystem traversal order, causing false preflight failures
    even though libN_demuxed.* existed.

    This function therefore prefers exact ``libN_demuxed.counts`` prefixes and
    rejects raw/species/audit/pass1 count files as core audit prefixes.
    """
    root = Path(demux_root)
    out: Dict[str, str] = {}

    if not root.exists():
        return out

    def lib_from_path(p: Path) -> str:
        m = re.search(r"(lib\d+)", p.name, re.I)
        if not m:
            m = re.search(r"(lib\d+)", str(p), re.I)
        return normalize_library(m.group(1)) if m else ""

    def is_rejected_counts_name(name: str) -> bool:
        lower = name.lower()
        rejected_markers = (
            "_raw.counts",
            "_species_counts.counts",
            ".species_counts",
            ".atac.counts",
            ".pass1.counts",
            "_audit_unconstrained.counts",
            "_audit_constrained.counts",
            "_audit.counts",
        )
        return any(marker in lower for marker in rejected_markers)

    # Group all count files by library so that each library can be resolved
    # deterministically.
    by_lib: Dict[str, List[Path]] = {}
    for p in root.rglob("*.counts"):
        lib = lib_from_path(p)
        if not lib:
            continue
        by_lib.setdefault(lib, []).append(p)

    for lib, candidates in by_lib.items():
        # 1) Exact preferred filtered demux prefix in the supplied demux root.
        preferred = root / f"{lib}_demuxed.counts"
        if preferred.exists():
            out[lib] = str(preferred)[:-len(".counts")]
            continue

        # 2) Recursive preferred filtered demux prefix.
        demuxed = sorted(p for p in candidates if p.name == f"{lib}_demuxed.counts")
        if demuxed:
            out[lib] = str(demuxed[0])[:-len(".counts")]
            continue

        # 3) Last-resort legacy support: accept only non-raw, non-species,
        # non-audit prefixes that still look like filtered demux output.
        usable = sorted(
            p for p in candidates
            if not is_rejected_counts_name(p.name)
            and p.name.endswith("_demuxed.counts")
        )
        if usable:
            out[lib] = str(usable[0])[:-len(".counts")]

    return out


def file_state(path: os.PathLike | str) -> str:
    return "present" if Path(path).exists() else "missing"


def safe_symlink(src: os.PathLike | str, dst: os.PathLike | str, overwrite: bool = False) -> None:
    srcp = Path(src).resolve()
    dstp = Path(dst)
    dstp.parent.mkdir(parents=True, exist_ok=True)
    if dstp.exists() or dstp.is_symlink():
        if overwrite:
            dstp.unlink()
        else:
            # If already points to the right target, leave it alone.
            if dstp.is_symlink() and Path(os.readlink(dstp)).resolve() == srcp:
                return
            return
    os.symlink(str(srcp), str(dstp))


def mean(vals: Sequence[float]) -> float:
    vals = [v for v in vals if v is not None and not math.isnan(v)]
    return sum(vals) / len(vals) if vals else math.nan


def median(vals: Sequence[float]) -> float:
    vals = sorted(v for v in vals if v is not None and not math.isnan(v))
    if not vals:
        return math.nan
    n = len(vals)
    mid = n // 2
    if n % 2:
        return vals[mid]
    return (vals[mid - 1] + vals[mid]) / 2.0


def percentile(vals: Sequence[float], pct: float) -> float:
    vals = sorted(v for v in vals if v is not None and not math.isnan(v))
    if not vals:
        return math.nan
    if len(vals) == 1:
        return vals[0]
    pos = (len(vals) - 1) * pct / 100.0
    lo = int(math.floor(pos))
    hi = int(math.ceil(pos))
    if lo == hi:
        return vals[lo]
    return vals[lo] * (hi - pos) + vals[hi] * (pos - lo)


def truthy(x: object) -> bool:
    return str(x).strip().lower() in {"1", "true", "t", "yes", "y"}


def read_panel_metadata(path: os.PathLike | str) -> Dict[str, str]:
    rows = read_tsv(path)
    if not rows:
        return {}
    keys = rows[0].keys()
    id_col = "indiv_id" if "indiv_id" in keys else ("VCF_ID" if "VCF_ID" in keys else list(keys)[0])
    species_col = "species" if "species" in keys else ("Species" if "Species" in keys else None)
    if species_col is None:
        raise ValueError(f"panel metadata lacks species column: {path}")
    return {r[id_col].strip(): r[species_col].strip() for r in rows if r.get(id_col, "").strip()}


HYBRID_SPECIES_EXPANSION: Dict[str, Tuple[str, ...]] = {
    # Chinobo-mCherry is a bonobo-chimp F1 hybrid.  Keep it as an individual
    # identity, but expand it to native species components for species-level QC.
    "Hy": ("B", "C"),
    "Chinobo-mCherry": ("B", "C"),
}


def expand_species_label(label: str) -> Tuple[str, ...]:
    label = str(label).strip()
    if not label or label == NA:
        return tuple()
    return HYBRID_SPECIES_EXPANSION.get(label, (label,))


def species_set_members(identity: str, panel: Mapping[str, str]) -> Set[str]:
    vals: Set[str] = set()
    for c in identity_components(identity):
        raw = panel.get(c, c if c in HYBRID_SPECIES_EXPANSION else "UNKNOWN")
        for sp in expand_species_label(raw):
            if sp and sp != NA:
                vals.add(sp)
    return vals


def format_species_members(vals: Iterable[str]) -> str:
    clean = sorted({str(v).strip() for v in vals if str(v).strip() and str(v).strip() != NA})
    return ",".join(clean) if clean else NA


def species_set(identity: str, panel: Mapping[str, str]) -> str:
    return format_species_members(species_set_members(identity, panel))


def parse_species_set_value(value: object) -> Set[str]:
    raw = str(value or "").strip()
    if not raw or raw == NA:
        return set()
    out: Set[str] = set()
    for tok in re.split(r"[+,;]", raw):
        tok = tok.strip()
        if not tok or tok == NA:
            continue
        out.update(expand_species_label(tok))
    return {x for x in out if x and x != NA}


def species_relation(expected: object, observed: object) -> str:
    exp = parse_species_set_value(expected)
    obs = parse_species_set_value(observed)
    if not exp or not obs:
        return "missing_species_evidence"
    if obs == exp:
        return "exact_match"
    if obs < exp:
        return "expected_subset_only_component_missing"
    if exp < obs:
        return "expected_superset_with_extra_species"
    if obs & exp:
        return "partial_overlap_with_extra_and_missing"
    return "disjoint_wrong_species"


def species_relation_booleans(expected: object, observed: object) -> Tuple[int, int, int]:
    rel = species_relation(expected, observed)
    missing_expected = 1 if rel in {"expected_subset_only_component_missing", "partial_overlap_with_extra_and_missing", "missing_species_evidence"} else 0
    has_unexpected = 1 if rel in {"expected_superset_with_extra_species", "partial_overlap_with_extra_and_missing", "disjoint_wrong_species"} else 0
    disjoint_wrong = 1 if rel == "disjoint_wrong_species" else 0
    return missing_expected, has_unexpected, disjoint_wrong


def load_simple_table(path: os.PathLike | str) -> List[List[str]]:
    rows = []
    with open_text(path, "rt") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line:
                rows.append(line.split("\t"))
    return rows
