#!/usr/bin/env python3
"""Shared helpers for the CellBouncer identity-reconciliation pipeline.

This module deliberately contains no ambient-RNA/contamination integration.
"""
from __future__ import annotations

import csv
import gzip
import hashlib
import json
import math
import os
import re
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Mapping, MutableMapping, Optional, Sequence, Set, Tuple

SCHEMA_VERSION = "identity_reconciliation_v1"
POLICY_VERSION = "identity_reconciliation_policy_v11_component_singlets_and_library_exchange"
NA_TOKENS = {"", ".", "na", "nan", "none", "null", "unavailable"}

DEFAULT_IDENTITY_POLICY = {
    "policy_version": POLICY_VERSION,
    "schema_version": SCHEMA_VERSION,
    # Conservative pilot thresholds. These are intentionally not treated as a
    # calibrated production policy; auto-apply is disabled by default below.
    "nuclear_strong_delta_ll": 100.0,
    "nuclear_suggestive_delta_ll": 20.0,
    "nuclear_strong_depth_normalized_delta": 0.01,
    "nuclear_min_informative_depth": 5000,
    "nuclear_low_information_depth": 500,
    "site_fold_min_evaluable": 4,
    "site_fold_support_fraction": 0.80,
    "mt_support_delta_ll": 4.0,
    "atac_support_delta_ll": 100.0,
    "nn_homotypic_prob_tet": 0.90,
    "event_min_cells": 10,
    # Unexpected non-homotet identities may be reported as coherent events at
    # event_min_cells, but production assignments are only rewritten once the
    # exact canonical identity has this many robust cells in the same library.
    # NN-driven A -> A+A homotet reclassification is exempt from this mass
    # threshold, but the occupancy guards below still apply.
    "unexpected_line_autoapply_min_cells": 100,
    # A repeated A+B donor-composition signal does not by itself distinguish an
    # intact biological T[A+B] cell from a technical M{D[A]|D[B]} droplet when
    # both donor singlets are independently present in the same library.
    "heterotypic_pair_occupancy_guard": True,
    # Likewise, a library-unexpected A+A call is not auto-applied from NN ploidy
    # evidence alone when diploid A is independently present; A+A and a same-
    # donor technical doublet are then observationally confounded here.
    "unexpected_homotet_occupancy_guard": True,
    # A positive quad-pattern score is genotype-composition evidence, not
    # independent proof that multiple cells occupied the droplet.  Only an
    # explicit non-single droplet flag is allowed to trigger the hard
    # technical-multiplet veto.
    "explicit_multiplet_evidence_source": "droplet_flag_only",
    "event_partner_preservation_fraction": 0.80,
    "event_primary_source_fraction": 0.80,
    "full_replacement_fraction": 0.90,
    "partial_replacement_fraction": 0.20,
    "foreign_carryover_max_fraction": 0.10,
    # Explicit opt-in until held-out calibration establishes the false-change rate.
    "auto_apply_decisive": False,
}
DEFAULT_IDENTITY_POLICY_SHA256 = hashlib.sha256(
    json.dumps(DEFAULT_IDENTITY_POLICY, sort_keys=True, separators=(",", ":")).encode("utf-8")
).hexdigest()


def clean(value) -> str:
    if value is None:
        return ""
    s = str(value).strip()
    return "" if s.lower() in NA_TOKENS else s


def natural_key(value: str):
    s = clean(value)
    if re.fullmatch(r"\d+", s):
        return (0, int(s), "")
    return (1, 0, s)


def canonical_uid_set(values: Iterable[str]) -> str:
    toks = set()
    for value in values:
        for token in re.split(r"[|,;]", clean(value)):
            token = clean(token)
            if token:
                toks.add(token)
    return "|".join(sorted(toks, key=natural_key))


def uid_count(value: str) -> int:
    return len([x for x in clean(value).split("|") if x])


def canonical_genotype(value: str, aliases: Optional[Mapping[str, str]] = None) -> str:
    aliases = aliases or {}
    raw = clean(value)
    if not raw:
        return ""
    # Remove state wrappers used by the reconciliation layer.
    m = re.fullmatch(r"(?:D|T|UNKNOWN_SINGLE_CELL)\[(.*)\]", raw)
    if m:
        raw = m.group(1)
    if raw.startswith("M{"):
        return raw
    parts = [clean(x) for x in raw.replace("x", "+").split("+")]
    parts = [aliases.get(x, x) for x in parts if x]
    return "+".join(sorted(parts))


def donor_components(genotype: str) -> List[str]:
    g = canonical_genotype(genotype)
    if not g or g.startswith("M{"):
        return []
    return [x for x in g.split("+") if x]


def expected_library_context(expected_genotypes_by_library: Mapping[str, object]):
    """Derive canonical per-library singlet/component/composite metadata views.

    ``expected_genotypes_by_library`` may map each library either to a mapping
    keyed by genotype or to an iterable of genotype strings.  All parsing flows
    through ``canonical_genotype``/``donor_components`` so this does not create a
    second genotype interpretation path.
    """
    explicit_singlets: Dict[str, Set[str]] = {}
    donor_universe: Dict[str, Set[str]] = {}
    composites_by_component: Dict[str, Dict[str, List[str]]] = {}

    for lib, value in expected_genotypes_by_library.items():
        genotypes = value.keys() if isinstance(value, Mapping) else value
        singlets: Set[str] = set()
        donors: Set[str] = set()
        composites = {}
        for raw in genotypes:
            genotype = canonical_genotype(raw)
            comps = donor_components(genotype)
            if not comps:
                continue
            donors.update(comps)
            if len(comps) == 1:
                singlets.add(comps[0])
            else:
                for donor in set(comps):
                    composites.setdefault(donor, set()).add(genotype)
        explicit_singlets[lib] = singlets
        donor_universe[lib] = donors
        composites_by_component[lib] = {
            donor: sorted(values, key=natural_key)
            for donor, values in composites.items()
        }
    return explicit_singlets, donor_universe, composites_by_component


def classify_singlet_library_relationship(
    genotype: str,
    explicit_expected_singlets: Iterable[str],
    expected_donor_components: Iterable[str],
    expected_composites_by_component: Mapping[str, Sequence[str]],
    globally_valid: bool = True,
) -> Tuple[str, str]:
    """Classify one single-donor identity relative to a library's biology."""
    comps = donor_components(genotype)
    if len(comps) != 1:
        return "NOT_SINGLET", ""
    donor = comps[0]
    if not globally_valid:
        return "UNRESOLVED", ""
    singlets = set(explicit_expected_singlets)
    components = set(expected_donor_components)
    if donor in singlets:
        return "EXPECTED_SINGLET", ""
    if donor in components:
        context = sorted(
            {canonical_genotype(x) for x in expected_composites_by_component.get(donor, []) if canonical_genotype(x)},
            key=natural_key,
        )
        return "EXPECTED_COMPOSITE_COMPONENT_SINGLET", ";".join(context)
    return "UNEXPECTED_SINGLET", ""


def pairwise_library_roster(
    library_a: str, library_b: str, expected_donor_components_by_library: Mapping[str, Iterable[str]]
) -> Dict[str, object]:
    """Return transparent donor-roster overlap and discriminability for a pair."""
    a = set(expected_donor_components_by_library.get(library_a, set()))
    b = set(expected_donor_components_by_library.get(library_b, set()))
    shared = a & b
    a_specific = a - b
    b_specific = b - a
    n_distinguishing = len(a_specific) + len(b_specific)
    if not a_specific and not b_specific:
        relation = "ROSTER_EQUIVALENT_NONDISCRIMINATING"
        discriminability = "NONE"
    else:
        if not shared:
            relation = "DISTINCT_DONOR_ROSTERS"
        elif not a_specific or not b_specific:
            relation = "NESTED_DONOR_ROSTERS"
        else:
            relation = "PARTIALLY_OVERLAPPING_DONOR_ROSTERS"
        if n_distinguishing <= 2:
            discriminability = "WEAK"
        elif n_distinguishing <= 5:
            discriminability = "MODERATE"
        else:
            discriminability = "STRONG"
    return {
        "roster_relation": relation,
        "pair_discriminability": discriminability,
        "shared": shared,
        "a_specific": a_specific,
        "b_specific": b_specific,
    }


def biological_state(genotype: str, ploidy: str = "") -> str:
    comps = donor_components(genotype)
    if not comps:
        return "UNKNOWN"
    pl = clean(ploidy).upper()
    if len(comps) == 1:
        if pl in {"T", "TET", "TETRAPLOID", "HOMOTYPIC"}:
            return f"T[{comps[0]}+{comps[0]}]"
        return f"D[{comps[0]}]"
    if len(comps) == 2:
        return f"T[{'+'.join(comps)}]"
    return "M{" + "+".join(comps) + "}"


def parse_library(value) -> str:
    s = clean(value)
    m = re.search(r"(?:lib|RNA[_-]?)(\d+)$", s, flags=re.I)
    if m:
        return f"lib{int(m.group(1))}"
    if re.fullmatch(r"\d+", s):
        return f"lib{int(s)}"
    return s


def library_number(value) -> Optional[int]:
    lib = parse_library(value)
    m = re.fullmatch(r"lib(\d+)", lib)
    return int(m.group(1)) if m else None


def parse_library_spec(values: Sequence[str]) -> List[int]:
    out = set()
    for item in values:
        for token in str(item).replace(",", " ").split():
            if "-" in token:
                a, b = token.split("-", 1)
                out.update(range(int(a), int(b) + 1))
            else:
                n = library_number(token)
                out.add(n if n is not None else int(token))
    return sorted(out)


def open_text(path: str, mode: str = "rt"):
    return gzip.open(path, mode, newline="") if str(path).endswith(".gz") else open(path, mode, newline="")


def read_tsv(path: str) -> List[Dict[str, str]]:
    with open_text(path, "rt") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        return [dict(row) for row in reader]


def iter_tsv(path: str) -> Iterator[Dict[str, str]]:
    with open_text(path, "rt") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            yield dict(row)


def write_tsv(path: str, rows: Iterable[Mapping[str, object]], fields: Sequence[str]) -> None:
    path = str(path)
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    tmp = path + ".tmp"
    opener = gzip.open if path.endswith(".gz") else open
    kwargs = {"newline": ""}
    if path.endswith(".gz"):
        kwargs["mode"] = "wt"
    else:
        kwargs["mode"] = "w"
    with opener(tmp, **kwargs) as fh:
        writer = csv.DictWriter(fh, fieldnames=list(fields), delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: format_value(row.get(field, "")) for field in fields})
    os.replace(tmp, path)


def write_headerless_tsv(path: str, rows: Iterable[Sequence[object]]) -> None:
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    tmp = str(path) + ".tmp"
    with open(tmp, "w", newline="") as fh:
        for row in rows:
            fh.write("\t".join(format_value(x) for x in row) + "\n")
    os.replace(tmp, path)


def format_value(value) -> str:
    if value is None:
        return "NA"
    if isinstance(value, bool):
        return "TRUE" if value else "FALSE"
    if isinstance(value, float):
        if not math.isfinite(value):
            return "NA"
        return f"{value:.10g}"
    s = str(value)
    return s if s else "NA"


def ffloat(value, default: float = math.nan) -> float:
    try:
        x = float(clean(value))
        return x if math.isfinite(x) else default
    except Exception:
        return default


def fint(value, default: int = 0) -> int:
    try:
        return int(float(clean(value)))
    except Exception:
        return default


def sha256_file(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def stable_site_fold(tid: int, pos: int, folds: int) -> int:
    if folds <= 1:
        return 0
    payload = f"{tid}:{pos}".encode()
    digest = hashlib.blake2b(payload, digest_size=8, person=b"cellbounc").digest()
    return int.from_bytes(digest, "little") % folds


def read_assignments(path: str) -> Dict[str, Dict[str, str]]:
    """Read the production headerless four-column .assignments convention."""
    out: Dict[str, Dict[str, str]] = {}
    with open_text(path, "rt") as fh:
        for line_no, raw in enumerate(fh, 1):
            if not raw.strip() or raw.startswith("#"):
                continue
            parts = raw.rstrip("\n\r").split("\t")
            if len(parts) < 2:
                continue
            bc = clean(parts[0])
            out[bc] = {
                "barcode": bc,
                "assignment": canonical_genotype(parts[1]),
                "type": clean(parts[2]) if len(parts) > 2 else "",
                "score": clean(parts[3]) if len(parts) > 3 else "",
                "source_line": str(line_no),
            }
    return out


def read_refined_assignments(path: str) -> Dict[str, Dict[str, str]]:
    if not path or not os.path.isfile(path):
        return {}
    rows = read_tsv(path)
    out = {}
    for row in rows:
        bc = clean(row.get("barcode"))
        if bc:
            out[bc] = row
    return out


def load_alias_table(path: str) -> Dict[str, str]:
    if not path or not os.path.isfile(path):
        return {}
    aliases: Dict[str, str] = {}
    for row in read_tsv(path):
        src = clean(row.get("library_label") or row.get("alias") or row.get("source"))
        dst = clean(row.get("canonical_vcf_id") or row.get("canonical") or row.get("target"))
        if src and dst:
            aliases[src] = dst
        if dst:
            aliases.setdefault(dst, dst)
        for eq in re.split(r"[,;|]", clean(row.get("equivalent_or_reporter_vcf_ids"))):
            eq = clean(eq)
            if eq and dst:
                aliases.setdefault(eq, dst)
    return aliases


def choose_first(row: Mapping[str, object], names: Sequence[str], default: str = "") -> str:
    normalized = {re.sub(r"[^a-z0-9]", "", str(k).lower()): k for k in row}
    for name in names:
        key = normalized.get(re.sub(r"[^a-z0-9]", "", name.lower()))
        if key is not None:
            val = clean(row.get(key))
            if val:
                return val
    return default


def json_dump_atomic(path: str, obj) -> None:
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    tmp = path + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(obj, fh, indent=2, sort_keys=True)
        fh.write("\n")
    os.replace(tmp, path)
