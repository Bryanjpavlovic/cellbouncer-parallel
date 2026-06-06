#!/usr/bin/env python3
"""Run tetra_refine for one Tet2025 library with optional NN ploidy calls.

This wrapper is intentionally thin: it does not reimplement tetra_refine logic.
It resolves the standard Tet2025 file layout, validates what is present, and
passes optional evidence layers only when available unless --require-* flags are
set.

Default output:
  /mnt/beegfs/tetmultiome_rna_mapped/QUANT3_CONTAM_V2/tetra_refine/lib<N>/lib<N>.refined_assignments

The output root is shaped to match swap_audit_prepare.py's existing
--refined_assignments_root search pattern:
  <root>/lib<N>/lib<N>.refined_assignments
"""
from __future__ import annotations

import argparse
import os
import shlex
import subprocess
import sys
from pathlib import Path
from typing import List, Optional

BEEGFS_ROOT = Path("/mnt/beegfs/tetmultiome_rna_mapped/mapping_output")
QUANT_ROOT = Path("/mnt/beegfs/tetmultiome_rna_mapped/QUANT3_CONTAM_V2")
EXPECTED_LINES_DIR = Path("/mnt/beegfs/tetmultiome_rna_mapped/Misc_Metadata")
SOFTWARE_BIN = Path("/nvme/software/packages/cellbouncer/dev/bin")
LIB_PREFIX = "Tet_2025_Multiome-RNA_"
DEMUX_SUBDIR = "demux_het_v2_test"


def positive_int(text: str) -> int:
    try:
        value = int(str(text).replace("lib", ""))
    except ValueError as exc:
        raise argparse.ArgumentTypeError(f"not a library number: {text!r}") from exc
    if value < 1:
        raise argparse.ArgumentTypeError("library number must be >= 1")
    return value


def default_demux_prefix(lib: int) -> Path:
    return BEEGFS_ROOT / f"{LIB_PREFIX}{lib}" / DEMUX_SUBDIR / f"lib{lib}_demuxed"


def default_expected_lines(lib: int) -> Path:
    """Resolve expected lines using the same search order as orchestrate_pipeline.py."""
    candidates = [
        BEEGFS_ROOT / f"{LIB_PREFIX}{lib}" / DEMUX_SUBDIR / f"lib{lib}_expected_lines.txt",
        BEEGFS_ROOT / f"{LIB_PREFIX}{lib}" / f"lib{lib}_expected_lines.txt",
        EXPECTED_LINES_DIR / f"lib{lib}_expected_lines.txt",
    ]
    for c in candidates:
        if c.is_file():
            return c
    # Return the first candidate for a helpful error message from require_file().
    return candidates[0]


def default_ploidy_calls(ploidy_root: Path, lib: int) -> Path:
    return ploidy_root / f"lib{lib}.ploidy_calls_nn.tsv"


def default_contam_rate(quant_root: Path, condition: str, lib: int) -> Path:
    return quant_root / condition / f"lib{lib}" / f"lib{lib}_demuxed.contam_rate"


def require_file(path: Path, label: str) -> None:
    if not path.is_file():
        raise FileNotFoundError(f"{label} missing: {path}")


def existing_optional(path: Optional[Path], label: str, required: bool) -> Optional[Path]:
    if path is None:
        return None
    if path.is_file():
        return path
    msg = f"{label} missing: {path}"
    if required:
        raise FileNotFoundError(msg)
    print(f"WARNING: {msg}; continuing without it", file=sys.stderr)
    return None


def build_command(args: argparse.Namespace) -> List[str]:
    lib = args.lib
    tetra_refine = Path(args.tetra_refine_bin)
    demux_prefix = Path(args.demux_prefix) if args.demux_prefix else default_demux_prefix(lib)
    expected = Path(args.expected) if args.expected else default_expected_lines(lib)
    out_dir = Path(args.output_root) / f"lib{lib}"
    out_prefix = out_dir / f"lib{lib}"

    require_file(tetra_refine, "tetra_refine binary")
    for suffix, label in [
        (".assignments", "demux assignments"),
        (".diagnostics.gz", "demux diagnostics"),
        (".runner_ups.gz", "demux runner-ups"),
    ]:
        require_file(Path(str(demux_prefix) + suffix), label)
    require_file(expected, "expected lines")

    ploidy_path: Optional[Path] = None
    if args.external_ploidy:
        ploidy_path = Path(args.external_ploidy)
    elif args.ploidy_calls_root:
        ploidy_path = default_ploidy_calls(Path(args.ploidy_calls_root), lib)
    ploidy_path = existing_optional(ploidy_path, "external ploidy calls", args.require_external_ploidy)

    contam_rate_path: Optional[Path] = None
    if args.contam_rate:
        contam_rate_path = Path(args.contam_rate)
    elif args.contam_condition:
        contam_rate_path = default_contam_rate(Path(args.quant_root), args.contam_condition, lib)
    contam_rate_path = existing_optional(contam_rate_path, "contam_rate", args.require_contam_rate)

    out_dir.mkdir(parents=True, exist_ok=True)

    cmd = [
        str(tetra_refine),
        "--assignments", str(demux_prefix) + ".assignments",
        "--diagnostics", str(demux_prefix) + ".diagnostics.gz",
        "--runner_ups", str(demux_prefix) + ".runner_ups.gz",
        "--expected", str(expected),
        "--output", str(out_prefix),
        "--write_changed_only",
        "--write_simple",
        "--verbose",
    ]
    if contam_rate_path:
        cmd += ["--contam_rate", str(contam_rate_path)]
    if ploidy_path:
        cmd += [
            "--external_ploidy", str(ploidy_path),
            "--external_ploidy_library", str(lib),
            "--external_ploidy_min_prob", str(args.external_ploidy_min_prob),
        ]
    if args.scoring_only:
        cmd += ["--scoring_only"]
    return cmd


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--lib", required=True, type=positive_int, help="Library number, e.g. 19 or lib19")
    p.add_argument("--tetra_refine_bin", default=str(SOFTWARE_BIN / "tetra_refine"))
    p.add_argument("--demux_prefix", default=None, help="Override full demux prefix without suffix")
    p.add_argument("--expected", default=None, help="Override expected lines file")
    p.add_argument("--output_root", default=str(QUANT_ROOT / "tetra_refine"))
    p.add_argument("--ploidy_calls_root", default=str(Path("/mnt/beegfs/tetmultiome_rna_mapped/mapping_output/ploidy_calls_unfiltered")),
                   help="Root containing lib<N>.ploidy_calls_nn.tsv")
    p.add_argument("--external_ploidy", default=None, help="Override external ploidy file path")
    p.add_argument("--external_ploidy_min_prob", type=float, default=0.90,
                   help="Minimum external ploidy confidence required to relabel A -> A+A")
    p.add_argument("--require_external_ploidy", action="store_true", help="Fail instead of warning if NN ploidy calls are missing")
    p.add_argument("--quant_root", default=str(QUANT_ROOT))
    p.add_argument("--contam_condition", default="IND_WS_RA_LOO", help="Condition used to auto-locate optional contam_rate; empty string disables")
    p.add_argument("--contam_rate", default=None, help="Override contam_rate path")
    p.add_argument("--require_contam_rate", action="store_true", help="Fail instead of warning if contam_rate is missing")
    p.add_argument("--scoring_only", action="store_true", help="Pass through assignments and only emit scoring columns")
    p.add_argument("--dry_run", action="store_true")
    return p.parse_args()


def main() -> int:
    args = parse_args()
    if args.contam_condition == "":
        args.contam_condition = None
    if not (0.0 <= args.external_ploidy_min_prob <= 1.0):
        print("ERROR: --external_ploidy_min_prob must be between 0 and 1", file=sys.stderr)
        return 2
    try:
        cmd = build_command(args)
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2

    print("Running tetra_refine command:")
    print(" ".join(shlex.quote(x) for x in cmd))
    if args.dry_run:
        return 0
    proc = subprocess.run(cmd)
    return proc.returncode


if __name__ == "__main__":
    raise SystemExit(main())
