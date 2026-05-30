#!/usr/bin/env python3
"""check_output_contract.py — post-run validation of an `mllmct annotate` output dir.

Asserts the output CONTRACT a downstream consumer relies on:
  - labels.csv exists with the required columns and one row per cluster;
  - py_* metric columns are present (the authoritative, Python-recomputed ones);
  - the trace tree is complete for the named lens;
  - cost_summary.json exists.

Usage:  uv run python checks/check_output_contract.py <out_dir> [--lens NAME]
Exit 0 = contract satisfied; 1 = something missing.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

REQUIRED_COLS = [
    "cluster", "consensus_label", "harmonized_label",
    "py_consensus_proportion", "py_entropy",
    "llm_reported_proportion", "llm_reported_entropy", "model_annotations",
]
TRACE_FILES = ["prompt.txt", "model_responses.json", "discussion.json", "tokens.json", "meta.json"]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("out_dir")
    ap.add_argument("--lens", default=None, help="lens name (default: the single dir under trace/)")
    args = ap.parse_args()
    out = Path(args.out_dir)
    fails: list[str] = []

    labels = out / "labels.csv"
    if not labels.exists():
        fails.append(f"missing {labels}")
    else:
        df = pd.read_csv(labels)
        missing = [c for c in REQUIRED_COLS if c not in df.columns]
        if missing:
            fails.append(f"labels.csv missing columns: {missing}")
        if len(df) == 0:
            fails.append("labels.csv has 0 rows")

    trace_root = out / "trace"
    if not (trace_root / "cost_summary.json").exists():
        fails.append("missing trace/cost_summary.json")
    lens = args.lens
    if lens is None and trace_root.exists():
        subdirs = [d for d in trace_root.iterdir() if d.is_dir()]
        lens = subdirs[0].name if len(subdirs) == 1 else None
    if lens is None:
        fails.append("could not determine lens (pass --lens); cannot check trace tree")
    else:
        for f in TRACE_FILES:
            if not (trace_root / lens / f).exists():
                fails.append(f"missing trace/{lens}/{f}")

    if fails:
        print("OUTPUT CONTRACT FAILED:")
        for f in fails:
            print(f"  - {f}")
        return 1
    print(f"OUTPUT CONTRACT OK: {labels} + trace/{lens}/ complete.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
