"""check_program_redundancy.py — flag cross-run pairs with r > 0.9 outside merged clusters.

A pair of programs from different variants with very high correlation (r > 0.9)
that did NOT end up in the same merged cluster is suspicious — usually a sign
that K was too high in one of the variants (the program was duplicated within
the variant) or that the merge threshold was too strict.

Usage:
    python check_program_redundancy.py <correlation_csv> <merged_csv>
        [--threshold 0.9]

Exits 0 on no flags, 1 on flagged pairs (so a CI step can fail fast).
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def check(
    correlation_csv: Path,
    merged_csv: Path,
    threshold: float = 0.9,
) -> int:
    corr = pd.read_csv(correlation_csv, index_col=0)
    merged = pd.read_csv(merged_csv)

    # Build program -> cluster_id map
    prog_to_cluster: dict[str, int] = {}
    for _, row in merged.iterrows():
        cluster_id = int(row["cluster_id"])
        for prog in str(row["source_programs"]).split(","):
            prog_to_cluster[prog.strip()] = cluster_id

    flags: list[tuple[str, str, float]] = []
    cols = list(corr.columns)
    for i, c1 in enumerate(cols):
        for c2 in cols[i + 1:]:
            r = corr.loc[c1, c2]
            if pd.isna(r):
                continue
            if abs(r) <= threshold:
                continue
            cluster1 = prog_to_cluster.get(c1)
            cluster2 = prog_to_cluster.get(c2)
            if cluster1 is None or cluster2 is None:
                continue
            if cluster1 != cluster2:
                flags.append((c1, c2, float(r)))

    if not flags:
        print(f"OK: no cross-cluster pairs with |r| > {threshold}.")
        return 0

    print(f"FAIL: {len(flags)} cross-cluster pairs with |r| > {threshold}:")
    for c1, c2, r in flags:
        print(f"  {c1} <-> {c2}: r={r:.3f} "
              f"(clusters {prog_to_cluster[c1]} vs {prog_to_cluster[c2]})")
    print("\n  Likely causes:")
    print("    1. K too high in one variant — duplicate programs within a variant")
    print("    2. Merge threshold too strict — programs that should have merged were split")
    print("    3. Correlation matrix has NaN-filled regions causing distance distortion")
    return 1


def _parse() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("correlation_csv", type=Path)
    p.add_argument("merged_csv", type=Path)
    p.add_argument("--threshold", type=float, default=0.9)
    return p.parse_args()


def main() -> int:
    a = _parse()
    return check(a.correlation_csv, a.merged_csv, a.threshold)


if __name__ == "__main__":
    raise SystemExit(main())
