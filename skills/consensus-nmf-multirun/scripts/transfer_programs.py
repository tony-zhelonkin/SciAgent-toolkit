"""transfer_programs.py — Stage 3, port of <ref-scrna>'s 06e_transfer_programs.py.

Score per-variant cNMF programs onto the full dataset's barcode space and
compute the cross-source correlation matrix that Stage 5 consumes.
"""

from __future__ import annotations

import argparse
import pickle
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc


def transfer(
    base_h5ad: str | Path,
    sources: dict[str, str | Path],
    out_h5ad: str | Path,
    correlation_csv: str | Path,
    redundancy_threshold: float = 0.7,
) -> tuple[ad.AnnData, pd.DataFrame]:
    """Join per-variant usage matrices into base AnnData; write correlation matrix.

    sources: dict mapping variant_name -> path to <variant>_results.pkl from Stage 1.
    """
    base_h5ad = Path(base_h5ad)
    out_h5ad = Path(out_h5ad)
    correlation_csv = Path(correlation_csv)
    out_h5ad.parent.mkdir(parents=True, exist_ok=True)
    correlation_csv.parent.mkdir(parents=True, exist_ok=True)

    adata = sc.read_h5ad(base_h5ad)
    print(f"Loaded base {base_h5ad}: {adata.n_obs:,} cells × {adata.n_vars:,} genes")

    summary: list[dict] = []

    for variant_name, pkl_path in sources.items():
        pkl_path = Path(pkl_path)
        if not pkl_path.exists():
            print(f"  skipping {variant_name}: {pkl_path} not found")
            continue

        with open(pkl_path, "rb") as f:
            res = pickle.load(f)
        usage: pd.DataFrame = res["usage"].copy()
        k = int(res.get("k", usage.shape[1]))

        # Sanity-check usage values
        if not ((usage.values >= 0).all() and (usage.values.sum(axis=1) <= 1.001).all()):
            print(f"  WARNING: {variant_name} usage values are outside [0, 1] simplex; "
                  "check whether spectra_scores were used by accident.")

        # Rename columns -> cNMF_<variant>_P<i>
        usage.columns = [f"cNMF_{variant_name}_P{i+1}" for i in range(usage.shape[1])]

        # Initialise NaN columns
        for col in usage.columns:
            adata.obs[col] = np.nan

        # Join on intersection of cells
        common = adata.obs_names.intersection(usage.index)
        coverage = len(common) / adata.n_obs
        adata.obs.loc[common, list(usage.columns)] = usage.loc[common].values

        print(f"  {variant_name}: K={k}, {len(common):,}/{adata.n_obs:,} cells "
              f"({coverage * 100:.1f}%) scored")
        summary.append({
            "variant": variant_name,
            "k": k,
            "n_scored": len(common),
            "coverage": round(coverage, 4),
            "programs": ",".join(usage.columns),
        })

    program_cols = [c for c in adata.obs.columns if c.startswith("cNMF_")]
    if not program_cols:
        raise RuntimeError("No cNMF columns produced; no variants succeeded.")

    # Cross-source correlation
    prog_data = adata.obs[program_cols].dropna(how="all")
    corr_matrix = prog_data.corr()
    corr_matrix.to_csv(correlation_csv)
    print(f"\n  wrote correlation matrix: {correlation_csv}")

    # Surface high-correlation pairs (potential merge candidates)
    print(f"\n  Cross-program pairs with |r| > {redundancy_threshold}:")
    n_flagged = 0
    for i, c1 in enumerate(corr_matrix.columns):
        for c2 in corr_matrix.columns[i + 1:]:
            r = corr_matrix.loc[c1, c2]
            if pd.notna(r) and abs(r) > redundancy_threshold:
                print(f"    {c1} <-> {c2}: r={r:.3f}")
                n_flagged += 1
    print(f"  {n_flagged} pairs flagged.")

    # Write checkpoint + summary
    adata.write_h5ad(out_h5ad)
    print(f"  wrote {out_h5ad}")
    summary_df = pd.DataFrame(summary)
    summary_csv = correlation_csv.parent / "program_sources_summary.csv"
    summary_df.to_csv(summary_csv, index=False)
    print(f"  wrote {summary_csv}")

    return adata, corr_matrix


def _parse() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Transfer per-variant cNMF programs onto base AnnData.")
    p.add_argument("--base", dest="base_h5ad", type=Path, required=True)
    p.add_argument("--source", action="append", required=True,
                   metavar="NAME=PATH",
                   help="Variant name and path to <variant>_results.pkl. "
                        "Pass multiple times: --source full=03_results/cnmf/full/full_results.pkl")
    p.add_argument("--out", dest="out_h5ad", type=Path, required=True)
    p.add_argument("--correlation-csv", type=Path, required=True)
    p.add_argument("--redundancy-threshold", type=float, default=0.7)
    return p.parse_args()


def main() -> int:
    a = _parse()
    sources: dict[str, str] = {}
    for entry in a.source:
        if "=" not in entry:
            print(f"--source must be NAME=PATH (got {entry!r})", file=sys.stderr)
            return 2
        name, path = entry.split("=", 1)
        sources[name] = path
    transfer(
        base_h5ad=a.base_h5ad,
        sources=sources,
        out_h5ad=a.out_h5ad,
        correlation_csv=a.correlation_csv,
        redundancy_threshold=a.redundancy_threshold,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
