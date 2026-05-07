"""run_cnmf_subset.py — Stage 1, parameterised port of 13403-YD's 06_run_cnmf.py family.

One call per (subset × QC) variant. Decision Pauses 1, 2 resolve the variant
list before this fires; Decision Pause 3 (K) fires after the call returns.

Default parameters lifted from the reference:
  - Full dataset:  K = range(8, 13)
  - Subset:        K = range(6, 12)
  - n_iter = 100, num_highvar_genes = 2000, density_threshold = 0.01, seed = 42
  - QC filter:     n_genes 200-6000, pct_mt ≤ 20, pct_ribo ≤ 50, no doublets
"""

from __future__ import annotations

import argparse
import pickle
import sys
from pathlib import Path
from typing import Iterable

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc


QC_DEFAULTS = {
    "min_genes": 200,
    "max_genes": 6000,
    "max_pct_mt": 20.0,
    "max_pct_ribo": 50.0,
    "remove_doublets": True,
}


def apply_qc_filter(adata: ad.AnnData, qc: dict | None = None) -> ad.AnnData:
    """Apply the reference QC threshold panel; return filtered copy."""
    cfg = {**QC_DEFAULTS, **(qc or {})}
    n0 = adata.n_obs

    masks: list[pd.Series] = []
    if "n_genes_by_counts" in adata.obs:
        masks.append(
            (adata.obs["n_genes_by_counts"] >= cfg["min_genes"])
            & (adata.obs["n_genes_by_counts"] <= cfg["max_genes"])
        )
    if "pct_counts_mt" in adata.obs:
        masks.append(adata.obs["pct_counts_mt"] <= cfg["max_pct_mt"])
    if "pct_counts_ribo" in adata.obs:
        masks.append(adata.obs["pct_counts_ribo"] <= cfg["max_pct_ribo"])
    if cfg["remove_doublets"] and "predicted_doublet" in adata.obs:
        masks.append(~adata.obs["predicted_doublet"].fillna(False))

    if not masks:
        print("  WARNING: no QC columns found in obs; QC filter is a no-op")
        return adata.copy()

    combined = masks[0]
    for m in masks[1:]:
        combined = combined & m
    out = adata[combined].copy()
    print(f"  QC filter: {n0:,} -> {out.n_obs:,} cells "
          f"({out.n_obs / n0 * 100:.1f}% retained)")
    return out


def apply_obs_filter(adata: ad.AnnData, expr: str | None) -> ad.AnnData:
    """Apply user-defined obs filter expression; return filtered copy."""
    if not expr:
        return adata
    mask = adata.obs.eval(expr)
    n_kept = int(mask.sum())
    if n_kept == 0:
        raise ValueError(f"obs filter {expr!r} matched 0 cells")
    if n_kept < 1000:
        print(f"  WARNING: subset has {n_kept:,} cells; cNMF stability may be degraded")
    return adata[mask].copy()


def run(
    in_h5ad: str | Path,
    output_dir: str | Path,
    name: str,
    obs_filter: str | None = None,
    qc_filter: bool = False,
    qc_overrides: dict | None = None,
    k_range: Iterable[int] = range(8, 13),
    n_iter: int = 100,
    n_top_genes: int = 2000,
    density_threshold: float = 0.01,
    seed: int = 42,
    counts_layer: str = "counts",
    consensus_k: int | None = None,
) -> Path:
    """Run cNMF prepare → factorize → combine → k_selection_plot.

    If `consensus_k` is set, also runs `consensus()` and dumps results.pkl.
    Returns the variant directory (`<output_dir>/<name>/`).
    """
    try:
        from cnmf import cNMF
    except ImportError as e:
        raise ImportError(
            "cnmf is not installed in this environment. "
            "Install via `pip install cnmf` or check the project's env."
        ) from e

    in_h5ad = Path(in_h5ad)
    output_dir = Path(output_dir)
    variant_dir = output_dir / name
    variant_dir.parent.mkdir(parents=True, exist_ok=True)

    adata = sc.read_h5ad(in_h5ad)
    print(f"Loaded {in_h5ad}: {adata.n_obs:,} cells × {adata.n_vars:,} genes")

    adata = apply_obs_filter(adata, obs_filter)
    if qc_filter:
        adata = apply_qc_filter(adata, qc_overrides)

    if counts_layer not in adata.layers:
        raise ValueError(
            f"layers[{counts_layer!r}] not present. cNMF needs raw counts; "
            "see SKILL.md Common Pitfalls 'Forgetting raw counts in layers'."
        )

    # cNMF wants the counts in .X (it ignores layers)
    counts_adata = adata.copy()
    counts_adata.X = adata.layers[counts_layer]
    counts_path = output_dir / f"{name}_counts.h5ad"
    counts_adata.write_h5ad(counts_path)
    print(f"  wrote counts file: {counts_path}")

    cnmf_obj = cNMF(output_dir=str(output_dir), name=name)
    print(f"cNMF prepare: K={list(k_range)}, n_iter={n_iter}, n_hvg={n_top_genes}")
    cnmf_obj.prepare(
        counts_fn=str(counts_path),
        components=np.array(list(k_range)),
        n_iter=n_iter,
        num_highvar_genes=n_top_genes,
        seed=seed,
    )
    print("cNMF factorize ...")
    cnmf_obj.factorize(worker_i=0, total_workers=1)
    print("cNMF combine ...")
    cnmf_obj.combine()
    print("cNMF k_selection_plot ...")
    cnmf_obj.k_selection_plot()
    print(f"  wrote {variant_dir / f'{name}.k_selection.png'}")

    if consensus_k is not None:
        print(f"cNMF consensus K={consensus_k}, density_threshold={density_threshold}")
        cnmf_obj.consensus(k=consensus_k, density_threshold=density_threshold,
                           show_clustering=True)
        usage, spectra_scores, spectra_tpm, top_genes = cnmf_obj.load_results(
            K=consensus_k, density_threshold=density_threshold,
        )
        results = {
            "usage": usage,
            "spectra_scores": spectra_scores,
            "spectra_tpm": spectra_tpm,
            "top_genes": top_genes,
            "k": consensus_k,
        }
        results_path = output_dir / f"{name}_results.pkl"
        with open(results_path, "wb") as f:
            pickle.dump(results, f)
        print(f"  wrote {results_path} (K={consensus_k}, "
              f"{usage.shape[0]:,} cells × {usage.shape[1]} programs)")

    return variant_dir


def _parse(argv: Iterable[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Run one cNMF variant.")
    p.add_argument("--in", dest="in_h5ad", type=Path, required=True)
    p.add_argument("--output-dir", type=Path, required=True)
    p.add_argument("--name", type=str, required=True,
                   help="Variant name (e.g., 'full', 'fullQC', 'T_cells', 'T_cellsQC').")
    p.add_argument("--obs-filter", type=str, default=None,
                   help="pandas eval expression to subset obs (e.g., \"celltype == 'T_cells'\")")
    p.add_argument("--qc-filter", action="store_true",
                   help="Apply project QC thresholds (see SKILL.md QC defaults).")
    p.add_argument("--k-min", type=int, default=8)
    p.add_argument("--k-max", type=int, default=12)
    p.add_argument("--n-iter", type=int, default=100)
    p.add_argument("--n-top-genes", type=int, default=2000)
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--consensus-k", type=int, default=None,
                   help="If set, also run consensus(k=...) and write <name>_results.pkl")
    return p.parse_args(argv)


def main(argv: Iterable[str] | None = None) -> int:
    a = _parse(argv)
    run(
        in_h5ad=a.in_h5ad,
        output_dir=a.output_dir,
        name=a.name,
        obs_filter=a.obs_filter,
        qc_filter=a.qc_filter,
        k_range=range(a.k_min, a.k_max + 1),
        n_iter=a.n_iter,
        n_top_genes=a.n_top_genes,
        seed=a.seed,
        consensus_k=a.consensus_k,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
