#!/usr/bin/env python
"""
export_explorers.py — build compact interactive-explorer tables (utility, not a stage).
=======================================================================================
SKELETON — copy this into `02_analysis/scripts/export_explorers.py` and fill the
`# TODO(project):` stubs. It is NOT symlinked: each project owns its exporter because the
per-stage checkpoints, marker genes, and signatures are project-specific.

Materializes small parquet tables under `03_results/interactive/` that the jscatter breakpoint
explorers (`02_analysis/notebooks/*_explore/`) load instantly — so the live interactive kernel
never opens a multi-GB checkpoint. This is a READ-ONLY projection of already-computed
checkpoints; it recomputes NO biology.

Re-run whenever an upstream checkpoint changes:
    python 02_analysis/scripts/export_explorers.py

Outputs (03_results/interactive/): one `<stage>_explore.parquet` per gated inflection point,
each = x,y coords + a few obs columns + a handful of marker/score columns, indexed by barcode.
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import yaml

# --- root + config -----------------------------------------------------------------------------
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "02_analysis"))
os.chdir(ROOT)

with open(ROOT / "02_analysis/config/analysis_config.yaml") as _fh:
    CONFIG = yaml.safe_load(_fh) or {}
ICFG = CONFIG.get("interactive", {}) or {}
PATHS = CONFIG.get("paths", {}) or {}

OUT = ROOT / PATHS.get("interactive", "03_results/interactive/")
OUT.mkdir(parents=True, exist_ok=True)

# Which embedding in adata.obsm becomes the x,y of the explorer (e.g. "X_umap").
OBSM = ICFG.get("load_coords_obsm", "X_umap")
# TODO(project): marker genes + named signatures to overlay (prefer reading them from config).
MARKERS = list(ICFG.get("marker_genes", []))                 # e.g. ["FOXP3", "CD8A", ...]
SIGNATURES = dict(ICFG.get("signatures", {}) or {})          # {"score_name": [gene, gene, ...]}


# =============================================================================================
# Generic builders — reusable across projects (no project-specific literals here)
# =============================================================================================
def _coords(adata) -> pd.DataFrame:
    """x,y from adata.obsm[OBSM], indexed by barcode."""
    xy = np.asarray(adata.obsm[OBSM])
    return pd.DataFrame({"x": xy[:, 0], "y": xy[:, 1]}, index=adata.obs_names.astype(str))


def _add_obs(df, adata, cols) -> None:
    """Add obs columns; numeric stay numeric, everything else -> str (jscatter needs str/category)."""
    for c in cols:
        if c not in adata.obs:
            continue
        col = adata.obs[c]
        num = (pd.api.types.is_numeric_dtype(col)
               and not isinstance(col.dtype, pd.CategoricalDtype)
               and not pd.api.types.is_bool_dtype(col))
        df[c] = col.to_numpy() if num else col.astype(str).to_numpy()


def _sym_map(adata) -> dict:
    """symbol -> var_name map. TODO(project): adapt if var_names are already gene symbols."""
    if "gene_symbol" not in adata.var:
        return {vn: vn for vn in adata.var_names.astype(str)}   # var_names ARE symbols
    m = {}
    for vn, sym in zip(adata.var_names.astype(str), adata.var["gene_symbol"].astype(str)):
        if sym and sym != "nan" and sym not in m:
            m[sym] = vn
    return m


def _add_genes(df, adata, genes) -> None:
    """Add per-cell expression for each gene present in the object."""
    sym_to_var = _sym_map(adata)
    for g in genes:
        vn = sym_to_var.get(g)
        if vn is None:
            continue
        x = adata[:, vn].X
        df[g] = np.asarray(x.todense()).ravel() if hasattr(x, "todense") else np.asarray(x).ravel()


def _add_modules(df, adata, modules) -> None:
    """Score each gene set via scanpy score_genes (on lognorm X) and add as a column."""
    sym_to_var = _sym_map(adata)
    for name, symbols in modules.items():
        present = [sym_to_var[s] for s in symbols if s in sym_to_var]
        if len(present) < 3:   # score_genes needs a few genes for a stable background
            continue
        sc.tl.score_genes(adata, gene_list=present, score_name=name, use_raw=False, random_state=0)
        df[name] = adata.obs[name].to_numpy()


def _slim(adata):
    """Free the heavy bits not needed for the explorer (raw snapshot + counts layers); keep X,
    var, obs, obsm. Halves peak memory so scoring does not OOM."""
    adata.raw = None
    for k in list(adata.layers.keys()):
        del adata.layers[k]
    return adata


def _stage_h5ad(stage: str) -> Path:
    """Resolve the checkpoint for a stage. TODO(project): map stage id -> its .h5ad path."""
    # Common convention: 03_results/objects/<stage>.h5ad. Adapt to your checkpoint layout.
    objects = ROOT / PATHS.get("objects", "03_results/objects/")
    return objects / f"{stage}.h5ad"


def _export_stage(stage: str, obs_cols, genes, modules) -> None:
    """Build ONE <stage>_explore.parquet: coords + obs + genes + module scores."""
    h5 = _stage_h5ad(stage)
    if not h5.exists():
        print(f"[export] {stage}: checkpoint {h5} not found — skipping")
        return
    adata = _slim(sc.read_h5ad(h5))
    df = _coords(adata)
    _add_obs(df, adata, obs_cols)
    _add_genes(df, adata, genes)
    _add_modules(df, adata, modules)
    out = OUT / f"{stage}_explore.parquet"
    df.to_parquet(out)
    print(f"[export] {out.name} {df.shape}")
    del adata


# =============================================================================================
# main — one _export_stage call per gated inflection point
# =============================================================================================
def main() -> None:
    # TODO(project): list the stages to export + the obs columns each explorer should carry.
    #   Prefer driving `stages` from config (interactive.export_stages) and keeping the per-stage
    #   obs subset here (it is genuinely stage-specific — QC metrics for 01_qc, labels for 02_*).
    stages = ICFG.get("export_stages", []) or ["01_qc"]

    # Per-stage obs subsets. TODO(project): fill for each stage you export.
    obs_by_stage = {
        # "01_qc": ["leiden", "pct_counts_mt", "total_counts", "n_genes_by_counts", "batch"],
        # "02_annotation": ["cell_type", "predicted_label", "batch"],
    }

    for stage in stages:
        obs_cols = obs_by_stage.get(stage, [])
        _export_stage(stage, obs_cols=obs_cols, genes=MARKERS, modules=SIGNATURES)

    print(f"[export] wrote explorer tables to {OUT}")


if __name__ == "__main__":
    main()
