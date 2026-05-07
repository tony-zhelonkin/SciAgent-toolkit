"""validate_cxg_h5ad.py — pre-deploy schema check.

Run before `docker compose up` to fail fast if Phase A was skipped or
incomplete. Usage:

    python checks/validate_cxg_h5ad.py 03_results/objects/08_explore_full.h5ad

Exits 0 on success, 1 on any failure. Prints a per-check verdict so the
agent (and the user) can see exactly what failed.
"""

from __future__ import annotations

import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd


def _check(name: str, ok: bool, detail: str = "") -> bool:
    tag = "OK  " if ok else "FAIL"
    suffix = f" — {detail}" if detail else ""
    print(f"  [{tag}] {name}{suffix}")
    return ok


def validate(h5ad_path: Path) -> int:
    if not h5ad_path.exists():
        print(f"  [FAIL] file does not exist: {h5ad_path}")
        return 1

    try:
        adata = ad.read_h5ad(h5ad_path)
    except Exception as e:
        print(f"  [FAIL] could not read AnnData: {e}")
        return 1

    print(f"Validating CXG schema for {h5ad_path}")
    print(f"  shape: {adata.n_obs} cells × {adata.n_vars} genes")

    results: list[bool] = []

    # 1. obs index unique + string-typed
    results.append(_check(
        "obs.index unique",
        adata.obs_names.is_unique,
        f"non-unique count: {adata.n_obs - adata.obs_names.nunique()}",
    ))
    results.append(_check(
        "obs.index inferred dtype is 'string'",
        pd.api.types.infer_dtype(adata.obs.index) == "string",
        f"got {pd.api.types.infer_dtype(adata.obs.index)}",
    ))

    # 2. var_names unique + string-typed
    results.append(_check(
        "var_names unique",
        adata.var_names.is_unique,
        f"non-unique count: {adata.n_vars - adata.var_names.nunique()}",
    ))
    results.append(_check(
        "var.index inferred dtype is 'string'",
        pd.api.types.infer_dtype(adata.var.index) == "string",
        f"got {pd.api.types.infer_dtype(adata.var.index)}",
    ))

    # 3. obs['barcode'] present and unique
    results.append(_check(
        "obs['barcode'] present",
        "barcode" in adata.obs,
    ))
    if "barcode" in adata.obs:
        results.append(_check(
            "obs['barcode'] unique",
            bool(adata.obs["barcode"].is_unique),
        ))

    # 4. X_umap shape
    has_umap = "X_umap" in adata.obsm and adata.obsm["X_umap"].shape == (adata.n_obs, 2)
    results.append(_check(
        "obsm['X_umap'] present, shape (n_obs, 2)",
        has_umap,
        f"shape={adata.obsm['X_umap'].shape}" if "X_umap" in adata.obsm else "missing",
    ))

    # 5. uns['title'] set
    results.append(_check(
        "uns['title'] set",
        "title" in adata.uns,
    ))

    # 6. all obsm entries are np.ndarray (no DataFrames left over)
    df_keys = [k for k, v in adata.obsm.items() if isinstance(v, pd.DataFrame)]
    results.append(_check(
        "all obsm entries are numpy arrays",
        not df_keys,
        f"DataFrames remain: {df_keys}" if df_keys else "",
    ))

    # 7. no NaN in numeric obsm entries
    nan_keys = []
    for k, v in adata.obsm.items():
        if isinstance(v, np.ndarray) and np.issubdtype(v.dtype, np.floating):
            if np.isnan(v).any():
                nan_keys.append(k)
    results.append(_check(
        "no NaN in numeric obsm entries",
        not nan_keys,
        f"NaN in: {nan_keys}" if nan_keys else "",
    ))

    # 8. all DataFrame obsm/varm entries' indices match
    drift_keys = []
    for k, v in adata.obsm.items():
        if hasattr(v, "index") and not isinstance(v, np.ndarray):
            if not v.index.equals(adata.obs.index):
                drift_keys.append(("obsm", k))
    for k, v in (getattr(adata, "varm", {}) or {}).items():
        if hasattr(v, "index") and not isinstance(v, np.ndarray):
            if not v.index.equals(adata.var.index):
                drift_keys.append(("varm", k))
    results.append(_check(
        "no obsm/varm DataFrame index drift",
        not drift_keys,
        f"drift in: {drift_keys}" if drift_keys else "",
    ))

    n_pass = sum(results)
    n_total = len(results)
    print(f"\n  {n_pass}/{n_total} checks passed")
    return 0 if n_pass == n_total else 1


def main() -> int:
    if len(sys.argv) != 2:
        print("Usage: python validate_cxg_h5ad.py <path-to-h5ad>")
        return 2
    return validate(Path(sys.argv[1]))


if __name__ == "__main__":
    raise SystemExit(main())
