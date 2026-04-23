#!/usr/bin/env python3
"""Pre-flight check: is this AnnData ready for any scvi-tools model?

Usage:
    python validate_scvi_adata.py path/to/adata.h5ad \
        --batch-key batch \
        [--layer counts] \
        [--labels-key cell_type] \
        [--sample-key donor]

Returns non-zero on any hard failure and prints a summary of warnings.
This is the shared check referenced by scvi-framework/checks/pre-train-checklist.md.
It does NOT train anything and does NOT modify the file.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("adata_path", type=Path, help=".h5ad file to check")
    parser.add_argument("--batch-key", default="batch")
    parser.add_argument("--layer", default=None, help="Layer name with raw counts (default: check adata.X)")
    parser.add_argument("--labels-key", default=None, help="For scANVI")
    parser.add_argument("--sample-key", default=None, help="For MrVI")
    args = parser.parse_args()

    try:
        import anndata as ad
        import numpy as np
        import scipy.sparse as sp
    except ImportError as exc:
        print(f"[FAIL] missing dependency: {exc}", file=sys.stderr)
        return 2

    if not args.adata_path.exists():
        print(f"[FAIL] file not found: {args.adata_path}", file=sys.stderr)
        return 2

    adata = ad.read_h5ad(args.adata_path)

    failures: list[str] = []
    warnings: list[str] = []

    # Raw counts
    X = adata.layers[args.layer] if args.layer else adata.X
    values = X.data if sp.issparse(X) else np.asarray(X).ravel()
    if values.size and not np.all(values == values.astype(int)):
        failures.append("Non-integer values found — scvi-tools models require raw counts.")

    # Batch key
    if args.batch_key not in adata.obs.columns:
        failures.append(f"batch_key='{args.batch_key}' not in adata.obs.columns.")
    else:
        dtype = adata.obs[args.batch_key].dtype.name
        if dtype not in ("category", "object"):
            failures.append(f"batch_key='{args.batch_key}' has dtype '{dtype}'. Cast to str or category.")

    # Sparse format
    if sp.issparse(adata.X) and not sp.isspmatrix_csr(adata.X):
        warnings.append("adata.X is sparse but not CSR. Convert with scipy.sparse.csr_matrix(adata.X).")

    # Optional keys
    for key, label in [(args.labels_key, "labels_key"), (args.sample_key, "sample_key")]:
        if key and key not in adata.obs.columns:
            failures.append(f"{label}='{key}' not in adata.obs.columns.")

    # Size sanity
    if adata.n_obs < 500:
        warnings.append(f"Only {adata.n_obs} cells — training may be unstable.")
    if adata.n_vars > 10_000:
        warnings.append(f"{adata.n_vars} features — consider HVG selection before training.")

    # Report
    print(f"AnnData: {adata.shape[0]} cells x {adata.shape[1]} features")
    for w in warnings:
        print(f"[WARN] {w}")
    for f in failures:
        print(f"[FAIL] {f}")

    if failures:
        return 1
    print("[OK] ready for setup_anndata")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
