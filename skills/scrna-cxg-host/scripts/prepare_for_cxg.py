"""prepare_for_cxg.py — Phase A (schema prepare) and Phase A-bis (subset re-embed).

Ported from <ref-scrna>'s 01_Scripts/Python_scripts/cxg_utils.py and
02_Analysis/04_subcluster.py. The Decision Pauses live in SKILL.md; this
script is the mechanical part — given resolved options, prepare the .h5ad.

Two entry points:

    prepare(in_path, out_path, title)
        # Phase A only — schema prepare on an annotated full dataset

    prepare_subsets(
        in_path, out_dir,
        celltype_column, celltypes,                # Decision Pause 1 resolves these
        leiden_resolution=0.8,
    )
        # Phase A-bis — subset by celltype, independent re-embed, then Phase A
        # on each subset; writes one .h5ad per celltype
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Iterable

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse


# --------------------------------------------------------------------------- #
# Phase A — seven steps (A.0 + A.1..A.6)
# --------------------------------------------------------------------------- #

# pandas extension dtype name -> numpy-native target dtype
# CellxGene 1.2.0 (pandas 1.5.3 + numpy 1.23.5) cannot decode pandas-nullable
# extension arrays at load time. Anything left as Int64/boolean/string survives
# the .h5ad write through the MaskedArray codec and crashes the cellxgene
# server with `TypeError: did not understand one of the types; 'None' not
# accepted`. The mapping below promotes/demotes each extension dtype to the
# closest numpy-native dtype that the older stack can read.
_EXT_DTYPE_TO_NUMPY = {
    "Int8": "float32", "Int16": "float32", "Int32": "float64", "Int64": "float64",
    "UInt8": "float32", "UInt16": "float32", "UInt32": "float64", "UInt64": "float64",
    "Float32": "float32", "Float64": "float64",
    "boolean": "bool",
    "string": "object",
}


def _demote_extension_columns(df: pd.DataFrame, label: str, verbose: bool) -> int:
    """Demote pandas extension dtypes in `df` in place. Returns count of changes."""
    n = 0
    for c in list(df.columns):
        dt = str(df[c].dtype)
        target = _EXT_DTYPE_TO_NUMPY.get(dt)
        if target is None:
            continue
        s = df[c]
        if target == "bool":
            df[c] = s.fillna(False).astype(bool)
        elif target == "object":
            df[c] = s.fillna("").astype(object)
        else:  # float32/float64 — pd.NA round-trips to NaN
            df[c] = s.astype(target)
        n += 1
        if verbose:
            print(f"  {label}[{c!r}]: {dt} -> {target}")
    return n


def _sanitize_dot_columns(df: pd.DataFrame, label: str, verbose: bool) -> int:
    """Replace `.` with `_` in column names of `df` in place. Returns count."""
    renames = {c: c.replace(".", "_") for c in df.columns if "." in c}
    if not renames:
        return 0
    df.rename(columns=renames, inplace=True)
    if verbose:
        print(f"  {label} columns renamed (dots): {renames}")
    return len(renames)


def enforce_cxg_dtypes(
    adata: ad.AnnData,
    sanitize_column_names: bool = True,
    verbose: bool = True,
) -> None:
    """A.0 — Demote pandas extension dtypes in obs/var (and raw.var) to numpy-native.

    CellxGene 1.2.0 ships with `pandas==1.5.3 + numpy==1.23.5` (per the skill's
    pinned Dockerfile). That stack cannot decode pandas-nullable extension
    arrays — `Int64`, `boolean`, `string` — that anndata serialises with the
    MaskedArray codec. The cellxgene loader hits `np.dtype(None)` while
    parsing the column header and the container crash-loops with::

        TypeError: did not understand one of the types; 'None' not accepted

    This step catches the failure mode BEFORE write. It must run before any
    other Phase A step because some downstream helpers (e.g. `final_checks`)
    will themselves choke on `pd.NA`.

    Sources of extension dtypes in the wild:
    - `anndataR` Seurat→AnnData conversion (integer columns with NA → Int64)
    - pandas >= 1.0 default for nullable string columns (`StringDtype`)
    - pyarrow-backed reads from `.parquet` / `.feather` provenance

    Also (when `sanitize_column_names`): replaces `.` with `_` in column
    names. Dotted column names (e.g. `var.features.rank`) have triggered
    cellxgene category-lookup edge cases on older versions.
    """
    n_changes = _demote_extension_columns(adata.obs, "obs", verbose)
    n_changes += _demote_extension_columns(adata.var, "var", verbose)

    if adata.raw is not None:
        raw = adata.raw.to_adata()
        n_raw = _demote_extension_columns(raw.var, "raw.var", verbose)
        if n_raw:
            adata.raw = raw
            n_changes += n_raw

    if sanitize_column_names:
        n_changes += _sanitize_dot_columns(adata.obs, "obs", verbose)
        n_changes += _sanitize_dot_columns(adata.var, "var", verbose)
        if adata.raw is not None:
            raw = adata.raw.to_adata()
            if _sanitize_dot_columns(raw.var, "raw.var", verbose):
                adata.raw = raw

    if verbose:
        print(f"  enforce_cxg_dtypes: {n_changes} column(s) modified")


def make_unique(values: Iterable) -> pd.Index:
    """Add __N suffix to duplicates, preserving order."""
    seen: dict[str, int] = {}
    out: list[str] = []
    for v in map(str, values):
        c = seen.get(v, 0)
        out.append(v if c == 0 else f"{v}__{c}")
        seen[v] = c + 1
    return pd.Index(out)


def convert_obsm_to_arrays(adata: ad.AnnData, verbose: bool = True) -> None:
    """Convert all DataFrame entries in obsm to float32 arrays.

    Column names are stashed in uns[f"{key}_columns"].
    """
    if not getattr(adata, "obsm", None):
        return
    for key in list(adata.obsm.keys()):
        v = adata.obsm[key]
        if isinstance(v, pd.DataFrame):
            if v.columns is not None and len(v.columns) > 0:
                adata.uns[f"{key}_columns"] = list(v.columns)
            adata.obsm[key] = v.values.astype("float32")
            if verbose:
                print(f"  obsm[{key!r}] DataFrame -> float32 array {adata.obsm[key].shape}")


def ensure_unique_varnames(adata: ad.AnnData) -> None:
    """Set var_names from gene_name → gene_id → existing index, with __N dedup.

    Applied to both ``.var`` AND ``.raw.var`` (when ``.raw`` is set). CellxGene
    queries dispatch through ``.raw.X`` whenever ``.raw`` is present, and the
    column lookup uses ``.raw.var.index`` — so failing to reindex ``.raw.var``
    leaves the gene autocomplete returning Ensembl ids even after ``.var`` is
    fixed. The original index is preserved in ``__orig_var_index`` on each side.

    ``.raw`` is read-only, so it must be rebuilt via ``adata.raw.to_adata()``,
    re-indexed in place, and reassigned to ``adata.raw``.
    """
    def _build_new_index(var_df: pd.DataFrame) -> pd.Index:
        if "gene_name" in var_df:
            base = var_df["gene_name"].astype(str)
        elif "gene_id" in var_df:
            base = var_df["gene_id"].astype(str)
        else:
            base = var_df.index.astype(str)
        return make_unique(base)

    # Main .var
    adata.var["__orig_var_index"] = adata.var.index.astype(str)
    adata.var_names = _build_new_index(adata.var)
    adata.var.index = pd.Index(adata.var.index.astype(str))

    # .raw.var — rebuild via to_adata() because .raw is read-only
    if adata.raw is None:
        return
    raw = adata.raw.to_adata()
    raw.var["__orig_var_index"] = raw.var.index.astype(str)
    raw.var_names = _build_new_index(raw.var)
    raw.var.index = pd.Index(raw.var.index.astype(str))
    adata.raw = raw


def ensure_unique_barcode_and_index(
    adata: ad.AnnData, joiner: str = "_", debug: bool = False,
) -> None:
    """Build obs['barcode'] (unique) and obs.index = barcode as string cell_id."""
    if "barcode" in adata.obs:
        adata.obs["barcode_raw"] = adata.obs["barcode"].astype(str)
        base_barcode = adata.obs["barcode_raw"]
    else:
        base_barcode = adata.obs.index.astype(str)
        adata.obs["barcode_raw"] = base_barcode

    if "sample_id" in adata.obs:
        sample_ids = adata.obs["sample_id"].astype(str)
        # Sniff first 10 to detect already-prefixed barcodes
        n_check = min(10, len(base_barcode))
        already = True
        for i in range(n_check):
            bc = base_barcode.iloc[i]
            sid = sample_ids.iloc[i]
            if not (bc.startswith(sid + "-") or bc.startswith(sid + "_") or bc.startswith(sid + ".")):
                already = False
                break
        base = base_barcode if already else (sample_ids + joiner + base_barcode)
    else:
        base = base_barcode

    adata.obs["barcode"] = make_unique(base)
    adata.obs["cell_id"] = adata.obs["barcode"].astype(str)
    adata.obs.index = pd.Index(adata.obs["cell_id"].values, dtype=str, name="cell_id")
    if debug:
        print(f"  cell_id index: dtype={adata.obs.index.dtype}, "
              f"first={adata.obs.index[0]}, unique={adata.obs.index.is_unique}")


def realign_aligned_mappings(adata: ad.AnnData, debug: bool = False) -> None:
    """Replace any DataFrame index in obsm/varm with the new string indices."""
    target_obs = pd.Index(adata.obs.index.astype(str))
    if hasattr(adata, "_obsm") and adata._obsm is not None:
        store = adata._obsm
        d = getattr(store, "_data", store)
        for k in list(d.keys()):
            v = d[k]
            if hasattr(v, "index"):  # DataFrame
                if len(v) != adata.n_obs:
                    print(f"  obsm[{k!r}] length {len(v)} != n_obs {adata.n_obs}; skipping")
                    continue
                df = v.copy()
                df.index = pd.Index(target_obs.values, dtype=str, name="cell_id")
                d[k] = df
            else:                     # ndarray
                if getattr(v, "shape", (0,))[0] != adata.n_obs:
                    print(f"  obsm[{k!r}] shape mismatch; skipping")

    target_var = pd.Index(adata.var.index.astype(str))
    if hasattr(adata, "_varm") and adata._varm is not None:
        store = adata._varm
        d = getattr(store, "_data", store)
        for k in list(d.keys()):
            v = d[k]
            if hasattr(v, "index"):
                if len(v) != adata.n_vars:
                    continue
                df = v.copy()
                df.index = pd.Index(target_var.values, dtype=str)
                d[k] = df


def ensure_umap(adata: ad.AnnData, n_pcs: int = 50) -> None:
    """Compute PCA, neighbors, and UMAP only if missing or mis-shaped."""
    if not ("X_pca" in adata.obsm and adata.obsm["X_pca"].shape[0] == adata.n_obs):
        sc.pp.pca(adata, n_comps=min(n_pcs, max(2, adata.n_vars - 1)))
    has_nbrs = ("neighbors" in adata.uns
                and all(k in adata.obsp for k in ("distances", "connectivities")))
    if not has_nbrs:
        sc.pp.neighbors(adata)
    if not any(k.startswith("X_umap") and adata.obsm[k].shape == (adata.n_obs, 2) for k in adata.obsm):
        sc.tl.umap(adata)


def final_checks(adata: ad.AnnData) -> None:
    """Hard assertions that must pass before write."""
    assert adata.obs_names.is_unique, "obs index not unique"
    assert adata.var_names.is_unique, "var names not unique"
    assert "barcode" in adata.obs and adata.obs["barcode"].is_unique, "obs['barcode'] missing or not unique"
    umap_keys = [k for k in adata.obsm if k.startswith("X_umap")]
    assert umap_keys, "no X_umap* embedding present — cellxgene will show an empty panel"
    for k in umap_keys:
        assert adata.obsm[k].shape == (adata.n_obs, 2), \
            f"obsm['{k}'] is shape {adata.obsm[k].shape}, expected ({adata.n_obs}, 2)"
    assert pd.api.types.infer_dtype(adata.obs.index) == "string", \
        f"obs index dtype must infer as 'string' (got {pd.api.types.infer_dtype(adata.obs.index)})"
    assert pd.api.types.infer_dtype(adata.var.index) == "string", \
        f"var index dtype must infer as 'string' (got {pd.api.types.infer_dtype(adata.var.index)})"
    for k, v in adata.obsm.items():
        if hasattr(v, "index"):
            assert v.index.equals(adata.obs.index), f"obsm[{k!r}] index drift"
            assert pd.api.types.infer_dtype(v.index) == "string", \
                f"obsm[{k!r}] index not string-typed"
        else:
            assert v.shape[0] == adata.n_obs, f"obsm[{k!r}] shape mismatch"


# --------------------------------------------------------------------------- #
# Phase A entry point
# --------------------------------------------------------------------------- #

def prepare(
    in_path: str | Path,
    out_path: str | Path,
    title: str,
    debug: bool = False,
) -> ad.AnnData:
    """Run all six Phase A steps, set title, write checkpoint."""
    in_path, out_path = Path(in_path), Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    adata = sc.read_h5ad(in_path)
    print(f"BEFORE PREPARE: n_obs={adata.n_obs}, n_vars={adata.n_vars}, "
          f"index_unique={adata.obs_names.is_unique}")

    enforce_cxg_dtypes(adata)
    convert_obsm_to_arrays(adata)
    ensure_unique_varnames(adata)
    ensure_unique_barcode_and_index(adata, joiner="_", debug=debug)
    realign_aligned_mappings(adata, debug=debug)
    ensure_umap(adata)
    final_checks(adata)

    adata.uns["title"] = title
    adata.write_h5ad(out_path, compression="gzip")
    print(f"AFTER PREPARE: wrote {out_path} ({adata.n_obs} × {adata.n_vars})")
    return adata


# --------------------------------------------------------------------------- #
# Phase A-bis — per-celltype subset + re-embed + prepare
# --------------------------------------------------------------------------- #

def reembed_subset(
    adata_sub: ad.AnnData,
    leiden_resolution: float = 0.8,
    n_top_genes: int = 2000,
    n_pcs: int = 30,
    n_neighbors: int = 15,
) -> ad.AnnData:
    """Subset-specific HVG → scale → PCA → neighbors → UMAP → leiden.

    Zero-variance genes are detected before scaling; any NaN that scaling
    produces is replaced with 0 (the correct scaled value for std==0).
    """
    # Preserve raw counts
    if "counts" not in adata_sub.layers:
        adata_sub.layers["counts"] = (
            adata_sub.raw.X.copy() if adata_sub.raw is not None else adata_sub.X.copy()
        )

    sc.pp.highly_variable_genes(adata_sub, n_top_genes=n_top_genes, flavor="seurat", subset=False)

    # Detect zero-variance genes BEFORE scaling
    if sparse.issparse(adata_sub.X):
        gene_vars = np.array(
            adata_sub.X.power(2).mean(axis=0) - np.power(adata_sub.X.mean(axis=0), 2)
        ).flatten()
    else:
        gene_vars = np.var(adata_sub.X, axis=0)
    n_zero = int(((gene_vars == 0) | (np.abs(gene_vars) < 1e-10)).sum())
    if n_zero:
        print(f"   {n_zero} zero-variance genes — will be NaN after scale, zero-filled")

    sc.pp.scale(adata_sub, max_value=10, zero_center=True)

    # NaN handling
    if sparse.issparse(adata_sub.X):
        n_nan = int(np.isnan(adata_sub.X.data).sum())
        if n_nan:
            adata_sub.X.data = np.nan_to_num(adata_sub.X.data, nan=0.0, posinf=0.0, neginf=0.0)
    else:
        n_nan = int(np.isnan(adata_sub.X).sum())
        if n_nan:
            adata_sub.X = np.nan_to_num(adata_sub.X, nan=0.0, posinf=0.0, neginf=0.0)
    if n_nan:
        print(f"   replaced {n_nan} NaN values with 0")

    sc.tl.pca(adata_sub, n_comps=n_pcs, mask_var="highly_variable")
    sc.pp.neighbors(adata_sub, n_pcs=n_pcs, n_neighbors=n_neighbors)
    sc.tl.umap(adata_sub)
    sc.tl.leiden(
        adata_sub,
        resolution=leiden_resolution,
        key_added=f"leiden_{leiden_resolution}",
        flavor="igraph",
        n_iterations=2,
        directed=False,
    )
    return adata_sub


def prepare_subsets(
    in_path: str | Path,
    out_dir: str | Path,
    celltype_column: str,
    celltypes: list[str],
    project_name: str = "subset",
    leiden_resolution: float = 0.8,
    out_name_template: str = "09_explore_{celltype}.h5ad",
) -> dict[str, Path]:
    """For each celltype, subset → reembed → Phase A prepare → write."""
    in_path = Path(in_path)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    adata = sc.read_h5ad(in_path)
    if celltype_column not in adata.obs:
        raise KeyError(f"obs has no column {celltype_column!r}; "
                       f"available: {list(adata.obs.columns)[:20]}...")

    out_paths: dict[str, Path] = {}
    for ct in celltypes:
        mask = adata.obs[celltype_column].astype(str) == str(ct)
        if mask.sum() == 0:
            print(f"  WARNING: 0 cells with {celltype_column}={ct!r}; skipping")
            continue
        adata_sub = adata[mask].copy()
        print(f"\n{ct}: {adata_sub.n_obs} cells × {adata_sub.n_vars} genes")
        adata_sub = reembed_subset(adata_sub, leiden_resolution=leiden_resolution)

        # Run Phase A on the subset
        enforce_cxg_dtypes(adata_sub)
        convert_obsm_to_arrays(adata_sub)
        ensure_unique_varnames(adata_sub)
        ensure_unique_barcode_and_index(adata_sub, joiner="_")
        realign_aligned_mappings(adata_sub)
        ensure_umap(adata_sub)
        final_checks(adata_sub)

        adata_sub.uns["title"] = f"{project_name} • {ct} Subset — {adata_sub.n_obs} cells"
        # Sanitise celltype name for filename
        safe_ct = str(ct).replace("/", "_").replace(" ", "_")
        out_path = out_dir / out_name_template.format(celltype=safe_ct)
        adata_sub.write_h5ad(out_path, compression="gzip")
        out_paths[str(ct)] = out_path
        print(f"  wrote {out_path} ({adata_sub.n_obs} × {adata_sub.n_vars})")
    return out_paths


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #

def _parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Prepare AnnData for CellxGene hosting.")
    sub = p.add_subparsers(dest="cmd", required=True)

    pa = sub.add_parser("prepare", help="Phase A — schema prep on full dataset")
    pa.add_argument("--in", dest="in_path", type=Path, required=True)
    pa.add_argument("--out", dest="out_path", type=Path, required=True)
    pa.add_argument("--title", type=str, required=True)
    pa.add_argument("--debug", action="store_true")

    pb = sub.add_parser("subsets", help="Phase A-bis — per-celltype re-embed + prepare")
    pb.add_argument("--in", dest="in_path", type=Path, required=True)
    pb.add_argument("--out-dir", type=Path, required=True)
    pb.add_argument("--celltype-column", type=str, required=True)
    pb.add_argument("--celltypes", nargs="+", required=True)
    pb.add_argument("--project-name", type=str, default="subset")
    pb.add_argument("--leiden-resolution", type=float, default=0.8)
    return p.parse_args()


def main() -> int:
    args = _parse_args()
    if args.cmd == "prepare":
        prepare(args.in_path, args.out_path, title=args.title, debug=args.debug)
    elif args.cmd == "subsets":
        prepare_subsets(
            in_path=args.in_path,
            out_dir=args.out_dir,
            celltype_column=args.celltype_column,
            celltypes=args.celltypes,
            project_name=args.project_name,
            leiden_resolution=args.leiden_resolution,
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
