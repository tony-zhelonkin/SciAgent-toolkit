# Phase A-bis — per-celltype subset re-embedding

The reason to subset: when the wet lab opens a CXG instance for "T cells only", the UMAP should reflect *within-T-cell* heterogeneity (memory vs naive vs effector states), not the full-dataset embedding where T cells form one tight blob.

## The pipeline

Six substantive steps. Parameters lifted from 13403-YD `02_Analysis/04_subcluster.py`.

```python
def reembed_subset(adata_sub: AnnData, leiden_resolution: float = 0.8) -> AnnData:
    # 0. Preserve raw counts
    if "counts" not in adata_sub.layers:
        adata_sub.layers["counts"] = (
            adata_sub.raw.X.copy() if adata_sub.raw is not None else adata_sub.X.copy()
        )

    # 1. Subset-specific HVGs
    sc.pp.highly_variable_genes(adata_sub, n_top_genes=2000, flavor="seurat", subset=False)

    # 2. Detect zero-variance genes BEFORE scaling
    if sparse.issparse(adata_sub.X):
        gene_vars = np.array(
            adata_sub.X.power(2).mean(axis=0) - np.power(adata_sub.X.mean(axis=0), 2)
        ).flatten()
    else:
        gene_vars = np.var(adata_sub.X, axis=0)
    zero_var = (gene_vars == 0) | (np.abs(gene_vars) < 1e-10)
    if zero_var.sum() > 0:
        print(f"   {zero_var.sum()} zero-variance genes — will be NaN after scale, zero-filled")

    # 3. Scale (NaN-tolerant)
    sc.pp.scale(adata_sub, max_value=10, zero_center=True)
    if sparse.issparse(adata_sub.X):
        n_nan = np.isnan(adata_sub.X.data).sum()
        if n_nan:
            adata_sub.X.data = np.nan_to_num(adata_sub.X.data, nan=0.0, posinf=0.0, neginf=0.0)
    else:
        n_nan = np.isnan(adata_sub.X).sum()
        if n_nan:
            adata_sub.X = np.nan_to_num(adata_sub.X, nan=0.0, posinf=0.0, neginf=0.0)
    if n_nan:
        print(f"   ✓ Replaced {n_nan} NaN values with 0 (zero-variance genes)")

    # 4. PCA on HVGs
    sc.tl.pca(adata_sub, n_comps=30, mask_var="highly_variable")

    # 5. Neighbours + UMAP
    sc.pp.neighbors(adata_sub, n_pcs=30, n_neighbors=15)
    sc.tl.umap(adata_sub)

    # 6. Leiden
    sc.tl.leiden(
        adata_sub,
        resolution=leiden_resolution,
        key_added=f"leiden_{leiden_resolution}",
        flavor="igraph",
        n_iterations=2,
        directed=False,
    )
    return adata_sub
```

## Why each parameter

- **`n_top_genes=2000, flavor='seurat'`** — matches the upstream HVG selection in standard scanpy workflows; more genes means slower PCA without obvious gains for the within-celltype UMAP.
- **`subset=False`** — flag HVGs but keep the full gene matrix in `adata_sub.X` for downstream marker plots in CXG.
- **`zero_center=True, max_value=10`** — `max_value` clips outliers (huge counts in low-mean genes); `zero_center` is required for PCA to behave well.
- **`mask_var='highly_variable'`** — modern scanpy API; deprecates the older `use_highly_variable=True` and avoids the corresponding `FutureWarning`.
- **`n_pcs=30, n_neighbors=15`** — defaults that work for 1k–50k cell subsets; tune up for larger.
- **`leiden flavor='igraph', n_iterations=2, directed=False`** — modern leiden flavour; reproducible with `seed=42` set globally.

## Why zero-variance handling is non-negotiable

`sc.pp.scale(adata, zero_center=True)` divides each gene by its standard deviation. A gene with zero variance has std=0 → division by zero → `NaN`. PCA's input check then raises:

```
ValueError: Input contains NaN, infinity or a value too large for dtype('float32').
```

The detection is cheap (one mean + power per gene). The post-scale `np.nan_to_num` converts NaN → 0, which is the *correct* scaled value for a zero-variance gene (no information, no contribution to PCA).

## After re-embed → Phase A on the subset

The subset is now a fresh AnnData with its own embedding; run `prepare()` (the six Phase A steps) on it before writing. This ensures the subset gets unique cell_ids in its own namespace, which avoids cross-subset collisions when the wet lab compares CXG instances side-by-side.

## Naming and write paths

Per `scrna-pipeline-conventions`:

```python
adata_sub.uns["title"] = f"<project> • {celltype_name} Subset — {adata_sub.n_obs} cells"
adata_sub.write_h5ad(f"03_results/objects/09_explore_{celltype_name}.h5ad", compression="gzip")
```

The 09_ prefix matches the reference; `cellranger-multi-to-anndata` writes 00_, QC writes 01–04_, integration writes 05–07_, full CXG-prep writes 08_, subset-CXG writes 09_.
