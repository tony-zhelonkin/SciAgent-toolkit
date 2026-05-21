# Robust gene-union concatenation

`anndata.concat(adatas, join='outer', fill_value=0)` does most of what we want — union-of-genes, zero-fill the missing — but it loses gene metadata when a gene is absent from the first input. For multi-pool CellRanger projects where different pools may have been counted against different reference releases, that metadata loss is the failure mode that produces "what is `ENSMUSG00000123456`?" hours later.

The robust pattern (ported from the <ref-scrna> reference `Python_scripts/anndata_utils.py::robust_concat_with_gene_union`) reconstructs a *union* `var` table first, then re-indexes each input to the full gene set before concat.

## The four steps

### 1. Collect the union of gene IDs and per-gene metadata

```python
all_genes = set()
gene_metadata = {}     # gene_id -> dict of var-row values (first-occurrence wins)

for adata in adatas:
    sample_genes = set(adata.var_names)
    all_genes.update(sample_genes)
    for gene in sample_genes:
        if gene not in gene_metadata:
            gene_metadata[gene] = adata.var.loc[gene].to_dict()

all_genes = sorted(all_genes)
```

First-occurrence-wins is acceptable when the sequence of inputs is stable (sorted file paths). For projects where input order is non-deterministic, prefer "highest-quality-source-wins" — usually the input with the most non-null `gene_name` values.

### 2. Build a consistent union `var` table

```python
import pandas as pd
union_var = pd.DataFrame(index=all_genes)
union_var["gene_id"] = union_var.index.astype(str)

# Symbol resolution with a preference order (handles 10x version drift)
gene_names = []
for gene in all_genes:
    meta = gene_metadata.get(gene, {})
    sym = next(
        (str(meta[c]) for c in ("gene_symbol", "gene_symbols", "gene_name", "feature_name")
         if c in meta and meta[c] and not pd.isna(meta[c]) and str(meta[c]) != str(gene)),
        gene,                                  # fallback: keep Ensembl ID
    )
    gene_names.append(sym)
union_var["gene_name"] = gene_names
```

Other metadata columns are ported from `gene_metadata` with a `False`/`"unknown"`/`np.nan` default depending on observed dtype.

### 3. Re-index each input to the full gene set (zero-fill)

For sparse inputs use `scipy.sparse.hstack`; for dense, `np.hstack`:

```python
from scipy import sparse
import numpy as np
import anndata as ad

standardized = []
for adata in adatas:
    missing = [g for g in all_genes if g not in adata.var_names]
    if missing:
        n_cells, n_missing = adata.n_obs, len(missing)
        if sparse.isspmatrix(adata.X):
            zero = sparse.csr_matrix((n_cells, n_missing))
            new_X = sparse.hstack([adata.X, zero])
        else:
            zero = np.zeros((n_cells, n_missing))
            new_X = np.hstack([adata.X, zero])
        new_var = pd.concat([adata.var, union_var.loc[missing]])
        new_adata = ad.AnnData(X=new_X, obs=adata.obs.copy(), var=new_var,
                               obsm=adata.obsm.copy() if adata.obsm else {},
                               uns=adata.uns.copy() if adata.uns else {})
        new_adata.var_names = list(adata.var_names) + missing
    else:
        new_adata = adata.copy()
    new_adata = new_adata[:, all_genes].copy()    # consistent order
    new_adata.var = union_var.copy()
    standardized.append(new_adata)
```

### 4. Concat with `join='outer'` for safety

```python
result = ad.concat(
    standardized,
    axis=0,
    join="outer",      # should be exact; safety net
    fill_value=0,
    merge="first",
    uns_merge="unique",
    label="orig_ident",
    keys=keys,         # list of sample_ids
    index_unique=None, # obs_names already prefixed with sample_id
)
```

## Sparse-vs-dense — the memory cliff

A 50-pool project with 25k–35k genes per pool and ~10% non-overlap means the union is ~40k genes. Zero-filling 250k cells × 4k missing genes = 1B float values. In dense `np.float32` that is 4 GB **per pool**, 200 GB total — enough to OOM a workstation. Sparse-fill is ~0.5 GB total. **Always check that `sparse.isspmatrix(adata.X) is True` before concat** — if any input is dense, convert it: `adata.X = sparse.csr_matrix(adata.X)`.

## When to prefer plain `ad.concat`

If all pools were counted against the *same* reference (same Ensembl release, same chromosome inclusion list), `ad.concat(..., join='inner')` keeps the ~99%-shared gene set without zero-fill — faster, smaller, lossless for shared genes. Use the robust gene-union path only when reference divergence is real or unknown.
