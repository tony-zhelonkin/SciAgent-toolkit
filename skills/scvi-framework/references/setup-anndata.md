# setup_anndata: the shared registration pattern

Every scvi-tools model begins with `setup_anndata` to register data structure. Run this **before** constructing the model. If `adata` is subset or its columns change, re-run `setup_anndata`.

## Canonical pattern

```python
import scvi
import scanpy as sc

# 1. Raw counts MUST be in .X or a named layer
adata.layers["counts"] = adata.X.copy()

# 2. HVG selection (for scRNA models)
sc.pp.highly_variable_genes(
    adata, n_top_genes=2000, flavor="seurat_v3",
    layer="counts", batch_key="batch", subset=True,
)

# 3. Register
scvi.model.SCVI.setup_anndata(
    adata,
    layer="counts",                                # raw counts location (CRITICAL)
    batch_key="batch",                             # technical batches, categorical str
    categorical_covariate_keys=["donor"],          # optional extra nuisance vars
    continuous_covariate_keys=["percent_mito"],
)
```

## Batch vs sample vs covariate

| Key | Meaning | Where it goes |
|---|---|---|
| `batch_key` | Technical batch to correct out | `setup_anndata(batch_key=...)` |
| `sample_key` | **Biological** donor/sample — modeled, not corrected | MrVI only (`setup_anndata(sample_key=...)`) |
| `categorical_covariate_keys` | Extra nuisance categorical vars | list of obs columns |
| `continuous_covariate_keys` | Continuous nuisance (e.g., %mito) | list of obs columns |
| `labels_key` | Cell-type labels for scANVI | scANVI only |

Batch keys must be **categorical strings**, not integers. Convert with `adata.obs["batch"] = adata.obs["batch"].astype(str)`.

## Required-layer rule

scvi-tools models expect **raw integer counts**. Do not feed normalized or log-transformed values. If you've already normalized in place, reload from source or keep a `.layers["counts"]` copy before normalizing.

## When to re-run setup_anndata

- After `adata = adata[mask].copy()` (subsetting)
- After adding/removing genes
- After changing dtype of a registered obs column
- After concatenating another AnnData

## Sparse matrix format

Ensure CSR: `adata.X = scipy.sparse.csr_matrix(adata.X)`. CSC or dense can cause silent slowdowns or OOM.

## Per-model extras

- **scANVI:** add `labels_key=...` and `unlabeled_category="Unknown"`.
- **MrVI:** add `sample_key=...`, keep `batch_key=...` for technical effects.
- **MultiVI:** requires concatenated RNA+ATAC genes; see the skill.
- **PeakVI / PoissonVI:** peaks-by-cells matrix; no HVG step.
- **totalVI / CITE-seq:** register a `protein_expression_obsm_key`.
