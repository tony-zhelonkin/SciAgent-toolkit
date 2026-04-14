---
name: scvi-linearscvi
description: "Trains LinearSCVI (scvi-tools variant with a LINEAR decoder) on raw scRNA-seq counts to produce continuous, signed gene-to-factor LOADINGS via model.get_loadings(), giving interpretable factors that identify which genes drive each latent dimension. Use to extract gene-program signatures as ranked gene lists, find drivers of population structure, or replace PCA-on-logcounts with a count-aware equivalent. Trades model capacity and batch-correction strength for interpretability; does NOT support scArches/transfer learning. Unlike scvi-lda (discrete Dirichlet topics, unsigned topic-gene probabilities), LinearSCVI yields continuous signed loadings closer to PCA. Unlike genenmf-metaprogram-discovery (per-sample NMF + cross-sample consensus meta-programs), LinearSCVI trains one joint model across all cells. For strong batch correction use scvi-basic; for discrete topics use scvi-lda; for multi-sample consensus programs use genenmf-metaprogram-discovery."
license: MIT
metadata:
  skill-author: SciAgent-toolkit
  last-reviewed: 2026-04-14
  category: analysis
  tier: simple
  tags:
  - scvi
  - scvi-tools
  - linearscvi
  - factor-analysis
  - interpretable
  - gene-loadings
  - scrna-seq
  complementary-skills:
  - scvi-framework
  - scvi-basic
  - scvi-lda
  - scanpy
  - anndata
  contraindications:
  - Do not use on normalized / log-transformed data. Raw counts required.
  - Do not use for strong multi-batch correction. Use scvi-basic.
  - Do not use with scArches or reference mapping. LinearSCVI does not support transfer learning.
  version: 1.0.0
  upstream-docs: https://docs.scvi-tools.org/en/stable/user_guide/models/linearscvi.html
---

# LinearSCVI: Interpretable Factor Analysis

**Foundation:** See `scvi-framework.md` for installation and core patterns.

## When to Use LinearSCVI

- Interpretable factors (gene programs, biological modules)
- Understanding which genes drive population structure
- Factor analysis on count data (avoiding log-transform biases)
- Extracting gene signatures for cell states

**Key difference from scVI:** Linear decoder → direct gene-to-factor mappings via `get_loadings()`.

**Limitations:**
- Less capacity than scVI (simpler decoder)
- Weaker batch correction for complex effects
- No scArches/transfer learning support

---

## Quick Start

```python
import scvi

scvi.model.LinearSCVI.setup_anndata(adata, layer="counts", batch_key="batch")

model = scvi.model.LinearSCVI(adata, n_latent=10)
model.train()

# Get interpretable loadings
loadings = model.get_loadings()  # genes × latent dimensions
adata.obsm["X_LinearSCVI"] = model.get_latent_representation()
```

---

## Loadings Matrix

The key output - direct gene weights for each latent factor:

```python
# Get loadings: DataFrame (genes × factors)
loadings = model.get_loadings()

# Top genes for factor 0
factor_0_genes = loadings[0].sort_values(ascending=False).head(20)
print(factor_0_genes)

# Interpret: high-weight genes define what factor 0 represents
# Moving high→low along factor 0 = decreasing expression of these genes
```

---

## Biological Interpretation

```python
import pandas as pd

loadings = model.get_loadings()

# For each factor, get top positive and negative genes
for i in range(model.n_latent):
    factor = loadings[i].sort_values()
    print(f"\n=== Factor {i} ===")
    print("Top positive:", list(factor.tail(10).index))
    print("Top negative:", list(factor.head(10).index))

# Use gene set enrichment on top genes per factor
# High positive weights → genes that increase with factor
# High negative weights → genes that decrease with factor
```

---

## Key Parameters

```python
model = scvi.model.LinearSCVI(
    adata,
    n_latent=10,            # Number of factors (start with 10-20)
    n_hidden=128,           # Encoder hidden size
    n_layers=1,             # Encoder layers
    use_observed_lib_size=True,
)

model.train(
    max_epochs=250,
    early_stopping=True,
)
```

---

## Comparison to PCA

| Aspect | PCA | LinearSCVI |
|--------|-----|------------|
| Input | Log-normalized | Raw counts |
| Likelihood | Gaussian | ZINB/NB |
| Loadings | Linear | Linear |
| Batch correction | No | Yes |
| Sparsity | No | Handles zeros |

---

## Typical Workflow

```python
import scvi
import scanpy as sc

# 1. Prep
adata.layers["counts"] = adata.X.copy()
sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor="seurat_v3",
                            layer="counts", subset=True)

# 2. Train
scvi.model.LinearSCVI.setup_anndata(adata, layer="counts", batch_key="batch")
model = scvi.model.LinearSCVI(adata, n_latent=15)
model.train()

# 3. Get outputs
adata.obsm["X_LinearSCVI"] = model.get_latent_representation()
loadings = model.get_loadings()

# 4. Visualize
sc.pp.neighbors(adata, use_rep="X_LinearSCVI")
sc.tl.umap(adata)

# 5. Color by individual factors
for i in range(model.n_latent):
    adata.obs[f"factor_{i}"] = adata.obsm["X_LinearSCVI"][:, i]
sc.pl.umap(adata, color=[f"factor_{i}" for i in range(5)])
```

---

## Critical Gotchas

| Issue | Solution |
|-------|----------|
| Factors not interpretable | Try different `n_latent`; check gene selection |
| Poor batch correction | Use scVI instead for complex batch effects |
| All factors similar | May need more factors; check for outlier cells |
| Want transfer learning | Use scVI + scArches instead |

---

## Resources

- **Docs:** https://docs.scvi-tools.org/en/stable/api/reference/scvi.model.LinearSCVI.html
- **Tutorial:** https://docs.scvi-tools.org/en/stable/tutorials/notebooks/scrna/linear_decoder.html
- **Paper:** https://academic.oup.com/bioinformatics/article/36/11/3418/5807606
