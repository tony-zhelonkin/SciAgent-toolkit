---
name: scvi-lda
description: "Runs AmortizedLDA (scvi-tools) for latent Dirichlet allocation topic modeling on scRNA-seq raw COUNT data, treating cells as documents and genes as words to discover discrete transcriptional programs (topics) and per-cell topic proportions. Use for unsupervised gene-module discovery on RNA, finding intermediate/transitional cell states, and interpretable dimensionality reduction before annotation. Requires raw counts with MT genes removed first. Unlike pycistopic-atac-topic-modeling, which runs LDA on scATAC-seq regions for the SCENIC+ pipeline, this skill is RNA-only. Unlike scvi-linearscvi, which gives continuous gene LOADINGS from a linear decoder, LDA gives discrete topic-gene distributions. Unlike genenmf-metaprogram-discovery (per-sample NMF + consensus meta-programs across patients), AmortizedLDA learns a single shared topic space jointly across all cells. For ATAC topics use pycistopic-atac-topic-modeling; for supervised label transfer use scvi-scanvi."
license: MIT
metadata:
  skill-author: SciAgent-toolkit
  last-reviewed: 2026-04-14
  category: analysis
  tier: simple
  tags:
  - scvi
  - scvi-tools
  - amortized-lda
  - topic-modeling
  - gene-programs
  - unsupervised
  - scrna-seq
  complementary-skills:
  - scvi-framework
  - scvi-linearscvi
  - scanpy
  - anndata
  - genenmf-metaprogram-discovery
  contraindications:
  - Do not use on normalized / log-transformed data. Raw counts required.
  - Do not keep MT genes in the input. Filter them first (they dominate topics).
  - Do not use for batch correction or integration. Use scvi-basic.
  version: 1.0.0
  upstream-docs: https://docs.scvi-tools.org/en/stable/user_guide/models/amortizedlda.html
---

# AmortizedLDA: Topic Modeling for Single-Cell

**Foundation:** inherits shared patterns from `scvi-framework`. Jump to:
- `scvi-framework/references/setup-anndata.md` — raw counts with MT genes removed
- `scvi-framework/references/training-and-gpu.md`
- `scvi-framework/references/gotchas.md`, `scvi-framework/checks/pre-train-checklist.md`

This file covers only the parts that differ for AmortizedLDA (discrete topic modeling).

## When to Use AmortizedLDA

- Discover shared transcriptional programs (topics) across cell types
- Identify gene modules without prior annotation
- Find intermediate/transitional states
- Exploratory analysis before supervised annotation
- Dimensionality reduction to interpretable topic space

**Key insight:** Treats cells as documents, genes as words → discovers latent topics.

---

## Quick Start

```python
import scvi

# CRITICAL: Remove MT genes and use raw counts
adata = adata[:, ~adata.var_names.str.startswith("mt-")]  # Mouse
adata.layers["counts"] = adata.X.copy()

scvi.model.AmortizedLDA.setup_anndata(adata, layer="counts")

model = scvi.model.AmortizedLDA(adata, n_topics=10)
model.train()

# Get topic proportions per cell
topic_prop = model.get_latent_representation()
adata.obsm["X_LDA"] = topic_prop
```

---

## Key Outputs

```python
# 1. Cell topic proportions (Monte Carlo estimate)
topic_prop = model.get_latent_representation()  # cells × topics

# 2. Gene-topic distribution
gene_topic = model.get_feature_by_topic()  # genes × topics

# Store in adata
adata.obsm["X_LDA"] = topic_prop
for i in range(n_topics):
    adata.obs[f"topic_{i}"] = topic_prop[f"topic_{i}"].values
```

---

## Topic Interpretation

```python
import pandas as pd

gene_topic = model.get_feature_by_topic()

# Top genes per topic
for i in range(model.n_topics):
    topic = gene_topic[f"topic_{i}"].sort_values(ascending=False)
    print(f"\n=== Topic {i} ===")
    print(topic.head(15))

# Topic enrichment by cell type
topic_cols = [f"topic_{i}" for i in range(model.n_topics)]
topic_by_celltype = adata.obs.groupby("cell_type")[topic_cols].mean()
print(topic_by_celltype)
```

---

## Visualization

```python
import scanpy as sc

# Option 1: Color expression UMAP by topics
sc.pl.umap(adata, color=[f"topic_{i}" for i in range(n_topics)])

# Option 2: Build UMAP in topic space (use Hellinger distance)
sc.pp.neighbors(adata, use_rep="X_LDA", metric="hellinger", n_neighbors=20)
sc.tl.umap(adata)
adata.obsm["topic_umap"] = adata.obsm["X_umap"].copy()
```

---

## Preprocessing Notes

**Critical:** Remove ubiquitous genes that create trivial topics.

```python
# Remove before LDA:
adata = adata[:, ~adata.var_names.str.startswith("mt-")]  # MT genes
adata = adata[:, ~adata.var_names.str.match("^Rp[sl]")]   # Ribosomal

# Select HVGs (500-2000 typical for LDA)
sc.pp.highly_variable_genes(adata, flavor="seurat_v3", layer="counts",
                            n_top_genes=1000, subset=True)
```

---

## Model Selection (Number of Topics)

```python
# Test range of topics
n_topics_range = [5, 10, 15, 20, 30]
elbos = []

for n in n_topics_range:
    scvi.model.AmortizedLDA.setup_anndata(adata, layer="counts")
    model = scvi.model.AmortizedLDA(adata, n_topics=n)
    model.train(max_epochs=500)
    elbos.append(model.get_elbo())
    print(f"n_topics={n}: ELBO={elbos[-1]:.2f}")

# Higher ELBO = better fit
# But also check biological interpretability of topics
```

---

## Key Parameters

```python
model = scvi.model.AmortizedLDA(
    adata,
    n_topics=10,
    n_hidden=128,
    n_layers=1,
    dropout_rate=0.1,
    cell_topic_prior=None,   # Dirichlet prior (default: 1/n_topics)
    topic_gene_prior=None,
)

model.train(
    max_epochs=1000,
    early_stopping=True,
    batch_size=128,
)
```

---

## Critical Gotchas

| Issue | Solution |
|-------|----------|
| Topics dominated by MT/ribosomal | Remove these genes before training |
| ELBO not decreasing | Normal due to KL annealing; check final value |
| One topic dominates | Decrease `n_topics`; filter genes more aggressively |
| Empty topics | Some topics may be unused; reduce `n_topics` |
| Different results each run | Set `scvi.settings.seed = 42` |

---

## Comparison: AmortizedLDA vs pycisTopic

| Feature | AmortizedLDA | pycisTopic LDA |
|---------|--------------|----------------|
| Data type | scRNA-seq | scATAC-seq |
| Inference | Variational (fast) | MCMC (thorough) |
| Scalability | Millions of cells | ~100k cells |
| GPU | Yes | No |

---

## Resources

- **Docs:** https://docs.scvi-tools.org/en/stable/user_guide/models/amortizedlda.html
- **Tutorial:** https://docs.scvi-tools.org/en/stable/tutorials/notebooks/scrna/amortized_lda.html
