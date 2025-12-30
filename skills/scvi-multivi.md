# MultiVI: RNA + ATAC Multimodal Integration

**Foundation:** See `scvi-framework.md` for installation and core patterns.

## When to Use MultiVI

- Integrating multiome (paired RNA+ATAC) with unimodal data
- Joint embedding of RNA-only, ATAC-only, and paired cells
- Cross-modality imputation
- Differential expression + accessibility testing

**Requirements:** Features must be ordered: genes first, then peaks. Shared peak set across datasets (requires re-calling peaks on merged fragments).

---

## Quick Start

```python
import scvi
import anndata

# Prepare: concatenate [genes | peaks] in feature axis
# adata.var has 'modality' column: "Gene Expression" or "Peaks"

scvi.model.MULTIVI.setup_anndata(
    adata,
    layer="counts",
    batch_key="batch",
    modality_key="modality",  # Distinguishes genes from peaks
)

model = scvi.model.MULTIVI(adata, n_latent=20)
model.train()

adata.obsm["X_MultiVI"] = model.get_latent_representation()
```

---

## Data Preparation

```python
import anndata

# RNA data: genes as features
rna_adata = sc.read_h5ad("rna.h5ad")
rna_adata.var["modality"] = "Gene Expression"

# ATAC data: peaks as features
atac_adata = sc.read_h5ad("atac.h5ad")
atac_adata.var["modality"] = "Peaks"

# Paired multiome: already has both
paired_adata = sc.read_h5ad("multiome.h5ad")
paired_adata.var["modality"] = ["Gene Expression"] * n_genes + ["Peaks"] * n_peaks

# Concatenate cells, union of features
# CRITICAL: genes must come before peaks in var_names
adata = anndata.concat([rna_adata, atac_adata, paired_adata], join="outer")
adata.layers["counts"] = adata.X.copy()
```

---

## Key Model Parameters

```python
model = scvi.model.MULTIVI(
    adata,
    n_hidden=128,
    n_latent=20,             # Latent dimensions
    region_factors=True,     # Learn region-specific scaling (recommended)
    fully_paired=False,      # Set True if ALL cells have both modalities
    n_layers_encoder=2,
    n_layers_decoder=2,
)

model.train(
    max_epochs=500,
    early_stopping=True,
    batch_size=256,
)
```

---

## Outputs

```python
# Joint latent representation
latent = model.get_latent_representation()

# Normalized gene expression (denoised)
norm_expr = model.get_normalized_expression()

# Accessibility probability per peak
accessibility = model.get_accessibility_estimates()

# Impute missing modality
# For RNA-only cells: impute ATAC
# For ATAC-only cells: impute RNA
imputed_rna = model.get_normalized_expression(
    adata=atac_only_adata,
    imputation=True
)
```

---

## Differential Analysis

```python
# Differential expression (genes)
de = model.differential_expression(
    groupby="cell_type",
    group1="cDC1A",
    group2="cDC1B",
    mode="change"
)

# Differential accessibility (peaks)
da = model.differential_accessibility(
    groupby="cell_type",
    group1="cDC1A",
    group2="cDC1B",
    mode="change"
)
```

---

## Critical Gotchas

| Issue | Solution |
|-------|----------|
| Feature order wrong | Genes MUST come before peaks in `var_names` |
| Peaks don't match | Re-call peaks on merged fragment files |
| Missing modality indicator | Each cell needs `modality_key` in `obs` or inferred from data |
| OOM errors | Reduce features; use `fully_paired=True` if applicable |

---

## Resources

- **Docs:** https://docs.scvi-tools.org/en/stable/user_guide/models/multivi.html
- **Tutorial:** https://docs.scvi-tools.org/en/stable/tutorials/notebooks/multimodal/MultiVI_tutorial.html
- **Paper:** https://www.biorxiv.org/content/10.1101/2021.08.11.455920
