# scvi-tools Framework Foundation

## Overview

scvi-tools provides probabilistic models for single-cell omics using PyTorch. Part of the scverse ecosystem (Scanpy, AnnData, MuData).

**Citation:** Cite [scvi-tools manuscript](https://www.nature.com/articles/s41587-021-01206-w) + model-specific paper.

---

## Installation

```bash
pip install -U scvi-tools                    # CPU
pip install -U "scvi-tools[cuda]"            # NVIDIA GPU
pip install -U "scvi-tools[metal]"           # Apple Silicon
pip install -U "scvi-tools[tutorials]"       # With scanpy, seaborn, etc.
```

Optional extras: `hub` (HuggingFace), `autotune` (Ray tuning), `jax`, `regseq`.

---

## Model Selection

```
Have pretrained reference?
├─ YES → scArches (scvi-scarches-reference-mapping.md)
└─ NO → Train new model:
    ├─ scRNA-seq only?
    │   ├─ Label transfer needed? → scANVI (scvi-scanvi.md)
    │   ├─ Multi-sample analysis? → MrVI (scvi-mrvi.md)
    │   ├─ Interpretable factors? → LinearSCVI (scvi-linearscvi.md)
    │   ├─ Topic modeling? → AmortizedLDA (scvi-lda.md)
    │   └─ Basic integration → scVI (scvi-basic.md)
    ├─ CITE-seq? → totalVI
    ├─ Multiome (RNA+ATAC)? → MultiVI (scvi-multivi.md)
    ├─ scATAC-seq only? → PeakVI (scvi-peakvi.md)
    └─ Perturbation study? → contrastiveVI (scvi-contrastivevi.md)
```

| Model | Modality | Key Feature | Skill |
|-------|----------|-------------|-------|
| scVI | scRNA-seq | Batch correction | `scvi-basic.md` |
| scANVI | scRNA-seq + labels | Semi-supervised label transfer | `scvi-scanvi.md` |
| MrVI | scRNA-seq + samples | Sample-level effects, covariate DE | `scvi-mrvi.md` |
| MultiVI | Multiome | RNA + ATAC joint modeling | `scvi-multivi.md` |
| PeakVI | scATAC-seq | Accessibility embedding & DA | `scvi-peakvi.md` |
| contrastiveVI | Perturbation | Salient vs background variation | `scvi-contrastivevi.md` |
| LinearSCVI | scRNA-seq | Interpretable linear decoder | `scvi-linearscvi.md` |
| AmortizedLDA | scRNA-seq | Topic modeling | `scvi-lda.md` |
| scArches | Any | Query-to-reference projection | `scvi-scarches-reference-mapping.md` |

---

## Core Pattern: setup_anndata → Model → train

All scvi-tools models follow the same workflow:

```python
import scvi
import scanpy as sc

# 1. Prepare data - MUST have raw counts
adata.layers["counts"] = adata.X.copy()
sc.pp.normalize_total(adata)  # For visualization only
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor="seurat_v3",
                            layer="counts", batch_key="batch", subset=True)

# 2. Register data structure
scvi.model.SCVI.setup_anndata(
    adata,
    layer="counts",           # Raw counts location (CRITICAL)
    batch_key="batch",        # Technical batches
    categorical_covariate_keys=["donor"],
    continuous_covariate_keys=["percent_mito"]
)

# 3. Create and train
model = scvi.model.SCVI(adata, n_latent=30, n_layers=2)
model.train(max_epochs=400, early_stopping=True)

# 4. Extract and use
adata.obsm["X_scVI"] = model.get_latent_representation()
sc.pp.neighbors(adata, use_rep="X_scVI")
sc.tl.umap(adata)

# 5. Save/load
model.save("my_model", save_anndata=True)
model = scvi.model.SCVI.load("my_model")
```

---

## Critical Gotchas

| Issue | Solution |
|-------|----------|
| **Using normalized data** | Models need RAW COUNTS in `.X` or specified layer |
| **var_names mismatch** | Query must have EXACT same genes as reference |
| **Batch key is numeric** | Must be categorical strings, not integers |
| **Sample vs batch confusion** | MrVI: `sample_key` (biological) ≠ `batch_key` (technical) |
| **setup_anndata after subsetting** | Must re-run `setup_anndata` if adata changes |
| **Sparse matrix issues** | Ensure CSR format: `adata.X = scipy.sparse.csr_matrix(adata.X)` |

---

## GPU & Performance

```python
import scvi
import torch

scvi.settings.seed = 42                          # Reproducibility
scvi.settings.device = "cuda"                    # Force GPU (auto-detected)
torch.set_float32_matmul_precision("high")       # Faster matmul
scvi.settings.dl_num_workers = 4                 # Parallel loading (0 if OOM)
```

---

## scvi-hub: Pretrained Models

Load models from [HuggingFace](https://huggingface.co/scvi-tools):

```python
from scvi.hub import HubModel

hmo = HubModel.pull_from_huggingface_hub(
    repo_name="scvi-tools/human-lung-atlas",
    revision="main"
)
model = hmo.model
adata = hmo.adata
```

Browse: https://huggingface.co/models?library=scvi-tools

---

## scArches: Query-to-Reference Mapping

Project new data onto existing atlas without full retraining:

```python
# 1. Train reference with scArches-compatible params
model = scvi.model.SCVI(
    ref_adata,
    use_layer_norm="both",       # REQUIRED for scArches
    use_batch_norm="none",       # REQUIRED for scArches
    encode_covariates=True       # REQUIRED for scArches
)
model.train()
model.save("ref_model/")

# 2. Map query
scvi.model.SCVI.prepare_query_anndata(query_adata, "ref_model/")
query_model = scvi.model.SCVI.load_query_data(query_adata, "ref_model/")
query_model.train(max_epochs=200, plan_kwargs={"weight_decay": 0.0})

query_adata.obsm["X_scVI"] = query_model.get_latent_representation()
```

See `scvi-scarches-reference-mapping.md` for advanced workflows.

---

## Common Outputs

```python
# Latent representation (for UMAP, clustering)
latent = model.get_latent_representation()

# Denoised expression (batch-corrected)
norm = model.get_normalized_expression(library_size=1e4)

# Differential expression (Bayesian)
de = model.differential_expression(groupby="cell_type", mode="change")
# Key columns: proba_de, bayes_factor, lfc_mean, is_de_fdr_0.05

# Model quality
elbo = model.get_elbo()  # Higher = better
```

---

## Hyperparameter Tuning

```python
from scvi import autotune
from ray import tune

results = autotune.run_autotune(
    model_cls=scvi.model.SCVI,
    data=adata,
    metrics="validation_loss",
    mode="min",
    search_space={
        "model_params": {"n_hidden": tune.choice([64, 128, 256])},
        "train_params": {"plan_kwargs": {"lr": tune.loguniform(1e-4, 1e-2)}}
    },
    num_samples=10
)
```

Can also optimize for scib-metrics: `"Batch correction"`, `"Bio conservation"`, etc.

---

## Count Distributions

| Distribution | Model Default | Use Case |
|--------------|---------------|----------|
| ZINB | scVI, scANVI | scRNA-seq with dropout |
| NB | scVI (`gene_likelihood="nb"`) | Simpler, sometimes better integration |
| Bernoulli | PeakVI | scATAC-seq binary accessibility |
| Poisson | PoissonVI | scATAC-seq counts |

---

## Resources

- **Docs:** https://docs.scvi-tools.org/
- **Tutorials:** https://docs.scvi-tools.org/en/stable/tutorials/
- **API:** https://docs.scvi-tools.org/en/stable/api/
- **Hub Models:** https://huggingface.co/scvi-tools
- **GitHub:** https://github.com/scverse/scvi-tools
