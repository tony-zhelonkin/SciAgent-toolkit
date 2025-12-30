# scANVI: Semi-Supervised Label Transfer

**Foundation:** See `scvi-framework.md` for installation and `scvi-basic.md` for scVI.

## When to Use scANVI

- Partial or complete cell type labels available
- Transfer labels from reference to query
- Better bio-conservation than scVI alone
- Propagate seed labels from marker gene scoring

**Critical:** Always initialize from trained scVI model.

---

## Core Workflow: scVI → scANVI

```python
import scvi

# 1. Train scVI first (unsupervised base)
scvi.model.SCVI.setup_anndata(adata, layer="counts", batch_key="batch")
scvi_model = scvi.model.SCVI(adata, n_layers=2, n_latent=30)
scvi_model.train()

# 2. Initialize scANVI from scVI
scanvi_model = scvi.model.SCANVI.from_scvi_model(
    scvi_model,
    labels_key="cell_type",
    unlabeled_category="Unknown"  # Label for cells without annotations
)

# 3. Train scANVI (few epochs needed)
scanvi_model.train(max_epochs=20, n_samples_per_label=100)

# 4. Get predictions
adata.obs["predicted"] = scanvi_model.predict()
adata.obsm["X_scANVI"] = scanvi_model.get_latent_representation()
```

---

## Seed Labeling (from Markers)

Propagate labels from confidently identified cells:

```python
import numpy as np

# Score cells with markers
markers = {"T_cell": ["CD3D", "CD3E"], "B_cell": ["CD19", "MS4A1"]}

# Create seed labels (Unknown for most cells)
seed_labels = np.array(["Unknown"] * adata.n_obs)
for ct, genes in markers.items():
    score = adata[:, genes].X.mean(axis=1).A1  # Mean expression
    top_idx = score.argsort()[-50:]  # Top 50 cells
    seed_labels[top_idx] = ct

adata.obs["seed_labels"] = seed_labels

# Train with seeds
scvi.model.SCVI.setup_anndata(adata, layer="counts", labels_key="seed_labels")
scvi_model = scvi.model.SCVI(adata)
scvi_model.train()

scanvi_model = scvi.model.SCANVI.from_scvi_model(scvi_model, unlabeled_category="Unknown")
scanvi_model.train(max_epochs=25)

adata.obs["predicted"] = scanvi_model.predict()
```

---

## Label Transfer Between Datasets

```python
import anndata

# Combine reference (labeled) and query (unlabeled)
ref_adata.obs["source"] = "reference"
query_adata.obs["source"] = "query"
adata = anndata.concat([ref_adata, query_adata])

# Create transfer labels
adata.obs["transfer_labels"] = "Unknown"
ref_mask = adata.obs["source"] == "reference"
adata.obs.loc[ref_mask, "transfer_labels"] = ref_adata.obs["cell_type"].values

# Train and predict
scvi.model.SCVI.setup_anndata(adata, layer="counts", batch_key="source")
scvi_model = scvi.model.SCVI(adata)
scvi_model.train()

scanvi_model = scvi.model.SCANVI.from_scvi_model(
    scvi_model, labels_key="transfer_labels", unlabeled_category="Unknown"
)
scanvi_model.train()

# Query predictions
query_preds = adata[adata.obs["source"] == "query"].obs["predicted"]
```

---

## Prediction Outputs

```python
# Hard predictions (most likely label)
predictions = scanvi_model.predict()

# Soft predictions (probabilities per class)
probs = scanvi_model.predict(soft=True)

# Uncertainty as entropy
entropy = -(probs * np.log(probs + 1e-10)).sum(axis=1)
adata.obs["uncertainty"] = entropy

# Confidence (max probability)
adata.obs["confidence"] = probs.max(axis=1)
```

---

## Training Parameters

```python
scanvi_model = scvi.model.SCANVI.from_scvi_model(
    scvi_model,
    labels_key="cell_type",
    unlabeled_category="Unknown",
    n_layers=1,  # Classifier layers
)

scanvi_model.train(
    max_epochs=20,
    n_samples_per_label=100,  # Balance rare types
    plan_kwargs={
        "classification_ratio": 50  # Weight of classification loss (default)
        # Low (1-10): prioritize reconstruction
        # High (100+): prioritize classification
    }
)
```

---

## Handling Novel Cell Types

Cells with low confidence may be novel types:

```python
probs = scanvi_model.predict(soft=True)
max_prob = probs.max(axis=1)

# Flag uncertain cells
uncertain = max_prob < 0.5
adata.obs["is_uncertain"] = uncertain

# Subcluster uncertain cells to check for novel populations
uncertain_adata = adata[uncertain]
sc.pp.neighbors(uncertain_adata, use_rep="X_scANVI")
sc.tl.leiden(uncertain_adata, resolution=0.3)
```

---

## Common Issues

| Issue | Solution |
|-------|----------|
| Poor predictions | Ensure scVI trained well first |
| All cells same type | Check label balance; use `n_samples_per_label` |
| Novel types missed | Lower confidence threshold; subcluster uncertain |
| Reference/query mismatch | Same genes, similar normalization |

---

## Resources

- **Docs:** https://docs.scvi-tools.org/en/stable/user_guide/models/scanvi.html
- **Seed Labeling:** https://docs.scvi-tools.org/en/stable/tutorials/notebooks/scrna/seed_labeling.html
- **Paper:** https://www.embopress.org/doi/full/10.15252/msb.20209620
