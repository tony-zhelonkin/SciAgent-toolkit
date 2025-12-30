# scArches: Query-to-Reference Mapping

**Foundation:** See `scvi-framework.md` for installation and core patterns.

## When to Use scArches

- Map query data to existing reference atlas
- Transfer labels without retraining from scratch
- Identify novel cell states not in reference
- Collaborative atlas building

**Two options:**
- `scvi-tools` (native): scVI, scANVI, totalVI, MultiVI
- `scarches` package: trVAE, scGen, expiMap, treeArches, scPoli

---

## Model Selection

```
Cell type annotation?
├─ With labels → scANVI
└─ Without labels → scVI + KNN transfer

Novel cell states?
└─ treeArches (see treearches-hierarchy-learning.md)

Multimodal?
├─ CITE-seq → totalVI
└─ Multiome → MultiVI
```

---

## Critical scArches Parameters

**MUST use these when training reference model:**

```python
model = scvi.model.SCVI(
    ref_adata,
    use_layer_norm="both",       # REQUIRED
    use_batch_norm="none",       # REQUIRED
    encode_covariates=True,      # REQUIRED
    n_layers=2,
    n_latent=30,
)
```

Without these, surgery will fail or produce poor results.

---

## Core Workflow (scvi-tools Native)

```python
import scvi
import anndata

# ========================================
# 1. TRAIN REFERENCE (with scArches params)
# ========================================
scvi.model.SCVI.setup_anndata(ref_adata, layer="counts", batch_key="batch")

ref_model = scvi.model.SCVI(
    ref_adata,
    use_layer_norm="both",
    use_batch_norm="none",
    encode_covariates=True,
    n_latent=30,
)
ref_model.train()
ref_model.save("ref_model/")

# ========================================
# 2. PREPARE QUERY
# ========================================
# Validates and pads missing genes
scvi.model.SCVI.prepare_query_anndata(query_adata, "ref_model/")

# ========================================
# 3. LOAD QUERY MODEL (surgery)
# ========================================
query_model = scvi.model.SCVI.load_query_data(query_adata, "ref_model/")

# ========================================
# 4. FINE-TUNE (weight_decay=0 preserves reference)
# ========================================
query_model.train(max_epochs=200, plan_kwargs={"weight_decay": 0.0})

# ========================================
# 5. GET JOINT EMBEDDING
# ========================================
query_adata.obsm["X_scVI"] = query_model.get_latent_representation()

# For combined visualization:
full_adata = anndata.concat([query_adata, ref_adata])
full_adata.obsm["X_scVI"] = query_model.get_latent_representation(full_adata)
```

---

## scANVI Reference Mapping (with Labels)

```python
# 1. Train reference scANVI
scvi.model.SCVI.setup_anndata(ref_adata, layer="counts", batch_key="batch")
ref_scvi = scvi.model.SCVI(ref_adata, use_layer_norm="both",
                            use_batch_norm="none", encode_covariates=True)
ref_scvi.train()

ref_scanvi = scvi.model.SCANVI.from_scvi_model(ref_scvi, labels_key="cell_type",
                                                unlabeled_category="Unknown")
ref_scanvi.train(max_epochs=20)
ref_scanvi.save("ref_scanvi/")

# 2. Map query
scvi.model.SCANVI.prepare_query_anndata(query_adata, "ref_scanvi/")
query_scanvi = scvi.model.SCANVI.load_query_data(query_adata, "ref_scanvi/")
query_scanvi.train(max_epochs=100, plan_kwargs={"weight_decay": 0.0})

# 3. Predict labels
query_adata.obs["predicted"] = query_scanvi.predict()
query_adata.obsm["X_scANVI"] = query_scanvi.get_latent_representation()
```

---

## Label Transfer via Weighted KNN

For scVI (no built-in classifier):

```python
from sklearn.neighbors import KNeighborsClassifier

# Get embeddings
ref_latent = ref_model.get_latent_representation(ref_adata)
query_latent = query_model.get_latent_representation(query_adata)

# Train KNN on reference
knn = KNeighborsClassifier(n_neighbors=50, weights="distance")
knn.fit(ref_latent, ref_adata.obs["cell_type"])

# Predict query
query_adata.obs["knn_predicted"] = knn.predict(query_latent)
query_adata.obs["knn_proba"] = knn.predict_proba(query_latent).max(axis=1)
```

---

## Detecting Novel Cell Types

Low prediction confidence may indicate novel populations:

```python
# From scANVI
probs = query_scanvi.predict(soft=True)
query_adata.obs["confidence"] = probs.max(axis=1)
query_adata.obs["is_novel"] = query_adata.obs["confidence"] < 0.5

# Subcluster novel cells
novel_adata = query_adata[query_adata.obs["is_novel"]]
sc.pp.neighbors(novel_adata, use_rep="X_scANVI")
sc.tl.leiden(novel_adata, resolution=0.3)
```

---

## scarches Package (Additional Models)

For models not in scvi-tools:

```python
import scarches as sca

# Example: trVAE for strong batch effects
sca.models.TRVAE.setup_anndata(ref_adata, batch_key="batch")
model = sca.models.TRVAE(ref_adata)
model.train()
model.save("trvae_ref/")

# Query mapping
query_model = sca.models.TRVAE.load_query_data(query_adata, "trvae_ref/")
query_model.train(max_epochs=200)
```

See `treearches-hierarchy-learning.md` for novel cell type detection with hierarchy.

---

## HLCA (Human Lung Cell Atlas) Example

```python
from scvi.hub import HubModel

# Load HLCA reference
hmo = HubModel.pull_from_huggingface_hub(
    repo_name="scvi-tools/hlca-human-lung-scvi",
    revision="main"
)
ref_model = hmo.model
ref_adata = hmo.adata

# Map your query
scvi.model.SCVI.prepare_query_anndata(query_adata, ref_model)
query_model = scvi.model.SCVI.load_query_data(query_adata, ref_model)
query_model.train(max_epochs=200, plan_kwargs={"weight_decay": 0.0})

query_adata.obsm["X_scVI"] = query_model.get_latent_representation()
```

---

## Critical Gotchas

| Issue | Solution |
|-------|----------|
| Surgery fails | Reference must use `use_layer_norm="both"`, `use_batch_norm="none"` |
| Gene mismatch | `prepare_query_anndata` handles this; genes padded with zeros |
| Poor mapping | Use `weight_decay=0.0` to preserve reference embedding |
| New batches wrong | Ensure `encode_covariates=True` in reference |
| Wrong var_names order | Query genes reordered automatically by `prepare_query_anndata` |

---

## Resources

- **scvi-tools:** https://docs.scvi-tools.org/en/stable/tutorials/notebooks/scrna/scarches_scvi_tools.html
- **scarches:** https://docs.scarches.org/
- **HLCA:** https://www.nature.com/articles/s41591-023-02327-2
- **Paper:** https://www.nature.com/articles/s41587-021-01133-0
