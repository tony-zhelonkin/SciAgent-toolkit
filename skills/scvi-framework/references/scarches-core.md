# scArches core: shared query-to-reference pattern

`scvi-scarches-reference-mapping` is the full skill for query mapping. This file captures the **shared substrate** every scvi-tools model uses.

## Reference training flags (mandatory)

```python
model = scvi.model.SCVI(
    ref_adata,
    use_layer_norm="both",         # REQUIRED
    use_batch_norm="none",         # REQUIRED
    encode_covariates=True,        # REQUIRED
)
model.train()
model.save("ref_model/")
```

Without these three flags, the reference cannot accept query fine-tuning.

## Query projection

```python
scvi.model.SCVI.prepare_query_anndata(query_adata, "ref_model/")
query_model = scvi.model.SCVI.load_query_data(query_adata, "ref_model/")

query_model.train(
    max_epochs=200,
    plan_kwargs={"weight_decay": 0.0},   # freeze bio variation, adapt only batch
)

query_adata.obsm["X_scVI"] = query_model.get_latent_representation()
```

`prepare_query_anndata` pads missing genes with zeros and aligns `var_names`. `load_query_data` attaches a fresh encoder layer ("architectural surgery") while freezing the reference backbone.

## Supported models

scVI, scANVI, totalVI, MultiVI, PeakVI, trVAE — all inherit the same contract. MrVI, contrastiveVI, AmortizedLDA do not currently support scArches.

## Label transfer

For scANVI-trained references, labels propagate automatically via `query_model.predict()`. For scVI references, use KNN or weighted-KNN in latent space against the reference.

## When not to use scArches

- Reference wasn't trained with the three required flags → retrain it.
- Query has a completely disjoint biology (no cell types overlap) → train from scratch, don't map.
- Reference was trained on genes that don't appear in query at all → need a shared gene panel.

See the `scvi-scarches-reference-mapping` skill for advanced surgery, weighted-KNN label transfer, and novel-state detection.
