# Model interoperability matrix

All scvi-tools models share the `setup_anndata → train → get_latent_representation → save/load` contract. A few useful conversions between them:

## Transitions

| From | To | How |
|---|---|---|
| scVI | scANVI | `scvi.model.SCANVI.from_scvi_model(scvi_model, labels_key=..., unlabeled_category=...)` |
| scVI | any scArches-compatible | Retrain with `use_layer_norm="both"`, `use_batch_norm="none"`, `encode_covariates=True` |
| MultiVI | PeakVI | Subset adata to ATAC features, re-run `setup_anndata`, train PeakVI |
| PeakVI | MultiVI | Concatenate RNA features, re-run `setup_anndata` with both modalities |
| Any | scArches query | Train reference with scArches flags → `prepare_query_anndata` → `load_query_data` |

## When to chain

- **scVI → scANVI:** any time you have partial labels. Initializing scANVI from scVI is strictly better than training scANVI from scratch.
- **scVI latent → SCENIC+ cell selection:** pass `adata.obsm["X_scVI"]` into pycisTopic workflows.
- **scVI latent → scGLUE RNA input:** use scVI-preprocessed RNA before guidance-graph alignment.
- **scVI → Hub upload:** save with `save_anndata=True`, minify, push with `HubModel`.

## Key rules that make interop work

1. **Raw counts always live in the same layer** across transitions.
2. **`var_names` stay consistent** across ref and query.
3. **Batch keys are strings.** Re-cast after any subset-copy-concat cycle.
4. **Re-register** with `setup_anndata` whenever the AnnData structure changes.
5. **scArches flags at reference training time** are not optional — cannot be retrofitted post hoc.

See `setup-anndata.md` for registration details and `scarches-core.md` for the transfer-learning pattern.
