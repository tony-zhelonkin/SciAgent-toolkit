# Pretrained models via scvi-hub

```python
from scvi.hub import HubModel

hmo = HubModel.pull_from_huggingface_hub(
    repo_name="scvi-tools/human-lung-atlas",
    revision="main",
    cache_dir="~/.cache/scvi-hub",
)
model = hmo.model
adata = hmo.adata
```

Browse atlases: https://huggingface.co/models?library=scvi-tools

## Typical flow

1. Pick a reference atlas whose organism, tissue, and modality match your query.
2. `HubModel.pull_from_huggingface_hub(...)` downloads checkpoint + minified adata.
3. Use the loaded model directly OR pass it to `load_query_data` for scArches mapping.
4. Optional: fine-tune with scANVI labels from the reference to propagate annotations.

## Minification

References are often "minified" — stripped of raw counts, keeping only what the model needs for inference. Check `adata.uns["_scvi_adata_minify_type"]`.

- Minified models can project queries but cannot retrain from scratch.
- If you need full retraining, pull the unminified companion or the raw `.h5ad`.

## Saving your own Hub model

```python
from scvi.hub import HubMetadata, HubModelCardHelper

model.save("ref_model/", save_anndata=True)
# Then use `scvi-hub` utilities to generate a model card and push to HuggingFace.
```

See the `scvi-hub-models` skill for browsing, selection heuristics, and upload recipes.
