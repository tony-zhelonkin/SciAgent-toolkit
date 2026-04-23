# Training, GPU, and reproducibility

## Global settings

```python
import scvi
import torch

scvi.settings.seed = 42                          # reproducibility
scvi.settings.device = "cuda"                    # auto-detects by default
torch.set_float32_matmul_precision("high")       # faster matmul on modern GPUs
scvi.settings.dl_num_workers = 4                 # parallel data loading; set 0 if OOM
```

## Typical train call

```python
model.train(
    max_epochs=400,
    early_stopping=True,
    early_stopping_patience=45,
    check_val_every_n_epoch=1,
    batch_size=128,
    plan_kwargs={"lr": 1e-3, "weight_decay": 1e-6},
)
```

## When to override defaults

| Situation | Change |
|---|---|
| Small dataset (<20k cells) | `max_epochs=200`, `batch_size=64` |
| Very large dataset (>1M cells) | `batch_size=512–1024`, `dl_num_workers=8` |
| GPU OOM | Lower `batch_size`, set `dl_num_workers=0`, reduce `n_hidden` |
| Unstable loss | Lower `lr` to `5e-4`, increase `weight_decay` |
| scArches query fine-tuning | `max_epochs=200`, `plan_kwargs={"weight_decay": 0.0}` |

## Save / load

```python
model.save("my_model/", save_anndata=True, overwrite=True)
model = scvi.model.SCVI.load("my_model/", adata=adata)   # pass adata if save_anndata=False
```

`save_anndata=True` is convenient but bloats disk; prefer `False` for sharable reference models and keep the AnnData separately.

## Validation loss sanity checks

```python
model.history["elbo_train"].plot()
model.history["elbo_validation"].plot()
```

Train ELBO should decrease monotonically; a widening train/val gap indicates overfitting — raise `weight_decay` or lower `n_latent`.
