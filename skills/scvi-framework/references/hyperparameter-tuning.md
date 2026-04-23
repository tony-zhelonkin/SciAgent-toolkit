# Hyperparameter tuning (autotune)

```python
from scvi import autotune
from ray import tune

results = autotune.run_autotune(
    model_cls=scvi.model.SCVI,
    data=adata,
    metrics="validation_loss",
    mode="min",
    search_space={
        "model_params": {
            "n_hidden": tune.choice([64, 128, 256]),
            "n_latent": tune.choice([10, 20, 30]),
            "n_layers": tune.choice([1, 2, 3]),
        },
        "train_params": {
            "plan_kwargs": {
                "lr": tune.loguniform(1e-4, 1e-2),
            }
        },
    },
    num_samples=10,
)
```

## Choosing a metric

| Metric | When to use |
|---|---|
| `validation_loss` (ELBO) | Default; proxies reconstruction + KL |
| `Batch correction` (scib) | When downstream goal is integration quality |
| `Bio conservation` (scib) | When downstream goal is preserving known cell-type structure |
| Weighted combination | Production atlases |

## Practical tips

- Start with `num_samples=10` on a small HVG subset to find a plausible neighborhood, then run a focused 5-sample search on full data.
- Fix `n_latent` first (sweep 10/20/30), then tune `n_hidden` and `lr` jointly.
- `plan_kwargs.weight_decay` is usually not worth tuning — `1e-6` is a reasonable default.
- Autotune requires Ray; install via `pip install -U "scvi-tools[autotune]"`.
