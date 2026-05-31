# Troubleshooting playbook (append-only)

This file grows as agents encounter novel failure modes across the MOFA family (`mofapy2`, `MOFA2`, `MOFAcellulaR`, `mofax`). The format is deliberately minimal: **symptom → cause → fix**, optionally with a code snippet.

> Append new entries at the bottom. Do not delete or reorder — this is a log, not a reference tree.

---

## Template

```
### <one-line symptom>

**Cause:** <what actually went wrong>
**Fix:** <concrete steps>
**Applies to:** <mofapy2 | MOFA2 | MOFAcellulaR | mofax | all>
**Seen:** <YYYY-MM-DD, skill that referred you here>
```

---

## Entries

### ELBO is non-monotone across iterations

**Cause:** Non-Gaussian views (`bernoulli`, `poisson`) use a Bohning / Jaakkola variational lower bound. Monotonicity of the ELBO is only guaranteed when the variational family is fully conjugate and updates are exact. With the bounds, small ELBO decreases are expected and harmless if the long-run trend is increasing.
**Fix:** Inspect `model@expectations$ELBO` or `model.hdf5/training_stats/elbo`. If the long-run trend is increasing, ignore the local non-monotonicity. If it is consistently decreasing or oscillating, raise `maxiter`, set `convergence_mode = "slow"`, and try `init = "pca"` (`init_model.py:78`). If still pathological, drop the non-Gaussian view temporarily to isolate the cause.
**Applies to:** mofapy2, MOFA2, MOFAcellulaR
**Seen:** 2026-05-28, initial seed (docs/vision/latent/mofa.md §1.1, §2.4)

### Factor pruning surprises me — K shrank from 15 to 3

**Cause:** `drop_factor_threshold` (`entry_point.py:867, 938-948`) drops factors whose total R² is below the threshold across all views; `removeInactiveFactors` (`BayesNet.py:178-220`) executes this at every iteration. If multiple views are small or noisy, ARD shrinks most of the latent dimensions to zero, and the pruner removes them.
**Fix:** If you wanted to keep K fixed, set `drop_factor_threshold = 0` (or `drop_factor_threshold = NULL` in R). If you wanted dynamic K but the pruning is too aggressive, lower the threshold (default is 0.01 in some versions, 0 in others). Check `model@cache$variance_explained$r2_per_factor` to see which factors were pruned and why.
**Applies to:** mofapy2, MOFA2
**Seen:** 2026-05-28, initial seed

### `weight_views = TRUE` produced different factors than `weight_views = FALSE`

**Cause:** The two settings rescale view contributions to Z's update differently (`Z_nodes.py:114-118`). With `weight_views = TRUE`, each view's loss is rescaled by `total_w / (M * Y[m].shape[1])` — the inverse of feature count. This changes the relative pull each view has on Z, and therefore changes which factors emerge as shared vs view-specific.
**Fix:** This is expected behaviour, not a bug. Document your choice in the methods section. Use `weight_views = TRUE` when view sizes are unbalanced (e.g. RNA + methylation) and you want balanced contribution; use `weight_views = FALSE` when one view is intentionally the dominant signal.
**Applies to:** mofapy2, MOFA2
**Seen:** 2026-05-28, initial seed

### NaNs propagate from one view to all factor values for some samples

**Cause:** If a sample has missing values in *all* views (no view has data for it), the variational posterior for that sample's Z row is undefined. MOFA's NaN-masking (`Y_nodes.py:46-50`) integrates out per-element missingness, but it cannot fabricate Z for samples with no observations anywhere.
**Fix:** Before fit, check `np.all(np.isnan(Y[m][i, :]))` for each sample i across all views. Drop samples whose data is entirely missing. If you intended such samples to be "projected later", drop them from the fit and use `MOFAcellulaR::project_data` after.
**Applies to:** mofapy2, MOFA2, MOFAcellulaR
**Seen:** 2026-05-28, initial seed

### GPU OOM at 10⁶ cells with `gpu_mode = True`

**Cause:** Full-batch VI at 10⁶ cells holds the entire Z (N × K) and all W lists on GPU. With K = 30 and float64, Z alone is ~250 MB; the per-view `tau · Y · W` products dominate.
**Fix:** Enable stochastic VI: `ent.set_stochastic_options(batch_size=0.1, learning_rate=0.5, forgetting_rate=0.25, start_stochastic=10)` (`entry_point.py:1033-1064`). Note that stochastic mode disables factor-dropping (`entry_point.py:1049-1051`) and is incompatible with MEFISTO smooth covariates (`entry_point.py:1037-1042`). For 10⁶+ cells consider pseudo-bulking via MOFAcellulaR instead.
**Applies to:** mofapy2 (and downstream)
**Seen:** 2026-05-28, initial seed (docs/vision/latent/mofa.md §3)

### Sign of Factor 1 flipped between two runs with different seeds

**Cause:** The bilinear likelihood `Y ≈ ZW^T` is sign-flip equivariant per factor: replacing `(Z[:, k], W[:, k])` with `(-Z[:, k], -W[:, k])` leaves the likelihood unchanged. Variational inference has no preferred sign; `save_model.py:115-119` reorders factors by total R² but does not pin sign.
**Fix:** Impose a canonical sign before downstream use. The simplest is largest-absolute-loading-positive across views. See `references/geometric-caveats.md` §5.1 for a code sketch. Always document the convention in the figure caption when the sign carries biological meaning.
**Applies to:** all
**Seen:** 2026-05-28, initial seed

### `MOFA2::load_model` fails with "basilisk environment not found"

**Cause:** R's `MOFA2` calls Python `mofapy2` via the `basilisk` package, which manages a private conda env (`MOFA2/R/basilisk.R`). On first use the env must be created; this can fail behind firewalls or when the conda channel is blocked.
**Fix:** Run `MOFA2::basilisk.utils::basiliskRunLocal(env = basilisk::BasiliskEnvironment("mofa_env", pkgname = "MOFA2", packages = c("mofapy2==X.Y.Z")), fun = function() invisible(NULL))` to force env creation. If that fails, install `mofapy2` system-wide via `pip` and set `reticulate::use_python("/path/to/python")` before calling `MOFA2`. Confirm with `reticulate::py_module_available("mofapy2")`.
**Applies to:** MOFA2, MOFAcellulaR
**Seen:** 2026-05-28, initial seed
