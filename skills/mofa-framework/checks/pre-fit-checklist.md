# Pre-fit checklist (walk before calling `ent.run()` / `run_mofa()`)

Run through these eight checks before training any MOFA model. Each item has a one-line verification you can paste directly.

1. **Likelihoods set per view and match data type.** Gaussian for log-normalised continuous (log-CPM, VST, beta-on-logit); bernoulli for 0/1 indicators (binarised ATAC peaks, methylation flags); poisson for raw integer counts.
   ```python
   # The auto-detector (build_model/utils.py:99-115) classifies all-0/1 as bernoulli,
   # all-integer as poisson, else gaussian. Set likelihoods explicitly when in doubt.
   assert set(likelihoods).issubset({"gaussian", "bernoulli", "poisson"})
   for m, lik in enumerate(likelihoods):
       if lik == "gaussian":
           assert not np.all((Y[m] == 0) | (Y[m] == 1) | np.isnan(Y[m])), \
               f"View {m}: gaussian likelihood on 0/1 data — did you mean bernoulli?"
   ```

2. **Per-view centring/scaling applied OR `process_data` will handle it.** Gaussian views only. Bernoulli and Poisson are not centred/scaled.
   ```python
   # If center_features=True and scale_features=True (the defaults) in process_data,
   # mofapy2 will centre and scale gaussian views to unit variance per-feature.
   # Otherwise pre-centre yourself and pass center_features=False, scale_features=True.
   ```

3. **Missing values represented as `NaN` / `NA`, not zero-imputed.** Zero is real signal in counts.
   ```python
   # NaN masking is handled at Y_nodes.py:46-50.
   for m in range(M):
       assert not np.all(Y[m] == 0, axis=1).any(), \
           f"View {m}: some rows are all zeros — confirm not zero-imputed missing"
   ```

4. **Batch / sample / group structure declared via `set_group` before `build()`.** Donor or batch as the group key gives `AlphaZ_Node` something to shrink against.
   ```python
   # In mofapy2:
   ent.set_data_options(use_float32=True)
   ent.set_data_matrix(data=Y_list, likelihoods=lik, groups_names=group_ids)
   # Group identity is recorded at set_data_matrix time — it is structural, not a covariate.
   ```

5. **Seed set explicitly and convergence_mode recorded.** Bayesian VI is seed-sensitive.
   ```python
   ent.set_train_options(seed=42, convergence_mode="medium")
   # convergence_mode ∈ {"fast", "medium", "slow"}. "fast" is for prototyping.
   # For publication-grade fits use "medium" or "slow" + maxiter ≥ 1000.
   ```

6. **GPU availability and `gpu_mode` flag in agreement.** `gpu_mode=True` requires CuPy.
   ```python
   try:
       import cupy as cp
       has_gpu = cp.cuda.runtime.getDeviceCount() > 0
   except (ImportError, RuntimeError):
       has_gpu = False
   assert (gpu_mode and has_gpu) or (not gpu_mode), \
       "gpu_mode=True but CuPy not available — set gpu_mode=False or install cupy"
   ```

7. **SVI options compatible with the rest of the model.** Stochastic mode is incompatible with MEFISTO smooth covariates (`entry_point.py:1037-1042`) and disables factor-dropping (`entry_point.py:1049-1051`).
   ```python
   if stochastic:
       assert not has_smooth_covariates, "SVI incompatible with MEFISTO"
       # drop_factor_threshold will be ignored — K stays fixed.
   ```

8. **Output path planned and write-permission confirmed.**
   ```python
   from pathlib import Path
   outfile = Path("model.hdf5")
   assert outfile.parent.exists() and os.access(outfile.parent, os.W_OK), \
       f"Cannot write to {outfile.parent}"
   ```

## Extending this checklist

When you discover a preventable failure that should have been caught before `.run()` / `run_mofa()`, add a numbered item here. Keep each item to one line + one verification block. Move anything bigger to `references/troubleshooting.md`.
