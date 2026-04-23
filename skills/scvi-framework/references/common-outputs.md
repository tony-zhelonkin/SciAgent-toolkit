# Common outputs across scvi-tools models

## Latent representation

```python
latent = model.get_latent_representation()          # posterior mean, shape (n_cells, n_latent)
latent_sample = model.get_latent_representation(give_mean=False)   # stochastic sample
adata.obsm["X_scVI"] = latent
sc.pp.neighbors(adata, use_rep="X_scVI")
sc.tl.umap(adata)
sc.tl.leiden(adata)
```

Use the latent as input to neighbors / UMAP / Leiden instead of PCA.

## Denoised expression

```python
norm = model.get_normalized_expression(library_size=1e4)

# Counterfactual: expression as if all cells were in one batch
norm_b1 = model.get_normalized_expression(transform_batch="batch1")
```

Useful for batch-corrected gene plots and pseudobulk.

## Differential expression (Bayesian)

```python
de = model.differential_expression(
    groupby="cell_type",
    group1="T_cells", group2="B_cells",   # omit for 1-vs-all
    mode="change", delta=0.25,            # tests LFC > delta
    batch_correction=True,
)
markers = de[(de["lfc_mean"] > 0) & (de["bayes_factor"] > 3)]
```

Key columns: `proba_de`, `bayes_factor`, `lfc_mean`, `is_de_fdr_0.05`.

Availability: scVI, scANVI, MultiVI, MrVI, contrastiveVI. Not for PeakVI (use DA).

## Model quality

```python
elbo = model.get_elbo()                             # lower = better (negative ELBO)
recon = model.get_reconstruction_error()
marginal_ll = model.get_marginal_ll(n_mc_samples=100)
```

Use ELBO for model selection only within the same architecture and dataset.
