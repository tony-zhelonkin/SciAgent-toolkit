# Critical gotchas (all scvi-tools models)

| Symptom | Cause | Fix |
|---|---|---|
| "RuntimeError: Expected integer tensor" | Model got normalized / log-transformed values | Restore raw counts to `.X` or `.layers["counts"]`; re-run `setup_anndata` |
| Query mapping fails with KeyError on genes | Query `var_names` differ from reference | Subset query to reference genes; pad missing with zeros via `scvi.model.SCVI.prepare_query_anndata` |
| Batch correction has no effect | `batch_key` is numeric | Cast to str: `adata.obs["batch"] = adata.obs["batch"].astype(str)` |
| Model ignores biological sample variation | Using `batch_key` for donor in MrVI | MrVI distinguishes `sample_key` (biological) from `batch_key` (technical) |
| Training crashes after `adata = adata[mask].copy()` | Registration is stale | Re-run `setup_anndata` |
| Silent CPU-only training on a GPU box | PyTorch CPU wheel or driver mismatch | Install `scvi-tools[cuda]`; verify with `torch.cuda.is_available()` |
| OOM on a normally-sized dataset | CSC or dense `.X` | Convert to CSR: `adata.X = scipy.sparse.csr_matrix(adata.X)` |
| DataLoader hang / zombie workers | `dl_num_workers > 0` in some containers | `scvi.settings.dl_num_workers = 0` |
| scArches query won't load | Reference was trained without scArches flags | Reference needs `use_layer_norm="both"`, `use_batch_norm="none"`, `encode_covariates=True` |
| Latent collapses to ~1–2 effective dimensions | `n_latent` too small or training too short | Raise `n_latent` to 20–30; ensure `max_epochs ≥ 200` |

## Extending this list

When you encounter a new failure pattern, add a row here and, if recurring across models, a short how-to in `troubleshooting.md`. Prefer extending this table over carving out a new skill.
