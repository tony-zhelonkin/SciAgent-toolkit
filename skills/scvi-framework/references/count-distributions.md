# Count distributions

| Distribution | Model default | Use case |
|---|---|---|
| ZINB (zero-inflated NB) | scVI, scANVI, MultiVI, MrVI, contrastiveVI | scRNA-seq with dropout |
| NB | `scVI(gene_likelihood="nb")` | Often better for batch integration; simpler |
| Bernoulli | PeakVI | Binary scATAC-seq accessibility |
| Poisson | PoissonVI | Fragment-count scATAC-seq |
| NB (protein) | totalVI | CITE-seq protein counts |

## Choosing between ZINB and NB for scRNA

- **Default to ZINB** for 10x 3'/5' and Smart-seq2 — explicit dropout term helps with sparse data.
- **Switch to NB** when integration looks poor, UMAP shows batch effects despite correction, or differential expression is unstable. The extra ZINB parameters can overfit.

## Dispersion options

```python
scvi.model.SCVI(adata, dispersion="gene")         # default; one dispersion per gene
scvi.model.SCVI(adata, dispersion="gene-batch")   # per-gene × per-batch
scvi.model.SCVI(adata, dispersion="gene-label")   # per-gene × per-cluster (if labels exist)
```

Use `"gene-batch"` when batch-specific technical variance is large (different protocols, different sequencers). Use `"gene-label"` rarely — it can leak label information into the latent space.
