# biomart Ensembl→symbol annotation

Two distinct calls into scanpy's biomart wrapper — one for symbols + biotype, one for the canonical mitochondrial gene list. Cache lives at `~/.cache/scanpy/biomart/`.

## The annotation call

```python
import scanpy as sc

annotations = sc.queries.biomart_annotations(
    organism="mmusculus",                          # "hsapiens" for human
    ["ensembl_gene_id", "external_gene_name", "gene_biotype"],
    use_cache=True,
).set_index("ensembl_gene_id")
```

Maps `var_names` (Ensembl IDs) → `external_gene_name` (symbol) and `gene_biotype` (e.g., `protein_coding`, `lncRNA`, `pseudogene`). The cache file is keyed by organism + attribute list; a different attribute list is a different cache.

## The mitochondrial-gene call

```python
mt_genes = sc.queries.mitochondrial_genes(organism="mmusculus", attrname="ensembl_gene_id")
mt_gene_ids = set(mt_genes["ensembl_gene_id"])
adata.var["mt_biomart"] = adata.var_names.isin(mt_gene_ids)
```

`mt_biomart` is the authoritative MT flag — it does not depend on the symbol prefix (`mt-`/`MT-`), so it survives projects where biomart returned IDs instead of symbols. Downstream QC should prefer `mt_biomart` over `gene_name.str.startswith("mt-")` when both are present.

## Cache invalidation

Symptoms of a stale cache: known protein-coding genes (e.g., `Pdcd1`, `Ifng`) come back as Ensembl IDs; the matched count is suspiciously low (<70%); a fresh wet-lab annotation that should have a known symbol does not.

Resolution:

```bash
ls ~/.cache/scanpy/biomart/                      # inspect
rm ~/.cache/scanpy/biomart/<organism>_<...>      # selective delete
```

Or call `use_cache=False` for one run; the result is written back to cache, so this is also a "force-refresh" pattern.

## Fallback when biomart is unreachable

The reference helper wraps the calls in a try/except; on failure it sets `gene_name = gene_id` (Ensembl IDs preserved as symbols) and `mt_biomart = False`. Downstream QC over MT% then falls back to symbol-prefix detection only — in mouse, `mt-` symbols persist even with old caches because they are stable in the Ensembl mouse reference.

## Cross-reference with downstream skills

- `single-cell-rna-qc` reads `adata.var["mt_biomart"]` if present, falling back to symbol prefix.
- `consensus-nmf-multirun` filters out genes where `gene_biotype == "pseudogene"` before HVG selection (configurable).
- `scvi-basic` is biotype-agnostic but expects `gene_name` non-null for the per-gene loadings tables it produces.

## Performance

`use_cache=True` resolves ~30k mouse genes in <2s after the first call; cold-cache calls take 10–30s depending on biomart load. For a 5-pool project the call happens once after concat — not per-pool — so the cache only matters across sessions. Across-project caches benefit from a shared `~/.cache/scanpy/biomart/`; in a Docker container, mount it to the host: `-v $HOME/.cache/scanpy:/home/devuser/.cache/scanpy`.
