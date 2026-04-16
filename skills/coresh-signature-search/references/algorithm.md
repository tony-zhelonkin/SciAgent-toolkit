# CORESH algorithm -- what `pctVar` actually measures

## The score

For a query gene set of size `k` and a dataset with expression matrix `E` (genes x samples, centered), CORESH computes

```
queryVar = sum( colSums(E[query, ])^2 )    # variance captured by the query direction
pctVar   = queryVar / k / totalVar * 100
```

where `totalVar = sum(E^2)` is pre-computed once per dataset and stored on the object.

### Why this is "PCA-inspired" but not PCA

Ordinary PCA finds the direction in sample space that **maximizes** `sum(colSums(E[subset, ])^2)` over all possible gene subsets. CORESH fixes the direction to "the sum of the query genes" and measures how much variance that fixed direction captures relative to the total.

The key intuition: if the query genes are strongly co-regulated in a dataset, they share a common axis of variation, and summing them produces a long (high-variance) sample vector. If they are uncorrelated, the sum averages out to short (low-variance). `pctVar` normalizes by `k` so signatures of different sizes can be compared.

### The `E1024` representation

On-disk, the expression matrix is stored as `E1024 = round(E * 1024)` (integer, for storage efficiency). Every use divides by 1024 to recover `E`. Additionally, `E` on these objects is **already centered and optionally PCA-reduced** -- the preprocessing step on the Sergushichev end collapses redundant sample directions before chunking. This has two consequences:

1. You never need to center again in your own code.
2. The CoReSh-to-GSEA bridge (`references/coresh-to-gsea-bridge.md`) is a simple projection of all genes onto the query direction -- not a fresh PCA -- because the data is already in PC space.

### GESECA p-values

Variance-only (`pctVar`) gives a **relative** ranking but no specificity test. For specificity, CORESH calls `fgsea:::gesecaCpp()`, the C++ core of the GESECA algorithm (Gene Set Co-regulation Analysis, part of `fgsea` from the same lab). GESECA constructs a null by sampling gene sets of the same size `k` from the dataset and returns an empirical p-value for the observed `queryVar`.

Call signature:

```r
fgsea:::gesecaCpp(E, queryVar, k, sampleSize = 21, seed = 1, eps = 1e-300)[[1]]$pval
```

- `sampleSize = 21`: the number of random gene-set draws per null iteration (the CORESH default).
- `eps`: the smallest representable p-value; 1e-300 is effectively "no floor."
- Cost: ~0.1-1 s per dataset, vs ~milliseconds for variance-only. Across 44,000 human datasets on 8 cores: variance-only ~10-20 s, with p-values ~2-5 min.

### When to use p-values

- **Screen first with `pctVar` only.** It's 10-100x faster and gives a usable ranking.
- **Re-rank the top ~200 hits by p-value** if `pctVar`-only has too many false positives (cancer cell lines, housekeeping-adjacent datasets -- see `references/interpretation-protocol.md`).
- **Always use p-values for formal reporting** -- the paper quotes GESECA p-values, not raw `pctVar`, for Figure 2.

### Relationship to fgsea

CORESH is roughly equivalent to running `fgsea::geseca()` on the user's query across all GEO datasets in the compendium. The package exists because (a) the preprocessed chunks are a specific `qs2`-packaged format optimized for repeated batch queries, (b) the scoring is simplified and fixed (no user tuning of sample size, normalization, etc.), and (c) the `pctVar` score is a cheap screening step GESECA by itself does not expose.

## Per-object attributes

Each object in a `*_full_objects.qs2` chunk has:

| Attribute | Type | Meaning |
|---|---|---|
| `gseId` | character | GEO Series accession, e.g. `"GSE12345"` |
| `gplId` | character | GEO Platform accession, e.g. `"GPL570"` |
| `E1024` | integer matrix (genes x samples) | `round(centered_expression * 1024)` |
| `rownames` | integer vector | Entrez Gene IDs, one per row of `E1024` |
| `totalVar` | numeric scalar | `sum((E1024/1024)^2)`, pre-computed |

Column names (sample IDs) are not preserved in the compendium -- CORESH only needs the variance structure, not sample identity. If you want per-sample metadata for a top hit, fetch it from GEO directly (e.g. via `GEOquery::getGEO` or the `acc.cgi` endpoint).
