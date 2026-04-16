# CoReSh-to-GSEA bridge -- deriving novel gene sets from top hits

This is CORESH's most underused output. Each top-ranking GEO dataset implicitly defines a gene set: the query plus every other gene that moves with the query in that dataset's particular biological context. Extract those co-moving genes and you have a data-driven gene set collection tailored to your signature.

## The math

For a query of size `k` on a dataset with expression matrix `E` (genes x samples, already centered):

```
profile       = colSums(E[query, ]) / 1024                      # direction in sample space
profile_unit  = profile / sqrt(sum(profile^2))                   # normalize
gene_loadings = (E / 1024) %*% profile_unit                      # length = n_genes
```

`gene_loadings[i]` is the projection of gene `i` onto the (normalized) query direction. Large |loading| = this gene varies together with the query genes in this dataset.

Normalizing `profile` to unit length is important: without it, loadings are scaled by `sqrt(queryVar)` which varies wildly across datasets and breaks comparability.

## Why it works without running PCA

CORESH's `E1024` is already centered and often PC-reduced during preprocessing. So the "direction defined by the query" is already expressible in PC space, and projecting all genes onto it is just a matrix-vector product -- no fresh eigendecomposition.

## The recipe (wrapped in `scripts/extract_gene_loadings.R`)

```r
extract_gene_loadings <- function(chunk_path, gse_id, query, n_top = 50) {
  chunk <- qs_read(chunk_path)
  obj <- Filter(function(o) o$gseId == gse_id, chunk)[[1]]
  if (is.null(obj)) stop(sprintf("GSE %s not in %s", gse_id, chunk_path))

  query_idxs <- na.omit(match(query, obj$rownames))
  stopifnot(length(query_idxs) >= 3)

  E <- obj$E1024 / 1024
  profile <- colSums(E[query_idxs, , drop = FALSE])
  profile <- profile / sqrt(sum(profile^2))
  gene_loadings <- as.numeric(E %*% profile)

  # rank all genes by absolute loading
  ord <- order(abs(gene_loadings), decreasing = TRUE)
  top_entrez <- obj$rownames[head(ord, n_top)]
  data.table(
    entrez   = top_entrez,
    loading  = gene_loadings[head(ord, n_top)],
    rank     = seq_len(min(n_top, length(ord)))
  )
}
```

## From top-hit list to GMT file

```r
# top_ranking is the head of a coresh_batch.R output, with one row per query x GSE
derived_sets <- lapply(seq_len(nrow(top_ranking)), function(i) {
  row <- top_ranking[i]
  # find which chunk this GSE lives in (indexed during coresh_batch.R)
  chunk_path <- resolve_chunk_for_gse(row$gse)
  loadings <- extract_gene_loadings(chunk_path, row$gse,
                                     queries[[row$query_name]])
  symbols <- ent2sym_human(loadings$entrez)     # see gene-id-conversion.md
  list(name = sprintf("CORESH_%s_%s", row$query_name, row$gse),
       genes = unique(na.omit(symbols)))
})

# filter and write GMT
derived_sets <- Filter(function(s) length(s$genes) >= 15 &&
                                   length(s$genes) <= 500, derived_sets)
derived_sets <- dedupe_jaccard(derived_sets, threshold = 0.8)   # see below

write_gmt <- function(sets, path) {
  lines <- sapply(sets, function(s) {
    paste(c(s$name, "-", s$genes), collapse = "\t")
  })
  writeLines(lines, path)
}
write_gmt(derived_sets, "03_results/cdc1_path/coresh/coresh_derived_sets.gmt")
```

The resulting `.gmt` file is directly consumable by `bulk-rnaseq-gsea-custom-db` via its GMT ingestion path.

## Caveats (also in SKILL.md pitfalls)

### Sign ambiguity

`gene_loadings` is signed -- positive means "moves with the query," negative means "moves against." BUT the overall sign is arbitrary (like PC1 loadings can flip between PCA runs). Treat derived sets as **unsigned** for GSEA by default (both `up` and `down` directions are tested). Don't try to merge signed leading edges across top hits.

If you *must* have a signed set, anchor on one "seed" dataset, flip others to match its query-gene-mean-loading sign, and accept that this is a heuristic.

### Deduplication

Top-ranking datasets often represent the same biology (two GSEs of PGE2-treated macrophages, say). Their derived sets overlap heavily. Dedupe by Jaccard similarity:

```r
dedupe_jaccard <- function(sets, threshold = 0.8) {
  keep <- rep(TRUE, length(sets))
  for (i in seq_along(sets)) {
    if (!keep[i]) next
    for (j in seq_len(i - 1)) {
      if (!keep[j]) next
      a <- sets[[i]]$genes; b <- sets[[j]]$genes
      jacc <- length(intersect(a, b)) / length(union(a, b))
      if (jacc > threshold) { keep[i] <- FALSE; break }
    }
  }
  sets[keep]
}
```

### Size filter

- Drop sets with `< 15` genes (too small for stable GSEA ES).
- Drop sets with `> 500` genes (too broad; diluted signal).
- Drop sets that are >90% overlap with the original query (trivial -- you already know the query).

### Provenance

Always keep a side-table mapping `CORESH_<query>_<GSE>` -> `(query_name, GSE, chunk_path, loading_cutoff, rank_in_coresh)`. When a derived set enriches in your dataset, you'll want to trace it back to the biological context that produced it.

## Downstream: handoff to `bulk-rnaseq-gsea-custom-db`

The GMT file fits the skill's T2G/T2N ingestion pattern:

```r
# (inside bulk-rnaseq-gsea-custom-db pipeline)
library(clusterProfiler)
T2G <- read.gmt("03_results/cdc1_path/coresh/coresh_derived_sets.gmt")
# T2G columns: term, gene
```

That skill handles the rest (GSEA execution, result normalization, master table assembly).

## When this bridge fails

- **Derived sets all look the same across queries:** your query genes are too housekeeping-like; the projection onto any dataset returns similar gene rankings. Filter query before the sweep (`references/query-design.md`).
- **Derived sets don't enrich in your own data:** either (a) your DE ranking is noisy (low power), (b) the biological context of the top GEO hit doesn't actually match your system, or (c) the query was too broad. Sanity-check the top GSE metadata before trusting its derived set.
