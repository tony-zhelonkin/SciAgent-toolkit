---
name: coresh-signature-search
description: "CORESH -- rank ~80,000 public GEO datasets by how strongly a user-supplied gene signature is co-regulated in each, returning a PCA-inspired pctVar score and optional GESECA p-value for hypothesis generation and interpretation. Use when you have a short gene signature (10-200 genes: top DE hits, cluster markers, leading edge of a pathway) and want to find public datasets where those genes move together -- to connect a signature to biological context (cell types, perturbations, diseases), find perturbation analogs for wet-lab follow-up, or derive novel coregulation-based gene sets for downstream GSEA. For canonical-pathway interpretation of a signature or running GSEA on CORESH-derived gene sets use bulk-rnaseq-gsea; for metadata-based GEO search (keyword, organism) use the GEO web interface or RummaGEO."
license: MIT
metadata:
  skill-author: SciAgent-toolkit
  last-reviewed: 2026-04-14
  version: 0.1.0
  upstream-docs: https://alserglab.wustl.edu/coresh/
  category: analysis
  tier: standard
  tags:
    - coresh
    - geo
    - signature-search
    - coregulation
    - geseca
    - fgsea
    - hypothesis-generation
    - public-data-reuse
    - gene-set-discovery
    - synapse
  complementary-skills:
    - bulk-rnaseq-gsea
    - bulk-rnaseq-pathway-explorer
    - gatom-metabolomic-predictions
  contraindications:
    - "Do not use for MSigDB canonical-pathway interpretation of a signature. Use bulk-rnaseq-gsea instead."
    - "Do not use for metadata-based GEO search (organism, disease keyword, platform). Use the GEO web interface or RummaGEO -- CORESH is data-driven only."
    - "Do not pass gene symbols. The chunk rownames are integer Entrez IDs and match() silently returns all NAs on character input."
    - "Do not use human queries against mmu/ chunks (or vice versa). Species mismatch returns near-zero sizes with no error."
---

# CORESH -- Gene Signature Search over Public GEO Data

## Overview

CORESH (COregulation REgulation SearcH) is a gene-signature-based search engine from the Sergushichev lab (authors of `fgsea`, `gatom`). Given a short gene list, it ranks ~80,000 preprocessed GEO datasets by how much variance those genes jointly explain in each dataset -- a PCA-inspired score the authors call `pctVar`, with an optional `GESECA` p-value for specificity. The compendium covers **44,253 human** and **42,224 mouse** datasets across microarray and RNA-seq.

Under the hood CORESH is a **thin R wrapper around `fgsea::geseca`**. The core algorithm is ~20 lines (see `references/algorithm.md`); the complexity comes from (a) the preprocessed chunk format hosted on Synapse, (b) the strict Entrez-integer ID requirement, and (c) deciding which top hits are signal versus artifact.

The headline differentiator: **coregulation, not metadata.** Classical GEO search (or tools like ReGEO, Gemma, RummaGEO) match on keywords or study design. CORESH asks "where in public data do *these genes* move together?" -- a data-driven question that finds perturbation analogs and unexpected biological contexts that no metadata search would surface. The original use case (Mehrotra et al., *Nature* 2024) was identifying PGE2 as the mediator of a pyroptosis secretome by noticing that macrophage-PGE2 GEO datasets topped the ranking for the Pyro-1 signature.

**When to use this skill:**
- You have a gene signature (10-200 genes) from DE, cluster markers, or a pathway's leading edge, and want to find public biological contexts where those genes co-vary
- You want to derive novel, data-driven gene sets from top-ranking datasets and feed them into GSEA as a custom database
- You want to cross-validate a mouse finding against human public data (or vice versa)
- You want to generate hypotheses about perturbations or conditions that reproduce your signature, as leads for wet-lab validation

**When NOT to use this skill:**
- You want canonical-pathway interpretation (Hallmark, KEGG, Reactome, GO) -> use `bulk-rnaseq-gsea`
- You want metadata search (organism, disease, platform) -> use the GEO web UI or RummaGEO
- You want a single-dataset DE analysis -> use `pydeseq2` / DESeq2 directly
- You want topology-aware metabolic module finding -> use `gatom-metabolomic-predictions`

---

## Decision Tree

```
Got a gene signature?
|
+-- One-off query, quick exploration? --> Web UI (https://alserglab.wustl.edu/coresh/)
|                                         No install, paste symbols, ranked hits in ~30 s.
|                                         Good for: sanity checks, single signatures, demo.
|
+-- Batch (>=3 queries) or programmatic pipeline? --> R package + local chunks
|                                                     Install, download ~20 GB Synapse chunks,
|                                                     run coresh_batch.R across hsa/ or mmu/.
|                                                     Good for: systematic sweeps, reproducibility,
|                                                     extracting derived gene sets.
|
+-- Cross-species validation? --> Run query against BOTH hsa/ and mmu/ with orthologs
|                                (symbols_to_entrez.R handles both directions).
|
+-- Derived gene sets for GSEA? --> extract_gene_loadings.R on top-ranking GSEs,
                                    then route output to bulk-rnaseq-gsea (custom-db reference).
```

*Web UI is always the right first stop for a new signature, even if you plan to go programmatic -- it tells you in under a minute whether the signature produces interpretable hits at all.*

---

## Quick Start

The fastest path to a visible, verifiable result uses the **web UI**. No install, no Synapse, no ID conversion.

1. Open <https://alserglab.wustl.edu/coresh/>.
2. Paste ~50 symbols from HALLMARK_HYPOXIA (a strong positive control -- see `references/interpretation-protocol.md`).
3. Select species (human or mouse), submit.
4. Inspect top 20 hits. Expected: hypoxia-treatment studies, tumor-hypoxia studies, HIF-related perturbations appear repeatedly.

If the positive control looks right, you can trust the tool. Move to the R-package path for anything beyond one query.

---

## Progressive Depth

### Basic: single query in R, one chunk

Install package, load one chunk (~500 datasets), run `coreshMatch()` on a query. Mirrors the upstream vignette (`01_modules/.ref/coresh/vignettes/coresh-local.Rmd`).

```r
# one-time install
remotes::install_github("alserglab/coresh")

library(coresh)
library(qs2)
library(data.table)
library(fgsea)

# load one chunk for a smoke test (chunks live in preprocessed_chunks/{hsa,mmu}/)
chunk <- qs_read("preprocessed_chunks/hsa/001_full_objects.qs2")
length(chunk)          # ~500 dataset objects
str(chunk[[1]])        # $E1024 (scaled expression), $rownames (Entrez int), $gseId, $gplId, $totalVar

# query: HALLMARK_HYPOXIA Entrez IDs (integers, not symbols!)
query <- c(5230L, 226L, 3484L, 5209L, 5163L)  # PGK1, ALDOA, IGFBP1, PFKP, PDK1 (truncated)

# score the first dataset
coreshMatch <- function(obj, query, calculatePvalues = FALSE) {
  E <- obj$E1024 / 1024
  queryIdxs <- na.omit(match(query, obj$rownames))
  k <- length(queryIdxs)
  curProfile <- colSums(E[queryIdxs, , drop = FALSE])
  queryVar <- sum(curProfile^2)
  pval <- if (calculatePvalues) {
    fgsea:::gesecaCpp(E, queryVar, k, sampleSize = 21, seed = 1, eps = 1e-300)[[1]]$pval
  } else NA
  data.table(gse = obj$gseId, gpl = obj$gplId,
             pctVar = queryVar / k / obj$totalVar * 100,
             pval = pval, size = k)
}

coreshMatch(chunk[[1]], query, calculatePvalues = TRUE)
```

The function lives in the package as `coreshMatch()`. Inlined here because the vignette defines it inline and its ~20-line body is the algorithm.

### Intermediate: full-compendium sweep with Entrez conversion

Three things change at scale: (1) symbol->Entrez conversion becomes mandatory, (2) you parallelize across chunks with `BiocParallel`, (3) you choose variance-only (seconds) vs p-value (minutes).

```r
library(BiocParallel)

# --- 1. symbol -> Entrez (see scripts/symbols_to_entrez.R for helper) ---
library(org.Hs.eg.db)
symbols <- c("TFRC", "SLC11A2", "STEAP3", "STEAP4")
query   <- as.integer(na.omit(mapIds(org.Hs.eg.db, keys = symbols,
                                     keytype = "SYMBOL", column = "ENTREZID")))
# If any symbols fail to map, log them -- silent NAs are the #1 failure mode.

# --- 2. parallel sweep across all hsa chunks ---
bpparam <- MulticoreParam(8, progressbar = TRUE)
chunkPaths <- list.files("preprocessed_chunks/hsa/",
                         pattern = "full_objects.qs2", full.names = TRUE)

varRanking <- rbindlist(bplapply(chunkPaths, function(p) {
  ch <- qs_read(p)
  rbindlist(lapply(ch, coreshMatch, query = query, calculatePvalues = FALSE))
}, BPPARAM = bpparam))
varRanking <- varRanking[order(-pctVar)]
head(varRanking, 20)
```

Timing on an 8-core machine: variance-only ~10-20 s across the full hsa compendium; with p-values ~2-5 min. For batch execution with query logging, error handling, and tidy CSV output, use `scripts/coresh_batch.R`.

**Verify observable output:**

```r
stopifnot(nrow(varRanking) > 40000,          # full compendium hit
          all(varRanking$size <= length(query)),
          any(varRanking$pctVar > 5))        # at least some datasets above baseline noise
```

### Advanced: CoReSh-to-GSEA bridge

The most underappreciated feature. Each top-ranking GEO dataset in CORESH implicitly defines a **coregulation module**: the query genes plus every other gene that moves with them in that dataset's particular biological context. Extracting those co-moving genes produces a data-driven gene set you can feed into GSEA on your own dataset.

The mechanism: `E1024` is already centered (and PC-reduced -- see `references/algorithm.md`), so the query defines a direction in sample space via `profile = colSums(E[query, ])`. Projecting every gene onto that direction gives gene-level loadings. Top-|loading| genes are the "context-specific coregulation partners" of your query.

```r
# for a top-ranking object `obj` and integer Entrez query
profile <- colSums(obj$E1024[match(query, obj$rownames), , drop = FALSE]) / 1024
profile <- profile / sqrt(sum(profile^2))           # normalize so loadings are comparable
gene_loadings <- (obj$E1024 / 1024) %*% profile      # length = n_genes
top_genes <- obj$rownames[order(abs(gene_loadings), decreasing = TRUE)[1:50]]
# convert back to symbols for downstream GSEA (see references/gene-id-conversion.md)
```

Wrapper in `scripts/extract_gene_loadings.R`. Emitted GMT records are consumed by `bulk-rnaseq-gsea` (custom-db reference) via its custom-database ingestion pattern. Caveats:

- **Sign ambiguity:** `gene_loadings` is signed but the overall direction flip is not biologically meaningful (PCA-loaded) -- treat as unsigned by default.
- **Filter size:** drop derived sets with fewer than 15 genes (after ID conversion dropouts) or more than 500.
- **Deduplicate:** multiple top hits often produce overlapping sets; Jaccard-dedupe at threshold 0.8 before GSEA.

For cross-species validation, run the same signature in parallel on `hsa/` and `mmu/` chunks (with orthologs on the mouse side) and look for concordant top hits.

---

## Verification Checklist

Run these checks every time you use CORESH on a new signature:

- [ ] **Positive control:** Hallmark hypoxia signature ranks hypoxia/HIF/tumor-hypoxia studies in top 20 (`references/interpretation-protocol.md` has the expected GSE pattern).
- [ ] **Size sanity:** `size` column on top hits is at or near `length(query)` -- if most top hits have `size < 3`, the platform coverage is wrong or you passed symbols instead of Entrez IDs.
- [ ] **Interpretability:** at least 5 of the top 10 hits belong to biology adjacent to your signature (expected category) -- if fewer, the signature is too generic or too narrow.
- [ ] **Negative control:** a random housekeeping set (ACTB, GAPDH, B2M, ...) does NOT rank the same datasets -- if it does, your signal is confounded by ubiquitous coregulation, not specific biology.
- [ ] **Cross-query convergence (if running multiple related queries):** at least two related signatures share at least one top-10 dataset.

---

## Common Pitfalls

### Pitfall: Symbols instead of Entrez IDs

- **Symptom:** `size = 0` or near-zero across the compendium; `pctVar` is noise.
- **Cause:** Chunk `rownames` are **integer Entrez IDs**. `match("TFRC", c(7037L, 4891L, ...))` silently returns `NA` for every element -- R does not warn on type mismatch.
- **Fix:** always convert symbols to Entrez integers via `scripts/symbols_to_entrez.R` and assert `length(query) == length(na.omit(match(query, chunk[[1]]$rownames)))` on one chunk before the full sweep.

### Pitfall: Species mismatch

- **Symptom:** human signature ranks mouse `mmu/` chunks with tiny `pctVar` and `size` near 0.
- **Cause:** human Entrez IDs do not match mouse Entrez IDs.
- **Fix:** use orthologs (`homologene`, `babelgene`, or Ensembl BioMart). `symbols_to_entrez.R` accepts a `species` argument for both directions.

### Pitfall: Signature too small or too ubiquitous

- **Symptom:** `size < 3` across most datasets, OR high `pctVar` everywhere (no selectivity).
- **Cause:** signatures of 1-5 genes have unstable variance estimates; housekeeping-heavy signatures coregulate everywhere.
- **Fix:** aim for 10-200 genes after filtering housekeeping. For DE-derived signatures, take top 20-100 by `|stat|`, filter out ribosomal (RPS*/RPL*), mitochondrial-encoded (MT-*), and cell-cycle genes unless those are the biology you're studying.

### Pitfall: Cancer cell line dominance

- **Symptom:** iron / proliferation / metabolism signatures pull mostly CCLE-like tumor-line studies, obscuring primary-tissue signal.
- **Cause:** broad housekeeping-proximate programs drive variance in every cycling-cell dataset.
- **Fix:** re-rank with `fgsea::geseca` p-values (slower, more specific), or post-filter top hits by manual dataset-type categorization. CORESH does not provide automatic cell-line filtering.

### Pitfall: PC-sign ambiguity in the bridge

- **Symptom:** when comparing derived gene sets across datasets, "up" and "down" sides appear flipped between hits.
- **Cause:** the projection direction has an arbitrary overall sign (like PC1 loadings).
- **Fix:** treat derived sets as unsigned for GSEA (both directions tested); don't try to merge signed leading-edge lists across datasets.

### Pitfall: Synapse download gap

- **Symptom:** `coresh_batch.R` fails with "no chunks found" or `qs_read` errors.
- **Cause:** Claude cannot perform the Synapse download -- it requires the user's Synapse access token.
- **Fix:** the user runs the one-time download per `references/synapse-data-setup.md`. `scripts/validate_coresh_install.R` checks for chunk presence and fails clearly if absent.

---

## Complementary Skills

| When you need... | Use skill | Relationship |
|---|---|---|
| Canonical-pathway interpretation (MSigDB Hallmark/KEGG/Reactome/GO) of the same signature | `bulk-rnaseq-gsea` | Alternative (prior-knowledge question) |
| Run GSEA with CORESH-derived gene sets as a custom database | `bulk-rnaseq-gsea` | Next step (consumes `extract_gene_loadings.R` output) |
| Interactive HTML dashboards over the resulting master table | `bulk-rnaseq-pathway-explorer` | Downstream (visualizes GSEA output) |
| Topology-aware metabolic module discovery from DE | `gatom-metabolomic-predictions` | Alternative (module-level, not signature-level) |

---

## Resources

- **Paper:** Sukhov V et al. *CORESH: a gene signature-based search engine for public gene expression datasets.* Nucleic Acids Res. 2025 Jul 7;53(W1):W187-W192. doi: [10.1093/nar/gkaf372](https://doi.org/10.1093/nar/gkaf372). PMID: 40322919; PMCID: PMC12230675.
- **Method precedent:** Mehrotra P et al. *Oxylipins and metabolites from pyroptotic cells act as promoters of tissue repair.* Nature 631, 207-215 (2024). doi: [10.1038/s41586-024-07585-9](https://doi.org/10.1038/s41586-024-07585-9).
- **Web UI:** <https://alserglab.wustl.edu/coresh/>
- **Local vignette:** `01_modules/.ref/coresh/vignettes/coresh-local.Rmd` (authoritative code pattern)
- **Public vignette:** <https://rpubs.com/asergushichev/coresh-local>
- **Repo:** <https://github.com/alserglab/coresh>
- **Synapse data:** <https://www.synapse.org/coresh> (project `syn66227307`, ~20 GB)
- **Related skills in this toolkit:** `bulk-rnaseq-gsea` (downstream), `gatom-metabolomic-predictions` (same lab)

---

## Deeper reference (load-on-demand)

- `references/algorithm.md` -- the `pctVar` math, GESECA relationship, E1024 scaling, when p-values are worth the cost
- `references/synapse-data-setup.md` -- one-time Synapse token + `synapse get` workflow; chunk layout; must be run by the user, not by Claude
- `references/gene-id-conversion.md` -- symbol<->Entrez recipes (R and Python), species handling, ortholog mapping
- `references/query-design.md` -- how to build good queries: size, filtering, DE-derived vs curated, project-specific guidance for cDC1/iron biology (Q1-Q12)
- `references/coresh-to-gsea-bridge.md` -- full recipe for deriving gene sets from top hits and feeding them to `bulk-rnaseq-gsea` (custom-db reference)
- `references/interpretation-protocol.md` -- read-the-top-20 protocol, positive/negative controls, expected-hit categories, red flags
