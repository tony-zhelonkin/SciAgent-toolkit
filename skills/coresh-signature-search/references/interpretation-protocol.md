# Interpretation protocol -- how to read CORESH output

CORESH gives you a ranked list of GEO datasets. Interpreting it well is the difference between a useful hypothesis and confused noise.

## Before you look at the results

Decide up front what the **expected** top-hit categories are. If you cannot write them down before running the query, your signature is too vague for CORESH -- go back and sharpen it.

| Query | Expected top-hit categories (positive control) |
|---|---|
| HALLMARK_HYPOXIA | HIF-1a activation, hypoxia chambers, tumor hypoxia, VHL knockdown |
| HALLMARK_INTERFERON_GAMMA_RESPONSE | IFN-gamma stimulation, viral infection, M1 macrophage polarization |
| Iron uptake (TFRC, STEAP3, SLC11A2, ...) | Erythroid differentiation, macrophage iron recycling, liver iron overload, iron chelation |
| DC immunogenic markers | LPS/poly(I:C) stimulation, TLR agonists, DC maturation time-courses |
| DC tolerogenic markers | IL-10/TGF-b treatment, tumor-associated DC profiles, Treg-inducing contexts |
| Cross-presentation | Antigen cross-priming, TAP-related studies, Batf3 DCs |
| Hypothesis generating? Open question? | Run the positive controls (iron uptake, iron storage) in parallel and see if they converge on the same top hits |

If you run the query and the top 20 contain **zero** of your expected categories, something is wrong: check Entrez conversion, species matching, and size column (see `references/gene-id-conversion.md`).

## The read-the-top-20 protocol

For every query, walk through the top 20 hits and categorize each:

1. **Retrieve metadata** for the GSE (not in the CORESH object; fetch from GEO):

   ```r
   library(GEOquery)
   meta <- getGEO(gse_id, GSEMatrix = FALSE)    # slower, full record
   # or, faster, just the summary line:
   meta_url <- sprintf(
     "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=%s&targ=self&form=text&view=brief",
     gse_id)
   head(readLines(meta_url), 20)
   ```

   Alternative: use the CORESH web UI to look up hits interactively once you have the GSE IDs.

2. **Categorize each hit** into one of:

   | Category | Example | Meaning |
   |---|---|---|
   | Expected-biology | "macrophage M1 polarization" for an iron-uptake query | Confirms signal is real |
   | Adjacent-biology | "granulocyte differentiation" for an iron query | Plausibly related; worth follow-up |
   | Novel | "ER stress in beta cells" for a cDC1 query | Hypothesis lead -- verify before trusting |
   | Cell-line / CCLE | "CCLE iron chelation across 200 lines" | Usually confounded with proliferation; low weight |
   | Noise | `size < 3` or random-looking metadata | Discard |

3. **Score the query overall:**

   - `>= 5/10 top hits` in Expected + Adjacent -> signature is robust for CORESH; proceed to derived-set extraction.
   - `2-4/10 top hits` interpretable -> signature is noisy; consider refining before the derived-set step.
   - `<= 1/10 top hits` interpretable -> CORESH is not giving you signal; try a different signature construction.

## Cross-query convergence

The strongest signal comes from datasets that rank high for **multiple related queries**. If Iron Uptake AND Iron Storage AND Heme Biosynthesis all have the same GSE in their top 20, that dataset captures iron biology broadly -- a high-confidence lead.

Protocol:

```r
# for a set of related queries run via coresh_batch.R
top10_by_query <- varRanking[, .SD[1:10], by = query_name]
convergence <- top10_by_query[, .(n_queries = uniqueN(query_name)),
                              by = gse]
convergence[order(-n_queries)][n_queries >= 2]
```

GSEs with `n_queries >= 2` are convergence hits; examine them first before pulling derived gene sets.

## Negative controls

Run these any time you're unsure whether CORESH ranking is specific to your biology or just finds "variable datasets":

| Control | Genes | Expected behavior |
|---|---|---|
| Housekeeping | ACTB, GAPDH, B2M, HPRT1, TBP, PPIA | Should rank uniformly across compendium; no clustering of top hits around a biology |
| Random | sample 50 genes with seed | Same -- no coherent top hits |
| X/Y-linked (sex-confounded) | XIST, RPS4Y1 | Top hits will cluster around sex-variable datasets (tumors, gametogenesis) -- useful sanity check that the machinery works |

If your real query's top hits look "like" the housekeeping control (no coherent biology), your signature is driven by housekeeping-adjacent genes. Filter and retry (`references/query-design.md` has the filter patterns).

## Red flags

| Flag | What it means |
|---|---|
| `size < 3` in most top hits | Platform coverage issue or ID conversion failure |
| Top hits all from one lab / one year / one platform | Batch-effect-driven -- deprioritize |
| Top hits are all tumor cell lines | Proliferation-adjacent signal, not your biology |
| Top hits span wildly unrelated biology | Signature is too broad or noisy |
| `pval` disagrees strongly with `pctVar` ranking | Use `pval` -- it's more specific. Large mismatch usually means the dataset has very high `totalVar` but low real specificity |
| One GSE dominates multiple queries with `pctVar > 50%` | Often a pathological dataset (single-gene overexpression, extreme outlier samples); check manually |

## Writing up CORESH results

Minimum provenance to include when a CORESH finding makes it into a manuscript:

- Exact query gene set (Entrez IDs, not just symbols)
- Species (hsa or mmu)
- Chunk compendium version / Synapse snapshot date
- Ranking metric used (`pctVar` or `pval`) and threshold
- Top hits by accession + their biological category
- Whether derived gene sets were used downstream, and how they were generated (link to `coresh-to-gsea-bridge.md` or cite the script version)

The published paper (Sukhov 2025) uses GESECA p-values for its main figures and reports accession + `padj` for each top hit -- a reasonable convention to follow.
