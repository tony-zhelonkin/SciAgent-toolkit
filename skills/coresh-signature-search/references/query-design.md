# Query design -- what makes a good CORESH signature

The `pctVar` score has a characteristic regime where it's informative, flanked by two failure modes:

| Signature size | Behavior |
|---|---|
| 1-5 genes | Very noisy; `queryVar` is dominated by a single highly-variable gene in each dataset; almost any cycling-cell or cell-line study ranks at the top |
| 10-200 genes | **Optimal** -- enough genes to average out per-gene noise, few enough that coregulation is specific |
| 200-1000 genes | Dilution; coregulation gets averaged over unrelated sub-programs; `pctVar` becomes a general "how variable is this dataset?" measure |
| >1000 genes | Essentially a proxy for `totalVar`; use GSEA on leading edges instead |

## Recipes by signature source

### From differential expression

```r
# after DESeq2::results() or edgeR::topTags()
top_up   <- head(res[order(-res$stat), ], 50)          # top 50 up-regulated
top_down <- head(res[order(res$stat), ], 50)           # top 50 down-regulated
signature <- rownames(top_up)                           # use up only (default)
# or combine both sides if the biology is bidirectional
```

Prefer `stat` (Wald/moderated t) over raw `log2FoldChange` + `padj` because `stat` is continuous and well-calibrated at the top.

### From single-cell cluster markers

```r
library(scran)  # or sc.tl.rank_genes_groups in Python
# pull top-20 markers per cluster with padj < 0.05 and logfc > 1
markers <- findMarkers(sce, groups = sce$cluster, direction = "up")
cluster_7_sig <- head(rownames(markers[[7]]), 20)
```

Cluster markers are typically 20-50 genes and dominate their cluster's biology -- good for CORESH.

### From curated pathways

Leading-edge genes from a prior GSEA are ideal. Alternatively, the MSigDB Hallmark sets (50-200 genes each) work directly.

```r
library(msigdbr)
hallmarks <- msigdbr(species = "Homo sapiens", collection = "H")
hypoxia <- hallmarks |>
  dplyr::filter(gs_name == "HALLMARK_HYPOXIA") |>
  dplyr::pull(ncbi_gene) |>
  as.integer()
# ready to use -- already Entrez IDs
```

## Filtering to remove confounders

Remove genes that coregulate everywhere regardless of biology:

| Category | Pattern | Why |
|---|---|---|
| Ribosomal | `^RP[SL]\d+` | Drive variance in every proliferating-cell study |
| Mitochondrial-encoded | `^MT-` | Correlated with QC/dissociation state |
| Hemoglobin | `^HB[ABDEG]\d*` | Dominate any blood-contaminated sample |
| Pseudogenes | `^[A-Z0-9]+P\d+$` | Noise |
| Cell cycle (unless that's the question) | MKI67, TOP2A, CCNB1/2, CDK1, CDC20, BIRC5, TYMS | Drive variance in every cycling-cell study |
| X/Y-linked (for sex-confounded comparisons) | XIST, RPS4Y1, DDX3Y, KDM5D | Sex effect, not biology |

Decision: remove them **before** building the signature, unless the removed category IS the biology (e.g. "I want to rank cell-cycle datasets" -> keep cell-cycle genes).

## Project-specific guidance: cDC1 / iron biology signatures

For the `docs/coresh-plan/03_coresh_integration.md` query set, all signatures come from the `02_analysis/modules/markers.py` IRON_PROGRAMS and DC_STATE_PROGRAMS dictionaries plus data-derived DE results. Guidance per query type:

| Query type | Size typical | Filter notes |
|---|---|---|
| Iron program (uptake, storage, export, regulatory) | 3-12 genes | Small; expect high variance in `pctVar` across datasets. Prefer running at least 4 programs in parallel and looking for CROSS-QUERY CONVERGENCE on shared top hits rather than trusting any single ranking. |
| DC state program (immunogenic, tolerogenic, antigen-presentation, cross-presentation) | 10-40 genes | Well-sized; filter HLA-* if cross-species (mouse compendium) |
| Cluster marker signature | top 20 by `stat` from scanpy `rank_genes_groups` | Remove mitochondrial and ribosomal automatically via the helper in `scripts/symbols_to_entrez.R` |
| TFRC-high vs TFRC-low signature | top 20-50 by `stat` from pseudobulk DE | Check positive control first -- expect iron-uptake/heme-biosynthesis studies near the top |

Run queries in matched pairs (healthy vs tumor, hi-TFRC vs lo-TFRC, cluster 7 vs cluster 6) and examine which datasets rank high for both members of the pair -- that's where the biology is shared.

## Asserting on the query before the sweep

```r
stopifnot(is.integer(query),
          length(query) >= 5,
          length(query) <= 500)
```

The `coresh_batch.R` wrapper bakes these assertions in. If you hand-roll the sweep, always assert.
