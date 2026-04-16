---
name: bulk-rnaseq-gsea-msigdb
description: "MSigDB Gene Set Enrichment Analysis via clusterProfiler/fgsea -- runs pre-ranked GSEA on limma-voom DE results against MSigDB collections (H, C2, C3, C5). Use when executing GSEA with MSigDB gene sets, configuring msigdbr databases, or building the checkpoint-cached GSEA pipeline. For custom gene sets (MitoCarta, TransportDB) use bulk-rnaseq-gsea-custom-db. For visualization use bulk-rnaseq-gsea-visualization."
license: MIT
metadata:
  skill-author: SciAgent-toolkit
  last-reviewed: 2026-04-14
  version: 1.0.0
  upstream-docs: https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html
  category: analysis
  tier: standard
  tags:
    - gsea
    - msigdb
    - clusterProfiler
    - fgsea
    - bulk-rnaseq
    - pathway-analysis
    - limma-voom
    - gene-set-enrichment
    - hallmark
    - go-terms
  complementary-skills:
    - bulk-rnaseq-gsea-custom-db
    - bulk-rnaseq-gsea-master-tables
    - bulk-rnaseq-gsea-visualization
    - gatom-metabolomic-predictions
  contraindications:
    - "Do not use for non-MSigDB gene sets (MitoCarta, TransportDB, mitoXplorer, GMT/GMX files). Use bulk-rnaseq-gsea-custom-db instead."
    - "Do not use for visualization or export. Use bulk-rnaseq-gsea-visualization or bulk-rnaseq-gsea-master-tables instead."
---

# MSigDB GSEA Execution with clusterProfiler/fgsea

## Overview

This skill covers running Gene Set Enrichment Analysis (GSEA) against Molecular Signatures Database (MSigDB) collections using the RNAseq-toolkit wrapper around `clusterProfiler::GSEA()` with the `fgsea` backend. It handles gene ranking from limma-voom DE results, MSigDB gene set loading via `msigdbr` (with v7.5/v8+ API compatibility), GSEA execution with recommended parameters, checkpoint caching, and database batching strategy.

**When to use this skill:**
- Running GSEA against MSigDB collections (Hallmark, KEGG, Reactome, GO, WikiPathways, TF targets)
- Configuring which MSigDB databases to include in a pipeline
- Setting up the `run_gsea()` wrapper with correct parameters
- Handling `msigdbr` API version differences

**When NOT to use this skill:**
- Custom/external gene sets (MitoCarta, TransportDB, GMT files) -> use `bulk-rnaseq-gsea-custom-db`
- Plotting GSEA results (dotplots, barplots, heatmaps) -> use `bulk-rnaseq-gsea-visualization`
- Creating master CSV tables from gseaResult objects -> use `bulk-rnaseq-gsea-master-tables`
- Active metabolic module discovery -> use `gatom-metabolomic-predictions`

---

## Decision Tree

```
Need pathway enrichment from DE results?
|
+- Gene sets from MSigDB (H, C2, C3, C5)?
|  |
|  +- Single database run?      -> Quick Start (this skill)
|  +- Full multi-database pipeline? -> Progressive Depth: Intermediate (this skill)
|  +- Need to visualize results?    -> bulk-rnaseq-gsea-visualization
|  +- Need master CSV export?       -> bulk-rnaseq-gsea-master-tables
|
+- Gene sets from external source (GMT, TransportDB, MitoCarta)?
|  -> bulk-rnaseq-gsea-custom-db
|
+- Need metabolic network modules, not pathway enrichment?
   -> gatom-metabolomic-predictions
```

---

## Quick Start

Minimal single-database GSEA run using the toolkit wrapper:

```r
# Source toolkit and config
source("02_analysis/config/config.R")
source(file.path(DIR_TOOLKIT, "GSEA/GSEA_processing/run_gsea.R"))

# Load DE results (gene symbols as rownames, 't' column present)
de_results <- readRDS("03_results/checkpoints/1.1_de_results.rds")
de_table <- de_results[["IL2RA_KO - NTC"]]

# Run GSEA for Hallmark gene sets
gsea_hallmark <- run_gsea(
    DE_results   = de_table,
    rank_metric  = "t",
    species      = "Mus musculus",
    collection   = "H",
    subcollection = "",
    pvalue_cutoff = 1.0,
    nperm        = 100000,
    seed         = 123
)
```

**Verify it worked:**

```r
# Check the result object
stopifnot(is(gsea_hallmark, "gseaResult"))
stopifnot(nrow(gsea_hallmark@result) > 0)
message("Hallmark GSEA returned ", nrow(gsea_hallmark@result), " gene sets")

# Inspect top results
head(gsea_hallmark@result[order(gsea_hallmark@result$p.adjust), c("ID", "NES", "p.adjust", "setSize")])
```

---

## Progressive Depth

### Basic Usage

#### Gene Ranking from limma-voom

Genes are ranked by the moderated t-statistic (default), which incorporates both fold-change magnitude and estimation precision. The ranked vector must be named with gene symbols and sorted descending:

```r
# Inside run_gsea(), the ranking is built automatically:
gene_vector <- DE_results[[rank_metric]]   # Extract t-statistics
names(gene_vector) <- rownames(DE_results) # Gene symbols as names
ranked_genes <- sort(gene_vector, decreasing = TRUE)
```

The DE results table must come from `limma::topTable()` with `number = Inf` and `sort.by = "none"`, with gene symbols as rownames.

#### MSigDB Gene Set Loading

Gene sets are loaded via `msigdbr` and converted to a two-column TERM2GENE data frame:

```r
# run_gsea() builds TERM2GENE internally:
msigdb_df <- msigdbr(species = "Mus musculus", category = "H")
term2gene_df <- msigdb_df[, c("gs_name", "gene_symbol")]
```

The TERM2GENE format is a long-format table mapping gene set names to individual gene symbols:

| gs_name | gene_symbol |
|---------|-------------|
| HALLMARK_MYC_TARGETS_V1 | Gnl3 |
| HALLMARK_MYC_TARGETS_V1 | Ddx21 |
| HALLMARK_E2F_TARGETS | Brms1l |

#### Core GSEA Call

```r
set.seed(123)
GSEA_result <- clusterProfiler::GSEA(
    geneList      = ranked_genes,
    TERM2GENE     = term2gene_df,
    pvalueCutoff  = 1.0,          # Retain ALL pathways
    pAdjustMethod = "fdr",
    eps           = 0,            # Exact p-values via adaptive algorithm
    by            = "fgsea",
    nPermSimple   = 100000,
    verbose       = FALSE
)
```

**Critical parameters:**

| Parameter | Value | Rationale |
|-----------|-------|-----------|
| `eps` | `0` | Enables exact p-value calculation via fgsea's adaptive multi-level splitting. Without this, extremely significant pathways get truncated p-values. |
| `nPermSimple` | `100000` | Permutation count for precision. 100K is sufficient for most analyses. |
| `pvalueCutoff` | `1.0` | Retains ALL pathways in the result. Filtering happens downstream at visualization/export, avoiding expensive re-runs at different thresholds. |
| `seed` | `123` | Set via `set.seed()` before the GSEA call for reproducibility. |

### Intermediate Usage

#### Database Configuration in config.R

The project-level configuration defines which MSigDB collections to analyze. Each entry is a two-element vector: `c(category, subcategory)`:

```r
# From 02_analysis/config/config.R
MSIGDB_DATABASES <- list(
  H              = c("H", ""),                # Hallmark gene sets (50 sets)
  C2_KEGG        = c("C2", "CP:KEGG"),        # KEGG pathways (~180 sets)
  C2_REACTOME    = c("C2", "CP:REACTOME"),    # Reactome (~1600 sets)
  C2_WIKIPATHWAYS= c("C2", "CP:WIKIPATHWAYS"),# WikiPathways (~600 sets)
  C3_TF          = c("C3", "TFT:GTRD"),       # TF targets (~500 sets)
  C5_BP          = c("C5", "GO:BP"),           # GO Biological Process (~7500 sets)
  C5_MF          = c("C5", "GO:MF"),           # GO Molecular Function (~1700 sets)
  C5_CC          = c("C5", "GO:CC")            # GO Cellular Component (~1000 sets)
)
```

Parse entries in the pipeline loop:

```r
db_params <- MSIGDB_DATABASES[[db_name]]
category <- db_params[1]
subcategory <- if (nzchar(db_params[2])) db_params[2] else NULL
```

#### Database Batching Strategy

GSEA is split into batches cached as separate checkpoints. This allows fast re-loading of smaller collections without waiting for the expensive GO term analysis:

**Batch 1:** Hallmark + C2 + C3 (faster, fewer gene sets)
```r
gsea_h_c2 <- load_or_compute(
  checkpoint_file = "1.1_gsea_H_C2.rds",
  description = "GSEA results (H, C2, C3)",
  compute_fn = function() {
    db_names <- c("H", "C2_KEGG", "C2_REACTOME", "C2_WIKIPATHWAYS", "C3_TF")
    results <- list()
    for (db_name in db_names) {
      db_params <- MSIGDB_DATABASES[[db_name]]
      category <- db_params[1]
      subcategory <- if (nzchar(db_params[2])) db_params[2] else NULL
      results[[db_name]] <- run_gsea(
        DE_results   = de_table,
        rank_metric  = RANK_METRIC,
        species      = SPECIES,
        collection   = category,
        subcollection = subcategory %||% "",
        pvalue_cutoff = GSEA_PVALUE_CUTOFF,
        nperm        = GSEA_NPERM,
        seed         = GSEA_SEED
      )
    }
    return(results)
  }
)
```

**Batch 2:** GO terms / C5 (thousands of gene sets, most expensive)
```r
gsea_c5 <- load_or_compute(
  checkpoint_file = "1.1_gsea_C5.rds",
  description = "GSEA results (C5 GO terms)",
  compute_fn = function() {
    db_names <- c("C5_BP", "C5_MF", "C5_CC")
    results <- list()
    for (db_name in db_names) {
      # ... same pattern as above ...
    }
    return(results)
  }
)
```

**Combination checkpoint:**
```r
all_gsea <- load_or_compute(
  checkpoint_file = "1.1_all_gsea_results.rds",
  description = "All GSEA results combined",
  compute_fn = function() {
    combined <- c(gsea_h_c2, gsea_c5)
    return(combined)
  }
)
```

The final `all_gsea` object is a named list:
```
all_gsea$H               -> gseaResult (Hallmark)
all_gsea$C2_KEGG         -> gseaResult (KEGG)
all_gsea$C2_REACTOME     -> gseaResult (Reactome)
all_gsea$C2_WIKIPATHWAYS -> gseaResult (WikiPathways)
all_gsea$C3_TF           -> gseaResult (TF targets)
all_gsea$C5_BP           -> gseaResult (GO BP)
all_gsea$C5_MF           -> gseaResult (GO MF)
all_gsea$C5_CC           -> gseaResult (GO CC)
```

#### Checkpoint Caching Pattern

Every expensive GSEA computation is wrapped in `load_or_compute()`:

```r
load_or_compute <- function(checkpoint_file, compute_fn,
                            force_recompute = FALSE, description = "Result") {
  if (!grepl("^/", checkpoint_file) && !grepl("^[A-Z]:", checkpoint_file)) {
    checkpoint_path <- file.path(DIR_CHECKPOINTS, checkpoint_file)
  } else {
    checkpoint_path <- checkpoint_file
  }
  if (file.exists(checkpoint_path) && !force_recompute) {
    message(sprintf("[CACHE] Loading %s from: %s", description, checkpoint_file))
    return(readRDS(checkpoint_path))
  }
  message(sprintf("[COMPUTE] Computing %s...", description))
  start_time <- Sys.time()
  result <- compute_fn()
  elapsed <- round(difftime(Sys.time(), start_time, units = "mins"), 2)
  message(sprintf("[SAVE] Saving %s to: %s (took %.2f min)", description, checkpoint_file, elapsed))
  saveRDS(result, checkpoint_path)
  return(result)
}
```

Force recomputation when parameters change: `load_or_compute(..., force_recompute = TRUE)`.

### Advanced Usage

#### msigdbr v7.5 vs v8+ API Compatibility

The `msigdbr` package underwent a breaking API change. The toolkit's `run_gsea()` detects the installed version at runtime:

```r
# From run_gsea.R -- runtime API detection
msigdbr_params <- names(formals(msigdbr::msigdbr))
use_new_api <- "collection" %in% msigdbr_params

if (use_new_api) {
  # v8+: uses collection/subcollection/db_species
  msigdb_df <- msigdbr(
    db_species    = db_species,
    species       = species,
    collection    = collection,
    subcollection = subcollection
  )
} else {
  # v7.5.x: uses category/subcategory
  msigdb_df <- msigdbr(
    species     = species,
    category    = collection,
    subcategory = subcollection
  )
}
```

Always use `collection`/`subcollection` as parameter names when calling `run_gsea()` -- the function translates internally.

#### The gseaResult Object Structure

`clusterProfiler::GSEA()` returns a `gseaResult` S4 object. Key slot is `@result`:

| Column | Description |
|--------|-------------|
| `ID` | Gene set identifier (e.g., `HALLMARK_MYC_TARGETS_V1`) |
| `setSize` | Genes in set found in ranked list |
| `NES` | Normalized Enrichment Score |
| `pvalue` | Nominal p-value |
| `p.adjust` | BH-adjusted p-value |
| `rank` | Position where running sum peaks |
| `core_enrichment` | Slash-separated leading edge genes |

#### Script Organization

| File | Role |
|------|------|
| `02_analysis/config/config.R` | Defines `MSIGDB_DATABASES`, parameters, `load_or_compute()` |
| `02_analysis/config/pipeline.yaml` | Shared YAML config (species, colors, schemas) |
| `01_scripts/RNAseq-toolkit/scripts/GSEA/GSEA_processing/run_gsea.R` | Core single-database wrapper |
| `01_scripts/RNAseq-toolkit/scripts/GSEA/GSEA_processing/run_gsea_analysis.R` | Multi-database pipeline with auto-plotting |
| `02_analysis/1.1.core_pipeline.R` | Executes GSEA across all configured databases |
| `02_analysis/1.5.create_master_tables.R` | Aggregates into master CSV (see bulk-rnaseq-gsea-master-tables) |

#### Alternative: run_gsea_analysis() for Standalone Use

For standalone analysis outside a project pipeline, the toolkit provides `run_gsea_analysis()` which auto-sources dependencies and generates plots:

```r
source(file.path(toolkit_dir, "GSEA/GSEA_processing/run_gsea_analysis.R"))

results <- run_gsea_analysis(
    de_table      = de_results,
    analysis_name = "MyAnalysis",
    rank_metric   = "t",
    species       = "Mus musculus",
    n_pathways    = 30,
    padj_cutoff   = 0.05,
    nperm         = 100000,
    pvalue_cutoff = 1,
    save_plots    = TRUE,
    output_dir    = "./GSEA_Plots"
)
```

---

## Verification Checklist

After running GSEA, confirm:

- [ ] **Result object type:** Each element of the results list `is(result, "gseaResult")` returns TRUE
- [ ] **Non-empty results:** `nrow(result@result) > 0` for every database (zero rows means no gene sets were tested)
- [ ] **Checkpoint files exist:** `file.exists("03_results/checkpoints/1.1_gsea_H_C2.rds")` and `1.1_gsea_C5.rds`
- [ ] **All databases present:** `names(all_gsea)` contains all expected keys (`H`, `C2_KEGG`, `C2_REACTOME`, etc.)
- [ ] **Biological plausibility:** Hallmark results include known biology (e.g., E2F targets, MYC targets should appear if proliferation is affected)
- [ ] **P-value precision:** Very significant pathways have p-values well below 1e-10 (confirms `eps=0` is working; if all p-values are truncated at ~1e-4, eps was not set to 0)

---

## Common Pitfalls

### Pitfall: Empty GSEA results (zero gene sets tested)

- **Symptom:** `nrow(gsea_result@result) == 0` or warning about no enrichment.
- **Cause:** Gene symbol mismatch between DE results rownames and msigdbr gene symbols. Most common when DE table uses Ensembl IDs instead of gene symbols.
- **Fix:** Ensure DE results have gene symbols as rownames. Check with `head(rownames(de_table))` -- should be symbols like `Il2ra`, `Cd3e`, not `ENSMUSG00000...`.

### Pitfall: Truncated p-values at ~1e-4

- **Symptom:** The most significant pathways all show the same minimum p-value (e.g., 9.99e-05).
- **Cause:** `eps` parameter not set to 0. The default fgsea eps truncates small p-values.
- **Fix:** Ensure `eps = 0` in the `clusterProfiler::GSEA()` call. This enables the adaptive multi-level splitting algorithm for exact p-value estimation.

### Pitfall: msigdbr API error on category/collection

- **Symptom:** Error like `unused argument (category = "H")` or `unused argument (collection = "H")`.
- **Cause:** Calling `msigdbr()` directly with the wrong parameter name for the installed version. v7.5 uses `category`/`subcategory`; v8+ uses `collection`/`subcollection`.
- **Fix:** Use the toolkit's `run_gsea()` wrapper, which auto-detects the API version. If calling `msigdbr` directly, check with `names(formals(msigdbr::msigdbr))`.

### Pitfall: Missing subcollection returns entire category

- **Symptom:** Running C2 without specifying subcollection returns all of C2 (thousands of gene sets from KEGG + Reactome + CGP + etc.) instead of just one subcollection.
- **Cause:** Subcollection was `""` or `NULL` when a specific subcollection was intended.
- **Fix:** Verify the config entry. For Hallmark (`H`), empty subcollection is correct (the entire category is one collection). For C2, C3, C5, always specify the subcollection (e.g., `"CP:KEGG"`, `"GO:BP"`).

### Pitfall: GSEA re-runs despite existing checkpoint

- **Symptom:** GSEA computation runs from scratch even though the checkpoint RDS file exists.
- **Cause:** Checkpoint path resolution failure -- typically because `DIR_CHECKPOINTS` is not set or the relative path does not match.
- **Fix:** Verify `config.R` is sourced before `load_or_compute()`. Check that `file.exists(file.path(DIR_CHECKPOINTS, "1.1_gsea_H_C2.rds"))` returns TRUE.

---

## Complementary Skills

| When you need... | Use skill | Relationship |
|---|---|---|
| GSEA with custom gene sets (MitoCarta, TransportDB, GMT) | `bulk-rnaseq-gsea-custom-db` | Alternative |
| Normalize gseaResult to CSV master tables | `bulk-rnaseq-gsea-master-tables` | Next step |
| Dotplots, barplots, heatmaps of GSEA results | `bulk-rnaseq-gsea-visualization` | Next step |
| Active metabolic module discovery from DE | `gatom-metabolomic-predictions` | Extension |

---

## Resources

- **clusterProfiler docs:** https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html
- **fgsea paper:** Korotkevich et al. (2021) bioRxiv, doi:10.1101/060012
- **msigdbr CRAN:** https://cran.r-project.org/package=msigdbr
- **MSigDB collections:** https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp
- **RNAseq-toolkit GSEA docs:** `01_scripts/RNAseq-toolkit/docs/GSEA-workflow/02-msigdb-gsea-pipeline.md`
