---
name: decoupler-tf-activity
description: "DecoupleR TF activity inference -- infers transcription factor activities from bulk RNA-seq DE results using ULM with CollecTRI regulatory network. Use when computing TF activity scores from limma/edgeR DE output, identifying differentially active transcription factors, or building master TF activity tables. For pathway-level activity (PROGENy) use progeny-pathway-activity. For GSEA-based enrichment use bulk-rnaseq-gsea-msigdb."
license: MIT
metadata:
  skill-author: SciAgent-toolkit
  last-reviewed: 2026-04-14
  version: 1.0.0
  upstream-docs: https://saezlab.github.io/decoupleR/
  category: analysis
  tier: standard
  tags:
    - decoupler
    - transcription-factor
    - ulm
    - collectri
    - tf-activity
    - bulk-rnaseq
    - omnipath
    - regulatory-network
    - differential-activity
    - limma
  complementary-skills:
    - progeny-pathway-activity
    - bulk-rnaseq-gsea-msigdb
    - bulk-rnaseq-gsea-master-tables
    - decoupler-activity-visualization
  contraindications:
    - "Do not use for pathway-level activity inference. Use progeny-pathway-activity instead."
    - "Do not use for gene set enrichment (GSEA). Use bulk-rnaseq-gsea-msigdb instead."
    - "Do not use for visualization. Use decoupler-activity-visualization instead."
---

# DecoupleR Transcription Factor Activity Inference

## Overview

DecoupleR infers transcription factor (TF) activities from gene expression data by leveraging prior-knowledge regulatory networks. The Univariate Linear Model (ULM) method regresses each TF's known target gene weights against the observed expression changes (t-statistics from DE analysis), producing a per-TF activity score and p-value. This approach is fundamentally different from GSEA -- ULM uses the continuous mode-of-regulation (MoR) weights from CollecTRI (activator=+1, repressor=-1, with intermediate values), while GSEA treats gene sets as unweighted membership lists.

**When to use this skill:**
- Inferring which transcription factors are differentially active between conditions
- Building a master TF activity table from limma/edgeR DE results
- Working with the CollecTRI regulatory network via OmnipathR
- Needing continuous activity scores (not just enrichment p-values) for TFs

**When NOT to use this skill:**
- Pathway-level activity inference (PROGENy/14 pathways) -> use `progeny-pathway-activity`
- Gene set enrichment analysis against MSigDB -> use `bulk-rnaseq-gsea-msigdb`
- Plotting TF activity heatmaps, barplots, or volcano plots -> use `decoupler-activity-visualization`

---

## Decision Tree

```
Need to identify active transcription factors from DE results?
|
+- Want TF-level activity scores with regulatory weights?
|  -> This skill (DecoupleR ULM + CollecTRI)
|
+- Want pathway-level activity scores (14 cancer pathways)?
|  -> progeny-pathway-activity (DecoupleR + PROGENy)
|
+- Want enrichment of curated TF target gene sets?
|  -> bulk-rnaseq-gsea-msigdb (C3 TFT:GTRD collection)
|
+- Want to visualize TF activity results?
   -> decoupler-activity-visualization
```

---

## Quick Start

Minimal TF activity inference from a limma DE results table:

```r
library(decoupleR)
library(OmnipathR)
library(dplyr)

# Load DE results (must have gene symbols as rownames and a 't' column)
de_results <- readRDS("03_results/checkpoints/1.1_de_results.rds")
deg_table <- de_results[["IL2RAKO_vs_NTC"]]

# Prepare input matrix: genes x contrasts, using t-statistics
mat <- deg_table %>% dplyr::select(t) %>% as.matrix()
mat[is.na(mat)] <- 0

# Fetch CollecTRI network for mouse
net <- decoupleR::get_collectri(organism = "mouse", split_complexes = FALSE)

# Run ULM
tf_acts <- decoupleR::run_ulm(
  mat = mat, net = net,
  .source = "source", .target = "target", .mor = "mor",
  minsize = 5
)

# Filter to the contrast column (named 't' because that was the column name in mat)
tf_acts <- tf_acts %>% filter(condition == "t")
```

**Verify it worked:**

```r
stopifnot(nrow(tf_acts) > 0)
stopifnot(all(c("source", "score", "p_value") %in% colnames(tf_acts)))
message("Inferred activities for ", nrow(tf_acts), " TFs")
# Expect 100-1000 TFs depending on species and minsize
stopifnot(nrow(tf_acts) >= 50 && nrow(tf_acts) <= 2000)
# Activity scores should be roughly centered near 0
stopifnot(abs(mean(tf_acts$score)) < 5)
```

---

## Progressive Depth

### Basic Usage

#### Input Matrix Preparation

ULM expects a numeric matrix with genes as rows and conditions/samples as columns. For DE-based TF activity, the t-statistic is the ideal input because it captures both direction and precision of expression change:

```r
# From a limma topTable (gene symbols as rownames, 't' column present)
mat <- deg_table %>% dplyr::select(t) %>% as.matrix()
mat[is.na(mat)] <- 0

# Dimensions: genes x 1 (single contrast)
# rownames: gene symbols (must match network target names)
# colnames: "t" (becomes the 'condition' in output)
```

Why t-statistics and not logFC? The t-statistic down-weights genes with high variance, giving more reliable TF activity estimates. LogFC can be used but produces noisier results.

#### Loading CollecTRI

CollecTRI is a comprehensive collection of TF-target interactions curated from literature and ChIP-seq data, accessed via OmnipathR. The network is a three-column data frame:

| Column | Description |
|--------|-------------|
| `source` | Transcription factor name (gene symbol) |
| `target` | Target gene name (gene symbol) |
| `mor` | Mode of regulation: +1 (activator), -1 (repressor), or intermediate values |

```r
# Fetch and cache the network
net <- decoupleR::get_collectri(organism = "mouse", split_complexes = FALSE)
# organism: "mouse" or "human" (not "Mus musculus")
# split_complexes: FALSE keeps TF complexes as single entries
```

Always cache the network to avoid repeated downloads:

```r
collectri_file <- file.path(DIR_CHECKPOINTS, "collectri_mouse.rds")
if (file.exists(collectri_file)) {
  net <- readRDS(collectri_file)
} else {
  net <- decoupleR::get_collectri(organism = "mouse", split_complexes = FALSE)
  saveRDS(net, collectri_file)
}
```

#### Running ULM

```r
tf_acts <- decoupleR::run_ulm(
  mat     = mat,
  net     = net,
  .source = "source",   # Column in net with TF names
  .target = "target",   # Column in net with target gene names
  .mor    = "mor",       # Column in net with mode of regulation weights
  minsize = 5            # Minimum regulon size (TFs with fewer targets are dropped)
)
```

The `minsize = 5` parameter filters out TFs with fewer than 5 targets present in the input matrix. This prevents unreliable activity estimates from tiny regulons. Increasing to 10 gives more robust but fewer TFs.

Output columns:

| Column | Description |
|--------|-------------|
| `statistic` | Always "ulm" |
| `source` | TF name |
| `condition` | Column name from input matrix ("t") |
| `score` | Activity score (t-value from the linear model) |
| `p_value` | P-value for the activity score |

### Intermediate Usage

#### Master Table Schema Mapping

To integrate TF activities into the pathway explorer and master table ecosystem, map ULM output columns to the standard GSEA-like schema:

```r
# Get all regulon targets present in our data for core_enrichment
targets_list <- net %>%
  filter(source %in% tf_acts$source) %>%
  filter(target %in% rownames(mat)) %>%
  group_by(source) %>%
  summarise(genes = paste(unique(target), collapse = "/"))

# Build master table
final_df <- tf_acts %>%
  left_join(targets_list, by = "source") %>%
  mutate(
    pathway_id       = source,
    pathway_name     = source,
    database         = "CollecTRI",
    contrast         = "IL2RAKO_vs_NTC",
    nes              = score,         # ULM score acts as NES equivalent
    padj             = p_value,       # See pitfall: no multiple-testing correction
    pvalue           = p_value,
    set_size         = NA_integer_,   # Can be computed from regulon size
    leading_edge_size = NA_integer_,
    gene_ratio       = NA_character_,
    core_enrichment  = genes,         # Slash-separated targets for Jaccard overlap
    direction        = ifelse(score > 0, "Up", "Down"),
    neg_log_padj     = -log10(p_value)
  ) %>%
  dplyr::select(
    pathway_id, pathway_name, nes, pvalue, padj, set_size,
    leading_edge_size, gene_ratio, core_enrichment, database,
    contrast, neg_log_padj, direction
  )

write_csv(final_df, file.path(DIR_TABLES, "master_tf_activities.csv"))
```

The `core_enrichment` field uses the full regulon (all targets in the dataset) rather than a leading-edge subset. This is intentional -- for Jaccard similarity calculations in the pathway explorer, broader gene lists produce more meaningful overlap scores.

#### Checkpoint Caching Pattern

Wrap the full TF analysis in `load_or_compute()`:

```r
source("02_analysis/config/config.R")

tf_results <- load_or_compute(
  checkpoint_file = "1.7_tf_results.rds",
  description     = "DecoupleR TF activity inference (ULM + CollecTRI)",
  compute_fn      = function() {
    net <- get_or_cache_collectri()
    mat <- prepare_tstat_matrix(de_results[[contrast_name]])
    tf_acts <- decoupleR::run_ulm(
      mat = mat, net = net,
      .source = "source", .target = "target", .mor = "mor",
      minsize = 5
    )
    tf_acts %>% filter(condition == "t")
  }
)
```

### Advanced Usage

#### Multiple Contrasts

When you have multiple contrasts, build a multi-column matrix:

```r
# Build matrix with one column per contrast
contrasts <- names(de_results)
mat <- sapply(contrasts, function(cname) {
  tvals <- de_results[[cname]]$t
  names(tvals) <- rownames(de_results[[cname]])
  tvals
})
mat[is.na(mat)] <- 0

# run_ulm handles multi-column matrices automatically
tf_acts <- decoupleR::run_ulm(
  mat = mat, net = net,
  .source = "source", .target = "target", .mor = "mor",
  minsize = 5
)

# The 'condition' column now contains the contrast names
tf_acts %>% filter(condition == "IL2RAKO_vs_NTC")
```

#### Alternative Networks: DoRothEA

DoRothEA provides confidence-level-filtered TF regulons (A-E, where A is highest confidence). It is an alternative to CollecTRI when you need stricter evidence filtering:

```r
# DoRothEA with confidence levels A and B only
net_dorothea <- decoupleR::get_dorothea(
  organism = "mouse",
  levels = c("A", "B")
)

tf_acts_dorothea <- decoupleR::run_ulm(
  mat = mat, net = net_dorothea,
  .source = "source", .target = "target", .mor = "mor",
  minsize = 5
)
```

CollecTRI is generally preferred over DoRothEA for bulk RNA-seq because it has broader coverage and uses continuous MoR weights. DoRothEA is more common in single-cell workflows where stricter filtering compensates for noisier data.

#### Species Handling

The `organism` parameter in `get_collectri()` and `get_dorothea()` accepts short-form species names, not binomial nomenclature:

| Species | Argument |
|---------|----------|
| Mouse | `organism = "mouse"` |
| Human | `organism = "human"` |
| Rat | `organism = "rat"` |

Gene symbols are automatically mapped to the correct species. If your DE results use human symbols for a mouse experiment (e.g., after ortholog mapping), use `organism = "human"`.

---

## Verification Checklist

After running TF activity inference, confirm:

- [ ] **TF count is reasonable:** Between 100-1000 TFs for mouse/human with `minsize = 5`. Fewer than 50 suggests gene symbol mismatch with the network. More than 2000 is unexpected.
- [ ] **Activity distribution centered near 0:** `mean(tf_acts$score)` should be close to 0 (within +/-5). A strong shift suggests input bias.
- [ ] **Known biology validation:** Check that TFs expected to change based on experimental design show activity. For example, in an IL2RA-KO experiment, STAT5 targets should be affected.
- [ ] **P-value distribution:** `hist(tf_acts$p_value)` should show a roughly uniform distribution with a spike near 0 (for truly active TFs), not all values near 1 or all near 0.
- [ ] **Checkpoint file exists:** `file.exists("03_results/checkpoints/1.7_tf_results.rds")` returns TRUE.
- [ ] **Master table written:** `file.exists("03_results/tables/master_tf_activities.csv")` returns TRUE and contains expected columns.

---

## Common Pitfalls

### Pitfall: Missing DE checkpoint

- **Symptom:** `Error: DE results checkpoint not found` at script startup.
- **Cause:** The DE analysis script (`1.1.core_pipeline.R`) has not been run, or `CHECKPOINT_DE_RESULTS` in `config.R` points to a non-existent file.
- **Fix:** Run `Rscript 02_analysis/1.1.core_pipeline.R` first, or verify that `file.path(DIR_CHECKPOINTS, CHECKPOINT_DE_RESULTS)` exists.

### Pitfall: Gene symbol mismatch with network

- **Symptom:** Very few TFs inferred (under 50), or `run_ulm()` returns empty results.
- **Cause:** The rownames of the input matrix do not match the `target` column in the CollecTRI network. Common when DE results use Ensembl IDs instead of gene symbols, or when human symbols are used with `organism = "mouse"`.
- **Fix:** Check overlap: `sum(rownames(mat) %in% net$target)` should be several thousand. If low, verify `head(rownames(mat))` shows gene symbols and that the organism argument matches.

### Pitfall: padj equals pvalue (no multiple-testing correction)

- **Symptom:** The `padj` column in the master table is identical to `pvalue`.
- **Cause:** `run_ulm()` returns raw p-values without multiple-testing correction. The current script maps `p_value` directly to `padj`, which is technically incorrect.
- **Fix:** Apply correction explicitly if needed:
  ```r
  tf_acts$padj <- p.adjust(tf_acts$p_value, method = "BH")
  ```
  Note: this is a known simplification in the current pipeline. For exploratory analysis, the raw p-values are often sufficient because downstream filtering uses both score magnitude and p-value.

### Pitfall: CollecTRI network download timeout

- **Symptom:** `get_collectri()` hangs or errors with a connection timeout, especially in restricted network environments.
- **Cause:** OmnipathR fetches the network from the OmniPath web service, which requires internet access.
- **Fix:** Cache the network on first successful download (`saveRDS(net, "collectri_mouse.rds")`). On subsequent runs, load from the cache file. The analysis script already implements this pattern.

### Pitfall: Condition column mismatch after filtering

- **Symptom:** `tf_acts %>% filter(condition == "t")` returns zero rows.
- **Cause:** The column name in the input matrix was not `"t"`. When using `dplyr::select(t)`, the column name is preserved. But if the matrix was built differently (e.g., renamed columns), the `condition` value in the output will differ.
- **Fix:** Check available conditions with `unique(tf_acts$condition)` and filter accordingly. The condition name matches the column name of the input matrix exactly.

---

## Complementary Skills

| When you need... | Use skill | Relationship |
|---|---|---|
| Pathway-level activity (PROGENy, 14 pathways) | `progeny-pathway-activity` | Alternative |
| GSEA against MSigDB collections | `bulk-rnaseq-gsea-msigdb` | Alternative |
| Master CSV tables from GSEA results | `bulk-rnaseq-gsea-master-tables` | Next step |
| Visualization of TF activity results | `decoupler-activity-visualization` | Next step |

---

## Resources

- **DecoupleR docs:** https://saezlab.github.io/decoupleR/
- **DecoupleR paper:** Badia-i-Mompel et al. (2022) Bioinformatics, doi:10.1093/bioinformatics/btac651
- **CollecTRI paper:** Muller-Dott et al. (2024) Nucleic Acids Research, doi:10.1093/nar/gkad841
- **OmnipathR docs:** https://saezlab.github.io/OmnipathR/
- **ULM method:** Described in DecoupleR vignette, adapted from Nichols et al. (2006)
