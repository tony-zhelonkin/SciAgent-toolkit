# DecoupleR Transcription Factor Activity Inference

> Part of the bulk-rnaseq-activity-inference pipeline. See `../SKILL.md` for the full routing decision tree.

## Overview

DecoupleR infers transcription factor (TF) activities from gene expression data by leveraging prior-knowledge regulatory networks. The Univariate Linear Model (ULM) method regresses each TF's known target gene weights against the observed expression changes (t-statistics from DE analysis), producing a per-TF activity score and p-value. This approach is fundamentally different from GSEA — ULM uses the continuous mode-of-regulation (MoR) weights from CollecTRI (activator=+1, repressor=-1, with intermediate values), while GSEA treats gene sets as unweighted membership lists.

**When to use this reference:**
- Inferring which transcription factors are differentially active between conditions
- Building `master_tf_activities.csv` from limma/edgeR DE results
- Working with the CollecTRI regulatory network via OmnipathR
- Needing continuous activity scores (not just enrichment p-values) for TFs

**Route elsewhere for:**
- Pathway-level activity (PROGENy/14 pathways) → `references/progeny.md`
- Visualization of TF activity results → `references/visualization.md`
- GSEA-based TF target enrichment (C3 TFT:GTRD collection) → use `bulk-rnaseq-gsea`

---

## Decision Tree

```
Need to identify active transcription factors from DE results?
|
+- Want TF-level activity scores with regulatory weights?
|  -> This reference (DecoupleR ULM + CollecTRI)
|
+- Want pathway-level activity scores (14 cancer pathways)?
|  -> references/progeny.md (DecoupleR + PROGENy)
|
+- Want enrichment of curated TF target gene sets (GTRD, ENCODE)?
|  -> bulk-rnaseq-gsea (C3 TFT:GTRD collection)
|
+- Want to visualize TF activity results?
   -> references/visualization.md
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
deg_table  <- de_results[["IL2RAKO_vs_NTC"]]

# Prepare input matrix: genes x contrasts, using t-statistics
mat        <- deg_table %>% dplyr::select(t) %>% as.matrix()
mat[is.na(mat)] <- 0

# Fetch CollecTRI network for mouse
net <- decoupleR::get_collectri(organism = "mouse", split_complexes = FALSE)

# Run ULM
tf_acts <- decoupleR::run_ulm(
  mat = mat, net = net,
  .source = "source", .target = "target", .mor = "mor",
  minsize = 5
)

# Filter to the contrast column
tf_acts <- tf_acts %>% filter(condition == "t")
```

**Verify it worked:**

```r
stopifnot(nrow(tf_acts) > 0)
stopifnot(all(c("source", "score", "p_value") %in% colnames(tf_acts)))
message("Inferred activities for ", nrow(tf_acts), " TFs")
stopifnot(nrow(tf_acts) >= 50 && nrow(tf_acts) <= 2000)
stopifnot(abs(mean(tf_acts$score)) < 5)
```

---

## Progressive Depth

### Basic Usage

#### Input Matrix Preparation

ULM expects a numeric matrix with genes as rows and conditions/samples as columns. For DE-based TF activity, the t-statistic is the ideal input because it captures both direction and precision of expression change:

```r
mat <- deg_table %>% dplyr::select(t) %>% as.matrix()
mat[is.na(mat)] <- 0
# Dimensions: genes x 1 (single contrast)
# rownames: gene symbols (must match network target names)
# colnames: "t" (becomes the 'condition' in output)
```

Why t-statistics and not logFC? The t-statistic down-weights genes with high variance, giving more reliable TF activity estimates.

#### Loading CollecTRI

CollecTRI is a comprehensive collection of TF-target interactions curated from literature and ChIP-seq data, accessed via OmnipathR. The network is a three-column data frame:

| Column | Description |
|--------|-------------|
| `source` | Transcription factor name (gene symbol) |
| `target` | Target gene name (gene symbol) |
| `mor` | Mode of regulation: +1 (activator), -1 (repressor), or intermediate values |

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
  .mor    = "mor",      # Column in net with mode of regulation weights
  minsize = 5           # Minimum regulon size (TFs with fewer targets are dropped)
)
```

The `minsize = 5` parameter filters out TFs with fewer than 5 targets present in the input matrix. Increasing to 10 gives more robust but fewer TFs.

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

To integrate TF activities into the pathway explorer and master table ecosystem, map ULM output columns to the standard schema:

```r
# Get all regulon targets present in our data for core_enrichment
targets_list <- net %>%
  filter(source %in% tf_acts$source) %>%
  filter(target %in% rownames(mat)) %>%
  group_by(source) %>%
  summarise(genes = paste(unique(target), collapse = "/"))

# Apply BH correction
tf_acts$padj <- p.adjust(tf_acts$p_value, method = "BH")

# Build master table
final_df <- tf_acts %>%
  left_join(targets_list, by = "source") %>%
  mutate(
    pathway_id        = source,
    pathway_name      = source,
    database          = "CollecTRI",
    contrast          = "IL2RAKO_vs_NTC",
    nes               = score,
    pvalue            = p_value,
    set_size          = NA_integer_,
    leading_edge_size = NA_integer_,
    gene_ratio        = NA_character_,
    core_enrichment   = genes,
    direction         = ifelse(score > 0, "Up", "Down"),
    neg_log_padj      = -log10(padj)
  ) %>%
  dplyr::select(pathway_id, pathway_name, nes, pvalue, padj, set_size,
                leading_edge_size, gene_ratio, core_enrichment, database,
                contrast, neg_log_padj, direction)

write_csv(final_df, file.path(DIR_TABLES, "master_tf_activities.csv"))
```

The `core_enrichment` field uses the full regulon (all targets in the dataset), not a leading-edge subset. This is intentional — for Jaccard similarity in the pathway explorer, broader gene lists produce more meaningful overlap scores.

#### Checkpoint Caching Pattern

```r
source("02_analysis/config/config.R")

tf_results <- load_or_compute(
  checkpoint_file = "1.7_tf_results.rds",
  description     = "DecoupleR TF activity inference (ULM + CollecTRI)",
  compute_fn      = function() {
    net     <- get_or_cache_collectri()
    mat     <- prepare_tstat_matrix(de_results[[contrast_name]])
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

```r
contrasts <- names(de_results)
mat <- sapply(contrasts, function(cname) {
  tvals         <- de_results[[cname]]$t
  names(tvals)  <- rownames(de_results[[cname]])
  tvals
})
mat[is.na(mat)] <- 0

# run_ulm handles multi-column matrices; 'condition' column contains contrast names
tf_acts <- decoupleR::run_ulm(
  mat = mat, net = net,
  .source = "source", .target = "target", .mor = "mor",
  minsize = 5
)
tf_acts %>% filter(condition == "IL2RAKO_vs_NTC")
```

#### Alternative Network: DoRothEA

DoRothEA provides confidence-level-filtered TF regulons (A-E). Use when stricter evidence filtering is needed (common in single-cell workflows):

```r
net_dorothea <- decoupleR::get_dorothea(organism = "mouse", levels = c("A", "B"))
tf_acts_dorothea <- decoupleR::run_ulm(
  mat = mat, net = net_dorothea,
  .source = "source", .target = "target", .mor = "mor",
  minsize = 5
)
```

CollecTRI is preferred for bulk RNA-seq: broader coverage and continuous MoR weights. DoRothEA is more common in single-cell workflows.

#### Species Handling

`get_collectri()` accepts short-form species names, not binomial nomenclature:

| Species | Argument |
|---------|----------|
| Mouse | `organism = "mouse"` |
| Human | `organism = "human"` |
| Rat | `organism = "rat"` |

---

## Verification Checklist

- [ ] **TF count reasonable:** 100-1500 TFs for mouse/human with `minsize = 5`. Fewer than 50 suggests gene symbol mismatch.
- [ ] **Activity distribution centered near 0:** `mean(tf_acts$score)` within +/-5.
- [ ] **Known biology validation:** Expected TFs show activity (e.g., STAT5 for IL-2 signaling experiments).
- [ ] **P-value distribution:** `hist(tf_acts$p_value)` shows roughly uniform distribution with spike near 0.
- [ ] **BH correction applied:** `padj` column is not identical to `pvalue`.
- [ ] **Checkpoint exists:** `file.exists("03_results/checkpoints/1.7_tf_results.rds")`.
- [ ] **Master table written:** `master_tf_activities.csv` exists with correct 13-column schema.

---

## Common Pitfalls

### Pitfall: Gene symbol mismatch with network

- **Symptom:** Very few TFs inferred (under 50), or `run_ulm()` returns empty results.
- **Cause:** Rownames of the input matrix do not match the `target` column in CollecTRI. Common when DE results use Ensembl IDs, or human symbols are used with `organism = "mouse"`.
- **Fix:** `sum(rownames(mat) %in% net$target)` should be several thousand. Verify `head(rownames(mat))` shows gene symbols.

### Pitfall: padj equals pvalue (no BH correction)

- **Symptom:** The `padj` column in the master table is identical to `pvalue`.
- **Cause:** `run_ulm()` returns raw p-values without multiple-testing correction.
- **Fix:** `tf_acts$padj <- p.adjust(tf_acts$p_value, method = "BH")`

### Pitfall: CollecTRI network download timeout

- **Symptom:** `get_collectri()` hangs or errors with a connection timeout.
- **Cause:** OmnipathR fetches from the OmniPath web service — requires internet access.
- **Fix:** Cache after first successful download. The analysis script implements this pattern.

### Pitfall: Condition column mismatch after filtering

- **Symptom:** `tf_acts %>% filter(condition == "t")` returns zero rows.
- **Cause:** The column name in the input matrix was not `"t"` — the `condition` value in output matches the column name exactly.
- **Fix:** Check `unique(tf_acts$condition)` and filter accordingly.

---

## Resources

- **DecoupleR docs:** https://saezlab.github.io/decoupleR/
- **DecoupleR paper:** Badia-i-Mompel et al. (2022) Bioinformatics, doi:10.1093/bioinformatics/btac651
- **CollecTRI paper:** Muller-Dott et al. (2024) Nucleic Acids Research, doi:10.1093/nar/gkad841
- **OmnipathR docs:** https://saezlab.github.io/OmnipathR/
