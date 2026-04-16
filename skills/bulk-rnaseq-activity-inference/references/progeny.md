# PROGENy Pathway Activity Inference via DecoupleR MLM

> Part of the bulk-rnaseq-activity-inference pipeline. See `../SKILL.md` for the full routing decision tree.

## Overview

PROGENy (Pathway RespOnsive GENes) infers the activity of 14 canonical signaling pathways from gene expression data using consensus footprint gene signatures. Unlike GSEA, which tests for gene set enrichment using binary membership, PROGENy assigns continuous weights to target genes reflecting how strongly and in which direction each gene responds to pathway perturbation. The DecoupleR Multivariate Linear Model (MLM) is the recommended inference method because it handles these continuous weights natively and accounts for genes that participate in multiple pathways simultaneously.

**The 14 PROGENy pathways:** Androgen, EGFR, Estrogen, Hypoxia, JAK-STAT, MAPK, NFkB, p53, PI3K, TGFb, TNFa, Trail, VEGF, WNT.

**When to use this reference:**
- Computing signaling pathway activity scores from limma-voom or edgeR DE results
- Comparing pathway activation/repression between experimental conditions
- Building `master_progeny_activities.csv` with the `PROGENY_` prefix convention
- Complementing GSEA and TF activity results with focused signaling estimates

**Route elsewhere for:**
- TF activity inference → `references/tf-activity.md`
- Broad gene set enrichment (hundreds/thousands of sets) → `bulk-rnaseq-gsea`
- Visualization → `references/visualization.md`

---

## Decision Tree

```
Need to assess pathway-level changes from DE results?
|
+- Want signaling pathway ACTIVITY scores (14 pathways, continuous weights)?
|  -> This reference (progeny-pathway-activity)
|
+- Want transcription factor ACTIVITY scores (hundreds of TFs)?
|  -> references/tf-activity.md
|
+- Want gene set ENRICHMENT across MSigDB (hundreds/thousands of sets)?
|  -> bulk-rnaseq-gsea (MSigDB collections)
|
+- Want enrichment with custom gene sets (MitoCarta, TransportDB)?
   -> bulk-rnaseq-gsea (custom-db reference)
```

---

## Quick Start

Minimal PROGENy pathway activity inference from DE results:

```r
library(tidyverse)
library(decoupleR)

# Load DE results (gene symbols as rownames, 't' column present)
de_results  <- readRDS("03_results/checkpoints/1.1_de_results.rds")
deg_table   <- de_results[[names(de_results)[1]]]

# Prepare input matrix
mat <- deg_table %>% dplyr::select(t) %>% as.matrix()
mat[is.na(mat)] <- 0

# Get PROGENy network (mouse, top 500 genes per pathway)
net <- decoupleR::get_progeny(organism = "mouse", top = 500)

# Run MLM inference
pathway_acts <- decoupleR::run_mlm(
  mat = mat, net = net,
  .source = "source", .target = "target", .mor = "weight",
  minsize = 5
) %>% filter(condition == "t")
```

**Verify it worked:**

```r
stopifnot(nrow(pathway_acts) == 14)
stopifnot(all(c("score", "p_value", "source") %in% colnames(pathway_acts)))
message("PROGENy returned ", nrow(pathway_acts), " pathway activity scores")
pathway_acts %>% filter(source == "JAK-STAT") %>% print()
```

---

## Progressive Depth

### Basic Usage

#### Input Matrix Preparation

MLM operates on a matrix of t-statistics (genes x contrasts). Gene symbols must be rownames:

```r
mat <- deg_table %>%
  dplyr::select(t) %>%
  as.matrix()
mat[is.na(mat)] <- 0
message("Input: ", nrow(mat), " genes x ", ncol(mat), " contrasts")
```

#### Loading the PROGENy Network

The PROGENy model is a table of gene-pathway associations with continuous weights. Use `top = 500` for bulk RNA-seq (top 500 most responsive genes per pathway):

```r
net <- decoupleR::get_progeny(organism = "mouse", top = 500)
```

The network has three key columns:

| Column | Description |
|--------|-------------|
| `source` | Pathway name (e.g., "JAK-STAT", "MAPK") |
| `target` | Gene symbol |
| `weight` | Continuous response weight (positive = activated by pathway, negative = repressed) |

**Critical:** The weight column is named `weight`, not `mor`. Always use `.mor = "weight"` when passing the PROGENy network to `run_mlm()`.

#### The 14 PROGENy Pathways

| Pathway | Biological Context |
|---------|-------------------|
| Androgen | Androgen receptor signaling |
| EGFR | Epidermal growth factor receptor |
| Estrogen | Estrogen receptor signaling |
| Hypoxia | Oxygen sensing (HIF pathway) |
| JAK-STAT | Cytokine signaling (IL-2, IL-6, IFN) |
| MAPK | Ras-Raf-MEK-ERK cascade |
| NFkB | Inflammatory/immune signaling |
| p53 | DNA damage response, apoptosis |
| PI3K | Growth/survival (AKT/mTOR) |
| TGFb | TGF-beta superfamily signaling |
| TNFa | Tumor necrosis factor signaling |
| Trail | TRAIL-induced apoptosis |
| VEGF | Angiogenesis |
| WNT | Wnt/beta-catenin signaling |

#### Caching the PROGENy Network

```r
progeny_file <- file.path(DIR_CHECKPOINTS, "progeny_mouse.rds")
if (file.exists(progeny_file)) {
  net <- readRDS(progeny_file)
} else {
  net <- decoupleR::get_progeny(organism = "mouse", top = 500)
  saveRDS(net, progeny_file)
}
```

### Intermediate Usage

#### Master Table Schema

PROGENy results are formatted to match the unified schema used by GSEA and TF activity tables:

```r
final_df <- pathway_acts %>%
  left_join(targets_list, by = c("source" = "source")) %>%
  mutate(
    pathway_id        = paste0("PROGENY_", source),   # PROGENY_ prefix convention
    pathway_name      = source,
    database          = "PROGENy",
    contrast          = contrast_name,
    nes               = score,
    pvalue            = p_value,
    padj              = p_value,                      # No MTC needed (only 14 tests)
    set_size          = n_genes,
    leading_edge_size = n_genes,
    direction         = ifelse(score > 0, "Up", "Down"),
    core_enrichment   = genes
  ) %>%
  dplyr::select(pathway_id, pathway_name, database, contrast,
                nes, pvalue, padj, set_size, leading_edge_size,
                direction, core_enrichment) %>%
  arrange(pvalue)

write_csv(final_df, file.path(DIR_TABLES, "master_progeny_activities.csv"))
```

**Key conventions:**
- **`PROGENY_` prefix**: All pathway IDs use `PROGENY_JAK-STAT`, `PROGENY_MAPK`, etc. to distinguish from GSEA pathway IDs in combined views.
- **No multiple testing correction**: With only 14 independent tests, `padj = pvalue` is acceptable. Bonferroni at 14 tests barely differs from nominal.
- **Direction**: `score > 0` means pathway is more active in the numerator group.

#### Significance Markers

```r
final_df %>%
  mutate(
    sig = case_when(
      pvalue < 0.001 ~ "***",
      pvalue < 0.01  ~ "**",
      pvalue < 0.05  ~ "*",
      TRUE           ~ ""
    )
  )
```

#### Checkpoint Caching

```r
pathway_acts <- load_or_compute(
  checkpoint_file = "1.8_progeny_results.rds",
  description     = "PROGENy pathway activities",
  compute_fn      = function() {
    mat <- deg_table %>% dplyr::select(t) %>% as.matrix()
    mat[is.na(mat)] <- 0
    net  <- decoupleR::get_progeny(organism = "mouse", top = 500)
    acts <- decoupleR::run_mlm(
      mat = mat, net = net,
      .source = "source", .target = "target", .mor = "weight",
      minsize = 5
    )
    acts %>% filter(condition == "t")
  }
)
```

### Advanced Usage

#### Relationship to GSEA

PROGENy and GSEA are complementary, not redundant:

| Aspect | PROGENy (MLM) | GSEA (fgsea) |
|--------|---------------|--------------|
| Gene weights | Continuous (signed weights) | Binary (member or not) |
| Gene sets | 14 curated pathways | Hundreds to thousands |
| Multi-pathway genes | Handled via multivariate model | Tested independently per set |
| Scoring | Linear model t-statistic | Normalized enrichment score |
| Interpretation | "Pathway activity change" | "Gene set enrichment" |
| Best for | Signaling pathway activity | Broad functional annotation |

When both are available, compare: a GSEA-significant `HALLMARK_IL2_STAT5_SIGNALING` should be concordant with PROGENy's JAK-STAT activity. Discordance may indicate the gene set captures a different aspect than the signaling footprint.

#### Target Gene Concordance Analysis

```r
jak_targets <- net %>%
  filter(source == "JAK-STAT", target %in% rownames(mat))

jak_de <- deg_table[jak_targets$target, ] %>%
  rownames_to_column("gene") %>%
  left_join(jak_targets %>% select(target, weight), by = c("gene" = "target")) %>%
  mutate(concordant = sign(logFC) == sign(weight))

message("Concordant: ", sum(jak_de$concordant, na.rm = TRUE), "/", nrow(jak_de))
```

#### Species Support

```r
net_mouse <- get_progeny(organism = "mouse", top = 500)  # This project
net_human <- get_progeny(organism = "human", top = 500)
# Use top = 100 for single-cell RNA-seq
```

---

## Verification Checklist

- [ ] **Exactly 14 pathways:** `nrow(pathway_acts) == 14`. Any other count indicates a network loading or filtering problem.
- [ ] **Score and p-value columns present:** `all(c("score", "p_value", "source") %in% colnames(pathway_acts))`
- [ ] **Biological plausibility:** For IL-2 signaling experiments, JAK-STAT should show activity change. For growth factor experiments, check EGFR, MAPK, PI3K.
- [ ] **Master table saved:** `master_progeny_activities.csv` exists and has 14 rows.
- [ ] **PROGENY_ prefix:** All `pathway_id` values start with `PROGENY_`.
- [ ] **No NA scores:** `sum(is.na(pathway_acts$score)) == 0`.

---

## Common Pitfalls

### Pitfall: Wrong .mor column name in run_mlm

- **Symptom:** Error or all-zero scores from `run_mlm()`.
- **Cause:** Using `.mor = "mor"` instead of `.mor = "weight"`. The PROGENy network column is named `weight`.
- **Fix:** Always use `.mor = "weight"` when passing the PROGENy network to `run_mlm()`.

### Pitfall: Only 14 pathways (user expecting hundreds)

- **Symptom:** User expects GSEA-like output with many pathways but gets only 14 rows.
- **Cause:** PROGENy by design covers exactly 14 canonical signaling pathways.
- **Fix:** This is correct behavior. For broad pathway coverage, use `bulk-rnaseq-gsea`. PROGENy complements GSEA with focused, high-confidence signaling activity estimates.

### Pitfall: Network download fails (OmniPath timeout)

- **Symptom:** `get_progeny()` hangs or errors with a connection timeout.
- **Cause:** OmniPath server unreachable or slow from current network.
- **Fix:** Cache the network locally after the first successful download.

### Pitfall: No multiple testing correction applied

- **Symptom:** User questions why `padj` equals `pvalue` in the master table.
- **Cause:** With only 14 tests, formal multiple testing correction is unnecessary.
- **Fix:** `padj = pvalue` is intentional. Bonferroni at n=14 gives threshold 0.05/14 = 0.0036, barely different from nominal for strong effects.

---

## Resources

- **DecoupleR docs:** https://saezlab.github.io/decoupleR/
- **PROGENy paper:** Schubert et al. (2018) Nature Communications, doi:10.1038/s41467-017-02391-6
- **DecoupleR paper:** Badia-i-Mompel et al. (2022) Bioinformatics, doi:10.1093/bioinformatics/btac651
- **OmniPath:** https://omnipathdb.org/
- **PROGENy model details:** https://saezlab.github.io/progeny/
