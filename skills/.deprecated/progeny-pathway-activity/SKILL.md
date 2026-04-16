---
name: progeny-pathway-activity
description: "PROGENy pathway activity inference -- infers signaling pathway activities from bulk RNA-seq DE results using DecoupleR MLM with PROGENy consensus gene signatures (14 canonical pathways). Use when computing pathway activity scores from limma/edgeR DE output, comparing signaling pathway activation between conditions, or building master pathway activity tables. For TF-level activity use decoupler-tf-activity. For gene set enrichment use bulk-rnaseq-gsea-msigdb."
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
    - progeny
    - pathway-activity
    - mlm
    - signaling-pathways
    - bulk-rnaseq
    - jak-stat
    - mapk
    - pi3k
    - nfkb
  complementary-skills:
    - decoupler-tf-activity
    - bulk-rnaseq-gsea-msigdb
    - bulk-rnaseq-gsea-master-tables
    - decoupler-activity-visualization
  contraindications:
    - "Do not use for TF activity inference. Use decoupler-tf-activity instead."
    - "Do not use for gene set enrichment (GSEA). Use bulk-rnaseq-gsea-msigdb instead."
    - "Do not use for visualization. Use decoupler-activity-visualization instead."
---

# PROGENy Pathway Activity Inference via DecoupleR MLM

## Overview

PROGENy (Pathway RespOnsive GENes) infers the activity of 14 canonical signaling pathways from gene expression data using consensus footprint gene signatures. Unlike GSEA, which tests for gene set enrichment using binary membership, PROGENy assigns continuous weights to target genes reflecting how strongly and in which direction each gene responds to pathway perturbation. The DecoupleR Multivariate Linear Model (MLM) is the recommended inference method because it handles these continuous weights natively and accounts for genes that participate in multiple pathways simultaneously.

The 14 PROGENy pathways are: **Androgen, EGFR, Estrogen, Hypoxia, JAK-STAT, MAPK, NFkB, p53, PI3K, TGFb, TNFa, Trail, VEGF, WNT**.

**When to use this skill:**
- Computing signaling pathway activity scores from limma-voom or edgeR DE results
- Comparing pathway activation/repression between experimental conditions
- Building a master pathway activity table with the `PROGENY_` prefix convention
- Complementing GSEA results with a footprint-based activity estimate

**When NOT to use this skill:**
- Transcription factor activity inference -> use `decoupler-tf-activity`
- Gene set enrichment analysis (hundreds/thousands of gene sets) -> use `bulk-rnaseq-gsea-msigdb`
- Plotting pathway activity results -> use `decoupler-activity-visualization`

---

## Decision Tree

```
Need to assess pathway-level changes from DE results?
|
+- Want signaling pathway ACTIVITY scores (14 pathways, continuous weights)?
|  -> This skill (progeny-pathway-activity)
|
+- Want transcription factor ACTIVITY scores (hundreds of TFs)?
|  -> decoupler-tf-activity
|
+- Want gene set ENRICHMENT across MSigDB (hundreds/thousands of sets)?
|  -> bulk-rnaseq-gsea-msigdb
|
+- Want enrichment with custom gene sets (MitoCarta, TransportDB)?
|  -> bulk-rnaseq-gsea-custom-db
```

---

## Quick Start

Minimal PROGENy pathway activity inference from DE results:

```r
library(tidyverse)
library(decoupleR)

# Load DE results (gene symbols as rownames, 't' column present)
de_results <- readRDS("03_results/checkpoints/1.1_de_results.rds")
deg_table <- de_results[[names(de_results)[1]]]

# Prepare input matrix: genes x contrasts (t-statistics)
mat <- deg_table %>% dplyr::select(t) %>% as.matrix()
mat[is.na(mat)] <- 0

# Get PROGENy network (mouse, top 500 genes per pathway)
net <- decoupleR::get_progeny(organism = "mouse", top = 500)

# Run MLM inference
pathway_acts <- decoupleR::run_mlm(
  mat = mat, net = net,
  .source = "source", .target = "target", .mor = "weight",
  minsize = 5
)

# Filter to contrast column
pathway_acts <- pathway_acts %>% filter(condition == "t")
```

**Verify it worked:**

```r
# Must produce exactly 14 pathway results
stopifnot(nrow(pathway_acts) == 14)
stopifnot(all(c("score", "p_value", "source") %in% colnames(pathway_acts)))
message("PROGENy returned ", nrow(pathway_acts), " pathway activity scores")

# Check known biology: JAK-STAT should appear for IL-2 signaling experiments
pathway_acts %>% filter(source == "JAK-STAT") %>% print()
```

---

## Progressive Depth

### Basic Usage

#### Input Matrix Preparation

MLM operates on a matrix of t-statistics (genes x contrasts). The t-statistic is preferred over logFC because it incorporates estimation precision. Gene symbols must be rownames:

```r
# DE results from limma::topTable(number=Inf, sort.by="none")
mat <- deg_table %>%
  dplyr::select(t) %>%
  as.matrix()

# Replace NA t-values with 0 (genes with no test result)
mat[is.na(mat)] <- 0

# Verify: should have thousands of genes, 1 column
message("Input: ", nrow(mat), " genes x ", ncol(mat), " contrasts")
```

#### Loading the PROGENy Network

The PROGENy model is a table of gene-pathway associations with continuous weights. Use `top = 500` for bulk RNA-seq (the top 500 most responsive genes per pathway):

```r
net <- decoupleR::get_progeny(organism = "mouse", top = 500)
```

The network has three key columns:

| Column | Description |
|--------|-------------|
| `source` | Pathway name (e.g., "JAK-STAT", "MAPK") |
| `target` | Gene symbol |
| `weight` | Continuous response weight (positive = activated by pathway, negative = repressed) |

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

The `get_progeny()` call downloads from OmniPath. Cache locally to avoid repeated network calls:

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

PROGENy results are formatted to match the same schema used by GSEA and TF activity master tables, enabling unified downstream analysis:

```r
final_df <- pathway_acts %>%
  left_join(targets_list, by = c("source" = "source")) %>%
  mutate(
    pathway_id = paste0("PROGENY_", source),    # PROGENY_ prefix convention
    pathway_name = source,
    database = "PROGENy",
    contrast = contrast_name,
    nes = score,                                 # MLM score as NES equivalent
    pvalue = p_value,
    padj = p_value,                              # No MTC needed (only 14 tests)
    set_size = n_genes,
    leading_edge_size = n_genes,                 # All target genes used in MLM
    direction = ifelse(score > 0, "Up", "Down"),
    core_enrichment = genes                      # Slash-separated target genes
  ) %>%
  dplyr::select(
    pathway_id, pathway_name, database, contrast,
    nes, pvalue, padj, set_size, leading_edge_size,
    direction, core_enrichment
  ) %>%
  arrange(pvalue)

write_csv(final_df, file.path(DIR_TABLES, "master_progeny_activities.csv"))
```

**Key conventions:**

- **`PROGENY_` prefix**: All pathway IDs use `PROGENY_JAK-STAT`, `PROGENY_MAPK`, etc. This distinguishes them from GSEA pathway IDs in combined views.
- **Direction assignment**: `score > 0` means pathway is more active in the numerator group (e.g., IL2RA_KO); `score < 0` means less active.
- **No multiple testing correction needed**: With only 14 independent tests, `padj = pvalue` is acceptable. Bonferroni at 14 tests barely differs from nominal.

#### Significance Markers

Apply star markers for display:

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

#### Checkpoint Caching with load_or_compute

Wrap the full analysis in `load_or_compute()` for pipeline integration:

```r
source("02_analysis/config/config.R")

pathway_acts <- load_or_compute(
  checkpoint_file = "1.8_progeny_results.rds",
  description = "PROGENy pathway activities",
  compute_fn = function() {
    mat <- deg_table %>% dplyr::select(t) %>% as.matrix()
    mat[is.na(mat)] <- 0
    net <- decoupleR::get_progeny(organism = "mouse", top = 500)
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

When both are available, compare: a GSEA-significant pathway (e.g., HALLMARK_IL2_STAT5_SIGNALING) should be concordant with the corresponding PROGENy activity (JAK-STAT). Discordance may indicate the gene set captures a different aspect of the pathway than the signaling footprint.

#### Target Gene Concordance Analysis

To check which PROGENy target genes are actually DE in your data:

```r
# Get target genes for a specific pathway present in our data
jak_targets <- net %>%
  filter(source == "JAK-STAT", target %in% rownames(mat))

# Check their DE status
jak_de <- deg_table[jak_targets$target, ] %>%
  rownames_to_column("gene") %>%
  left_join(jak_targets %>% select(target, weight), by = c("gene" = "target")) %>%
  mutate(concordant = sign(logFC) == sign(weight))

# Proportion of concordant genes supports the activity call
message("Concordant: ", sum(jak_de$concordant, na.rm = TRUE), "/", nrow(jak_de))
```

#### Species Support

PROGENy supports human and mouse via `get_progeny()`:

```r
# Mouse (this project)
net_mouse <- get_progeny(organism = "mouse", top = 500)

# Human
net_human <- get_progeny(organism = "human", top = 500)
```

The `top` parameter controls how many footprint genes per pathway to include. Use `top = 500` for bulk RNA-seq and `top = 100` for single-cell RNA-seq (where fewer genes are reliably detected per cell).

#### Script Organization

| File | Role |
|------|------|
| `02_analysis/config/config.R` | Defines directories, `load_or_compute()`, checkpoints |
| `02_analysis/1.8.progeny_analysis.R` | Full PROGENy pipeline (load, infer, format, export) |
| `03_results/checkpoints/1.8_progeny_results.rds` | Cached raw MLM results |
| `03_results/tables/master_progeny_activities.csv` | Final master table |

---

## Verification Checklist

After running PROGENy analysis, confirm:

- [ ] **Exactly 14 pathways**: `nrow(pathway_acts) == 14` -- any other count indicates a network loading or filtering problem
- [ ] **Score and p-value columns present**: `all(c("score", "p_value", "source") %in% colnames(pathway_acts))`
- [ ] **Biological plausibility**: For IL-2 signaling experiments, JAK-STAT should show activity change. For growth factor experiments, check EGFR, MAPK, PI3K.
- [ ] **Master table saved**: `file.exists("03_results/tables/master_progeny_activities.csv")` and has 14 rows
- [ ] **PROGENY_ prefix**: All `pathway_id` values start with `PROGENY_`
- [ ] **No NA scores**: `sum(is.na(pathway_acts$score)) == 0`

---

## Common Pitfalls

### Pitfall: Missing DE checkpoint

- **Symptom:** `Error: DE results checkpoint not found` when running `1.8.progeny_analysis.R`.
- **Cause:** The core pipeline (`1.1.core_pipeline.R`) has not been run, so `1.1_de_results.rds` does not exist.
- **Fix:** Run the core pipeline first: `Rscript 02_analysis/1.1.core_pipeline.R`. Verify with `file.exists(file.path(DIR_CHECKPOINTS, "1.1_de_results.rds"))`.

### Pitfall: Only 14 pathways (expecting hundreds)

- **Symptom:** User expects many pathway results like GSEA but gets only 14 rows.
- **Cause:** PROGENy by design covers exactly 14 canonical signaling pathways. It is not a gene set enrichment method.
- **Fix:** This is correct behavior. For broad pathway coverage, use GSEA (`bulk-rnaseq-gsea-msigdb`). PROGENy complements GSEA with focused, high-confidence signaling activity estimates.

### Pitfall: Wrong .mor column name in run_mlm

- **Symptom:** Error or all-zero scores from `run_mlm()`.
- **Cause:** Using `.mor = "mor"` instead of `.mor = "weight"`. The PROGENy network from `get_progeny()` stores weights in a column named `weight`, not `mor`.
- **Fix:** Always use `.mor = "weight"` when passing the PROGENy network to `run_mlm()`.

### Pitfall: No multiple testing correction applied

- **Symptom:** User applies BH or Bonferroni correction to 14 p-values and questions why `padj` equals `pvalue` in the master table.
- **Cause:** With only 14 tests, formal multiple testing correction is unnecessary. The Bonferroni threshold at alpha=0.05 is 0.05/14 = 0.0036, which barely differs from nominal significance for strong effects.
- **Fix:** Setting `padj = pvalue` is intentional. If strict correction is desired, `p.adjust(pvalue, method = "bonferroni")` can be applied, but it changes few conclusions at n=14.

### Pitfall: Network download fails (OmniPath timeout)

- **Symptom:** `get_progeny()` hangs or errors with a connection timeout.
- **Cause:** OmniPath server is unreachable or slow from the current network.
- **Fix:** Cache the network locally after the first successful download (see Basic Usage: Caching the PROGENy Network). If the cached RDS exists, the download is skipped entirely.

---

## Complementary Skills

| When you need... | Use skill | Relationship |
|---|---|---|
| TF activity inference (CollecTRI/DoRothEA) | `decoupler-tf-activity` | Alternative (TF-level, not pathway-level) |
| Broad gene set enrichment (MSigDB) | `bulk-rnaseq-gsea-msigdb` | Extension (complementary analysis) |
| Master table creation/export | `bulk-rnaseq-gsea-master-tables` | Next step |
| Plotting pathway activities | `decoupler-activity-visualization` | Next step |

---

## Resources

- **DecoupleR docs:** https://saezlab.github.io/decoupleR/
- **PROGENy paper:** Schubert et al. (2018) Nature Communications, doi:10.1038/s41467-017-02391-6
- **DecoupleR paper:** Badia-i-Mompel et al. (2022) Bioinformatics, doi:10.1093/bioinformatics/btac651
- **OmniPath:** https://omnipathdb.org/
- **PROGENy model details:** https://saezlab.github.io/progeny/
