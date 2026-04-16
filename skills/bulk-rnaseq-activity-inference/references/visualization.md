# DecoupleR / PROGENy Activity Visualization

> Part of the bulk-rnaseq-activity-inference pipeline. See `../SKILL.md` for the full routing decision tree.

## Overview

This reference covers publication-quality static visualization for two types of activity inference results: transcription factor (TF) activities from DecoupleR ULM and signaling pathway activities from PROGENy MLM. Both pipelines produce master CSV tables (`master_tf_activities.csv`, `master_progeny_activities.csv`) that share a common schema, enabling a unified set of visualization patterns — diverging barplots, volcano plots, and direction summaries — with tool-specific extensions (TF distribution histograms, PROGENy target gene scatter plots). All plots use a colorblind-safe blue-white-orange diverging palette and are saved as paired PDF + PNG outputs.

**When to use this reference:**
- Creating barplots of top differentially active TFs or all 14 PROGENy pathways
- Generating volcano plots of activity score vs significance
- Building TF direction summary charts (4-category: Up/Down × Significant/Not)
- Creating PROGENy target gene scatter plots (weight vs t-value with concordance coloring)
- Writing interpretation README files

**Route elsewhere for:**
- Running TF activity inference → `references/tf-activity.md`
- Running PROGENy inference → `references/progeny.md`
- Visualizing GSEA enrichment results (dotplots, running sum) → `bulk-rnaseq-gsea`
- Interactive pathway UMAP explorer → `bulk-rnaseq-pathway-explorer`

---

## Decision Tree

```
Need to visualize activity inference results?
|
+-- What type of activity?
|   |
|   +-- Transcription factor activities (DecoupleR ULM)
|   |   +-- "Which TFs are most differentially active?"    -> Activity barplot (top N by p-value)
|   |   +-- "Which TFs are significant AND large effect?"  -> Volcano plot
|   |   +-- "Overall balance of up vs down TFs?"           -> Direction summary (4-category)
|   |   +-- "Full activity distribution?"                  -> Distribution histogram
|   |
|   +-- PROGENy signaling pathway activities
|       +-- "Which pathways are active/suppressed?"        -> Barplot (all 14, with stars)
|       +-- "How significant are pathway activities?"      -> Volcano plot (all 14 labeled)
|       +-- "Which target genes drive a pathway?"          -> Target gene scatter
|
+-- What format?
    +-- Publication figure (PDF + PNG)  -->  R scripts (2.3 / 2.4)
    +-- Interactive exploration          -->  bulk-rnaseq-pathway-explorer
```

---

## Quick Start

```r
source("02_analysis/config/config.R")
source("02_analysis/config/color_config.R")

library(tidyverse)
library(ggplot2)

# Load TF activity results
tf_results <- read_csv(file.path(DIR_TABLES, "master_tf_activities.csv"),
                       show_col_types = FALSE)

# Quick barplot of top 25 TFs by p-value
top_tfs <- tf_results %>% arrange(pvalue) %>% head(25)

p <- ggplot(top_tfs, aes(x = reorder(pathway_name, nes), y = nes, fill = nes)) +
  geom_bar(stat = "identity", width = 0.8) +
  scale_fill_gradient2(
    low      = DIVERGING_COLORS$negative,   # #2166AC (blue)
    mid      = DIVERGING_COLORS$neutral,    # #F7F7F7 (white)
    high     = DIVERGING_COLORS$positive,   # #B35806 (orange)
    midpoint = 0, name = "Activity\nScore"
  ) +
  coord_flip() +
  labs(title = "Top 25 Differentially Active TFs",
       x = "Transcription Factor", y = "Activity Score (ULM)") +
  theme_minimal(base_size = 12)

ggsave("tf_activity_barplot.pdf", p, width = 10, height = 8, dpi = 300)
ggsave("tf_activity_barplot.png", p, width = 10, height = 8, dpi = 300)
```

**Verify:** PDF and PNG are non-empty; bars are blue (negative) through orange (positive); ordered by activity score.

---

## Progressive Depth

### Basic Usage

#### Shared Schema and Patterns

Both `master_tf_activities.csv` and `master_progeny_activities.csv` share a common column schema:

| Column | Type | Description |
|--------|------|-------------|
| `pathway_id` | character | Unique identifier |
| `pathway_name` | character | Display name (TF symbol or pathway name) |
| `database` | character | Source database (CollecTRI or PROGENy) |
| `nes` | numeric | Activity score (ULM t-value or MLM t-value) |
| `pvalue` | numeric | Raw p-value |
| `padj` | numeric | Adjusted p-value (BH for TF; = pvalue for PROGENy) |
| `direction` | character | "Up" or "Down" (title case, case-sensitive) |
| `core_enrichment` | character | Slash-separated target gene list |
| `contrast` | character | Comparison label |

#### Diverging Color Palette

All activity plots use the colorblind-safe diverging scale from `config/color_config.R`:

| Role | Hex | Color |
|------|-----|-------|
| Negative (suppressed) | `#2166AC` | Blue |
| Neutral (zero) | `#F7F7F7` | White |
| Positive (activated) | `#B35806` | Orange |

```r
# Fill gradient for continuous NES barplots
scale_fill_gradient2(
  low = DIVERGING_COLORS$negative, mid = DIVERGING_COLORS$neutral,
  high = DIVERGING_COLORS$positive, midpoint = 0
)

# Discrete direction coloring for volcano/scatter
scale_color_manual(
  values = c(Up = DIVERGING_COLORS$positive, Down = DIVERGING_COLORS$negative)
)
```

#### Activity Barplot Pattern (TF and PROGENy)

The barplot is the primary ranking visualization. Pattern is identical for both:

```r
plot_data <- results %>%
  arrange(pvalue) %>%
  head(N_TOP) %>%
  mutate(pathway_name = factor(pathway_name, levels = pathway_name[order(nes)]))

p <- ggplot(plot_data, aes(x = pathway_name, y = nes, fill = nes)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_hline(yintercept = 0, linetype = "solid", color = "gray40") +
  scale_fill_gradient2(
    low = DIVERGING_COLORS$negative, mid = DIVERGING_COLORS$neutral,
    high = DIVERGING_COLORS$positive, midpoint = 0
  ) +
  coord_flip() +
  theme_minimal(base_size = 12) +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank())
```

#### Volcano Plot Pattern

```r
results <- results %>%
  mutate(
    neg_log_padj  = -log10(padj),
    is_significant = padj < FDR_CUTOFF,
    label         = if_else(is_significant & abs(nes) >= ACTIVITY_CUTOFF, pathway_name, NA_character_)
  )

p <- ggplot(results, aes(x = nes, y = neg_log_padj, color = direction)) +
  geom_point(aes(alpha = is_significant), size = 2) +
  geom_hline(yintercept = -log10(FDR_CUTOFF), linetype = "dashed", color = "gray50") +
  geom_vline(xintercept = 0, linetype = "solid", color = "gray70") +
  geom_text_repel(aes(label = label), size = 3, max.overlaps = 20, na.rm = TRUE) +
  scale_color_manual(values = c(Up = DIVERGING_COLORS$positive, Down = DIVERGING_COLORS$negative)) +
  scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.3), guide = "none")
```

### Intermediate Usage

#### TF-Specific: Direction Summary (4-Category Bar Chart)

Groups TFs into four categories based on significance and direction:

```r
direction_summary <- tf_results %>%
  group_by(is_significant, direction) %>%
  summarize(count = n(), .groups = "drop") %>%
  mutate(
    category = case_when(
      is_significant & direction == "Up"    ~ "Up (Significant)",
      is_significant & direction == "Down"  ~ "Down (Significant)",
      !is_significant & direction == "Up"   ~ "Up (Not Sig.)",
      !is_significant & direction == "Down" ~ "Down (Not Sig.)"
    ),
    fill_color = case_when(
      direction == "Up"   & is_significant  ~ DIVERGING_COLORS$positive,
      direction == "Down" & is_significant  ~ DIVERGING_COLORS$negative,
      direction == "Up"                     ~ adjustcolor(DIVERGING_COLORS$positive, alpha.f = 0.4),
      direction == "Down"                   ~ adjustcolor(DIVERGING_COLORS$negative, alpha.f = 0.4)
    )
  )
```

#### TF-Specific: Activity Distribution Histogram

```r
p_dist <- ggplot(tf_results, aes(x = nes, fill = direction)) +
  geom_histogram(bins = 50, alpha = 0.7, position = "identity") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  scale_fill_manual(values = c(Up = DIVERGING_COLORS$positive, Down = DIVERGING_COLORS$negative))
```

**Top N selection and volcano label thresholds:**

```r
TF_FDR_CUTOFF      <- 0.05
TF_ACTIVITY_CUTOFF <- 2.0  # Minimum |activity score| to label on volcano

label = case_when(
  is_significant & abs(nes) >= TF_ACTIVITY_CUTOFF ~ pathway_name,
  TRUE ~ NA_character_
)
```

#### PROGENy-Specific: All 14 Pathways (No Filtering)

Unlike TF visualization (top N), PROGENy shows all 14 pathways — no filtering step:

```r
progeny_ordered <- progeny_results %>%
  mutate(pathway_name = factor(pathway_name, levels = pathway_name[order(nes)]))
```

#### PROGENy-Specific: Significance Stars at Bar Tips

```r
progeny_ordered <- progeny_ordered %>%
  mutate(
    sig_label = case_when(
      pvalue < 0.001 ~ "***",
      pvalue < 0.01  ~ "**",
      pvalue < 0.05  ~ "*",
      TRUE           ~ ""
    )
  )

# In ggplot:
geom_text(aes(label = sig_label, y = nes + sign(nes) * 0.8), size = 5)
```

Note: PROGENy uses raw `pvalue` for stars (not `padj`) — with 14 pathways the multiple testing burden is minimal.

#### PROGENy-Specific: Target Gene Scatter (Weight vs T-Value)

For each top pathway, a scatter shows PROGENy gene weights (x) vs DE t-statistics (y), colored by concordance:

```r
df_genes <- pathway_genes %>%
  inner_join(deg_table %>% rownames_to_column("gene") %>% select(gene, t, logFC, adj.P.Val),
             by = c("target" = "gene")) %>%
  mutate(
    effect = case_when(
      weight > 0 & t > 0 ~ "Concordant (++)",
      weight < 0 & t < 0 ~ "Concordant (--)",
      weight > 0 & t < 0 ~ "Discordant (+-)",
      weight < 0 & t > 0 ~ "Discordant (-+)",
      TRUE ~ "Neutral"
    )
  )

concordant_colors <- c(
  "Concordant (++)" = "#2ca02c",
  "Concordant (--)" = "#2ca02c",
  "Discordant (+-)" = "#d62728",
  "Discordant (-+)" = "#d62728",
  "Neutral"         = "gray50"
)
```

Top 20 genes by absolute t-value are labeled with `ggrepel::geom_text_repel()`.

#### Output Organization

```
03_results/plots/
  TF/
    tf_activity_barplot.pdf / .png     # Top 25 TFs by p-value
    tf_activity_volcano.pdf / .png     # Activity vs significance
    tf_direction_summary.pdf / .png    # 4-category direction chart
    tf_activity_distribution.pdf / .png
    README.md                          # Key statistics and interpretation
  PROGENy/
    progeny_activity_barplot.pdf / .png    # All 14 pathways with stars
    progeny_activity_volcano.pdf / .png    # Activity vs significance
    progeny_targets_*.pdf / .png           # Target gene scatter (top 4 pathways)
    README.md                              # Results table and interpretation
```

Publication standards: 300 DPI PNG, PDF vector, default 10 x 8 inches (PROGENy barplot 9 x 7).

#### README Generation

Both scripts generate README.md files containing:
- File manifest with descriptions
- Key statistics (total tested, significant count, up/down breakdown)
- Top upregulated and downregulated entries
- Method description (ULM for TF, MLM for PROGENy)
- For PROGENy: markdown table of all 14 pathways with scores, p-values, and directions

### Advanced Usage

#### Color Helper Functions

`config/color_config.R` provides reusable helpers:

```r
nes_ggplot_fill_scale(limits = c(-3.5, 3.5), name = "NES")  # fill aesthetic
nes_ggplot_scale(limits = c(-3.5, 3.5), name = "NES")       # color aesthetic
col_fun <- nes_color_scale(limits = c(-3.5, 3.5), n_colors = 3)  # ComplexHeatmap
```

#### Integration with Pathway Explorer

Both master activity CSVs feed into `bulk-rnaseq-pathway-explorer`, where TF and PROGENy entries appear alongside GSEA pathways in the UMAP embedding:
- Shape coding: diamond for TFs, square for PROGENy, circle for GSEA pathways
- Same blue-white-orange NES coloring
- Hybrid similarity: TF-Pathway uses Overlap coefficient, PROGENy-Pathway uses Jaccard

No extra visualization code needed — the explorer auto-detects entries from the master CSVs.

---

## Verification Checklist

- [ ] **TF plot files exist:** `tf_activity_barplot.pdf`, `tf_activity_volcano.pdf`, `tf_direction_summary.pdf`, `tf_activity_distribution.pdf` in `plots/TF/`
- [ ] **PROGENy plot files exist:** `progeny_activity_barplot.pdf`, `progeny_activity_volcano.pdf`, `progeny_targets_*.pdf` in `plots/PROGENy/`
- [ ] **PDF/PNG paired:** Every `.pdf` has a corresponding `.png` file
- [ ] **READMEs generated:** `plots/TF/README.md` and `plots/PROGENy/README.md` exist with correct statistics
- [ ] **NES color scale centered:** Barplot fill gradient symmetric around zero (blue-white-orange)
- [ ] **PROGENy shows all 14:** Barplot includes all 14 signaling pathways, not a filtered subset
- [ ] **Significance markers consistent:** PROGENy uses stars (* p<0.05, ** p<0.01, *** p<0.001)

---

## Common Pitfalls

### Pitfall: ggplot2 4.0+ NA Color Removes Points

- **Symptom:** Shape 21 points disappear; console warns "removed rows containing missing values."
- **Cause:** ggplot2 4.0+ treats `color = NA` as a missing aesthetic, removing the geom entirely.
- **Fix:**
  ```r
  # WRONG (ggplot2 4.0+)
  geom_point(shape = 21, color = NA, fill = "blue")
  # CORRECT
  geom_point(shape = 21, stroke = 0, fill = "blue")
  ```

### Pitfall: Missing Master CSV

- **Symptom:** Script exits with "TF results file not found" or "PROGENy results not found."
- **Cause:** Activity inference scripts have not been run.
- **Fix:**
  ```bash
  Rscript 02_analysis/1.7.decoupler_tf_analysis.R   # produces master_tf_activities.csv
  Rscript 02_analysis/1.8.progeny_analysis.R         # produces master_progeny_activities.csv
  ```

### Pitfall: Significance Threshold Inconsistency

- **Symptom:** TF volcano labels do not match expected count; README statistics disagree with plot.
- **Cause:** `TF_FDR_CUTOFF` in the visualization script diverges from threshold in the analysis script.
- **Fix:** Ensure `TF_FDR_CUTOFF` in `2.3.tf_viz.R` matches the threshold in `1.7.decoupler_tf_analysis.R` and `config/pipeline.yaml`.

### Pitfall: Direction Column Case Sensitivity

- **Symptom:** `scale_color_manual()` warns about unused values; points have unexpected colors.
- **Cause:** `direction` column must contain exactly "Up" or "Down" (title case).
- **Fix:**
  ```r
  stopifnot(all(results$direction %in% c("Up", "Down")))
  # If values differ: mutate(direction = if_else(nes > 0, "Up", "Down"))
  ```

### Pitfall: PROGENy Barplot Color Scale Asymmetry

- **Symptom:** Barplot gradient appears skewed — one direction more saturated.
- **Cause:** `scale_fill_gradient2()` limits are set asymmetrically.
- **Fix:** Use symmetric limits based on the maximum absolute value:
  ```r
  nes_max <- max(abs(progeny_ordered$nes)) + 1
  scale_fill_gradient2(..., limits = c(-nes_max, nes_max))
  ```

---

## Resources

- **DecoupleR:** https://saezlab.github.io/decoupleR/
- **PROGENy:** https://saezlab.github.io/progeny/
- **ggrepel:** https://ggrepel.slowkow.com/
- **Color config:** `02_analysis/config/color_config.R`
- **TF viz script:** `02_analysis/2.3.tf_viz.R`
- **PROGENy viz script:** `02_analysis/2.4.progeny_viz.R`
