---
name: decoupler-activity-visualization
description: "DecoupleR activity visualization -- creates publication-quality plots for TF and pathway activity results from DecoupleR/PROGENy analysis. Use when visualizing TF activity scores (barplots, volcano, distribution), PROGENy pathway activities (barplots, volcano, target gene scatter), or generating interpretation README files. For running the activity inference use decoupler-tf-activity or progeny-pathway-activity."
license: MIT
metadata:
  skill-author: SciAgent-toolkit
  last-reviewed: 2026-04-14
  version: 1.0.0
  upstream-docs: https://saezlab.github.io/decoupleR/

  category: workflow
  tier: standard

  tags:
    - decoupler
    - visualization
    - tf-activity
    - progeny
    - barplot
    - volcano-plot
    - publication-quality
    - ggplot2
    - colorblind-safe
    - bulk-rnaseq

  complementary-skills:
    - decoupler-tf-activity
    - progeny-pathway-activity
    - bulk-rnaseq-gsea-visualization

  contraindications:
    - "Do not use for running activity inference. Use decoupler-tf-activity or progeny-pathway-activity instead."
    - "Do not use for GSEA visualization. Use bulk-rnaseq-gsea-visualization instead."
---

# DecoupleR / PROGENy Activity Visualization

## Overview

This skill covers publication-quality static visualization for two types of activity inference results: transcription factor (TF) activities from DecoupleR's ULM method and signaling pathway activities from PROGENy's MLM method. Both pipelines produce master CSV tables (`master_tf_activities.csv`, `master_progeny_activities.csv`) that follow a shared schema, enabling a unified set of visualization patterns -- diverging barplots, volcano plots, and direction summaries -- with tool-specific extensions (TF distribution histograms, PROGENy target gene scatter plots). All plots use a colorblind-safe blue-white-orange diverging palette and are saved as paired PDF + PNG outputs.

**When to use this skill:**
- Creating barplots of top differentially active TFs or all 14 PROGENy pathways
- Generating volcano plots of activity score vs significance for TFs or pathways
- Building TF direction summary charts (4-category: Up/Down x Significant/Not)
- Creating PROGENy target gene scatter plots (weight vs t-value with concordance coloring)
- Writing interpretation README files with key statistics

**When NOT to use this skill:**
- Running TF activity inference (DecoupleR ULM) -> use `decoupler-tf-activity`
- Running PROGENy pathway activity inference -> use `progeny-pathway-activity`
- Visualizing GSEA enrichment results (dotplots, running sum) -> use `bulk-rnaseq-gsea-visualization`

---

## Decision Tree

```
Need to visualize activity inference results?
|
+-- What type of activity?
|   |
|   +-- Transcription factor activities (DecoupleR ULM)
|   |   |
|   |   +-- "Which TFs are most differentially active?"
|   |   |   --> Activity barplot (top N by p-value, fill by NES)
|   |   |
|   |   +-- "Which TFs are significant AND have large effect?"
|   |   |   --> Volcano plot (NES vs -log10(FDR), label |activity| >= 2.0)
|   |   |
|   |   +-- "What is the overall balance of up vs down TFs?"
|   |   |   --> Direction summary bar chart (4-category)
|   |   |
|   |   +-- "What does the full activity distribution look like?"
|   |       --> Activity distribution histogram (overlaid by direction)
|   |
|   +-- PROGENy signaling pathway activities
|       |
|       +-- "Which pathways are active/suppressed?"
|       |   --> Activity barplot (all 14 pathways, significance stars)
|       |
|       +-- "How significant are the pathway activities?"
|       |   --> Volcano plot (NES vs -log10(pvalue), label all 14)
|       |
|       +-- "Which target genes drive a pathway's activity?"
|           --> Target gene scatter (PROGENy weight vs DE t-value)
|
+-- What format?
    +-- Publication figure (PDF + PNG)  -->  R scripts (2.3 / 2.4)
    +-- Interactive exploration          -->  Pathway explorer (see bulk-rnaseq-gsea-visualization)
```

---

## Quick Start

```r
# Source configuration and color palette
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
    low = DIVERGING_COLORS$negative,   # #2166AC (blue)
    mid = DIVERGING_COLORS$neutral,    # #F7F7F7 (white)
    high = DIVERGING_COLORS$positive,  # #B35806 (orange)
    midpoint = 0, name = "Activity\nScore"
  ) +
  coord_flip() +
  labs(title = "Top 25 Differentially Active TFs",
       x = "Transcription Factor", y = "Activity Score (ULM)") +
  theme_minimal(base_size = 12)

ggsave("tf_activity_barplot.pdf", p, width = 10, height = 8, dpi = 300)
ggsave("tf_activity_barplot.png", p, width = 10, height = 8, dpi = 300)
```

**Verify it worked:**
- PDF and PNG files are non-empty and contain bars colored blue (negative) through orange (positive)
- Bars are ordered by activity score (most negative at bottom, most positive at top)
- The top TFs by p-value are visible regardless of significance status

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
| `padj` | numeric | FDR-adjusted p-value |
| `direction` | character | "Up" or "Down" (case-sensitive) |
| `core_enrichment` | character | Slash-separated target gene list |
| `contrast` | character | Comparison label |

This shared schema enables a unified visualization pattern for both TF and PROGENy results.

#### Diverging Color Palette

All activity plots use the same colorblind-safe diverging scale, defined in `config/color_config.R`:

| Role | Hex | Color |
|------|-----|-------|
| Negative (suppressed) | `#2166AC` | Blue |
| Neutral (zero) | `#F7F7F7` | White |
| Positive (activated) | `#B35806` | Orange |

Apply via `scale_fill_gradient2()` for barplots or `scale_color_manual()` for direction-based coloring:

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

#### Activity Barplot Pattern

The barplot is the primary ranking visualization. The pattern is identical for both TFs and PROGENy:

```r
# 1. Select and order data
plot_data <- results %>%
  arrange(pvalue) %>%         # or arrange(desc(abs(nes)))
  head(N_TOP) %>%
  mutate(pathway_name = factor(pathway_name, levels = pathway_name[order(nes)]))

# 2. Build barplot with diverging fill
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

The volcano plot maps activity score (x) against significance (y), with direction-based coloring and selective labeling:

```r
# 1. Add derived columns
results <- results %>%
  mutate(
    neg_log_padj = -log10(padj),
    is_significant = padj < FDR_CUTOFF,
    label = if_else(is_significant & abs(nes) >= ACTIVITY_CUTOFF, pathway_name, NA_character_)
  )

# 2. Build volcano
p <- ggplot(results, aes(x = nes, y = neg_log_padj, color = direction)) +
  geom_point(aes(alpha = is_significant), size = 2) +
  geom_hline(yintercept = -log10(FDR_CUTOFF), linetype = "dashed", color = "gray50") +
  geom_vline(xintercept = 0, linetype = "solid", color = "gray70") +
  geom_text_repel(aes(label = label), size = 3, max.overlaps = 20, na.rm = TRUE) +
  scale_color_manual(values = c(Up = DIVERGING_COLORS$positive, Down = DIVERGING_COLORS$negative)) +
  scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.3), guide = "none")
```

### Intermediate Usage

#### TF-Specific Visualizations

**Direction Summary (4-Category Bar Chart):**

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
      direction == "Up" & is_significant  ~ DIVERGING_COLORS$positive,
      direction == "Down" & is_significant ~ DIVERGING_COLORS$negative,
      direction == "Up"   ~ adjustcolor(DIVERGING_COLORS$positive, alpha.f = 0.4),
      direction == "Down" ~ adjustcolor(DIVERGING_COLORS$negative, alpha.f = 0.4)
    )
  )
```

**Activity Distribution Histogram:**

Overlaid histograms by direction show the full score distribution:

```r
p_dist <- ggplot(tf_results, aes(x = nes, fill = direction)) +
  geom_histogram(bins = 50, alpha = 0.7, position = "identity") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  scale_fill_manual(values = c(Up = DIVERGING_COLORS$positive, Down = DIVERGING_COLORS$negative))
```

**Top N Selection and Label Threshold:**

TFs are selected by p-value (`arrange(pvalue) %>% head(N_TOP_TFS)` where `N_TOP_TFS = 25`), but volcano labels require both significance AND a minimum activity magnitude:

```r
TF_FDR_CUTOFF <- 0.05
TF_ACTIVITY_CUTOFF <- 2.0  # Minimum |activity score| to label on volcano

label = case_when(
  is_significant & abs(nes) >= TF_ACTIVITY_CUTOFF ~ pathway_name,
  TRUE ~ NA_character_
)
```

#### PROGENy-Specific Visualizations

**All 14 Pathways (No Filtering):**

Unlike TF visualization which selects top N, PROGENy shows all 14 signaling pathways. There is no filtering step:

```r
progeny_ordered <- progeny_results %>%
  mutate(pathway_name = factor(pathway_name, levels = pathway_name[order(nes)]))
```

**Significance Stars:**

PROGENy barplots use significance markers placed at bar tips:

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

Note: PROGENy uses raw `pvalue` for stars (not `padj`), since with only 14 pathways the multiple testing burden is minimal.

**Target Gene Scatter Plot (Weight vs T-Value with Concordance):**

For each top pathway, a scatter plot shows PROGENy gene weights (x) against DE t-statistics (y), colored by concordance:

```r
# Join PROGENy network weights with DE results
df_genes <- pathway_genes %>%
  inner_join(deg_table %>% rownames_to_column("gene") %>% select(gene, t, logFC, adj.P.Val),
             by = c("target" = "gene"))

# Classify concordance
df_genes <- df_genes %>%
  mutate(
    effect = case_when(
      weight > 0 & t > 0 ~ "Concordant (++)",
      weight < 0 & t < 0 ~ "Concordant (--)",
      weight > 0 & t < 0 ~ "Discordant (+-)",
      weight < 0 & t > 0 ~ "Discordant (-+)",
      TRUE ~ "Neutral"
    )
  )

# Concordance color scheme
concordant_colors <- c(
  "Concordant (++)" = "#2ca02c",   # Green
  "Concordant (--)" = "#2ca02c",
  "Discordant (+-)" = "#d62728",   # Red
  "Discordant (-+)" = "#d62728",
  "Neutral"         = "gray50"
)
```

The top 20 genes by absolute t-value are labeled with `ggrepel::geom_text_repel()`. The subtitle reports the pathway activity score and p-value.

#### Output Organization

```
03_results/plots/
  TF/
    tf_activity_barplot.pdf         # Top 25 TFs by p-value
    tf_activity_barplot.png
    tf_activity_volcano.pdf         # Activity vs significance
    tf_activity_volcano.png
    tf_direction_summary.pdf        # 4-category direction chart
    tf_direction_summary.png
    tf_activity_distribution.pdf    # Score histogram by direction
    tf_activity_distribution.png
    README.md                       # Key statistics and interpretation
  PROGENy/
    progeny_activity_barplot.pdf    # All 14 pathways with stars
    progeny_activity_barplot.png
    progeny_activity_volcano.pdf    # Activity vs significance
    progeny_activity_volcano.png
    progeny_targets_*.pdf           # Target gene scatter (top 4 pathways)
    progeny_targets_*.png
    README.md                       # Results table and interpretation
```

Publication standards: 300 DPI PNG, PDF vector, default 10 x 8 inches (barplot 9 x 7 for PROGENy).

#### README Generation

Both scripts generate README.md files containing:
- File manifest with descriptions
- Key statistics (total tested, significant count, up/down breakdown)
- Top upregulated and downregulated entries
- Method description (ULM for TF, MLM for PROGENy)
- Interpretation notes (positive score = targets upregulated in KO vs control)

For PROGENy, the README also includes a markdown table of all 14 pathways with activity scores, p-values, and directions.

### Advanced Usage

#### Color Helper Functions

The `config/color_config.R` file provides reusable helper functions:

```r
# For ggplot2 fill aesthetic (barplots, shape 21 points)
nes_ggplot_fill_scale(limits = c(-3.5, 3.5), name = "NES")

# For ggplot2 color aesthetic (points, lines)
nes_ggplot_scale(limits = c(-3.5, 3.5), name = "NES")

# For ComplexHeatmap
col_fun <- nes_color_scale(limits = c(-3.5, 3.5), n_colors = 3)
```

These helpers centralize the diverging palette so changes propagate automatically to all plots.

#### Publication Standards

All visualization scripts follow these conventions:

| Parameter | Value |
|-----------|-------|
| DPI | 300 (from `PLOT_DPI` in config) |
| Default size | 10 x 8 inches |
| PROGENy barplot | 9 x 7 inches |
| Font base size | 12 pt (`theme_minimal(base_size = 12)`) |
| Output formats | PDF (vector) + PNG (raster), always paired |

#### Integration with Pathway Explorer

Both master activity CSVs feed into the interactive pathway explorer (`3.1.pathway_explorer.py`), where TF and PROGENy entries appear alongside GSEA pathways in the UMAP embedding. The pathway explorer uses:
- Shape coding: diamond for TFs, square for PROGENy (circle for GSEA pathways)
- Same blue-white-orange NES coloring
- Hybrid similarity: TF-Pathway uses Overlap coefficient, PROGENy-Pathway uses Jaccard

No additional visualization code is needed -- the explorer auto-detects entries from the master CSVs.

---

## Verification Checklist

After running TF and PROGENy visualization scripts, confirm:

- [ ] **TF plot files exist:** `tf_activity_barplot.pdf`, `tf_activity_volcano.pdf`, `tf_direction_summary.pdf`, `tf_activity_distribution.pdf` in `plots/TF/`
- [ ] **PROGENy plot files exist:** `progeny_activity_barplot.pdf`, `progeny_activity_volcano.pdf`, `progeny_targets_*.pdf` in `plots/PROGENy/`
- [ ] **PDF/PNG paired:** Every `.pdf` has a corresponding `.png` file
- [ ] **README generated:** Both `plots/TF/README.md` and `plots/PROGENy/README.md` exist and contain correct statistics (total TFs tested, significant count, up/down split)
- [ ] **NES color scale centered:** Barplot fill gradient is symmetric around zero (blue-white-orange)
- [ ] **PROGENy shows all 14:** The barplot includes all 14 signaling pathways, not a filtered subset
- [ ] **Significance markers consistent:** PROGENy uses stars (* p<0.05, ** p<0.01, *** p<0.001); TF uses FDR < 0.05 bold face on axis labels

---

## Common Pitfalls

### Pitfall: ggplot2 4.0+ NA Color Removes Points

- **Symptom:** Shape 21 (filled circle with border) points disappear entirely from volcano or scatter plots. Console warns "removed rows containing missing values."
- **Cause:** ggplot2 4.0+ treats `color = NA` as a missing aesthetic, causing the entire geom to be dropped. In older versions, `color = NA` meant "no border."
- **Fix:** Use `color = "transparent"` or `stroke = 0` instead of `color = NA`:

```r
# WRONG (ggplot2 4.0+)
geom_point(shape = 21, color = NA, fill = "blue")

# CORRECT
geom_point(shape = 21, stroke = 0, fill = "blue")
```

### Pitfall: Missing Master CSV

- **Symptom:** Script exits with "TF results file not found" or "PROGENy results not found."
- **Cause:** The activity inference scripts (`1.7.decoupler_tf_analysis.R`, `1.8.progeny_analysis.R`) have not been run, or the master CSV was not saved to `DIR_TABLES`.
- **Fix:** Run the upstream analysis scripts first. The visualization scripts are Phase 2 and depend on Phase 1 outputs:

```bash
Rscript 02_analysis/1.7.decoupler_tf_analysis.R   # produces master_tf_activities.csv
Rscript 02_analysis/1.8.progeny_analysis.R         # produces master_progeny_activities.csv
```

### Pitfall: Significance Threshold Inconsistency

- **Symptom:** TF volcano labels do not match expected count, or README statistics disagree with plot highlighting.
- **Cause:** The TF significance threshold (`TF_FDR_CUTOFF = 0.05`) is defined locally in the visualization script. If it diverges from the threshold used in the analysis script or the pathway explorer, results appear inconsistent.
- **Fix:** Ensure `TF_FDR_CUTOFF` in `2.3.tf_viz.R` matches the FDR threshold in `1.7.decoupler_tf_analysis.R` and `config/pipeline.yaml`. For PROGENy, raw p-value thresholds (0.05/0.01/0.001) are used for stars since only 14 pathways are tested.

### Pitfall: Direction Column Case Sensitivity

- **Symptom:** `scale_color_manual()` warns about unused values, or points appear with unexpected default colors.
- **Cause:** The `direction` column must contain exactly "Up" or "Down" (title case). Values like "up", "UP", or "upregulated" will not match the manual color mapping.
- **Fix:** Verify the direction column values before plotting:

```r
stopifnot(all(results$direction %in% c("Up", "Down")))
```

If values differ, recode them:

```r
results <- results %>%
  mutate(direction = if_else(nes > 0, "Up", "Down"))
```

### Pitfall: PROGENy Barplot Color Scale Asymmetry

- **Symptom:** The barplot gradient appears skewed -- one direction is more saturated than the other.
- **Cause:** `scale_fill_gradient2()` limits are set asymmetrically (e.g., `limits = c(min(nes) - 1, max(nes) + 1)`) and the data range is not symmetric around zero.
- **Fix:** Set symmetric limits based on the maximum absolute value:

```r
nes_max <- max(abs(progeny_ordered$nes)) + 1
scale_fill_gradient2(..., limits = c(-nes_max, nes_max))
```

---

## Complementary Skills

| When you need... | Use skill | Relationship |
|---|---|---|
| Run TF activity inference (DecoupleR ULM) | `decoupler-tf-activity` | Prerequisite |
| Run PROGENy pathway activity inference | `progeny-pathway-activity` | Prerequisite |
| Visualize GSEA enrichment results (dotplots, running sum) | `bulk-rnaseq-gsea-visualization` | Alternative (for GSEA-type viz) |
| Interactive pathway exploration with TF/PROGENy overlays | `bulk-rnaseq-gsea-visualization` | Extension |

---

## Resources

- **DecoupleR:** https://saezlab.github.io/decoupleR/
- **PROGENy:** https://saezlab.github.io/progeny/
- **CollecTRI (TF regulons):** https://github.com/saezlab/CollecTRI
- **ggrepel:** https://ggrepel.slowkow.com/
- **Color config:** `02_analysis/config/color_config.R`
- **Pipeline config:** `02_analysis/config/pipeline.yaml`
- **TF viz script:** `02_analysis/2.3.tf_viz.R`
- **PROGENy viz script:** `02_analysis/2.4.progeny_viz.R`
