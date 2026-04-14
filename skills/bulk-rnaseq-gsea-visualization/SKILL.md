---
name: bulk-rnaseq-gsea-visualization
description: "RNAseq-toolkit GSEA visualization and interactive pathway exploration -- publication-quality dotplots, barplots, running sum plots, and interactive HTML dashboards from master GSEA tables. Use when generating GSEA figures from completed GSEA analysis, creating cross-database summary plots, or building interactive pathway explorers. For running GSEA itself use bulk-rnaseq-gsea-msigdb or bulk-rnaseq-gsea-custom-db. For master table assembly use bulk-rnaseq-gsea-master-tables."
license: MIT
metadata:
  skill-author: SciAgent-toolkit
  last-reviewed: 2026-04-14
  version: 1.0.0
  upstream-docs: https://github.com/tony-zhelonkin/RNAseq-toolkit
  category: workflow
  tier: standard
  tags:
    - gsea
    - visualization
    - dotplot
    - running-sum
    - pathway-explorer
    - ggplot2
    - plotly
    - publication-figures
    - interactive-dashboard
    - bulk-rnaseq
  complementary-skills:
    - bulk-rnaseq-gsea-msigdb
    - bulk-rnaseq-gsea-custom-db
    - bulk-rnaseq-gsea-master-tables
  contraindications:
    - "Do not use for running GSEA or preparing gene sets. Use bulk-rnaseq-gsea-msigdb or bulk-rnaseq-gsea-custom-db instead."
    - "Do not use for master table assembly. Use bulk-rnaseq-gsea-master-tables instead."
    - "Do not use for GATOM network visualization. Use gatom-metabolomic-predictions instead."
---

# Bulk RNA-seq GSEA Visualization and Interactive Pathway Exploration

## Overview

This skill covers publication-quality GSEA visualization (R) and interactive pathway exploration (Python) from completed GSEA analyses. It bridges the gap between cached GSEA results (RDS checkpoints) and final figures by using a two-track system: R scripts (2.x.*.R) consume gseaResult objects directly for static plots, while a Python dashboard (3.x.*.py) consumes master CSV tables for interactive exploration. The master table CSV is the single bridge between R and Python.

**When to use this skill:**
- Creating per-database GSEA visualizations (dotplots, barplots, running sum plots)
- Building cross-database pooled summary figures
- Generating an interactive pathway explorer HTML dashboard
- Adding a new database's results to the existing visualization pipeline
- Troubleshooting ggplot2 or enrichplot rendering issues with GSEA objects

**When NOT to use this skill:**
- Running GSEA analysis itself -> use `bulk-rnaseq-gsea-msigdb`
- Building custom gene sets -> use `bulk-rnaseq-gsea-custom-db`
- Assembling master_gsea_table.csv -> use `bulk-rnaseq-gsea-master-tables`
- GATOM metabolic network visualization -> use `gatom-metabolomic-predictions`

---

## Decision Tree

```
Need GSEA figures?
|
+-- Which biological question?
|   |
|   +-- "Which pathways are significant?"
|   |   +-- Few databases (1-3)     -->  gsea_dotplot() per database
|   |   +-- Many databases (4+)     -->  Cross-database pooled dotplot
|   |
|   +-- "What direction are pathways enriched?"
|   |   +-- Compact overview         -->  gsea_barplot() (NES bars with direction coloring)
|   |   +-- Separated up/down        -->  gsea_dotplot_facet() or filterBy="NES_positive"/"NES_negative"
|   |
|   +-- "How does enrichment accumulate across the ranked list?"
|   |   +-- Top 5 pathways together  -->  gsea_running_sum_plot() stacked
|   |   +-- Single pathway detail    -->  gsea_running_sum_plot() individual
|   |
|   +-- "How do pathways relate to each other by gene overlap?"
|   |   -->  Pathway explorer (Python interactive HTML)
|   |
|   +-- "What is the expression pattern of pathway genes?"
|       -->  gsea_heatmap() (pathway x sample z-score)
|
+-- What format?
    +-- Publication figure (PDF/PNG)  -->  R scripts (2.x.*.R)
    +-- Interactive exploration       -->  Python (3.x.*.py) -> HTML
```

---

## Quick Start

### R Static Visualizations (from RDS checkpoints)

```r
# Source configuration and toolkit
source("02_analysis/config/config.R")
source("02_analysis/config/color_config.R")
source_toolkit(verbose = TRUE)

library(clusterProfiler)
library(ggplot2)

# Load GSEA results
all_gsea <- readRDS(file.path(DIR_CHECKPOINTS, "1.1_all_gsea_results.rds"))
hallmark <- all_gsea[["H"]]

# Basic dotplot: top 20 by FDR, outline FDR < 0.05
p <- gsea_dotplot(
  hallmark,
  filterBy = "p.adjust",
  showCategory = 20,
  padj_cutoff = 0.05,
  title = "Hallmark: IL2RAKO vs NTC"
)
ggsave("hallmark_dotplot.pdf", p, width = 12, height = 10)
```

### Python Interactive Dashboard (from master CSV)

```bash
cd 02_analysis
python 3.1.pathway_explorer.py
# Output: 03_results/interactive/pathway_explorer.html
```

**Verify it worked:**
- R plots: Check that PDF files are non-empty and contain visible pathways
- Pathway explorer: Open `pathway_explorer.html` in a browser; verify scatter plot loads, database filters work, and clicking a pathway shows its gene table and running sum

---

## Progressive Depth

### Basic Usage

#### Per-Database Visualization (7-Step Pattern)

Every database follows an identical 7-step procedure in `2.2.gsea_viz.R`:

```r
# For each database (e.g., Hallmark, KEGG, Reactome, MitoPathways, TransportDB):

# Step 1: Combined dotplot (all significant, up + down)
p_dotplot <- gsea_dotplot(
  gsea_res,
  filterBy = "p.adjust",
  showCategory = min(20, n_sig),
  padj_cutoff = 0.05,
  title = sprintf("%s: %s", db_display, contrast_name)
)
ggsave(file.path(db_dir, paste0(prefix, "_dotplot.pdf")), p_dotplot, width = 12, height = 10)

# Step 2: Upregulated-only dotplot
p_up <- gsea_dotplot(gsea_res, filterBy = "NES_positive",
                     showCategory = min(20, n_up), padj_cutoff = 0.05)
ggsave(file.path(db_dir, paste0(prefix, "_up_dot.pdf")), p_up, width = 12, height = 10)

# Step 3: Downregulated-only dotplot
p_down <- gsea_dotplot(gsea_res, filterBy = "NES_negative",
                       showCategory = min(20, n_down), padj_cutoff = 0.05)
ggsave(file.path(db_dir, paste0(prefix, "_down_dot.pdf")), p_down, width = 12, height = 10)

# Step 4: NES barplot
p_bar <- gsea_barplot(gsea_res, top_n = min(20, n_sig), padj_cutoff = 0.05,
                      title = sprintf("%s NES", db_display))
ggsave(file.path(db_dir, paste0(prefix, "_nes_bar.pdf")), p_bar, width = 10, height = 8)

# Step 5: Stacked running sum (top 5 by |NES|)
top_ids <- gsea_res@result %>%
  filter(p.adjust < 0.05) %>%
  arrange(desc(abs(NES))) %>%
  head(5) %>%
  pull(ID)
p_running <- gsea_running_sum_plot(gsea_res, gene_set_ids = top_ids,
                                   palette = RUNNING_SUM_PALETTE[seq_along(top_ids)])
ggsave(file.path(db_dir, paste0(prefix, "_running_sum.pdf")), p_running, width = 12, height = 10)

# Step 6: Individual running sums (top 10)
for (pathway in top10_ids) {
  p_single <- gsea_running_sum_plot(gsea_res, gene_set_ids = pathway,
                                    palette = RUNNING_SUM_PALETTE[1])
  ggsave(file.path(running_dir, paste0("running_", safe_filename(pathway), ".pdf")),
         p_single, width = 10, height = 8)
}

# Step 7: Text summary
save_gsea_summary(gsea_res, file.path(db_dir, paste0(prefix, "_results.txt")),
                  contrast_name, db_display, 0.05)
```

#### Key Concept: Selection vs Highlighting in gsea_dotplot()

The dotplot separates **which pathways to show** from **which to highlight**:

```r
gsea_dotplot(
  gsea_obj,
  # SELECTION: which pathways to display
  filterBy = "NES",           # Sort criterion: "NES", "p.adjust", "NES_positive", "NES_negative"
  showCategory = 20,          # How many to show

  # HIGHLIGHTING: which get a black outline
  padj_cutoff = 0.10,         # FDR threshold for outline
  highlight_sig = TRUE,       # Enable/disable outlines
  highlight_threshold = 0.01, # Optional stricter outline threshold

  # APPEARANCE
  use_gradient = TRUE         # Continuous NES color gradient
)
```

This "show all, highlight significant" pattern ensures top-ranked pathways are always visible even if not all pass FDR, while significant ones are marked by a black border.

### Intermediate Usage

#### Cross-Database Pooled Dotplot

Collects top significant pathways from all databases into a single figure:

```r
# Collect top 10 significant from each database
all_pathways <- list()
for (db_name in names(all_gsea)) {
  gsea_res <- all_gsea[[db_name]]
  top_paths <- gsea_res@result %>%
    filter(p.adjust < 0.05) %>%
    arrange(desc(abs(NES))) %>%
    head(10) %>%
    mutate(
      GeneRatio = sapply(seq_len(n()), function(i) {
        length(strsplit(core_enrichment[i], "/")[[1]]) / setSize[i]
      }),
      negLogPadj = -log10(p.adjust),
      Database = db_display,
      PathwayClean = format_pathway_name(Description)
    )
  all_pathways[[db_name]] <- top_paths
}

combined_df <- bind_rows(all_pathways) %>%
  mutate(
    PathwayID = paste0(Database, ": ", PathwayClean),
    is_highly_sig = p.adjust < 0.01
  )

# Dual-layer geom_point for conditional black outline
pooled_plot <- ggplot(combined_df, aes(x = GeneRatio, y = PathwayID)) +
  # Layer 1: Base points (no outline)
  geom_point(aes(size = negLogPadj, fill = NES), shape = 21, stroke = 0) +
  # Layer 2: Outline on highly significant only
  geom_point(data = filter(combined_df, is_highly_sig),
             aes(size = negLogPadj, fill = NES),
             shape = 21, stroke = 0.7, color = "black") +
  nes_ggplot_fill_scale(limits = c(-nes_max, nes_max), name = "NES") +
  scale_size_continuous(name = expression(-log[10](FDR)), range = c(2, 8))
```

Three variants are produced:
- `pooled_{contrast}.pdf` -- all databases
- `focused_top5_{contrast}.pdf` -- key databases, top 5 per db, faceted
- `focused_top10_{contrast}.pdf` -- key databases, top 10 per db, faceted

#### Color System

All colors are defined in `config/pipeline.yaml` (shared YAML, single source of truth) and loaded by both R (`config/color_config.R`) and Python (`pathway_explorer/config.py`).

**NES Diverging Scale:**

| Role | Hex | Color |
|------|-----|-------|
| Negative (downregulated) | `#2166AC` | Blue |
| Neutral | `#F7F7F7` | White |
| Positive (upregulated) | `#B35806` | Orange |

R helpers: `nes_ggplot_scale()`, `nes_ggplot_fill_scale()`, `nes_color_scale()`

**Database Colors (Okabe-Ito Colorblind-Safe, 13 databases):**

| Database | Hex |
|----------|-----|
| Hallmark | `#E69F00` |
| KEGG | `#56B4E9` |
| Reactome | `#009E73` |
| WikiPathways | `#F0E442` |
| GO_BP | `#0072B2` |
| GO_MF | `#D55E00` |
| GO_CC | `#CC79A7` |
| MitoPathways | `#332288` |
| MitoXplorer | `#999999` |
| CollecTRI | `#882255` |
| PROGENy | `#44AA99` |
| TransportDB | `#DDCC77` |
| GATOM | `#117733` |

**Running Sum Palette** (9 colors for multi-pathway overlay):

```r
RUNNING_SUM_PALETTE <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3",
                         "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999")
```

#### Output Organization

```
03_results/plots/GSEA/
  {Contrast}/
    {Database}/
      {contrast}_{database}_dotplot.pdf
      {contrast}_{database}_up_dot.pdf
      {contrast}_{database}_down_dot.pdf
      {contrast}_{database}_nes_bar.pdf
      {contrast}_{database}_running_sum.pdf
      {contrast}_{database}_results.txt
      running_sum/
        running_{PATHWAY_ID}.pdf
    cross_database_pooled/
      pooled_{contrast}.pdf
      focused_top5_{contrast}.pdf
      focused_top10_{contrast}.pdf
      pooled_pathways_data.csv

03_results/interactive/
  pathway_explorer.html
```

Publication standards: 300 DPI PNG, PDF vector, default 10 x 8 inches (dotplots 12 x 10, running sum stacked 12 x 10).

### Advanced Usage

#### Pathway Explorer Architecture

The interactive pathway explorer (`3.1.pathway_explorer.py`) wraps a modular Python package:

```
pathway_explorer/
  config.py        # Paths, YAML loading, schema validation
  data_loader.py   # CSV loading, score standardization
  similarity.py    # Jaccard/Overlap similarity matrices
  embedding.py     # UMAP dimensionality reduction
  html_generator.py# Self-contained HTML dashboard
```

**Data flow:**
1. Load `master_gsea_table.csv`, `master_tf_activities.csv`, `master_progeny_activities.csv`, `master_de_table.csv`
2. Validate schemas against `pipeline.yaml` definitions
3. Compute pairwise gene-set similarity (hybrid strategy):
   - Same-type (Pathway-Pathway, TF-TF, PROGENy-PROGENy): **Jaccard** index
   - TF cross-type (TF-Pathway, TF-PROGENy): **Overlap** coefficient
   - PROGENy-Pathway: **Jaccard** (both are pathway-like)
4. UMAP embedding (metric=precomputed, n_neighbors=15, min_dist=0.1)
5. Generate self-contained HTML (~7.8 MB) with embedded Plotly.js

**Interactive features:**
- 2D UMAP scatter colored by NES (blue-white-orange diverging)
- Size = -log10(FDR), Shape = entity type (circle=Pathway, diamond=TF, square=PROGENy)
- Click: edge overlay for connected pathways, running sum enrichment plot (client-side), gene table
- Sidebar: database checkboxes, FDR slider, NES threshold, text search

#### Adding a New Database (7-Step Guide)

**Step 1: Create GSEA checkpoint**
```r
mycustomdb_gsea <- run_gsea(DE_results = de_table, rank_metric = "t",
                             custom_gene_sets = my_gene_sets)
saveRDS(mycustomdb_gsea, file.path(DIR_CHECKPOINTS, "1.X_gsea_mycustomdb.rds"))
```

**Step 2: Add to master_gsea_table.csv**
```r
new_rows <- normalize_gsea_results(mycustomdb_gsea, database = "MyCustomDB",
                                    contrast = contrast_name)
master <- bind_rows(read_csv("03_results/tables/master_gsea_table.csv"), new_rows)
write_csv(master, "03_results/tables/master_gsea_table.csv")
```

**Step 3: Add processing block to 2.2.gsea_viz.R** (follow TransportDB pattern)

**Step 4: Register display name and color**
- `2.2.gsea_viz.R`: Add to `DB_DISPLAY_NAMES`
- `config/color_config.R` and `config/pipeline.yaml`: Add to `DATABASE_COLORS`
- `pathway_explorer/config.py`: Add to `DB_COLORS`

**Step 5: Add to cross-database pooled plot** (V7 section of 2.2.gsea_viz.R)

**Step 6: Pathway explorer auto-detects** from `master_gsea_table.csv` -- no code change needed if Step 2 is done. Unregistered databases get default gray.

**Step 7: (Optional) Update `pipeline.yaml` schema** if new columns are introduced.

---

## Verification Checklist

After running visualizations, confirm:

- [ ] **Per-database plots exist:** Every database with significant results has `_dotplot.pdf`, `_nes_bar.pdf`, `_running_sum.pdf` in its directory
- [ ] **Running sum plots render correctly:** Three-panel layout (enrichment score curve, tick marks, ranked metric) with colored gene set lines
- [ ] **Pooled dotplot includes all databases:** Check `pooled_pathways_data.csv` has rows from every database with significant results
- [ ] **NES color scale is symmetric:** Blue-white-orange centered at zero (verify limits are +/- max |NES|)
- [ ] **Pathway explorer loads:** HTML opens in browser, scatter plot renders, database filters toggle visibility
- [ ] **Schema validation passes:** Pathway explorer reports no schema validation warnings on startup
- [ ] **No NA color warnings:** ggplot2 does not warn about "removed rows containing missing values" for shape 21 points

---

## Common Pitfalls

### Pitfall: Description/ID Mismatch in Custom DB Running Sum Plots

- **Symptom:** Running sum plots for custom databases (MitoPathways, MitoXplorer, TransportDB) have wrong colors, missing labels, or fail silently. `enrichplot::gseaplot2()` shows default colors instead of the specified palette.
- **Cause:** `enrichplot::gseaplot2()` internally uses `Description` for color mapping, but `gsea_running_sum_plot()` names palette entries by `ID`. For MSigDB databases, Description and ID match. For custom databases, they diverge (e.g., ID=`MITOPATHWAYS_OXPHOS`, Description=`Oxidative Phosphorylation`).
- **Fix:** Create a plotting copy with Description set to ID, pass original descriptions via the `labels` parameter:

```r
# 1. Save original descriptions
original_desc <- gsea_obj@result$Description
names(original_desc) <- gsea_obj@result$ID

# 2. Set Description = ID for internal color mapping
plot_obj <- gsea_obj
plot_obj@result$Description <- plot_obj@result$ID

# 3. Pass original descriptions via labels parameter
p <- gsea_running_sum_plot(
  plot_obj,
  gene_set_ids = top_ids,
  palette = RUNNING_SUM_PALETTE[seq_along(top_ids)],
  labels = original_desc[top_ids]
)
```

### Pitfall: ggplot2 4.0+ NA Color Causes Point Removal

- **Symptom:** Shape 21 (filled circle with border) points disappear entirely. Console shows "removed rows containing missing values" warning.
- **Cause:** ggplot2 4.0+ treats `color = NA` as a missing aesthetic, which causes the entire point to be removed. In older ggplot2, `color = NA` meant "no border".
- **Fix:** Use `color = "transparent"` instead of `color = NA` for shape 21 points:

```r
# WRONG (ggplot2 4.0+): points are removed
geom_point(shape = 21, color = NA, fill = "blue")

# CORRECT: transparent border
geom_point(shape = 21, color = "transparent", fill = "blue")

# ALSO CORRECT: zero stroke width (avoids the issue entirely)
geom_point(shape = 21, stroke = 0, fill = "blue")
```

The `2.2.gsea_viz.R` script uses `stroke = 0` on base-layer points to avoid this issue.

### Pitfall: Pathway Explorer Schema Validation Failure

- **Symptom:** `pathway_explorer.py` exits with a schema validation error listing missing columns.
- **Cause:** The master CSV table is missing required columns defined in `config/pipeline.yaml` under `schemas.master_gsea_table.required_columns` (e.g., `pathway_id`, `pathway_name`, `database`, `nes`, `pvalue`, `padj`, `core_enrichment`).
- **Fix:** Regenerate the master table using `1.5.create_master_tables.R` or verify that `normalize_gsea_results()` was used correctly when adding new databases.

### Pitfall: Hybrid Similarity Confusion (Jaccard vs Overlap)

- **Symptom:** In the pathway explorer, transcription factors cluster too tightly with pathways, or PROGENy pathways appear disconnected from related GSEA pathways.
- **Cause:** The similarity computation uses different metrics for different entity-pair types. TF-Pathway uses Overlap coefficient (sensitive to small TF target sets overlapping large pathways), while PROGENy-Pathway uses Jaccard (treats both as pathway-like).
- **Fix:** This is by design. If you need to adjust, modify `pathway_explorer/similarity.py`. The key decision: Overlap coefficient inflates similarity when one set is much smaller than the other (common for TFs), while Jaccard penalizes size asymmetry.

### Pitfall: Cross-Database Pooled Plot Has Duplicate Factor Levels

- **Symptom:** Error: `duplicated levels in factors are deprecated` or pathways from different databases with the same name collapse into one row.
- **Cause:** Multiple databases can contain pathways with identical names. The PathwayID must include the database prefix.
- **Fix:** Prefix pathway names with database: `PathwayID = paste0(Database, ": ", PathwayClean)` and use `factor(PathwayID, levels = unique(PathwayID))`.

---

## Complementary Skills

| When you need... | Use skill | Relationship |
|---|---|---|
| Run GSEA on MSigDB collections | `bulk-rnaseq-gsea-msigdb` | Prerequisite |
| Run GSEA on custom gene sets (MitoCarta, TransportDB) | `bulk-rnaseq-gsea-custom-db` | Prerequisite |
| Build master CSV tables from checkpoints | `bulk-rnaseq-gsea-master-tables` | Prerequisite |
| Visualize GATOM metabolic networks | `gatom-metabolomic-predictions` | Alternative (for network viz) |

---

## Resources

- **RNAseq-toolkit:** `01_scripts/RNAseq-toolkit/scripts/GSEA/GSEA_plotting/`
- **Visualization workflow docs:** `01_scripts/RNAseq-toolkit/docs/GSEA-workflow/04-output-artifacts-and-visualization.md`
- **Core pipeline docs:** `01_scripts/RNAseq-toolkit/docs/GSEA-workflow/01-core-pipeline-and-toolkit.md`
- **clusterProfiler:** https://bioconductor.org/packages/clusterProfiler/
- **enrichplot:** https://bioconductor.org/packages/enrichplot/
- **Plotly.js:** https://plotly.com/javascript/
- **UMAP:** https://umap-learn.readthedocs.io/
