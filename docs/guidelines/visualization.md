# Visualization Standards

**Module:** `visualization.md`
**Purpose:** Colors, themes, publication standards, and plot quality

---

## 1. Colorblind-Safe Palettes

### 1.1 Okabe-Ito Palette (Default for Categories)

```yaml
# In analysis_config.yaml
colors:
  okabe_ito:
    orange: "#E69F00"
    sky_blue: "#56B4E9"
    bluish_green: "#009E73"
    yellow: "#F0E442"
    blue: "#0072B2"
    vermillion: "#D55E00"
    reddish_purple: "#CC79A7"
    black: "#000000"
```

### 1.2 Database Colors

```yaml
colors:
  databases:
    Hallmark: "#E69F00"       # Orange
    KEGG: "#56B4E9"           # Sky Blue
    Reactome: "#009E73"       # Bluish Green
    WikiPathways: "#F0E442"   # Yellow
    GO_BP: "#0072B2"          # Blue
    GO_MF: "#D55E00"          # Vermillion
    GO_CC: "#CC79A7"          # Reddish Purple
    MitoCarta: "#117733"      # Dark Green
    Custom: "#999999"         # Gray
```

### 1.3 Diverging Scale (NES, logFC)

```yaml
colors:
  diverging:
    down: "#2166AC"           # Blue (negative)
    neutral: "#F7F7F7"        # White/light gray
    up: "#B35806"             # Orange (positive)
```

### 1.4 Sequential Scale (p-values, counts)

```yaml
colors:
  sequential:
    low: "#FEE0D2"            # Light
    mid: "#FC9272"            # Medium
    high: "#DE2D26"           # Dark (significant)
```

---

## 2. R Implementation

### 2.1 Color Configuration

```r
# In 02_analysis/config/color_config.R

# Load from YAML
colors_config <- yaml::read_yaml("analysis_config.yaml")$colors

# Database colors
DB_COLORS <- unlist(colors_config$databases)

# Diverging scale
DIVERGING_COLORS <- c(
  colors_config$diverging$down,
  colors_config$diverging$neutral,
  colors_config$diverging$up
)

# ggplot2 scale functions
scale_fill_database <- function() {
  scale_fill_manual(values = DB_COLORS)
}

scale_color_database <- function() {
  scale_color_manual(values = DB_COLORS)
}

scale_fill_diverging <- function(limits = NULL, midpoint = 0) {
  scale_fill_gradient2(
    low = DIVERGING_COLORS[1],
    mid = DIVERGING_COLORS[2],
    high = DIVERGING_COLORS[3],
    midpoint = midpoint,
    limits = limits
  )
}
```

### 2.2 Custom Theme

```r
# Minimal theme for publication
theme_publication <- function(base_size = 12) {
  theme_minimal(base_size = base_size) +
    theme(
      # Clean panel
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),

      # Axis
      axis.line = element_line(color = "black", linewidth = 0.5),
      axis.ticks = element_line(color = "black", linewidth = 0.5),
      axis.text = element_text(color = "black"),
      axis.title = element_text(face = "bold"),

      # Legend
      legend.background = element_blank(),
      legend.key = element_blank(),

      # Title
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),

      # Margins
      plot.margin = margin(10, 10, 10, 10)
    )
}
```

---

## 3. Plot Types and Standards

### 3.1 Volcano Plot

```r
create_volcano_plot <- function(de_results,
                                fdr_cutoff = 0.05,
                                logfc_cutoff = 1,
                                label_top = 10) {
  # Prepare data
  plot_data <- de_results %>%
    mutate(
      significance = case_when(
        adj.P.Val >= fdr_cutoff ~ "NS",
        abs(logFC) < logfc_cutoff ~ "Sig only",
        logFC > logfc_cutoff ~ "Up",
        logFC < -logfc_cutoff ~ "Down"
      ),
      neg_log_p = -log10(P.Value)
    )

  # Identify top genes for labeling
  top_genes <- plot_data %>%
    filter(significance %in% c("Up", "Down")) %>%
    slice_max(abs(logFC), n = label_top)

  # Plot
  ggplot(plot_data, aes(x = logFC, y = neg_log_p)) +
    geom_point(aes(color = significance), alpha = 0.6, size = 1.5) +
    geom_vline(xintercept = c(-logfc_cutoff, logfc_cutoff),
               linetype = "dashed", color = "gray50") +
    geom_hline(yintercept = -log10(fdr_cutoff),
               linetype = "dashed", color = "gray50") +
    ggrepel::geom_text_repel(
      data = top_genes,
      aes(label = gene_symbol),
      size = 3, max.overlaps = 20
    ) +
    scale_color_manual(
      values = c("Up" = "#B35806", "Down" = "#2166AC",
                 "Sig only" = "gray50", "NS" = "gray80")
    ) +
    labs(
      x = expression(log[2]~"Fold Change"),
      y = expression(-log[10]~"P-value"),
      color = "Significance"
    ) +
    theme_publication()
}
```

### 3.2 GSEA Dotplot

```r
create_gsea_dotplot <- function(gsea_data,
                                n_pathways = 20,
                                fdr_cutoff = 0.05) {
  # Filter and prepare
  plot_data <- gsea_data %>%
    filter(padj < fdr_cutoff) %>%
    slice_max(abs(nes), n = n_pathways) %>%
    mutate(pathway_name = fct_reorder(pathway_name, nes))

  # Plot
  ggplot(plot_data, aes(x = nes, y = pathway_name)) +
    geom_point(aes(size = set_size, color = padj)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    scale_color_gradient(
      low = "#DE2D26",    # Significant (low p)
      high = "#FEE0D2",   # Less significant
      trans = "log10",
      guide = guide_colorbar(reverse = TRUE)
    ) +
    scale_size_continuous(range = c(2, 8)) +
    labs(
      x = "Normalized Enrichment Score (NES)",
      y = NULL,
      color = "Adjusted\np-value",
      size = "Gene Set\nSize"
    ) +
    theme_publication() +
    theme(axis.text.y = element_text(size = 10))
}
```

### 3.3 Heatmap

```r
create_expression_heatmap <- function(expr_matrix,
                                      sample_annotation,
                                      genes = NULL,
                                      cluster_rows = TRUE,
                                      cluster_cols = TRUE) {
  library(ComplexHeatmap)

  # Subset genes if specified
  if (!is.null(genes)) {
    expr_matrix <- expr_matrix[rownames(expr_matrix) %in% genes, ]
  }

  # Scale by row (z-score)
  scaled_matrix <- t(scale(t(expr_matrix)))

  # Annotation
  col_annotation <- HeatmapAnnotation(
    Group = sample_annotation$Group,
    col = list(Group = c("Control" = "#56B4E9", "Treatment" = "#E69F00"))
  )

  # Heatmap
  Heatmap(
    scaled_matrix,
    name = "Z-score",
    col = colorRamp2(c(-2, 0, 2), c("#2166AC", "#F7F7F7", "#B35806")),
    top_annotation = col_annotation,
    cluster_rows = cluster_rows,
    cluster_columns = cluster_cols,
    show_row_names = nrow(scaled_matrix) <= 50,
    row_names_gp = gpar(fontsize = 8)
  )
}
```

---

## 4. Publication Standards

### 4.1 Resolution and Size

```yaml
# In analysis_config.yaml
visualization:
  dpi: 300               # Print quality
  width: 10              # inches (for most plots)
  height: 8              # inches
  width_narrow: 6        # Single column
  width_wide: 14         # Two columns

  # Journal-specific
  nature_width: 89       # mm (single column)
  nature_width_double: 183  # mm (double column)
```

### 4.2 Saving Plots

```r
# Save both PDF and PNG
save_publication_plot <- function(plot, filename, width = 10, height = 8) {
  # Create directory if needed
  dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)

  # PDF for publication
  ggsave(
    paste0(filename, ".pdf"),
    plot = plot,
    width = width,
    height = height,
    device = cairo_pdf
  )

  # PNG for quick viewing
  ggsave(
    paste0(filename, ".png"),
    plot = plot,
    width = width,
    height = height,
    dpi = 300
  )

  message("Saved: ", filename, ".{pdf,png}")
}
```

### 4.3 Figure Panels

```r
# Combining plots for figures
library(patchwork)

# Side by side
combined <- plot_a + plot_b

# Stacked
combined <- plot_a / plot_b

# With labels
combined <- (plot_a + plot_b) / (plot_c + plot_d) +
  plot_annotation(tag_levels = 'A')

# Save
save_publication_plot(combined, "03_results/plots/Publication/Figure1",
                      width = 14, height = 10)
```

---

## 5. Python Visualization

### 5.1 Matplotlib Configuration

```python
import matplotlib.pyplot as plt
import seaborn as sns

# Set style
plt.style.use('seaborn-v0_8-whitegrid')
plt.rcParams.update({
    'figure.figsize': (10, 8),
    'figure.dpi': 100,
    'savefig.dpi': 300,
    'font.size': 12,
    'axes.labelsize': 14,
    'axes.titlesize': 14,
    'legend.fontsize': 11,
    'xtick.labelsize': 11,
    'ytick.labelsize': 11,
})

# Color palette
COLORS = {
    'diverging': ['#2166AC', '#F7F7F7', '#B35806'],
    'databases': {
        'Hallmark': '#E69F00',
        'KEGG': '#56B4E9',
        'Reactome': '#009E73',
    }
}
```

### 5.2 Plotly Interactive

```python
import plotly.express as px
import plotly.graph_objects as go

def create_interactive_volcano(de_data):
    """Create interactive volcano plot."""
    fig = px.scatter(
        de_data,
        x='logFC',
        y=-np.log10(de_data['P.Value']),
        color='direction',
        hover_data=['gene_symbol', 'adj.P.Val'],
        color_discrete_map={
            'Up': '#B35806',
            'Down': '#2166AC',
            'NS': '#CCCCCC'
        }
    )

    fig.update_layout(
        xaxis_title='log₂ Fold Change',
        yaxis_title='-log₁₀ P-value',
        template='plotly_white'
    )

    return fig
```

---

## 6. Quality Checklist

### 6.1 Before Saving

- [ ] Colors are colorblind-safe
- [ ] Font sizes are readable at target size
- [ ] Axis labels include units
- [ ] Legend is clear and complete
- [ ] No overlapping text
- [ ] Significance thresholds are marked

### 6.2 File Requirements

- [ ] PDF for vector graphics (publication)
- [ ] PNG for raster (presentations, web)
- [ ] DPI >= 300 for print
- [ ] File size reasonable (<10MB per plot)

### 6.3 Documentation

- [ ] README in plot directory
- [ ] Script path noted
- [ ] Data source documented
- [ ] Statistical methods noted
