# Code Style

**Module:** `code_style.md`
**Purpose:** R and Python coding conventions, function patterns, and documentation

---

## 1. R Style Guide

### 1.1 Naming Conventions

```r
# Variables: snake_case
de_results <- topTable(fit, coef = contrast_name)
gsea_hallmark <- run_gsea(de_results, category = "H")
n_significant <- sum(results$adj.P.Val < 0.05)

# Functions: snake_case (verbs preferred)
run_gsea <- function(...) { }
create_volcano_plot <- function(...) { }
load_or_compute <- function(...) { }

# Constants: SCREAMING_SNAKE_CASE
DE_FDR_CUTOFF <- 0.05
GSEA_NPERM <- 100000
DIR_CHECKPOINTS <- "03_results/checkpoints"

# Avoid:
# - camelCase (Java/JavaScript style)
# - dot.names (old R style, conflicts with S3)
```

### 1.2 Function Documentation (roxygen2)

```r
#' Run GSEA analysis on DE results
#'
#' @description
#' Performs Gene Set Enrichment Analysis using clusterProfiler or fgsea
#' on differential expression results.
#'
#' @param de_results Data frame with columns: gene (rownames), t, logFC, P.Value
#' @param rank_metric Column to use for gene ranking. Default: "t"
#' @param species Species for msigdbr. Default: "Mus musculus"
#' @param category MSigDB category (H, C2, C5, etc.)
#' @param subcategory MSigDB subcategory (e.g., "CP:KEGG")
#' @param nperm Number of permutations. Default: 100000
#' @param seed Random seed for reproducibility. Default: 123
#'
#' @return gseaResult object from clusterProfiler
#'
#' @examples
#' gsea_res <- run_gsea(de_results, category = "H")
#'
#' @export
run_gsea <- function(de_results,
                     rank_metric = "t",
                     species = "Mus musculus",
                     category = "H",
                     subcategory = "",
                     nperm = 100000,
                     seed = 123) {
  # Implementation
}
```

### 1.3 Progress Messages

```r
# Use message() for progress, not print()
message("[STEP 1] Loading data...")
message("[STEP 2] Running GSEA for ", db_name, "...")
message("[DONE] Saved ", n_files, " plots to ", output_dir)

# Structured format
message("=================================================================")
message("PROJECT-ID: Analysis Pipeline")
message("=================================================================")

# Warning for recoverable issues
warning("Gene ", gene_name, " not found, skipping...")

# Stop for fatal errors
stop("Required file not found: ", filepath)
```

### 1.4 Tidyverse Style

```r
# Pipes for data transformation
results <- de_data %>%
  filter(adj.P.Val < 0.05) %>%
  arrange(desc(abs(logFC))) %>%
  select(gene_symbol, logFC, adj.P.Val) %>%
  head(20)

# Not this:
results <- head(select(arrange(filter(de_data, adj.P.Val < 0.05),
                                desc(abs(logFC))),
                       gene_symbol, logFC, adj.P.Val), 20)

# Line breaks for readability
long_result <- data %>%
  filter(
    condition1 == TRUE,
    condition2 > threshold,
    category %in% valid_categories
  ) %>%
  mutate(
    new_col1 = complex_calculation(col1),
    new_col2 = another_calculation(col2, col3)
  )
```

### 1.5 Function Length

```r
# Functions should be < 50 lines
# If longer, extract helper functions

# BAD: Monolithic function
run_full_analysis <- function(data) {
  # 200 lines of code...
}

# GOOD: Modular functions
run_full_analysis <- function(data) {
  data <- preprocess_data(data)
  data <- filter_data(data)
  results <- run_de_analysis(data)
  results <- annotate_results(results)
  return(results)
}

preprocess_data <- function(data) { ... }
filter_data <- function(data) { ... }
run_de_analysis <- function(data) { ... }
annotate_results <- function(results) { ... }
```

---

## 2. Python Style Guide

### 2.1 Naming Conventions

```python
# Variables: snake_case
de_results = load_de_data(filepath)
gsea_hallmark = run_gsea(de_results)
n_significant = sum(results['adj_p_val'] < 0.05)

# Functions: snake_case
def run_gsea(de_results: pd.DataFrame) -> pd.DataFrame:
    pass

# Classes: PascalCase
class PathwayExplorer:
    pass

# Constants: SCREAMING_SNAKE_CASE
DE_FDR_CUTOFF = 0.05
GSEA_NPERM = 100000
DIR_CHECKPOINTS = Path("03_results/checkpoints")

# Private: leading underscore
def _internal_helper():
    pass
```

### 2.2 Type Hints

```python
from typing import Dict, List, Optional, Union
from pathlib import Path
import pandas as pd

def load_gsea_data(
    data_dir: Path,
    database: str = "all",
    fdr_cutoff: float = 0.05
) -> pd.DataFrame:
    """Load and filter GSEA results.

    Args:
        data_dir: Path to tables directory
        database: Filter by database name, or "all"
        fdr_cutoff: FDR significance threshold

    Returns:
        DataFrame with filtered GSEA results
    """
    df = pd.read_csv(data_dir / "master_gsea_table.csv")

    if database != "all":
        df = df[df['database'] == database]

    df = df[df['padj'] < fdr_cutoff]

    return df
```

### 2.3 Path Handling

```python
from pathlib import Path

# Always use pathlib, not os.path
PROJECT_ROOT = Path(__file__).parent.parent.parent
DATA_DIR = PROJECT_ROOT / "03_results" / "tables"
PLOT_DIR = PROJECT_ROOT / "03_results" / "plots"

# Create directories safely
PLOT_DIR.mkdir(parents=True, exist_ok=True)

# Iterate over files
for csv_file in DATA_DIR.glob("*.csv"):
    print(csv_file.name)
```

### 2.4 Docstrings

```python
def create_volcano_plot(
    de_data: pd.DataFrame,
    fdr_cutoff: float = 0.05,
    logfc_cutoff: float = 1.0,
    output_path: Optional[Path] = None
) -> go.Figure:
    """Create interactive volcano plot from DE results.

    Creates a Plotly scatter plot with log2 fold change on x-axis
    and -log10(p-value) on y-axis. Points are colored by significance
    and direction.

    Args:
        de_data: DataFrame with columns: gene_symbol, logFC, P.Value, adj.P.Val
        fdr_cutoff: FDR threshold for significance (default: 0.05)
        logfc_cutoff: Absolute logFC threshold (default: 1.0)
        output_path: If provided, save HTML to this path

    Returns:
        Plotly Figure object

    Raises:
        ValueError: If required columns are missing

    Example:
        >>> de_data = pd.read_csv("master_de_table.csv")
        >>> fig = create_volcano_plot(de_data)
        >>> fig.show()
    """
    # Implementation
```

---

## 3. Script Structure

### 3.1 R Script Template

```r
#!/usr/bin/env Rscript
# 1.3.gsea_analysis.R - Multi-database GSEA analysis
# Project: PROJECT-ID
# Phase: 1 (Core Analysis)
#
# Description:
#   Runs GSEA across multiple MSigDB databases and saves checkpoints.
#
# Dependencies:
#   - 1.1_dge_normalized.rds
#   - 1.1_de_results.rds
#
# Outputs:
#   - 1.3_gsea_hallmark.rds
#   - 1.3_gsea_kegg.rds
#   - 1.3_gsea_reactome.rds
#
# Updated: 2025-12-15
# ============================================================================

# ============================================================================
# SETUP
# ============================================================================

message("=================================================================")
message("PROJECT-ID: GSEA Analysis")
message("=================================================================")

# Source configuration (ALWAYS FIRST)
source("02_analysis/config/config.R")
load_packages()
source_toolkit()

# ============================================================================
# LOAD DEPENDENCIES
# ============================================================================

message("[STEP] Loading dependencies...")
dge <- readRDS(file.path(DIR_CHECKPOINTS, "1.1_dge_normalized.rds"))
de_results <- readRDS(file.path(DIR_CHECKPOINTS, "1.1_de_results.rds"))

# ============================================================================
# MAIN ANALYSIS
# ============================================================================

message("[STEP] Running GSEA...")

# ... analysis code ...

# ============================================================================
# COMPLETE
# ============================================================================

message("=================================================================")
message("Script complete!")
message("Outputs saved to: ", DIR_CHECKPOINTS)
message("=================================================================")
```

### 3.2 Python Script Template

```python
#!/usr/bin/env python3
"""
3.1.pathway_explorer.py - Interactive pathway visualization dashboard

Project: PROJECT-ID
Phase: 4-5 (Python Visualization / Interactive)

Description:
    Creates an interactive HTML dashboard for exploring GSEA results
    with pathway similarity clustering.

Dependencies:
    - master_gsea_table.csv
    - master_de_table.csv

Outputs:
    - interactive/pathway_explorer.html

Updated: 2025-12-15
"""

# ============================================================================
# IMPORTS
# ============================================================================

from pathlib import Path
from typing import Dict, Optional
import pandas as pd
import plotly.express as px

# ============================================================================
# CONFIGURATION
# ============================================================================

PROJECT_ROOT = Path(__file__).parent.parent
DATA_DIR = PROJECT_ROOT / "03_results" / "tables"
OUTPUT_DIR = PROJECT_ROOT / "03_results" / "interactive"

# ============================================================================
# FUNCTIONS
# ============================================================================

def main():
    """Main entry point."""
    print("=" * 65)
    print("PROJECT-ID: Pathway Explorer")
    print("=" * 65)

    # ... implementation ...

    print("=" * 65)
    print("Complete! Output:", OUTPUT_DIR / "pathway_explorer.html")
    print("=" * 65)


# ============================================================================
# ENTRY POINT
# ============================================================================

if __name__ == "__main__":
    main()
```

---

## 4. Common Patterns

### 4.1 Configuration Loading

```r
# R: Source config at script start
source("02_analysis/config/config.R")

# Access config values
fdr_cutoff <- DE_FDR_CUTOFF
output_dir <- DIR_PLOTS
```

```python
# Python: Import config
import yaml
from pathlib import Path

with open("analysis_config.yaml") as f:
    config = yaml.safe_load(f)

fdr_cutoff = config['thresholds']['de_fdr']
```

### 4.2 Error Handling

```r
# R: tryCatch for recoverable errors
result <- tryCatch(
  run_analysis(data),
  error = function(e) {
    warning("Analysis failed: ", e$message)
    return(NULL)
  }
)

if (is.null(result)) {
  message("Skipping downstream analysis due to error")
}
```

```python
# Python: try/except
try:
    result = run_analysis(data)
except ValueError as e:
    print(f"Warning: Analysis failed: {e}")
    result = None

if result is None:
    print("Skipping downstream analysis due to error")
```

### 4.3 Directory Creation

```r
# R: Always ensure directory exists
output_dir <- file.path(DIR_PLOTS, "GSEA", db_name)
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}
```

```python
# Python: Use pathlib
output_dir = Path("03_results/plots/GSEA") / db_name
output_dir.mkdir(parents=True, exist_ok=True)
```

---

## 5. What NOT to Do

### 5.1 Avoid

```r
# Don't use T/F for TRUE/FALSE
keep <- T  # BAD
keep <- TRUE  # GOOD

# Don't use attach()
attach(data)  # BAD - pollutes namespace

# Don't use <<- for global assignment
result <<- value  # BAD - confusing scope

# Don't hardcode paths
path <- "/home/user/project/data.csv"  # BAD
path <- file.path(DATA_DIR, "data.csv")  # GOOD

# Don't suppress warnings blindly
suppressWarnings(result <- risky_function())  # BAD
```

### 5.2 R-Specific Issues

```r
# Don't use 1:length(x) in loops
for (i in 1:length(vec)) { }  # BAD - breaks if vec is empty
for (i in seq_along(vec)) { }  # GOOD

# Don't grow vectors in loops
result <- c()
for (x in data) {
  result <- c(result, process(x))  # BAD - O(n²)
}

# Use pre-allocation or lapply
result <- lapply(data, process)  # GOOD
```

### 5.3 Python-Specific Issues

```python
# Don't use mutable default arguments
def bad_function(data, result=[]):  # BAD
    result.append(data)
    return result

def good_function(data, result=None):  # GOOD
    if result is None:
        result = []
    result.append(data)
    return result

# Don't use bare except
try:
    risky_operation()
except:  # BAD - catches everything including KeyboardInterrupt
    pass

try:
    risky_operation()
except ValueError:  # GOOD - specific exception
    pass
```
