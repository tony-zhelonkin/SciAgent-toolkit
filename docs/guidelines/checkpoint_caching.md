# Checkpoint Caching

**Module:** `checkpoint_caching.md`
**Purpose:** The `load_or_compute` pattern and caching strategy

---

## 1. The Core Pattern

### 1.1 `load_or_compute()` Function

This is the **single most important pattern** for efficient analysis:

```r
load_or_compute <- function(checkpoint_file, compute_fn,
                            force_recompute = FALSE,
                            description = "Result") {
  # Handle relative paths
  if (!grepl("^/", checkpoint_file)) {
    checkpoint_path <- file.path(DIR_CHECKPOINTS, checkpoint_file)
  } else {
    checkpoint_path <- checkpoint_file
  }

  # Check for cached result
  if (file.exists(checkpoint_path) && !force_recompute) {
    message("[CACHE] Loading ", description, " from: ", basename(checkpoint_file))
    return(readRDS(checkpoint_path))
  }

  # Compute and cache
  message("[COMPUTE] Computing ", description, "...")
  start_time <- Sys.time()
  result <- compute_fn()
  elapsed <- round(difftime(Sys.time(), start_time, units = "mins"), 2)

  message("[SAVE] Saving ", description, " (took ", elapsed, " min)")
  saveRDS(result, checkpoint_path)
  return(result)
}
```

### 1.2 Usage Pattern

```r
# Always wrap expensive operations
gsea_result <- load_or_compute(
  checkpoint_file = "1.3_gsea_hallmark.rds",
  description = "Hallmark GSEA",
  compute_fn = function() {
    run_gsea(
      DE_results = de_table,
      species = "Mus musculus",
      collection = "H",
      nperm = 100000
    )
  }
)
```

---

## 2. What to Cache

### 2.1 Always Cache (>1 minute)

| Object | Typical Time | Checkpoint Name |
|--------|--------------|-----------------|
| DGEList (normalized) | ~1 min | `1.1_dge_normalized.rds` |
| limma fit object | 2-5 min | `1.1_fit_object.rds` |
| DE results (all genes) | Part of fit | `1.1_de_results.rds` |
| GSEA per database | 5-15 min each | `1.3_gsea_{db}.rds` |
| TF activities | 5-10 min | `1.7_tf_activities.rds` |
| Pathway activities | 3-5 min | `1.8_progeny_activities.rds` |
| Gene annotations (biomaRt) | 2-5 min | `1.1_gene_annotations.rds` |

### 2.2 Sometimes Cache (30s - 1 min)

| Object | When to Cache |
|--------|---------------|
| Custom gene sets | If parsing is complex |
| Intermediate DE steps | If debugging |
| Filtered DGEList | If filtering logic is complex |

### 2.3 Never Cache

| Object | Why Not |
|--------|---------|
| Plot objects | Regenerate quickly, large files |
| Temporary variables | Not needed across sessions |
| Simple transformations | Too fast to matter |
| Configuration objects | Should be reloaded fresh |

---

## 3. Checkpoint Naming Convention

### 3.1 Standard Format

```
{phase}.{step}_{description}.rds

Components:
- phase: Which pipeline phase (1, 2, etc.)
- step: Order within phase (1, 2, 3...)
- description: What's cached (snake_case)
```

### 3.2 Examples

```
1.1_dge_raw.rds           # Initial DGEList before processing
1.1_dge_normalized.rds    # After filterByExpr + TMM
1.1_fit_object.rds        # limma fit with all contrasts
1.1_de_results.rds        # topTable output for all genes

1.3_gsea_hallmark.rds     # GSEA results for Hallmark
1.3_gsea_kegg.rds         # GSEA results for KEGG
1.3_gsea_reactome.rds     # GSEA results for Reactome
1.3_gsea_go_bp.rds        # GSEA results for GO:BP

1.4_gsea_mito.rds         # Custom mitochondrial GSEA
1.5_gsea_combined.rds     # All GSEA results combined

1.7_tf_activities.rds     # decoupleR TF activities
1.8_progeny_activities.rds # PROGENy pathway activities
```

---

## 4. Force Recompute

### 4.1 During Development

```r
# Set TRUE to regenerate during development
result <- load_or_compute(
  checkpoint_file = "1.3_gsea_hallmark.rds",
  force_recompute = TRUE,  # DEVELOPMENT MODE
  ...
)
```

### 4.2 Production Mode

```r
# Default FALSE for normal runs
result <- load_or_compute(
  checkpoint_file = "1.3_gsea_hallmark.rds",
  force_recompute = FALSE,  # PRODUCTION (default)
  ...
)
```

### 4.3 Global Force Flag

```r
# Set at script level for full recompute
FORCE_ALL <- FALSE  # Set TRUE to regenerate everything

result <- load_or_compute(
  checkpoint_file = "1.3_gsea_hallmark.rds",
  force_recompute = FORCE_ALL,
  ...
)
```

---

## 5. Checkpoint Directory Management

### 5.1 Standard Setup

```r
# In config.R
DIR_CHECKPOINTS <- "03_results/checkpoints"

# Ensure directory exists
if (!dir.exists(DIR_CHECKPOINTS)) {
  dir.create(DIR_CHECKPOINTS, recursive = TRUE)
}
```

### 5.2 Listing Checkpoints

```r
# See all cached objects
list.files(DIR_CHECKPOINTS, pattern = "\\.rds$")

# With file sizes and dates
file.info(list.files(DIR_CHECKPOINTS, pattern = "\\.rds$", full.names = TRUE))[, c("size", "mtime")]
```

### 5.3 Selective Cleanup

```r
# Remove specific checkpoint to force recompute
unlink(file.path(DIR_CHECKPOINTS, "1.3_gsea_hallmark.rds"))

# Remove all GSEA checkpoints
unlink(list.files(DIR_CHECKPOINTS, pattern = "gsea", full.names = TRUE))

# Nuclear option: remove all (use carefully!)
# unlink(DIR_CHECKPOINTS, recursive = TRUE)
```

---

## 6. Checkpoint Dependencies

### 6.1 Documenting Dependencies

```r
# In script header, document what's needed
# Dependencies: 1.1_dge_normalized.rds, 1.1_fit_object.rds
# Outputs: 1.3_gsea_hallmark.rds, 1.3_gsea_kegg.rds
```

### 6.2 Loading Dependencies

```r
# Load required checkpoints at script start
dge <- readRDS(file.path(DIR_CHECKPOINTS, "1.1_dge_normalized.rds"))
fit <- readRDS(file.path(DIR_CHECKPOINTS, "1.1_fit_object.rds"))
de_results <- readRDS(file.path(DIR_CHECKPOINTS, "1.1_de_results.rds"))
```

### 6.3 Dependency Chain

```
1.1_dge_normalized.rds
    ↓
1.1_fit_object.rds
    ↓
1.1_de_results.rds
    ↓
1.3_gsea_*.rds
    ↓
1.5_gsea_combined.rds
```

**Rule:** If upstream checkpoint changes, downstream must be regenerated.

---

## 7. Python Equivalent

### 7.1 Python Implementation

```python
import pickle
from pathlib import Path
from datetime import datetime

def load_or_compute(checkpoint_file, compute_fn, force_recompute=False, description="Result"):
    """Python equivalent of R load_or_compute pattern."""
    checkpoint_path = Path(DIR_CHECKPOINTS) / checkpoint_file

    if checkpoint_path.exists() and not force_recompute:
        print(f"[CACHE] Loading {description} from: {checkpoint_file}")
        with open(checkpoint_path, 'rb') as f:
            return pickle.load(f)

    print(f"[COMPUTE] Computing {description}...")
    start_time = datetime.now()
    result = compute_fn()
    elapsed = (datetime.now() - start_time).total_seconds() / 60

    print(f"[SAVE] Saving {description} (took {elapsed:.2f} min)")
    with open(checkpoint_path, 'wb') as f:
        pickle.dump(result, f)

    return result
```

### 7.2 Usage

```python
gsea_result = load_or_compute(
    checkpoint_file="1.3_gsea_hallmark.pkl",
    description="Hallmark GSEA",
    compute_fn=lambda: run_gsea(de_results, database="H")
)
```

---

## 8. Troubleshooting

### 8.1 Corrupted Checkpoint

```r
# If checkpoint fails to load
tryCatch(
  readRDS(checkpoint_path),
  error = function(e) {
    warning("Corrupted checkpoint, removing: ", checkpoint_path)
    unlink(checkpoint_path)
    return(NULL)
  }
)
```

### 8.2 Stale Checkpoints

If upstream data changes but checkpoints exist:

```r
# Option 1: Delete specific checkpoint
unlink(file.path(DIR_CHECKPOINTS, "1.3_gsea_hallmark.rds"))

# Option 2: Force recompute in code
result <- load_or_compute(..., force_recompute = TRUE)

# Option 3: Check file dates
checkpoint_mtime <- file.mtime(checkpoint_path)
source_mtime <- file.mtime(source_data_path)
if (source_mtime > checkpoint_mtime) {
  message("Source newer than checkpoint, recomputing...")
  force_recompute <- TRUE
}
```

### 8.3 Disk Space Management

```r
# Check checkpoint sizes
checkpoint_files <- list.files(DIR_CHECKPOINTS, pattern = "\\.rds$", full.names = TRUE)
sizes <- file.size(checkpoint_files) / 1e6  # MB
names(sizes) <- basename(checkpoint_files)
sort(sizes, decreasing = TRUE)
```
