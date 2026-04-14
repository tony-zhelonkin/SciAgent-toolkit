---
name: bulk-rnaseq-gsea-master-tables
description: "GSEA result normalization and master table assembly -- converts heterogeneous gseaResult objects (MSigDB, custom databases, pseudo-GSEA) into a unified CSV schema for downstream visualization. Use when assembling master_gsea_table.csv, normalizing clusterProfiler output, appending custom database results idempotently, or debugging column mismatches between R compute and Python visualization layers. For running GSEA use bulk-rnaseq-gsea-msigdb or bulk-rnaseq-gsea-custom-db. For visualization use bulk-rnaseq-gsea-visualization."
license: MIT
metadata:
  skill-author: SciAgent-toolkit
  last-reviewed: 2026-04-14
  version: 1.0.0
  upstream-docs: https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html
  category: workflow
  tier: standard
  tags:
    - gsea
    - master-table
    - normalization
    - csv-schema
    - bulk-rnaseq
    - clusterProfiler
    - fgsea
    - data-wrangling
    - pipeline
  complementary-skills:
    - bulk-rnaseq-gsea-msigdb
    - bulk-rnaseq-gsea-custom-db
    - bulk-rnaseq-gsea-visualization
  contraindications:
    - "Do not use for running GSEA. Use bulk-rnaseq-gsea-msigdb or bulk-rnaseq-gsea-custom-db instead."
    - "Do not use for creating visualizations. Use bulk-rnaseq-gsea-visualization instead."
    - "Do not use for differential expression analysis. This skill operates on GSEA output objects, not raw counts."
---

# GSEA Result Normalization and Master Table Assembly

## Overview

This skill describes how to convert heterogeneous GSEA result objects -- from MSigDB collections, custom mitochondrial databases, TransportDB, GATOM active modules, and any future gene set database -- into a single unified CSV schema (`master_gsea_table.csv`). The master table is the bridge between the R computation layer and the Python visualization layer: all downstream plots, interactive explorers, and summary statistics read from this CSV, never from raw RDS checkpoint objects.

The core pattern is **normalize-then-visualize**: each GSEA analysis produces a `gseaResult` S4 object (or a manually constructed data frame for non-GSEA sources like GATOM), which is normalized to a 13-column tibble, then appended to the master table using an idempotent filter-then-bind pattern.

**When to use this skill:**
- Assembling `master_gsea_table.csv` from checkpoint RDS files
- Adding a new database's GSEA results to an existing master table
- Debugging column name mismatches between toolkit and project normalization functions
- Understanding the master table schema for Python consumers
- Creating derived tables (`master_gsea_significant.csv`, `gsea_summary_stats.csv`)

**When NOT to use this skill:**
- Running GSEA itself (gene ranking, permutation testing) -> use `bulk-rnaseq-gsea-msigdb`
- Parsing raw gene set files (GMT, GMX, TSV) -> use `bulk-rnaseq-gsea-custom-db`
- Creating dotplots, heatmaps, or running sum plots -> use `bulk-rnaseq-gsea-visualization`

---

## Decision Tree

```
Need to work with GSEA results in tabular form?
|
+-- Have gseaResult objects in RDS checkpoints?
|   |
|   +-- First time creating master table?
|   |   --> Use 1.5.create_master_tables.R pattern (Section: Basic Usage)
|   |
|   +-- Appending a new database to existing master table?
|       --> Use idempotent append pattern (Section: Intermediate Usage)
|
+-- Have non-GSEA results (GATOM modules, TF activities)?
|   --> Manually construct schema-compliant rows (Section: Advanced Usage)
|
+-- Master table exists but Python code reports missing columns?
|   --> Check NES column casing: toolkit uses "NES", project uses "nes"
|       (Section: Common Pitfalls)
|
+-- Need filtered/summary tables?
    --> Derive from master table (Section: Derived Tables)
```

---

## Quick Start

Normalize a single gseaResult object and write to CSV:

```r
source("02_analysis/helpers/normalize_gsea.R")

# Load a checkpoint
gsea_hallmark <- readRDS("03_results/checkpoints/1.1_gsea_H_C2.rds")[["H"]]

# Normalize to standard tibble
df <- normalize_gsea_results(
  gsea_obj = gsea_hallmark,
  database = "Hallmark",
  contrast = "IL2RAKO_vs_NTC"
)

# Inspect
str(df)
# tibble [50 x 13]
#  $ pathway_id       : chr  "HALLMARK_OXIDATIVE_PHOSPHORYLATION" ...
#  $ pathway_name     : chr  "Oxidative Phosphorylation" ...
#  $ nes              : num  -2.41 -2.18 ...
#  $ pvalue           : num  1e-10 1e-08 ...
#  $ padj             : num  5e-09 3e-07 ...
#  $ set_size         : int  200 197 ...
#  $ leading_edge_size: int  89 72 ...
#  $ gene_ratio       : num  0.445 0.365 ...
#  $ core_enrichment  : chr  "Ndufs1/Ndufs2/Sdha/..." ...
#  $ database         : chr  "Hallmark" ...
#  $ contrast         : chr  "IL2RAKO_vs_NTC" ...
#  $ neg_log_padj     : num  8.30 6.52 ...
#  $ direction        : chr  "Down" "Down" ...
```

**Verify it worked:**

```r
stopifnot(ncol(df) == 13)
stopifnot(all(c("pathway_id", "pathway_name", "nes", "padj",
                "database", "contrast", "direction") %in% colnames(df)))
stopifnot(all(df$direction %in% c("Up", "Down")))
stopifnot(!any(is.na(df$nes)))
```

---

## Progressive Depth

### Basic Usage: Initial Master Table Creation

The `1.5.create_master_tables.R` script creates the master table by normalizing all MSigDB and custom GSEA checkpoints into a single CSV. This is the initial assembly step that runs after all Phase 1 GSEA scripts complete.

```r
source("02_analysis/config/config.R")
source("02_analysis/helpers/normalize_gsea.R")

# Load all GSEA checkpoints
all_gsea <- readRDS(file.path(DIR_CHECKPOINTS, "1.1_all_gsea_results.rds"))

# Database name mapping (config key -> display name)
DB_DISPLAY_NAMES <- list(
  H              = "Hallmark",
  C2_KEGG        = "KEGG",
  C2_REACTOME    = "Reactome",
  C2_WIKIPATHWAYS = "WikiPathways",
  C3_TF          = "TF_Targets",
  C5_BP          = "GO_BP",
  C5_MF          = "GO_MF",
  C5_CC          = "GO_CC"
)

# Normalize each database
all_normalized <- list()
for (db_name in names(all_gsea)) {
  display_name <- DB_DISPLAY_NAMES[[db_name]] %||% db_name
  all_normalized[[db_name]] <- normalize_gsea_results(
    gsea_obj = all_gsea[[db_name]],
    database = display_name,
    contrast = names(CONTRASTS)[1]
  )
}

# Combine and write
master_gsea <- dplyr::bind_rows(all_normalized)
readr::write_csv(master_gsea, file.path(DIR_TABLES, "master_gsea_table.csv"))
message(sprintf("Wrote master_gsea_table.csv: %d rows, %d databases",
                nrow(master_gsea), length(unique(master_gsea$database))))
```

**The `normalize_gsea_results()` function does the following internally:**

1. Extracts `@result` slot from the gseaResult S4 object (or accepts a plain data frame)
2. Detects the adjusted p-value column (`p.adjust`, `qvalue`, or `padj`)
3. Computes `leading_edge_size` by splitting `core_enrichment` on `/` and counting genes
4. Computes `gene_ratio = leading_edge_size / set_size`
5. Cleans pathway names (strips prefixes like `HALLMARK_`, title-cases, truncates)
6. Assigns `direction` as `"Up"` or `"Down"` from NES sign
7. Computes `neg_log_padj = -log10(padj)`, capped at 16

### Intermediate Usage: Idempotent Append Pattern

Scripts that run after master table creation (1.9, 1.10, etc.) must append their results without creating duplicates on re-run. The pattern uses the `database` column as a partition key:

```r
# From 1.9.transportdb_gsea.R -- idempotent append
master_file <- file.path(DIR_TABLES, "master_gsea_table.csv")

if (file.exists(master_file)) {
  master_df <- readr::read_csv(master_file, show_col_types = FALSE)

  # DELETE existing rows from this database (idempotent)
  master_df <- master_df %>%
    dplyr::filter(database != "TransportDB")

  # APPEND fresh results
  master_df <- dplyr::bind_rows(master_df, export_df)
  readr::write_csv(master_df, master_file)
  message(sprintf("Updated master table: removed old TransportDB, added %d rows", nrow(export_df)))
} else {
  # First write -- no existing file
  readr::write_csv(export_df, master_file)
}
```

**Key rule:** The string passed to `filter(database != ...)` must exactly match the `database` value used when normalizing. A mismatch (e.g., `"transportdb"` vs `"TransportDB"`) causes silent duplicates.

**Execution order matters for initial creation only:**

```
1.1 core_pipeline.R         -> MSigDB GSEA checkpoints
1.3 mito_db_prepare.R       -> Parsed mito databases (RDS)
1.4 mito_gsea.R             -> Mito GSEA checkpoints
1.5 create_master_tables.R  -> master_gsea_table.csv (MSigDB + Mito)
1.9 transportdb_gsea.R      -> Appends TransportDB rows (idempotent)
1.10 gatom_to_gsea.R        -> Appends GATOM rows (idempotent)
```

Scripts 1.9 and 1.10 can run in any order after 1.5. They can also be safely re-run.

### Advanced Usage: Manual Schema Construction (Non-GSEA Sources)

Some data sources (GATOM modules, DecoupleR TF activities, PROGENy pathways) do not produce `gseaResult` objects. Their results must be manually constructed to match the master table schema.

**GATOM modules (pseudo-NES from average logFC):**

```r
# From 1.10.gatom_to_gsea.R
entry <- data.frame(
  pathway_id        = naming_result$pathway_id,     # e.g. "GATOM_GLUTAMINE_DEPENDENT_REDOX___ANAPLEROSIS"
  pathway_name      = naming_result$pathway_name,   # e.g. "Glutamine-Dependent Redox & Anaplerosis"
  nes               = ifelse(is.na(avg_logfc), 0, avg_logfc),  # pseudo-NES from edge logFC
  pvalue            = geom_mean_pval,               # geometric mean of edge p-values
  padj              = geom_mean_pval,               # no multiple testing correction for single module
  set_size          = length(module_genes),
  leading_edge_size = length(module_genes),          # all genes are "leading edge"
  gene_ratio        = NA_real_,                      # not meaningful for topology-derived sets
  core_enrichment   = paste(module_genes, collapse = "/"),
  database          = "GATOM",
  contrast          = contrast_name,
  neg_log_padj      = -log10(max(geom_mean_pval, 1e-50)),
  direction         = ifelse(avg_logfc > 0, "Up", "Down"),
  stringsAsFactors  = FALSE
)
```

**Important caveats for pseudo-NES values:**
- GATOM pseudo-NES is the average logFC of module edges, not a permutation-derived enrichment score
- It should not be compared to MSigDB NES values on the same scale
- The pathway explorer handles this by grouping databases separately in visualizations

### Derived Tables

After the master table is assembled, several derived tables are created:

```r
# master_gsea_significant.csv -- FDR-filtered subset
sig_gsea <- master_gsea %>%
  dplyr::filter(padj < DE_FDR_CUTOFF)
readr::write_csv(sig_gsea, file.path(DIR_TABLES, "master_gsea_significant.csv"))

# gsea_summary_stats.csv -- per-database counts
summary_stats <- master_gsea %>%
  dplyr::group_by(database, contrast, direction) %>%
  dplyr::summarise(
    n_total       = dplyr::n(),
    n_significant = sum(padj < DE_FDR_CUTOFF, na.rm = TRUE),
    mean_abs_nes  = mean(abs(nes[padj < DE_FDR_CUTOFF]), na.rm = TRUE),
    .groups       = "drop"
  )
readr::write_csv(summary_stats, file.path(DIR_TABLES, "gsea_summary_stats.csv"))

# de_summary.csv -- one-row DE summary
de_summary <- data.frame(
  contrast    = contrast_name,
  total_genes = nrow(de_table),
  sig_genes   = sum(de_table$adj.P.Val < DE_FDR_CUTOFF),
  sig_up      = sum(de_table$adj.P.Val < DE_FDR_CUTOFF & de_table$logFC > 0),
  sig_down    = sum(de_table$adj.P.Val < DE_FDR_CUTOFF & de_table$logFC < 0)
)
readr::write_csv(de_summary, file.path(DIR_TABLES, "de_summary.csv"))
```

### Schema Validation (Config-Driven)

The `pipeline.yaml` configuration file defines the required columns for the master table. The Python pathway explorer validates against this schema at load time:

```yaml
# From config/pipeline.yaml
schemas:
  master_gsea_table:
    required_columns:
      - pathway_id
      - pathway_name
      - database
      - nes
      - pvalue
      - padj
      - core_enrichment
```

Python validation pattern:

```python
import yaml
import pandas as pd

with open("02_analysis/config/pipeline.yaml") as f:
    config = yaml.safe_load(f)

df = pd.read_csv("03_results/tables/master_gsea_table.csv")
required = config["schemas"]["master_gsea_table"]["required_columns"]
missing = [c for c in required if c not in df.columns]
if missing:
    raise ValueError(f"Master table missing required columns: {missing}")
```

---

## Master Table Schemas

### master_gsea_table.csv (13 core columns + optional extras)

| Column | Type | Description | Source |
|--------|------|-------------|--------|
| `pathway_id` | string | Original gene set ID with database prefix | `gseaResult@result$ID` |
| `pathway_name` | string | Cleaned human-readable name | `clean_pathway_name(Description)` |
| `nes` | float | Normalized Enrichment Score | `gseaResult@result$NES` |
| `pvalue` | float | Raw p-value from fgsea | `gseaResult@result$pvalue` |
| `padj` | float | BH-adjusted p-value | `gseaResult@result$p.adjust` or `$qvalue` |
| `set_size` | int | Genes in gene set | `gseaResult@result$setSize` |
| `leading_edge_size` | int | Genes in leading edge | Computed: `length(strsplit(core_enrichment, "/"))` |
| `gene_ratio` | float | `leading_edge_size / set_size` | Computed |
| `core_enrichment` | string | Slash-separated leading edge genes | `gseaResult@result$core_enrichment` |
| `database` | string | Source database display name | User-provided (e.g., `"Hallmark"`, `"KEGG"`) |
| `contrast` | string | Contrast name | User-provided (e.g., `"IL2RAKO_vs_NTC"`) |
| `neg_log_padj` | float | `-log10(padj)`, capped at 16 | Computed |
| `direction` | string | `"Up"` or `"Down"` from NES sign | Computed: `ifelse(nes > 0, "Up", "Down")` |

**Optional GATOM-specific columns** (present only for GATOM rows, NA for others):

| Column | Type | Description |
|--------|------|-------------|
| `module_theme` | string | GATOM module theme annotation |
| `n_metabolites` | int | Metabolites in GATOM module |
| `n_enzymes` | int | Enzymes in GATOM module |

### master_de_table.csv

| Column | Type | Description |
|--------|------|-------------|
| `contrast` | string | Contrast name |
| `gene_symbol` | string | Gene symbol |
| `ensembl_id` | string | Ensembl gene ID with version |
| `ensembl_id_base` | string | Ensembl gene ID without version |
| `logFC` | float | Log2 fold change |
| `AveExpr` | float | Average expression (log2-CPM) |
| `t` | float | Moderated t-statistic |
| `P.Value` | float | Raw p-value |
| `adj.P.Val` | float | FDR-adjusted p-value |
| `B` | float | B-statistic (log-odds of DE) |

---

## Verification Checklist

After assembling the master table, confirm:

- [ ] **Column count:** `master_gsea_table.csv` has exactly 13 core columns (or 13 + GATOM extras)
- [ ] **No duplicates:** `nrow(distinct(df, pathway_id, database, contrast)) == nrow(df)`
- [ ] **All databases present:** `unique(df$database)` includes all expected databases (Hallmark, KEGG, Reactome, WikiPathways, TF_Targets, GO_BP, GO_MF, GO_CC, MitoPathways, MitoXplorer, TransportDB, GATOM)
- [ ] **Direction validity:** All values in `direction` are `"Up"` or `"Down"`; no NAs
- [ ] **NES-direction consistency:** All rows with `nes > 0` have `direction == "Up"` and vice versa
- [ ] **neg_log_padj bounded:** No `Inf` values; maximum is 16
- [ ] **core_enrichment format:** Genes are slash-separated (`Hk2/Ldha/Slc2a1`), not comma or space
- [ ] **Schema validation passes:** Required columns from `pipeline.yaml` are all present
- [ ] **Idempotent re-run:** Running any append script twice produces the same row count

```r
# Quick validation snippet
df <- readr::read_csv("03_results/tables/master_gsea_table.csv", show_col_types = FALSE)
stopifnot(ncol(df) >= 13)
stopifnot(!any(duplicated(paste(df$pathway_id, df$database, df$contrast))))
stopifnot(all(df$direction %in% c("Up", "Down")))
stopifnot(all((df$nes > 0) == (df$direction == "Up")))
stopifnot(all(df$neg_log_padj <= 16))
stopifnot(all(grepl("/", df$core_enrichment[nchar(df$core_enrichment) > 5])))
message("All checks passed")
```

---

## Common Pitfalls

### Pitfall: NES Column Casing Mismatch (Toolkit vs Project)

- **Symptom:** Python code or downstream R scripts fail with `KeyError: 'nes'` or `Error: Column 'NES' not found`. Alternatively, `bind_rows()` produces both `NES` and `nes` columns with NAs.
- **Cause:** The toolkit's `normalize_gsea.R` outputs the column as uppercase `NES`, while the project-local `normalize_gsea.R` outputs it as lowercase `nes`. If both functions are sourced in the same session, the last one loaded wins.
- **Fix:** Use the project-local version (`02_analysis/helpers/normalize_gsea.R`) for master table assembly. It uses lowercase `nes` consistently, which matches the `pipeline.yaml` schema and the Python consumer. If you must use the toolkit version, rename after: `df <- df %>% rename(nes = NES)`.

### Pitfall: Silent Duplicates from Database Name Mismatch

- **Symptom:** The master table grows each time an append script is re-run. Row counts increase by the number of new rows without removing old ones.
- **Cause:** The `filter(database != "X")` step uses a string that does not exactly match the `database` values already in the CSV. For example, the initial creation used `"TransportDB"` but the append script filters on `"transportdb"` (case mismatch) or `"Transport_DB"` (underscore mismatch).
- **Fix:** Always use the exact same string for both the `database` parameter in `normalize_gsea_results()` and the `filter()` call. Define database display names in `config.R` as constants and reference them everywhere.

### Pitfall: Missing Database in Master Table After Append

- **Symptom:** After running an append script, the master table is missing a previously present database (e.g., Hallmark rows disappear after TransportDB append).
- **Cause:** The append script accidentally calls `filter(database != "Hallmark")` instead of `filter(database != "TransportDB")` due to copy-paste error.
- **Fix:** Always parameterize the database name. Use a variable: `DB_NAME <- "TransportDB"` and reference it in both the normalization call and the filter.

### Pitfall: Inf Values in neg_log_padj

- **Symptom:** Downstream visualizations crash or produce extreme axis ranges. Python `plotly` may fail silently.
- **Cause:** Some pathways have `padj = 0` (exact zero from fgsea when the enrichment is extreme), producing `-log10(0) = Inf`.
- **Fix:** The normalization function caps `neg_log_padj` at 16: `ifelse(is.infinite(neg_log_padj), 16, neg_log_padj)`. Verify this cap is applied. If manually constructing rows, use `neg_log_padj = -log10(max(padj, 1e-16))` or `-log10(max(padj, 1e-50))` for GATOM pseudo-p-values.

### Pitfall: Leading Edge Size Parsed as Zero

- **Symptom:** `leading_edge_size` is 0 for all rows, or `gene_ratio` is 0.
- **Cause:** The `core_enrichment` column contains comma-separated genes instead of slash-separated, or the column is NA/empty. This can happen when using a non-standard GSEA implementation or when the `@result` slot was manually constructed.
- **Fix:** Check the separator: `head(df$core_enrichment, 1)`. If commas are used, pre-convert: `df$core_enrichment <- gsub(",", "/", df$core_enrichment)`. Ensure NAs are handled: the normalization function returns 0 for NA/empty strings.

### Pitfall: padj Column Ambiguity (p.adjust vs qvalue)

- **Symptom:** Different databases show different p-value scales in the master table, or filtering by `padj < 0.05` gives unexpected results.
- **Cause:** clusterProfiler's `gseaResult@result` contains both `p.adjust` and `qvalue` columns. Some GSEA runs populate one, some populate the other. The project-local normalization prefers `qvalue` over `p.adjust`, while the toolkit version checks in order: `p.adjust` > `qvalue` > `padj`.
- **Fix:** Be aware of which column is used. The project-local function (`02_analysis/helpers/normalize_gsea.R`) uses: `if "qvalue" in colnames -> use qvalue, else use p.adjust`. This is consistent within a project but may differ if you switch between toolkit and project normalization functions.

---

## Complementary Skills

| When you need... | Use skill | Relationship |
|---|---|---|
| Run MSigDB GSEA (Hallmark, KEGG, GO, etc.) | `bulk-rnaseq-gsea-msigdb` | Prerequisite |
| Run GSEA with custom databases (MitoPathways, TransportDB) | `bulk-rnaseq-gsea-custom-db` | Prerequisite |
| Create dotplots, heatmaps, running sum plots from master table | `bulk-rnaseq-gsea-visualization` | Next step |

---

## Resources

- **clusterProfiler docs:** https://yulab-smu.top/biomedical-knowledge-mining-book/enrichment-overview.html
- **fgsea paper:** Korotkevich et al., bioRxiv 2021, doi:10.1101/060012
- **RNAseq-toolkit normalize_gsea.R:** `01_scripts/RNAseq-toolkit/scripts/GSEA/GSEA_processing/normalize_gsea.R`
- **Project-local normalize_gsea.R:** `02_analysis/helpers/normalize_gsea.R`
- **GSEA workflow documentation:** `01_scripts/RNAseq-toolkit/docs/GSEA-workflow/`
