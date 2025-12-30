# Core Architecture

**Module:** `core_architecture.md`
**Purpose:** Phased workflow, directory structure, and fundamental principles

---

## 1. Fundamental Principles

### 1.1 The Five Rules

1. **Normalize Once, Visualize Many**
   All data processing happens in Phase 1. Visualization scripts only read pre-computed results.

2. **Single Source of Truth**
   Configuration lives in one place, colors defined once, schemas validated everywhere.

3. **Checkpoint Everything Expensive**
   Any computation taking >1 minute must be cached. First run: 45-60 min. Subsequent runs: 5-10 min.

4. **Master Tables as Bridges**
   CSV exports connect R computation to Python visualization. No language lock-in.

5. **Separation of Concerns**
   Data processing scripts never plot. Visualization scripts never compute. Config files never contain logic.

### 1.2 Guiding Questions Before Writing Code

- Does a checkpoint already exist for this computation?
- Is this the right phase for this code (processing vs visualization)?
- Is there a single source of truth for this parameter/color/threshold?
- Can this be done with existing toolkit functions?
- Will the output be reusable by both R and Python?

---

## 2. Phase-Based Execution Model

### 2.1 The Five Phases

| Phase | Scripts | Purpose | Output |
|-------|---------|---------|--------|
| **1** | `1.x.*.R` | Core analysis (DE, GSEA) | Checkpoints (RDS) |
| **2** | `1.x.*.R` | Master tables | CSV exports |
| **3** | `2.x.*.R` | R visualizations | PDF/PNG plots |
| **4** | `3.x.*.py` | Python visualizations | Publication figures |
| **5** | `3.x.*.py` | Interactive dashboards | HTML files |

### 2.2 Data Flow

```
Stage 0: Raw Data → Counts Matrix
├── Input: FASTQ files (external)
├── Process: QC, alignment, quantification
└── Output: counts matrix (genes × samples)

Stage 1: Counts → Normalized Data
├── Input: Raw counts + sample metadata
├── Process: Filtering + TMM/VST normalization
└── Output: DGEList object (checkpoint)

Stage 2: Normalized → Differential Expression
├── Input: Normalized counts + design formula
├── Process: limma-voom/DESeq2 modeling
└── Output: DE results per contrast (checkpoint)

Stage 3: DE Results → Pathway Enrichment
├── Input: Ranked gene list from DE
├── Process: GSEA across multiple databases
└── Output: Enrichment tables (checkpoint)

Stage 4: Results → Master Tables
├── Input: All checkpoints
├── Process: Normalize schemas, combine
└── Output: CSV files (language-agnostic)

Stage 5: Master Tables → Visualization
├── Input: CSV tables (Python) or checkpoints (R)
├── Process: Plotting only
└── Output: Publication figures
```

### 2.3 Script Dependencies

```
1.1.core_pipeline.R ─────────────────┐
    ↓                                │
1.2.annotate_de.R                    │
    ↓                                │
1.3.custom_gsea.R                    │
    ↓                                │
1.5.create_master_tables.R ←─────────┘
    ↓
┌───┴───┐
│       │
2.x.R   3.x.py
(R viz) (Python viz)
```

---

## 3. Directory Structure

### 3.1 Standard Layout

```
project-root/
├── 00_data/                    # READ-ONLY input data
│   ├── processed/              # Count matrices
│   └── references/             # Gene sets, annotations
│
├── 01_scripts/                 # Shared tools
│   ├── RNAseq-toolkit/         # Git submodule (reusable functions)
│   └── SciAgent-toolkit/       # AI agent infrastructure
│
├── 02_analysis/                # Project-specific code
│   ├── config/                 # Configuration files
│   │   ├── pipeline.yaml       # Single source of truth
│   │   ├── config.R            # R-side config loader
│   │   └── color_config.R      # Color definitions
│   ├── helpers/                # Project-specific utilities
│   ├── 1.x.*.R                 # Phase 1-2: Data processing
│   ├── 2.x.*.R                 # Phase 3: R visualizations
│   └── 3.x.*.py                # Phase 4-5: Python/interactive
│
├── 03_results/                 # Generated outputs
│   ├── checkpoints/            # Cached RDS objects
│   ├── tables/                 # Master CSV exports
│   ├── plots/                  # Publication figures
│   │   ├── QC/                 # Quality control
│   │   ├── Volcano/            # Volcano plots
│   │   └── GSEA/               # Pathway visualizations
│   └── interactive/            # HTML dashboards
│
├── CLAUDE.md                   # AI assistant context (thin wrapper)
├── GEMINI.md                   # Gemini context (thin wrapper)
├── analysis_config.yaml        # Project-specific values
├── PROJECT_CONTEXT.md          # Biological question, datasets
├── plan.md                     # Analysis strategy
├── tasks.md                    # Execution tracker
└── notes.md                    # Research findings
```

### 3.2 Data Flow Rules

1. **Immutable Inputs** - Raw data in `00_data/` is read-only. Never modify source files.
2. **Checkpoints are Intermediate** - RDS files in `checkpoints/` can be regenerated. They're cache, not archive.
3. **Master Tables are the Contract** - CSVs in `tables/` define the interface between R and Python.
4. **Outputs are Disposable** - Everything in `plots/` can be regenerated from checkpoints.

---

## 4. Naming Conventions

### 4.1 Scripts

```
{phase}.{order}.{description}.{ext}

Examples:
1.1.core_pipeline.R          # Phase 1, first script
1.5.create_master_tables.R   # Phase 1, fifth script
2.1.visualizations.R         # Phase 3 (viz), first script
3.1.pathway_explorer.py      # Phase 4-5, first script
```

### 4.2 Checkpoints

```
{phase}.{step}_{description}.rds

Examples:
1.1_dge_raw.rds              # Initial DGEList
1.1_dge_normalized.rds       # After TMM
1.1_fit_object.rds           # limma fit
1.1_de_results.rds           # topTable output
1.4_gsea_mito.rds            # Mitochondrial GSEA
```

### 4.3 Master Tables

```
master_{type}.csv

Examples:
master_gsea_table.csv        # All GSEA results
master_de_table.csv          # All DE results
master_tf_activities.csv     # TF analysis
```

### 4.4 Directories

**Convention:** lowercase with underscores (`00_data`, `checkpoints`)

---

## 5. Configuration Architecture

### 5.1 Hierarchy

```
analysis_config.yaml          # Project root - project-specific values
    ↓ imports into
02_analysis/config/config.R   # R configuration loader
02_analysis/config/config.py  # Python configuration loader
```

### 5.2 What Goes Where

| Config Type | Location | Example |
|-------------|----------|---------|
| Project metadata | `analysis_config.yaml` | Project ID, species |
| Analysis thresholds | `analysis_config.yaml` | FDR cutoffs, logFC cutoffs |
| Color palettes | `analysis_config.yaml` | Database colors, diverging scale |
| Directory paths | `config.R` / `config.py` | Derived from project structure |
| Toolkit paths | `config.R` / `config.py` | Submodule locations |

---

## 6. Anti-Patterns to Avoid

### 6.1 Mixing Computation and Visualization

```r
# BAD: Computing in viz function
volcano_plot <- function(dge, contrast) {
  fit <- lmFit(dge, design)        # WRONG - computation in viz
  fit2 <- eBayes(fit)
  results <- topTable(fit2, coef = contrast)
  ggplot(results, ...)
}

# GOOD: Visualization only
create_volcano_plot <- function(de_results, contrast_name) {
  # de_results already computed and loaded
  ggplot(de_results, aes(logFC, -log10(P.Value))) + ...
}
```

### 6.2 Hardcoding Values

```r
# BAD: Magic numbers
results %>% filter(padj < 0.05, abs(logFC) > 1)

# GOOD: Config-driven
results %>% filter(padj < DE_FDR_CUTOFF, abs(logFC) > DE_LOGFC_CUTOFF)
```

### 6.3 Skipping Checkpoints

```r
# BAD: Recomputing every time
gsea_result <- run_gsea(de_results, ...)  # Takes 10+ minutes

# GOOD: Using checkpoint pattern
gsea_result <- load_or_compute(
  checkpoint_file = "1.3_gsea_hallmark.rds",
  description = "Hallmark GSEA",
  compute_fn = function() run_gsea(de_results, ...)
)
```
