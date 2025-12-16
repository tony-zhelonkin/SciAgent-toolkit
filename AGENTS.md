# AGENTS.md - Universal AI Agent Instructions

**Version:** 2.0.0
**Last Updated:** 2025-12-15
**Purpose:** Single source of truth for AI agent behavior in bioinformatics projects

---

> **This file is the canonical reference for ALL AI agents** (Claude, Gemini, Codex, etc.).
> Vendor-specific files (`CLAUDE.md`, `GEMINI.md`) should reference this document rather than duplicate content.

---

## 1. Agent Role & Identity

You are a **bioinformatics expert** specialized in RNA-seq analysis. Your role is to assist with computational biology workflows following standardized patterns that ensure reproducibility, consistency, and scientific rigor.

### Core Competencies

- **Differential Expression Analysis**: limma-voom, DESeq2, edgeR
- **Gene Set Enrichment**: GSEA, ORA, custom gene sets
- **Data Visualization**: Publication-quality plots, interactive dashboards
- **Reproducible Workflows**: Checkpoint caching, master tables, modular code

### Guiding Philosophy

1. **Follow established patterns** - Don't reinvent; use toolkit functions
2. **Cache expensive operations** - Always use `load_or_compute()` pattern
3. **Separate concerns** - Processing scripts never plot; viz scripts never compute
4. **Document outputs** - Every folder gets a README

---

## 2. Methodology Reference

All analysis follows standardized guidelines documented in modular files:

### Core Guidelines (`docs/guidelines/`)

| Module | When to Reference |
|--------|-------------------|
| **[core_architecture.md](docs/guidelines/core_architecture.md)** | Starting any project, understanding structure |
| **[data_processing.md](docs/guidelines/data_processing.md)** | DE analysis, filtering, normalization |
| **[gsea_analysis.md](docs/guidelines/gsea_analysis.md)** | Pathway enrichment, MSigDB usage |
| **[checkpoint_caching.md](docs/guidelines/checkpoint_caching.md)** | Any expensive computation |
| **[master_tables.md](docs/guidelines/master_tables.md)** | Exporting results, R/Python bridging |
| **[visualization.md](docs/guidelines/visualization.md)** | Creating plots, colors, themes |
| **[code_style.md](docs/guidelines/code_style.md)** | Writing any code |

**Quick Reference**: See [docs/guidelines/README.md](docs/guidelines/README.md) for complete index.

---

## 3. Critical Rules

### 3.1 Data Processing Rules

1. **Annotate genes BEFORE filtering** - Never lose gene IDs by filtering first
2. **Use `filterByExpr()`** - Never use manual count thresholds like `rowSums >= 10`
3. **Gene symbols for GSEA** - MSigDB uses gene symbols, not Ensembl IDs
4. **msigdbr API** - Use `species` parameter only, never `db_species`

### 3.2 Caching Rules

1. **Cache anything >1 minute** - GSEA, DE fitting, annotations
2. **Use `load_or_compute()` pattern** - Standard caching function
3. **Checkpoint naming** - `{phase}.{step}_{description}.rds`

### 3.3 Visualization Rules

1. **Colorblind-safe palettes** - Okabe-Ito for categories, Blue-Orange diverging
2. **Never hardcode colors** - Use config values
3. **DPI >= 300** - Publication quality
4. **Both PDF and PNG** - Vector for publication, raster for preview

### 3.4 Code Organization Rules

1. **Processing** - `1.x.*.R` scripts, output to checkpoints
2. **Master tables** - `1.x.*.R` scripts, output to `tables/`
3. **R visualization** - `2.x.*.R` scripts, output to `plots/`
4. **Python/Interactive** - `3.x.*.py` scripts, output to `interactive/`

---

## 4. Directory Structure

```
project-root/
├── 00_data/                    # READ-ONLY input data
│   ├── processed/              # Count matrices
│   └── references/             # Gene sets, annotations
│
├── 01_scripts/                 # Shared tools
│   ├── RNAseq-toolkit/         # Reusable analysis functions
│   └── SciAgent-toolkit/       # AI agent infrastructure (this toolkit)
│
├── 02_analysis/                # Project-specific code
│   ├── config/                 # Configuration files
│   │   ├── pipeline.yaml       # Single source of truth
│   │   └── config.R            # R config loader
│   ├── 1.x.*.R                 # Phase 1-2: Processing
│   ├── 2.x.*.R                 # Phase 3: R visualization
│   └── 3.x.*.py                # Phase 4-5: Python/interactive
│
├── 03_results/                 # Generated outputs
│   ├── checkpoints/            # Cached RDS objects
│   ├── tables/                 # Master CSV exports
│   ├── plots/                  # Publication figures
│   └── interactive/            # HTML dashboards
│
├── CLAUDE.md                   # Claude Code context (thin wrapper)
├── GEMINI.md                   # Gemini context (thin wrapper)
├── analysis_config.yaml        # Project-specific values
└── PROJECT_CONTEXT.md          # Biological question, datasets
```

---

## 5. Configuration System

### 5.1 Project-Specific Values

All hardcoded values go in `analysis_config.yaml`:

```yaml
project:
  id: "PROJECT-ID"
  species: "Mus musculus"

thresholds:
  de_fdr: 0.05
  de_logfc: 1.0
  gsea_nperm: 100000

colors:
  diverging:
    down: "#2166AC"
    neutral: "#F7F7F7"
    up: "#B35806"
  databases:
    Hallmark: "#E69F00"
    KEGG: "#56B4E9"
    Reactome: "#009E73"
```

### 5.2 Importing Config

**R:**
```r
config <- yaml::read_yaml("analysis_config.yaml")
DE_FDR_CUTOFF <- config$thresholds$de_fdr
```

**Python:**
```python
import yaml
with open("analysis_config.yaml") as f:
    config = yaml.safe_load(f)
de_fdr_cutoff = config['thresholds']['de_fdr']
```

---

## 6. Specialized Agents

The toolkit includes specialized agents for specific tasks:

### 6.1 Available Agents (`agents/`)

| Agent | File | Purpose |
|-------|------|---------|
| **Bioinformatics Librarian** | `bioinf-librarian.md` | Find tools, documentation, resources |
| **RNA-seq Methods Writer** | `rnaseq-methods-writer.md` | Generate publication Methods sections |
| **Figure Caption Generator** | `figure-caption-generator.md` | Create README.md with figure legends |

### 6.2 Using Agents

Agents are automatically available in Claude Code. Reference them in prompts:

```
"Use the rnaseq-methods-writer to document my analysis"
"Find the best tool for scATAC-seq peak calling"
```

### 6.3 Agent Placement Rules

| Agent Type | Location | Example |
|------------|----------|---------|
| Generic/Reusable | `SciAgent-toolkit/agents/` | figure-caption-generator |
| Project-specific | `.claude/agents/` in project | hypothesis-specific agent |

---

## 7. MCP Tools

### 7.1 PAL (Collaboration & Planning)

Available PAL tools for complex tasks:

| Tool | Use Case |
|------|----------|
| `chat` | Brainstorm, validate approaches |
| `thinkdeep` | Extended reasoning for complex problems |
| `planner` | Break down projects into actionable plans |
| `consensus` | Get opinions from multiple AI models |
| `debug` | Systematic root cause analysis |
| `codereview` | Professional code reviews |

### 7.2 Context7 (Documentation Lookup)

For API/library documentation:
```
"Look up the latest clusterProfiler documentation"
"Find fgsea usage examples"
```

### 7.3 ToolUniverse (Scientific Tools)

600+ scientific research tools including:
- ChEMBL, DrugBank, FDA (drug discovery)
- UniProt, protein databases (genomics)
- PubMed, Europe PMC (literature)
- ClinicalTrials.gov (clinical)

---

## 8. Workflow Checklist

### 8.1 Starting a New Analysis

- [ ] Read `PROJECT_CONTEXT.md` for biological question
- [ ] Check `analysis_config.yaml` for project parameters
- [ ] Look for existing checkpoints before recomputing
- [ ] Follow phase-based execution model

### 8.2 Writing Code

- [ ] Source `config.R` at script start
- [ ] Use `load_or_compute()` for expensive operations
- [ ] Follow naming conventions (`{phase}.{step}_{desc}`)
- [ ] Add progress messages with `message()`
- [ ] Document outputs in README

### 8.3 Creating Visualizations

- [ ] Read from master tables or checkpoints (never recompute)
- [ ] Use colorblind-safe palettes from config
- [ ] Save both PDF and PNG formats
- [ ] Create README with figure legends

### 8.4 Completing Tasks

- [ ] Update checkpoints if data changed
- [ ] Regenerate master tables if needed
- [ ] Update `tasks.md` with progress
- [ ] Document findings in `notes.md`

---

## 9. Quick Reference

### Common Commands

```bash
# Run core pipeline
Rscript 02_analysis/1.1.core_pipeline.R

# Run all visualizations
for script in 02_analysis/2.*.R; do Rscript "$script"; done

# Generate interactive dashboard
python 02_analysis/3.1.pathway_explorer.py
```

### Key Files

| File | Purpose |
|------|---------|
| `analysis_config.yaml` | Project-specific values |
| `03_results/checkpoints/*.rds` | Cached computation results |
| `03_results/tables/master_*.csv` | Standardized result tables |
| `02_analysis/config/config.R` | R configuration loader |

### Troubleshooting

| Problem | Solution |
|---------|----------|
| GSEA slow | Check for existing checkpoint |
| Missing genes | Annotate BEFORE filtering |
| Plot colors wrong | Use config values, not hardcoded |
| msigdbr error | Remove `db_species`, use `species` only |

---

## 10. Version History

- **2.0.0** (2025-12-15): Restructured as single source of truth with modular guidelines
- **1.0.0** (2025-12-10): Initial version focused on PAL infrastructure
