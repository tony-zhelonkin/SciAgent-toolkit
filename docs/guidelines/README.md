# RNA-seq Analysis Guidelines for AI Agents

**Version:** 2.0.0
**Last Updated:** 2025-12-15
**Purpose:** Modular methodology guidelines for AI agents performing bulk RNA-seq analysis

---

## Overview

This directory contains **universal, vendor-agnostic** methodology guidelines that ANY AI agent (Claude, Gemini, Codex, or future models) should follow when performing RNA-seq analysis.

These guidelines are the **single source of truth** for analysis patterns, coding conventions, and workflow standards. Vendor-specific context files (`CLAUDE.md`, `GEMINI.md`) should reference these guidelines rather than duplicate content.

---

## Module Index

| Module | Purpose | When to Reference |
|--------|---------|-------------------|
| [core_architecture.md](core_architecture.md) | Phased workflow, directory structure, philosophy | Starting any project |
| [data_processing.md](data_processing.md) | filterByExpr, normalization, DE workflow | Running DE analysis |
| [gsea_analysis.md](gsea_analysis.md) | GSEA patterns, msigdbr, pathway databases | Pathway enrichment |
| [checkpoint_caching.md](checkpoint_caching.md) | `load_or_compute` pattern, what to cache | Any expensive computation |
| [master_tables.md](master_tables.md) | Schema standardization, CSV contracts | Exporting results |
| [visualization.md](visualization.md) | Colors, themes, publication standards | Creating plots |
| [code_style.md](code_style.md) | R/Python conventions, function patterns | Writing any code |

---

## Quick Reference

### Core Principles

1. **Normalize Once, Visualize Many** - All data processing in Phase 1. Visualization scripts only read pre-computed results.
2. **Single Source of Truth** - Configuration in YAML, colors defined once, schemas validated everywhere.
3. **Checkpoint Everything Expensive** - Any computation >1 minute must be cached.
4. **Master Tables as Bridges** - CSV exports connect R computation to Python visualization.
5. **Separation of Concerns** - Processing scripts never plot. Visualization scripts never compute.

### Directory Convention

```
project-root/
├── 00_data/           # READ-ONLY input data
├── 01_scripts/        # Shared tools (submodules)
├── 02_analysis/       # Project-specific code
│   ├── config/        # Configuration files
│   ├── 1.x.*.R        # Phase 1: Data preparation, filtering, normalisation
│   ├── 2.x.*.R        # Phase 2: Data analysis and checkpoint caching
│   ├── 3.x.*.R        # Phase 3: e.g. R visualizations
│   └── 4.x.*.py       # Phase 4: e.g. Python/interactive viz
└── 03_results/        # Generated outputs
    ├── checkpoints/   # Cached RDS objects
    ├── tables/        # Master CSV exports
    └── plots/         # Publication figures
```

### Script Naming

```
{phase}.{order}.{description}.{ext}

Examples:
1.1.core_pipeline.R       # Phase 1, step 1
2.3.gsea_viz.R            # Phase 2, step 1
3.1.pathway_explorer.py   # Phase 4-5, step 1
```

---

## How to Use These Guidelines

### For AI Agents

When starting work on an RNA-seq project:
1. Read `core_architecture.md` for project structure understanding
2. Reference task-specific modules as needed (e.g., `gsea_analysis.md` for GSEA work)
3. Follow `code_style.md` for all code you write
4. Use patterns from `checkpoint_caching.md` for expensive operations

### For Humans

These guidelines encode best practices learned across multiple projects. They ensure:
- Consistent analysis patterns across AI assistants
- Reproducible workflows
- Maintainable code

To customize for a specific project, create a `PROJECT_CONTEXT.md` file with project-specific details (biological question, datasets, contrasts) rather than modifying these universal guidelines.

---

## Related Files

| File | Location | Purpose |
|------|----------|---------|
| `AGENTS.md` | Toolkit root | Master agent instructions (single source of truth) |
| `templates/vendor/*.template` | Toolkit templates | Thin wrappers for vendor-specific files |
| `agents/*.md` | Toolkit agents | Specialized agent definitions |

---

## Version History

- **2.0.0** (2025-12-15): Modularized from monolithic `rnaseq-analysis-guidelines.md`
- **1.0.0** (2025-12-12): Initial monolithic version
