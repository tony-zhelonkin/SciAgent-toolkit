---
name: bulk-rnaseq-gsea
description: "Bulk RNA-seq GSEA pipeline router -- covers the complete workflow: MSigDB execution via clusterProfiler/fgsea (Hallmark/KEGG/Reactome/GO), custom database integration (MitoCarta/TransportDB/GMT/igraph), master table assembly (13-column CSV schema, idempotent append), and visualization (dotplots/barplots/running sum/interactive HTML). Use when running GSEA, adding a custom gene set database, normalizing gseaResult objects to CSV, or generating publication figures and interactive pathway dashboards. Routes internally to references/msigdb.md (MSigDB execution), references/custom-db.md (custom databases), references/master-tables.md (CSV schema), references/visualization.md (plots + explorer). For interactive UMAP-based pathway explorer use bulk-rnaseq-pathway-explorer. For metabolic modules use gatom-metabolomic-predictions."
license: MIT
metadata:
  skill-author: SciAgent-toolkit
  last-reviewed: 2026-04-14
  version: 1.0.0
  upstream-docs: https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html
  category: analysis
  tier: standard
  tags:
    - gsea
    - msigdb
    - clusterProfiler
    - fgsea
    - bulk-rnaseq
    - pathway-analysis
    - limma-voom
    - gene-set-enrichment
    - custom-gene-sets
    - master-table
    - visualization
    - dotplot
    - running-sum
  complementary-skills:
    - gatom-metabolomic-predictions
    - progeny-pathway-activity
    - decoupler-tf-activity
    - coresh-signature-search
    - bulk-rnaseq-pathway-explorer
  contraindications:
    - "Do not use for topology-aware metabolic module discovery. Use gatom-metabolomic-predictions instead."
    - "Do not use for TF activity inference. Use decoupler-tf-activity instead."
    - "Do not use for PROGENy pathway activity scores. Use progeny-pathway-activity instead."
    - "Do not use for signature-based GEO dataset search. Use coresh-signature-search instead."
---

# Bulk RNA-seq GSEA Pipeline

## Overview

This skill covers the complete GSEA pipeline for bulk RNA-seq analyses: from limma-voom DE results through MSigDB GSEA execution, custom gene set database integration, master table normalization, and publication-quality visualization. The four stages are sequential — each stage produces outputs consumed by the next — but each can be entered independently when prior checkpoints already exist.

**Stage map:**

| Stage | Reference | Key output |
|-------|-----------|------------|
| 1. MSigDB GSEA | `references/msigdb.md` | `gseaResult` S4 objects (one per database) |
| 2. Custom DB GSEA | `references/custom-db.md` | Additional `gseaResult` objects + T2G/T2N frames |
| 3. Master Tables | `references/master-tables.md` | `master_gsea_table.csv` (13-column, all databases) |
| 4. Visualization | `references/visualization.md` | PDFs per database + `pathway_explorer.html` |

**When to use this skill:**
- Running GSEA against MSigDB collections from limma-voom DE results
- Adding a new gene set database (MitoCarta, TransportDB, custom GMT/GMX)
- Assembling or debugging `master_gsea_table.csv`
- Creating per-database dotplots/barplots/running sum plots
- Building the interactive pathway explorer HTML dashboard

**When NOT to use this skill:**
- Topology-aware metabolic module discovery → use `gatom-metabolomic-predictions`
- TF activity inference (CollecTRI) → use `decoupler-tf-activity`
- PROGENy signaling pathway scores → use `progeny-pathway-activity`
- Interactive UMAP explorer as a standalone tool → use `bulk-rnaseq-pathway-explorer`
- Signature-based GEO dataset ranking → use `coresh-signature-search`

---

## Decision Tree

```
Where are you in the GSEA pipeline?
|
+-- Starting from DE results, need MSigDB GSEA
|   (Hallmark, KEGG, Reactome, WikiPathways, TF targets, GO)?
|   --> references/msigdb.md
|       Covers: run_gsea() wrapper, gene ranking, msigdbr v7.5/v8+ compat,
|               checkpoint batching (H+C2+C3 vs C5), eps=0 requirement
|
+-- Need to add custom databases (MitoCarta, TransportDB, GMT/GMX, GATOM)?
|   --> references/custom-db.md
|       Covers: T2G/T2N data model, bundled databases (load_reference_db()),
|               GMX/TSV/CSV/igraph parsing patterns, homologene + AnnotationDbi
|               ID mapping, Jaccard dedup, @geneSets slot fix, idempotent append
|
+-- Have gseaResult checkpoints, need to build or append master CSV?
|   --> references/master-tables.md
|       Covers: normalize_gsea_results(), 13-column schema, idempotent
|               filter+bind, padj column detection, derived tables, Python
|               schema validation, NES casing pitfall
|
+-- Have master_gsea_table.csv, need figures or interactive HTML?
    --> references/visualization.md
        Covers: 7-step per-database R pattern, cross-database pooled dotplot,
                color system (NES diverging + database palette), pathway
                explorer architecture, running sum ID/Description pitfall,
                ggplot2 4.0 NA color fix
```

---

## Quick Start

The most common path: single-database MSigDB GSEA from limma-voom results.

```r
source("02_analysis/config/config.R")
source(file.path(DIR_TOOLKIT, "GSEA/GSEA_processing/run_gsea.R"))

# Load DE results (gene symbols as rownames, 't' column present)
de_table <- readRDS("03_results/checkpoints/1.1_de_results.rds")[["IL2RA_KO - NTC"]]

# Run Hallmark GSEA
gsea_hallmark <- run_gsea(
    DE_results    = de_table,
    rank_metric   = "t",
    species       = "Mus musculus",
    collection    = "H",
    subcollection = "",
    pvalue_cutoff = 1.0,   # Retain ALL; filter downstream
    nperm         = 100000,
    seed          = 123
)

stopifnot(is(gsea_hallmark, "gseaResult"))
stopifnot(nrow(gsea_hallmark@result) > 0)
message("Hallmark GSEA: ", nrow(gsea_hallmark@result), " gene sets")
```

For the full multi-database pipeline with checkpoint batching and msigdbr version compatibility → see `references/msigdb.md`.

---

## Pipeline Architecture

### Execution order

```
1.1  core_pipeline.R        -> DE results + MSigDB GSEA checkpoints
1.3  mito_db_prepare.R      -> Parsed mito database RDS files
1.4  mito_gsea.R            -> Mito GSEA checkpoints
1.5  create_master_tables.R -> master_gsea_table.csv (MSigDB + Mito)
1.9  transportdb_gsea.R     -> Appends TransportDB rows (idempotent)
1.10 gatom_to_gsea.R        -> Appends GATOM rows (idempotent)
2.2  gsea_viz.R             -> Per-database PDFs
3.1  pathway_explorer.py    -> Interactive HTML dashboard
```

Scripts 1.9 and 1.10 can run in any order after 1.5.

### Data flow

```
limma-voom DE results (gene symbols as rownames, t-statistic)
    |
    v  run_gsea() [MSigDB via msigdbr]
gseaResult S4 objects  <-- plus T2G/T2N from custom databases
    |                           (MitoPathways, mitoXplorer, TransportDB, GATOM)
    v  normalize_gsea_results() [13-column tibble per database]
master_gsea_table.csv  [single bridge to Python]
    |
    +-- R:      gsea_dotplot / gsea_barplot / gsea_running_sum_plot
    |           -> per-database PDFs in 03_results/plots/GSEA/
    |
    +-- Python: pathway_explorer.py (Jaccard/Overlap UMAP + Plotly)
                -> pathway_explorer.html in 03_results/interactive/
```

### Checkpoint sizes

| Checkpoint | Contents | Approx size |
|-----------|----------|-------------|
| `1.1_de_results.rds` | Named list of limma topTable data frames | ~1 MB |
| `1.1_gsea_H_C2.rds` | gseaResult objects: H, C2 (KEGG/Reactome/WikiPathways), C3 | ~5 MB |
| `1.1_gsea_C5.rds` | gseaResult objects: GO BP, MF, CC | ~30 MB |
| `1.1_all_gsea_results.rds` | Combined (all 8 MSigDB databases) | ~35 MB |
| `mito_mitopathways.rds` | T2G/T2N for MitoPathways | ~0.5 MB |
| `1.4_mito_gsea.rds` | gseaResult objects (3 mito variants) | ~2 MB |
| `tables/master_gsea_table.csv` | Normalized rows, all databases | ~5 MB |

---

## Verification Checklist

The full pipeline passes when:

- [ ] **DE results loaded:** `nrow(de_table) > 5000`, `"t" %in% colnames(de_table)`, gene symbols as rownames
- [ ] **MSigDB GSEA complete:** `names(all_gsea)` contains `H, C2_KEGG, C2_REACTOME, C2_WIKIPATHWAYS, C3_TF, C5_BP, C5_MF, C5_CC`
- [ ] **Custom DB GSEA complete:** gseaResult objects present for mito databases, TransportDB, GATOM
- [ ] **Master table assembled:** `ncol(master_gsea_table) == 13`, all databases present in `unique(df$database)`
- [ ] **No duplicates:** `nrow(distinct(df, pathway_id, database, contrast)) == nrow(df)`
- [ ] **Visualization complete:** Per-database PDFs exist in `03_results/plots/GSEA/`
- [ ] **Pathway explorer built:** `pathway_explorer.html` opens in browser and renders UMAP scatter

---

## Common Pitfalls

The most critical cross-stage failure modes:

| Symptom | Stage | Cause | Fix |
|---------|-------|-------|-----|
| `nrow(gsea_result@result) == 0` | MSigDB | Ensembl IDs instead of gene symbols as rownames | `head(rownames(de_table))` should show `Il2ra`, not `ENSMUSG...` |
| p-values truncated at ~1e-4 | MSigDB | `eps` not set to 0 | Set `eps = 0` in `clusterProfiler::GSEA()` call |
| Zero overlap between gene sets and ranked list | Custom DB | Gene ID mismatch (human vs mouse, RefSeq vs symbol) | Check `sum(db$T2G$gene_symbol %in% names(ranked_genes))` |
| `gseaplot2()` fails after custom GSEA | Custom DB | Empty `@geneSets` slot | `gsea_result@geneSets <- split(T2G$gene, T2G$term)` |
| Master table doubles on re-run | Master Tables | Missing idempotent filter | `filter(database != "MyDB")` before `bind_rows()` |
| Python `KeyError: 'nes'` | Master Tables | NES column casing mismatch | Use `rename(nes = NES)` or project-local normalizer |
| Running sum plots have wrong colors | Visualization | ID/Description divergence in custom DBs | Set `Description = ID` on plotting copy; pass labels separately |
| Shape-21 points disappear in ggplot | Visualization | `color = NA` removed in ggplot2 4.0+ | Use `color = "transparent"` or `stroke = 0` |

For detailed walkthroughs of each pitfall → see the relevant reference document.

---

## Complementary Skills

| When you need... | Use skill | Relationship |
|---|---|---|
| Topology-aware metabolic module discovery (KEGG/Rhea networks) | `gatom-metabolomic-predictions` | Extension (produces modules → feed to custom-db.md) |
| TF activity inference (CollecTRI) | `decoupler-tf-activity` | Parallel analysis (results appear in pathway explorer) |
| PROGENy signaling pathway activities (14 canonical) | `progeny-pathway-activity` | Parallel analysis (results appear in pathway explorer) |
| Standalone interactive UMAP pathway explorer | `bulk-rnaseq-pathway-explorer` | Next step (consumes master_gsea_table.csv) |
| Signature-based GEO dataset ranking and hypothesis generation | `coresh-signature-search` | Alternative question (data-driven, not prior-knowledge) |

---

## Resources

- **clusterProfiler:** https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html
- **fgsea paper:** Korotkevich et al. (2021) bioRxiv, doi:10.1101/060012
- **msigdbr:** https://cran.r-project.org/package=msigdbr
- **MSigDB collections:** https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp
- **RNAseq-toolkit GSEA docs:** `01_scripts/RNAseq-toolkit/docs/GSEA-workflow/`

---

## Deeper Reference (load on demand)

- `references/msigdb.md` — MSigDB GSEA execution: `run_gsea()` wrapper, gene ranking from limma-voom, msigdbr v7.5/v8+ API compatibility, checkpoint batching strategy, `gseaResult` object structure
- `references/custom-db.md` — Custom database integration: T2G/T2N data model, bundled databases (`load_reference_db()`), four parsing patterns (GMX/TSV/CSV/igraph), three ID mapping strategies (homologene/AnnotationDbi/passthrough), Jaccard dedup, `@geneSets` fix, master table append
- `references/master-tables.md` — Master table assembly: `normalize_gsea_results()`, 13-column CSV schema, idempotent filter+bind, non-GSEA sources (GATOM pseudo-NES), derived tables, Python schema validation, column casing pitfalls
- `references/visualization.md` — Visualization: 7-step per-database R pattern, cross-database pooled dotplot, NES diverging color system, pathway explorer architecture (Jaccard/Overlap UMAP), adding a new database (7-step guide), running sum pitfalls
