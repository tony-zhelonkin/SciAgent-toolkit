# Custom Database GSEA Integration

> Stage 2 of the bulk-rnaseq-gsea pipeline. See `../SKILL.md` for the full pipeline map and routing decision tree.

## Overview

This reference covers the complete workflow for adding non-MSigDB gene set databases to a clusterProfiler GSEA pipeline. It handles four database-specific parsing patterns (GMX, TSV, CSV, igraph), three gene ID mapping strategies (homologene orthologs, AnnotationDbi ID conversion, direct passthrough), quality control (size filtering, Jaccard deduplication), GSEA execution with custom T2G/T2N frames, and result normalization to the master table schema.

The core abstraction is **T2G/T2N** — every external database, regardless of its native format, is converted into two standardized data frames that clusterProfiler accepts via its `TERM2GENE` and `TERM2NAME` arguments. Once in T2G/T2N form, GSEA execution is identical to MSigDB.

**When to use this reference:**
- Adding a new gene set database (GMT, GMX, TSV, CSV, or igraph modules) to an existing GSEA pipeline
- Converting gene IDs from RefSeq, Ensembl, or Entrez to gene symbols for GSEA
- Converting human gene sets to mouse orthologs (or vice versa)
- Merging multiple custom databases with deduplication

**Route elsewhere for:**
- MSigDB collections (Hallmark, C2, C3, C5, etc.) → `references/msigdb.md`
- Visualizing GSEA results (dotplots, heatmaps, running sum) → `references/visualization.md`
- Building or appending to master CSV tables → `references/master-tables.md`
- Active metabolic module discovery (topology-aware) → use `gatom-metabolomic-predictions`

---

## Bundled Databases (Pre-packaged)

The RNAseq-toolkit ships these databases ready for immediate use via `load_reference_db()`:

| Database | Name | Description |
|----------|------|-------------|
| `mitopathways` | MitoPathways 3.0 | Hierarchical mitochondrial pathways (human, converted to mouse) |
| `mitoxplorer` | mitoXplorer 3.0 | Mitochondrial functional classifications (native mouse) |
| `mito_unified` | Unified Mito | MitoPathways + mitoXplorer merged with Jaccard deduplication |
| `transportdb` | TransportDB 2.0 | Membrane transporter families (TCDB) |

```r
source(file.path(toolkit_dir, "scripts/GSEA/GSEA_processing/load_reference_db.R"))
db <- load_reference_db("mito_unified")
# Returns list(T2G, T2N, source, created) -- ready for clusterProfiler::GSEA()
```

For GATOM network files (not bundled, ~24 MB), use `download_gatom_references()`.

Each bundled database includes `CITATIONS.bib` for proper attribution. Use `get_reference_db_info("mitopathways")` to retrieve citation info.

---

## Decision Tree

```
Have a non-MSigDB gene set database to run GSEA on?
|
+-- Is this a bundled database (mito, transportdb)?
|   +-- YES --> load_reference_db("name") -- done, skip to GSEA execution
|   +-- NO  --> continue below (BYOD workflow)
|
+-- What format is the raw file?
|   |
|   +-- GMT/GMX file
|   |   +-- Human genes, mouse analysis? --> parse_gmx() + homologene (Section: GMX)
|   |   +-- Same species?                --> parse_gmx() or parse_gmt() directly
|   |
|   +-- TSV/CSV with gene column + category column
|   |   +-- Gene symbols match your species?  --> Group by category, direct passthrough
|   |   +-- Gene IDs are RefSeq/Ensembl?      --> AnnotationDbi conversion first
|   |   +-- Gene IDs are other species?       --> homologene first, then group
|   |
|   +-- igraph object (e.g. GATOM module)
|   |   +-- Extract genes from EDGES, not vertices (Section: igraph)
|   |   +-- Lossy conversion -- topology is lost
|   |
|   +-- Two-column file (gene_set, gene) --> Already T2G format, just add prefix
|
+-- Gene IDs need conversion?
|   |
|   +-- Cross-species (human<->mouse) --> homologene (Strategy 1)
|   +-- Same-species, different ID type --> AnnotationDbi (Strategy 2)
|   +-- Already gene symbols           --> Direct passthrough (Strategy 3)
|
+-- Ready to run GSEA?
    +-- Filter gene sets by size (5-500)
    +-- Run clusterProfiler::GSEA() with TERM2GENE + TERM2NAME
    +-- Fix @geneSets slot post-GSEA
    +-- Normalize results to master table schema
    +-- Idempotent append to master_gsea_table.csv
```

---

## Quick Start: Add a New Database (Generic Recipe)

Minimal end-to-end example for a hypothetical two-column TSV database.

```r
# 1. Parse raw file into T2G/T2N
source("01_scripts/RNAseq-toolkit/scripts/GSEA/GSEA_processing/parse_external_genesets.R")

db <- parse_geneset_file(
  file     = "00_data/references/mydb/MyDatabase.tsv",
  gs_col   = "pathway",
  gene_col = "gene_symbol",
  desc_col = "description",   # Optional
  prefix   = "MYDB"
)

# 2. Filter by gene set size
source("02_analysis/helpers/pathway_utils.R")
filtered <- filter_pathways_by_size(db$T2G, db$T2N, min_size = 5, max_size = 500)
db$T2G   <- filtered$T2G
db$T2N   <- filtered$T2N

# 3. Save processed database
saveRDS(db, "00_data/references/mydb/mydb_genesets.rds")

# 4. Prepare for GSEA
T2G <- db$T2G %>% dplyr::rename(term = gs_name, gene = gene_symbol)
T2N <- db$T2N %>% dplyr::rename(term = gs_name, name = description)
T2G <- T2G %>% dplyr::filter(gene %in% names(ranked_genes))

# 5. Run GSEA
set.seed(GSEA_SEED)
gsea_result <- clusterProfiler::GSEA(
  geneList      = ranked_genes,
  TERM2GENE     = T2G,
  TERM2NAME     = T2N,
  pvalueCutoff  = 1,
  pAdjustMethod = "BH",
  eps           = 0,
  by            = "fgsea",
  nPermSimple   = 10000,
  verbose       = FALSE
)

# 6. Fix geneSets slot (required for enrichplot functions)
gsea_result@geneSets <- split(T2G$gene, T2G$term)

# 7. Normalize to master table schema and append (see master-tables.md)
```

**Verify it worked:**

```r
stopifnot(ncol(db$T2G) == 2)
stopifnot(all(colnames(db$T2G) == c("gs_name", "gene_symbol")))
stopifnot(all(grepl("^MYDB_", db$T2G$gs_name)))
sizes <- table(db$T2G$gs_name)
stopifnot(all(sizes >= 5 & sizes <= 500))
stopifnot(nrow(gsea_result@result) > 0)
```

---

## Progressive Depth

### Basic Usage: The T2G/T2N Data Model

Every custom database must be converted to two data frames before GSEA.

**T2G (TERM2GENE)** — long format, one row per gene-per-set:

| gs_name | gene_symbol |
|---------|-------------|
| `MITOPATHWAYS_Oxidative_Phosphorylation` | `Cox5a` |
| `MITOPATHWAYS_Oxidative_Phosphorylation` | `Ndufs1` |
| `TRANSPORTDB_AAAP` | `Slc1a1` |

**T2N (TERM2NAME)** — one row per gene set:

| gs_name | description |
|---------|-------------|
| `MITOPATHWAYS_Oxidative_Phosphorylation` | `Oxidative Phosphorylation` |
| `TRANSPORTDB_AAAP` | `Amino Acid/Auxin Permease` |

**Naming convention:** Gene set identifiers always carry a SCREAMING_SNAKE_CASE database prefix:

```
{DATABASE}_{pathway_name}

MITOPATHWAYS_Oxidative_Phosphorylation.Complex_I
MITOXPLORER_TRANSLATION_INITIATION
TRANSPORTDB_AAAP
GATOM_GLUTAMINE_GLUTATHIONE_AXIS
```

This prevents name collisions and allows downstream code to determine the source database from the ID alone.

**Storage format:** Parsed databases are persisted as RDS lists:

```r
db <- list(
  T2G    = data.frame(gs_name, gene_symbol),
  T2N    = data.frame(gs_name, description),
  source = "MitoPathways3.0",
  created = Sys.time()
)
saveRDS(db, "00_data/references/mitochonria/mito_mitopathways.rds")
```

### Intermediate Usage: Database-Specific Parsing Patterns

#### Pattern 1: GMX format (MitoPathways — human-to-mouse)

GMX is the transpose of GMT: row 1 = short names, row 2 = full names, rows 3+ = genes (one per row per column, variable length).

```r
mp_data <- parse_mitopathways_gmx(
  gmx_file         = "00_data/references/mitochonria/MitoPathways3.0.gmx",
  prefix           = "MITOPATHWAYS",
  convert_to_mouse = TRUE,
  org_db           = org.Mm.eg.db
)
```

The parser reads column-by-column, extracts genes from rows 3+, converts human symbols to mouse orthologs via homologene, and emits T2G/T2N. The full hierarchical name from row 2 becomes `gs_name`, the short name from row 1 becomes `description`.

#### Pattern 2: TSV format (mitoXplorer — native mouse)

Gene function file with `MGI_symbol` and `mito_process` columns. No species conversion needed.

```r
mx_data <- parse_mitoxplorer(
  gene_function_file = "00_data/references/mitochonria/mitoXplorer3.0/mouse_gene_function.txt",
  prefix             = "MITOXPLORER"
)
```

Groups genes by the `mito_process` column. Special characters are replaced with underscores:

```r
gs_name <- paste0(prefix, "_", toupper(gsub("[^A-Za-z0-9_]", "_", proc)))
```

#### Pattern 3: CSV format (TransportDB — RefSeq protein IDs)

Headerless CSV with 7 columns. Challenge: RefSeq protein accessions (`NP_*`) must be mapped to gene symbols through `org.Mm.eg.db`.

```r
transport_db <- parse_transportdb(
  transportdb_file = "00_data/references/TransportDB2.0.csv",
  prefix           = "TRANSPORTDB",
  org_db           = org.Mm.eg.db,
  group_by         = "family"
)
```

The ID conversion pipeline:
1. Strip version suffixes (`NP_001009950.1` → `NP_001009950`)
2. Primary: `AnnotationDbi::select(org_db, keys, columns = "SYMBOL", keytype = "REFSEQ")`
3. Fallback: REFSEQ → ENTREZID → SYMBOL (two-step)

Gene sets are formed by grouping mapped symbols by the `family_abbrev` column.

#### Pattern 4: igraph modules (GATOM — edge extraction)

GATOM produces igraph objects where **metabolites are nodes** and **enzymes/reactions are edges**. Extracting a gene list is lossy by design (topology and metabolites are lost).

```r
module  <- gatom_result$module
edge_df <- igraph::as_data_frame(module, "edges")

# Find the gene symbol column (varies by GATOM version)
gene_col <- NULL
for (col in c("Symbol", "symbol", "gene", "gene_symbol", "label")) {
  if (col %in% colnames(edge_df)) { gene_col <- col; break }
}
module_genes <- unique(na.omit(edge_df[[gene_col]]))
```

Biological naming uses metabolite and enzyme signature matching to generate informative names rather than "Active Module 1".

### Advanced Usage: Gene ID Mapping Strategies

#### Strategy 1: Homologene (cross-species ortholog conversion)

Used for human-to-mouse (or vice versa) conversion. Achieves 85-95% mapping with hybrid approach:

```r
# Primary: homologene package
h_res      <- homologene::homologene(human_genes, inTax = 9606, outTax = 10090)
mouse_map  <- setNames(h_res$`10090`, h_res$`9606`)

# Fallback: Title-case capitalization + org.Mm.eg.db validation
manual_mouse <- tools::toTitleCase(tolower(unmapped_genes))
valid_keys   <- AnnotationDbi::keys(org.Mm.eg.db, keytype = "SYMBOL")
manual_mouse <- manual_mouse[manual_mouse %in% valid_keys]
```

The fallback catches genes like `MT-CO1` that may be absent from homologene but have predictable capitalization patterns between species.

#### Strategy 2: AnnotationDbi (same-species ID type conversion)

Used for RefSeq, Ensembl, or Entrez to Symbol within the same species:

```r
mapping <- AnnotationDbi::select(
  org.Mm.eg.db,
  keys    = refseq_ids,
  columns = "SYMBOL",
  keytype = "REFSEQ"
)
```

Generic utility in the toolkit:

```r
T2G <- convert_geneset_ids(
  T2G,
  org_db        = org.Mm.eg.db,
  from_type     = "ENSEMBL",    # or "REFSEQ", "ENTREZID"
  to_type       = "SYMBOL",
  drop_unmapped = TRUE
)
```

#### Strategy 3: Direct passthrough (no conversion)

Used when the database already provides gene symbols matching your species (mitoXplorer mouse symbols, GATOM edge symbols). Just clean whitespace and remove NAs.

### Advanced Usage: Quality Control

#### Gene set size filtering (5-500)

```r
filter_pathways_by_size <- function(t2g, t2n, min_size = 5, max_size = 500) {
  sizes <- t2g %>%
    dplyr::group_by(gs_name) %>%
    dplyr::summarize(n = dplyr::n(), .groups = "drop")
  valid_pathways <- sizes$gs_name[sizes$n >= min_size & sizes$n <= max_size]
  list(
    T2G = t2g[t2g$gs_name %in% valid_pathways, ],
    T2N = t2n[t2n$gs_name %in% valid_pathways, ]
  )
}
```

**Rationale:** Fewer than 5 genes = unreliable enrichment statistics. More than 500 genes = too broad to be informative.

#### Jaccard deduplication (0.99 threshold)

When combining databases (e.g., MitoPathways + mitoXplorer), nearly identical gene sets are collapsed:

```r
merged_res <- merge_redundant_pathways(
  combined_t2g,
  combined_t2n,
  similarity_cutoff = 0.99,
  preference_order  = c("MITOPATHWAYS")  # keep this name on collision
)
```

This prevents inflated significance from testing the same biology twice under different names.

#### Pre-GSEA filtering to ranked list

Before running GSEA, filter T2G to genes actually present in the ranked list:

```r
T2G_filtered <- T2G %>% dplyr::filter(gene %in% names(ranked_genes))
```

### Advanced Usage: Result Normalization and Master Table Append

All results are normalized to the 13-column master table schema. See `references/master-tables.md` for the full normalization pattern. The idempotent append pattern:

```r
master_file <- file.path(DIR_TABLES, "master_gsea_table.csv")

if (file.exists(master_file)) {
  master_df <- readr::read_csv(master_file, show_col_types = FALSE)
  master_df <- master_df %>% dplyr::filter(database != "MyDB")  # Remove old rows
  master_df <- dplyr::bind_rows(master_df, export_df)           # Add fresh rows
  readr::write_csv(master_df, master_file)
} else {
  readr::write_csv(export_df, master_file)
}
```

### Advanced Usage: The @geneSets Slot Fix

After running GSEA with custom databases, the `@geneSets` slot may be empty. This must be populated for `enrichplot` functions (running sum plots, etc.) to work:

```r
T2G_list               <- split(T2G$gene, T2G$term)
gsea_result@geneSets   <- T2G_list
```

Without this fix, `gseaplot2()` and `enrichplot::gseaplot()` will fail silently or throw cryptic errors.

---

## Script Organization

| Script | Phase | Purpose |
|--------|-------|---------|
| `1.3.mito_db_prepare.R` | DB Prep | Parse and merge MitoPathways + mitoXplorer |
| `1.4.mito_gsea.R` | Execution | Run GSEA on 3 mito database variants |
| `1.9.transportdb_gsea.R` | Execution | TransportDB GSEA + master table append |
| `1.10.gatom_to_gsea.R` | Execution | GATOM module extraction + master table append |

**Execution order:**

```
1.1  core_pipeline.R         -> DE results + MSigDB GSEA checkpoints
1.3  mito_db_prepare.R       -> Parsed mito databases (RDS)
1.4  mito_gsea.R             -> Mito GSEA checkpoints
1.5  create_master_tables.R  -> master_gsea_table.csv (MSigDB + Mito)
1.9  transportdb_gsea.R      -> Appends TransportDB rows
1.10 gatom_to_gsea.R         -> Appends GATOM rows
```

Scripts 1.9 and 1.10 can run in any order after 1.5.

---

## Verification Checklist

After adding a new custom database, confirm:

- [ ] **T2G structure:** Exactly 2 columns named `gs_name` and `gene_symbol`
- [ ] **T2N structure:** Exactly 2 columns named `gs_name` and `description`
- [ ] **Prefix present:** All `gs_name` values start with `{DATABASE}_` in SCREAMING_SNAKE_CASE
- [ ] **Gene symbols match DE:** `sum(db$T2G$gene_symbol %in% rownames(de_table))` is substantial (>50% of unique symbols)
- [ ] **Size filtering applied:** No gene set has fewer than 5 or more than 500 genes
- [ ] **No duplicate gene set names:** `length(unique(db$T2G$gs_name)) == nrow(db$T2N)`
- [ ] **GSEA produced results:** `nrow(gsea_result@result) > 0`
- [ ] **geneSets slot populated:** `length(gsea_result@geneSets) > 0`
- [ ] **Master table schema:** Output has exactly 13 columns matching the schema
- [ ] **Idempotent re-run:** Running the script twice produces the same number of rows in master table
- [ ] **Checkpoint cached:** GSEA wrapped in `load_or_compute()` so re-runs skip permutation

---

## Common Pitfalls

### Pitfall: Zero overlap between gene sets and ranked list

- **Symptom:** GSEA returns 0 results. Warning: "No gene sets have size between X and Y."
- **Cause:** Gene IDs in T2G do not match the rownames of the DE table. Most common: database has human symbols but DE uses mouse, or database has Ensembl IDs but DE uses symbols.
- **Fix:** Check overlap before running GSEA:
  ```r
  overlap <- sum(db$T2G$gene_symbol %in% names(ranked_genes))
  message(sprintf("Overlap: %d / %d genes (%.0f%%)",
    overlap, length(unique(db$T2G$gene_symbol)),
    100 * overlap / length(unique(db$T2G$gene_symbol))))
  ```
  If overlap is <30%, you need ID conversion (see Gene ID Mapping Strategies).

### Pitfall: Homologene misses mitochondrial genes

- **Symptom:** MitoPathways conversion achieves only 60-70% mapping. Mitochondrially-encoded genes (MT-CO1, MT-ND1) are all missing.
- **Cause:** Mitochondrial genes often absent from homologene database. Mouse equivalents have predictable `mt-Co1` / `mt-Nd1` naming.
- **Fix:** Use the hybrid approach — homologene first, then validated title-case fallback:
  ```r
  manual_mouse <- tools::toTitleCase(tolower(unmapped_genes))
  valid_keys   <- AnnotationDbi::keys(org.Mm.eg.db, keytype = "SYMBOL")
  manual_mouse <- manual_mouse[manual_mouse %in% valid_keys]
  ```

### Pitfall: RefSeq version suffixes block AnnotationDbi lookup

- **Symptom:** `AnnotationDbi::select()` returns NA for most keys when converting RefSeq IDs.
- **Cause:** TransportDB provides `NP_001009950.1` but `org.Mm.eg.db` keys are `NP_001009950` (no version).
- **Fix:** Strip version suffixes before lookup:
  ```r
  clean_ids <- gsub("\\.\\d+$", "", refseq_ids)
  ```

### Pitfall: GATOM genes extracted from vertices instead of edges

- **Symptom:** Gene list contains metabolite labels (e.g., "pyruvate", "glutathione") instead of gene symbols.
- **Cause:** GATOM graphs have an inverted structure — metabolites are nodes, enzymes/reactions are edges. Extracting from `V(module)` returns metabolites.
- **Fix:** Always extract from edges:
  ```r
  edge_df      <- igraph::as_data_frame(module, "edges")
  module_genes <- unique(na.omit(edge_df$Symbol))
  ```

### Pitfall: enrichplot functions fail after custom GSEA

- **Symptom:** `gseaplot2()` throws "subscript out of bounds" or produces empty plots after running GSEA with custom T2G/T2N.
- **Cause:** The `@geneSets` slot of the gseaResult object is empty when using custom gene sets (clusterProfiler does not auto-populate it).
- **Fix:** Manually populate after GSEA:
  ```r
  gsea_result@geneSets <- split(T2G$gene, T2G$term)
  ```

### Pitfall: Duplicate rows in master table after re-running script

- **Symptom:** Master table row count doubles each time the script runs. Downstream plots show duplicate entries.
- **Cause:** Appending without first removing old rows from the same database.
- **Fix:** Use the idempotent pattern — filter out existing rows by `database` column before appending:
  ```r
  master_df <- master_df %>% dplyr::filter(database != "MyDb")
  master_df <- dplyr::bind_rows(master_df, export_df)
  ```

---

## Resources

- **clusterProfiler GSEA with custom gene sets:** https://yulab-smu.top/biomedical-knowledge-mining-book/universal-api.html
- **homologene R package:** https://cran.r-project.org/package=homologene
- **AnnotationDbi:** https://bioconductor.org/packages/AnnotationDbi/
- **MitoCarta 3.0:** https://www.broadinstitute.org/mitocarta
- **TransportDB 2.0:** http://www.membranetransport.org/
- **mitoXplorer:** https://mitoxplorer.ibdm.univ-mrs.fr/
