# Gene-Symbol Annotation (the usual case)

Annotate a bulk RNA-seq gene count matrix: map Ensembl IDs to gene symbols via
Ensembl/biomaRt, then assemble an edgeR `DGEList` for downstream limma/edgeR DE.
This is the common path — every preprocessed dataset needs it.

> **Shared principle (also in SKILL.md):** Always annotate BEFORE filtering. Never
> drop low-count rows first — you lose Ensembl IDs irreversibly.

## Inputs

- A gene count matrix (featureCounts or Salmon gene-level output, Ensembl IDs as rownames), e.g.
  `00_data/processed/fc_genes/count_matrices_fc/sorted_counts_matrix.txt`.
- Sample metadata matching the count columns, e.g. `00_data/metadata/Metadata.csv` or `.xlsx`.
- Project config: `02_analysis/config/analysis_config.yaml` (keys: `project.id`,
  `project.genome_build`). `reference.ensembl_version` is **not required** — when absent,
  `biomart_version` stays `NULL` (floating; repro via provenance file only — see Risk 4 in the
  refactor plan).

## Helpers (SSoT — do not duplicate code)

From **RNAseq-toolkit v0.2.0**:

- `read_counts_matrix()` — `scripts/General/io_helpers.R` — load count matrix from file;
  auto-detects Salmon gene-level shape (`gene_id`+`gene_name`+samples) and attaches
  `attr(mat,"input_gene_name")` (NA for featureCounts)
- `aggregate_duplicate_ids()` — `scripts/General/io_helpers.R` — SUM-collapse duplicate
  stripped Ensembl IDs; output rownames in `sort(unique(ids))` order — re-match annotation
  downstream
- `read_metadata()` — `scripts/General/io_helpers.R` — load sample metadata (CSV or xlsx)
- `align_metadata_to_counts()` — `scripts/General/io_helpers.R` — match metadata rows to
  count columns
- `write_annotated_matrix()` — `scripts/General/io_helpers.R` — export counts with annotation
  columns and embedded factor rows (table B — transposed format)
- `annotate_genes_from_ensembl()` — `scripts/General/annotate_genes.R` — add gene symbols via
  biomaRt/org.db; emits `Symbol, Ensembl, mgi_symbol, ENTREZID, gene_biotype, input_gene_name`
- `build_dge()` — `scripts/General/dge_helpers.R` — create the edgeR `DGEList`
- `write_gene_index()` — `scripts/General/provenance.R` — write `*_gene_index.tsv` with
  renamed cols (`ensembl_gene_id, SYMBOL, mgi_symbol, ENTREZID, gene_biotype, input_gene_name`)
- `write_session_provenance()` — `scripts/General/provenance.R` — write `sessionInfo()` +
  genome build + resolved biomaRt archive release

## How-to

```r
suppressPackageStartupMessages({
  library(data.table); library(tibble); library(dplyr); library(stringr); library(yaml)
})

# ---- config (drive all paths and IDs from here) ----
cfg        <- yaml::read_yaml("02_analysis/config/analysis_config.yaml")
project_id <- cfg$project$id
genome_build <- cfg$project$genome_build
# reference.ensembl_version is absent from this config — biomart_version stays NULL
# (floating; reproducibility via provenance file only)
biomart_ver <- NULL   # set to e.g. 111 if cfg$reference$ensembl_version exists

outdir <- "03_results/annotated_outputs"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

counts_gene_fp <- "00_data/processed/fc_genes/count_matrices_fc/sorted_counts_matrix.txt"
metadata_fp    <- "00_data/metadata/Metadata.csv"   # or .xlsx; read_metadata() dispatches

# ---- source helpers (RNAseq-toolkit v0.2.0) ----
source("01_modules/RNAseq-toolkit/scripts/General/io_helpers.R")
source("01_modules/RNAseq-toolkit/scripts/General/annotate_genes.R")
source("01_modules/RNAseq-toolkit/scripts/General/dge_helpers.R")
source("01_modules/RNAseq-toolkit/scripts/General/provenance.R")

# ---- read counts ----
message("[info] reading counts...")
mat_gene <- read_counts_matrix(counts_gene_fp)   # Ensembl rows (versioned IDs stripped)
samples  <- colnames(mat_gene)

# ---- aggregate duplicate stripped Ensembl IDs (SUM-collapse) ----
# Output rownames are in sort(unique(ids)) order — NOT input order.
# Annotation must be re-matched to rownames(mat_gene) AFTER this step.
message("[info] aggregating duplicate Ensembl IDs...")
mat_gene <- aggregate_duplicate_ids(mat_gene)

# ---- read metadata ----
message("[info] reading metadata...")
md <- read_metadata(metadata_fp)
md_aligned <- align_metadata_to_counts(md, samples)
samp_df <- md_aligned
rownames(samp_df) <- samp_df$Sample_ID

# ---- annotate genes (BEFORE filtering) ----
# Pass attr(mat_gene,"input_gene_name") so Salmon gene_name is threaded through.
# Pass biomart_version= for pinned Ensembl release (NULL = floating).
message("[info] annotating genes...")
gene_ids <- rownames(mat_gene)
ann_gene <- annotate_genes_from_ensembl(
  gene_ids,
  try_biomart     = TRUE,
  input_gene_name = attr(mat_gene, "input_gene_name"),
  biomart_version = biomart_ver
)

# ---- RE-ORDER annotation to rownames(mat_gene) via match() ----
# aggregate_duplicate_ids() outputs rownames in alphabetical order; the annotation
# tibble must be aligned back to the matrix before any downstream step.
ann_df <- as.data.frame(ann_gene)
rownames(ann_df) <- ann_df$Ensembl
ann_df <- ann_df[match(rownames(mat_gene), rownames(ann_df)), , drop = FALSE]
stopifnot(identical(rownames(mat_gene), rownames(ann_df)))

# write_gene_index expects the 6-col tibble with Symbol/Ensembl casing (it renames internally)
# Build write_gene_index-ready frame; note: ann_df still has Symbol/Ensembl casing here.
ann_for_index <- ann_df  # write_gene_index renames Ensembl->ensembl_gene_id, Symbol->SYMBOL

# For write_annotated_matrix + alignment check we need the renamed key
ann_df$ensembl_gene_id <- ann_df$Ensembl
stopifnot(identical(ann_df$ensembl_gene_id, rownames(mat_gene)))

# ---- derive factor columns for the transposed matrix (table B) ----
# Adapt factor derivation to your experimental design (G4: factor_cols is caller responsibility).
# Example for a 3-factor experiment:
factor_cols <- c("LPS_treat", "IFNg_treat", "genotype")
# md_aligned must contain these columns (derive them from Group/Treatment if needed):
# md_aligned$LPS_treat  <- if_else(grepl("LPS",  md_aligned$Treatment), "LPS",  "ctrl")
# md_aligned$IFNg_treat <- if_else(grepl("IFNg", md_aligned$Treatment), "IFNg", "ctrl")
# md_aligned$genotype   <- md_aligned$Genotype

# Annotation columns for the matrix body
add_cols_gene <- ann_df %>%
  dplyr::transmute(
    Symbol,
    Ensembl,
    gene_biotype = ifelse(is.na(gene_biotype), "", gene_biotype)
  )

# ---- write table B: transposed annotated matrix with factor rows at TOP ----
# Rows = factor rows (one per factor_col) + gene rows
# Format: factor name in Symbol col, "" in Ensembl col, factor values across sample columns.
message("[info] writing counts_transposed (table B)...")
write_annotated_matrix(
  mat_gene, md_aligned, add_cols_gene,
  file.path(outdir, paste0(project_id, "_counts_transposed.tsv")),
  factor_cols = factor_cols
)

# ---- write table A: wide-by-sample matrix ----
# Rows = samples, columns = factor columns + gene columns.
message("[info] writing counts_wide_by_sample (table A)...")
wide_mat <- cbind(
  md_aligned[, factor_cols, drop = FALSE],
  as.data.frame(t(mat_gene))
)
data.table::fwrite(
  wide_mat,
  file.path(outdir, paste0(project_id, "_counts_wide_by_sample.tsv")),
  sep = "\t", row.names = FALSE
)

# ---- build DGEList ----
message("[info] building DGEList: genes...")
dge_genes <- build_dge(mat_gene, samples_df = samp_df, genes_df = ann_df, round_nonint = TRUE)
saveRDS(dge_genes, file.path(outdir, paste0(project_id, "_DGEList.rds")))

# ---- write gene index (table C) ----
# write_gene_index renames Symbol->SYMBOL and Ensembl->ensembl_gene_id internally.
message("[info] writing gene index (table C)...")
write_gene_index(
  ann_for_index,
  file.path(outdir, paste0(project_id, "_gene_index.tsv"))
)

# ---- write session provenance ----
message("[info] writing session provenance...")
write_session_provenance(
  file.path(outdir, paste0(project_id, "_provenance.txt")),
  genome_build    = genome_build,
  ensembl_version = biomart_ver
)

message("[done] Outputs written under: ", normalizePath(outdir))
```

## Output Files

| File | Description |
|------|-------------|
| `<id>_counts_wide_by_sample.tsv` | Table A — wide-by-sample matrix: rows = samples, columns = factor cols + gene cols |
| `<id>_counts_transposed.tsv` | Table B — transposed annotated matrix: factor rows at TOP, then gene rows (Symbol, Ensembl, biotype) |
| `<id>_gene_index.tsv` | Table C — gene dictionary: `ensembl_gene_id, SYMBOL, mgi_symbol, ENTREZID, gene_biotype, input_gene_name` |
| `<id>_DGEList.rds` | edgeR DGEList for genes (counts + genes + samples + TMM) |
| `<id>_provenance.txt` | Session provenance: `sessionInfo()` + genome build + resolved biomaRt archive release |

## Customization Points

1. **Config keys**: `project.id` and `project.genome_build` drive all output filenames and
   provenance. Add `reference.ensembl_version: 111` (or your release) to the config and set
   `biomart_ver <- cfg$reference$ensembl_version` for pinned Ensembl reproducibility.
2. **File paths**: Update `counts_gene_fp`, `metadata_fp`.
3. **Factor derivation**: derive `factor_cols` columns in `md_aligned` before calling
   `write_annotated_matrix()`. The toolkit does not derive factors — it only embeds named columns.
4. **Salmon vs featureCounts**: `read_counts_matrix()` auto-detects both shapes.
   For Salmon gene-level, `attr(mat_gene,"input_gene_name")` is populated automatically.

## Notes on the `aggregate_duplicate_ids` + `match()` pattern

After SUM-collapsing duplicate stripped Ensembl IDs, `rownames(mat_gene)` are in
`sort(unique(ids))` order — NOT the original input order. The `annotate_genes_from_ensembl()`
call can return rows in a different order. The `match()` + `stopifnot` block above is the
mandatory alignment step (mirrors reference script lines 219/267): it ensures the annotation
data frame and the count matrix share the same row order before building the DGEList or writing
any output.

## Combining with TEs

If this dataset also has a TE count matrix and you want a single combined gene+TE
DGEList, read `references/te-annotation.md` — the gene annotation produced here feeds
`create_combined_dge()` there.
