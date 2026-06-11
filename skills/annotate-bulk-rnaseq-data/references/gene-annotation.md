# Gene-Symbol Annotation (the usual case)

Annotate a bulk RNA-seq gene count matrix: map Ensembl IDs to gene symbols via
Ensembl/biomaRt, then assemble an edgeR `DGEList` for downstream limma/edgeR DE.
This is the common path — every preprocessed dataset needs it.

> **Shared principle (also in SKILL.md):** Always annotate BEFORE filtering. Never
> drop low-count rows first — you lose Ensembl IDs irreversibly.

## Inputs

- A gene count matrix (featureCounts output, Ensembl IDs as rownames), e.g.
  `00_data/processed/fc_genes/count_matrices_fc/sorted_counts_matrix.txt`.
- Sample metadata matching the count columns, e.g. `00_data/metadata/Metadata.xlsx`.

## Helpers (SSoT — do not duplicate code)

From **RNAseq-toolkit v2.0.0**:

- `read_counts_matrix()` — `scripts/General/io_helpers.R` — load count matrix from file
- `read_metadata()` — `scripts/General/io_helpers.R` — load sample metadata
- `align_metadata_to_counts()` — `scripts/General/io_helpers.R` — match metadata rows to count columns
- `aggregate_duplicate_ids()` — handle duplicate row IDs
- `annotate_genes_from_ensembl()` — `scripts/General/annotate_genes.R` — add gene symbols via biomaRt/org.db
- `write_annotated_matrix()` — export counts with annotation columns
- `build_dge()` — `scripts/General/dge_helpers.R` — create the edgeR `DGEList`

## How-to

```r
suppressPackageStartupMessages({
  library(data.table); library(tibble); library(dplyr); library(stringr)
})

# ---- paths (customize per project) ----
counts_gene_fp <- "00_data/processed/fc_genes/count_matrices_fc/sorted_counts_matrix.txt"
metadata_xlsx  <- "00_data/metadata/Metadata.xlsx"
outdir         <- "03_results/annotated_outputs"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ---- source helpers (RNAseq-toolkit v2.0.0) ----
source("01_modules/RNAseq-toolkit/scripts/General/io_helpers.R")
source("01_modules/RNAseq-toolkit/scripts/General/annotate_genes.R")
source("01_modules/RNAseq-toolkit/scripts/General/dge_helpers.R")

# ---- read counts ----
message("[info] reading counts...")
mat_gene <- read_counts_matrix(counts_gene_fp)   # Ensembl rows
samples  <- colnames(mat_gene)

# ---- read metadata ----
message("[info] reading metadata...")
md <- read_metadata(metadata_xlsx)
md_aligned <- align_metadata_to_counts(md, samples)
samp_df <- md_aligned
rownames(samp_df) <- samp_df$Sample_ID

# ---- annotate genes (BEFORE filtering) ----
message("[info] annotating genes...")
mat_gene <- aggregate_duplicate_ids(mat_gene)
gene_ids <- rownames(mat_gene)
ann_gene <- annotate_genes_from_ensembl(gene_ids, try_biomart = TRUE)

# Ensure row order matches counts
ann_gene_df <- as.data.frame(ann_gene)
rownames(ann_gene_df) <- ann_gene_df$Ensembl
ann_gene_df <- ann_gene_df[match(rownames(mat_gene), rownames(ann_gene_df)), , drop = FALSE]
stopifnot(identical(rownames(mat_gene), rownames(ann_gene_df)))

# Columns for gene CSVs
add_cols_gene <- ann_gene_df %>%
  dplyr::transmute(
    Symbol,
    Ensembl,
    gene_biotype = ifelse(is.na(gene_biotype), "", gene_biotype)
  )

# ---- write annotated matrix ----
message("[info] writing annotated matrix...")
write_annotated_matrix(mat_gene, md_aligned, add_cols_gene,
                       file.path(outdir, "genes_annotated.tsv"))

# ---- build DGEList ----
message("[info] building DGEList: genes...")
dge_genes <- build_dge(mat_gene, samples_df = samp_df, genes_df = ann_gene_df, round_nonint = TRUE)
saveRDS(dge_genes, file.path(outdir, "DGEList_genes.rds"))

message("[done] Outputs written under: ", normalizePath(outdir))
```

## Output Files

| File | Description |
|------|-------------|
| `genes_annotated.tsv` | Gene counts with Symbol, Ensembl, biotype |
| `DGEList_genes.rds` | edgeR DGEList for genes |

## Customization Points

1. **File paths**: Update `counts_gene_fp`, `metadata_xlsx`.
2. **Output directory**: Change `outdir` as needed.
3. **Metadata columns**: Adjust `samp_df` selection based on experimental design.

## Combining with TEs

If this dataset also has a TE count matrix and you want a single combined gene+TE
annotated matrix and DGEList, read `references/te-annotation.md` — the gene block
produced here (`add_cols_gene`) is row-bound with the TE block there.
