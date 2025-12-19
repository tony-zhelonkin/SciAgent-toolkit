# Annotate RNAseq Data

Template for annotating count matrices with gene symbols and transposable element classifications.

## When to Use

Use this skill when you need to:
- Annotate Ensembl IDs with gene symbols before filtering
- Process transposable element (TE) count matrices
- Build DGEList objects for downstream analysis
- Combine gene and TE counts into unified matrices

## Key Principle

**Always annotate BEFORE filtering** - Never lose gene IDs by filtering first.

## Code Template

```r
# RNA-seq Annotation Template
# Annotates count matrices with gene symbols and TE families/subfamilies

suppressPackageStartupMessages({
  library(data.table); library(tibble); library(dplyr); library(stringr)
})

# ---- paths (customize per project) ----
counts_gene_fp <- "00_data/processed/fc_genes/count_matrices_fc/sorted_counts_matrix.txt"
counts_te_fp   <- "00_data/processed/featurecounts_TE/te_counts_matrix.txt"
metadata_xlsx  <- "00_data/metadata/Metadata.xlsx"
outdir         <- "03_results/annotated_outputs"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ---- source helpers ----
source("01_modules/RNAseq-toolkit/scripts/General/io_helpers.R")
source("01_modules/RNAseq-toolkit/scripts/General/annotate_genes.R")
source("01_modules/RNAseq-toolkit/scripts/General/dge_helpers.R")
source("01_modules/TE-RNAseq-toolkit/scripts/te_utils.R")

# ---- read counts ----
message("[info] reading counts...")
mat_gene <- read_counts_matrix(counts_gene_fp)   # Ensembl rows
mat_te   <- read_counts_matrix(counts_te_fp)     # HAL1:L1:LINE rows

# Align samples between matrices
stopifnot(identical(colnames(mat_gene), colnames(mat_te)))
samples <- colnames(mat_gene)

# ---- read metadata ----
message("[info] reading metadata...")
md <- read_metadata(metadata_xlsx)
md_aligned <- align_metadata_to_counts(md, samples)
samp_df <- md_aligned
rownames(samp_df) <- samp_df$Sample_ID

# ---- annotate genes ----
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

# ---- annotate TEs ----
message("[info] parsing TE labels...")
te_ids <- rownames(mat_te)
ann_te <- build_te_annotation(te_ids)
add_cols_te <- ann_te %>% dplyr::select(Symbol, Ensembl, subfamily, family, class)

# ---- build combined table ----
message("[info] building combined table...")
gene_block_cols <- c("Symbol", "Ensembl", "gene_biotype", "subfamily", "family", "class", "type")

gene_block <- add_cols_gene %>%
  dplyr::mutate(subfamily = "", family = "", class = "", type = "gene") %>%
  dplyr::select(all_of(gene_block_cols))

te_block <- add_cols_te %>%
  dplyr::mutate(gene_biotype = "", type = "TE") %>%
  dplyr::select(all_of(gene_block_cols))

combined_counts <- rbind(mat_gene, mat_te)
combined_annot  <- rbind(gene_block, te_block)
stopifnot(nrow(combined_counts) == nrow(combined_annot))

# ---- write annotated matrices ----
message("[info] writing annotated matrices...")
write_annotated_matrix(mat_gene, md_aligned, add_cols_gene,
                       file.path(outdir, "genes_annotated.tsv"))

write_annotated_matrix(mat_te, md_aligned, add_cols_te,
                       file.path(outdir, "TEs_annotated.tsv"))

write_annotated_matrix(combined_counts, md_aligned, combined_annot,
                       file.path(outdir, "genes_TEs_combined_annotated.tsv"))

# ---- build DGEList objects ----
message("[info] building DGEList: genes...")
dge_genes <- build_dge(mat_gene, samples_df = samp_df, genes_df = ann_gene_df, round_nonint = TRUE)
saveRDS(dge_genes, file.path(outdir, "DGEList_genes.rds"))

message("[info] building DGEList: TEs...")
te_genes_df <- ann_te %>%
  dplyr::select(Symbol, Ensembl = Symbol) %>%
  as.data.frame()
rownames(te_genes_df) <- rownames(mat_te)
dge_te <- build_dge(mat_te, samples_df = samp_df, genes_df = te_genes_df, round_nonint = TRUE)
saveRDS(dge_te, file.path(outdir, "DGEList_TEs.rds"))

message("[done] Outputs written under: ", normalizePath(outdir))
```

## Output Files

| File | Description |
|------|-------------|
| `genes_annotated.tsv` | Gene counts with Symbol, Ensembl, biotype |
| `TEs_annotated.tsv` | TE counts with subfamily, family, class |
| `genes_TEs_combined_annotated.tsv` | Combined matrix |
| `DGEList_genes.rds` | edgeR DGEList for genes |
| `DGEList_TEs.rds` | edgeR DGEList for TEs |

## Required Helper Functions

From `RNAseq-toolkit`:
- `read_counts_matrix()` - Load count matrix from file
- `read_metadata()` - Load sample metadata
- `align_metadata_to_counts()` - Match metadata rows to count columns
- `aggregate_duplicate_ids()` - Handle duplicate row IDs
- `annotate_genes_from_ensembl()` - Add gene symbols via biomaRt
- `write_annotated_matrix()` - Export with annotations
- `build_dge()` - Create edgeR DGEList

From `TE-RNAseq-toolkit`:
- `build_te_annotation()` - Parse TE IDs into subfamily/family/class

## Customization Points

1. **File paths**: Update `counts_gene_fp`, `counts_te_fp`, `metadata_xlsx`
2. **Output directory**: Change `outdir` as needed
3. **Metadata columns**: Adjust `samp_df` selection based on experimental design
4. **TE format**: TE IDs expected in `name:subfamily:family` format (e.g., `HAL1:L1:LINE`)
