# Transposable-Element Annotation (the rare case)

Annotate a transposable-element count matrix: parse TE IDs into
`Subfamily:Family:Class`, build the TE `DGEList`, and (optionally) row-bind genes +
TEs into a unified annotated matrix. This is the rarer path — only datasets carrying
a TE count matrix from `star-te-preprocessing` need it.

> **Shared principle (also in SKILL.md):** Always annotate BEFORE filtering. Never
> drop low-count rows first — you lose IDs irreversibly.

## TE-ID label

TE matrix rows use the toolkit-canonical 3-field label **`Subfamily:Family:Class`**
(e.g. `L1Md_A:L1:LINE`). This is the same label `star-te-preprocessing` writes into
the SAF `GeneID`, so it parses cleanly downstream.

> **Authoritative parser:** TE-RNAseq-toolkit **v2.0.1** —
> `R/te_utils.R::parse_te_id` (version-pinned). Do not restate the parser logic here;
> the code is the spec.

## Inputs

- A TE count matrix (featureCounts on the grouped TE SAF; `Subfamily:Family:Class`
  rownames), e.g. `00_data/processed/featurecounts_TE/te_counts_matrix.txt`.
- Sample metadata matching the count columns (the same metadata used on the gene
  path; columns must align with the gene matrix when building a combined matrix).

## Helpers (SSoT — do not duplicate code)

From **TE-RNAseq-toolkit v2.0.1**:

- `build_te_annotation()` — `scripts/te_utils.R` — parse TE IDs into subfamily/family/class
- `parse_te_id` — `R/te_utils.R` — authoritative `Subfamily:Family:Class` parser

From **RNAseq-toolkit v2.0.0** (shared with the gene path):

- `read_counts_matrix()`, `read_metadata()`, `align_metadata_to_counts()` — `scripts/General/io_helpers.R`
- `write_annotated_matrix()` — export with annotations
- `build_dge()` — `scripts/General/dge_helpers.R` — create the edgeR `DGEList`

## How-to

```r
suppressPackageStartupMessages({
  library(data.table); library(tibble); library(dplyr); library(stringr)
})

# ---- paths (customize per project) ----
counts_te_fp  <- "00_data/processed/featurecounts_TE/te_counts_matrix.txt"
metadata_xlsx <- "00_data/metadata/Metadata.xlsx"
outdir        <- "03_results/annotated_outputs"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ---- source helpers ----
source("01_modules/RNAseq-toolkit/scripts/General/io_helpers.R")
source("01_modules/RNAseq-toolkit/scripts/General/dge_helpers.R")
source("01_modules/TE-RNAseq-toolkit/scripts/te_utils.R")   # TE-RNAseq-toolkit v2.0.1

# ---- read counts ----
message("[info] reading counts...")
mat_te  <- read_counts_matrix(counts_te_fp)     # Subfamily:Family:Class rows, e.g. L1Md_A:L1:LINE
samples <- colnames(mat_te)

# ---- read metadata ----
md <- read_metadata(metadata_xlsx)
md_aligned <- align_metadata_to_counts(md, samples)
samp_df <- md_aligned
rownames(samp_df) <- samp_df$Sample_ID

# ---- annotate TEs (BEFORE filtering) ----
message("[info] parsing TE labels...")
te_ids <- rownames(mat_te)
ann_te <- build_te_annotation(te_ids)   # uses parse_te_id internally
add_cols_te <- ann_te %>% dplyr::select(Symbol, Ensembl, subfamily, family, class)

# ---- write annotated matrix ----
message("[info] writing annotated TE matrix...")
write_annotated_matrix(mat_te, md_aligned, add_cols_te,
                       file.path(outdir, "TEs_annotated.tsv"))

# ---- build DGEList: TEs ----
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
| `TEs_annotated.tsv` | TE counts with subfamily, family, class |
| `DGEList_TEs.rds` | edgeR DGEList for TEs |

## Combined gene + TE matrix

When you want a single unified matrix (combined-mode DE), annotate genes first per
`references/gene-annotation.md` (producing `mat_gene` + `add_cols_gene`), then row-bind
the gene and TE blocks. Sample columns of the two matrices must be identical:

```r
# Align samples between matrices
stopifnot(identical(colnames(mat_gene), colnames(mat_te)))

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

write_annotated_matrix(combined_counts, md_aligned, combined_annot,
                       file.path(outdir, "genes_TEs_combined_annotated.tsv"))
```

| File | Description |
|------|-------------|
| `genes_TEs_combined_annotated.tsv` | Combined gene + TE matrix |

> The combined matrix is valid only because exonic TE loci were subtracted upstream
> (`bedtools subtract`) in `star-te-preprocessing`, so no read is double-counted. Exon
> subtraction stops double-counting (necessary) but does **not** equalize the gene/TE
> measurement bases — mutual exclusivity is necessary, not sufficient.

> **Joint normalization/DE caveats (graded options, not mandates — grades + gaps in
> `te-gene-featurecounts/SKILL.md` "Evidence & open questions"; the joint matrix itself is grade
> A).** When this combined matrix goes into joint normalization/DE:
> - **Size factors from genes only — grade B / contested:**
>   `estimateSizeFactors(dds, controlGenes = which(type == "gene"))`. TE-Seq advocates it (the
>   `-M`-inflated, long-tailed TE minority can drag gene LFCs); TEtranscripts **pools** instead.
>   Sanity-check against pooled factors.
> - **Valid for within-feature-type, across-sample DE only — grade C / inference** (gene-vs-sample,
>   TE-vs-sample).
> - **Never compare gene-vs-TE magnitude within a sample — grade C / inference**, and treat "TE %"
>   as a QC band, not biology — genes and TEs sit on different strandedness / multimapper / length
>   bases.
> - **No TPM/FPKM for TE rows — grade C / inference** — a summed multi-locus subfamily has no single
>   length; use model-normalized counts / logCPM / DESeq2 LFCs only.
>
> If TE rows were counted stranded with a sense/antisense split upstream (`--te-strand
> sense_antisense`), keep sense and antisense as separate TE features — the more principled
> best-practice (**grade B / SQuIRE-specific**) for preserving bidirectional TE biology without
> breaking gene-comparability. The field is split (TEtranscripts defaults `--stranded no`), so
> `-s 0` standalone matrices remain valid and mode-switching is not required.

## Customization Points

1. **File paths**: Update `counts_te_fp`, `metadata_xlsx`.
2. **Output directory**: Change `outdir` as needed.
3. **Metadata columns**: Adjust `samp_df` selection based on experimental design.
4. **TE format**: TE IDs are expected in `Subfamily:Family:Class` format (e.g.
   `L1Md_A:L1:LINE`), as produced by `star-te-preprocessing` and parsed by
   TE-RNAseq-toolkit v2.0.1 `R/te_utils.R::parse_te_id`.
