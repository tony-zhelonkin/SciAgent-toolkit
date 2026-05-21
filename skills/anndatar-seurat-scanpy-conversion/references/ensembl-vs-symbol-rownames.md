# Quirk: Ensembl IDs vs Gene Symbols on h5ad → Seurat round-trip

> Discovered during <ref-multi> Phase C (2026-04-30). Captured here for the next person who hits a `FeaturePlot("Xbp1")` that returns "feature not found".

## TL;DR

`anndataR::read_h5ad(path, as = "Seurat")` preserves whatever `var_names` the upstream h5ad has — **including Ensembl IDs**. If the scanpy pipeline that produced the h5ad used `var_names_strategy = "ensembl"` (a common pattern when symbols aren't unique enough to be safe row keys), the resulting Seurat object will have `ENSMUSG…` (or `ENSG…`) rownames and biologist-facing functions (`FeaturePlot`, `VlnPlot`, `DotPlot`, marker lookups) fail until you rename rows to gene symbols.

The `var_names` carry through. The other `var/*` columns (including `gene_name` / `feature_name` / `gene_symbol`) **do not** carry into the resulting Seurat assay's `@meta.data` — that slot only gets `var.features` and `var.features.rank`.

## Detection

```r
obj <- anndataR::read_h5ad("processed.h5ad", as = "Seurat")

# Symptom 1: rownames look like database IDs, not symbols
head(rownames(obj))            # "ENSMUSG00000000001" "ENSMUSG00000000003" …

# Symptom 2: gene_name column is missing from feature meta.data
colnames(obj[["RNA"]]@meta.data)  # only "var.features" + "var.features.rank"

# Symptom 3: lookup by symbol fails
"Xbp1" %in% rownames(obj)      # FALSE
```

If all three symptoms appear, you have the Ensembl-rowname round-trip — fix below.

## Fix

The `gene_name` column lives in `var` and is exposed by `anndataR`'s **SingleCellExperiment** form (which does carry `var` into `rowData`), even though the `as="Seurat"` form drops it. Read the same h5ad twice — once as Seurat for the assay, once as SCE for the var metadata — and apply a symbol rename.

When two Ensembl IDs map to the same symbol (~17 / 26,336 in mouse, ~50 / 38,000 in human — typically pseudogenes, paralogs, or annotation-table artifacts), **sum-collapse** on the counts layer. Don't rename `data` (log-normalized) by collapsing in place: `log(a + b) ≠ log(a) + log(b)`, so you'd corrupt the normalization. Re-derive `data` from the collapsed counts via `NormalizeData()`.

```r
library(anndataR)
library(Seurat)
library(Matrix)

H5AD <- "processed.h5ad"

# 1. Read the assay once as Seurat, and read var separately (lazily — much cheaper)
obj <- read_h5ad(H5AD, as = "Seurat")
adata_var   <- read_h5ad(H5AD)                # default = InMemoryAnnData (lazy on X)
var_df      <- as.data.frame(adata_var$var)

# Note: `read_h5ad(H5AD, as = "SingleCellExperiment")` also gives you var (in
# rowData) but materializes X a second time — for a 145k × 26k h5ad that's ~2 min
# of redundant I/O. The InMemoryAnnData form skips X materialization since we only
# touch $var.

# 2. Pull gene_name
ensembl_ids <- rownames(var_df)
symbols     <- as.character(var_df$gene_name) # adjust if column differs
empty       <- is.na(symbols) | nchar(symbols) == 0L
symbols[empty] <- ensembl_ids[empty]          # fall back to Ensembl when no symbol
stopifnot(identical(ensembl_ids, rownames(obj[["RNA"]])))

# 3. Sum-collapse counts on duplicate symbols (degenerates to a permutation if unique)
sym_factor <- factor(symbols, levels = unique(symbols))
C <- sparseMatrix(
  i = as.integer(sym_factor),
  j = seq_along(symbols),
  x = 1,
  dims = c(nlevels(sym_factor), length(symbols))
)
counts_old      <- LayerData(obj[["RNA"]], layer = "counts")
counts_collapsed <- as(C %*% counts_old, "CsparseMatrix")
rownames(counts_collapsed) <- levels(sym_factor)

# 4. Re-derive `data` (log-normalized) from collapsed counts — DO NOT collapse 'data' in place
tmp <- CreateSeuratObject(counts = CreateAssay5Object(counts = counts_collapsed))
tmp <- NormalizeData(tmp, normalization.method = "LogNormalize", scale.factor = 1e4,
                      verbose = FALSE)
data_collapsed <- LayerData(tmp[["RNA"]], layer = "data")

# 5. Rebuild the assay with symbol rownames
new_assay <- CreateAssay5Object(counts = counts_collapsed)
LayerData(new_assay, layer = "data") <- data_collapsed
# scale.data / HVGs / pca@feature.loadings: rename rows by Ensembl→symbol lookup;
# first-wins drop on the rare HVG-set collisions (see project's 04_seurat_export.R
# section 4b for the full pattern). DO NOT recompute PCA — that would change cell
# embeddings and invalidate downstream UMAP / clustering.

# 6. Verify
"Xbp1" %in% rownames(counts_collapsed)        # TRUE — biologists can now FeaturePlot
```

## Why not just use `var_names_strategy = "symbol"` upstream?

Because **symbols are not guaranteed unique**. In mouse Ensembl ~0.06% of features share a symbol; in human Ensembl ~0.13%. scanpy will silently `make_unique` on duplicate `var_names` (turning `Xbp1` into `Xbp1`, `Xbp1-1`, `Xbp1-2`), which is worse than the Ensembl round-trip — biologists query `Xbp1` and get the first one without realizing the others exist.

The right pattern is: **store Ensembl IDs as `var_names` upstream, ship symbols in `var$gene_name`**, and rename to symbols only at the analysis-output boundary (Seurat .rds for collaborators, plotting layers, etc.). The fix above is what runs at that boundary.

## Why not collapse `data` directly instead of re-deriving?

Because `LogNormalize` is `log1p(x / sum_per_cell * scale.factor)`, and `log1p(a) + log1p(b) ≠ log1p(a + b)`. Summing log-normalized values across collapsed-symbol genes produces a value that's not on the log-normalized scale anymore — it overestimates the merged feature's true `log1p` of summed counts. Re-deriving from collapsed counts is the only correct path. (Same reason 05_contrasts_de.R's symbol-collapse pattern only sums *raw counts*, never log-values.)

## Why not recompute PCA on the symbol-rowed scale.data?

For the rare HVG-set collisions, recomputing PCA would give numerically different cell embeddings. The downstream UMAP was built from the original PCA — recomputing would invalidate it (and the leiden clusterings already on disk). The pragmatic choice: keep `pca@cell.embeddings` as-is (it's a cell × dim matrix with no gene rownames), apply first-wins drop to `pca@feature.loadings` for the rare collisions, accept that the loadings are slightly inconsistent with the new `scale.data` for ~0–2 HVGs, and document. The cell embeddings remain authoritative for downstream analysis; loadings serve as a gene-contribution annotation.

## Cross-references

- `02_analysis/scripts/04_seurat_export.R` — full implementation in section 4b.
- `02_analysis/scripts/05_contrasts_de.R` — same sparse-matmul pattern applied to bulk DE pseudobulk counts (lines ~209–224).
- scverse AnnData docs on `var_names`: https://anndata.readthedocs.io/en/stable/generated/anndata.AnnData.var_names.html
