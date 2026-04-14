---
name: louper-seurat-conversion
description: "Convert AnnData (.h5ad) to 10x Loupe Browser (.cloupe) via Seurat, using\
  \ loupeR and anndataR. Use when you need to share scRNA-seq data with collaborators\
  \ who prefer Loupe Browser, or when generating .cloupe files from non-10x data (requires\
  \ cellid_ barcode format). For R\u2194Python round-tripping without Loupe export,\
  \ use anndatar-seurat-scanpy-conversion."
license: MIT
metadata:
  skill-author: SciAgent-toolkit
  last-reviewed: 2026-04-14
  category: foundation
  tier: rich
  tags:
  - louper
  - loupe
  - seurat
  - r
  - conversion
  - cloupe
  - 10x
  complementary-skills:
  - anndatar-seurat-scanpy-conversion
  - anndata
  contraindications:
  - "Do not use create_loupe() with argument name 'counts' \u2014 v1.1.5 expects 'count_mat'."
  - Do not pass non-10x barcodes unchanged. Use cellid_XXXXXXXXX format.
  - "Do not rely on force=TRUE to bypass barcode validation \u2014 it only skips output-file\
    \ overwrite."
  version: 1.1.5
  upstream-docs: https://github.com/10XGenomics/loupeR
---

# loupeR — h5ad → Seurat → Loupe Browser

**Tested against:** loupeR 1.1.5 · anndataR 1.1.0 · Seurat 5.2.1 · R 4.5.0

---

## Purpose

Convert an AnnData `.h5ad` file to a 10x Loupe Browser `.cloupe` file via anndataR
(h5ad → Seurat) and loupeR (Seurat → .cloupe).

Use `create_loupe_from_seurat()` for standard 10x data with ACGT-16 barcodes.
Use `create_loupe()` directly for custom data (non-10x barcodes, explicit projections).
See `scripts/h5ad_to_seurat_to_loupe.R` in this skill for a generalized reference
script covering both barcode paths.

---

## Quick Reference

```r
library(anndataR); library(Seurat); library(loupeR)

# h5ad → Seurat
adata  <- read_h5ad("input.h5ad")
seurat <- adata$as_Seurat(layers_mapping = c("counts", "scvi_normalized"))

# Seurat → .cloupe (simple path, 10x ACGT barcodes only)
create_loupe_from_seurat(seurat, output_dir = "results/loupe", output_name = "export")

# Seurat → .cloupe (full control)
create_loupe(
  count_mat   = counts_mat,      # NOTE: argument is count_mat, NOT counts
  clusters    = clusters_list,
  projections = projections_list,
  output_dir  = "results/loupe",
  output_name = "export"
)
```

---

## Critical Gotchas (all verified v1.1.5)

### 1. `create_loupe()` argument is `count_mat`, not `counts`

The raw count matrix argument changed between versions. In v1.1.5 it is `count_mat`.
Using `counts =` silently fails or errors.

```r
# WRONG (old API / wrong guess)
create_loupe(counts = mat, ...)

# CORRECT
create_loupe(count_mat = mat, ...)
```

### 2. `validate_barcodes()` returns a list, not a boolean

The return value is `list($success, $msg)`. Using it as a boolean crashes.

```r
# WRONG
if (validate_barcodes(barcodes)) { ... }

# CORRECT
check <- validate_barcodes(barcodes)
if (isTRUE(check$success)) {
  cat("Valid\n")
} else {
  stop("Barcodes invalid: ", check$msg)
}
```

### 3. `force = TRUE` does NOT bypass the barcode whitelist

`force = TRUE` only skips the output-file overwrite check. It does **not** bypass the
louper binary's internal barcode validation. If your barcodes are unrecognised, the
binary will hard-fail regardless of `force = TRUE`.

### 4. Non-10x barcodes require `cellid_XXXXXXXXX` format

Custom / project-internal barcodes (e.g. `dcverse_RHP4930_2`) are not in any format
loupeR accepts. Use the Loupe-native `cellid_` format for all non-10x data:

```r
loupe_barcodes <- sprintf("cellid_%09d", seq_len(ncol(seurat)))
# e.g. "cellid_000000001", "cellid_000000002", ...

# Validate before calling create_loupe()
check <- validate_barcodes(loupe_barcodes)
if (!isTRUE(check$success)) stop(check$msg)
```

Always preserve the original IDs as a cluster metadata column so cells remain
traceable in the Loupe Browser.

---

## Barcode Formats

| Format | Pattern | Example | Notes |
|--------|---------|---------|-------|
| 10x ACGT-16 | `[ACGT]{16}(-1)?` | `ACGTACGTACGTACGT-1` | Standard 10x output |
| 10x with prefix | `\w+_[ACGT]{16}(-1)?` | `sample1_ACGTACGT...` | Multi-sample |
| 10x with suffix | `[ACGT]{16}(-1)?_\w+` | `ACGTACGT..._lib` | Multi-library |
| cellid (custom) | `cellid_[0-9]+` | `cellid_000000001` | **Use for non-10x data** |

Unrecognised formats cause a hard fail in the louper binary.

---

## Layer Mapping: AnnData → Seurat

```r
# Pass AnnData layer NAMES only (not "X")
seurat <- adata$as_Seurat(layers_mapping = c("counts", "scvi_normalized"))
```

| AnnData slot | Seurat layer name | Content |
|---|---|---|
| `layers["counts"]` | `counts` | Raw integer UMI |
| `layers["scvi_normalized"]` | `scvi_normalized` | scVI log-normalized floats |
| `X` (float32) | `X` | Auto-added; same values as scvi_normalized |

**Note:** anndataR auto-converts **all** `obsm` entries to Seurat reductions, including
internal scVI metadata keys like `_scvi_extra_categorical_covs`. Expect extra
reductions you didn't ask for. Add clean `umap` / `scvi` keys explicitly if needed.

---

## Dependencies

```r
install.packages("remotes")
remotes::install_github("10XGenomics/loupeR")

# System requirement: HDF5 must be installed at the OS level
# Ubuntu/Debian: sudo apt-get install libhdf5-dev
# macOS:         brew install hdf5
```

loupeR ≥1.1.5 requires Loupe Browser ≥8.1.

---

## Ensembl IDs

Seurat drops Ensembl IDs during import; gene links in Loupe Browser won't work without
them. To preserve:

```r
feature_ids <- read_feature_ids_from_tsv("path/to/features.tsv.gz")
create_loupe_from_seurat(seurat, feature_ids = feature_ids, ...)
```

---

## Reference Script

`scripts/h5ad_to_seurat_to_loupe.R` — generalized, configurable version of a
production export script. Covers:
- Configurable input/output paths
- Barcode auto-detection with `cellid_` fallback
- Parameterised layer names and cluster columns
- Configurable UMAP obsm key
- Cell ID traceability pattern
