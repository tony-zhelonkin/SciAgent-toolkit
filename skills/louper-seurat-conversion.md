# loupeR - Seurat to Loupe Browser Converter skill

## Purpose
Converts Seurat objects to .cloupe files for 10x Loupe Browser visualization. Single-cell gene expression only.

## Quick Reference
**Simple**: `create_loupe_from_seurat(seurat_obj)`  
**Custom clusters/projections**: `create_loupe(counts, clusters, projections)` with `counts_matrix_from_assay(assay)` helper

## Critical Gotcha: Ensembl IDs
Seurat drops Ensembl IDs during import (keeps only feature names as rownames). Gene links in Loupe Browser won't work without manually importing them:
```r
feature_ids <- read_feature_ids_from_tsv("path/to/features.tsv.gz")
create_loupe_from_seurat(seurat_obj, feature_ids = feature_ids)
```

## Barcode Format
10x-specific: 16 ACGT characters + optional `-1` GEM well suffix. Optional prefix/suffix delimited by `_` or `:`.  
**Invalid barcodes → conversion fails silently**

## Dependencies
**HDF5 must be pre-installed** (system-level, not R package)

## Version Compatibility
LoupeR v1.1.2 requires Loupe Browser ≥8.1. Check versions if .cloupe file won't open.

---

**Installation**: Save as `/mnt/skills/user/loupeR.md`

**Test command**: 
```r
library(loupeR)
create_loupe_from_seurat(seurat_obj, 
                         feature_ids = read_feature_ids_from_tsv("features.tsv.gz"))
```