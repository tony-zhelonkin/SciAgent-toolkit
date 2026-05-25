---
name: anndatar-seurat-scanpy-conversion
description: Convert between AnnData (.h5ad, Python/scverse) and Seurat (.rds, R/Seurat5) using the modern anndataR API. Use when bridging R and Python workflows, loading .h5ad in R, or exporting a Seurat object to scverse. For multi-assay / multimodal containers (CITE-seq, multiome) use multimodal-anndata-mudata. For Seurat→Loupe use louper-seurat-conversion.
license: MIT
metadata:
  scope: implementation
  requires: []
  skill-author: SciAgent-toolkit
  last-reviewed: 2026-04-30
  category: foundation
  tier: standard
  tags: []
  complementary-skills:
  - anndata
  - multimodal-anndata-mudata
  - louper-seurat-conversion
  contraindications:
  - Do not use for multi-assay Seurat objects (CITE-seq, multiome). Use multimodal-anndata-mudata.
  - Do not use the older SeuratDisk API. This skill uses modern anndataR.
  - Do not pass AnnData X layer name when calling as_Seurat — pass layer NAMES only (e.g., c('counts')).
  version: 1.0.0
  upstream-docs: https://scverse.org/anndataR/
---

# anndataR: Seurat ↔ h5ad Conversion Skill

## Purpose
Convert h5ad <-> Seurat using modern anndataR API

## Quick Reference

### Read h5ad → Seurat
```r
library(anndataR)
library(Seurat)

# Direct read
seurat <- read_h5ad("file.h5ad", as = "Seurat")

# Two-step (for customization)
adata <- read_h5ad("file.h5ad")
seurat <- adata$as_Seurat()
```

### Write Seurat → h5ad
```r
# Direct write
write_h5ad(seurat, "output.h5ad")

# Two-step (for customization)
adata <- as_AnnData(seurat)
adata$write_h5ad("output.h5ad")
```

## Custom Mapping
```r
# Read with mapping
seurat <- adata$as_Seurat(
  layers_mapping = c("counts", "dense_counts"),  # which layers to include
  object_metadata_mapping = c(new_name = "obs_column"),
  reduction_mapping = list(
    pca = c(key = "PC_", embeddings = "X_pca", loadings = "PCs")
  ),
  graph_mapping = TRUE  # or FALSE to skip
)

# Write with mapping
adata <- as_AnnData(seurat,
  x_mapping = "counts",  # main assay layer
  layers_mapping = c("data"),  # additional layers
  obs_mapping = c(cell_type = "celltype"),  # metadata
  obsm_mapping = list(X_pca = "pca", X_umap = "umap")  # reductions
)
```

## Structure Mapping
**Key difference: matrices are TRANSPOSED**
- AnnData: cells × genes (obs × var)
- Seurat: genes × cells (features × samples)

```
AnnData          →  Seurat
X                →  Assay layer (data/counts)
obs              →  meta.data (columns)
var              →  feature metadata
obsm (X_pca)     →  reductions (pca@cell.embeddings)
varm (PCs)       →  reductions (pca@feature.loadings)
obsp             →  graphs
uns              →  misc
layers           →  assay layers
```

## Common Issues & Solutions

### 1. Memory Issues (large files)
```r
# Use HDF5-backed for reading
adata <- read_h5ad("big.h5ad", as = "HDF5AnnData")
# Process in chunks/subsets before converting
subset <- adata[1:1000, ]
seurat <- subset$as_Seurat()
```

### 2. Missing Metadata After Conversion
```r
# Check what's available
names(adata$obs)  # cell metadata
names(adata$var)  # gene metadata
names(adata$uns)  # unstructured data

# Explicitly map
seurat <- adata$as_Seurat(
  object_metadata_mapping = c("leiden", "cell_type", "batch")
)
```

### 3. R Version Error
```r
# Requires R >= 4.5.0
getRversion()  # check version
# If < 4.5.0, install dependencies may fail
```

### 4. Layer/Assay Confusion
```r
# List available layers
names(adata$layers)
# ["counts", "data", "scale.data"]

# Specify which to convert
seurat <- adata$as_Seurat(layers_mapping = c("counts", "data"))
```

### 5. Sparse Matrix Handling
```r
# anndataR handles sparse automatically
# But check if preserved:
class(seurat[["RNA"]]$counts)  # dgCMatrix for sparse
```

### 6. Ensembl IDs as rownames (var_names round-trip)

`read_h5ad(as = "Seurat")` preserves whatever `var_names` the upstream h5ad has — **including Ensembl IDs** when the producing scanpy pipeline used `var_names_strategy = "ensembl"`. The symptom: `FeaturePlot("Xbp1")` returns "feature not found" because rownames are `ENSMUSG…`. The `var/gene_name` column is in the h5ad but NOT carried into the Seurat assay's `@meta.data`.

**Quick detection:**
```r
head(rownames(obj))                   # ENSMUSG... → quirk present
"Xbp1" %in% rownames(obj)             # FALSE → biologist queries fail
colnames(obj[["RNA"]]@meta.data)      # only var.features / var.features.rank → no gene_name
```

**Fix:** read the same h5ad as a `SingleCellExperiment` (which DOES carry `var` into `rowData`), pull `gene_name`, sum-collapse counts on duplicate symbols, re-derive `data` via `NormalizeData()` (do not collapse log-normalized values in place — `log(a+b) ≠ log(a) + log(b)`).

Full recipe — including HVG handling, scale.data + PCA-loading rename, and the reasons not to recompute PCA — at [`references/ensembl-vs-symbol-rownames.md`](references/ensembl-vs-symbol-rownames.md).

## Inspection
```r
# AnnData structure
adata
dim(adata)
names(adata$obs)      # cell metadata columns
names(adata$var)      # gene metadata columns
names(adata$obsm)     # embeddings (PCA, UMAP)
names(adata$layers)   # data layers

# Seurat structure
seurat
dim(seurat)
colnames(seurat@meta.data)  # cell metadata
Reductions(seurat)          # PCA, UMAP
Layers(seurat)              # data layers
```

## Subsetting (creates views)
```r
# Subset by condition (lazy evaluation)
subset <- adata[adata$obs$cell_type == "T cell", ]

# Convert view to concrete
concrete <- subset$as_InMemoryAnnData()

# Then convert to Seurat
seurat <- concrete$as_Seurat()
```

## Installation
```r
# Install anndataR
BiocManager::install("anndataR")

# With all dependencies
install.packages("pak")
pak::pak("scverse/anndataR", dependencies = TRUE)
```

## Multi-Assay / Multimodal Conversion

anndataR converts **one assay at a time**. For multi-assay Seurat objects (RNA + ATAC, CITE-seq), export each assay separately and assemble in Python.

### Quick Path
```r
# Export each assay
rna_adata <- as_AnnData(seurat_obj, assay_name = "RNA", x_mapping = "counts")
rna_adata$write_h5ad("rna.h5ad")

atac_adata <- as_AnnData(seurat_obj, assay_name = "ATAC")
atac_adata$write_h5ad("atac.h5ad")
```

```python
# Assemble in Python
import mudata as md
rna = sc.read_h5ad("rna.h5ad")
atac = sc.read_h5ad("atac.h5ad")
mdata = md.MuData({"rna": rna, "atac": atac})
```

### ChromatinAssay Warning

Signac's ChromatinAssay has extra slots (peak GRanges, fragment file paths) that are **lost** during anndataR conversion. Extract these manually before converting.

**For full multi-assay conversion workflow, see [multimodal-anndata-mudata.md](multimodal-anndata-mudata.md).**

---

## Debug Checklist
1. Check R version: `getRversion() >= "4.5.0"`
2. Verify file exists: `file.exists("file.h5ad")`
3. Check structure: `adata` (print object)
4. List available slots: `names(adata$obs)`, `names(adata$obsm)`
5. Test small subset first: `adata[1:100, 1:500]`
6. Check memory: `pryr::object_size(seurat)`
7. Validate conversion: compare dimensions, metadata columns