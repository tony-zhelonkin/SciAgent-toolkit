# Seurat Multimodal Analysis Skill

## Overview

Comprehensive skill for multimodal single-cell analysis using the Seurat ecosystem (R). Covers CITE-seq (RNA+Protein), 10x Multiome (RNA+ATAC), weighted nearest neighbors integration, and downstream chromatin analysis with Signac.

**Key Features:**
- Multi-assay Seurat object management (RNA, ADT, ATAC, peaks)
- Weighted Nearest Neighbors (WNN) for multimodal integration
- CITE-seq protein quantification and clustering
- 10x Multiome joint RNA+ATAC analysis
- Signac for chromatin accessibility analysis
- ChromVAR for motif enrichment scoring
- Peak-gene linkage analysis
- Reference mapping with Azimuth/bridge integration

**Documentation:** https://satijalab.org/seurat/

---

## Installation

```r
# Core Seurat + Signac
install.packages("Seurat")
install.packages("Signac")

# For motif analysis
BiocManager::install(c("JASPAR2020", "TFBSTools", "BSgenome.Hsapiens.UCSC.hg38",
                        "BSgenome.Mmusculus.UCSC.mm10", "chromVAR", "motifmatchr"))

# For reference mapping
install.packages("SeuratData")
install.packages("Azimuth")

# Optional: BPCells for large datasets
install.packages("BPCells")
```

---

## Part 1: CITE-seq Analysis (RNA + Protein)

### 1.1 Create Multimodal Seurat Object

```r
library(Seurat)
library(SeuratData)

# Load CITE-seq data (example: PBMC dataset)
# InstallData("bmcite")
bm <- LoadData(ds = "bmcite")

# Object structure: multiple assays in one object
Assays(bm)
# [1] "RNA" "ADT"

# Access different assays
DefaultAssay(bm) <- "RNA"  # Switch to RNA
DefaultAssay(bm) <- "ADT"  # Switch to protein
```

### 1.2 From Raw Data (10x Cell Ranger)

```r
library(Seurat)

# Read 10x CITE-seq outputs
counts <- Read10X("filtered_feature_bc_matrix/")
# Returns list: Gene Expression, Antibody Capture

# Create Seurat object with RNA
seurat_obj <- CreateSeuratObject(
  counts = counts$`Gene Expression`,
  project = "CITE-seq",
  min.cells = 3,
  min.features = 200
)

# Add protein (ADT) as separate assay
seurat_obj[["ADT"]] <- CreateAssay5Object(
  counts = counts$`Antibody Capture`
)
```

### 1.3 CITE-seq QC and Normalization

```r
# QC for both modalities
DefaultAssay(bm) <- "RNA"
bm <- PercentageFeatureSet(bm, pattern = "^MT-", col.name = "percent.mt")

# Filter cells
bm <- subset(bm,
  nFeature_RNA > 200 &
  nFeature_RNA < 5000 &
  percent.mt < 10
)

# Normalize RNA (SCTransform recommended)
DefaultAssay(bm) <- "RNA"
bm <- SCTransform(bm, verbose = FALSE)

# Normalize ADT (CLR recommended for protein)
DefaultAssay(bm) <- "ADT"
bm <- NormalizeData(bm, normalization.method = "CLR", margin = 2)

# Variable features and scaling for ADT
# Note: For ADT, typically use all antibodies (not HVG selection)
VariableFeatures(bm) <- rownames(bm[["ADT"]])
bm <- ScaleData(bm)
```

### 1.4 Dimensionality Reduction (Each Modality)

```r
# RNA: PCA
DefaultAssay(bm) <- "SCT"
bm <- RunPCA(bm, verbose = FALSE)

# ADT: PCA on protein space
DefaultAssay(bm) <- "ADT"
bm <- RunPCA(bm, reduction.name = "apca", verbose = FALSE)
```

### 1.5 Weighted Nearest Neighbors (WNN)

```r
# Find multimodal neighbors using WNN
# Combines RNA (PCA) and ADT (apca) information
bm <- FindMultiModalNeighbors(
  bm,
  reduction.list = list("pca", "apca"),
  dims.list = list(1:30, 1:18),  # Dimensions per modality
  modality.weight.name = "RNA.weight"  # Store weights
)

# WNN UMAP (joint embedding)
bm <- RunUMAP(
  bm,
  nn.name = "weighted.nn",  # Use WNN graph
  reduction.name = "wnn.umap",
  reduction.key = "wnnUMAP_"
)

# WNN clustering
bm <- FindClusters(
  bm,
  graph.name = "wsnn",  # WNN shared nearest neighbor graph
  algorithm = 3,  # SLM algorithm
  resolution = 0.5,
  verbose = FALSE
)

# Visualize
DimPlot(bm, reduction = "wnn.umap", label = TRUE)
```

### 1.6 Modality Weights Analysis

```r
# Check how much each modality contributes per cell
VlnPlot(bm, features = "RNA.weight", group.by = "celltype", pt.size = 0)

# High RNA weight = RNA more informative for that cell type
# High ADT weight (1 - RNA.weight) = protein more informative
```

### 1.7 Protein Marker Visualization

```r
# Protein feature plots
DefaultAssay(bm) <- "ADT"
FeaturePlot(bm,
  features = c("CD4", "CD8", "CD14", "CD19"),
  reduction = "wnn.umap",
  cols = c("lightgrey", "darkgreen"),
  ncol = 2
)

# Dotplot comparing RNA vs protein
DefaultAssay(bm) <- "ADT"
p1 <- DotPlot(bm, features = rownames(bm), group.by = "celltype") +
  RotatedAxis() + ggtitle("Protein")

DefaultAssay(bm) <- "SCT"
# Use RNA gene names that match ADT
p2 <- DotPlot(bm, features = c("CD4", "CD8A", "CD14", "CD19"), group.by = "celltype") +
  RotatedAxis() + ggtitle("RNA")

p1 | p2
```

---

## Part 2: 10x Multiome Analysis (RNA + ATAC)

### 2.1 Load 10x Multiome Data

```r
library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)  # or EnsDb.Mmusculus.v79 for mouse

# Load Cell Ranger ARC outputs
counts <- Read10X_h5("filtered_feature_bc_matrix.h5")
fragpath <- "atac_fragments.tsv.gz"

# Get gene annotations for ATAC
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- "UCSC"

# Create Seurat with RNA
seurat_obj <- CreateSeuratObject(
  counts = counts$`Gene Expression`,
  assay = "RNA"
)

# Create ChromatinAssay for ATAC
seurat_obj[["ATAC"]] <- CreateChromatinAssay(
  counts = counts$Peaks,
  sep = c(":", "-"),
  fragments = fragpath,
  annotation = annotations
)
```

### 2.2 Multiome QC

```r
# RNA QC
DefaultAssay(seurat_obj) <- "RNA"
seurat_obj <- PercentageFeatureSet(seurat_obj, pattern = "^MT-", col.name = "percent.mt")

# ATAC QC (requires fragment file)
DefaultAssay(seurat_obj) <- "ATAC"
seurat_obj <- NucleosomeSignal(seurat_obj)
seurat_obj <- TSSEnrichment(seurat_obj)

# Visualize QC metrics
VlnPlot(
  seurat_obj,
  features = c("nCount_RNA", "nFeature_RNA", "percent.mt",
               "nCount_ATAC", "TSS.enrichment", "nucleosome_signal"),
  ncol = 3,
  pt.size = 0
)

# Filter cells
seurat_obj <- subset(
  seurat_obj,
  nCount_ATAC < 100000 &
  nCount_RNA < 25000 &
  nCount_ATAC > 1000 &
  nCount_RNA > 1000 &
  nucleosome_signal < 2 &
  TSS.enrichment > 1
)
```

### 2.3 Process RNA Modality

```r
DefaultAssay(seurat_obj) <- "RNA"

# Normalize and find variable features
seurat_obj <- SCTransform(seurat_obj, verbose = FALSE)

# PCA
seurat_obj <- RunPCA(seurat_obj, verbose = FALSE)
```

### 2.4 Process ATAC Modality

```r
DefaultAssay(seurat_obj) <- "ATAC"

# TF-IDF normalization (standard for ATAC)
seurat_obj <- RunTFIDF(seurat_obj)

# Feature selection for ATAC
seurat_obj <- FindTopFeatures(seurat_obj, min.cutoff = "q0")

# LSI (Latent Semantic Indexing) - ATAC equivalent of PCA
seurat_obj <- RunSVD(seurat_obj)

# Check for correlation with sequencing depth
DepthCor(seurat_obj)
# If component 1 correlates with depth, exclude from downstream
```

### 2.5 WNN Integration

```r
# Multimodal neighbors using RNA PCA + ATAC LSI
seurat_obj <- FindMultiModalNeighbors(
  seurat_obj,
  reduction.list = list("pca", "lsi"),
  dims.list = list(1:50, 2:40),  # LSI: skip component 1 if correlated with depth
  modality.weight.name = "RNA.weight"
)

# Joint UMAP
seurat_obj <- RunUMAP(
  seurat_obj,
  nn.name = "weighted.nn",
  reduction.name = "wnn.umap",
  reduction.key = "wnnUMAP_"
)

# Clustering
seurat_obj <- FindClusters(
  seurat_obj,
  graph.name = "wsnn",
  algorithm = 3,
  resolution = 0.5
)

DimPlot(seurat_obj, reduction = "wnn.umap", label = TRUE)
```

### 2.6 Gene Activity Scores

```r
# Compute gene activity (accessibility near gene body)
DefaultAssay(seurat_obj) <- "ATAC"
gene.activities <- GeneActivity(seurat_obj)

# Add as new assay
seurat_obj[["GeneActivity"]] <- CreateAssay5Object(counts = gene.activities)

# Normalize
seurat_obj <- NormalizeData(
  seurat_obj,
  assay = "GeneActivity",
  normalization.method = "LogNormalize"
)

# Compare RNA expression vs gene activity
DefaultAssay(seurat_obj) <- "SCT"
p1 <- FeaturePlot(seurat_obj, features = "MS4A1", reduction = "wnn.umap") +
  ggtitle("MS4A1 RNA")

DefaultAssay(seurat_obj) <- "GeneActivity"
p2 <- FeaturePlot(seurat_obj, features = "MS4A1", reduction = "wnn.umap") +
  ggtitle("MS4A1 Activity")

p1 | p2
```

---

## Part 3: Signac Chromatin Analysis

### 3.1 Peak Calling and Visualization

```r
library(Signac)

DefaultAssay(seurat_obj) <- "ATAC"

# Coverage plot for a region
CoveragePlot(
  seurat_obj,
  region = "chr14-99700000-99760000",  # Example: BCL11B locus
  features = "BCL11B",
  assay = "ATAC",
  extend.upstream = 1000,
  extend.downstream = 1000
)

# Coverage by group
CoveragePlot(
  seurat_obj,
  region = "MS4A1",  # Gene name works too
  group.by = "seurat_clusters"
)
```

### 3.2 Differential Accessibility

```r
DefaultAssay(seurat_obj) <- "ATAC"

# Find DA peaks between clusters
da_peaks <- FindMarkers(
  seurat_obj,
  ident.1 = "0",
  ident.2 = "1",
  min.pct = 0.05,
  test.use = "LR",  # Logistic regression recommended for ATAC
  latent.vars = "nCount_ATAC"  # Control for sequencing depth
)

# Get closest gene for top DA peaks
closest_genes <- ClosestFeature(seurat_obj, regions = rownames(da_peaks))
da_peaks$gene <- closest_genes$gene_name
```

### 3.3 Motif Analysis with ChromVAR

```r
library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Hsapiens.UCSC.hg38)  # Match your genome

DefaultAssay(seurat_obj) <- "ATAC"

# Get JASPAR motifs
pwm_set <- getMatrixSet(
  JASPAR2020,
  opts = list(collection = "CORE", tax_group = "vertebrates", all_versions = FALSE)
)

# Add motif information to assay
seurat_obj <- AddMotifs(
  seurat_obj,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  pfm = pwm_set
)

# Run ChromVAR
seurat_obj <- RunChromVAR(
  seurat_obj,
  genome = BSgenome.Hsapiens.UCSC.hg38
)

# Access chromVAR scores
DefaultAssay(seurat_obj) <- "chromvar"

# Plot motif activity
FeaturePlot(seurat_obj, features = "MA0139.1", reduction = "wnn.umap")  # CTCF

# Differential motif activity
differential.activity <- FindMarkers(
  seurat_obj,
  ident.1 = "B_cells",
  ident.2 = "T_cells",
  only.pos = TRUE,
  mean.fxn = rowMeans,
  fc.name = "avg_diff"
)

# Get motif names
motif.names <- differential.activity %>%
  rownames() %>%
  ConvertMotifID(seurat_obj)
```

### 3.4 Peak-Gene Linkage

```r
DefaultAssay(seurat_obj) <- "ATAC"

# Link peaks to genes (correlation-based)
seurat_obj <- RegionStats(seurat_obj, genome = BSgenome.Hsapiens.UCSC.hg38)
seurat_obj <- LinkPeaks(
  seurat_obj,
  peak.assay = "ATAC",
  expression.assay = "SCT",
  genes.use = c("MS4A1", "CD3D", "CD8A")  # Or use all variable genes
)

# Visualize links
CoveragePlot(
  seurat_obj,
  region = "MS4A1",
  features = "MS4A1",
  expression.assay = "SCT",
  links = TRUE
)

# Get all links
links <- Links(seurat_obj)
```

---

## Part 4: Reference Mapping with Bridge Integration

### 4.1 Query to Reference Mapping (RNA)

```r
library(Azimuth)
library(SeuratData)

# Load reference
# InstallData("pbmcref")
reference <- LoadData("pbmcref")

# Map query to reference
query <- RunAzimuth(
  query = seurat_obj,
  reference = "pbmcref"
)

# Predicted labels transferred
head(query$predicted.celltype.l2)
head(query$predicted.celltype.l2.score)

DimPlot(query, group.by = "predicted.celltype.l2", reduction = "wnn.umap")
```

### 4.2 Bridge Integration (ATAC to RNA Reference)

```r
# For mapping ATAC-only data to RNA reference
library(Seurat)
library(Signac)
library(SeuratData)

# Load ATAC query
atac_query <- LoadData("pbmc_multiome_atac")

# Load multiome bridge (has both modalities)
bridge <- LoadData("pbmc_multiome")

# Load RNA reference
reference <- LoadData("pbmcref")

# Preprocess bridge
DefaultAssay(bridge) <- "RNA"
bridge <- SCTransform(bridge, verbose = FALSE)
bridge <- RunPCA(bridge, verbose = FALSE)

DefaultAssay(bridge) <- "ATAC"
bridge <- RunTFIDF(bridge)
bridge <- FindTopFeatures(bridge, min.cutoff = "q0")
bridge <- RunSVD(bridge)

# Find transfer anchors using bridge
anchors <- FindBridgeTransferAnchors(
  reference = reference,
  bridge = bridge,
  query = atac_query,
  reduction = "lsiproject",  # Project ATAC query into bridge LSI
  dims = 2:30
)

# Transfer labels
atac_query <- MapQuery(
  anchorset = anchors,
  reference = reference,
  query = atac_query,
  refdata = list(
    l1 = "celltype.l1",
    l2 = "celltype.l2"
  ),
  reduction.model = "wnn.umap"
)

DimPlot(atac_query, group.by = "predicted.l2")
```

---

## Part 5: Multi-Sample Integration

### 5.1 Integration Workflow

```r
# For multiple samples/batches
library(Seurat)

# Load multiple samples
sample_list <- list(
  sample1 = Read10X_h5("sample1/filtered_feature_bc_matrix.h5"),
  sample2 = Read10X_h5("sample2/filtered_feature_bc_matrix.h5"),
  sample3 = Read10X_h5("sample3/filtered_feature_bc_matrix.h5")
)

# Create Seurat objects
seurat_list <- lapply(names(sample_list), function(x) {
  CreateSeuratObject(
    counts = sample_list[[x]]$`Gene Expression`,
    project = x,
    min.cells = 3
  )
})
names(seurat_list) <- names(sample_list)

# Add ATAC assay to each
for (name in names(seurat_list)) {
  seurat_list[[name]][["ATAC"]] <- CreateChromatinAssay(
    counts = sample_list[[name]]$Peaks,
    sep = c(":", "-"),
    fragments = paste0(name, "/atac_fragments.tsv.gz")
  )
}

# Merge objects
merged <- merge(
  seurat_list[[1]],
  y = seurat_list[2:length(seurat_list)],
  add.cell.ids = names(seurat_list)
)

# Standard integration using rpca or CCA
# Then proceed with WNN on integrated data
```

---

## Common Pitfalls and Solutions

### Data Preparation
| Pitfall | Solution |
|---------|----------|
| ADT rownames have prefixes | Strip prefixes before analysis |
| Fragment file not indexed | Run `tabix` to index fragment file |
| Wrong genome annotation | Match EnsDb version to your alignment reference |
| Missing peaks for some cells | Check Cell Ranger QC, may need to rerun with lower thresholds |

### WNN Integration
| Pitfall | Solution |
|---------|----------|
| One modality dominates | Check modality weights; adjust `dims.list` |
| LSI component 1 correlates with depth | Exclude component 1 from `dims.list` |
| Poor WNN clustering | Try different resolutions; check per-modality UMAPs first |

### Motif Analysis
| Pitfall | Solution |
|---------|----------|
| ChromVAR fails | Check genome package matches annotation style |
| No motif hits | Peaks may be too narrow; check peak calling |
| Memory error | Subset to fewer peaks or use parallel processing |

### Bridge Integration
| Pitfall | Solution |
|---------|----------|
| Low mapping scores | Check that bridge shares cell types with both query and reference |
| Missing features | Ensure consistent peak set between bridge and query |

---

## Key API Reference

### Assay Management
```r
DefaultAssay(obj) <- "RNA"                    # Switch active assay
Assays(obj)                                    # List all assays
obj[["NewAssay"]] <- CreateAssay5Object()     # Add new assay
GetAssayData(obj, assay = "RNA", layer = "counts")  # Access data
```

### WNN Functions
```r
FindMultiModalNeighbors()   # Build multimodal graph
RunUMAP(nn.name = "weighted.nn")  # Use WNN for UMAP
FindClusters(graph.name = "wsnn")  # Cluster on WNN graph
```

### Signac Functions
```r
CreateChromatinAssay()      # Create ATAC assay
RunTFIDF() + RunSVD()       # Standard ATAC processing
GeneActivity()              # Compute gene activity scores
AddMotifs() + RunChromVAR() # Motif enrichment
LinkPeaks()                 # Peak-gene correlation
CoveragePlot()              # Visualize accessibility
```

### Reference Mapping
```r
RunAzimuth()                        # Automated reference mapping
FindBridgeTransferAnchors()         # Bridge integration anchors
MapQuery()                          # Project query onto reference
```

---

## Resources

- **Seurat v5:** https://satijalab.org/seurat/
- **Signac:** https://stuartlab.org/signac/
- **WNN Tutorial:** https://satijalab.org/seurat/articles/weighted_nearest_neighbor_analysis
- **Multiome Tutorial:** https://satijalab.org/seurat/articles/seurat5_multimodal_vignette
- **Bridge Integration:** https://satijalab.org/seurat/articles/seurat5_integration_bridge
- **ChromVAR:** https://greenleaflab.github.io/chromVAR/
- **JASPAR:** https://jaspar.genereg.net/
