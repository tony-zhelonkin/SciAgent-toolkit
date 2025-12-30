# CREscendo Skill

## Purpose
Analyze differential CRE usage in scATAC-seq data via Tn5 cleavage frequency chi-square tests.

## Biology Problem
```
Peak-based DA:        |-------- PEAK --------|  → "No difference" (similar total accessibility)
CREscendo:            |CRE1|CRE2|gap|CRE3|CRE4|  → CRE3 active in cDC1A, CRE4 in cDC1B ✓
```

**Peak-level methods miss cell-type-specific signals when CREs within the same peak have opposite patterns.**

## Installation
```r
devtools::install_github("ChenMengjie/CREscendo")
library(CREscendo)
library(Signac)
library(Seurat)
```

## Prerequisites

**Fragment files**: Must be bgzip-compressed + tabix-indexed
```bash
# If not indexed:
bgzip fragments.tsv && tabix -p bed fragments.tsv.gz
```

**From Signac/Seurat object**:
```r
# Extract peaks as GRanges
peaks_gr <- granges(seurat_obj@assays$peaks)
peaks_gr <- peaks_gr[grep("chr", seqnames(peaks_gr))]  # Remove scaffolds

# Get fragment file path (already in Signac object)
frag_path <- Fragments(seurat_obj)[[1]]@path
```

## Core Workflow (~1hr total)

### 1. Load CRE Annotations (ENCODE Registry v3)
```r
# Human: 1,063,878 CREs | Mouse: 368,121 CREs
CRE_gr <- CREscendo:::load_CRE_annotations(genome = "hg38")  # or "mm10"
```

### 2. Prepare Cell Annotations
```r
# Two-column data.frame: cellname (barcode), celltype (label)
cells <- colnames(seurat_obj@assays$peaks)
cell_type_list <- c("cDC1A", "cDC1B", "cDC2")  # Cell types to compare

cell_annotation <- data.frame(
  celltype = seurat_obj@active.ident[seurat_obj@active.ident %in% cell_type_list],
  cellname = cells[seurat_obj@active.ident %in% cell_type_list]
)
```

### 3. Open TabixFile Connection
```r
tabix_file <- Rsamtools::TabixFile(
  file = "fragments.tsv.gz",
  index = "fragments.tsv.gz.tbi"
)
open(tabix_file)
```

### 4. Create CREscendo Object (~25-30 min)
```r
# Human
cre_obj <- PrepareCREscendoObject(
  peaks_gr, tabix_file, cell_annotation, cell_type_list, CRE_gr,
  tssRegion = c(-1000, 1000),
  TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene,
  annoDb = "org.Hs.eg.db",
  multiCREs = TRUE  # Include peaks with multiple CREs only
)

# Mouse
cre_obj <- PrepareCREscendoObject(
  peaks_gr, tabix_file, cell_annotation, cell_type_list, CRE_gr,
  tssRegion = c(-1000, 1000),
  TxDb = TxDb.Mmusculus.UCSC.mm10.knownGene,
  annoDb = "org.Mm.eg.db"
)
```

### 5. Chi-Square Test (~12-25 min)
```r
# Compare two cell types
results <- CREtest(cre_obj, celltypes = c("cDC1A", "cDC1B"))

# Or compare one vs all others
results <- CREtest(cre_obj, celltypes = c("cDC1A"))
```

### 6. Filter & Summarize
```r
summary <- Summarize(results)  # Full results

# Filter high-confidence results
filtered <- Filter(results,
  fdr_cutoff = 0.05,
  coverage_cutoff = 100,   # Min fragments per peak
  abs_diff_cutoff = 0.2    # Min frequency difference
)
```

### 7. Visualize
```r
# By peak ID
Visualize(results, which_peak = "chr1-207335620-207337449")

# By gene symbol
Visualize(results, which_gene = "CD55")

# Extract detailed info + contribution plot
info <- Extract(results, which_gene = "CD55")
info$peak_annotation    # All peaks near gene
info$plot_contribution  # CRE contribution to chi-square
```

## Key Output Columns

| Column | Interpretation |
|--------|----------------|
| `chisq_stat` | Overall differential usage (>2000 = very strong) |
| `fdr` | FDR-adjusted p-value |
| `max_contribution` | % of chi-sq from top CRE |
| `max_refID` | Index of top-contributing CRE |
| `max_dif` | Frequency difference in top CRE |
| `abs_dif` | Absolute frequency difference |
| `n_1`, `n_2` | Fragment counts per cell type |
| `CRE_type` | ENCODE annotations (dELS, pELS, PLS, CTCF) |
| `SYMBOL` | Gene symbol |

## Chi-Square Interpretation

| chisq_stat | Meaning |
|------------|---------|
| >2000 | Very strong differential CRE usage |
| 500-2000 | Moderate differential usage |
| <500 | Weak or no differential usage |

**max_contribution >80%**: Single CRE drives signal (clear cell-type-specific enhancer)
**max_contribution ~50%**: Multiple CREs contribute (complex regulation)

## CRE Type Annotations

| Type | ENCODE Category |
|------|-----------------|
| `PLS` | Promoter-like signature (H3K4me3 high) |
| `pELS` | Proximal enhancer-like (<2kb from TSS) |
| `dELS` | Distal enhancer-like (>2kb from TSS) |
| `CTCF-bound` | CTCF insulator binding |

## Comparison with FindMarkers (DA)

```r
# Standard DA test (for comparison)
da_peaks <- FindMarkers(seurat_obj,
  ident.1 = "cDC1A", ident.2 = "cDC1B",
  test.use = 'LR', latent.vars = 'nCount_peaks'
)  # ~30-60 min

# Compare results
commons <- intersect(rownames(summary), rownames(da_peaks))
combined <- cbind(summary[commons,], da_peaks[commons,])

# Peaks with low DA log2FC but high chi-square = CRE-level signals missed by DA
missed_by_da <- combined[abs(combined$avg_log2FC) < 1 & combined$chisq_stat > 1000, ]
```

## Troubleshooting

| Issue | Solution |
|-------|----------|
| No CREs found | Check genome version matches (hg38/mm10) |
| Empty results | Verify cell_annotation format (celltype, cellname columns) |
| TabixFile error | Ensure .tbi index exists and path is correct |
| Memory error | Reduce cell_type_list, process in chunks |
| Few significant hits | Lower coverage_cutoff, check cell type sample sizes |

## Quick Checklist

- [ ] Fragment files: bgzip + tabix indexed
- [ ] Genome version matches peaks (hg38/mm10)
- [ ] cell_annotation: 2 columns (celltype, cellname)
- [ ] Sufficient cells per type (recommend >100)
- [ ] TxDb + annoDb match genome
- [ ] CRE annotations loaded for correct genome

## Advantages Over Peak-Based DA

1. **Base-resolution**: Analyzes CRE-level patterns within peaks
2. **Recovers hidden signals**: Detects opposing CRE usage that cancels at peak level
3. **Interpretable**: Maps to ENCODE regulatory elements
4. **Reproducible**: Standardized CRE coordinates across studies