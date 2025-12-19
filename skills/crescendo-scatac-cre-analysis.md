# CREscendo Skill

## Purpose
Analyze differential usage of cis-regulatory elements (CREs) in scATAC-seq data using Tn5 cleavage frequencies.

## Biology Context
- **CREs**: Enhancers, promoters, insulators that control gene expression
- **Problem**: Peak-based methods miss cell-type-specific regulatory signals because peaks span multiple CREs
- **Solution**: Segment peaks by ENCODE CRE annotations, test differential Tn5 cleavage patterns between cell types

## Installation
```r
devtools::install_github("ChenMengjie/CREscendo")
library(CREscendo)
```

## Core Workflow
```r
# 1. Test differential CRE usage (chi-square test on cleavage frequencies)
results <- CREscendo_test(
  peak_data,           # Peak regions containing multiple CREs
  fragments,           # Fragment files with Tn5 cleavage sites
  cre_annotations,     # ENCODE CRE coordinates
  cell_type_A,         # First cell type
  cell_type_B          # Second cell type
)

# 2. Filter high-confidence results
filtered <- Filter(
  results,
  fdr_cutoff = 0.05,          # FDR adjusted p-value
  coverage_cutoff = 100,       # Min fragments per peak
  abs_diff_cutoff = 0.2        # Min absolute difference in CRE frequencies
)

# 3. Visualize specific peaks/genes
Visualize(
  results,
  gene = "CD248",              # Or specify peak directly
  show_encode = TRUE           # Include ENCODE annotations
)
```

## Key Outputs

### Test Statistics
- **chi_squared**: Overall differential usage (higher = stronger cell-type difference)
- **p_value**: Significance of differential usage
- **fdr**: FDR-adjusted p-value
- **partial_chi_squared**: Per-CRE contribution to signal

### CRE Frequency Interpretation
- **Frequency vector**: Proportion of cleavage sites in each CRE segment per cell type
- **High freq in CRE**: Active regulatory element in that cell type
- **Differential freq**: Cell-type-specific enhancer/promoter activity

## Visualization Best Practices
1. **Track plot**: Show Tn5 cleavage counts across peak region by cell type
2. **CRE segments**: Overlay ENCODE CRE annotations
3. **Frequency barplot**: Compare CRE usage proportions between cell types
4. **ENCODE validation**: Include bulk ATAC/ChIP tracks if available

## Interpretation Logic

**Peak-level vs CRE-level**:
- Peak-level DA (e.g., Signac): Can miss signals when overall accessibility similar but specific CREs differ
- CREscendo: Detects when CRE3 active in CD8+ but CRE4 active in CD14+ within same peak

**Chi-square > 2000**: Very strong differential CRE usage
**Chi-square 500-2000**: Moderate differential usage
**Partial chi-square**: Identifies which specific CRE(s) drive the signal

**Biological meaning**:
- High CRE activity in cell type → likely enhancer/promoter for genes in that cell type
- CTCF-bound CREs → insulators regulating chromatin structure
- Cell-type-specific CRE patterns → regulatory mechanisms controlling cell identity

## Common Patterns

**Example 1 - CD248 gene**:
- CRE3 highly active in CD8+ T cells → T cell-specific enhancer
- CRE4 highly active in CD14+ monocytes → monocyte-specific enhancer
- Overall peak accessibility similar → missed by peak-based methods

**Example 2 - VAV1 gene**:
- CRE9 accessible in CD14+ monocytes, linked to platelet-related GWAS SNPs
- Shows monocyte-platelet interaction regulation

## Data Requirements
- scATAC-seq fragment files (Tn5 cleavage positions)
- Peak calls from MACS2, Cell Ranger, or MACS3
- ENCODE CRE annotations (Registry v3: 1M+ human CREs, 368K+ mouse CREs)
- Cell type labels for cells

## Advantages Over Peak-Based Methods
1. **Precision**: Base-resolution CRE analysis vs broad peaks
2. **Interpretability**: Identifies specific regulatory elements, not just accessible regions
3. **Reproducibility**: Standardized ENCODE CRE coordinates across studies
4. **Sensitivity**: Recovers signals when differential usage within peak cancels out at peak level