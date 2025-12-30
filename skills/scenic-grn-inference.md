# SCENIC+ GRN Inference Skill

## Overview

SCENIC+ builds enhancer-driven gene regulatory networks (eRegulons: TF→region→gene links) from paired scRNA-seq + scATAC-seq data. Integrates motif enrichment with expression correlation to identify high-confidence regulatory relationships.

## SCENIC+ Ecosystem

```
┌─────────────────┐     ┌─────────────────┐
│   scATAC-seq    │     │   scRNA-seq     │
│   (fragments)   │     │   (counts)      │
└────────┬────────┘     └────────┬────────┘
         │                       │
         ▼                       ▼
┌─────────────────┐     ┌─────────────────┐
│   pycisTopic    │     │  Save .raw      │
│  (topic model)  │     │  BEFORE norm    │
└────────┬────────┘     └────────┬────────┘
         │                       │
         ▼                       │
┌─────────────────┐              │
│   pycistarget   │              │
│ (motif enrich)  │              │
└────────┬────────┘              │
         │                       │
         ▼                       ▼
┌─────────────────────────────────────────┐
│              SCENIC+ Snakemake          │
│  (integrates motifs + expression → GRN) │
└────────────────────┬────────────────────┘
                     │
                     ▼
            ┌─────────────────┐
            │   eRegulons     │
            │ (TF→region→gene)│
            └─────────────────┘
```

**Prerequisites:**
1. **pycisTopic output**: `cistopic_obj.pkl` + `region_sets/*.bed` — see [pycisTopic skill](pycistopic-atac-topic-modeling.md)
2. **pycistarget** (optional standalone): For detailed motif enrichment outside Snakemake — see [pycistarget skill](pycistarget-motif-enrichment.md)
3. **scRNA-seq AnnData**: With `.raw` saved BEFORE normalization

## Core Concepts

**eRegulons**: TF-region-gene modules combining:
- GBM importance scores (TF→gene prediction strength)
- Correlation (rho) between TF and target gene expression
- Motif enrichment linking TF to accessible regions

**Direct vs Extended**:
- Direct = database motif directly maps to TF
- Extended = includes orthologs and motif similarity matches

## Critical Requirements

```python
# CRITICAL: Save raw BEFORE any normalization
adata.raw = adata  # Required for SCENIC+ GBM models

# Barcodes MUST match between RNA and ATAC (for multiome)
# Use bc_transform_func in config if they differ

# cisTarget databases (download from Aertslab)
ctx_db = '*.rankings.feather'  # Recovery curves
dem_db = '*.scores.feather'    # Wilcoxon tests
motif_anno = 'motifs-v10-nr.{species}.tbl'
```

## Data Modes: Multiome vs Non-Multiome (Unpaired)

SCENIC+ supports two data modes:

| Mode | Description | Key Requirement |
|------|-------------|-----------------|
| **Multiome** | Same cells have both RNA + ATAC | Matching barcodes |
| **Non-multiome** | Separate RNA and ATAC experiments | Matching cell-type annotations |

### Non-Multiome / Unpaired Data (Metacell Approach)

For unpaired data (different cells for RNA vs ATAC), SCENIC+ creates **pseudo-multiome metacells**:

> "Sampling a predefined number of cells from each data modality within the same cell-type annotation label and averaging the raw gene expression and imputed chromatin-accessibility data across these cells to create a multiome meta-cell." — [Nature Methods 2023](https://www.nature.com/articles/s41592-023-01938-4)

**Prerequisites for unpaired data:**
1. Good cell-type annotations for BOTH RNA and ATAC
2. Annotations must match between modalities (same label vocabulary)
3. Sufficient cells per cell-type (recommend ≥50 per group)

**Config for unpaired data:**
```yaml
params_data_preparation:
  is_multiome: False                    # CRITICAL: Set to False for unpaired
  key_to_group_by: "cell_type"          # Metadata column for metacell grouping
  nr_cells_per_metacells: 10            # Cells to sample per metacell (default: 10)
```

**Expected output:** ~100-300 pseudo-multiome metacells (depends on cell-type granularity)

## Workflow

### 1. Prepare scRNA-seq Data

```python
import scanpy as sc

# Standard preprocessing
adata = sc.read_h5ad('rna_data.h5ad')
adata.raw = adata  # CRITICAL: Before normalize_total/log1p/scale
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata)
sc.pp.scale(adata)
sc.pp.pca(adata)
sc.pp.neighbors(adata)
sc.tl.umap(adata)

# Save with raw preserved
adata.write('adata_for_scenic.h5ad')
```

### 2. Custom cisTarget Database (Recommended, 2-3x better)

```bash
# Download cbust
wget https://resources.aertslab.org/cistarget/programs/cbust

# Create FASTA from peaks
create_fasta_with_padded_bg_from_bed.sh ${GENOME} ${CHROMSIZES} peaks.bed out.fa 1000 yes

# Build database
create_cistarget_motif_databases.py -f out.fa -M motifs/ -o db_prefix \
    --bgpadding 1000 -t 20
```

**Precomputed databases**: https://resources.aertslab.org/cistarget/ (hg38/mm10/dm6)

### 3. Initialize SCENIC+ Pipeline

```bash
# Initialize Snakemake workflow
scenicplus init_snakemake --out_dir scplus_pipeline
```

### 4. Configure Pipeline

Edit `scplus_pipeline/Snakemake/config.yaml`:

```yaml
input_data:
  cisTopic_obj_fname: "/path/to/cistopic_obj.pkl"      # From pycisTopic
  GEX_anndata_fname: "/path/to/adata_for_scenic.h5ad"  # With .raw saved
  region_set_folder: "/path/to/region_sets/"           # From pycisTopic

params_data_preparation:
  bc_transform_func: "lambda x: f'{x}-sample'"  # If barcodes don't match
  is_multiome: True                              # True for paired/multiome data
  # --- OR for unpaired/non-multiome data: ---
  # is_multiome: False
  # key_to_group_by: "cell_type"                 # Metadata column for metacell grouping
  # nr_cells_per_metacells: 10                   # Cells per metacell (default: 10)
  species: "mmusculus"                           # or "hsapiens"

params_motif_enrichment:
  dem_adj_pval_thr: 0.05     # DEM significance threshold
  dem_log2fc_thr: 1.0        # DEM fold-change threshold
  ctx_auc_threshold: 0.005   # cisTarget AUC threshold
  ctx_nes_threshold: 3.0     # cisTarget NES threshold (context specificity)
  ctx_rank_threshold: 0.05   # cisTarget rank threshold
  # For detailed motif enrichment parameters, see pycistarget skill

cistarget_databases:
  ctx_db: "/path/to/*.rankings.feather"
  dem_db: "/path/to/*.scores.feather"
  motif_annotation: "/path/to/motifs-v10-nr.mmus-m0.001-o0.0.tbl"
```

### 5. Run Pipeline

```bash
cd scplus_pipeline/Snakemake/
snakemake --cores 20
```

### 6. Analyze Results

```python
import mudata
from scenicplus.RSS import regulon_specificity_scores
from scenicplus.plotting.dotplot import heatmap_dotplot

# Load results
scplus_mdata = mudata.read('outs/scplusmdata.h5mu')

# eRegulons metadata
eregulons = scplus_mdata.uns['direct_e_regulon_metadata']
# Columns: TF, Gene, Region, importance_TF2G, rho_TF2G, regulation, triplet_rank

# Regulon Specificity Score (RSS)
rss = regulon_specificity_scores(scplus_mdata, variable='scRNA_counts:cell_type',
    modalities=['direct_gene_based_AUC'])

# Visualization: gene AUC (color) + region AUC (size)
heatmap_dotplot(scplus_mdata, color_modality='direct_gene_based_AUC',
    size_modality='direct_region_based_AUC', group_variable='scRNA_counts:cell_type')
```

## Quality Metrics

| Metric | Good | Interpretation |
|--------|------|----------------|
| triplet_rank | <1000 | High-confidence TF-region-gene link |
| importance_TF2G | >1.0 | Strong GBM feature importance |
| \|rho\| | >0.3 | Strong TF-gene correlation |
| RSS | >0.5 | Cell-type specific regulon |
| N eRegulons | >50 | Sufficient regulatory coverage |

## Filtering eRegulons

```python
# High-confidence eRegulons
hc_eregulons = eregulons[
    (eregulons['triplet_rank'] < 500) &
    (eregulons['importance_TF2G'] > 1.5) &
    (abs(eregulons['rho_TF2G']) > 0.3)
]

# Cell-type specific (high RSS)
tf_rss = rss.loc[rss.max(axis=1) > 0.5]
```

## Standalone pycistarget vs Snakemake

| Approach | When to Use |
|----------|-------------|
| **SCENIC+ Snakemake** | Full GRN workflow, standard parameters |
| **Standalone pycistarget** | Custom motif analysis, parameter tuning, exploratory cistrome analysis |

For detailed motif enrichment outside the Snakemake pipeline (cisTarget, DEM, cistrome extraction), see [pycistarget skill](pycistarget-motif-enrichment.md).

## Performance Tips

- **Custom databases**: 2-3x better motif recovery than precomputed
- **Memory**: Use `backed='r'` mode for large MuData files
- **Parallel**: Set `--cores 20` in Snakemake
- **Species**: Must match between annotation, databases, and data

## Troubleshooting

| Issue | Solution |
|-------|----------|
| Barcode mismatch | Set `bc_transform_func` in config |
| No eRegulons found | Check `.raw` was saved before normalization |
| Memory errors | Use custom database subset or increase RAM |
| Low motif recovery | Build custom cisTarget database (see [pycistarget skill](pycistarget-motif-enrichment.md)) |
| Unpaired data fails | Set `is_multiome: False` + `key_to_group_by` |
| Too few metacells | Increase `nr_cells_per_metacells` or use coarser annotations |
| Metacell mismatch | Ensure RNA and ATAC annotations use identical label vocabulary |

## Quick Checklist

**For all data:**
- [ ] pycisTopic completed: `cistopic_obj.pkl` + `region_sets/` ready
- [ ] RNA: `.raw` saved BEFORE normalize_total/log1p
- [ ] cisTarget database matches genome version + species
- [ ] Config.yaml paths are absolute
- [ ] Sufficient eRegulons (>50) with good metrics

**For multiome (paired) data:**
- [ ] Barcodes match (or bc_transform_func configured)
- [ ] `is_multiome: True` in config

**For non-multiome (unpaired) data:**
- [ ] `is_multiome: False` in config
- [ ] `key_to_group_by` set to cell-type annotation column
- [ ] RNA and ATAC annotations use identical label vocabulary
- [ ] ≥50 cells per group for robust metacells
