# SCENIC+ Skill

## Overview
SCENIC+ builds enhancer-driven gene regulatory networks (eRegulons: TF→region→gene links) from scRNA-seq + scATAC-seq. Integrates: pycisTopic (topic modeling), pycisTarget (motif enrichment), SCENIC+ (GRN inference).

## Core Concepts

**Topics**: Co-accessible region sets from LDA representing regulatory programs  
**eRegulons**: TF-region-gene modules with GBM importance + correlation scores  
**Cistromes**: TF→region links from motif enrichment (not ChIP-seq)  
**Direct/Extended**: Direct = database motif→TF; Extended = +orthologs/similarity

## Critical Requirements

```python
# CRITICAL: Save raw BEFORE normalization
adata.raw = adata  # Required for SCENIC+ GBM models

# Fragments: sorted, bgzipped, tabix-indexed
fragments_dict = {'sample_id': '/path/fragments.tsv.gz'}

# Barcodes MUST match between RNA/ATAC
# Multiome: exact match
# Separate data: use bc_transform_func or metacells

# cisTarget databases
ctx_db = '*.rankings.feather'  # Recovery curves
dem_db = '*.scores.feather'    # Wilcoxon tests
motif_anno = 'motifs-v10-nr.{species}.tbl'
```

## Workflow (SCENIC+-Specific Steps)

### 1. scRNA-seq Preprocessing
```python
# Standard Scanpy workflow - see tutorials
# ONLY CRITICAL SCENIC+ REQUIREMENT:
adata.raw = adata  # BEFORE normalize_total/log1p/scale
```

### 2. scATAC-seq with pycisTopic
```python
from pycisTopic.cistopic_class import *

# Consensus peaks (recommended over bulk peaks)
consensus_peaks = peak_calling(
    fragments_dict=fragments_dict,
    path_to_blacklist=blacklist_bed,
    n_cpu=10
)

# Create object + run LDA models
cistopic_obj = create_cistopic_object_from_fragments(...)
models = run_cgs_models(
    cistopic_obj,
    n_topics=[20, 25, 30, 35, 40, 45, 50],
    n_iter=500, alpha=50, alpha_by_topic=True, eta=0.1
)

# Select model (use coherence + log-likelihood)
cistopic_obj.add_LDA_model(evaluate_models(models, select_model=40))

# Binarize topics for motif analysis
binarize_topics(cistopic_obj, method='otsu')  # or ntop=3000

# Impute for DARs
impute_accessibility(cistopic_obj, scale_factor=10**6)
markers_dict = find_diff_features(
    cistopic_obj, 
    variable='cell_type',
    adjpval_thr=0.05, log2fc_thr=0.5
)
```

### 3. Custom cisTarget Database (Recommended)
```bash
# Download cluster-buster + motif collection
wget https://resources.aertslab.org/cistarget/programs/cbust
wget https://resources.aertslab.org/.../v10nr_clust_public.zip

# Create FASTA with 1kb background padding
create_fasta_with_padded_bg_from_bed.sh \
    ${GENOME_FASTA} ${CHROMSIZES} ${CONSENSUS_PEAKS_BED} \
    output.fa 1000 yes

# Create databases (~hours, use compute node)
create_cistarget_motif_databases.py \
    -f output.fa \
    -M path/to/motifs/singletons/ \
    -m motifs_list.txt \
    -o output_prefix \
    --bgpadding 1000 -t 20
# Outputs: .rankings.feather + .scores.feather
```
Precomputed databases available at https://resources.aertslab.org/cistarget/ for hg38/mm10/dm6

### 4. SCENIC+ Pipeline via Snakemake
```bash
scenicplus init_snakemake --out_dir scplus_pipeline
# Edit config.yaml with your paths and parameters:
```

**Critical config.yaml parameters:**
```yaml
input_data:
  cisTopic_obj_fname: "cistopic_obj.pkl"
  GEX_anndata_fname: "adata.h5ad"
  region_set_folder: "region_sets/"  # Topics/DARs as .bed
  ctx_db_fname: "*.rankings.feather"
  dem_db_fname: "*.scores.feather"
  path_to_motif_annotations: "motifs-v10-nr.hgnc.tbl"

params_data_preparation:
  bc_transform_func: "lambda x: f'{x}-sample'"  # If barcodes don't match
  is_multiome: True
  species: "hsapiens"
  search_space_upstream: "1000 150000"  # 1kb-150kb from TSS
  search_space_downstream: "1000 150000"

params_motif_enrichment:
  dem_adj_pval_thr: 0.05
  dem_log2fc_thr: 1.0
  ctx_nes_threshold: 3.0
  fraction_overlap: 0.4  # 40% overlap required

params_inference:
  min_target_genes: 10  # Filter low-confidence TFs
  rho_threshold: 0.05   # Min correlation
```

Run: `cd scplus_pipeline/Snakemake/ && snakemake --cores 20`

**Pipeline outputs:**
1. TF/region/gene importance scores (GBM)
2. eRegulons (TF-region-gene links with rankings)
3. AUCell enrichment scores per cell
4. Final MuData: `scplusmdata.h5mu`

### 5. Analyze Results
```python
import mudata
scplus_mdata = mudata.read('outs/scplusmdata.h5mu')

# eRegulons: TF-region-gene links with quality metrics
eregulons = scplus_mdata.uns['direct_e_regulon_metadata']
# Key columns: TF, Gene, Region, importance_TF2G, rho_TF2G, 
#              importance_R2G, rho_R2G, regulation (+1/-1), triplet_rank

# Regulon Specificity Score (cell type specificity)
from scenicplus.RSS import regulon_specificity_scores, plot_rss
rss = regulon_specificity_scores(
    scplus_mudata=scplus_mdata,
    variable='scRNA_counts:cell_type',
    modalities=['direct_gene_based_AUC', 'extended_gene_based_AUC']
)
plot_rss(rss, top_n=3, num_columns=5)

# SCENIC+-specific heatmap: gene AUC (color) + region AUC (size)
from scenicplus.plotting.dotplot import heatmap_dotplot
heatmap_dotplot(
    scplus_mudata=scplus_mdata,
    color_modality='direct_gene_based_AUC',
    size_modality='direct_region_based_AUC',
    group_variable='scRNA_counts:cell_type',
    sort_data_by='direct_gene_based_AUC'
)

# Standard: UMAP on eRegulon activity (concat gene AUCs + use scanpy)
```

## Interpretation & Quality Metrics

**eRegulon Quality:**
- `triplet_rank`: <1000 = high confidence
- `importance_TF2G`: >1.0 = strong TF→gene link (GBM feature importance)
- `|rho|`: >0.3 = strong correlation
- `regulation`: +1 (activator), -1 (repressor)

**RSS (Regulon Specificity):**
- RSS >0.5: highly cell-type specific
- RSS <0.2: broadly active/housekeeping

**Topic Quality:**
- Coherence: higher = better (check with `evaluate_models`)
- Marginal distribution: contribution to model
- Clear nucleosome pattern in insertion size plot = good data

## Integration & Export

**Seurat/Signac**: Export cell-topic matrix as CSV → `CreateDimReducObject()` in R

**ArchR**: Start from ArchR peaks/fragments

**Perturbation simulation**: `train_gene_expression_models()` + `simulate_perturbation()` for TF knockout effects

**Export formats**:
- Cytoscape: eRegulons as TSV (TF, Gene, Region, scores)
- IGV: regions as BED, topics as BigWig
- MuData→old SCENIC+ object: `mudata_to_scenicplus()` for legacy functions

## Performance Tips

- Custom databases: 2-3x better than precomputed
- Binarization: 'ntop' for balanced region sets across topics
- Parallel: Set `n_cpu=20` in all steps
- Batch correction: Use `harmony` in pycisTopic
- Memory: MuData files large; use `backed='r'` mode

## Key Visualizations

**Standard plots** (use scanpy/seaborn): Cell-topic heatmaps, eRegulon UMAP, topic score distributions

**SCENIC+-Specific:**

```python
# 1. Topic-specific TF cistrome comparison (context-dependent activity)
# Shows same TF (e.g., SPI1) has different spatial patterns in B-cells vs monocytes
topic6_regions = eregulons[(eregulons['TF']=='SPI1') & (eregulons['eRegulon_name'].str.contains('Topic6'))]['Region']
topic10_regions = eregulons[(eregulons['TF']=='SPI1') & (eregulons['eRegulon_name'].str.contains('Topic10'))]['Region']
# Score using region AUC, visualize on UMAP

# 2. Gene accessibility scores (when no paired RNA-seq)
# Sum region accessibility around genes to create pseudo-gene expression
region_to_gene = cistopic_obj.region_data[['SYMBOL']].groupby('SYMBOL').apply(lambda x: x.index.tolist()).to_dict()
gene_scores = {gene: imputed_acc[:, regions].sum(1) for gene, regions in region_to_gene.items()}
# Visualize marker genes: CD3D, CD14, MS4A1

# 3. Bulk signature validation
# Compare topics to bulk ATAC-seq from purified cell types using AUC
# Proves topics aren't artifacts

# 4. Export for genome browsers
# BigWig: topic scores per region for IGV/UCSC
# BED: eRegulon regions for locus-specific validation
```

## Quick Checklist

- [ ] RNA: raw counts saved in `.raw` before normalization
- [ ] ATAC: fragments file is bgzipped and tabix-indexed  
- [ ] Barcodes match exactly between RNA and ATAC (or transformed)
- [ ] Cell annotations present in both modalities
- [ ] Custom cisTarget database created (or precomputed downloaded)
- [ ] Motif annotations match species (hgnc/mgi/flybase)
- [ ] Config.yaml paths all absolute and verified
- [ ] Sufficient disk space (~50GB+ for full pipeline)
- [ ] Topics show clear nucleosome pattern in QC
- [ ] eRegulons contain known cell-type markers