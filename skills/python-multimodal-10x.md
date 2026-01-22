# Python Multimodal 10x Analysis Skill

## Overview

Comprehensive skill for multimodal single-cell analysis using the Python/scverse ecosystem. Covers 10x Multiome (RNA+ATAC), CITE-seq (RNA+Protein), multi-sample integration, and downstream chromatin analysis.

**Key Tools:**
- **muon/MuData**: Multimodal data container (analogous to Seurat multi-assay)
- **SnapATAC2**: Fast and scalable scATAC-seq processing
- **scanpy**: Core single-cell analysis toolkit
- **scvi-tools**: Deep learning models (see `scvi-multivi.md` for detailed MultiVI usage)
- **Harmony**: Multi-sample batch correction
- **episcanpy**: Epigenomic single-cell analysis

**Key Features:**
- MuData multimodal object management
- SnapATAC2 for efficient ATAC-seq preprocessing
- TF-IDF + LSI for ATAC dimensionality reduction
- Harmony integration across samples
- Differential accessibility testing
- Joint RNA+ATAC visualization

**Documentation:** https://muon.scverse.org/

---

## Installation

```bash
# Core tools
pip install muon mudata scanpy

# ATAC-seq analysis
pip install snapatac2
pip install episcanpy  # Alternative for ATAC

# Multi-sample integration
pip install harmonypy

# Deep learning (optional)
pip install scvi-tools

# Visualization
pip install matplotlib seaborn plotly
```

---

## Part 1: MuData Multimodal Container

### 1.1 Understanding MuData Structure

```python
import muon as mu
import mudata as md
import scanpy as sc
import anndata as ad
import numpy as np

# MuData is a container for multiple AnnData objects (modalities)
# Similar to multi-assay Seurat object

# Structure:
# mdata.mod["rna"]  -> AnnData for RNA
# mdata.mod["atac"] -> AnnData for ATAC
# mdata.obs         -> Cell metadata (shared)
# mdata.var         -> Feature metadata (combined)
# mdata.obsm        -> Cell embeddings (shared)
```

### 1.2 Create MuData from AnnData Objects

```python
# Load individual modalities
rna = sc.read_h5ad("rna.h5ad")
atac = sc.read_h5ad("atac.h5ad")

# Create MuData
mdata = md.MuData({"rna": rna, "atac": atac})

# Access modalities
print(mdata.mod["rna"])
print(mdata.mod["atac"])

# Shared observations (intersection by default)
print(mdata.n_obs)  # Number of cells with both modalities
print(mdata.obs)    # Shared cell metadata

# Update shared metadata
mdata.update()
```

### 1.3 Load 10x Multiome Data

```python
import muon as mu
import scanpy as sc

# Read 10x Cell Ranger ARC output
mdata = mu.read_10x_h5("filtered_feature_bc_matrix.h5")
# Returns MuData with 'rna' and 'atac' modalities

# Or read individual modalities
rna = sc.read_10x_mtx("filtered_feature_bc_matrix/", gex_only=True)
atac = mu.atac.read_10x_mtx("filtered_feature_bc_matrix/")

mdata = md.MuData({"rna": rna, "atac": atac})

# Add fragment file for ATAC
mdata.mod["atac"].uns["files"] = {"fragments": "atac_fragments.tsv.gz"}
```

### 1.4 MuData I/O

```python
# Save
mdata.write("multiome.h5mu")

# Load
mdata = mu.read("multiome.h5mu")

# Export modality to AnnData
rna = mdata.mod["rna"].copy()
rna.write_h5ad("rna_only.h5ad")

# Subset MuData
mdata_subset = mdata[mdata.obs["cell_type"] == "B_cell"].copy()
```

---

## Part 2: RNA Preprocessing with Scanpy

### 2.1 Standard RNA QC and Normalization

```python
import scanpy as sc

rna = mdata.mod["rna"]

# QC metrics
rna.var["mt"] = rna.var_names.str.startswith("MT-")  # or "mt-" for mouse
sc.pp.calculate_qc_metrics(
    rna,
    qc_vars=["mt"],
    percent_top=None,
    log1p=False,
    inplace=True
)

# Visualize QC
sc.pl.violin(rna, ["n_genes_by_counts", "total_counts", "pct_counts_mt"])

# Filter cells
sc.pp.filter_cells(rna, min_genes=200)
sc.pp.filter_genes(rna, min_cells=3)
rna = rna[rna.obs["pct_counts_mt"] < 10].copy()

# Store raw counts
rna.layers["counts"] = rna.X.copy()

# Normalize
sc.pp.normalize_total(rna, target_sum=1e4)
sc.pp.log1p(rna)

# Highly variable genes
sc.pp.highly_variable_genes(rna, n_top_genes=3000, flavor="seurat_v3",
                             layer="counts")  # Use raw counts for HVG

# Scale and PCA
rna.raw = rna  # Store normalized data
sc.pp.scale(rna, max_value=10)
sc.pp.pca(rna, n_comps=50)
```

### 2.2 RNA Clustering and Visualization

```python
# Neighbors and UMAP
sc.pp.neighbors(rna, n_pcs=30)
sc.tl.umap(rna)
sc.tl.leiden(rna, resolution=0.5)

# Visualize
sc.pl.umap(rna, color=["leiden", "n_genes_by_counts"])
```

---

## Part 3: ATAC Preprocessing with SnapATAC2

### 3.1 SnapATAC2 Quick Start

```python
import snapatac2 as snap

# Read fragment file directly
adata_atac = snap.pp.import_data(
    fragment_file="atac_fragments.tsv.gz",
    chrom_sizes=snap.genome.hg38,  # or mm10
    sorted_by_barcode=False,
    min_num_fragments=500,
    max_num_fragments=100000,
)

# Or read existing AnnData
adata_atac = mdata.mod["atac"]
```

### 3.2 ATAC QC Metrics

```python
# Calculate QC metrics
snap.metrics.tsse(adata_atac, snap.genome.hg38)  # TSS enrichment
snap.metrics.frag_size_distr(adata_atac)  # Fragment size distribution

# Visualize QC
snap.pl.frag_size_distr(adata_atac)

# Filter cells
snap.pp.filter_cells(adata_atac, min_counts=1000, max_counts=100000,
                     min_tsse=2)  # TSS enrichment > 2
```

### 3.3 Peak Calling and Matrix

```python
# Call peaks with MACS2 (requires macs2 installed)
snap.pp.make_peak_matrix(
    adata_atac,
    file="peaks.h5ad",  # Output file
    peak_file=None,  # Call peaks de novo
    use_rep="X_spectral"  # Use spectral embedding for grouping
)

# Or use existing peaks
snap.pp.make_peak_matrix(
    adata_atac,
    file="peaks.h5ad",
    peak_file="peaks.bed"  # BED file with peak coordinates
)
```

### 3.4 TF-IDF and LSI (Alternative to SnapATAC2)

```python
import episcanpy as epi
import scanpy as sc

atac = mdata.mod["atac"]

# Binary matrix (presence/absence)
atac.X = (atac.X > 0).astype(float)

# TF-IDF normalization
epi.pp.tfidf(atac)

# Or manual TF-IDF
from sklearn.feature_extraction.text import TfidfTransformer
tfidf = TfidfTransformer(norm="l2", use_idf=True)
atac.X = tfidf.fit_transform(atac.X)

# LSI (Latent Semantic Indexing) - essentially TruncatedSVD
from sklearn.decomposition import TruncatedSVD
lsi = TruncatedSVD(n_components=50, random_state=0)
atac.obsm["X_lsi"] = lsi.fit_transform(atac.X)

# Note: LSI component 1 often correlates with sequencing depth
# Check correlation and exclude if necessary
import pandas as pd
corr = pd.Series(atac.obsm["X_lsi"][:, 0]).corr(
    pd.Series(np.array(atac.X.sum(axis=1)).flatten())
)
print(f"LSI1 correlation with depth: {corr:.3f}")
# If |corr| > 0.5, use X_lsi[:, 1:] for downstream
```

### 3.5 SnapATAC2 Spectral Embedding

```python
# Feature selection
snap.pp.select_features(adata_atac, n_features=50000)

# Spectral embedding (preferred for ATAC)
snap.pp.spectral(adata_atac, n_comps=50)

# UMAP
snap.pp.umap(adata_atac)

# Clustering
snap.pp.knn(adata_atac)
snap.tl.leiden(adata_atac, resolution=0.5)

# Visualize
snap.pl.umap(adata_atac, color="leiden")
```

---

## Part 4: Multi-Sample Integration with Harmony

### 4.1 Harmony Integration for RNA

```python
import scanpy as sc
import scanpy.external as sce

# Assuming rna has 'sample' column in obs
rna = mdata.mod["rna"]

# Standard preprocessing
sc.pp.normalize_total(rna, target_sum=1e4)
sc.pp.log1p(rna)
sc.pp.highly_variable_genes(rna, n_top_genes=3000)
rna.raw = rna
sc.pp.scale(rna, max_value=10)
sc.pp.pca(rna, n_comps=50)

# Harmony integration
sce.pp.harmony_integrate(
    rna,
    key="sample",  # Batch key
    basis="X_pca",  # Input embedding
    adjusted_basis="X_pca_harmony"  # Output embedding
)

# Use harmonized PCA for downstream
sc.pp.neighbors(rna, use_rep="X_pca_harmony")
sc.tl.umap(rna)
sc.tl.leiden(rna)

# Check integration
sc.pl.umap(rna, color=["sample", "leiden"])
```

### 4.2 Harmony Integration for ATAC

```python
import snapatac2 as snap
import scanpy.external as sce

atac = mdata.mod["atac"]

# SnapATAC2 preprocessing
snap.pp.select_features(atac, n_features=50000)
snap.pp.spectral(atac, n_comps=50)

# Harmony on spectral embedding
sce.pp.harmony_integrate(
    atac,
    key="sample",
    basis="X_spectral",
    adjusted_basis="X_spectral_harmony"
)

# Downstream
snap.pp.umap(atac, use_rep="X_spectral_harmony")
snap.pp.knn(atac, use_rep="X_spectral_harmony")
snap.tl.leiden(atac)
```

### 4.3 Joint Integration (WNN-style)

```python
import muon as mu
import numpy as np

# Method 1: Simple concatenation of embeddings
rna = mdata.mod["rna"]
atac = mdata.mod["atac"]

# Align cells (intersection)
common_cells = list(set(rna.obs_names) & set(atac.obs_names))
rna_sub = rna[common_cells].copy()
atac_sub = atac[common_cells].copy()

# Concatenate embeddings
joint_embedding = np.concatenate([
    rna_sub.obsm["X_pca"][:, :30],    # First 30 RNA PCs
    atac_sub.obsm["X_lsi"][:, 1:31]   # LSI 2-31 (skip depth-correlated)
], axis=1)

rna_sub.obsm["X_joint"] = joint_embedding

# UMAP on joint embedding
sc.pp.neighbors(rna_sub, use_rep="X_joint")
sc.tl.umap(rna_sub)

# Method 2: muon's WNN implementation
mu.pp.neighbors(mdata, key_added="wnn")
mu.tl.umap(mdata, neighbors_key="wnn")
```

---

## Part 5: Differential Accessibility Analysis

### 5.1 Marker Peaks with Scanpy

```python
import scanpy as sc

atac = mdata.mod["atac"]

# Ensure raw counts available
if "counts" in atac.layers:
    atac.X = atac.layers["counts"].copy()

# Normalize for DE testing
sc.pp.normalize_total(atac, target_sum=1e4)
sc.pp.log1p(atac)

# Find marker peaks per cluster
sc.tl.rank_genes_groups(
    atac,
    groupby="leiden",
    method="wilcoxon",  # or "t-test", "logreg"
    pts=True  # Compute fraction of cells expressing
)

# Get results
markers = sc.get.rank_genes_groups_df(atac, group="0")  # Cluster 0
markers = markers[markers["pvals_adj"] < 0.05]
print(markers.head(20))

# Visualize
sc.pl.rank_genes_groups(atac, n_genes=10)
sc.pl.rank_genes_groups_dotplot(atac, n_genes=5)
```

### 5.2 Differential Accessibility with SnapATAC2

```python
import snapatac2 as snap

atac = mdata.mod["atac"]

# Marker peaks (one-vs-rest)
snap.tl.marker_regions(
    atac,
    groupby="leiden",
    pvalue=0.01
)

# Access results
markers = atac.uns["marker_regions"]

# Pairwise comparison
da_results = snap.tl.diff_test(
    atac,
    groupby="condition",
    groups=["treated", "control"],
    method="t-test"  # or "wilcoxon"
)
```

### 5.3 Pseudobulk Differential Accessibility

```python
import scanpy as sc
import numpy as np
import pandas as pd
from scipy import stats

atac = mdata.mod["atac"]

# Create pseudobulk by summing counts per sample/condition
def create_pseudobulk(adata, groupby):
    """Create pseudobulk matrix"""
    groups = adata.obs[groupby].unique()
    pseudobulk = {}
    for g in groups:
        mask = adata.obs[groupby] == g
        pseudobulk[g] = np.array(adata[mask].X.sum(axis=0)).flatten()
    return pd.DataFrame(pseudobulk, index=adata.var_names)

pb = create_pseudobulk(atac, groupby="sample")

# Then use DESeq2 or edgeR in R, or pyDESeq2 in Python
# pip install pydeseq2
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

# Prepare metadata
metadata = pd.DataFrame({
    "sample": pb.columns,
    "condition": ["treated", "treated", "control", "control"]  # Example
}).set_index("sample")

dds = DeseqDataSet(
    counts=pb.T,
    metadata=metadata,
    design_factors="condition"
)
dds.deseq2()

stat_res = DeseqStats(dds, contrast=["condition", "treated", "control"])
stat_res.summary()
results = stat_res.results_df
```

---

## Part 6: Gene Activity and Peak-Gene Links

### 6.1 Gene Activity Scores

```python
import snapatac2 as snap
import scanpy as sc

atac = mdata.mod["atac"]

# SnapATAC2 gene activity
gene_matrix = snap.pp.make_gene_matrix(
    atac,
    gene_anno=snap.genome.hg38,  # or custom GTF
    upstream=2000,
    downstream=0
)

# Add as new modality
mdata.mod["gene_activity"] = gene_matrix

# Or with muon
import muon as mu
mu.atac.pp.gene_activity(mdata, gene_anno=snap.genome.hg38)
```

### 6.2 Peak-Gene Correlation

```python
import numpy as np
import pandas as pd
from scipy import stats

rna = mdata.mod["rna"]
atac = mdata.mod["atac"]

# Get common cells
common = list(set(rna.obs_names) & set(atac.obs_names))
rna_sub = rna[common]
atac_sub = atac[common]

def correlate_peak_gene(peak_idx, gene_idx):
    """Correlate peak accessibility with gene expression"""
    peak_vals = np.array(atac_sub.X[:, peak_idx].todense()).flatten()
    gene_vals = np.array(rna_sub.X[:, gene_idx].todense()).flatten()
    r, p = stats.pearsonr(peak_vals, gene_vals)
    return r, p

# Example: correlate all peaks within 100kb of a gene
gene = "MS4A1"
gene_idx = list(rna.var_names).index(gene)

# Get peaks near gene (requires genomic coordinates in atac.var)
# atac.var should have 'chrom', 'start', 'end' columns
gene_chrom = "chr11"  # MS4A1 location
gene_start = 60223282

nearby_peaks = atac.var[
    (atac.var["chrom"] == gene_chrom) &
    (abs(atac.var["start"] - gene_start) < 100000)
]

correlations = []
for peak in nearby_peaks.index:
    peak_idx = list(atac.var_names).index(peak)
    r, p = correlate_peak_gene(peak_idx, gene_idx)
    correlations.append({"peak": peak, "gene": gene, "r": r, "p": p})

corr_df = pd.DataFrame(correlations)
corr_df = corr_df[corr_df["p"] < 0.05].sort_values("r", ascending=False)
```

---

## Part 7: Visualization

### 7.1 Joint UMAP Visualization

```python
import muon as mu
import matplotlib.pyplot as plt

# Plot both modalities
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

sc.pl.umap(mdata.mod["rna"], color="cell_type", ax=axes[0], show=False)
axes[0].set_title("RNA")

sc.pl.umap(mdata.mod["atac"], color="cell_type", ax=axes[1], show=False)
axes[1].set_title("ATAC")

plt.tight_layout()
plt.show()
```

### 7.2 Coverage Plots (pyGenomeTracks)

```python
# Install: pip install pygenometracks

# Create tracks.ini file
tracks_config = """
[bigwig]
file = coverage.bw
title = ATAC Coverage
height = 3
color = blue

[genes]
file = genes.gtf
title = Genes
fontsize = 10
height = 5
"""

with open("tracks.ini", "w") as f:
    f.write(tracks_config)

# Run from command line:
# pyGenomeTracks --tracks tracks.ini --region chr11:60000000-60500000 -o coverage.png
```

### 7.3 Heatmaps

```python
import scanpy as sc
import matplotlib.pyplot as plt

# RNA expression heatmap
sc.pl.heatmap(
    mdata.mod["rna"],
    var_names=["CD3D", "CD4", "CD8A", "MS4A1", "CD14"],
    groupby="cell_type",
    standard_scale="var"
)

# ATAC accessibility heatmap (top DA peaks)
top_peaks = markers.head(50)["names"].tolist()
sc.pl.heatmap(
    mdata.mod["atac"],
    var_names=top_peaks,
    groupby="cell_type",
    standard_scale="var"
)
```

---

## Part 8: Integration with scvi-tools

For deep learning-based multimodal integration, see `scvi-multivi.md`.

```python
import scvi

# Quick MultiVI setup
adata = prepare_multiome_for_multivi(mdata)  # Combine modalities
scvi.model.MULTIVI.setup_anndata(
    adata,
    layer="counts",
    batch_key="sample",
    modality_key="modality"
)

model = scvi.model.MULTIVI(adata)
model.train()

adata.obsm["X_multivi"] = model.get_latent_representation()
```

---

## Common Pitfalls and Solutions

### MuData
| Pitfall | Solution |
|---------|----------|
| Cell name mismatch between modalities | Use `mdata.update()` or manually align obs_names |
| Memory issues with large data | Use backed mode: `mu.read("file.h5mu", backed="r")` |
| Lost metadata after subsetting | Call `mdata.update()` after operations |

### ATAC Preprocessing
| Pitfall | Solution |
|---------|----------|
| LSI component 1 correlated with depth | Exclude first component from downstream |
| Low TSS enrichment | Check fragment file quality; adjust thresholds |
| Too few/many peaks | Adjust peak calling parameters or use different peak set |
| SnapATAC2 memory error | Use disk-backed mode or process in chunks |

### Integration
| Pitfall | Solution |
|---------|----------|
| Harmony overcorrection | Reduce `theta` parameter |
| Poor batch mixing | Increase number of PCs/LSI components |
| Cell type mixing | Use more conservative resolution |

### Differential Accessibility
| Pitfall | Solution |
|---------|----------|
| No significant peaks | Use more lenient thresholds; check normalization |
| Batch effects confounding | Include batch as covariate or use pseudobulk |
| Sparse data issues | Aggregate into larger genomic windows |

---

## Key API Reference

### muon
```python
mu.read_10x_h5()              # Read 10x multiome
mu.read() / mdata.write()     # I/O
mu.pp.neighbors()             # Joint neighbors
mu.tl.umap()                  # Joint UMAP
mu.atac.pp.gene_activity()    # Gene activity scores
```

### SnapATAC2
```python
snap.pp.import_data()         # Import fragments
snap.pp.filter_cells()        # QC filtering
snap.pp.select_features()     # Feature selection
snap.pp.spectral()            # Spectral embedding
snap.pp.make_peak_matrix()    # Peak matrix
snap.tl.marker_regions()      # Find markers
snap.tl.diff_test()           # Differential test
```

### scanpy
```python
sc.pp.normalize_total()       # Normalization
sc.pp.highly_variable_genes() # Feature selection
sc.pp.pca()                   # PCA
sc.pp.neighbors()             # KNN graph
sc.tl.leiden()                # Clustering
sc.tl.rank_genes_groups()     # DE/DA testing
```

### harmonypy
```python
sce.pp.harmony_integrate()    # Harmony batch correction
```

---

## Resources

- **muon:** https://muon.scverse.org/
- **SnapATAC2:** https://kzhang.org/SnapATAC2/
- **scanpy:** https://scanpy.readthedocs.io/
- **episcanpy:** https://episcanpy.readthedocs.io/
- **Harmony:** https://github.com/immunogenomics/harmony
- **scvi-tools MultiVI:** https://docs.scvi-tools.org/en/stable/user_guide/models/multivi.html
- **scverse ecosystem:** https://scverse.org/
- **pyGenomeTracks:** https://pygenometracks.readthedocs.io/
