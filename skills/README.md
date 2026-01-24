# SciAgent Skills Library

This directory contains modular skill files for single-cell and multi-omics analysis. Each skill provides structured guidance for a specific tool or workflow.

---

## Skill Categories

### scvi-tools Family (Deep Learning Models)

The scvi-tools ecosystem provides probabilistic deep learning models for single-cell analysis. All skills follow a consistent pattern and share installation/configuration via the framework document.

| Skill | Purpose | Modality |
|-------|---------|----------|
| `scvi-framework.md` | **Foundation** - installation, patterns, model selection | All |
| `scvi-basic.md` | Batch correction, DE | scRNA-seq |
| `scvi-scanvi.md` | Semi-supervised label transfer | scRNA-seq + labels |
| `scvi-multivi.md` | Joint RNA+ATAC modeling | Multiome |
| `scvi-peakvi.md` | Accessibility embedding, DA | scATAC-seq |
| `scvi-mrvi.md` | Multi-sample effects, covariate DE | scRNA-seq + samples |
| `scvi-contrastivevi.md` | Salient vs background variation | Perturbation |
| `scvi-linearscvi.md` | Interpretable linear decoder | scRNA-seq |
| `scvi-lda.md` | Topic modeling (AmortizedLDA) | scRNA-seq |
| `scvi-scarches-reference-mapping.md` | Query-to-reference projection | Any scvi model |

**Decision Tree:**
```
scRNA-seq only?
├─ Need label transfer? → scANVI
├─ Multi-sample analysis? → MrVI
├─ Interpretable factors? → LinearSCVI
├─ Topic modeling? → AmortizedLDA
└─ Basic integration → scVI

Multiome (RNA+ATAC)? → MultiVI
scATAC-seq only? → PeakVI
Perturbation study? → ContrastiveVI
Have pretrained reference? → scArches
```

---

### SCENIC+ Ecosystem (Gene Regulatory Networks)

The SCENIC+ pipeline builds enhancer-driven GRNs from paired scRNA-seq + scATAC-seq data.

| Skill | Purpose | Input |
|-------|---------|-------|
| `pycistopic-atac-topic-modeling.md` | Topic models for ATAC regions | scATAC-seq |
| `pycistarget-motif-enrichment.md` | TF motif enrichment in regions | Region sets |
| `scenic-grn-inference.md` | Full eRegulon inference | Both modalities |

**Pipeline Flow:**
```
scATAC-seq → pycisTopic → pycistarget
                  ↓            ↓
             region sets    motif scores
                  ↓            ↓
scRNA-seq → ─────────────────────────────→ SCENIC+ → eRegulons
```

---

### Multimodal Analysis (Paired Data)

| Skill | Purpose | Key Feature |
|-------|---------|-------------|
| `seurat-multimodal-analysis.md` | **R/Seurat** - CITE-seq, 10x Multiome, WNN, Signac, ChromVAR | Full R ecosystem |
| `python-multimodal-10x.md` | **Python** - muon, SnapATAC2, Harmony, scanpy for multiome | Full Python ecosystem |
| `scglue-unpaired-multiomics-integration.md` | Unpaired integration | Guidance graph |
| `treearches-hierarchy-learning.md` | Reference mapping with hierarchy | Tree-structured latent |

**Multiome Decision Tree:**
```
10x Multiome (RNA+ATAC paired)?
├─ R ecosystem preferred? → seurat-multimodal-analysis.md (WNN, Signac)
├─ Python preferred? → python-multimodal-10x.md (muon, SnapATAC2)
└─ Deep learning integration? → scvi-multivi.md

CITE-seq (RNA+Protein)?
├─ R/Seurat → seurat-multimodal-analysis.md (WNN for ADT)
└─ Python → python-multimodal-10x.md (muon)

Unpaired modalities (separate experiments)?
└─ → scglue-unpaired-multiomics-integration.md
```

---

### 10x Genomics Pipelines & QC

| Skill | Purpose |
|-------|---------|
| `cellranger-arc-multiome.md` | Cell Ranger ARC pipeline for 10x Multiome ATAC+GEX |
| `iterative-peak-merging.md` | Iterative summit-preserving peak set construction |

---

### Differential Accessibility & TF Activity

| Skill | Purpose | Key Feature |
|-------|---------|-------------|
| `crescendo-scatac-cre-analysis.md` | CRE-level DA testing | Sub-peak resolution |
| `chromvar-motif-accessibility.md` | TF motif activity from chromatin accessibility | Bias-corrected deviation scores |

---

### Trajectory & RNA Velocity

| Skill | Purpose | Key Tools |
|-------|---------|-----------|
| `rna-velocity-trajectory.md` | Infer cell dynamics from spliced/unspliced RNA | scVelo, VeloVI, Chronocell, CellRank |

**Tool Decision Tree:**
```
Quick exploration? → scVelo steady-state
Publication-quality velocity? → scVelo dynamical
Need uncertainty estimates? → VeloVI
Explicit topology/pseudotime? → Chronocell
Fate probabilities & terminal states? → CellRank
```

**Data Requirements:**
- Spliced/unspliced count quantification (velocyto, kallisto|bustools, alevin-fry)
- Process timescale must match RNA half-life (~hours)

---

### Cell Type Annotation & Transfer Learning

| Skill | Purpose | Method |
|-------|---------|--------|
| `cellxgene-census-annotation.md` | Annotation transfer via CZ CELLxGENE Census (33M+ cells) | KNN in pre-computed embedding space |
| `scembed-atac-annotation.md` | scATAC-seq cell-type annotation via transfer learning | Region embedding + vector DB search |

---

### Data Utilities & Conversion

| Skill | Purpose |
|-------|---------|
| `anndatar-seurat-scanpy-conversion.md` | R/Python data format conversion |
| `louper-seurat-conversion.md` | Seurat to 10x Loupe Browser (.cloupe) |
| `single-cell-vector-search.md` | Semantic search over cell embeddings |
| `annotate-te-rnaseq-data.md` | Transposable element quantification |

---

### Other Analysis

| Skill | Purpose |
|-------|---------|
| `genenmf-metaprogram-discovery.md` | Meta-program discovery via NMF |
| `gatom-metabolomic-predictions.md` | Metabolomics pathway analysis |

---

## scvi-tools Interoperability Framework

All scvi-tools models share a common architecture enabling seamless interoperability:

### Shared Core Pattern

```python
# 1. Register data structure
Model.setup_anndata(adata, layer="counts", batch_key="batch", ...)

# 2. Create and train
model = Model(adata, n_latent=30, ...)
model.train()

# 3. Extract latent representation
adata.obsm["X_model"] = model.get_latent_representation()

# 4. Save/load
model.save("path/", save_anndata=True)
model = Model.load("path/")
```

### Model Compatibility Matrix

| From \ To | scVI | scANVI | MultiVI | PeakVI | MrVI |
|-----------|------|--------|---------|--------|------|
| scVI | - | `from_scvi_model()` | Setup same adata | - | - |
| scANVI | - | - | - | - | - |
| MultiVI | - | - | - | Subset to ATAC | - |
| PeakVI | - | - | Add RNA | - | - |

### scArches Transfer Learning

All scvi-tools models support scArches-style query mapping:

```python
# Reference model (train once)
Model.setup_anndata(ref_adata, ...)
ref_model = Model(ref_adata, use_layer_norm="both", use_batch_norm="none")
ref_model.train()
ref_model.save("ref/")

# Query projection (fast)
Model.prepare_query_anndata(query_adata, "ref/")
query_model = Model.load_query_data(query_adata, "ref/")
query_model.train(max_epochs=200)  # Fine-tune only
```

### Data Flow Conventions

```
Raw counts (adata.X or layer)
         ↓
  setup_anndata() registers structure
         ↓
  Model training (learns latent space)
         ↓
  get_latent_representation() → adata.obsm["X_model"]
         ↓
  scanpy downstream: neighbors → umap → leiden → DE
```

### Key Interoperability Rules

1. **Raw Counts Required**: All models need raw counts, not normalized data
2. **Consistent var_names**: Query data must have EXACT same genes/peaks as reference
3. **Batch Keys**: Must be categorical strings, not integers
4. **Re-register After Subset**: Run `setup_anndata()` again if adata changes
5. **GPU Recommended**: Most models benefit significantly from GPU acceleration

---

## Integration Points

### scvi-tools ↔ SCENIC+

```python
# Use scVI latent for SCENIC+ cell selection
adata.obsm["X_scVI"] = scvi_model.get_latent_representation()
# Export to pycisTopic for ATAC topic modeling
```

### scvi-tools ↔ scGLUE

```python
# scGLUE can use scVI-preprocessed RNA data
rna.obsm["X_scVI"] = scvi_model.get_latent_representation()
# Then align with ATAC via guidance graph
```

### PeakVI ↔ MultiVI

```python
# PeakVI for ATAC-only analysis
peakvi_model = scvi.model.PEAKVI(atac_adata)

# MultiVI for joint RNA+ATAC
multivi_model = scvi.model.MULTIVI(multiome_adata)

# Compare: MultiVI ATAC latent vs PeakVI latent
# Both should recover similar cell populations
```

---

## Adding New Skills

1. Follow naming convention: `tool-name-purpose.md` or `scvi-modelname.md`
2. Include sections: Overview, Installation, Quick Start, Key Parameters, Common Issues
3. Reference related skills (e.g., "See `scvi-framework.md` for installation")
4. Add to appropriate category in this README

---

## Resources

- **scvi-tools**: https://docs.scvi-tools.org/
- **SCENIC+**: https://scenicplus.readthedocs.io/
- **scGLUE**: https://scglue.readthedocs.io/
- **scVelo**: https://scvelo.readthedocs.io/
- **CellRank**: https://cellrank.readthedocs.io/
- **Chronocell**: https://github.com/pachterlab/Chronocell
- **scverse ecosystem**: https://scverse.org/
