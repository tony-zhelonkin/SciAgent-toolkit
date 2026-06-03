# SciAgent Skills Library

This directory contains modular skills for single-cell and multi-omics analysis. Each skill provides structured guidance for a specific tool or workflow.

## Skill Format

Skills use the **canonical Claude Code directory format**:

```
skills/<skill-name>/
├── SKILL.md              # REQUIRED — entry point with YAML frontmatter
├── references/           # Deep-dive topic guides (loaded on demand)
├── scripts/              # Executable Python/Bash/R helpers
├── checks/               # Verification scripts and checklists
└── assets/               # Templates, configs, example data
```

The role activator (`sciagent activate`) resolves skills in this order:
1. **Directory format** (canonical): `skills/<name>/SKILL.md` → symlinks the directory
2. **Flat format** (legacy): `skills/<name>.md` → symlinks the single file

Both coexist during migration. Authoring new skills uses the directory format only.

See `_TEMPLATE/` for a canonical starter and `skill-creator/` for a full reference implementation.

---

## Taxonomy: Scope Classes & `requires:` Inheritance

Every skill's `metadata.scope:` field declares which of three classes it belongs to. The activator (`sciagent activate`) treats them identically — the distinction is purely curatorial, intended to keep individual SKILL.md files focused.

| Scope class | Line cap (target) | Role | Example |
|---|---|---|---|
| `atomic` | ~300 LOC | One tool, one workflow stage. No nested dependencies. | `harmonypy-batch-integration`, `pygenometracks-coverage-plots` |
| `foundation` | ~600 LOC | Shared container / methodology referenced by many other skills. | `anndata`, `scanpy`, `scvi-framework`, `multimodal-anndata-mudata` |
| `orchestrator` | ~250 LOC | Entry point that routes between leaves; contains decision tree + cross-cutting pitfalls but **no per-tool tutorial content** — that lives in the leaves. | `muon-multimodal-analysis`, `seurat-multimodal-analysis`, `scvi-hub-models` |

### `requires:` mechanism (added in ADR-0002)

A skill may declare other skills it transitively depends on via `metadata.requires:`. When a role activates a skill, the resolver walks its `requires:` graph (post-order DFS) and symlinks every transitive dependency into `.claude/skills/` and `.agents/skills/`. Inherited entries appear in the `AGENTS.md` block under a dedicated subsection:

```markdown
## Skills (inherited via requires:)
- `snapatac2-atac-preprocessing` — (via `muon-multimodal-analysis`)
- `harmonypy-batch-integration`   — (via `muon-multimodal-analysis`)
- ...
```

**Worked example** — adding `muon-multimodal-analysis` to a role pulls in:

```yaml
# skills/muon-multimodal-analysis/SKILL.md
metadata:
  scope: orchestrator
  requires:
    - multimodal-anndata-mudata
    - scanpy
    - snapatac2-atac-preprocessing
    - harmonypy-batch-integration
    - atac-differential-accessibility
    - pyranges-peak-gene-linkage
    - pygenometracks-coverage-plots
    - scvi-multivi
```

→ the role YAML stays terse; activating muon brings the whole multimodal leaf set with it. Cyclic and missing dependencies are caught **before** any filesystem mutation (the previous active stack is left intact on resolver failure).

---

## Skill Categories

### scverse Foundation (Core Data Structures & QC)

These skills are prerequisites for all other scverse workflows. Start here if you're new to the ecosystem or working with `.h5ad` files.

| Skill | Purpose |
|-------|---------|
| `anndata.md` | AnnData structure — the shared data container for all scverse tools |
| `scanpy.md` | Standard scRNA-seq QC, clustering, visualization, and DE |
| `single-cell-rna-qc.md` | MAD-based QC filtering — data-driven outlier detection for scRNA-seq |
| `cellranger-multi-to-anndata` | Build a single pooled AnnData from CellRanger Multi per-sample outputs (gene-union concat + biomart symbols + metadata join) |

**Starting a new scRNA-seq project?**
```
New to single-cell or running a standard workflow?
└─ scanpy.md (QC → normalize → UMAP → Leiden → markers)

Need to integrate multiple batches or samples with deep learning?
└─ scvi-framework.md → then scanpy for downstream

Have spliced/unspliced RNA data for trajectory direction?
└─ rna-velocity-trajectory.md (start after scanpy clustering)

Working with format questions (.h5ad reading, subsetting, concat)?
└─ anndata.md
```

---

### scvi-tools Family (Deep Learning Models)

The scvi-tools ecosystem provides probabilistic deep learning models for single-cell analysis. All skills follow a consistent pattern and share installation/configuration via the framework document.

**Architecture (2026-04-20):** `scvi-framework/` is the **Rich-tier router** for this cluster. Shared patterns (`setup_anndata`, batch vs sample keys, scArches flags, hub loading, gotchas, count distributions, GPU, troubleshooting) live in `scvi-framework/references/` and are referenced by every child SKILL.md instead of being duplicated. The cross-reference manifest at `scvi-framework/references/cross-refs.yaml` is the single source of truth for each child's `complementary-skills:` graph. Self-evolution hooks: `scvi-framework/references/troubleshooting.md` (append-only failure log) and `scvi-framework/checks/pre-train-checklist.md` (extend on encounter).

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

### Multimodal Analysis (Paired & Unpaired Data)

The two monolith multimodal skills were dissected in ADR-0002 into orchestrator + leaf form. The orchestrators stay as the entry points; the leaves are pulled in automatically via `requires:`.

#### Orchestrators

| Skill | Purpose | Key Feature |
|-------|---------|-------------|
| `multimodal-anndata-mudata` | **R→Python conversion** for multi-assay Seurat → MuData | Foundation container skill |
| `seurat-multimodal-analysis` | **R/Seurat orchestrator** — chains 4 leaves for CITE-seq, multiome, bridge, unpaired | Decision tree + cross-cutting pitfalls |
| `muon-multimodal-analysis`   | **Python orchestrator** — chains 5 leaves for muon/SnapATAC2/Harmony/DA/peak-gene/coverage | Decision tree + workflow choreography |
| `scglue-unpaired-multiomics-integration` | Unpaired integration + enhancer-gene inference | Guidance graph |
| `treearches-hierarchy-learning` | Reference mapping with hierarchy | Tree-structured latent |

#### Leaves (Python / scverse — pulled in by `muon-multimodal-analysis`)

| Leaf | Purpose |
|---|---|
| `snapatac2-atac-preprocessing` | Fragment import, QC, TF-IDF+LSI, spectral embedding, peak calling, gene activity |
| `harmonypy-batch-integration`  | Multi-sample batch correction on PCA / LSI / spectral embeddings |
| `atac-differential-accessibility` | Marker peaks via scanpy / SnapATAC2 and pseudobulk DA via pyDESeq2 |
| `pyranges-peak-gene-linkage`   | Peak-to-gene cis-regulatory correlations with configurable windows |
| `pygenometracks-coverage-plots`| INI-driven genome-browser figure panels from bigWig / BED / GTF |

#### Leaves (R / Seurat — pulled in by `seurat-multimodal-analysis`)

| Leaf | Purpose |
|---|---|
| `seurat-citeseq-wnn`            | Multi-assay RNA + ADT, CLR + SCTransform, WNN joint embedding |
| `signac-chromatin-analysis`     | 10x Multiome RNA+ATAC, TF-IDF+LSI, gene activity, ChromVAR, LinkPeaks, CoveragePlot |
| `seurat-bridge-integration`     | Azimuth + bridge-integration label transfer (ATAC query → RNA reference) |
| `seurat-unpaired-cross-modality`| CCA anchoring on gene activity for unpaired RNA+ATAC; label transfer + expression imputation |

**Multiome Decision Tree:**
```
10x Multiome (RNA+ATAC paired)?
├─ R ecosystem preferred? → seurat-multimodal-analysis.md (WNN, Signac)
├─ Python preferred? → muon-multimodal-analysis.md (muon, SnapATAC2)
└─ Deep learning integration? → scvi-multivi.md

CITE-seq (RNA+Protein)?
├─ R/Seurat → seurat-multimodal-analysis.md (WNN for ADT)
└─ Python → muon-multimodal-analysis.md (muon)

Unpaired modalities (separate experiments)?
├─ Need to convert Seurat→Python? → multimodal-anndata-mudata.md (conversion recipes)
├─ Integration → scglue-unpaired-multiomics-integration.md
├─ GRN inference → scenic-grn-inference.md (metacell mode)
└─ Stay in R (anchor imputation) → seurat-multimodal-analysis.md (Part 6)
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
| `tf-footprint-differential-analysis.md` | TF footprinting & differential binding | TOBIAS + HINT-ATAC + Signac |

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
| `genenmf-metaprogram-discovery.md` | Meta-program discovery via NMF (per-sample → cross-donor consensus) |
| `consensus-nmf-multirun` | Multi-run consensus cNMF on a single dataset (full ± QC × subset ± QC, merge at r > 0.7, g:Profiler annotation, per-celltype ANOVA) |
| `gatom-metabolomic-predictions.md` | Metabolomics pathway analysis |

---

### scRNA-seq Workflow (project-shared scaffolding)

These three skills share the house style documented by `scrna-pipeline-conventions`. They were authored together as the canonical pipeline shape for scRNA-seq projects matching the user's `scbio-docker` template.

| Skill | Tier | Purpose |
|-------|------|---------|
| `scrna-pipeline-conventions` | simple | House style: numbered scripts (00_build, 01_qc, …), multi-checkpoint, central `config.py` (PATHS / PARAMS), `03_results/{checkpoints,tables,plots,interactive,objects,annotation}/` |
| `cellranger-multi-to-anndata` | standard | Build pooled AnnData from `cellranger multi` outputs |
| `scrna-cxg-host` | rich | Two-phase CellxGene hosting: schema prep + Docker Compose deploy (nginx + htpasswd) |
| `consensus-nmf-multirun` | rich | Multi-run consensus cNMF (run × N variants → K-selection → merge → annotate → ANOVA) |

All four embody the **Decision Pause Contract** (documented per-skill in each SKILL.md's "Decision Pauses" section): each judgment moment surfaces as a `### DECISION PAUSE — <topic>` block that the agent stops at, with named defaults a user can clear with one word ("default"). Recorded choices land in `analysis_config.yaml::decisions::<skill>::<key>` for replay on re-run.

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

### 1. Copy the template

```bash
cp -r skills/_TEMPLATE skills/<your-skill-name>
cd skills/<your-skill-name>
```

Directory name MUST be kebab-case and match the `name:` frontmatter field.

The template ships with `name: SKILL_IDENTIFIER` as a placeholder. That value intentionally fails `quick_validate.py`'s kebab-case check — rename it and every other field should already validate cleanly.

### 2. Fill the SKILL.md frontmatter

The **canonical Claude Code skill schema** only allows these top-level keys:
`name`, `description`, `license`, `metadata`, `allowed-tools`, `compatibility`.

Anything else causes canonical validators (e.g. `skill-creator/scripts/quick_validate.py`,
Claude Desktop plugin loader) to fail. All v3 taxonomy fields must therefore be
**nested under `metadata:`**.

Required (top-level):
- `name` — matches directory name; kebab-case; `^[a-z0-9-]+$`
- `description` — 1-3 sentences, include trigger keywords so the AI activates on intent. If your description contains a raw colon (e.g. `"scVI: integration..."`), wrap it in quotes or use YAML block scalar (`|`).

Recommended (nested under `metadata:`):
- `category` — `foundation` | `integration` | `annotation` | `analysis` | `workflow` | `practice`
- `tier` — `simple` | `standard` | `rich` (see below)
- `tags` — free-form discovery tags (YAML list)
- `complementary-skills` — related skill names (YAML list)
- `contraindications` — when NOT to use this skill (YAML list of strings)
- `version` — semver of the skill itself
- `upstream-docs` — canonical documentation URL
- `skill-author`, `last-reviewed` — provenance

**Canonical frontmatter example:**

```yaml
---
name: my-skill
description: One sentence. Use when <trigger>. For <other case> use other-skill.
license: MIT
metadata:
  skill-author: SciAgent-toolkit
  last-reviewed: 2026-04-14
  version: 1.0.0
  upstream-docs: https://example.org/docs
  category: foundation
  tier: standard
  tags: [tag-a, tag-b]
  complementary-skills:
    - related-skill-a
  contraindications:
    - "Do not use for X. Use other-skill instead."
---
```

Validate with `python skill-creator/scripts/quick_validate.py skills/<your-skill>/`.

### 3. Choose a tier

| Tier | Contents | When |
|---|---|---|
| **Simple** | `SKILL.md` only | One primary use case, few parameters, fits in <200 lines |
| **Standard** | `SKILL.md` + `references/` | Multiple use cases with deep-dive guides |
| **Rich** | `SKILL.md` + `references/` + `scripts/` + `checks/` | Complex workflow with executable verification |

Skills are promoted from Simple → Standard → Rich additively. Start simple; add subdirectories only when needed.

### 4. Required SKILL.md sections (per tier)

| Section | Simple | Standard | Rich |
|---|---|---|---|
| Frontmatter | REQUIRED | REQUIRED | REQUIRED |
| Overview | REQUIRED | REQUIRED | REQUIRED |
| Decision Tree | Optional | REQUIRED | REQUIRED |
| Quick Start (with verification block) | REQUIRED | REQUIRED | REQUIRED |
| Progressive Depth (Basic/Intermediate/Advanced) | Optional | REQUIRED | REQUIRED |
| Verification Checklist | REQUIRED | REQUIRED | REQUIRED |
| Common Pitfalls (Symptom/Cause/Fix) | REQUIRED | REQUIRED | REQUIRED |
| Complementary Skills | Optional | RECOMMENDED | REQUIRED |

### 5. Keep SKILL.md provider-agnostic

The body prose must work for any agentic harness (Claude Code, Gemini CLI, Codex CLI). If a skill genuinely needs provider-specific instructions, factor them into `references/claude.md`, `references/gemini.md`, etc.

### 6. Factor long content into `references/`

SKILL.md should stay under ~400 lines. If a section grows beyond ~300 lines, move it to `references/<topic>.md` and reference from SKILL.md.

### 7. Activate and test

Add the skill name (without `.md`) to the role YAML (`roles/base.yaml` or a new role), then:

```bash
sciagent activate base
ls -la /path/to/project/.claude/skills/<your-skill-name>
# Should be a symlink pointing to 01_modules/SciAgent-toolkit/skills/<your-skill-name>
```

### 8. Add to the category table above

Place your skill in the appropriate category section in this README.

---

## Resources

- **scvi-tools**: https://docs.scvi-tools.org/
- **SCENIC+**: https://scenicplus.readthedocs.io/
- **scGLUE**: https://scglue.readthedocs.io/
- **scVelo**: https://scvelo.readthedocs.io/
- **CellRank**: https://cellrank.readthedocs.io/
- **Chronocell**: https://github.com/pachterlab/Chronocell
- **scverse ecosystem**: https://scverse.org/
