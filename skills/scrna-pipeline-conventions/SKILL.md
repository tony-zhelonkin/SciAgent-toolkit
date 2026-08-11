---
name: scrna-pipeline-conventions
description: "House style for scRNA-seq projects: numbered analysis stages, multi-checkpoint discipline in long-running stages, a central config sourced from analysis_config.yaml, the 03_results/<stage>/{figures,tables}/ layout, and per-utility helper modules. Use when starting a project that should match the established workflow, or retrofitting one."
license: MIT
---

# scRNA-seq Pipeline Conventions

A small set of patterns that make scRNA-seq pipelines restartable, debuggable, and reproducible across projects. These are the shared scaffolding the other three workflow skills (`cellranger-multi-to-anndata`, `scrna-cxg-host`, `consensus-nmf-multirun`) reference rather than re-derive. Adopt them when a project will run for more than a week or be revisited by a future session.

**When to use this skill:**
- Starting a new scRNA-seq project and want the established directory + script layout
- Retrofitting an existing project to a shared house style before adding new pipeline stages
- Author a SKILL.md or analysis script that must point to a stable convention rather than redefine one

**When NOT to use this skill:**
- Throwaway exploratory notebooks → adopt none of this; just hack
- Software-design discipline (ADRs, design reviews) → use `architecture-first-dev`
- Lint/style enforcement → out of scope; this is a convention, not a checker

---

## Conventions

### 1. Numbered scripts

Analysis scripts live at the project root (or `02_analysis/`) with a two-digit numeric prefix that encodes pipeline order:

```
00_build_anndata.py          # ingestion (CellRanger Multi → AnnData)
01_qc_anndata.py             # MAD QC + cell cycle + doublets + ambient
02_annotation.py             # AUCell + CellTypist + leiden
03_cxg_prepare.py            # CXG schema preparation
04_subcluster.py             # per-celltype HVG/PCA/UMAP/leiden re-embed
05_compute_scores.py         # signature scoring, pseudo-bulk DE
06_run_cnmf.py               # cNMF on full dataset
07_merge_programs.py         # cross-run consolidation
08_analyze_programs.py       # condition association
09_visualize.py              # publication figures
```

Variants (a, b, c suffixes) hold parallel-equivalent scripts that differ only in input slice or QC variant — for example `06b_cnmf_qc_filtered.py`, `06c_cnmf_subset_X.py`, `06d_cnmf_subset_Y.py`. The base number stays the stage; the letter is the variant.

**Why.** The prefix is the topological order. Anyone reading the directory listing knows the pipeline shape in five seconds. `git diff` and file-tree views sort meaningfully. Scripts are restartable per stage instead of buried in import graphs.

**How to apply.** New stage gets the next free number; new variant of an existing stage gets the next free letter. If you find yourself wanting `00.5_*`, that is signal the new step is structurally distinct and deserves its own number — re-number downstream files (rare; cheap on `git mv`).

---

### 2. Multi-checkpoint within long-running scripts

A long script that runs multiple substantive stages writes a checkpoint after each stage. AnnData checkpoints (`.h5ad`) live under `03_results/objects/` (or `03_results/objects/<stage>/`), keyed by stage name:

```
03_results/objects/00_raw.h5ad            # after ingestion
03_results/objects/01_qc.h5ad             # after MAD QC + filter
03_results/objects/01_qc_cycle.h5ad       # after cell-cycle scoring
03_results/objects/01_qc_doublets.h5ad    # after scrublet
03_results/objects/01_qc_ambient.h5ad     # after SoupX / CellBender
```

Every script reads from one checkpoint and writes to the next. The pattern at the top of the script:

```python
import scanpy as sc
from config import PATHS

# Resume point — read the latest checkpoint produced by the previous stage
adata = sc.read_h5ad(PATHS.object("00_raw"))
# ... work ...
adata.write_h5ad(PATHS.object("01_qc"), compression="gzip")
```

**Why.** A two-hour run that fails on step 4 should not lose steps 1–3. Each checkpoint is a Pareto-honest moment to inspect intermediate state. Sessions resume from the latest checkpoint, never from raw.

**How to apply.** Any script that takes more than ~10 minutes to run, or executes more than two substantive stages, splits its output into multiple checkpoints. Name them with the stage prefix the script implements (`01_qc.h5ad`, not `qc_step1.h5ad`). Place all `.h5ad` / `.rds` checkpoints under `objects/`; this directory is gitignored and never deleted between runs.

---

### 3. Central `config.py` sourced from `analysis_config.yaml`

A single `config.py` at the project root imports parameters from `02_analysis/config/analysis_config.yaml` and exposes two top-level frozen dataclasses:

- **`PATHS`** — path resolution. Properties (not attributes) so directories are auto-created on first access. Resolves stage dirs via `stage_dir(id)` / `figures(stage)` / `tables(stage)` rather than flat paths.
- **`PARAMS`** — thresholds and parameters. All numeric defaults (DE FDR cutoff, MAD scale, K range for cNMF, batch key, etc.) live here.

```python
# config.py — sketch
import yaml
from dataclasses import dataclass
from pathlib import Path

_cfg = yaml.safe_load(open("02_analysis/config/analysis_config.yaml"))

@dataclass(frozen=True)
class _Paths:
    root: Path = Path(_cfg["paths"]["root"])

    def _results(self) -> Path:
        return self.root / "03_results"

    def stage_dir(self, stage_id: str) -> Path:
        p = self._results() / stage_id
        p.mkdir(parents=True, exist_ok=True)
        return p

    def figures(self, stage_id: str) -> Path:
        p = self.stage_dir(stage_id) / "figures"
        p.mkdir(parents=True, exist_ok=True)
        return p

    def tables(self, stage_id: str) -> Path:
        p = self.stage_dir(stage_id) / "tables"
        p.mkdir(parents=True, exist_ok=True)
        return p

    @property
    def objects(self) -> Path:
        p = self._results() / "objects"
        p.mkdir(parents=True, exist_ok=True)
        return p

    def object(self, name: str) -> Path:
        return self.objects / f"{name}.h5ad"

    @property
    def master(self) -> Path:
        p = self._results() / "master"
        p.mkdir(parents=True, exist_ok=True)
        return p

    @property
    def interactive(self) -> Path:
        p = self._results() / "interactive"
        p.mkdir(parents=True, exist_ok=True)
        return p

    @property
    def scratch(self) -> Path:
        p = self._results() / "_scratch"
        p.mkdir(parents=True, exist_ok=True)
        return p

PATHS = _Paths()
PARAMS = _cfg["thresholds"]   # plain dict is fine; freeze if you prefer
```

**Why.** Changing a threshold is one edit. Switching between container and host filesystems is one edit. Skills downstream (e.g., `consensus-nmf-multirun`) read `analysis_config.yaml::decisions::*` for replay-from-config; the convention guarantees that file exists and has a known shape.

**How to apply.** First script in a project creates `analysis_config.yaml` and `config.py`. Every subsequent script imports `from config import PATHS, PARAMS` and never hard-codes a path or threshold.

---

### 4. Output directory layout `03_results/`

The canonical layout is stage-based, matching `analysis_config.yaml::stages:` and the CRAFT results-placement rule. Stages are registered in `analysis_config.yaml` before any script writes to them.

```
03_results/
├── <stage>/                 # one dir per stages: entry (e.g. 01_qc, 02_eda, 05_annotation)
│   ├── figures/
│   │   ├── _overview/       # cross-contrast figures (+ same-stem source tables under tables/_overview/)
│   │   └── by_contrast/<c>/ # per-contrast figures
│   ├── tables/
│   │   ├── _overview/
│   │   └── by_contrast/<c>/
│   └── README.md            # captions (how-to-read) for THIS stage's artifacts
├── objects/                 # finalized milestone deliverables (e.g. <project>_annotated.h5ad); source of truth — .rds / .cloupe / CXG derive from this
├── checkpoints/             # intermediate stage .h5ad files keyed by stage (gitignored, never deleted)
├── master/                  # cross-stage accumulator tables (append_master_table)
├── interactive/             # HTML dashboards
└── _scratch/                # sanctioned ephemeral zone (gitignored)
```

Figures, source-table adjacency, and captions follow the figure-style contract (see the `figure-style` skill) and the CRAFT results-placement rule — use `save_overview()`; do not hand-build `03_results/` paths.

Lower-case `03_results/` (not `03_Results/`) is the convention. Existing projects with capital-R `03_Results/` may keep it; new projects use lower-case.

**Why.** Each stage owns its artifacts. Anyone walking into the project sees what kind of output is where — per stage — without a flat namespace collision. `objects/` and `interactive/` at the results root are project-wide resources, not tied to a single stage.

**How to apply.** `PATHS` (Convention 3) auto-creates these directories. A new analysis output goes into the stage directory by purpose; never at the project root. Register the stage in `analysis_config.yaml::stages:` first.

#### `objects/` vs `checkpoints/`

`objects/` holds **finalized milestone deliverables** — a single file named `<project>_annotated.h5ad`, not stage-prefixed — from which all derived formats are generated. Cell order is guaranteed identical across derived formats because they all originate from this file.

`checkpoints/` holds **intermediate stage files** keyed by stage (`01_qc.h5ad`, `02_embedding.h5ad`, etc.); consumed by the next stage and never delivered externally.

Derivation chain: `<project>_annotated.h5ad` → `.rds` (via `anndatar-seurat-scanpy-conversion`) → `.cloupe` (via `louper-seurat-conversion`) → CXG (via `scrna-cxg-host`).

**Rule:** Change the `.h5ad` first; regenerate downstream formats from it. Never patch a derived format directly — the next regeneration overwrites the patch.

**Evidence:** `03_results/objects/README.md`; `02_analysis/scripts/08a_finalize_annotated_h5ad.py:1–9` docstring.

---

### 5. Per-utility custom modules

A `Python_scripts/` (or `src/`) directory holds focused utility modules with clear single-purpose names:

```
Python_scripts/
├── anndata_utils.py     # read_10x_h5_robust, robust_concat_with_gene_union
├── cxg_utils.py         # CXG schema validation, obsm-DataFrame → array
├── geneset_utils.py     # signature loading, Ensembl ↔ symbol mapping
├── biomart_utils.py     # cached biomart lookups
└── soupx_wrapper.py     # Python ↔ R SoupX bridge
```

Each module has ≤300 lines and covers one concern. Functions exported from these modules are imported by the numbered scripts; logic that recurs in two scripts is moved into a module.

**Why.** Numbered scripts stay procedural and readable end-to-end. Recurring logic is named, testable, and reusable. The skill `cellranger-multi-to-anndata` ships a `build_anndata.py` whose helpers came from this layer; `scrna-cxg-host` ships `prepare_for_cxg.py` from `cxg_utils.py`. The convention makes skills portable.

**How to apply.** First time a piece of logic is needed in a numbered script, it can live inline. Second time, extract into a `*_utils.py` module. Name modules by *what they touch* (anndata, cxg, geneset, biomart), not by the stage that called them first.

---

### 6. Milestone gene naming — symbols active, Ensembl preserved

Build scripts keep Ensembl IDs as `var_names` throughout (stable identifiers, no duplicates, safe for biomart lookups and gene-set matching). At the **packaging milestone** the convention flips: `var_names` becomes gene symbols (made unique with `-1`/`-2` suffixes), and Ensembl is preserved in `var['gene_id']`.

```python
a.var["gene_id"] = a.var_names.copy()          # preserve before switch
a.var_names      = a.var["gene_name"].values    # symbols become active
a.var_names_make_unique(join="-")               # Symbol-1, Symbol-2 for dups
a.uns["var_names_note"] = (
    "var_names are gene symbols (make-unique). Ensembl IDs preserved in "
    "var['gene_id']. To switch back: adata.var_names = adata.var['gene_id']."
)
```

R/Seurat equivalent at packaging: `rownames(obj) <- obj[["RNA"]][[]]$gene_name`; `gene_id` stays in feature metadata.

**Why.** Collaborators query `adata[:, ['Gzmk','Cd3e']]` and gene-set tools use symbols. Ensembl stays as the machine-readable backup. Flipping at build time would silently break intermediate scripts that join on Ensembl; flipping once at the finalize script is safe and discoverable (the `uns` note documents the round-trip).

**How to apply.** The switch happens ONCE, in the finalize script (e.g., `08a_finalize_annotated_h5ad.py`). Never flip mid-pipeline. The `analysis_config.yaml` key `var_name_strategy: ensembl` governs build stages only; packaging always delivers symbols. Confirm symbols carry through when exporting to `.rds` or `.cloupe` — see `anndatar-seurat-scanpy-conversion` and `louper-seurat-conversion`.

**Evidence:** `02_analysis/scripts/08a_finalize_annotated_h5ad.py:104–113`; `02_analysis/config/analysis_config.yaml:78`.

---

### 7. Dual-embedding deliverable — unsupervised + integrated

When a pipeline produces both a de-novo unsupervised embedding (state discovery) and a supervised/batch-corrected integrated embedding (label transfer), carry **both** in the packaged object, renamed at packaging:

| Key | Basis | Primary question |
|-----|-------|-----------------|
| `X_umap_unsupervised` | `X_pca` | What **cell states** exist? (sub-states, continua) |
| `X_umap_integrated` | `X_scANVI` | What **cell types** are these? (label transfer, condition comparison) |

```python
# verify before writing — all cells including query/treated
for emb in ("X_umap_unsupervised", "X_umap_integrated"):
    assert not np.isnan(a.obsm[emb]).any(), f"NaN in {emb}"

a.uns["embeddings"] = {
    "X_umap_unsupervised": "de novo PCA->UMAP; cell-STATE discovery",
    "X_umap_integrated":   "scANVI-corrected UMAP; cell-TYPE comparison across conditions",
}
```

Carry the naming into derived formats: Seurat `.rds` as `umap.unsup`/`umap.integrated`; Loupe as `umap_unsupervised`/`umap_integrated` (see `louper-seurat-conversion`).

**Why.** The supervised latent space is optimized to separate coarse training labels, which compresses within-label substructure. The unsupervised space keeps all variance so sub-states stay resolved. Empirical basis: GZMK⁺ effector kNN-purity is 0.817 in `X_pca` vs 0.767 in `X_scANVI` (base 0.213) — the population is present and distinguishable in both, but crisper in the unsupervised space. Dropping either embedding destroys one axis of interpretability.

**How to apply.** Do NOT drop the unsupervised embedding after annotation. Assert all cells — including query or treated cohort — carry valid (non-NaN) coordinates in both embeddings before writing. Document both in `uns['embeddings']`.

**Evidence:** `02_analysis/scripts/08a_finalize_annotated_h5ad.py:130–147`; `docs/_internal/reasoning/2026-06-27_gzmk_oldspace_validation_and_dual_embedding.md`.

---

### 8. Milestone packaging pass — one deliberate obs cleanup

Do **not** prune `obs` columns mid-pipeline. Intermediate stages join on original column names; a rename in script `02c` breaks `02d`. Instead, perform ONE cleanup pass at the finalize script:

1. **Dedup identity columns** to canonical `sample_id` / `mouse_id` / `pool`; drop aliases
2. **Drop derivable QC columns** (`log1p_*`, `outlier_*`); keep base metrics + `qc_pass`
3. **Move single-valued columns** (`Project`, `Organ`, `Sex`) to `uns`
4. **Apply collaborator-agreed naming** here and nowhere else (e.g., `Group == "Young"` → `age_group == "Adult"`; condition label `"Y_C"` → `"Adult ctrl"`)

The delivered obs schema must be documented in `objects/README.md`.

**Why.** Each stage safely joins on original names. Cleanup is correct only at the known-final point; it also halves the in-memory obs footprint for the delivered object.

**How to apply.** A single `cleanup_obs(adata)` call at the top of the finalize script. Never rename obs columns in numbered scripts before the finalize stage. If a downstream step needs a renamed column, add an alias column — do not rename the source.

**Evidence:** `02_analysis/scripts/08a_finalize_annotated_h5ad.py:62–99`; `docs/_internal/reasoning/2026-06-27_metadata_naming_and_cleanup.md`.

---

## Practical adoption checklist

- [ ] Project root has `config.py` and `02_analysis/config/analysis_config.yaml`
- [ ] `analysis_config.yaml` has a `stages:` block; each stage entry exists before any script writes to it
- [ ] `PATHS.stage_dir(id)`, `PATHS.figures(id)`, `PATHS.tables(id)` resolve and auto-create directories
- [ ] `PATHS.objects`, `PATHS.master`, `PATHS.interactive`, `PATHS.scratch` are the four root-level resource dirs
- [ ] First analysis script is named `00_<verb>_<noun>.py` (e.g., `00_build_anndata.py`)
- [ ] Long-running scripts (≥10 min, ≥2 substantive stages) write multiple `.h5ad` checkpoints under `objects/`
- [ ] Figures and tables use `save_overview()` (from the `figure-style` skill); no hand-built `03_results/` paths
- [ ] Recurring logic lives in `Python_scripts/<topic>_utils.py`, not at the top of an analysis script
- [ ] `analysis_config.yaml` has a `decisions:` section (empty is fine; skills append to it)

---

## DECISION PAUSE — Adoption scope

> **Stop here.** Present the options below to the user and wait for their selection before proceeding. Do not pick a default silently.

**Question for the user:** This project already has analysis files. How should the conventions be adopted?

| Option | What happens | When to choose | Default |
|--------|--------------|-----------------|---------|
| **A — Reference only** | Conventions are documented in `CLAUDE.md` / `AGENTS.md`; no project files are created or renamed | Established project; conventions apply to new files only | yes (safest) |
| **B — Scaffold new directories** | Stage dirs, `objects/`, `master/`, `interactive/`, `_scratch/` created; `config.py` + `analysis_config.yaml` template added; existing scripts left alone | New project; want the scaffolding without touching pre-existing analysis files | |
| **C — Full retrofit** | Skill proposes renames for existing scripts to match `NN_verb_noun.py`, reorganises outputs into the stage-based layout, and produces a single batch-rename script the user reviews before running | Project owner explicitly wants the retrofit and accepts the disruption | |

**After the user chooses:** proceed with the chosen option. Append the choice to `analysis_config.yaml` under `decisions.scrna-pipeline-conventions.adoption_scope` so re-runs replay the decision.

For a fresh project (no analysis files yet), Option B is the implicit default and the pause may be skipped silently — but only when the directory truly contains no `*.py` analysis files.

---

### Decision Pause anti-patterns to avoid

- Picking a default silently because "it's faster"
- Trying one option without naming the tradeoff to the user
- Skipping the pause because a previous session's choice is in config — but the input directory has changed since then (e.g., new scripts appeared)
- Promising end-to-end automation while a pause is upcoming

---

## Verification Checklist

After adopting the conventions, confirm:

- [ ] **`PATHS` resolves and creates directories.** Running `from config import PATHS; print(PATHS.stage_dir("01_qc"))` prints an absolute path and the directory now exists on disk.
- [ ] **`PARAMS` loads from YAML.** `from config import PARAMS; PARAMS["de_fdr"]` returns the expected numeric without `KeyError`.
- [ ] **At least one numbered script reads from a checkpoint and writes another.** `grep -l "PATHS.object" 02_analysis/*.py` returns ≥1 file (or analogous for `Python_scripts/`).
- [ ] **`analysis_config.yaml` has a `decisions:` section.** Even if empty (`decisions: {}`), the key exists for downstream skills to append to.

---

## Common Pitfalls

### Pitfall: PATHS as plain attributes instead of properties

- **Symptom:** Running a fresh `git clone` of the project, the first script crashes with `FileNotFoundError: 03_results/objects/`.
- **Cause:** `PATHS.objects` was defined as a `Path` attribute computed once at import time, before the directory existed. The auto-create behaviour requires `@property`.
- **Fix:** Convert each path to a `@property` that calls `mkdir(parents=True, exist_ok=True)` before returning.

### Pitfall: Hard-coded paths inside numbered scripts

- **Symptom:** Switching from container to host filesystem requires editing five scripts; one gets missed and silently writes to the wrong place.
- **Cause:** A script wrote `pd.read_csv("<project-root>/03_results/01_qc/tables/x.csv")` directly instead of `pd.read_csv(PATHS.tables("01_qc") / "x.csv")`.
- **Fix:** `grep -nE "/scratch|/data|03_results/" 02_analysis/*.py` — every match should be inside `config.py` only.

### Pitfall: Numbered script "00.5"

- **Symptom:** A new step is added between `00_build_anndata.py` and `01_qc_anndata.py`; the author hesitates to renumber and creates `00b_*` or `00.5_*`.
- **Cause:** Confusion between *stage* (the integer) and *variant* (the letter). `00b_*` should mean "alternate variant of stage 0", not "an interstitial new stage".
- **Fix:** If the new step is genuinely a new stage, renumber downstream scripts (`git mv`). If it is a variant of an existing stage, use the letter. If it is one-off exploration, it does not belong as a numbered script — put it in a notebook.

### Pitfall: `Python_scripts/` becomes a junk drawer

- **Symptom:** `Python_scripts/utils.py` grows to 1500 lines mixing AnnData, CXG, biomart, and plotting helpers.
- **Cause:** Logic was extracted into a single `utils.py` instead of a topic-named module.
- **Fix:** Split by *what it touches*. `anndata_utils.py`, `cxg_utils.py`, `biomart_utils.py`. Each ≤300 lines.

### Pitfall: Capital-R `03_Results/` clashes with new lower-case convention

- **Symptom:** A skill writes to `03_results/` but the project already has `03_Results/` from a prior workflow; outputs land in two places.
- **Cause:** The prior reference used capital-R; the convention is lower-case (per project `CLAUDE.md`).
- **Fix:** New projects use lower-case from day one. Existing projects keep their capital-R; do not migrate mid-stream. `PATHS` reads the actual directory name from `analysis_config.yaml::paths::results_dir` so a project can opt into either.

---

## Appendix: Annotation consensus — flag, don't force

When fusing multiple label-transfer methods (CellTypist, scANVI, Seurat), bias toward **flagging** uncertain cells rather than coercing them into the nearest coarse label.

```yaml
# analysis_config.yaml — annotation block
annotation:
  confidence_min: 0.5          # Seurat prediction.score.max floor; below this, abstain
  consensus_policy: conservative  # 'conservative' | 'balanced' | 'permissive'
  treearches_threshold: 0.7    # scHPL posterior threshold
```

Use `consensus_policy: conservative` whenever the perturbation may induce novel cell states. A weak `prediction.score.max` below `confidence_min` should **abstain** (`Uncertain`) or propagate `Novel:<subcluster>` rather than absorbing the cell into the nearest training label.

Report the confidently-unplaceable fraction per condition as a QC deliverable. Reference recall on known coarse types is **not** the success metric — it only confirms the classifier reproduces its training categories. The guard's job is to protect genuinely distinct treated-condition states from silent absorption.

**Evidence:** `docs/_internal/reasoning/2026-06-27_rare-type-reframe-flag-not-force.md`; `02_analysis/scripts/02f_consensus_finalize.py:12–19`; `02_analysis/config/analysis_config.yaml:167–170`.

---

## Resources

This SKILL.md is the canonical reference. Other skills point here. There is no upstream documentation — these conventions are local to this skill library.

---

## When not to use

- Do not use for software-design discipline (ADRs, design reviews, refactor planning). Use architecture-first-dev.
- Do not enforce on throwaway exploration notebooks or one-off scratch analyses. The conventions are for production analysis pipelines that will run more than once.
- Do not use as a code-review checklist. The skill documents conventions; lint/style enforcement is out of scope.

---

## See also

- `cellranger-multi-to-anndata` — Next step (Stage 0); build the first AnnData from CellRanger output
- `scrna-cxg-host` — Downstream consumer; host the prepared `.h5ad` on CellxGene for the wet lab (uses `PATHS.objects`)
- `consensus-nmf-multirun` — Downstream consumer; discover gene programs after annotation (uses `PATHS.objects`, `PATHS.tables`)
- `architecture-first-dev` — Different scope; discipline for software-design changes (ADRs, reviews), not data-pipeline scaffolding
- `anndatar-seurat-scanpy-conversion` — Downstream of milestone packaging; export the packaged `.h5ad` to Seurat `.rds` (Convention 6, 7), preserves `gene_id` and both embeddings
- `louper-seurat-conversion` — Export to Loupe `.cloupe` (Convention 7); carries `umap_unsupervised`/`umap_integrated` embedding names
- `single-cell-rna-qc` — Adjacent; run MAD-based QC filtering on the produced AnnData (consumes `00_raw.h5ad` from this convention's Stage 0)
- `figure-style` — Governs all `03_results/<stage>/figures/` output; this skill defers to it
