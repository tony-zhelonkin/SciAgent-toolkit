---
name: scrna-pipeline-conventions
description: "scrna-pipeline-conventions — house style for scRNA-seq analysis projects: numbered analysis scripts (00_build, 01_qc, 02_annotation, ...), multi-checkpoint discipline within long-running scripts, central config.py with PATHS and PARAMS dataclasses sourced from analysis_config.yaml, output-directory layout 03_results/[checkpoints,tables,plots,interactive,objects,annotation]/, and per-utility custom modules (anndata_utils, cxg_utils, geneset_utils). Use when starting a new scRNA-seq project that should match the established workflow shared by cellranger-multi-to-anndata, scrna-cxg-host, and consensus-nmf-multirun, or when retrofitting an existing project to the convention. Practice-tier skill — no executable code; documents the style other skills reference. For software-design discipline (ADRs, design reviews) use architecture-first-dev. For MAD-based QC filtering use single-cell-rna-qc."
license: MIT
metadata:
  skill-author: SciAgent-toolkit
  last-reviewed: 2026-04-29
  category: practice
  tier: simple
  version: 0.1.0
  upstream-docs: ""
  tags:
    - conventions
    - project-scaffold
    - reproducibility
    - checkpoints
    - config
    - scrna-seq
    - workflow-style
    - house-style
  complementary-skills:
    - cellranger-multi-to-anndata
    - scrna-cxg-host
    - consensus-nmf-multirun
    - architecture-first-dev
  contraindications:
    - "Do not use for software-design discipline (ADRs, design reviews, refactor planning). Use architecture-first-dev."
    - "Do not enforce on throwaway exploration notebooks or one-off scratch analyses. The conventions are for production analysis pipelines that will run more than once."
    - "Do not use as a code-review checklist. The skill documents conventions; lint/style enforcement is out of scope."
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

## The five conventions

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

A long script that runs multiple substantive stages writes a checkpoint after each stage. For example, `01_qc_anndata.py` produces:

```
03_results/checkpoints/01_qc.h5ad           # after MAD QC + filter
03_results/checkpoints/02_qc_cycle.h5ad     # after cell-cycle scoring
03_results/checkpoints/03_qc_doublets.h5ad  # after scrublet
03_results/checkpoints/04_qc_ambient.h5ad   # after SoupX / CellBender
```

Every script reads from one checkpoint and writes to the next. The pattern at the top of the script:

```python
import scanpy as sc
from config import PATHS

# Resume point — read the latest checkpoint produced by the previous stage
adata = sc.read_h5ad(PATHS.checkpoint("00_raw"))
# ... work ...
adata.write_h5ad(PATHS.checkpoint("01_qc"), compression="gzip")
```

**Why.** A two-hour run that fails on step 4 should not lose steps 1–3. Each checkpoint is a Pareto-honest moment to inspect intermediate state. Sessions resume from the latest checkpoint, never from raw.

**How to apply.** Any script that takes more than ~10 minutes to run, or executes more than two substantive stages, splits its output into multiple checkpoints. Name them with the stage prefix the script implements (`01_qc.h5ad`, not `qc_step1.h5ad`).

---

### 3. Central `config.py` sourced from `analysis_config.yaml`

A single `config.py` at the project root imports parameters from `02_analysis/config/analysis_config.yaml` and exposes two top-level frozen dataclasses:

- **`PATHS`** — path resolution. Properties (not attributes) so directories are auto-created on first access. Top-level keys are `raw`, `checkpoints`, `tables`, `plots`, `interactive`, `objects`, `annotation`.
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

    @property
    def checkpoints(self) -> Path:
        p = self.root / "03_results" / "checkpoints"
        p.mkdir(parents=True, exist_ok=True)
        return p

    def checkpoint(self, name: str) -> Path:
        return self.checkpoints / f"{name}.h5ad"

    # ... tables, plots, interactive, objects, annotation analogously ...

PATHS = _Paths()
PARAMS = _cfg["thresholds"]   # plain dict is fine; freeze if you prefer
```

**Why.** Changing a threshold is one edit. Switching between container and host filesystems is one edit. Skills downstream (e.g., `consensus-nmf-multirun`) read `analysis_config.yaml::decisions::*` for replay-from-config; the convention guarantees that file exists and has a known shape.

**How to apply.** First script in a project creates `analysis_config.yaml` and `config.py`. Every subsequent script imports `from config import PATHS, PARAMS` and never hard-codes a path or threshold.

---

### 4. Output directory layout `03_results/`

```
03_results/
├── checkpoints/    # .h5ad / .rds / .pkl — analysis state at each stage
├── tables/         # .csv — DE results, master tables, signatures, program metadata
├── plots/          # .png / .pdf — static figures
├── interactive/    # .html — pathway-explorer dashboards, plotly
├── objects/        # .h5ad ready for CellxGene hosting (post schema-prep)
└── annotation/     # CXG annotation autosave per dataset (read-write from container)
```

Lower-case `03_results/` (not `03_Results/`) is the new convention. Existing projects with capital-R `03_Results/` may keep it; new projects use lower-case.

**Why.** Deterministic mental model. Anyone walking into the project sees what kind of output is where. `objects/` and `annotation/` are mounted into the CXG container by `scrna-cxg-host`; the layout encodes the deploy contract.

**How to apply.** `PATHS` (Convention 3) auto-creates these directories. A new analysis output goes into the matching subdirectory by extension and purpose; never at the project root.

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

## Practical adoption checklist

- [ ] Project root has `config.py` and `02_analysis/config/analysis_config.yaml`
- [ ] `03_results/{checkpoints,tables,plots,interactive,objects,annotation}/` exist or are auto-created on `PATHS.<key>` access
- [ ] First analysis script is named `00_<verb>_<noun>.py` (e.g., `00_build_anndata.py`)
- [ ] Long-running scripts (≥10 min, ≥2 substantive stages) write multiple checkpoints
- [ ] Recurring logic lives in `Python_scripts/<topic>_utils.py`, not at the top of an analysis script
- [ ] `analysis_config.yaml` has a `decisions:` section (empty is fine; skills append to it)

---

## DECISION PAUSE — Adoption scope

> **Stop here.** Present the options below to the user and wait for their selection before proceeding. Do not pick a default silently.

**Question for the user:** This project already has analysis files. How should the conventions be adopted?

| Option | What happens | When to choose | Default |
|--------|--------------|-----------------|---------|
| **A — Reference only** | Conventions are documented in `CLAUDE.md` / `AGENTS.md`; no project files are created or renamed | Established project; conventions apply to new files only | yes (safest) |
| **B — Scaffold new directories** | `03_results/{checkpoints,tables,plots,interactive,objects,annotation}/` created; `config.py` + `analysis_config.yaml` template added; existing scripts left alone | New project; want the scaffolding without touching pre-existing analysis files | |
| **C — Full retrofit** | Skill proposes renames for existing scripts to match `NN_verb_noun.py`, reorganises outputs into `03_results/<subdir>/`, and produces a single batch-rename script the user reviews before running | Project owner explicitly wants the retrofit and accepts the disruption | |

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

- [ ] **`PATHS` resolves and creates directories.** Running `from config import PATHS; print(PATHS.checkpoints)` prints an absolute path and the directory now exists on disk.
- [ ] **`PARAMS` loads from YAML.** `from config import PARAMS; PARAMS["de_fdr"]` returns the expected numeric without `KeyError`.
- [ ] **At least one numbered script reads from a checkpoint and writes another.** `grep -l "PATHS.checkpoint" 02_analysis/*.py` returns ≥1 file (or analogous for `Python_scripts/`).
- [ ] **`analysis_config.yaml` has a `decisions:` section.** Even if empty (`decisions: {}`), the key exists for downstream skills to append to.

---

## Common Pitfalls

### Pitfall: PATHS as plain attributes instead of properties

- **Symptom:** Running a fresh `git clone` of the project, the first script crashes with `FileNotFoundError: 03_results/checkpoints/`.
- **Cause:** `PATHS.checkpoints` was defined as a `Path` attribute computed once at import time, before the directory existed. The auto-create behaviour requires `@property`.
- **Fix:** Convert each path to a `@property` that calls `mkdir(parents=True, exist_ok=True)` before returning.

### Pitfall: Hard-coded paths inside numbered scripts

- **Symptom:** Switching from container to host filesystem requires editing five scripts; one gets missed and silently writes to the wrong place.
- **Cause:** A script wrote `pd.read_csv("/scratch/14616-DM/03_results/tables/x.csv")` directly instead of `pd.read_csv(PATHS.tables / "x.csv")`.
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
- **Cause:** The 13403-YD reference used capital-R; the new convention is lower-case (per project `CLAUDE.md`).
- **Fix:** New projects use lower-case from day one. Existing projects keep their capital-R; do not migrate mid-stream. `PATHS` reads the actual directory name from `analysis_config.yaml::paths::results_dir` so a project can opt into either.

---

## Complementary Skills

| When you need... | Use skill | Relationship |
|---|---|---|
| Build the first AnnData from CellRanger output | `cellranger-multi-to-anndata` | Next step (Stage 0) |
| Host the prepared `.h5ad` on CellxGene for the wet lab | `scrna-cxg-host` | Downstream consumer (uses `PATHS.objects`, `PATHS.annotation`) |
| Discover gene programs after annotation | `consensus-nmf-multirun` | Downstream consumer (uses `PATHS.checkpoints`, `PATHS.tables`) |
| Discipline for software-design changes (ADRs, reviews) | `architecture-first-dev` | Different scope (software design, not data-pipeline scaffolding) |
| Run MAD-based QC filtering on the produced AnnData | `single-cell-rna-qc` | Adjacent (consumes `00_raw.h5ad` from this convention's Stage 0) |

---

## Resources

This SKILL.md is the canonical reference. Other skills point here. There is no upstream documentation — these conventions are local to this skill library.
