# Phase 4: Phase-Based Artifact Structure

## Summary

This phase replaces the flat `03_results/{checkpoints,plots,tables}` layout with a
**phase-based artifact grouping** in which each analysis stage gets its own top-level
folder (`01_qc/`, `02_programs/`, …) that contains *both* the tables and the figures it
produced. Checkpoint objects (`.h5ad`, `.rds`, `.parquet`) are pulled out into a single
shared `objects/` directory because they are pipeline *state*, not *results* — they are
recomputable, large, and shared across stages, whereas tables and figures are the
human-facing deliverables that get read together when writing a figure legend or
assembling supplementary material.

The change is grounded in two real projects. `13403-YD_Christina` (scRNA-seq, scanpy)
already grew an organic phase-based layout under `figures/02_eda/`, `tables/03_scvi/`,
`validation/02_qc2_validation.json`, with objects isolated in `objects/`. `AdaW_eWAT_WL_2025`
(bulk RNA-seq, R) kept a flatter `plots/` + `tables/` split but added `checkpoints/`,
`interactive/` (HTML explorers), and accumulated `master_de.csv` / `master_gsea.csv`
cross-stage tables. This phase canonicalizes the best of both into a single scaffold that
SciAgent-toolkit ships, plus the config-path constants and script-naming conventions that
make scripts write to the right place by default.

## Current State

The project scaffold and its results layout currently live in **scbio-docker**, not in
SciAgent-toolkit. The relevant artifacts:

- `scbio-docker/templates/base/03_results/plots/.gitkeep`
- `scbio-docker/templates/base/03_results/tables/.gitkeep`
- `scbio-docker/templates/base/03_results/checkpoints/.gitkeep`

  i.e. the flat three-way split: `03_results/{plots,tables,checkpoints}/`.

- `scbio-docker/templates/base/02_analysis/config/config.R.template`
  (~290 lines) hardcodes the path constants:
  ```r
  DIR_RESULTS     <- file.path(PROJECT_ROOT, "03_results")
  DIR_CHECKPOINTS <- file.path(DIR_RESULTS, "checkpoints")
  DIR_TABLES      <- file.path(DIR_RESULTS, "tables")
  DIR_PLOTS       <- file.path(DIR_RESULTS, "plots")
  for (dir in c(DIR_CHECKPOINTS, DIR_TABLES, DIR_PLOTS)) { ... dir.create ... }
  ```
  There is no Python equivalent (`config.py`) in the template tree today, even though the
  observed scanpy project is Python-first.

- `scbio-docker/templates/base/02_analysis/config/pipeline.yaml.template`
  already contains a `schemas:` block describing `master_gsea_table`, `master_de_table`,
  and `master_tf_activities` with `required_columns` — i.e. the master-table concept exists
  in config but has no home directory in the scaffold.

- `scbio-docker/templates/base/.gitignore` ignores results contents but keeps `.gitkeep`:
  ```
  03_results/*/*
  !03_results/*/.gitkeep
  ```
  This glob is **exactly one level deep** (`*/*`), so it does *not* correctly cover a
  nested phase layout like `03_results/01_qc/figures/foo.png`.

- `scbio-docker/templates/base/docs/...` and the whole `templates/base/` tree are slated to
  move to SciAgent-toolkit under Phase 1/2 of this refactor. This phase assumes that move
  has happened (see Dependencies) and edits the scaffold in its new SciAgent-toolkit home.

- `SciAgent-toolkit/templates/AGENTS.md.template` currently documents an **idealized and
  partly wrong** structure: it shows `03_Results/{objects,tables,figures}/` (capital R,
  flat) and a critical rule "Normalize, then visualize. Compute all scores/embeddings and
  save to checkpoint files before running any visualization script." There is no
  formalized `0X_` / `0Xg_` compute/viz script-naming convention written down anywhere.

## Changes

### 4.1 Define the canonical `03_results/` phase-based layout

**Files:**
- new scaffold dirs under `SciAgent-toolkit/templates/<project-type>/03_results/`
- (replaces) `scbio-docker/templates/base/03_results/{plots,tables,checkpoints}/.gitkeep`

**What:** Adopt the following canonical structure as the single documented standard:

```
03_results/
├── objects/                 # STATE, not results: checkpoints (.h5ad, .rds, .parquet)
│   └── .gitkeep
├── master/                  # cross-stage accumulator tables (master_de.csv, …)
│   └── .gitkeep
├── interactive/             # standalone HTML dashboards / explorers
│   └── .gitkeep
├── 01_<slug>/               # first analysis stage
│   ├── tables/
│   │   └── .gitkeep
│   └── figures/
│       └── .gitkeep
├── 02_<slug>/               # second stage
│   ├── tables/
│   │   └── .gitkeep
│   └── figures/
│       └── .gitkeep
└── ...
```

Rules that define the standard:

1. **Stage folders** are `NN_<slug>/` where `NN` is a zero-padded two-digit numeric prefix
   (`01_`, `02_`, …) and `<slug>` is a laconic lowercase snake_case stage name
   (`qc`, `programs`, `pan_t`, `scvi`, `dge`, `gsea`).
2. Each stage folder contains exactly two child dirs: `tables/` and `figures/`. Both
   belong to the *same* phase so a legend author finds the figure and its source table
   side by side.
3. **`objects/`** holds checkpoints (recomputable state). It is *not* a stage folder and
   has no `tables/`/`figures/` children. Per-tool run directories that emit their own
   nested state (observed: `cnmf/` in the scanpy project) live under `objects/<tool>/`.
4. **`master/`** holds cross-stage accumulator tables (see 4.5).
5. **`interactive/`** holds self-contained HTML dashboards (see 4.6).
6. **`validation/`** is *optional* and added by project type, not baked into every scaffold
   (see ADR-4.3).

**Why:** Tables-with-figures-by-phase matches how results are actually consumed (paper
prep, supplementary tables, legend writing) and matches the layout `13403-YD_Christina`
arrived at organically. Isolating `objects/` keeps multi-GB recomputable state out of the
human-readable deliverable tree and lets `.gitignore` / data-sync rules treat it
differently.

**How:**
1. In the SciAgent-toolkit scaffold for the analysis project type, create the directory
   skeleton above with `.gitkeep` in every leaf dir.
2. Pre-populate stage folders `01_qc/` and `02_eda/` only (see ADR-4.1), each with empty
   `tables/.gitkeep` and `figures/.gitkeep`.
3. Always create `objects/.gitkeep`, `master/.gitkeep`, `interactive/.gitkeep`.
4. Delete the old `templates/base/03_results/{plots,tables,checkpoints}/` skeleton when the
   scaffold lands in SciAgent-toolkit.

### 4.2 Fix `.gitignore` for the nested layout

**Files:** scaffold `.gitignore` (moved from `scbio-docker/templates/base/.gitignore`)

**What:** Replace the one-level results glob with a recursive ignore that keeps directory
skeletons via `.gitkeep`, and never commits objects.

**Why:** The current `03_results/*/*` + `!03_results/*/.gitkeep` only matches one nesting
level. With `03_results/01_qc/figures/plot.png` the figure sits *two* levels deep and would
be silently committed. We also want `objects/` contents (large checkpoints) hard-ignored
regardless of depth.

**How:**
1. Replace the results section of `.gitignore` with:
   ```gitignore
   # Results: track directory skeleton via .gitkeep, never track contents
   03_results/**
   !03_results/
   !03_results/**/
   !03_results/**/.gitkeep

   # Checkpoints / state objects are always large + recomputable: never commit
   03_results/objects/**
   ```
2. Confirm with `git check-ignore -v 03_results/01_qc/figures/x.png` that figures are
   ignored while `03_results/01_qc/figures/.gitkeep` is tracked.

### 4.3 Rewrite config path constants (R and Python) to be phase-aware

**Files:**
- `config.R.template` (scaffold)
- new `config.py.template` (scaffold) — Python parity for scanpy-first projects

**What:** Replace the flat `DIR_CHECKPOINTS / DIR_TABLES / DIR_PLOTS` triple with:
- a single `DIR_OBJECTS` (was `DIR_CHECKPOINTS`),
- `DIR_MASTER`, `DIR_INTERACTIVE`,
- and helper functions that *resolve a stage's* `tables/` and `figures/` dir on demand and
  create it lazily, rather than pre-creating a fixed flat set.

**Why:** Stage folders are open-ended and named per project; the config cannot enumerate
them ahead of time. A resolver (`stage_dir("01_qc", "figures")`) keeps "config not
hardcoding" (AGENTS.md rule 1) intact while supporting the nested layout. Python parity is
required because the canonical scanpy project has no R config today.

**How (R):**
1. Replace the path block with:
   ```r
   DIR_RESULTS     <- file.path(PROJECT_ROOT, "03_results")
   DIR_OBJECTS     <- file.path(DIR_RESULTS, "objects")       # state, not results
   DIR_MASTER      <- file.path(DIR_RESULTS, "master")
   DIR_INTERACTIVE <- file.path(DIR_RESULTS, "interactive")

   #' Resolve (and lazily create) a per-stage output dir.
   #' kind = "tables" | "figures"
   stage_dir <- function(stage, kind = c("tables", "figures")) {
     kind <- match.arg(kind)
     d <- file.path(DIR_RESULTS, stage, kind)
     if (!dir.exists(d)) dir.create(d, recursive = TRUE)
     d
   }
   ```
2. Keep `DIR_CHECKPOINTS <- DIR_OBJECTS` as a one-line deprecated alias *inside the template
   comment block only* if desired — but since breaking changes are fine (single user), prefer
   to drop it and rename call sites. Update `load_or_compute()` to default its checkpoint
   root to `DIR_OBJECTS`.
3. Only auto-create `DIR_OBJECTS`, `DIR_MASTER`, `DIR_INTERACTIVE` at config load. Stage dirs
   are created on first `stage_dir()` call.

**How (Python):** Ship `config.py.template` mirroring the same constants and a
`stage_dir(stage, kind="figures") -> pathlib.Path` helper that `mkdir(parents=True,
exist_ok=True)` on demand, plus `DIR_OBJECTS`, `DIR_MASTER`, `DIR_INTERACTIVE`.

### 4.4 Define the stage slug and how stages are created

**Files:** scaffold scripts; `pipeline.yaml.template` (new `stages:` list); `AGENTS.md.template`

**What:** Formalize that a stage slug = `NN_<laconic_snake_case>`. Stages are **generic in
the shipped scaffold** (only `01_qc`, `02_eda` pre-populated) and are **added by the author
or agent on demand** as the analysis grows, with an optional declared list in
`pipeline.yaml`.

**Why:** Slugs are inherently project-specific (`programs`, `pan_t`, `scvi` came from one
real project's biology). Hardcoding a long fixed list at init time would be wrong for most
projects and would create empty noise. But a *declared* list in `pipeline.yaml` gives the
agent a canonical, validatable enumeration so it doesn't invent inconsistent names
(`02_eda` vs `02_EDA` vs `2_eda`).

**How:**
1. Add to `pipeline.yaml.template`:
   ```yaml
   # Analysis stages. Add entries as the analysis grows.
   # Each stage gets 03_results/<id>/{tables,figures}/ created on first write.
   stages:
     - id: "01_qc"
       title: "Quality control"
     - id: "02_eda"
       title: "Exploratory analysis"
   ```
2. In `AGENTS.md.template`, document the slug rule: two-digit zero-padded prefix, lowercase
   snake_case, ≤ ~2 words, no spaces; new stages must be appended to `pipeline.yaml:stages`
   before scripts write to them.
3. `config.{R,py}` may optionally read `stages` to validate that a `stage_dir()` argument is
   a declared stage (warn, do not hard-fail).

### 4.5 Give master / cross-stage tables a canonical home

**Files:** scaffold `03_results/master/`; `pipeline.yaml.template` `schemas:` block; `config.{R,py}.template`

**What:** Cross-stage accumulator tables (observed: `master_de.csv`, `master_gsea.csv` in
`AdaW_eWAT_WL_2025`) live in `03_results/master/`, *not* inside any numbered stage, and are
validated against the existing `schemas:` block in `pipeline.yaml`.

**Why:** A master table is by definition not owned by one phase — it concatenates results
across stages/contrasts. Putting it under `02_dge/tables/` would lie about its scope.
`pipeline.yaml` already declares `master_de_table` / `master_gsea_table` schemas, so the
directory and the schema now line up.

**How:**
1. Create `03_results/master/.gitkeep`.
2. Expose `DIR_MASTER` in both configs (done in 4.3).
3. Document in `AGENTS.md.template` that any table whose rows span more than one stage goes
   to `master/` and must conform to a `schemas:` entry in `pipeline.yaml`; per-stage tables
   stay in their stage `tables/` dir.

### 4.6 Place interactive HTML dashboards

**Files:** scaffold `03_results/interactive/`; `AGENTS.md.template`

**What:** Self-contained HTML explorers/dashboards (observed: 15+ files in
`AdaW_eWAT_WL_2025/03_results/interactive/`) live in a single top-level
`03_results/interactive/` directory.

**Why:** HTML dashboards are usually cross-cutting (one explorer spans all DE results), are
large self-contained bundles, and are opened directly in a browser rather than dropped into
a paper. They don't belong inside a numbered stage's `figures/`, which is reserved for
static publication figures. Keeping them separate also lets `.gitignore` / sync rules treat
heavy HTML bundles distinctly.

**How:**
1. Create `03_results/interactive/.gitkeep`.
2. Expose `DIR_INTERACTIVE` (done in 4.3).
3. If a dashboard is genuinely single-stage, the author *may* place it in that stage's
   `figures/`; default and documented home is `interactive/`.

### 4.7 Formalize the `0X_` / `0Xg_` compute/viz script convention

**Files:** `AGENTS.md.template`; scaffold `02_analysis/scripts/.gitkeep` (add a README stub)

**What:** Document the script-naming split observed in `13403-YD_Christina`:
- `0X_<slug>.py` / `0X_<slug>.R` — **compute** scripts (produce checkpoints in `objects/`,
  tables in the stage `tables/`, and a validation gate).
- `0Xg_<slug>.py` / `0Xg_<slug>.R` — **viz** scripts (read checkpoints + tables, write to the
  stage `figures/`). The `g` suffix = "graphics".

This operationalizes the existing AGENTS.md rule "Normalize, then visualize": compute
scripts must finish (and pass validation) before the matching `…g` viz script runs.

**Why:** The naming makes the compute→viz dependency visible in `ls` ordering
(`03_scvi.py` sorts immediately before `03g_scvi.py`), enforces the normalize-then-visualize
discipline mechanically, and is already battle-tested in a real project. Writing it into
`AGENTS.md.template` means every agent in every new project inherits the convention instead
of rediscovering it.

**How:**
1. Add a "Script naming" subsection to `AGENTS.md.template`:
   - `NN_<slug>` compute, `NNg_<slug>` viz, both keyed to the same stage `NN_<slug>` in
     `03_results/` and the same `pipeline.yaml:stages` id.
   - Compute writes checkpoints to `objects/`, per-stage tables to `NN_<slug>/tables/`, and a
     validation artifact (see ADR-4.3). Viz writes only to `NN_<slug>/figures/`.
2. Replace the stale `03_Results/{objects,tables,figures}` block in `AGENTS.md.template`
   (capital R, flat) with the canonical 4.1 tree (lowercase `03_results`, phase-based).

### 4.9 Define the `03_results/{phase}/README.md` caption convention

**Files:** `templates/03_results/01_qc/README.md` (scaffold stub); `templates/03_results/02_eda/README.md` (scaffold stub); `templates/AGENTS.md.template` (new "Artifact captions" section); `agents/analysis-base/doc-curator.md` (completeness check); `agents/analysis-base/handoff.md` (pre-write verification, owned by Phase 3 / 3.4 but caption logic added here).

**What:** Every `03_results/{phase}/` directory ships with a `README.md` that serves as a **figure legend sheet**: a per-artifact inventory where each entry states the scientific finding, traces the exact script and function that produced the artifact, and records the active config parameter values. This is the layer that makes the repo self-documenting and reproducible without opening any other file.

**The canonical caption format** (baked into every scaffolded README stub):

```markdown
# NN_<slug> — <Phase Title>

One sentence: what decision does this phase make, what does it produce, what comes next.

---

## <filename.ext>

**<Finding: one sentence, scientific. What does this artifact show — not describe.>**

[Optional: 1–3 sentences of context if the finding alone is insufficient.]

| | |
|---|---|
| Script | `02_analysis/<phase>/<NN_script>.R` or `.py` |
| Function | `exact_function_name()` |
| Config | `path.to.key = value` |
| Input | `<primary input path>` |

---

## <next-filename.ext>
...
```

**Invariants agents must follow:**

1. **Finding is a scientific statement, not a label.** Write "Cells predicted as doublets form a distinct peripheral cluster isolated from the main manifold" not "UMAP colored by doublet score." The finding sentence should belong in a paper's figure legend.
2. **Function is exact** — must name the actual function/method/class that produced the artifact. Agents must read the script to fill this in. Invented function names are worse than a blank.
3. **Config key + value** — the value active when the artifact was produced, not a pointer to the config file. A reader must be able to reproduce the output from the README alone.
4. **No `docs/_internal/` references** — these READMEs are public-facing. The one-way rule applies (Phase 3, section 3.2). State outcomes, never cite internal reasoning docs.
5. **Artifact entries are written as artifacts are produced**, not batched at the end of a session. The handoff agent flags any `03_results/{phase}/` file not yet captioned.

**Why:** Without this convention, arriving cold at a results directory means opening every figure to understand what it shows, grepping scripts to find what produced it, and re-reading config to find what parameters were used. This friction compounds across sessions and is exactly the kind of observability loss the user cited ("I need to have 100% traceability"). The caption README eliminates that friction and makes the repo itself the communication layer between scientist and agent — no ad-hoc explanation needed before "go write captions for phase 01_qc."

**How:**
1. Write `templates/03_results/01_qc/README.md` and `templates/03_results/02_eda/README.md` as stubs containing the format spec and one empty artifact block to show the pattern.
2. Add an "Artifact captions" section to `AGENTS.md.template`:
   ```markdown
   ## Artifact captions

   Every file in `03_results/{phase}/` must have a caption entry in that phase's `README.md`
   before the session closes. Use the format in the README stub.

   Caption entry rules:
   - Finding: one scientific sentence (what it shows, not what it is)
   - Function: exact name from the script — read the script, do not guess
   - Config: active key = value at time of production
   - No references to `docs/_internal/`
   ```
3. Add to `doc-curator.md`: list every file in `03_results/**/*.{pdf,png,svg,csv,tsv,html}` that has no matching `## <filename>` heading in the sibling `README.md`. Report as "uncaptioned artifacts" with their paths.
4. The scaffold's `README.md` stubs are committed (they define the contract); the generated caption content is also committed (it is the human-facing deliverable documentation).

### 4.8 Correct the AGENTS.md.template directory diagram

**Files:** `SciAgent-toolkit/templates/AGENTS.md.template`

**What:** The "Directory structure" block currently shows `00_Data/`, `01_Scripts/`,
`02_Analysis/`, `03_Results/{objects,tables,figures}`. Replace with the real scaffold:
`00_data/`, `01_modules/`, `02_analysis/{config,scripts,notebooks,helpers}/`, and the Phase-4
`03_results/` tree from 4.1. Fix the casing (`03_Results` → `03_results`) so the documented
rule paths match the actual filesystem the config constants point at.

**Why:** The template currently lies about the layout; an agent reading it will write to
`03_Results/tables/` (wrong case, wrong shape) and miss every file. This must match 4.1
exactly.

**How:** Edit the fenced block; cross-check every path against the scaffold produced by 4.1.

## Open ADRs

### ADR-4.1: Pre-populate stage folders, or create entirely on demand?
**Options:**
- A — Ship an empty `03_results/` with only `objects/`, `master/`, `interactive/`; create all
  `NN_<slug>/` lazily via `stage_dir()`.
- B — Pre-populate a small generic seed set (`01_qc/`, `02_eda/`) plus the three shared dirs.
- C — Pre-populate a long fixed list (`01_qc`…`06_export`).
**Recommended:** B
**Why:** A leaves a near-empty `03_results/` that gives the author no example of the
tables+figures pattern and no `.gitkeep` to make the dir survive clone. C imposes biology
the project may not do (the scanpy project's `pan_t` is domain-specific) and creates empty
noise. B shows the pattern with the two near-universal first stages (QC always happens; EDA
almost always) while leaving the rest to grow on demand via `stage_dir()` + `pipeline.yaml`.
**Blocking implementation:** no (B is the assumed default in 4.1).

### ADR-4.2: Is the stage slug fixed at init time or always generic?
**Options:**
- A — Generic seed slugs only (`01_qc`, `02_eda`); author/agent adds the rest.
- B — Prompt for slugs at `sciagent new project` and bake them in.
- C — Project-type-specific seed lists (analysis vs software-tool).
**Recommended:** A, refined by C
**Why:** Slugs encode biology not known at init time, so prompting (B) front-loads decisions
the author can't yet make and produces stale names. A keeps init dumb and lets stages
accrete. The one refinement worth taking from C: the *software-tool* project type (the
`pathway-explorer` child-toolkit pattern, Phase-3 multi-type work) should not ship a
`03_results/` at all — it uses `src/`+`tests/`. So "generic seed" is correct *for the
analysis type*; the software-tool type omits this scaffold entirely.
**Blocking implementation:** no — but coordinate with the Phase-3 multi-project-type work so
the software-tool scaffold doesn't inherit `03_results/`.

### ADR-4.3: Bake a `validation/` dir + gate convention into the scaffold?
**Options:**
- A — Ship `03_results/validation/` + document the "every compute script writes
  `validation/<stage>.json`, exit on FAIL" gate in `AGENTS.md.template`.
- B — Document the gate convention but do not ship the dir; create it lazily like stage dirs.
- C — Leave validation gates entirely out of the template (too opinionated).
**Recommended:** B
**Why:** The validation-gate pattern (`validation/02_qc2_validation.json`, compute script
exits on FAIL) is one of the strongest, most transferable practices observed in
`13403-YD_Christina` — worth documenting and giving a canonical path (`03_results/validation/`
or per-stage `NN_<slug>/validation.json`; recommend a single top-level `validation/` keyed by
stage id, matching the observed project). But forcing an empty `validation/.gitkeep` on
every project (including ones with no compute gates yet) is mildly opinionated noise, so
create it lazily on first gate write and expose `DIR_VALIDATION` in config. This documents
the discipline without imposing empty dirs. C is rejected: the gate pattern is exactly the
kind of methodology SciAgent-toolkit should propagate.
**Blocking implementation:** no.

### ADR-4.4: Where do master / cross-stage tables live?
**Options:**
- A — `03_results/master/`.
- B — `03_results/00_master/` (sorts first as a pseudo-stage).
- C — Inside the latest stage that produces them.
**Recommended:** A
**Why:** C lies about scope (the table spans stages). B abuses the numeric-stage convention
for something that isn't a pipeline stage and would imply it has `tables/`+`figures/`
children, which it doesn't. A is an unambiguous, non-numeric sibling that pairs cleanly with
the existing `pipeline.yaml:schemas` master-table definitions.
**Blocking implementation:** no.

## Dependencies

- **Depends on:** Phase 1/2 (move `scbio-docker/templates/base/` → SciAgent-toolkit and strip
  `init-project.sh` to container-only). This phase edits the scaffold in its new
  SciAgent-toolkit home; doing it before the move would put the changes in the wrong repo.
- **Depends on (soft):** Phase 3 multi-project-type support — ADR-4.2/C requires that the
  software-tool project type exists so it can opt out of `03_results/`. If Phase 3 lands
  after this, the analysis-type scaffold here is still correct; only the software-type
  carve-out is deferred.
- **Enables:** consistent artifact paths that the AI harness (roles/skills/agents) can rely
  on — any "write the QC figure" command resolves deterministically to
  `03_results/01_qc/figures/`. Enables the validation-gate methodology to be referenced from
  guidelines.

## Breaking Changes

- `03_results/{plots,tables,checkpoints}/` flat layout is **removed**; replaced by
  `03_results/{objects,master,interactive}/` + per-stage `NN_<slug>/{tables,figures}/`.
- Config constants renamed: `DIR_CHECKPOINTS` → `DIR_OBJECTS`; `DIR_TABLES` / `DIR_PLOTS`
  removed in favor of `stage_dir(stage, kind)`. New: `DIR_MASTER`, `DIR_INTERACTIVE`,
  `DIR_VALIDATION`. Existing scripts referencing the old constants must be updated.
- `.gitignore` results rules changed from one-level (`03_results/*/*`) to recursive.
- `AGENTS.md.template` directory diagram and rules rewritten (casing + shape change:
  `03_Results` → `03_results`, flat → phase-based).
- New required `config.py.template` (Python parity) added to the analysis scaffold.

Acceptable per project ground rules: single user (Anton), no downstream consumers; clean
modularity preferred over backwards compatibility. No migration shims shipped.

## Estimated Scope

| Item | Files | Approx. line delta |
|---|---|---|
| Scaffold dir skeleton + `.gitkeep`s (4.1) | ~10 new `.gitkeep`, -3 old | +10 / -3 files |
| `.gitignore` results rules (4.2) | 1 | ~ -3 / +7 lines |
| `config.R.template` rewrite (4.3) | 1 | ~ -25 / +35 lines |
| `config.py.template` new (4.3) | 1 new | ~ +120 lines |
| `pipeline.yaml.template` `stages:` (4.4) | 1 | +8 lines |
| Master/interactive/validation dirs + config (4.5/4.6, ADR-4.3) | within above | +6 lines |
| `AGENTS.md.template` script naming + diagram + caption section (4.7/4.8/4.9) | 1 | ~ -15 / +65 lines |
| Scaffold `README.md` stubs for `01_qc/`, `02_eda/` (4.9) | 2 new | ~ +40 lines each |
| `doc-curator.md` uncaptioned-artifact check (4.9) | 1 | +15 lines |

**Total:** ~6 files edited, ~15 created (`.gitkeep`s + README stubs), net ≈ +350 / −45 lines. No
runtime code paths in scbio-docker touched (scaffold + templates only).
