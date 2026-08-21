# Phase 2: `sciagent new project` — Full Scaffold Expansion

## Summary

Today `sciagent new project` is a three-file stub: it copies `AGENTS.md`, `CLAUDE.md`,
and `context.md` into a target directory and tells the user to run `sciagent activate`.
Meanwhile, scbio-docker's `init-project.sh` (602 lines) is the real scaffolder — it builds
the entire `00_data/ … 03_results/` tree, renders R/Python config, generates the
devcontainer + docker-compose, optionally inits git, and registers submodules. Phase 1 of
this refactor establishes the clean seam: scbio-docker owns the *container substrate* only,
and SciAgent-toolkit owns *all project-level context*. This phase is where that decision
becomes concrete code.

Phase 2 expands `_new_project()` in `lib/sciagent/new.sh` so that `sciagent new project`
produces a complete, opinionated project scaffold: the numbered analysis tree, the new
`docs/_internal/` reasoning namespace, the phase-based `03_results/` artifact layout, and
decomposed context docs (replacing the `context.md` monolith). It introduces a
`--type {analysis,software}` flag so the scaffold can branch between the current
analysis-repo pattern and a `src/`/`tests/` software-tool pattern. It does **not** generate
docker-compose or devcontainer files — that is the container substrate's job and is invoked
as a separate, optional step (Phase 3). The templates that encode project structure migrate
out of `scbio-docker/templates/base/` into SciAgent-toolkit.

## Current State

**`lib/sciagent/new.sh` — `_new_project()` (the stub):**

```
/data1/users/antonz/pipeline/scbio-docker/toolkits/SciAgent-toolkit/lib/sciagent/new.sh
```

`_new_project()` (≈35 lines of the 180-line file) does exactly this:

1. `mkdir -p "$dir"`; derive `pid` from `basename`.
2. Loop over three templates in `$SCIAGENT_TOOLKIT/templates/`:
   `AGENTS.md.template`, `CLAUDE.md.template`, `context.md.template`.
3. For each, skip if the output already exists, else call `_subst` and write.
4. Print "Next: cd $dir && sciagent activate <role>".

`_subst()` is the only variable-substitution machinery:

```bash
_subst() {
    local src="$1" dst="$2" pid="$3"
    local date_str; date_str=$(date +%Y-%m-%d)
    sed -e "s|{{PROJECT_ID}}|$pid|g" \
        -e "s|{{DATE}}|$date_str|g" \
        -e "s|{{SKILL_NAME}}|$pid|g" \
        "$src" > "$dst"
}
```

So today's substitution vocabulary is `{{PROJECT_ID}}`, `{{DATE}}`, `{{SKILL_NAME}}`.

**Templates that exist in SciAgent-toolkit today:**

```
toolkits/SciAgent-toolkit/templates/
├── AGENTS.md.template            (1766 B — describes a 00_Data/…03_Results/ tree inline)
├── CLAUDE.md.template            (11 B — just "@AGENTS.md")
├── context.md.template           (735 B — the monolith to decompose)
└── analysis_config.yaml.template (5726 B — RNA-seq GSEA/colors/schemas; recently added)
```

Note `AGENTS.md.template` hardcodes a *different, capitalized* tree
(`00_Data/`, `01_Scripts/`, `02_Analysis/`, `03_Results/{objects,tables,figures}`) than
both `init-project.sh` (`00_data/ … 03_results/{checkpoints,plots,tables}`) and
`analysis_config.yaml.template` (`03_results/{checkpoints,tables,plots,interactive}`).
**There are three conflicting directory conventions in the repo right now.** Phase 2 must
pick one.

**`scbio-docker/scripts/init-project.sh` (602 lines) — the structure it currently builds:**

```
00_data/{raw,processed,references}
01_modules/.ref
02_analysis/{config,helpers,scripts,notebooks}
03_results/{checkpoints,plots,tables}     # flat
docs/{raw, ai-generated/{vignettes,research}, plan/phase-1}
logs/
.vscode/settings.json
.devcontainer/{devcontainer.json, docker-compose.yml, .env, scripts/}
README.md, notes.md
02_analysis/config/{config.R, pipeline.yaml, color_config.R}
```

…plus `.gitkeep` seeding, interactive prompts (species/mounts/resources), `--git-init`,
and `--with-submodules` (registers `RNAseq-toolkit` and `SciAgent-toolkit` at
`01_modules/`). Templates live under `scbio-docker/templates/base/`.

**Per Phase 1**, everything in `init-project.sh` that touches project structure, docs, or
config moves here; only docker-compose/devcontainer generation (≈150 lines) stays in
scbio-docker.

## Changes

### 2.1 Migrate project-structure templates into SciAgent-toolkit

**Files:**
- New: `toolkits/SciAgent-toolkit/templates/project/` (template root for project scaffolds)
- Source (to be deleted in Phase 1/this phase): `scbio-docker/templates/base/` (everything
  except `.devcontainer/` and `.vscode/`, which scbio-docker keeps)

**What:** Establish `templates/project/` as the canonical home for everything that defines a
project's directory structure, docs namespace, and language-agnostic config. Move the
following from `scbio-docker/templates/base/` into it:

- `docs/README.md.template` → `templates/project/docs/README.md.template` (rewritten, see 2.3)
- `docs/notes.md.template` → `templates/project/docs/notes.md.template`
- `docs/plan/README.md.template` → `templates/project/docs/plan/README.md.template`
- `README.md` (project root readme) → `templates/project/README.md.template`
- `.gitignore` → `templates/project/.gitignore.template`

**Why:** Phase 1's seam says scbio-docker "must know NOTHING about project directory
structure." These templates *are* the directory structure. They cannot stay in the
container repo.

**How:**
1. `git mv` the listed files from `scbio-docker/templates/base/` into
   `toolkits/SciAgent-toolkit/templates/project/` (the toolkit is a submodule, so these are
   two separate commits — one `git rm` in scbio-docker, one `git add` in the toolkit).
2. Keep `.devcontainer/` and `.vscode/` template trees in scbio-docker — those are
   substrate (consumed by Phase 3 docker-compose generation), not project structure.
3. Delete the now-orphaned directory placeholders in `scbio-docker/templates/base/`.

### 2.2 Decompose `context.md` into `docs/_internal/`

**Files:**
- Remove from scaffold output: `context.md` (top-level monolith)
- New: `templates/project/docs/_internal/scientific-context.md.template`
- New: `templates/project/docs/_internal/sessions/.gitkeep`
- New: `templates/project/docs/_internal/reasoning/.gitkeep`
- New: `templates/project/docs/_internal/scratch/.gitkeep`
- Edit: `templates/AGENTS.md.template` (point "Project context" at the new path)

**What:** Replace the single `context.md` file with the `docs/_internal/` namespace agreed
in the architecture:

```
docs/_internal/
├── scientific-context.md   # was context.md (question, datasets, hypotheses, goals)
├── reasoning/              # decision traces, why-not logs, agent thinking
├── sessions/              # dated handoffs: YYYY-MM-DD_*.md
└── scratch/               # exploratory writeups before they crystallize
docs/                      # public-facing: guides, plan, notes, README
```

**Why:** The `context.md` monolith conflates the stable scientific framing (rarely changes)
with running session notes and decision logs (change constantly). Splitting them lets the
`handoff` agent write dated session files without churning the canonical context, and gives
the `doc-curator`/`reasoning` agents a stable namespace. `_internal/` signals
"agent-facing, not the public deliverable" — distinct from `docs/` which is the published
plan/guides.

**How:**
1. Rename `templates/context.md.template` →
   `templates/project/docs/_internal/scientific-context.md.template`. Keep its existing
   sections (Scientific Question / Datasets / Hypotheses / Analysis Goals / References).
2. Create `.gitkeep` files under `reasoning/`, `sessions/`, `scratch/`.
3. Update `AGENTS.md.template`'s "Project context" block to reference
   `docs/_internal/scientific-context.md` as the canonical project framing, and to point
   handoffs at `docs/_internal/sessions/`.
4. Update the `handoff` agent's output convention (Phase 5 territory; note the dependency)
   to write `docs/_internal/sessions/YYYY-MM-DD_<slug>.md`.

### 2.3 Phase-based `03_results/` artifact layout

**Files:**
- New: `templates/project/docs/README.md.template` (documents the new layout)
- Edit: `_new_project()` directory-creation list
- Edit (later phases): `analysis_config.yaml.template` / R config `paths:` block (see ADR-2.1)

**What:** Replace the flat `03_results/{checkpoints,plots,tables}` with phase-based artifact
grouping. Each phase folder holds **both** the tables and figures it produces; checkpoint
state objects are kept separate because they are state, not deliverables:

```
03_results/
├── objects/                # checkpoint state (.h5ad, .rds) — NOT a deliverable
│   └── .gitkeep
├── 01_qc/                  # tables + figures for the QC phase
│   └── .gitkeep
├── 02_programs/            # tables + figures for program discovery, etc.
│   └── .gitkeep
└── _scratch/              # throwaway plots/tables
    └── .gitkeep
```

The scaffold seeds `objects/` and a single starter phase folder (`01_qc/`); analysts add
`02_*`, `03_*` as the plan progresses. `docs/plan/` phase numbering and `03_results/` phase
numbering are intended to track each other.

**Why:** The flat layout forces a plot and its source table into two unrelated trees
(`plots/` vs `tables/`), so reviewing "everything phase 2 produced" means cross-referencing
filename prefixes across directories. Phase-folder grouping co-locates the artifact set a
reader actually wants together. Keeping `objects/` separate prevents large reusable state
files from being mistaken for results and lets `.gitignore` treat them differently.

**How:**
1. In `_new_project()`, generate `03_results/objects/`, `03_results/01_qc/`,
   `03_results/_scratch/` with `.gitkeep`.
2. Document the convention in `docs/README.md.template`: "each `NN_phase/` folder under
   `03_results/` contains the tables and figures for that phase; reusable checkpoint objects
   live in `03_results/objects/`."
3. Update the R config / yaml `paths:` block in a follow-up (ADR-2.1 decides where that
   config lives, which gates whether Phase 2 edits it or a toolkit does).

### 2.4 Expand `_new_project()` to build the full tree

**Files:** `lib/sciagent/new.sh`

**What:** Rewrite `_new_project()` from a 3-file copy into a full scaffolder driven by a
template manifest. Pseudocode of the new flow:

```bash
_new_project() {
    # parse: <dir> [--type analysis|software-tool] [--species ...] [--genome ...]
    #               [--title ...] [--git] [--with-submodules] [--force]
    # 1. resolve dir, derive PROJECT_ID
    # 2. select template root:  templates/project/<type>/  (+ shared templates/project/_common/)
    # 3. create directory tree for that type (see 2.3 / 2.7)
    # 4. render each *.template through _subst (expanded vars, see 2.5)
    # 5. seed .gitkeep into empty dirs
    # 6. optionally: git init + submodule registration (2.6)
    # 7. print next-steps incl. the scbio-docker compose hook (2.8)
}
```

The directory list for `--type analysis` (default):

```
00_data/{raw,processed,references}
01_modules/.ref
02_analysis/{config,helpers,scripts,notebooks}
03_results/{objects,01_qc,_scratch}
docs/{plan,_internal/{reasoning,sessions,scratch}}
logs/
```

**Why:** This is the core of the phase — making `sciagent new project` do what
`init-project.sh` did for structure, plus the new namespaces. Driving it from a manifest
(rather than the long imperative `for dir in …; do mkdir` block in init-project.sh) keeps
the two project types (2.7) from duplicating logic.

**How:**
1. Replace the three-template `for` loop with two passes: (a) `mkdir -p` over a
   type-specific dir list; (b) walk the template root with `find … -name '*.template'`,
   compute each output path by stripping `.template` and re-rooting under `$dir`, render via
   `_subst`.
2. Honor a `--force` flag for the existing "skip (exists)" behavior; default stays
   non-destructive (skip + warn).
3. Keep the function under ~120 lines; push git/submodule logic into a helper `_new_git()`
   (2.6) so the structure pass stays readable.

### 2.5 Expand the template-variable vocabulary

**Files:** `lib/sciagent/new.sh` (`_subst`)

**What:** Grow `_subst` from 3 tokens to the set the migrated templates actually need.
init-project.sh substitutes: `{{PROJECT_NAME}}`, `{{DATE}}`, `{{TEMPLATE_TYPE}}`,
`{{SPECIES}}`, `{{SPECIES_DB}}`, `{{GENOME_BUILD}}`, `{{PROJECT_TITLE}}`, plus
compose-only tokens (`{{IMAGE_VERSION}}`, `{{MAX_CPUS}}`, etc.) that stay in scbio-docker.

Proposed `sciagent`-owned vocabulary (see ADR-2.2):

| Token | Source | Default |
|-------|--------|---------|
| `{{PROJECT_ID}}` | `basename $dir` | — |
| `{{PROJECT_TITLE}}` | `--title` | `<PROJECT_ID> Analysis` |
| `{{DATE}}` | `date +%Y-%m-%d` | — |
| `{{PROJECT_TYPE}}` | `--type` | `analysis` |
| `{{SPECIES}}` | `--species` | `Mus musculus` |
| `{{SPECIES_DB}}` | derived from species | `MM` |
| `{{GENOME_BUILD}}` | `--genome` | `mm10` |

Note: keep `{{PROJECT_NAME}}` as a *deprecated alias* mapped to the same value as
`{{PROJECT_ID}}` during migration so the migrated init-project.sh templates render without
a mass find/replace; remove the alias once templates are normalized to `{{PROJECT_ID}}`.

**Why:** The migrated `config.R.template`, `pipeline.yaml.template`,
`analysis_config.yaml.template`, `notes.md.template`, and `README.md.template` all reference
species/genome/title tokens. The current 3-token `_subst` would leave those literally
`{{SPECIES}}` in output. We also drop `{{SKILL_NAME}}` from the *project* substitution path
(it belongs only to `_new_skill`).

**How:**
1. Refactor `_subst` to take an associative array (or build the `sed -e` list from a
   `declare -A vars`) instead of positional args.
2. Add a `_derive_species_db()` helper mirroring init-project.sh's `case "$SPECIES"`
   (mouse→MM/mm10, human→HS/hg38, else prompt/blank).
3. Validate: any unresolved `{{…}}` token in rendered output prints a warning (cheap
   `grep -l '{{' "$dir"` post-pass) so missing-variable bugs surface immediately.

### 2.6 Git init + submodule registration (optional, behind flags)

**Files:** `lib/sciagent/new.sh` (new `_new_git()` helper)

**What:** Port init-project.sh's `--git-init` and `--with-submodules` behavior, gated behind
`--git` and `--with-submodules` flags (default off). `--with-submodules` registers, at
`01_modules/`, the toolkits relevant to the project type:

- analysis: `RNAseq-toolkit` (branch `dev`) + `SciAgent-toolkit` (branch `main`)
- software-tool: `SciAgent-toolkit` only (no analysis toolkit)

Preserve the SSH-then-`gh`-HTTPS fallback from init-project.sh.

**Why:** This is project-level concern (it provisions the project's module set), so per the
Phase 1 seam it belongs in SciAgent-toolkit, not scbio-docker. It is kept *optional* and
*off by default* because `sciagent new project` is often run *inside* an already-cloned
repo where git already exists (see ADR-2.3).

**How:**
1. Add `_new_git()` that runs only when `--git`/`--with-submodules` is passed.
2. `git init` only if `.git` absent; initial commit message references project type + date.
3. Lift the `add_submodule_with_fallback()` function verbatim from init-project.sh; make the
   toolkit list a per-type array.
4. If neither flag is passed, print a hint: "run with `--with-submodules` to attach
   RNAseq-toolkit / SciAgent-toolkit under 01_modules/."

### 2.7 `--type {analysis,software-tool}` branching

**Files:**
- `lib/sciagent/new.sh`
- New: `templates/project/analysis/` and `templates/project/software-tool/`
- New: `templates/project/_common/` (shared: AGENTS/CLAUDE/docs/_internal)

**What:** Two first-class project types:

- **`analysis`** (default): the numbered tree above. The DC_hum_verse pattern — an analysis
  repo — is this type.
- **`software-tool`**: a Python/R package layout:

  ```
  src/<pkg>/            # or R/  for an R package
  tests/
  docs/{_internal/{reasoning,sessions,scratch}}
  pyproject.toml.template  (or DESCRIPTION.template)
  README.md.template
  ```

  The pathway-explorer / child-toolkit pattern (an analysis repo with a software toolkit
  submodule) is expressed as: parent = `analysis`, child = `software-tool` registered under
  `01_modules/`.

`docs/_internal/`, `AGENTS.md`, `CLAUDE.md`, and the substitution machinery are **shared**
across both types (live in `templates/project/_common/`); only the code/data trees differ.

**Why:** The architecture explicitly requires multi-project-type support and names both
patterns as first-class. Sharing `_common/` keeps the AI-context and reasoning namespace
identical regardless of type, so the harness behaves the same whether you're in an analysis
repo or a tool repo.

**How:**
1. Template selection: render `_common/` first, then overlay `<type>/`.
2. `--type` validates against `{analysis,software-tool}`; unknown → error listing valid
   types.
3. Stamp the chosen type into `{{PROJECT_TYPE}}` and into `AGENTS.md` so agents know which
   conventions apply (analysis house-style vs `src/`/`tests/` engineering style).
4. Software-tool details (pyproject vs DESCRIPTION, test runner) are deferred to Phase 6;
   Phase 2 ships the directory skeleton + `_common/` and a TODO-stub `pyproject.toml.template`.

### 2.8 Container-substrate hand-off (compose generation stays in scbio-docker)

**Files:** `lib/sciagent/new.sh` (next-steps output only)

**What:** `sciagent new project` does **not** write `.devcontainer/` or
`docker-compose.yml`. After scaffolding it prints the exact scbio-docker command to layer
the container substrate on top:

```
Next steps:
  1. cd <dir>
  2. Add container substrate:
       <scbio-docker>/scripts/init-container.sh <dir> \
         --type <type> --data-mount raw:/path/to/data:ro
  3. Open in VS Code → Reopen in Container
  4. sciagent activate base
```

(`init-container.sh` is the ≈150-line stripped descendant of `init-project.sh` defined in
Phase 1/Phase 3; it consumes `--type` and the project dir and emits only compose +
devcontainer + `.env`.)

**Why:** Enforces the seam. The two tools compose cleanly: `sciagent new project` makes the
*project*, `init-container.sh` makes the *container around it*. Neither imports the other's
templates.

**How:**
1. Detect whether scbio-docker is reachable (e.g. via `01_modules/SciAgent-toolkit` sibling
   or an env var) and print an absolute path if known, else a relative hint.
2. Do not hard-fail if scbio-docker isn't found — the project scaffold is valid on its own.

## Open ADRs

### ADR-2.1: Where do R/Python config templates live — SciAgent-toolkit or RNAseq-toolkit?
**Options:**
- **A.** Keep `analysis_config.yaml.template` + `config.R.template` + `pipeline.yaml.template`
  + `color_config.R.template` in SciAgent-toolkit `templates/project/analysis/`.
- **B.** Move all R/language-specific config into `RNAseq-toolkit` (language toolkit owns
  language config); SciAgent-toolkit scaffolds only the language-agnostic
  `analysis_config.yaml` (the shared YAML both R and Python read).
- **C.** Split: the language-agnostic `analysis_config.yaml.template` stays in
  SciAgent-toolkit (it's the contract); `config.R.template`/`color_config.R.template`
  (R-specific loaders) move to RNAseq-toolkit and are dropped in by its own scaffolder when
  registered as a submodule.

**Recommended:** C
**Why:** The shared `analysis_config.yaml` is a cross-language data contract (both R and
Python `read_yaml`/`safe_load` it) and defines the project's experimental design and
schemas — that is project context, so SciAgent-toolkit owns it. But `config.R` is an
R-only loader (`yaml::read_yaml`, `.libPaths`, `load_or_compute`) tightly coupled to the
RNAseq-toolkit's helper functions; baking it into the AI-harness toolkit re-couples
SciAgent-toolkit to R. Option C keeps SciAgent-toolkit language-neutral while still giving
analysts a working config out of the box once they attach RNAseq-toolkit. It also resolves
the path-block edit in 2.3: SciAgent-toolkit edits only the yaml `paths:`; RNAseq-toolkit
owns the R mirror.
**Blocking implementation:** yes — determines which templates 2.1 migrates into
SciAgent-toolkit vs leaves for the language toolkit.

### ADR-2.2: How large should the template-variable set be?
**Options:**
- **A.** Minimal — keep `{{PROJECT_ID}}`, `{{DATE}}` only; require users to fill species/genome
  by hand in the rendered config.
- **B.** Full init-project.sh parity — `{{PROJECT_ID/TITLE}}`, `{{DATE}}`, `{{PROJECT_TYPE}}`,
  `{{SPECIES}}`, `{{SPECIES_DB}}`, `{{GENOME_BUILD}}` with `--species/--genome/--title` flags
  and species→DB derivation.
- **C.** Full set, but sourced from a `project.yaml` answer-file instead of CLI flags.
**Recommended:** B
**Why:** The migrated templates *already* contain `{{SPECIES}}`/`{{GENOME_BUILD}}` tokens;
Option A would render broken config. Option C is cleaner long-term but adds an answer-file
format and parser this phase doesn't need — defer it. B matches what users expect from
init-project.sh and keeps the diff small (flags + a `case` for species derivation).
**Blocking implementation:** yes — 2.5 cannot be written without fixing the token set.

### ADR-2.3: Does `sciagent new project` own git init + submodules, or is that separate?
**Options:**
- **A.** `sciagent new project` always inits git and registers submodules (init-project.sh
  default-ish behavior).
- **B.** Git/submodules are **optional flags** (`--git`, `--with-submodules`), default off.
- **C.** Git is entirely a separate command (`sciagent new git`/out of scope); scaffold never
  touches git.
**Recommended:** B
**Why:** `sciagent new project` is frequently run inside an *already-initialized* repo (you
clone your project, then scaffold into it), where an automatic `git init` is wrong and a
forced initial commit is hostile. But submodule registration is genuinely project-level work
that belongs to this tool (per the seam), so excising it entirely (C) would lose
functionality init-project.sh had. Optional flags give both: clean default, full
init-project.sh parity on request. This is also why the next-steps text (2.8) mentions the
flags rather than running them.
**Blocking implementation:** no — can ship the structure scaffold first and add `_new_git()`
behind flags after.

### ADR-2.4: One canonical directory convention — resolve the three-way conflict?
**Options:**
- **A.** lowercase `00_data/ … 03_results/` (init-project.sh + analysis_config.yaml).
- **B.** Capitalized `00_Data/ … 03_Results/` (current AGENTS.md.template).
**Recommended:** A
**Why:** Two of three sources (init-project.sh, analysis_config.yaml) and all the R config
path constants already use lowercase; only the AGENTS.md template uses capitals. Standardize
on lowercase and fix AGENTS.md.template's inline "Directory structure" block to match the
phase-based layout from 2.3. Cheaper migration, matches existing real projects.
**Blocking implementation:** no — but must be settled before 2.2/2.3 write the AGENTS.md
"Directory structure" section.

## Dependencies
- **Depends on:** Phase 1 (defines the scbio-docker ↔ SciAgent-toolkit seam and the
  `init-container.sh` split that 2.8 hands off to).
- **Enables:**
  - Phase 3 (container substrate / `init-container.sh` consuming `--type` and the scaffolded
    dir).
  - Phase 5 (`handoff`/reasoning agents writing into `docs/_internal/sessions/` and
    `reasoning/` — the namespace this phase creates).
  - Phase 6 (software-tool project type — this phase ships its skeleton; Phase 6 fleshes out
    pyproject/test scaffolding).

## Breaking Changes
- `context.md` no longer produced; replaced by `docs/_internal/scientific-context.md`. Any
  agent/doc referencing `context.md` must be updated (Phase 5).
- `03_results/{checkpoints,plots,tables}` flat layout replaced by phase-based
  `03_results/{objects,NN_phase/…}`. Existing projects are not migrated automatically.
- `sciagent new project` signature grows flags (`--type`, `--species`, `--genome`,
  `--title`, `--git`, `--with-submodules`); output changes from 3 files to a full tree.
- `scbio-docker/templates/base/` loses its project-structure templates (moved here);
  `init-project.sh` is superseded by `sciagent new project` + `init-container.sh`.
- `_subst` drops `{{SKILL_NAME}}` from the project path; `{{PROJECT_NAME}}` becomes a
  temporary deprecated alias for `{{PROJECT_ID}}`.

These are all acceptable — single user, no downstream consumers; clean modularity wins.

## Estimated Scope
- `lib/sciagent/new.sh`: `_new_project()` ~35 → ~120 lines; `_subst` ~10 → ~30 lines; new
  `_new_git()` + `_derive_species_db()` ~90 lines. Net ≈ **+200 lines** in new.sh.
- Template migration: ~10 files `git mv`'d scbio-docker → SciAgent-toolkit; ~6 new
  `.template`/`.gitkeep` files for `docs/_internal/` and the phase-based results tree.
- `templates/AGENTS.md.template`: rewrite "Project context" + "Directory structure" sections
  (~30 lines changed).
- New template trees: `templates/project/_common/`, `templates/project/analysis/`,
  `templates/project/software-tool/` (skeleton).
- scbio-docker side: net **−450 lines** as `init-project.sh` (602) collapses toward
  `init-container.sh` (~150) — counted in Phase 1/3, noted here for the seam.
- **Total Phase-2 delta:** roughly +250 / −60 lines within SciAgent-toolkit, plus the
  cross-repo template moves.
