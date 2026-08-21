# Phase 1: scbio-docker Scope Reduction

## Summary

This phase strips `scbio-docker` down to its single legitimate responsibility: being a
**container substrate**. After this phase, the repo provides Docker images, a
docker-compose template, a devcontainer.json template, R/Python environment definitions,
and build scripts — and nothing else. Everything that knows about *project directory
structure*, *analysis config*, *AI context*, or *docs namespaces* moves out (to
SciAgent-toolkit) or is deleted outright (teaching-platform → owned by scbio-instruct).

Concretely, the two big movers are `scripts/init-project.sh` (601 lines, ~75% of which is
project scaffolding) and `templates/base/` (the entire analysis-project tree). The
container-only residue of `init-project.sh` shrinks to ~150 lines whose sole job is to
render `docker-compose.yml` and `devcontainer.json` from templates with data-mount and
UID substitution. The analysis tree, config templates, results layout, and docs structure
are handed to SciAgent-toolkit, which already owns project context and will own scaffolding
after this refactor. Breaking changes are explicitly acceptable (single user, no downstream
consumers); we optimize for a clean seam, not backwards compatibility.

## Current State

Repo root: `/data1/users/antonz/pipeline/scbio-docker` (branch `teaching-platform`).

**`scripts/init-project.sh` — 601 lines.** Symlinked from repo root as `init-project.sh`.
Mixed responsibilities, line ranges (approx., from the 601-line file):
- Lines 1–155: arg parsing, usage, interactive prompts (data mounts, git, submodules,
  resource limits, **species/genome configuration**), template validation.
- Lines 156–175: project-dir existence check + `mkdir`.
- Lines ~180–200: **creates the analysis directory tree** (`00_data/{raw,processed,references}`,
  `01_modules/.ref`, `02_analysis/{config,helpers,scripts,notebooks}`,
  `03_results/{checkpoints,plots,tables}`, `docs/{raw,ai-generated/...,plan/phase-1}`, `logs`).
- Lines ~200–215: copies `.vscode/settings.json` and project `.gitignore`.
- Lines ~215–290: **renders project docs and config** via `sed` token substitution
  (`docs/README.md`, `docs/plan/README.md`, root `README.md`, `notes.md`,
  `02_analysis/config/{config.R,pipeline.yaml,color_config.R}`) — all carry
  `{{SPECIES}}`, `{{SPECIES_DB}}`, `{{GENOME_BUILD}}`, `{{PROJECT_NAME}}` tokens.
- Lines ~292–340: **container config rendering** — `devcontainer.json` (sed) and
  `docker-compose.yml` (inline Python heredoc for the multi-line `{{DATA_MOUNTS}}` block).
  **This is the only part that is genuinely container-domain.**
- Lines ~340–375: copy devcontainer `scripts/`, poststart sanity fallback, copy
  `MCP_AUTH_SETUP.md` / `PYTHON_VENV_GUIDE.md`.
- Lines ~375–420: create/update `.devcontainer/.env` (LOCAL_UID/GID + MCP API keys).
- Lines ~420–435: `.gitkeep` seeding for the analysis tree.
- Lines ~435–520: git init + **submodule attachment** (`RNAseq-toolkit`, `SciAgent-toolkit`)
  with SSH→gh-CLI fallback (`add_submodule_with_fallback`).
- Lines ~520–601: success summary + next-steps text referencing `setup-ai.sh`,
  `context.md`, `pipeline.yaml`.

**`templates/base/` — 31 files** (`find templates/base -type f`). Audit by domain:
- Container-domain (STAYS, possibly relocated to `templates/devcontainer/`):
  - `.devcontainer/devcontainer.json.template`
  - `.devcontainer/docker-compose.yml.template`
  - `.devcontainer/.env.example`
  - `.devcontainer/scripts/{setup_claude_mcp.sh,source_env.sh}` — *ambiguous, see ADR-1.3*
  - `.devcontainer/{MCP_AUTH_SETUP.md,PYTHON_VENV_GUIDE.md}` — *ambiguous, see ADR-1.3*
  - `.vscode/settings.json` — *ambiguous: R/Python interpreter paths are container-domain,
    but file lives in project root; see ADR-1.3*
- Project-domain (MOVES to SciAgent-toolkit):
  - `00_data/{raw,processed,references}/.gitkeep`
  - `01_modules/.ref/.gitkeep`
  - `02_analysis/config/{color_config.R.template,config.R.template,pipeline.yaml.template}`
  - `02_analysis/{helpers,notebooks,scripts}/.gitkeep`
  - `03_results/{checkpoints,plots,tables}/.gitkeep`
  - `docs/README.md.template`, `docs/plan/README.md.template`, `docs/plan/phase-1/.gitkeep`,
    `docs/notes.md.template`, `docs/raw/.gitkeep`, `docs/ai-generated/{research,vignettes}/.gitkeep`
  - `README.md` (project root README with `{{PROJECT_NAME}}` etc.)
  - `.gitignore` (project-level, excludes `00_data/`, `03_results/`)
  - `logs/.gitkeep`

**`teaching-platform/` — untracked dir, 33 files** (`01_base/`, `01_host/`, `02_image/`,
`03_hub/`, `04_provision/`, `docs/`, plus build logs). NOT in git index (`git ls-files
teaching-platform/` is empty). The *docs* counterpart `docs/teaching-platform/` (13 `.md`
files) IS tracked and is staged for deletion in the current working tree (see `git status`
`D docs/teaching-platform/...`).

**`QUICKSTART.md` — 66 lines.** Half single-cell quick start (stale: references nonexistent
`basic-rna` template), half teaching-platform operations runbook (JupyterHub, NVMe mounts,
provisioning). **`GEMINI-PROMPT.md` — 158 lines.** Entirely a teaching-platform plan-revision
prompt; zero container-substrate relevance.

**`README.md` — 84 lines.** Already mostly aligned with the strip ("containerization-only")
but has bugs: line 8 references `basic-rna` template (doesn't exist), line 29 same, line 60
typo "comb bio tooling", line 77 "setup-ai deprecated" dangling. References `init-project.sh
--with-submodules` (line 67) which will change.

**`.gitmodules` — 1 submodule:**
```
[submodule "toolkits/SciAgent-toolkit"]
	path = toolkits/SciAgent-toolkit
	url = https://github.com/tony-zhelonkin/SciAgent-toolkit
```
So scbio-docker currently vendors SciAgent-toolkit at `toolkits/SciAgent-toolkit/`.

**`CLAUDE.md` — 27,690 bytes.** Documents v0.5.4. Already claims "containerization-only" but
its "Project Setup Pattern" section (the whole `init-project.sh creates ...` tree, `03_results`
layout, `docs/_internal` non-existence, `setup-ai.sh` workflow) describes the to-be-removed
behavior and must be rewritten.

## Changes

### 1.1 Remove teaching-platform from scbio-docker entirely

**Files:**
- `teaching-platform/` (untracked dir, 33 files) — delete from working tree
- `docs/teaching-platform/` (13 tracked `.md` files) — already staged `D`; finalize
- `QUICKSTART.md` sections 2 & 3 (teaching-platform runbook)
- `GEMINI-PROMPT.md` (entire file)

**What:** Excise all teaching-platform code, docs, and prompts. scbio-instruct owns this
domain; scbio-docker keeps zero knowledge of JupyterHub/DockerSpawner/student provisioning.

**Why:** Teaching-platform is a *deployment* of a teaching image, not container-substrate
definition. It belongs in scbio-instruct (per ground truth). Keeping it here violates the
clean seam and bloats the repo with NVMe-mount runbooks and Pushover monitoring scripts.

**How:**
1. Confirm nothing in `teaching-platform/` is needed by image builds:
   `grep -rl "teaching-platform" docker/ scripts/ .devcontainer/` (expect no hits in the
   build path; only docs/runbooks reference it).
2. If any artifact is worth preserving for scbio-instruct, copy it there *first*
   (out of scope for this phase — coordinate with the scbio-instruct migration), then:
   `git rm -r --cached docs/teaching-platform/` (already staged; verify with `git status`).
3. Delete untracked code dir: `rm -rf teaching-platform/`.
4. Delete the teaching build logs at repo root referencing it if any
   (`build-archr-wrapper.log`, `build.log` are image logs — keep gitignored, not teaching).
5. Stage deletions; the commit message should state scbio-instruct now owns the platform.

### 1.2 Strip `scripts/init-project.sh` to container-only (~601 → ~150 lines)

**Files:** `scripts/init-project.sh`, root symlink `init-project.sh`.

**What:** Reduce `init-project.sh` to a single responsibility: **render the devcontainer +
compose config for a project directory** (data mounts, image version, resource limits,
UID/GID env). Remove: the analysis directory tree creation, all `02_analysis/config`
rendering, all docs/README/notes rendering, species/genome configuration, `.gitkeep`
seeding, submodule attachment, and the project-level `.gitignore`/`.vscode` copy *if*
ADR-1.3 routes those out.

**Why:** ~75% of the script is project scaffolding that belongs to SciAgent-toolkit. The
container repo must not encode `03_results/`, `pipeline.yaml`, or `docs/` layout. Keeping
only compose/devcontainer rendering makes the script's contract obvious: "given a project
path and mount spec, write a working dev container into `.devcontainer/`".

**How:**
1. **Keep** arg parsing for the container-relevant flags only:
   `<project-dir>`, `--data-mount KEY:PATH[:ro]` (repeatable), `--image-version vX.Y.Z`
   (default from `VERSION`), `--service dev-core|dev-archr` (new explicit flag; default
   `dev-core`), `--max-cpus`, `--max-memory`. **Remove** `--interactive` species prompts,
   `--git-init`, `--with-submodules`, and the `TEMPLATE` positional arg (only one template,
   and it's now SciAgent's concern). See ADR-1.1 for whether the command survives at all.
2. **Remove** the directory-tree loop (`for dir in 00_data/... 03_results/... docs/...`).
3. **Remove** all `sed`/heredoc blocks rendering `docs/README.md`, `docs/plan/README.md`,
   root `README.md`, `notes.md`, `02_analysis/config/{config.R,pipeline.yaml,color_config.R}`.
4. **Remove** `SPECIES`/`SPECIES_DB`/`GENOME_BUILD` derivation and all `{{SPECIES*}}` tokens.
5. **Remove** the `.gitkeep` seeding loop.
6. **Remove** `add_submodule_with_fallback`, `RNASEQ_TOOLKIT_URL`, `SCIAGENT_TOOLKIT_URL`,
   and the entire git-init + submodule block.
7. **Keep & keep working**: the `.devcontainer/` creation, `devcontainer.json` sed render
   (`{{PROJECT_NAME}}`, `{{SERVICE}}`), the Python heredoc that renders `docker-compose.yml`
   with `{{IMAGE_VERSION}}`, `{{PROJECT_NAME}}`, `{{MAX_CPUS}}`, `{{MAX_MEMORY}}`,
   `{{DATA_MOUNTS}}`, the `.devcontainer/scripts/` copy + poststart sanity fallback, and the
   `.devcontainer/.env` creation (LOCAL_UID/GID + MAX_CPUS/MAX_MEMORY). The MCP API-key lines
   in `.env` are container-runtime env, not AI scaffolding — keep them but documented as
   "consumed by setup-ai.sh later".
8. **Rewrite** the success/next-steps footer to drop `context.md`/`pipeline.yaml`/`setup-ai.sh`
   references; instead point at SciAgent-toolkit's `sciagent new project` for the rest.
9. Update `TEMPLATES_DIR` to point at the surviving container template dir (see 1.3 / ADR-1.3).
10. Keep the symlink `init-project.sh -> scripts/init-project.sh` (or remove if ADR-1.1
    deprecates the command).

**Resulting interface (container-only):**
```
init-project.sh <project-dir> [OPTIONS]
  --data-mount KEY:PATH[:ro]   (repeatable)
  --image-version vX.Y.Z       (default: VERSION file)
  --service dev-core|dev-archr (default: dev-core)
  --max-cpus N                 (default: 50)
  --max-memory NG              (default: 450G)
Outputs (into <project-dir>):
  .devcontainer/devcontainer.json
  .devcontainer/docker-compose.yml
  .devcontainer/.env
  .devcontainer/scripts/poststart_sanity.sh
```

### 1.3 Audit & relocate `templates/base/` — keep container, move project

**Files:** all of `templates/base/` (31 files).

**What:** Split `templates/base/` into:
- **STAY** in scbio-docker, relocated to `templates/devcontainer/`:
  `.devcontainer/devcontainer.json.template`, `.devcontainer/docker-compose.yml.template`,
  `.devcontainer/.env.example`, `.devcontainer/scripts/{source_env.sh}`, and the poststart
  sanity script.
- **MOVE** to SciAgent-toolkit `templates/`:
  `02_analysis/config/{config.R.template,pipeline.yaml.template,color_config.R.template}`,
  `docs/{README.md.template,plan/README.md.template,notes.md.template}`, project root
  `README.md`, the project `.gitignore`, and the `.gitkeep` skeleton for
  `00_data/`, `02_analysis/`, `03_results/`, `docs/`, `logs/`.
- **DECIDE** via ADR-1.3: `.vscode/settings.json`, `.devcontainer/scripts/setup_claude_mcp.sh`,
  `.devcontainer/{MCP_AUTH_SETUP.md,PYTHON_VENV_GUIDE.md}`.

**Why:** The container repo must not own the analysis tree, config schema (`pipeline.yaml`
with GSEA params, MSigDB collections), or docs namespace. Those are project-context, owned by
SciAgent-toolkit. Only the devcontainer/compose templates and the env stub are substrate.

**How:**
1. Create `templates/devcontainer/` in scbio-docker.
2. `git mv templates/base/.devcontainer/* templates/devcontainer/` (devcontainer.json,
   docker-compose.yml, .env.example, scripts/source_env.sh).
3. Move the project-domain files into SciAgent-toolkit `templates/` preserving subpaths
   (e.g. `templates/project/02_analysis/config/pipeline.yaml.template`). Coordinate exact
   destination layout with the SciAgent-toolkit scaffolding phase (Phase 2/3 of this refactor).
4. `git rm -r templates/base/` once everything is relocated.
5. Update `scripts/init-project.sh` `TEMPLATE_PATH` references to `templates/devcontainer/`.
6. Note for the receiving repo: `config.R.template` and `pipeline.yaml.template` carry the
   `03_results/{checkpoints,plots,tables}` flat layout (`DIR_CHECKPOINTS`, `DIR_TABLES`,
   `DIR_PLOTS` in `config.R.template`). Per ground truth, SciAgent-toolkit will rework these
   into phase-based artifact grouping (`01_qc/`, `02_programs/` each holding tables+figures)
   with checkpoints separated as state — that rework is **not** done here; this phase only
   *moves* the files so the rework happens in their new home.

### 1.4 Resolve the SciAgent-toolkit submodule relationship

**Files:** `.gitmodules`, `toolkits/SciAgent-toolkit/`.

**What:** Decide and implement whether scbio-docker keeps SciAgent-toolkit as a submodule at
all after the split. See ADR-1.2.

**Why:** Today the submodule exists so `init-project.sh --with-submodules` can re-attach
SciAgent-toolkit into projects and so this repo can reference toolkit scripts. After the
strip, init-project does no submodule attachment and the container repo has no functional
dependency on toolkit code. The only remaining justification is the *refactor itself* (this
plan lives under `toolkits/SciAgent-toolkit/docs/_internal/`).

**How (recommended path — ADR-1.2 option B):**
1. During the refactor, **keep** the submodule (the plan docs and the moved templates land in
   it; convenient to edit both repos from one checkout).
2. At the *end* of the refactor (a later phase), `git submodule deinit toolkits/SciAgent-toolkit`,
   `git rm toolkits/SciAgent-toolkit`, remove the stanza from `.gitmodules`, and delete the
   `[submodule ...]` entry from `.git/config`.
3. Document in CLAUDE.md that SciAgent-toolkit is a **sibling** repo, attached *per-project*
   (at `01_modules/SciAgent-toolkit/`), never vendored into the container repo.

### 1.5 Clean up root-level docs (QUICKSTART, GEMINI-PROMPT, README)

**Files:** `QUICKSTART.md`, `GEMINI-PROMPT.md`, `README.md`.

**What:**
- Delete `GEMINI-PROMPT.md` (teaching-platform-only).
- Rewrite `QUICKSTART.md` to cover *only* building the image and rendering a dev container;
  drop sections 2 (teaching platform) and 3 (maintenance/JupyterHub). Fix the stale
  `basic-rna` template reference.
- Fix `README.md`: remove `basic-rna` (lines 8/29 in the rendered README → in repo
  `README.md` lines 29, 67), fix the "comb bio tooling" typo (line 60), replace the dangling
  "setup-ai deprecated" (line 77) with a one-liner pointing at SciAgent-toolkit, and update
  the `init-project.sh --with-submodules` example to the new container-only interface.

**Why:** Root docs are the first thing a reader sees; they currently advertise a template that
doesn't exist and a teaching platform that's leaving. After the strip they must describe a
pure container substrate.

**How:**
1. `git rm GEMINI-PROMPT.md`.
2. Rewrite `QUICKSTART.md` to ~25 lines: build image → render container → reopen in VS Code →
   (optional) `sciagent new project` for project scaffolding.
3. Edit `README.md` per the list above. Verify no remaining `basic-rna`, `teaching-platform`,
   or `setup-ai deprecated` strings: `grep -nE 'basic-rna|teaching-platform|setup-ai deprecated' README.md QUICKSTART.md`.

### 1.6 Rewrite scbio-docker `CLAUDE.md` for the post-strip reality

**Files:** `CLAUDE.md` (repo root).

**What:** Remove the entire "Project Setup Pattern" section that documents the analysis tree,
`03_results` layout, `init-project.sh creates ...`, `setup-ai.sh` workflow, and the
`docs/_internal` namespace. Replace with a short "Project Scaffolding" section that says:
container config comes from `init-project.sh` (devcontainer/compose only); *all* project
structure, config, docs, and AI context come from SciAgent-toolkit (`sciagent new project`).

**Why:** CLAUDE.md is load-bearing for future agents working in this repo. If it keeps
describing the old monolith behavior, agents will recreate it. It must encode the seam.

**How:**
1. In CLAUDE.md, delete the "Project Setup Pattern" subsections that enumerate the
   `my-project/` tree, the `init-project.sh creates` block, and the `setup-ai.sh creates`
   block.
2. Update the "Build" / overview top matter only where it claims init-project does project
   scaffolding.
3. Add a "Boundary with SciAgent-toolkit" subsection stating the clean seam verbatim:
   scbio-docker = Dockerfiles, compose/devcontainer templates, env definitions, build
   scripts; SciAgent-toolkit = project scaffold, config templates, docs namespace, AI harness.
4. Update the `init-project.sh` documentation to the new flag set from 1.2.
5. Leave image-build, venv, and R-library sections untouched (still accurate and in-domain).

## Open ADRs

### ADR-1.1: Does scbio-docker keep a project-init command at all, or only provide templates that SciAgent-toolkit renders?
**Options:**
- **A — Keep a thin `init-project.sh`** in scbio-docker that renders devcontainer + compose,
  callable standalone.
- **B — Delete `init-project.sh` entirely**; ship only `templates/devcontainer/*.template`
  and have SciAgent-toolkit's `sciagent new project` call into scbio-docker (read its
  templates by path / `SCBIO_DOCKER` env) to render the container config.
- **C — Library function**: scbio-docker exposes a tiny `render-devcontainer.sh` *function*
  (not a project-init command) that SciAgent sources; scbio-docker never "inits a project".
**Recommended:** **C** (with B as the pragmatic interim).
**Why:** The clean seam says scbio-docker must not know about "projects". A *project* is a
SciAgent concept. scbio-docker should own one verb only: "render container config for a
directory". Naming it `init-project` leaks the project concept. A `render-devcontainer.sh`
that takes `(target_dir, image_version, service, mounts...)` and writes `.devcontainer/` is
the honest container-domain primitive; `sciagent new project` orchestrates it alongside the
analysis tree + AI harness. Interim: keep the thin script (B) working so nothing breaks
mid-refactor, rename/refactor to C in a later phase.
**Blocking implementation:** No (1.2 can land as a thin script; rename is cosmetic and can
follow once SciAgent's `new project` is wired up in the next phase).

### ADR-1.2: After the split, does scbio-docker still vendor SciAgent-toolkit as a submodule?
**Options:**
- **A — Keep submodule** at `toolkits/SciAgent-toolkit/` permanently.
- **B — Keep during refactor, drop at the end.**
- **C — Drop immediately.**
**Recommended:** **B.**
**Why:** Permanent vendoring (A) re-creates the coupling we're trying to remove — the
container repo would carry the entire AI toolkit. Dropping immediately (C) is wrong *right
now* because this very refactor plan and the relocated templates live inside the toolkit
checkout; editing both from one working tree is convenient until the moves are committed.
B keeps the seam clean in the final state while not blocking the migration. After drop,
SciAgent-toolkit is a sibling repo attached per-project at `01_modules/SciAgent-toolkit/`,
never inside scbio-docker.
**Blocking implementation:** No (the deinit is the *last* step of the overall refactor, not
this phase).

### ADR-1.3: Where do `.vscode/settings.json`, `setup_claude_mcp.sh`, and the MCP/venv docs live?
**Options (per artifact):**
- `.vscode/settings.json` — A: scbio-docker (it pins R/Python interpreter paths
  `/opt/venvs/base/bin/...` that are *image* facts) / B: SciAgent-toolkit (it's project-root
  editor config).
- `setup_claude_mcp.sh` + `MCP_AUTH_SETUP.md` + `PYTHON_VENV_GUIDE.md` — A: scbio-docker /
  B: SciAgent-toolkit.
**Recommended:** `.vscode/settings.json` → **A (scbio-docker)** because its content is
100% derived from image internals (interpreter paths, radian path, httpgd) and changes when
the image changes; it ships from `templates/devcontainer/`. `setup_claude_mcp.sh`,
`MCP_AUTH_SETUP.md`, `PYTHON_VENV_GUIDE.md` → **B (SciAgent-toolkit)** because MCP wiring is
AI-harness concern (already SciAgent's domain). `PYTHON_VENV_GUIDE.md` is borderline (it
documents the *image's* venvs) — recommend it move to scbio-docker `docs/` instead, not a
template.
**Why:** Decide by "does the content change when the *image* changes (scbio-docker) or when
the *AI harness/project* changes (SciAgent)?" `settings.json` tracks the image; MCP setup
tracks the harness.
**Blocking implementation:** No (1.2/1.3 can keep `.vscode/settings.json` flowing from the
relocated `templates/devcontainer/.vscode/` while the MCP files are simply dropped from the
container template copy step).

## Dependencies
- **Depends on:** none (this is the first phase; it intentionally lands before SciAgent
  gains its richer `sciagent new project`).
- **Enables:**
  - Phase 2 (SciAgent project scaffolding) — receives the moved `templates/base/` project
    tree, config templates, docs namespace, and owns the `03_results` → phase-based rework
    plus the `docs/_internal/{reasoning,sessions,scratch}` baking.
  - Phase(s) on `stack.sh`/`inject.sh`/`collisions.sh` debt — unblocked once scbio-docker is
    out of the project-scaffolding business and all scaffolding logic lives in one repo.
  - Final cleanup phase — submodule deinit (ADR-1.2) and `init-project.sh` → `render-
    devcontainer.sh` rename (ADR-1.1).

## Breaking Changes
- `init-project.sh` no longer creates the analysis tree, config files, docs, or `.gitkeep`s.
  Existing workflows that relied on a single command producing a full project break; project
  scaffolding now requires a second step (`sciagent new project`). **Acceptable** — single
  user, no downstream consumers.
- `--git-init`, `--with-submodules`, `--interactive` (species prompts), and the `TEMPLATE`
  positional argument are **removed** from `init-project.sh`.
- `templates/base/` is deleted; the only surviving templates are `templates/devcontainer/`.
- `GEMINI-PROMPT.md` and `teaching-platform/` (+ `docs/teaching-platform/`) are removed from
  the repo.
- `basic-rna` template references (already nonfunctional) are removed from README/QUICKSTART.

## Estimated Scope
- **`scripts/init-project.sh`**: ~601 → ~150 lines (−450).
- **`templates/base/`**: −31 files from scbio-docker (most re-homed in SciAgent-toolkit;
  ~5–6 container files relocated to `templates/devcontainer/`).
- **`teaching-platform/`**: −33 untracked files; **`docs/teaching-platform/`**: −13 tracked
  files; **`GEMINI-PROMPT.md`**: −158 lines.
- **`QUICKSTART.md`**: 66 → ~25 lines (−40).
- **`README.md`**: ~10 lines edited.
- **`CLAUDE.md`**: ~−150 lines (delete "Project Setup Pattern" section) / +30 lines
  (new "Boundary with SciAgent-toolkit" + updated init-project interface).
- **`.gitmodules`**: 1 stanza removed (deferred to final phase per ADR-1.2).
- **Net for scbio-docker:** roughly −700 lines of script/docs and −70+ files, leaving a
  pure container-substrate repo.
