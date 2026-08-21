# Phase 6: Multi-Project-Type Support

## Summary

Today the toolkit assumes one shape of project:
- a single-cell / bulk RNA-seq **analysis** repo with the `00_data/ 01_modules/ 02_analysis/ 03_results/` layout.
That assumption is baked into the scaffold templates (`AGENTS.md.template` hardcodes `00_Data/ 01_Scripts/ 03_Results/`),
into `analysis_config.yaml.template` (GSEA thresholds, MSigDB databases, contrast schemas),
and into every role in `roles/`, which only ship analysis skills.

But the real working pattern has two genuine structural shapes. An **analysis** project runs
a biological question from raw data to published results (`00_data/ 01_modules/ 02_analysis/ 03_results/`).
A **software-tool** project is a standalone, packageable library or CLI (`src/ tests/ docs/ pyproject.toml`).
The user switches between these modes routinely — `pathway-explorer` lives inside `DC_hum_verse/01_modules/`
and is a fully-fledged Python package, not an analysis tree.

A third pattern exists — the **umbrella** — but it does not require a separate `--type`.
An umbrella is an analysis project whose `01_modules/` may include peer *analysis projects*
(pinned as root-level submodules) in addition to the usual software-tool toolkits, and whose
integration layer (`integration/`, the functional equivalent of `02_analysis/`) consumes
only the published `03_results/` surfaces of its children. DC-nexus is the canonical instance.
See §6.6 for the structural plan.

This phase makes **project type** a first-class concept with exactly two `--type` values:
`analysis` and `software-tool`. It threads a `--type` selector through `sciagent new`,
ships per-type scaffold template sets, adds a `software-tool` role, resolves nested-toolkit
activation, and extends `docs/_internal/` per type. It also documents the umbrella as a named
layout variant of `analysis` (no separate scaffold needed). This phase depends on Phase 5
(init-project.sh ownership moving to SciAgent-toolkit) and Phase 3 (the `docs/_internal/`
namespace).

## Current State

**Scaffolding is monolithic and analysis-only.**

- `lib/sciagent/new.sh` (141 lines). `_new_project()` (lines 40–66) copies exactly three files —
  `AGENTS.md`, `CLAUDE.md`, `context.md` — with `{{PROJECT_ID}}`/`{{DATE}}` substitution.
  No directory structure is created. No `--type` flag exists.
- `templates/` holds four templates total: `AGENTS.md.template`, `CLAUDE.md.template`,
  `context.md.template`, `analysis_config.yaml.template`. All four assume an analysis project.
- There is no template namespace for project *type*.

**Roles are all analysis roles.**

`base.yaml` ships ~50 analysis skills; `architect.yaml` is the only software-oriented role
but explicitly carries no day-to-day software-eng skills (no testing, packaging, or API-design
skill) and is meant as an overlay over an analysis base, not a standalone tool-dev role.

**Stack depth is hard-capped at 2.** `lib/sciagent/activate.sh` lines 66–69: three+ positional
args is a hard error. `docs/architecture.md` codifies this as a deliberate invariant.

**Activation is single-context per checkout.** There is no notion of "activate context X for
the parent and context Y for a child subdir."

**`docs/_internal/` (from Phase 3) is analysis-shaped.** No provision for software-tool
internal docs (design records, benchmark logs).

## Changes

### 6.1 Define an explicit project-type taxonomy

**Files:** `docs/architecture.md` (new "Project types" section), `templates/project/` (new tree)

**What:** Two first-class project types with canonical directory trees, one documented layout
variant, and three explicitly-catalogued-and-deferred types.

#### First-class types (scaffolded by `sciagent new project --type <t>`)

| Type | One-liner | Canonical top-level layout |
|------|-----------|----------------------------|
| `analysis` | A scientific analysis project. | `00_data/ 01_modules/ 02_analysis/ 03_results/<NN_phase>/ docs/ docs/_internal/` |
| `software-tool` | A standalone, packageable software library or CLI. | `src/ tests/ docs/ examples/ docs/_internal/` + `pyproject.toml`/`package.json` |

#### Umbrella — layout variant of `analysis` (no separate `--type`)

An **umbrella** project is an `analysis` project whose integration scope is cross-project:

- Root-level submodules are other *analysis* projects (e.g. `DC_Dictionary/`, `DC_hum_verse/`).
  They are gitlinked at the repo root, **not** inside `01_modules/`.
- `01_modules/` holds the usual software-tool toolkits (e.g. `SciAgent-toolkit`) that serve
  the umbrella's own integration work.
- The `02_analysis/` role is played by an `integration/` directory whose scripts may only
  consume `<child>/03_results/` (the published surface); they never import from child
  `02_analysis/` or intermediate checkpoints.
- No separate `--type umbrella` exists. Initialize with `--type analysis`, then add child
  submodules and create `integration/` by hand. See §6.6 for the DC-nexus structural plan.

#### Deferred / catalogued-only (ADR-6.4)

`pipeline` (Nextflow/Snakemake), `paper` (manuscript + figures), `data-package` (curated dataset
+ loader). Real categories, but no concrete present-day instance demands a scaffold template now.
Named here to reserve the vocabulary; no `templates/project/pipeline/` directory is created.

**Why:** Two `--type` values cover the full real portfolio: every analysis project and every
packaged toolkit. The umbrella is a structural variant of `analysis` — its distinction is
*what it puts in `01_modules/` and what it calls `02_analysis/`*, not a different scaffold
shape. Naming it explicitly as a variant gives it a canonical description without bloating the
`--type` flag surface.

**How:**
1. Add a "Project types" section to `docs/architecture.md` with the taxonomy table above,
   the umbrella variant note, and the deferred-types catalogue.
2. Establish `templates/project/{_common,analysis,software-tool}/` as the template root
   (ADR-X1 resolution from the Phase 0 audit). `_common/` holds type-agnostic files
   (`CLAUDE.md.template`, the `docs/_internal/` skeleton, shared `.gitignore` seed).
3. State the invariant: `_new_project()` renders `_common/` first, then overlays `<type>/`.
   Umbrella projects have no dedicated template — they start from `analysis` and diverge by hand.

### 6.2 Thread `--type` through `sciagent new project`

**Files:** `lib/sciagent/new.sh`, `templates/project/analysis/`, `templates/project/software-tool/`

**What:** `sciagent new project [<dir>] --type <analysis|software-tool>` (default `analysis`)
selects which template set is materialized. `_new_project()` creates the canonical directory
tree and copies the type's templates with substitution. Re-running in an existing project is
non-destructive (existing files are skipped).

**How:**
1. Add `--type <t>` argument parsing to `_new_project()`; validate against `{analysis,software-tool}`;
   default `analysis`; reject unknown types with the list of valid values. (ADR-6.1: use flag,
   not subcommand.)
2. Reorganize templates into `templates/project/_common/` and `templates/project/<type>/`:
   - `_common/`: `CLAUDE.md.template` (`@AGENTS.md`), `docs/_internal/reasoning/.gitkeep`,
     shared `.gitignore` seed.
   - `analysis/`: `AGENTS.md.template`, `context.md.template` (pointer form, per ADR-X2),
     `docs/_internal/scientific-context.md.template`, `analysis_config.yaml.template`.
   - `software-tool/`: `AGENTS.md.template` (tool guidance, no analysis idioms),
     `tool_config.yaml.template` (~30 lines: package name, language, version, lint/test commands),
     `README.md.template` (Features / Installation / Usage / Development / License).
3. Directory skeletons per type:
   - `analysis`: `00_data/{raw,processed,references}`, `01_modules/`, `02_analysis/config/`,
     `03_results/`, `docs/_internal/{reasoning,sessions,scratch}/`, `logs/`.
   - `software-tool`: `src/`, `tests/`, `docs/`, `examples/`, `docs/_internal/{reasoning,design,benchmarks}/`.
4. Extend `_new_project()`'s copy loop to walk `templates/project/_common/` then
   `templates/project/<type>/` rather than the hardcoded three-file list.
5. Update `cmd_new()` usage text to document `--type`.

### 6.3 Nested-toolkit activation model

**Files:** `lib/sciagent/activate.sh` (doc only), `docs/architecture.md`

**What:** A child `software-tool` under `01_modules/<tool>/` gets its own independent
activation — its own `.claude/`, `.agents/`, and `AGENTS.md` managed block rooted at the child
directory — rather than a third stack tier on top of the parent. Activation is per-directory-root;
the two roots are siblings in context, not nested in the stack.

**Why:** The depth-2 cap is a clean invariant worth keeping (ADR-6.3). The nested-toolkit case
does not require a third tier: parent analysis and child tool are *different working contexts*.
When working on `pathway-explorer`, you want the software-tool role, not "analysis base + tool overlay."

**How:**
1. Document in `docs/architecture.md`: "Activation root = CWD. A child toolkit is activated
   by `cd`-ing into it and running `sciagent activate` there; it gets its own depth-≤2 stack
   independent of the parent."
2. When `sciagent new project --type software-tool` is run inside an existing project's
   `01_modules/`, emit an informational note: "Scaffolded a child software-tool under <parent>.
   Activate its context by `cd <dir> && sciagent activate <role>`."
3. Do not add cross-root awareness; auto-switching on `cd` is a shell/harness concern, not
   a sciagent-core concern (catalogued in ADR-6.3 option C for a future shell-hook phase).

### 6.4 Add software-tool role and align roles to types

**Files:** `roles/software-tool.yaml` (new), `roles/architect.yaml` (doc note only),
`docs/architecture.md`

**What:** Add a `software-tool` role as a standalone base for tool development. Map types to
default roles (used by `_new_project()`'s "Next:" hint).

**How:**
1. Author `roles/software-tool.yaml` as a standalone base role:
   - Agents: `code-reviewer`, `docs-librarian`, `doc-curator`, `handoff` (domain-neutral).
   - Skills: `skill-creator` (meta). Dedicated software-eng skills (packaging, testing, api,
     cli design) are a catalogued follow-on (ADR-6.5/A — ship now, author skills later).
   - Commands: `commit`.
   - Output style: unset (compose `architect` overlay for mentor style).
2. Type → default role mapping:
   - `analysis` → `base`
   - `software-tool` → `software-tool`
3. No `meta-project` entry — umbrella projects use `base` as their role (same as analysis).

### 6.5 Extend `docs/_internal/` for software-tool projects

**Files:** `templates/project/software-tool/` scaffold, `docs/architecture.md`

**What:** Software-tool projects get a `docs/_internal/` namespace tuned to software work:
- `docs/_internal/reasoning/` — decision traces / why-not logs (shared with analysis).
- `docs/_internal/design/` — design records, API drafts, ADRs for the tool itself
  (replaces `sessions/` from the analysis variant).
- `docs/_internal/benchmarks/` — benchmark results, profiling logs.
- `docs/` — public-facing API docs, usage guides.

Analysis projects keep the Phase 3 shape (`reasoning/ sessions/ scratch/`). `reasoning/` is
universal. The handoff agent targets `sessions/` for analysis and `design/` for software-tool
(type-conditional target — implement as a content edit to the agent's output-path instruction).

**How:**
1. In 6.2's directory skeleton for `software-tool`, create
   `docs/_internal/{reasoning,design,benchmarks}` with a one-line `README.md` seed in each.
2. Update Phase 3's namespace description in `docs/architecture.md` to note the type-conditional
   shape: `reasoning/` universal; `{sessions,scratch}` analysis; `{design,benchmarks}` tool.

### 6.6 Umbrella layout — DC-nexus structural plan

**Files:** `docs/architecture.md` (umbrella variant section), DC-nexus `.gitmodules` and
`.devcontainer/` (deferred — implement after SciAgent refactor phases complete)

**What:** DC-nexus is the canonical umbrella instance. It currently has no `01_modules/` at
the umbrella root; this section plans the structural changes to bring it into conformance with
the umbrella variant layout without disrupting the ongoing analysis work.

**Why:** The `integration/` scripts need toolkit access (SciAgent roles, `load_or_compute`,
visualization helpers) without coupling to DC_Dictionary's or DC_hum_verse's own submodule
copies. The umbrella should be self-contained: its own `01_modules/SciAgent-toolkit` pin,
its own toolkit activation, and clean env-var injection for all mounted paths.

#### Planned structural changes to DC-nexus

**1. Add `01_modules/` at the umbrella root**

```
DC-nexus/
├── DC_Dictionary/          ← gitlinked submodule (analysis child)
├── DC_hum_verse/           ← gitlinked submodule (analysis child)
├── DC_mouse_cancer/        ← gitlinked submodule (pending Opus audit)
├── 01_modules/
│   └── SciAgent-toolkit/  ← NEW: gitlinked submodule (software-tool toolkit)
├── integration/            ← umbrella's 02_analysis/ equivalent
│   ├── 00_data/
│   ├── 01_scripts/
│   ├── 02_analysis/
│   └── 03_results/
├── .devcontainer/
├── docs/
└── ...
```

Naming discipline: analysis children live at the root (not inside `01_modules/`);
`01_modules/` holds only software-tool toolkits serving the umbrella's own work.

**2. Register umbrella SciAgent-toolkit in `.gitmodules`**

```gitconfig
[submodule "01_modules/SciAgent-toolkit"]
    path = 01_modules/SciAgent-toolkit
    url = git@github.com:tony-zhelonkin/SciAgent-toolkit.git
    branch = main
```

This is an independent pin at the umbrella level — bumped separately from the children's
own pins. No cross-submodule coupling.

**3. Update devcontainer environment**

Add to both `dev-core` and `dev-archr` service `environment:` blocks in
`.devcontainer/docker-compose.yml`:

```yaml
- SCIAGENT_ROOT=/workspaces/DC-nexus/01_modules/SciAgent-toolkit
```

This lets `integration/` scripts call `${SCIAGENT_ROOT}/bin/sciagent` without hardcoding paths,
consistent with the `DC_DICTIONARY_ROOT` / `DC_HUM_VERSE_ROOT` pattern already in place.

**4. `integration/` stays as-is (no rename)**

`integration/` is the functional `02_analysis/` for the umbrella. Renaming it to `02_analysis/`
would obscure its cross-project identity and break the consumption contract documented in
`integration/README.md`. Keep the name; document the equivalence in `docs/architecture.md`.

**5. Implementation timing**

Do **not** implement these changes until:
1. SciAgent refactor phases P1–P5 are complete (SciAgent must own the scaffold before we add
   a new instance of it).
2. DC_mouse_cancer receives its dedicated Opus audit and wiring decision.
3. A fresh `sciagent activate <role>` is run at the DC-nexus root to initialize the
   umbrella's own `.claude/` context from the new `01_modules/SciAgent-toolkit`.

---

## Open ADRs

### ADR-6.1: `--type` flag vs separate `sciagent new <type>` commands?
**Options:** A. `new project --type <t>` · B. `new analysis` / `new tool` · C. Both.
**Resolved: A.** `cmd_new()` dispatches on `<kind>` (project/role/skill/agent); `--type` is
an orthogonal modifier on `project`, not a new kind. Consistent with P2's `_new_project()`
signature. **Blocking 6.2: yes — accept A.**

### ADR-6.2: Shared `reasoning/` vs fully per-type `docs/_internal/` namespaces?
**Options:** A. `reasoning/` universal; rest per-type · B. Identical shape for all types ·
C. Fully disjoint.
**Resolved: A.** `reasoning/` is genuinely type-agnostic (agents write there regardless of
type). `sessions/scratch` vs `design/benchmarks` is a real semantic difference. **Not blocking.**

### ADR-6.3: Does the nested-toolkit pattern require lifting the depth-2 stack cap?
**Options:** A. Keep depth-2, independent per-root activation · B. Lift to depth-3 ·
C. Keep depth-2 + shell-hook auto-switch on `cd`.
**Resolved: A.** Depth-2 is load-bearing; child tool does not want parent's analysis skills
in its stack. Two independent activations at two roots covers the case cleanly without touching
`stack_walk()`. Shell-hook auto-switch (C) is a future harness concern, not sciagent-core.
**Blocking 6.3: yes — accept A.**

### ADR-6.4: Which types are first-class?
**Options:** A. `analysis` + `software-tool` + `meta-project` live; rest catalogued ·
B. Also ship `pipeline` now · C. Only `analysis` + `software-tool`; umbrella as informal.
**Resolved: C (revised).** `analysis` and `software-tool` are the only `--type` values.
`meta-project` (as originally defined — "analysis hosting software-tool children") is simply
what every analysis project with `01_modules/` toolkits already is; it adds no distinct scaffold.
The *umbrella* pattern (analysis hosting *analysis* children) is a real structural variant,
documented explicitly in §6.1 and §6.6 but not a `--type`. This is a deliberate reduction from
the original ADR-6.4/A recommendation. **Not blocking.**

### ADR-6.5: Does `software-tool` ship its own software-eng skills now, or compose `architect`?
**Options:** A. Ship now with `skill-creator` + architect composition; author eng skills later ·
B. Block until skills are authored · C. Don't add a role; use `architect` over `base`.
**Resolved: A.** B inflates scope. C fails because `architect` has no day-to-day software skills
and is an overlay, not a base. **Not blocking.**

---

## Dependencies

- **Depends on:**
  - **Phase 5** (init-project.sh ownership moved to SciAgent-toolkit): 6.2 assumes all
    scaffolding lives in SciAgent.
  - **Phase 3** (`docs/_internal/` namespace + context.md decomposition): 6.5 extends the
    namespace Phase 3 defines; 6.2's analysis skeleton uses the decomposed form.
  - **ADR-X1** (P0 audit): template root layout must be `templates/project/{_common,<type>/}`
    — resolved before 6.2 lands.

- **Enables:**
  - A future software-eng skills content phase (packaging/testing/api/cli skills from ADR-6.5).
  - Shell-hook auto-context-switching (ADR-6.3 option C, deferred).
  - Onboarding deferred types (`pipeline`/`paper`/`data-package`) — the `--type` plumbing and
    `templates/project/<type>/` convention are the extension points.
  - DC-nexus umbrella structural changes (§6.6) — gated on P1–P5 completion.

## Breaking Changes

- **`templates/` reorganized into `templates/project/{_common,analysis,software-tool}/`.**
  The flat `templates/AGENTS.md.template` etc. move. External callers with hardcoded paths
  break. Acceptable (single user, no downstream consumers).
- **`sciagent new project` now creates a full directory tree**, not three files.
- **`analysis_config.yaml.template` is no longer the universal config.** Software-tool projects
  get `tool_config.yaml.template` instead.
- **`meta-project` is removed as a taxonomy entry.** No template existed for it anyway.
  Projects that would have been called meta-projects are just `analysis` projects. No migration
  needed.
- No change to `stack_walk()` / activation semantics (ADR-6.3/A).

## Estimated Scope

- `lib/sciagent/new.sh`: +~90/−~25 lines. Net ~+65.
- `templates/project/analysis/`: moved files + context.md pointer split (Phase 3 payload).
- `templates/project/software-tool/`: ~3 new template files (~130 lines) + `_internal` seeds.
- `templates/project/_common/`: ~2 shared files.
- `roles/software-tool.yaml`: new file, ~40 lines.
- `docs/architecture.md`: new "Project types" + umbrella variant + nested-activation rule +
  Phase 3 namespace amendment. ~+90 lines.
- **Total:** ~8 files touched/created, net ~+330 lines. No changes to the stack/activation
  engine. Software-eng skills (ADR-6.5) are explicitly out of scope.
