# Phase 0: Architecture Audit — Refactor 2026-06-02

Audit of the six-phase plan to cut a clean seam between **scbio-docker** (container
substrate) and **SciAgent-toolkit** (project context management). This document is the
synthesis and gatekeeper for the refactor: it states the target system, orders the work,
flags inter-phase conflicts, ranks the open decisions, names what the plans forgot, and
records the highest-risk changes.

Phase docs audited:
- `phase-1-scbio-docker-strip.md`
- `phase-2-new-project-scaffold.md`
- `phase-3-docs-namespace.md`
- `phase-4-artifact-structure.md`
- `phase-5-internals-refactor.md`
- `phase-6-project-types.md`

---

## 1. Executive Summary

This refactor severs a single overgrown responsibility — "set up a project" — that today
is wrongly co-located in `scbio-docker/scripts/init-project.sh` (601 lines), and splits it
along a clean seam. After the refactor, **scbio-docker** is a pure container substrate:
Dockerfiles, image definitions, R/Python environment specs, build scripts, and a thin
`init-container.sh` (~150 lines) that renders only `docker-compose.yml`, `devcontainer.json`,
and `.env` for a target directory. It knows nothing about analysis trees, results layouts,
docs namespaces, config schemas, or AI context. **SciAgent-toolkit** becomes the single
authority for everything project-level: it owns the directory scaffold, the `docs/_internal/`
reasoning namespace, the phase-based artifact layout, language config templates, project-type
selection, and the AI harness (roles/skills/agents/commands). The teaching-platform domain
leaves scbio-docker entirely (it belongs to scbio-instruct).

The clean system the user gets is a two-tool composition with no overlap: `sciagent new
project --type {analysis,software-tool}` materializes a complete, opinionated project — a
numbered analysis tree *or* a `src/`+`tests/` tool layout, the `docs/_internal/{reasoning,
sessions,scratch}` namespace, a decomposed scientific-context pointer, and phase-based
`03_results/NN_stage/{tables,figures}` — and then prints the one scbio-docker command that
layers the dev container on top. Neither tool imports the other's templates. Internally,
SciAgent-toolkit's brittle bits (`stack_walk()`'s triplicated arms, `inject.sh`'s post-hoc
collision detection, the hand-maintained module load graph) are de-duplicated and made
auditable, so the tool that now carries all the project-setup load is sound enough to carry it.

Concretely, the user moves from "one 601-line script does everything, mixing container and
project concerns, hardcoding three conflicting directory conventions, dumping every working
note into `context.md`/`AGENTS.md`/`CLAUDE.md` until they hit 500–755 lines" to "two
composable tools with crisp boundaries, one canonical lowercase phase-based layout, a
sanctioned home for every kind of agent output, and first-class support for both analysis
repos and the software-tool sub-toolkits they host." Breaking changes are accepted
throughout — single user (Anton), no downstream consumers — so every phase optimizes for
clean modularity over compatibility.

---

## 2. Phase Dependency Graph

### Declared dependencies (from each phase's own `Dependencies` section)

| Phase | Depends on | Enables | Hard blockers |
|-------|-----------|---------|---------------|
| **P1** scbio-docker strip | none | P2, P5, final cleanup | — |
| **P2** new-project scaffold | P1 (seam + `init-container.sh` handoff) | P3, P5, P6 | ADR-2.1, ADR-2.2 |
| **P3** docs namespace | P4 (results layout, *soft*), P2 (full scaffold) | doc-curator enforcement, lean role files | none |
| **P4** artifact structure | P1/P2 (templates moved to toolkit), P3 multi-type (*soft*) | deterministic harness paths, validation methodology | none |
| **P5** internals refactor | P1 only (mandate); parallel to P2–P4 | future verb/inject work | none |
| **P6** project types | P5 (init-project ownership moved), P3 (`_internal` namespace) | software-eng skills phase, shell-hook future, deferred types | ADR-6.1, ADR-6.3 |

### Narrative

**P1 is the keystone.** It is the only phase with no dependencies and it unblocks every
other phase: it cuts the seam and physically moves the templates. Nothing project-level can
land correctly until P1 has established *which repo owns what*.

**P2 is the second keystone.** It turns the seam into code: `_new_project()` grows from a
3-file stub into the real scaffolder. P3, P4, and P6 all assume P2's "full scaffold copies a
template tree" mechanism exists. In practice P2 must land before P3/P4/P6 can be *wired*
(they supply template payloads that P2's copy loop materializes), even though P3 and P4
describe their payloads as if independent.

**P3 and P4 are tightly coupled and mutually referential** (see Coherence Check §3). P4
defines the canonical `03_results/` layout; P3 *consumes* that layout in `AGENTS.md.template`
and the `context.md` pointer. They should be treated as one work unit or P4-then-P3.

**P5 is the independent track.** It depends only on P1 for *mandate* (SciAgent now carries
the load, so its internals must be sound), not for code. It touches `lib/sciagent/stack.sh`,
`inject.sh`, `collisions.sh`, `block.sh`, `bin/sciagent` — none of which P2/P3/P4 edit except
`new.sh`. It can run fully in parallel.

**P6 is last.** It depends on P5 ("init-project ownership moved" — though strictly that is P1)
and P3 (`_internal` namespace). It reorganizes `templates/` into `templates/types/<type>/`,
which **collides directly with P2's `templates/project/<type>/` layout** (see §3, the
single most important conflict in this audit).

### Recommended execution order

```
P1  (strip scbio-docker, move templates)           ── must be first
 │
 ├─► P5  (internals refactor)                       ── parallel track, P1-mandate only
 │
 └─► P2  (scaffold engine in new.sh)                ── second keystone
        │
        ├─► P4  (artifact structure / payload)      ── do P4 before P3
        │     │
        │     └─► P3  (docs namespace / payload)     ── consumes P4's layout
        │
        └─► P6  (project types)                      ── must reconcile templates/ layout with P2
```

Rationale for **P4 before P3**: P3's `AGENTS.md.template` directory section and `context.md`
pointer both *reference* the phase-based layout (`01_qc/`, `02_programs/`, `objects/`). If P3
lands first it must forward-declare paths P4 hasn't finalized. P4 explicitly owns the
`03_results/` tree and the `AGENTS.md.template` directory-diagram rewrite (P4 change 4.8) —
and so does P3 (change 3.2 step 2). One of them must go first and own that diagram; P4 is the
natural owner because it defines the shape. See §3 for this overlap.

Rationale for **P2 before P6**: both restructure `templates/`. P2 proposes
`templates/project/{_common,analysis,software-tool}/`; P6 proposes
`templates/types/{analysis,software-tool}/`. These are two names for the same thing. Whichever
lands first sets the convention; the second must conform, not re-move the files. **Decide the
template-root layout once, in P2, and have P6 inherit it.**

---

## 3. Coherence Check

### 3.1 P2 ↔ P6 — template directory layout (DIRECT CONFLICT)

**Conflict.** P2 §2.1/§2.7 establish the template root as
`templates/project/_common/`, `templates/project/analysis/`,
`templates/project/software-tool/`. P6 §6.1/§6.2 establish it as
`templates/types/analysis/`, `templates/types/software-tool/`, with the shared
`CLAUDE.md.template` left flat at `templates/CLAUDE.md.template` (no `_common/`).

These are incompatible directory conventions for the *same files*. P2 introduces a `_common/`
overlay tier ("render `_common/` first, then overlay `<type>/`"); P6 has no `_common/` and
keeps only `CLAUDE.md.template` shared at top level. If implemented as written, P6 would
re-move files P2 just placed, and the `_common/`-overlay rendering logic P2 builds would be
orphaned by P6's flat-per-type model.

**Resolution (must decide before P2 lands):** Pick one root name and one sharing model.
Recommend **`templates/project/<type>/` + `templates/project/_common/`** (P2's form) because
the `_common/` overlay is the cleaner mechanism for the genuinely-shared payload
(`AGENTS.md`, `CLAUDE.md`, the `docs/_internal/` skeleton, `_subst` machinery) that both
P2 and P3 and P6 agree is type-agnostic. Rewrite P6 §6.2 to target `templates/project/`
not `templates/types/`, and to use `_common/` for shared files instead of leaving
`CLAUDE.md.template` flat. This is a naming reconciliation, not a design change. **New ADR
needed — see §4 (ADR-X1).**

### 3.2 P3 ↔ P4 — who owns the `AGENTS.md.template` directory diagram (OVERLAP)

**Overlap, low risk if sequenced.** Both phases rewrite the `## Directory structure` block
of `AGENTS.md.template`:
- P4 §4.8 ("Correct the AGENTS.md.template directory diagram") replaces the capitalized flat
  `03_Results/{objects,tables,figures}` with the lowercase phase-based 4.1 tree.
- P3 §3.2 step 2 ("replace the hardcoded `00_Data/01_Scripts/...` tree with the
  namespace-aware structure") also rewrites it, and P3 §3.2 step 2 says to "remove the
  duplicated flat `objects/tables/figures` listing (defer authoritative tree to
  `analysis_config.yaml` + the results phase)."

These are *consistent in intent* (both want lowercase, phase-based, deferring the
authoritative tree to config) but both edit the same fenced block. If done independently they
will conflict in the file. **Resolution:** P4 owns the authoritative directory diagram (it
defines the shape); P3 owns the *namespace* additions (the `docs/_internal/` routing table,
the "never grow this file" rules). Sequence P4→P3 and have P3 add only the `docs/` namespace
content, treating the `03_results/` block as already-correct from P4.

### 3.3 P2 ↔ P3 — context.md decomposition: two slightly different target shapes (MINOR INCONSISTENCY)

P2 §2.2 renames `context.md.template` → `docs/_internal/scientific-context.md.template` and
**removes `context.md` from scaffold output entirely** ("Remove from scaffold output:
`context.md` (top-level monolith)").

P3 §3.3 **keeps `context.md` at project root** as a ~30-line pointer ("`context.md` stays at
the project root (agents and tools look for it there) but becomes a ~30-line pointer
document").

**Conflict:** P2 deletes `context.md`; P3 retains it as a pointer. P3's reasoning is stronger
("tools/agents/users expect `context.md` at root"; ADR-3.1 recommends the pointer over
deletion). **Resolution:** adopt P3's model — `context.md` survives as a pointer, plus
`docs/_internal/scientific-context.md` holds the body. Update P2 §2.2 to *not* remove
`context.md` from scaffold output but to render it as the pointer P3 defines. This is a real
contradiction that must be reconciled, not a phrasing nuance: the scaffold's file list
differs between the two phases.

### 3.4 P3 ↔ P6 — `docs/_internal/` shape diverges by type (CONSISTENT, needs handoff alignment)

P3 defines `docs/_internal/{reasoning,sessions,scratch}` (analysis-shaped). P6 §6.5 extends
this for software-tool to `docs/_internal/{reasoning,design,benchmarks}` — `reasoning/`
shared, `{sessions,scratch}` analysis-only, `{design,benchmarks}` tool-only. This is
*consistent and additive* (P6 ADR-6.2 option A). One alignment point: the **`handoff` agent**
is retargeted by P3 §3.4 to write `docs/_internal/sessions/`, but P6 §6.5 step 4 says the
handoff agent must write to `design/` for software-tool projects. The handoff agent therefore
needs **type-awareness** that P3 alone does not give it. P3 hardcodes `sessions/`; P6 makes it
conditional. **Resolution:** P3 should retarget handoff to `sessions/` as the analysis
default but leave a note that P6 will make the target type-conditional; P6 must actually
implement that branch. Flag: P6 §6.5 step 4 calls this "a doc/convention change, not a code
change" — but if the handoff agent's output path is hardcoded in `agents/analysis-base/
handoff.md`, making it type-conditional *is* a content edit to that agent file. Slight
under-scoping in P6.

### 3.5 P1 ↔ P2 ↔ P5 — the `init-container.sh` name (CONSISTENT, naming drift)

P1 §1.2/ADR-1.1 produces a stripped `init-project.sh` and floats renaming it to
`render-devcontainer.sh` (ADR-1.1 recommends option C). P2 §2.8 calls the stripped script
`init-container.sh`. P5 dependencies reference "Phase 5 (init-project.sh ownership moving to
SciAgent-toolkit)" — but P1 is what moves it, not P5; P6 repeats "Phase 5 (init-project.sh
split)." **Three names in play** (`init-project.sh` stripped / `render-devcontainer.sh` /
`init-container.sh`) and **a mis-attribution** (P5/P6 credit P5 with the init-project split
that P1 actually owns). Not a design conflict, but the name and the owning-phase must be
pinned so P2's next-steps text and P6's dependency claims point at the right artifact.
**Resolution:** pick `init-container.sh` (P2's name, most descriptive of the new contract),
and correct P5/P6 dependency prose to credit **Phase 1** with the split.

### 3.6 P2 ↔ P4 — config-template ownership is unresolved across both (SHARED OPEN QUESTION)

P2 ADR-2.1 (blocking) asks where `config.R.template`/`config.py`/`pipeline.yaml.template`
live (SciAgent vs RNAseq-toolkit), recommending option C (split: yaml in SciAgent, R loaders
in RNAseq-toolkit). P4 §4.3 then **writes new `config.R.template` and `config.py.template`
content** (the `stage_dir()` resolvers) — but if ADR-2.1/C is accepted, `config.R.template`
**doesn't live in SciAgent-toolkit at all** (it moves to RNAseq-toolkit). So P4's edits to
`config.R.template` would land in the wrong repo. **Resolution:** ADR-2.1 must be resolved
*before* P4 §4.3, and P4 §4.3's R-config edits must follow it to whichever repo wins. The
language-agnostic `pipeline.yaml.template` (`stages:`, `schemas:`, `paths:`) stays in
SciAgent per all options; the R/Python loaders' home is gated by ADR-2.1. This is the most
consequential cross-phase coupling on the config axis.

### 3.7 Directory-convention three-way conflict — resolved consistently

P2 ADR-2.4, P3, P4 §4.8, and P6 all independently flag and resolve the same defect: three
conflicting directory conventions exist today (`00_Data/` capitalized in
`AGENTS.md.template`; `00_data/` lowercase in init-project.sh and `analysis_config.yaml`;
`03_results/{checkpoints,plots,tables}` vs `{checkpoints,tables,plots,interactive}`). All four
phases converge on **lowercase `00_data/`** and **phase-based `03_results/`**. This is
coherent — no conflict — but it means **four phases all edit `AGENTS.md.template`'s directory
block** (P2 §2.5/§2.7, P3 §3.2, P4 §4.8, P6 §6.2). That file is a hot contention point;
single-ownership must be assigned (recommend: P4 owns the `03_results/` tree, P3 owns the
`docs/` namespace, P6 owns the per-type variants, P2 owns nothing in this file beyond the
`_subst` tokens).

### 3.8 Summary of coherence findings

| # | Phases | Type | Severity | Must resolve before |
|---|--------|------|----------|---------------------|
| 3.1 | P2 ↔ P6 | template root layout conflict | **High** | P2 lands |
| 3.2 | P3 ↔ P4 | AGENTS.md diagram double-edit | Medium | P3 lands |
| 3.3 | P2 ↔ P3 | context.md kept vs deleted | **High** | P2 lands |
| 3.4 | P3 ↔ P6 | handoff agent type-awareness | Medium | P6 lands |
| 3.5 | P1/P2/P5/P6 | script name + owning-phase drift | Low | P2 next-steps text |
| 3.6 | P2 ↔ P4 | config-template home (ADR-2.1) | **High** | P4 §4.3 |
| 3.7 | P2/P3/P4/P6 | AGENTS.md ownership (converged) | Low | implementation |

---

## 4. ADR Cross-Reference

Across the six phases there are 22 ADRs plus the new ones this audit surfaces. Grouped by
when they must be settled.

### 4.1 Must-resolve BEFORE any implementation

These gate file moves or the scaffold engine; resolving them wrong means re-doing work.

**ADR-X1 (NEW — from §3.1): Single template-root layout for P2 and P6.**
*Question:* `templates/project/<type>/` + `_common/` (P2) **or** `templates/types/<type>/`
(P6)?
**Recommendation: `templates/project/{_common,analysis,software-tool}/` (P2's form).**
The `_common/` overlay is the right mechanism for the type-agnostic payload all three of
P2/P3/P6 share. Rewrite P6 to target this root. Blocking: yes — it determines where every
migrated template lands.

**ADR-X2 (NEW — from §3.3): Does `context.md` survive as a root pointer, or is it deleted?**
*Question:* P2 deletes it; P3 keeps it as a pointer.
**Recommendation: keep it as the P3 pointer.** Discoverability convention is strong, cost is
~30 lines, ADR-3.1 already recommends this. Update P2 §2.2 accordingly. Blocking: yes — it
changes the scaffold's output file list.

**ADR-2.1 (P2, declared blocking): Where do R/Python config templates live?**
*Question:* SciAgent vs RNAseq-toolkit for `config.R.template`/`config.py`.
**Recommendation: option C (split).** Language-agnostic `pipeline.yaml`/`analysis_config.yaml`
(the cross-language contract) stays in SciAgent; R/Python loaders move to RNAseq-toolkit.
Keeps SciAgent language-neutral. **Critical:** this gates P4 §4.3 (§3.6) — P4 must write its
`stage_dir()` config edits to whichever repo wins. Blocking: yes.

**ADR-2.2 (P2, declared blocking): Template-variable vocabulary size.**
*Question:* minimal vs full init-project.sh parity vs answer-file.
**Recommendation: option B (full parity).** Migrated templates already contain
`{{SPECIES}}`/`{{GENOME_BUILD}}`; minimal would render broken config. Blocking: yes — `_subst`
cannot be written without the token set fixed.

**ADR-6.1 (P6, declared blocking): `--type` flag vs `sciagent new <type>` subcommands.**
*Question:* `new project --type analysis` vs `new analysis`.
**Recommendation: option A (flag).** `cmd_new()` already dispatches on `<kind>`; `--type` is
an orthogonal modifier on `project`. Blocking: yes — 6.2 cannot be written until fixed, but
it also affects P2's `_new_project()` signature (P2 already assumes `--type`, so A is
consistent with P2). Resolve jointly with P2.

**ADR-6.3 (P6, declared blocking): Lift the depth-2 stack cap for nested toolkits?**
*Question:* depth-2 independent roots vs depth-3 additive vs cd-autoswitch.
**Recommendation: option A (keep depth-2, independent per-root activation).** Lifting to
depth-3 would force a rewrite of `stack_walk()` and reopen the injected-state question P5
deliberately closes. A is also the no-code-change-to-stack option and keeps P5 and P6 from
fighting over `stack.sh`. Blocking: yes — sets the entire nesting model.

### 4.2 Can be deferred to DURING implementation

Settle as the relevant change is written; wrong choice is cheap to revisit.

- **ADR-1.1 (P1): keep `init-project.sh` / rename to `render-devcontainer.sh` / library
  function.** Recommend the descriptive name `init-container.sh` (per §3.5); the
  library-function form (C) is a later cosmetic. Not blocking.
- **ADR-1.3 (P1): homes for `.vscode/settings.json`, `setup_claude_mcp.sh`, MCP/venv docs.**
  Recommend: `settings.json` → scbio-docker (tracks the image); MCP files → SciAgent
  (track the harness). Not blocking; decide at move time.
- **ADR-2.3 (P2): git init + submodules — flags or automatic.** Recommend option B (optional
  flags, default off) — scaffold is often run inside an already-cloned repo. Not blocking.
- **ADR-3.4 (P3): is `docs/_internal/` committed or gitignored?** Recommend option C
  (`reasoning/`+`sessions/` committed, `scratch/` gitignored). Not blocking.
- **ADR-3.5 (P3): handoff filename date-only vs date+time.** Recommend A (date-only, time on
  collision). Not blocking.
- **ADR-4.1 (P4): pre-populate stage folders or fully on-demand.** Recommend B (seed
  `01_qc/`+`02_eda/`). Not blocking.
- **ADR-4.2 (P4): fixed vs generic stage slugs.** Recommend A refined by C (generic for
  analysis; software-tool ships no `03_results/`). Coordinate with ADR-X1. Not blocking.
- **ADR-4.3 (P4): bake a `validation/` dir.** Recommend B (lazy-create, document the gate).
  Not blocking.
- **ADR-4.4 (P4): master-table home.** Recommend A (`03_results/master/`). Not blocking.
- **ADR-5.1 (P5): inject collision behavior.** Recommend C (warn default, `--force`,
  `--strict`). Not blocking — warn-and-proceed core lands first.
- **ADR-5.4 (P5): declarative verb→modules table.** Recommend B (declarative table). Not
  blocking.
- **ADR-6.2 (P6): shared `reasoning/` vs fully per-type internal namespaces.** Recommend A
  (`reasoning/` universal, rest per-type). Not blocking.
- **ADR-6.5 (P6): software-tool ships own skills now or composes architect.** Recommend A
  (ship now with `skill-creator` + architect composition; author eng skills later). Not
  blocking.

### 4.3 Nice-to-have / purely mechanical

- **ADR-1.2 (P1): drop SciAgent submodule from scbio-docker.** Recommend B (keep during
  refactor, drop at end). It's literally the last step; the refactor lives inside the
  submodule checkout. Mechanical.
- **ADR-3.1 (P3): context.md pointer vs full vs delete.** Subsumed by ADR-X2; recommend
  pointer.
- **ADR-3.2 (P3): golden path in templates vs roles/skills.** Recommend B (keep in
  guideline + skill, template points only). Strong recommendation; proceed.
- **ADR-3.3 (P3): AGENTS.md vs CLAUDE.md thinness.** Recommend A (CLAUDE.md = `@AGENTS.md`
  only). Proceed.
- **ADR-5.2 (P5): nameref vs eval vs globals for `_sw_record()`.** Recommend A (namerefs).
  Mechanical; bash 4.x already assumed.
- **ADR-5.3 (P5): `block_render_and_write` in stack.sh vs block.sh.** Recommend A (stack.sh;
  keep block.sh a pure leaf). Mechanical.
- **ADR-6.4 (P6): which types first-class.** ~~Recommend A~~ **Revised: option C.**
  `analysis` + `software-tool` only. `meta-project` removed — every analysis project with
  `01_modules/` toolkits already is one; a separate type adds no distinct scaffold. The
  *umbrella* pattern (analysis hosting peer *analysis* children, e.g. DC-nexus) is documented
  as a named layout variant of `analysis` in P6 §6.1 and §6.6, but carries no `--type` value.
  `pipeline`/`paper`/`data-package` remain catalogued-only.

---

## 5. Missing Pieces

What the six phases collectively do **not** address but the full refactor needs.

1. **No migration path for existing real projects.** P3 explicitly names two live projects
   (`13403-YD_Christina`, `AdaW_eWAT_WL_2025`) and P6 names two more (`DC_hum_verse`,
   `pathway-explorer`) as evidence, but every phase says "existing projects are not migrated
   automatically." That is fine as a *policy* but there is **no documented procedure** for the
   user to manually port those four projects to the new layout (flat `03_results/` → phase-
   based, `context.md` monolith → pointer + scientific-context, root `handoff_*.md` →
   `sessions/`). At minimum a `docs/_internal/migration-existing-projects.md` checklist is
   needed, or those four projects diverge permanently from the scaffold the tooling assumes.

2. **No test coverage for the new scaffold engine.** P5 adds tests for its *own* changes
   (shadow-CSV, marker-boundary grep, no-`exit`-in-libs, per-verb smoke). But P2's rewritten
   `_new_project()` — the single biggest behavioral change in the refactor — ships with **no
   test plan**. There is no test that `sciagent new project --type analysis` produces the
   expected tree, that `--type software-tool` produces `src/`+`tests/`, that `_subst` leaves
   no unresolved `{{…}}` tokens (P2 §2.5 step 3 mentions a warning grep but not a test), or
   that re-running is non-destructive. This is a coverage gap for the highest-risk new code.

3. **No documentation of the new `sciagent new project` interface for end users.** P6 §6.2
   step 6 updates `cmd_new()` *usage text*, but the user-facing docs — scbio-docker's
   `README.md`/`QUICKSTART.md` (P1 rewrites these but to point at "`sciagent new project`"
   generically) and SciAgent's `CLAUDE.md`/`docs/architecture.md` — are not given a worked
   end-to-end example of the new two-step flow (`sciagent new project` → `init-container.sh`).
   The seam is described; the *workflow a human follows* is not written down in one place.

4. **`pipeline.yaml` vs `analysis_config.yaml` — two names, unclear which is canonical.** P2
   and P4 reference `pipeline.yaml.template`; P3 and P6 reference `analysis_config.yaml`. The
   current repo has *both* a `pipeline.yaml.template` (in scbio-docker's `templates/base/`)
   and an `analysis_config.yaml.template` (in SciAgent's `templates/`). No phase states
   whether these are the same file, whether one supersedes the other, or how `stages:`/
   `schemas:` (P4 adds to `pipeline.yaml`) relate to the GSEA/MSigDB content (in
   `analysis_config.yaml`). This must be unified or the scaffold ships two overlapping config
   files.

5. **The `.devcontainer/.env` MCP-key responsibility is split with no owner test.** P1 §1.2
   step 7 keeps MCP API-key lines in `.env` "documented as consumed by setup-ai.sh later,"
   but `setup-ai.sh` itself is never mentioned as surviving or being updated — and the CLAUDE.md
   workflow references it heavily. No phase clarifies whether `setup-ai.sh` still exists post-
   refactor, or whether `sciagent activate` subsumes it. The AI-tooling bootstrap path is left
   ambiguous.

6. **No rollback / partial-state story for `inject`'s new inline collision check.** P5 §5.2
   asserts the check is read-only and runs before any side effect, so "no rollback needed."
   But P5 §5.5 simultaneously hardens *all* side-effecting calls to fail loudly — meaning a
   `block_write` or `manifest_*` failure mid-inject is now a hard error. No phase describes
   what state the project is left in if inject fails *after* writing symlinks but *before*
   updating the manifest/block (the symlinks are orphaned). The half-mounted state P5.5 warns
   about is identified but not given a recovery command (`sciagent repair`?).

7. **Cross-repo commit choreography is hand-waved.** P1/P2 note that moving files from
   scbio-docker to SciAgent-toolkit is "two separate commits (one `git rm`, one `git add`)
   because the toolkit is a submodule." But the submodule is dropped at the *end* (ADR-1.2/B),
   so during the refactor scbio-docker pins a specific SciAgent SHA. No phase specifies the
   submodule-bump ordering: if SciAgent gains the templates in commit A, scbio-docker must
   bump its submodule pointer to include A before its own `git rm` is coherent. This
   choreography is real work and unspecified.

8. **No definition of "done" / acceptance criteria for the seam.** There is no phase that
   says "the refactor is complete when: (a) `grep -r '03_results\|00_data\|context.md' scbio-
   docker/` returns only container-irrelevant hits; (b) `sciagent new project` + `init-
   container.sh` produces a working VS Code dev container end-to-end; (c) all four real
   projects have a documented migration status." A final verification phase is missing.

---

## 6. Risk Assessment

The five highest-risk changes, ranked.

### Risk 1 — Template-root layout conflict between P2 and P6 (§3.1)

**What could go wrong:** P2 builds the `_common/`-overlay rendering engine against
`templates/project/`; P6 independently moves everything to `templates/types/` with no
`_common/`. Result: the overlay logic is orphaned, files are moved twice, and the scaffold
silently renders the wrong (or missing) shared payload. Because both phases touch the same
files, a merge of independently-implemented P2 and P6 produces a broken `_new_project()` that
may still *appear* to run.
**Mitigation:** Resolve ADR-X1 before P2 lands. Pin one root name (`templates/project/`) and
one sharing model (`_common/`). Make P6 a *conformance* change, not a re-layout. Add a test
(see Missing Piece 2) asserting the rendered tree for each `--type`.

### Risk 2 — Config-template home (ADR-2.1) lands P4's edits in the wrong repo (§3.6)

**What could go wrong:** If P4 §4.3 rewrites `config.R.template` in SciAgent-toolkit but
ADR-2.1/C later moves R config to RNAseq-toolkit, the `stage_dir()` resolver work is stranded
in the wrong repo and the cross-language `paths:` contract drifts between two homes. Scripts
that `source config.R` then write to the old flat paths and every file silently lands in the
wrong directory — a *data-loss-shaped* failure (results written where nothing looks for them).
**Mitigation:** Resolve ADR-2.1 *before* P4 §4.3. Make `pipeline.yaml`'s `paths:` block the
single contract both languages read; treat R/Python loaders as thin consumers placed per
ADR-2.1. Add a test that a freshly-scaffolded project's `stage_dir("01_qc","figures")`
resolves to `03_results/01_qc/figures/` in both R and Python.

### Risk 3 — The `.gitignore` recursive rule for nested results (P4 §4.2)

**What could go wrong:** P4 correctly identifies that the current `03_results/*/*` glob only
matches one level and silently commits two-level-deep figures. But the *replacement* recursive
pattern (`03_results/**` + `!03_results/**/` + `!03_results/**/.gitkeep` + hard-ignore
`03_results/objects/**`) is exactly the class of `.gitignore` rule that is easy to get subtly
wrong — negated globs interacting with directory-vs-file matching is a notorious footgun. A
wrong rule either commits multi-GB `.h5ad` checkpoints (repo bloat / push failure) or ignores
the `.gitkeep`s (scaffold dirs vanish on clone).
**Mitigation:** P4 §4.2 step 2 already mandates `git check-ignore -v` verification — make this
a *test*, not a manual step: assert a figure is ignored, its `.gitkeep` tracked, an
`objects/*.h5ad` hard-ignored, across at least two nesting depths.

### Risk 4 — `_sw_record()` nameref extraction silently changes shadow-CSV semantics (P5 §5.1)

**What could go wrong:** The triplicated `stack_walk()` arms encode subtle "append previous
provider to shadow CSV in accumulation order" logic. Bash namerefs have a self-aliasing
footgun (a nameref local colliding with the caller's variable name produces silent wrong
behavior, not an error). If the extraction gets the accumulation order or the alias hygiene
wrong, the managed AGENTS.md block renders with wrong "shadowed by" provenance — a *quiet*
correctness regression that won't crash and may not be noticed until provenance is audited.
**Mitigation:** P5 §5.1 step 5 already specifies the locking test (base+overlay shadowing a
skill, agent, and command simultaneously, asserting the multi-shadow CSV). Make that test a
**precondition** — write it against the *current* triplicated code first (capture the golden
output), then refactor and assert byte-identical. P5 §5.5's "byte-identical managed block for
the same stack" claim should be a CI guard.

### Risk 5 — Cross-repo / submodule commit choreography (Missing Piece 7)

**What could go wrong:** During the refactor scbio-docker vendors SciAgent as a submodule
pinned to a SHA. Moving templates out of scbio-docker requires SciAgent to gain them *first*,
then scbio-docker to bump its submodule pointer, then `git rm` its copies. Done out of order,
there is a window where the templates exist in *neither* committed state (scbio-docker removed
them, submodule pointer predates their addition) — a clone at that commit can't scaffold
anything. Compounded by ADR-1.2's "drop submodule at the end," the bump/rm/deinit sequence is
four-way ordered and unspecified.
**Mitigation:** Write the choreography explicitly as a Missing-Piece deliverable: (1) commit
template additions in SciAgent; (2) bump scbio-docker submodule pointer; (3) `git rm` in
scbio-docker referencing the bumped pointer; (4) at refactor end, deinit submodule. Never `git
rm` before the pointer bump. Consider doing all moves while the submodule is *present* (per
ADR-1.2/B) precisely to keep both trees in one working checkout.

---

## 7. Clean Architecture Statement

> **scbio-docker** is a container substrate. It provides versioned Docker images (R + Bioconductor + Python single-cell stacks), the R/Python environment definitions, build scripts, and one thin command — `init-container.sh` — that renders `docker-compose.yml`, `devcontainer.json`, and `.env` into any target directory. It knows nothing about project structure, results layouts, documentation, configuration schemas, or AI context.
>
> **SciAgent-toolkit** owns everything project-level. `sciagent new project [--type analysis|software-tool]` materializes the full project scaffold: the directory tree, the phase-based `03_results/NN_stage/{tables,figures}` artifact layout, the `docs/_internal/{reasoning,sessions,scratch}` reasoning namespace, the decomposed scientific-context pointer, and the language-agnostic config contract. `sciagent activate <role>` mounts a role — a bundle of skills, sub-agents, and commands — into `.claude/` and `.agents/`. Two first-class types exist: `analysis` and `software-tool`. The **umbrella** pattern (an analysis project hosting peer analysis projects as root-level submodules) is a documented layout variant of `analysis`, not a separate type.
>
> The two tools compose and never import each other. To start a project, a user runs `sciagent new project --type <type> <dir>`, then `init-container.sh <dir> --data-mount <...>`, opens the directory in VS Code, reopens in the dev container, and runs `sciagent activate <role>`. The project is the unit; the container wraps it; the harness reasons within it.
