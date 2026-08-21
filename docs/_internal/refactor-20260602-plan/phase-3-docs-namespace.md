# Phase 3: docs/_internal Namespace + Context Decomposition

## Summary

This phase gives agents and humans a **sanctioned place for working documents**. Today there is no home for decision traces, dated handoffs, or exploratory writeups, so they leak into `context.md`, `CLAUDE.md`, and a litter of root-level files (`ideas.md`, `research_notes.md`, `handoff_*.md`, `docs_notes.md`). The observed result: `AGENTS.md` ballooned to 634 lines, `context.md` to 755, and `CLAUDE.md` to 532 — each absorbing whatever had nowhere else to go. We fix the root cause by baking a `docs/_internal/` namespace into every project scaffold, declaring authoritatively in `AGENTS.md.template` where each kind of output belongs, and demoting `context.md` from a 477–755-line monolith to a ~30-line pointer.

The namespace is organized along **two independent axes**:

1. **Privacy axis** — `docs/_internal/` is the outer wrapper for everything not ready to be shared publicly. Both AI-generated and human-authored private docs live here. Public-facing material lives at `docs/` root and is opt-in (you consciously move or create things there).
2. **Authorship axis** — within `_internal/`, everything an agent writes goes into `ai-generated/` so you can audit what agent teams actually produced. Human-authored private notes sit directly in `_internal/` without any further namespace (free-form, project-specific).

The biological framework (hypothesis, gene sets, contrast families, methodology rationale) moves out of `context.md` into `docs/_internal/scientific-context.md`. After this phase, an agent doing a literature search knows it goes in `docs/_internal/ai-generated/research/`, a handoff goes in `docs/_internal/ai-generated/sessions/YYYY-MM-DD_slug.md`, a phase plan goes in `docs/_internal/ai-generated/plans/YYYY-MM-DD__name/` — never inline in a context or harness file. `AGENTS.md` stays lean because the namespace gives growth somewhere to go.

## Current State

Templates live at `toolkits/SciAgent-toolkit/templates/`:

- `AGENTS.md.template` — 1766 bytes, ~70 lines. Has: project context block, a `SCIAGENT:ROLES` managed block, four "Critical rules", a hardcoded `00_Data/01_Scripts/02_Analysis/03_Results` directory tree, and a Notes section. Has **no** guidance on where agents should write working documents.
- `CLAUDE.md.template` — 11 bytes. Single line: `@AGENTS.md`.
- `context.md.template` — 735 bytes. Full document: Scientific Question, Datasets table, Hypotheses, Analysis Goals checklist, Broad Direction, Key References, Notes. This is the seed that grew to 477–755 lines in real projects.
- `analysis_config.yaml.template` — 5726 bytes. Defines `paths:` with the **flat** results layout: `checkpoints/`, `tables/`, `plots/`, `interactive/` under `03_results/`.

Real-project evidence (ground truth for the problem):

- **13403-YD_Christina**: `AGENTS.md` = 634 lines (became a project-specific monolith), `context.md` = 477 lines, `CLAUDE.md` = 216 lines. Already independently invented `docs/_internal/reasoning/`, `docs/stages/`, dated session folders (`docs/2026-05-25/`), `ANALYSIS_STAGES.md`, `TROUBLESHOOTING.md` — convergent evidence that the namespace is needed, but ad-hoc per project.
- **AdaW_eWAT_WL_2025**: `context.md` = 755 lines (full TE-inflammation hypothesis, gene sets, contrast families). Working docs (`ideas.md`, `docs_notes.md`, `research_notes.md`, a handoff file) scattered at project root. `AGENTS.md` = 66 lines (good — minimal), `CLAUDE.md` = 532 lines (golden-path patterns, checkpoint cache pattern, master-tables schema).

The `handoff` agent (`agents/analysis-base/handoff.md`) currently writes `handoff_YYYYMMDD_HHMMSS.md` **to project root** and archives old ones to `.handoff_archive/`. This is part of the root-clutter problem and must be retargeted to `docs/_internal/sessions/`.

The golden-path / checkpoint pattern is already documented canonically at `docs/guidelines/checkpoint_caching.md` (the `load_or_compute()` function) — but in real projects it was copy-pasted into each `CLAUDE.md`, which is why those files grew to 500+ lines.

Note on dependency: this phase consumes the `03_results/` phase-based artifact layout defined in the results-restructure phase. Where this doc references `01_qc/`, `02_programs/`, `objects/` it assumes that layout is already chosen.

## Changes

### 3.1 Define the docs namespace authoritatively

**Files:** new `templates/docs/` skeleton (scaffolded by `sciagent new project`); documented in `AGENTS.md.template` (see 3.2) and `docs/architecture.md`.

**What:** Establish one canonical directory contract, baked into every new project as empty dirs with `.gitkeep` plus a top-level `docs/_internal/README.md` that states the contract. The namespace:

```
docs/
├── _internal/                 # private by default; gitignore entry in public repos
│   ├── README.md              # states the contract below
│   ├── scientific-context.md  # SSoT for biology (see 3.3)
│   ├── ai-generated/          # ALL agent output — never human prose here
│   │   ├── research/          # lit search, web crawl, paper synthesis
│   │   │   └── YYYY-MM-DD_topic.md
│   │   ├── explorations/      # codebase audits, data crawls, agent-team findings
│   │   │   └── YYYY-MM-DD_what.md
│   │   ├── sessions/          # dated handoffs, RESUME.md
│   │   │   └── YYYY-MM-DD_desc.md
│   │   ├── decisions/         # ADR-style, slug.md, referenced forward from code
│   │   │   └── slug.md
│   │   └── plans/             # devops / phase plans (AI-generated, internal by default)
│   │       └── YYYY-MM-DD__phase-name/
│   └── [free-form human notes]  # your private intelligence, audit notes — no fixed schema
└── [public-facing, opt-in; free-form, project-specific]
    architecture/, vignettes/, deck/, notes.md, README.md, plan/, guides/, methods.md
```

The public area is intentionally unstructured. You move or create things there deliberately when ready to share. Anything you are not sure about stays in `_internal/`.

Authoritative "what goes where" contract:

| Dir | Holds | Lifecycle | Example filename |
|-----|-------|-----------|------------------|
| `docs/_internal/ai-generated/research/` | Lit search, web crawl, paper synthesis, AI-compiled reading notes | Dated, append-only | `2026-05-25_te-family-survey.md` |
| `docs/_internal/ai-generated/explorations/` | Codebase audits, data crawls, agent-team findings ("what we found") | Dated, append-only | `2026-05-25_codebase-audit.md` |
| `docs/_internal/ai-generated/sessions/` | Dated session handoffs, RESUME.md (replaces root `handoff_*.md`) | Dated; no archive dir needed | `2026-05-25_integration-working.md` |
| `docs/_internal/ai-generated/decisions/` | ADR-style decision records, why-not logs, "considered X, rejected because Y" | Slug-named, stable | `ADR-003-checkpoint-format.md` |
| `docs/_internal/ai-generated/plans/` | Devops / phase plans generated by agents (internal reasoning, not for public sharing by default) | Dated subdir per plan | `2026-05-25__phase-2-qc/` |
| `docs/_internal/scientific-context.md` | Full biological framework, hypotheses, gene sets, contrast families | Curated SSoT | (single file) |
| `docs/_internal/[free-form]` | Human-authored private notes, intelligence, audit notes | Free-form, project-specific | `strategy.md`, `field-notes.md` |
| `docs/[public area]` | Anything ready to share; free-form | Opt-in, curated | `architecture/`, `vignettes/`, `deck/` |

**Why:** Two independent axes drove the design. First, **privacy**: not all docs should be public-facing by default — devops plans contain internal reasoning, research notes are working intelligence, session handoffs are kitchen. A single `_internal/` gitignore entry handles the entire private footprint in public repos. Second, **authorship traceability**: in computational biology you need to be able to audit exactly what an agent team found, not just trust the orchestrator's summary. Separating `ai-generated/` from free-form human notes makes that audit path unambiguous — `docs/_internal/ai-generated/explorations/YYYY-MM-DD_what.md` is a first-class artifact, not a hidden tmp file.

Both observed real projects independently grew `docs/ai-generated/` (14616-DM) and `docs/ai-research/` (DC_hum_verse) — the namespace is a discovered need. Baking it in stops per-project reinvention and gives agents somewhere to park output so `context.md`/`AGENTS.md`/`CLAUDE.md` stop absorbing everything.

**How:**
1. Create `templates/docs/` mirroring the tree above, each leaf dir holding a `.gitkeep`.
2. Write `templates/docs/_internal/README.md` containing the contract table verbatim.
3. Write `templates/docs/_internal/scientific-context.md.template` (the long-form seed; see 3.3).
4. Wire `sciagent new project` to copy `templates/docs/` into the new project (this is part of the broader "new project copies only 3 files → full scaffold" fix tracked in the scaffolding phase; here we only own the `docs/` subtree).
5. Add the `_internal/` convention to `.gitignore` policy decision (ADR-3.4): by default `_internal/` is **committed** (it is the project's reasoning record), not ignored.

### 3.2 Add a "Documentation Namespace" section to AGENTS.md.template

**Files:** `templates/AGENTS.md.template`

**What:** Insert a new section that tells agents WHERE to write outputs and what NOT to inline. Also replace the hardcoded `00_Data/01_Scripts/...` tree with the namespace-aware structure and a pointer to `scientific-context.md`.

Proposed section text:

```markdown
## Documentation namespace

All agent output must be saved as persistent, dated artifacts — never ephemeral
tmp files. Do NOT inline output into this file, CLAUDE.md, or context.md.

| You are writing… | It goes in… |
|------------------|-------------|
| Literature / web research | `docs/_internal/ai-generated/research/YYYY-MM-DD_<topic>.md` |
| Codebase / data exploration | `docs/_internal/ai-generated/explorations/YYYY-MM-DD_<what>.md` |
| Session handoff | `docs/_internal/ai-generated/sessions/YYYY-MM-DD_<slug>.md` |
| Decision trace / ADR / why-not log | `docs/_internal/ai-generated/decisions/<slug>.md` |
| Devops / phase plan | `docs/_internal/ai-generated/plans/YYYY-MM-DD__<phase-name>/` |
| Biological framework (hypothesis, gene sets) | `docs/_internal/scientific-context.md` |
| Public how-to-reproduce, stage runbook | `docs/guides/<stage>.md` (opt-in) |
| Public milestone tracking | `docs/plan/` (opt-in) |

Rules:
- **Never write agent output outside `docs/_internal/`** unless explicitly producing
  a public-facing artifact the human requested.
- **Never grow context.md.** It is a fixed-size pointer. New biology goes in
  `docs/_internal/scientific-context.md`.
- **Never grow this file (AGENTS.md) with project narrative.** Project-specific
  findings go in `docs/_internal/`. This file stays under ~120 lines (see "When
  AGENTS.md gets fat" below).
- **One current handoff per day.** Leave prior dated handoffs in
  `docs/_internal/ai-generated/sessions/` — they are the continuity record, not clutter.

**One-way reference rule:**

`docs/_internal/` CAN reference any public artifact — scripts, result files, config keys,
function names. The reverse is FORBIDDEN: `03_results/*/README.md`, `docs/guides/`, and
every other public-facing file MUST NOT link to or mention `docs/_internal/`. Public-facing
artifacts must be fully comprehensible to a reader who has never seen the kitchen. If you
need to explain a decision in a public README, state the outcome — not the internal doc
that records the reasoning.
```

Plus a "When AGENTS.md gets fat" subsection (see 3.6).

**Why:** The core finding is that agents had *no guidance* on where to park outputs, so they defaulted to the nearest always-open file (`context.md`, `CLAUDE.md`). Declaring the routing table inside the file agents always read closes that gap directly. Putting it in `AGENTS.md` (not a separate doc) guarantees every agent sees it via `@AGENTS.md`.

**How:**
1. Add the `## Documentation namespace` section after "Active role", before "Critical rules".
2. Replace the existing `## Directory structure` block's results subtree to reference the phase-based layout and remove the duplicated flat `objects/tables/figures` listing (defer authoritative tree to `analysis_config.yaml` + the results phase).
3. Add explicit "Never grow context.md / Never grow this file" rules to the existing "Critical rules" list.

### 3.3 Decompose context.md into a pointer + scientific-context.md

**Files:** `templates/context.md.template` (rewrite to pointer), new `templates/docs/_internal/scientific-context.md.template` (absorbs the long-form content).

**What:** `context.md` stays at the project root (agents and tools look for it there) but becomes a ~30-line **pointer document**. The long biological hypothesis, gene sets, contrast families, methodology, and references move into `docs/_internal/scientific-context.md`.

New `context.md.template` (full content):

```markdown
# Context: {{PROJECT_ID}}

**Created:** {{DATE}}

This is a pointer. The full scientific framework lives in
[`docs/_internal/scientific-context.md`](docs/_internal/scientific-context.md).
Keep this file short — when in doubt, edit scientific-context.md, not this file.

## One-line question
[What are you investigating? One sentence.]

## System
- **Species / system:** [e.g. Mus musculus, eWAT]
- **Design:** [e.g. 2 diet × 3 timepoints × 4 donors]

## Where things are
- Biological framework → `docs/_internal/scientific-context.md`
- Research traces → `docs/_internal/ai-generated/research/`
- Codebase / data explorations → `docs/_internal/ai-generated/explorations/`
- Decision records (ADRs) → `docs/_internal/ai-generated/decisions/`
- Latest handoff → `docs/_internal/ai-generated/sessions/` (newest dated file)
- Phase plans → `docs/_internal/ai-generated/plans/`
- Config / thresholds → `02_analysis/config/analysis_config.yaml`

## Status
- **Current stage:** [e.g. 02_programs — cNMF k-selection]
- **Blocked on:** [or "nothing"]
```

New `scientific-context.md.template` absorbs the old `context.md` body: Scientific Question (long form), Datasets table, Hypotheses, Analysis Goals, Broad Direction, Gene sets / signatures, Contrast families, Key References, interpretation Notes. This is where 477–755-line content is *allowed* to live, because it is curated SSoT, not a catch-all.

**Why:** `context.md` is loaded into agent context constantly; a 755-line file is expensive and gets stale. Splitting concerns means the always-loaded file is cheap and orienting, while the deep biology is one click away and edited deliberately. The pointer's "Where things are" block doubles as a navigation index that reinforces the namespace.

**How:**
1. Move the existing `context.md.template` body into `scientific-context.md.template`, expanding with `## Gene sets` and `## Contrast families` headers (observed real-project content that had no home).
2. Replace `context.md.template` with the pointer above.
3. Update any tooling/agents that grep `context.md` for biology to follow the link (search `agents/` and `skills/` for `context.md` references).

### 3.4 Retarget the handoff agent to docs/_internal/sessions/

**Files:** `agents/analysis-base/handoff.md`

**What:** Change the handoff agent's output target from project-root `handoff_YYYYMMDD_HHMMSS.md` + `.handoff_archive/` to `docs/_internal/ai-generated/sessions/YYYY-MM-DD_<slug>.md`, keeping all dated files in the same folder (newest = current; no separate archive dir).

**Why:** The handoff agent is a primary source of root clutter and directly contradicts the new namespace. Dated files in one folder are self-archiving — `ls docs/_internal/ai-generated/sessions/` is the manifest. The entire `ai-generated/` subtree is committed (it is the project's durable continuity record), so handoffs accumulate as the audit trail rather than being periodically archived or deleted.

**How:**
1. Change the MANDATORY filename format to `YYYY-MM-DD_<slug>.md` (keep a time suffix only if same-day collisions are expected: `YYYY-MM-DD_HHMM_<slug>.md`).
2. Change all write paths to `docs/_internal/ai-generated/sessions/`.
3. Delete the Step 3 "archive to `.handoff_archive/`" logic; replace with "leave prior dated handoffs in place — they are the continuity record".
4. Update the "Gather Current Context" step to read `docs/_internal/scientific-context.md` and `context.md` (pointer) instead of `plan.md` at root; read `docs/_internal/ai-generated/sessions/` for most recent prior handoff.
5. Add a pre-write verification step: before producing the session file, scan `03_results/` for any artifact files that lack a corresponding caption entry in their phase `README.md`. List missing entries in the handoff under a "## Uncaptioned artifacts" section so the next session can open with caption writing rather than silently losing provenance.

### 3.5 Golden path lives in a guideline, referenced — not inlined in templates

**Files:** `docs/guidelines/checkpoint_caching.md` (already canonical), `roles/base.yaml` (already activates `scrna-pipeline-conventions` skill), `templates/AGENTS.md.template` (reference only).

**What:** The compute/visualize separation + `load_or_compute()` checkpoint caching ("golden path") stays in `docs/guidelines/checkpoint_caching.md` and the `scrna-pipeline-conventions` skill. `AGENTS.md.template` carries a **one-line pointer**, not the pattern body. We do NOT copy the pattern into the template.

**Why:** In real projects the golden path was copy-pasted into `CLAUDE.md`, which is exactly why those files hit 500+ lines. The pattern is already version-controlled canonically in `docs/guidelines/` and surfaced through the `base` role's skills. Inlining it into a template would re-seed the bloat we are removing. Pointer + role-activated skill = single source of truth, no duplication. (See ADR-3.2.)

**How:**
1. In `AGENTS.md.template`, under "Critical rules", keep the existing rules 1–4 (config-not-hardcoding, normalize-then-visualize, read-only input, cache-expensive) since they are *principles*, and add a one-liner: "Full checkpoint pattern: see the `scrna-pipeline-conventions` skill and `docs/guidelines/checkpoint_caching.md`."
2. Do not add the `load_or_compute()` code body to any template.

### 3.6 AGENTS.md size discipline + split rule

**Files:** `templates/AGENTS.md.template` (add subsection), `agents/analysis-base/doc-curator.md` (add a check).

**What:** Add explicit size guidance so `AGENTS.md` cannot silently become the next monolith.

Subsection text for the template:

```markdown
### When AGENTS.md gets fat

AGENTS.md should stay under ~120 lines. It holds: project one-liner, active role
block, critical rules, namespace routing, directory pointer. Nothing else.

If it is growing, the content belongs elsewhere:
- Agent findings / research notes → `docs/_internal/ai-generated/research/` or `explorations/`
- Stage-by-stage instructions → `docs/guides/<stage>.md` (public) or a decision record
- Biology → `docs/_internal/scientific-context.md`
- Reusable patterns → a skill in SciAgent-toolkit, not here

Split trigger: if a section exceeds ~25 lines or describes a single stage in
detail, move it to docs/ and leave a one-line pointer.
```

**Why:** 634-line `AGENTS.md` was the worst observed monolith. A hard-ish budget plus a "where does this belong instead" table makes the right move obvious. Giving `doc-curator` a check operationalizes enforcement.

**How:**
1. Add the subsection to `AGENTS.md.template` (after the namespace section).
2. Add to `doc-curator.md` a discovery check: "Flag `AGENTS.md` / `CLAUDE.md` > ~150 lines and `context.md` that is not a pointer; recommend decomposition into `docs/_internal/`."
3. Add a second curator check: scan all public-facing files (`03_results/**/README.md`, `docs/guides/**`, `docs/*.md` outside `_internal/`) for any link or path reference containing `_internal/`. Report violations of the one-way reference rule and suggest rewriting them to state the outcome rather than cite the internal doc.

### 3.7 Clarify AGENTS.md vs CLAUDE.md division of labor

**Files:** `templates/CLAUDE.md.template` (keep as `@AGENTS.md`), documented in `docs/architecture.md`.

**What:** Lock the contract: `AGENTS.md` is the **single source of truth** for agent behavior (harness-agnostic); `CLAUDE.md` is a **thin Claude-specific import** (`@AGENTS.md`) plus, if ever needed, Claude-only directives. All project narrative that crept into `CLAUDE.md` (golden paths, schemas) moves to `docs/_internal/` or skills.

**Why:** Real `CLAUDE.md` files reached 216–532 lines by absorbing golden-path patterns and master-table schemas — content that is either harness-agnostic (belongs in `AGENTS.md`/skills/guidelines) or project-specific (belongs in `docs/_internal/`). Keeping `CLAUDE.md` as a one-line import enforces "edit AGENTS.md, not CLAUDE.md," which the template already states but reality violated. (See ADR-3.3.)

**How:**
1. Keep `CLAUDE.md.template` = `@AGENTS.md` (no change to content).
2. Add a comment line below it: `<!-- Claude-specific only. Shared rules live in AGENTS.md. Project narrative lives in docs/_internal/. -->`
3. In `doc-curator` check (3.6), flag any `CLAUDE.md` longer than 2 substantive lines as drift.

## Open ADRs

### ADR-3.1: Does context.md become a pointer, or keep full content?
**Options:**
A — Pointer (~30 lines) referencing `docs/_internal/scientific-context.md`.
B — Keep `context.md` as the full biological doc; no separate file.
C — Delete `context.md` entirely; agents read `scientific-context.md` directly.
**Recommended:** A
**Why:** `context.md` is loaded into agent context frequently; a 755-line always-loaded file is costly and goes stale. A pointer keeps the hot path cheap while the deep, deliberately-edited biology lives one click away. C is cleaner conceptually but breaks the strong convention (tools/agents/users expect `context.md` at root); a 30-line pointer preserves discoverability for near-zero cost.
**Blocking implementation:** no (recommendation is strong; proceed unless overruled)

### ADR-3.2: Golden path in templates or in roles/skills?
**Options:**
A — Inline the `load_or_compute()` pattern body into `AGENTS.md.template`.
B — Keep it in `docs/guidelines/checkpoint_caching.md` + `scrna-pipeline-conventions` skill; template only points.
C — New dedicated skill that every analysis role activates.
**Recommended:** B
**Why:** The pattern is already canonical in `docs/guidelines/` and surfaced via the `base` role's `scrna-pipeline-conventions` skill. Inlining (A) is precisely what bloated real `CLAUDE.md` files to 500+ lines. B = zero duplication, single SoT, already wired. C is viable later if more pipeline conventions accrete, but is redundant with the existing skill today.
**Blocking implementation:** no

### ADR-3.3: AGENTS.md vs CLAUDE.md division of labor — enforce thinness?
**Options:**
A — `CLAUDE.md` stays `@AGENTS.md` only; everything shared in `AGENTS.md`, everything project-specific in `docs/_internal/`.
B — Allow `CLAUDE.md` to carry Claude-specific golden paths/schemas.
**Recommended:** A
**Why:** Every byte of project narrative that landed in `CLAUDE.md` was either harness-agnostic (belongs in `AGENTS.md`/skills) or project-specific (belongs in `docs/_internal/`). There is no legitimate "Claude-only project narrative" category. A makes the boundary crisp and curator-checkable; the template already claims A — we just enforce it.
**Blocking implementation:** no

### ADR-3.4: Is docs/_internal/ committed or gitignored?
**Options:**
A — Committed entirely (full audit trail travels with the repo).
B — Gitignored entirely (treat as local scratch).
C — Split: `ai-generated/` committed, human free-form notes gitignored.
D — Neither committed nor gitignored by SciAgent-toolkit; per-project decision only.
**Recommended:** A for private/analysis repos; D for repos you intend to publish publicly (gitignore `docs/_internal/` at publish time).
**Why:** `docs/_internal/` is the project's reasoning record — that is exactly the content you want to preserve for future-self auditability ("why was this decision made in 2026?"). There is no `scratch/` subdirectory anymore (it was dropped; unclassified ephemera should not be saved), so there is no longer a case for splitting. For public repos: a single `echo 'docs/_internal/' >> .gitignore` at publication time strips the private kitchen from the public view. SciAgent-toolkit should not pre-configure the gitignore since it cannot know the repo's eventual publicity status.
**Blocking implementation:** no

### ADR-3.5: Handoff filename — date-only or date+time?
**Options:**
A — `YYYY-MM-DD_<slug>.md` (date-only, human-friendly, sorts well).
B — `YYYY-MM-DD_HHMM_<slug>.md` (collision-proof for multiple handoffs/day).
**Recommended:** A, with B as documented fallback when same-day collisions occur.
**Why:** Most days have one handoff; date-only reads cleanly and sorts chronologically. The slug carries meaning. Add the time component only when a same-day collision is actually hit, keeping the common case clean.
**Blocking implementation:** no

## Dependencies
- Depends on: results-restructure phase (phase-based `03_results/` layout: `01_qc/`, `02_programs/`, `objects/`) — referenced by `AGENTS.md.template` directory section and `context.md` pointer.
- Depends on: scaffolding phase (`sciagent new project` becomes a full scaffold rather than 3-file copy) — this phase owns only the `docs/` subtree it must copy.
- Enables: doc-curator enforcement (size/pointer checks) becomes meaningful only once the namespace exists.
- Enables: cleaner role-specific instructions, since project narrative now has a home outside role/harness files.

## Breaking Changes
- `context.md` semantics change: it is now a pointer, not the full biological doc. Any agent/skill/tool that reads the full hypothesis from `context.md` must follow the link to `docs/_internal/scientific-context.md`. (Acceptable — single user, no downstream consumers.)
- `handoff` agent output path moves from project root + `.handoff_archive/` to `docs/_internal/sessions/`. Old root `handoff_*.md` files in existing projects are not auto-migrated.
- `AGENTS.md.template` directory tree changes from flat `objects/tables/figures` to phase-based references; templates only — no migration of existing projects.
- New projects gain a `docs/` skeleton they did not have before (additive, non-breaking for new scaffolds).

## Estimated Scope
- `templates/context.md.template`: rewrite, ~735 B → ~30 lines (pointer + updated "Where things are" table).
- `templates/docs/_internal/scientific-context.md.template`: new, ~80–120 lines (absorbs old context body + gene-sets/contrast headers).
- `templates/AGENTS.md.template`: +~50 lines (namespace section with updated routing table, size-discipline subsection, rule additions, directory pointer rewrite).
- `templates/CLAUDE.md.template`: +1 comment line.
- `templates/docs/` skeleton: ~9 new files (`.gitkeep` × 6 for `research/`, `explorations/`, `sessions/`, `decisions/`, `plans/` under `ai-generated/`, plus `_internal/README.md`, top-level `README.md`/`notes.md` stubs).
- `agents/analysis-base/handoff.md`: ~30 lines changed (path updated to `ai-generated/sessions/`, archive logic removed, read-targets updated).
- `agents/analysis-base/doc-curator.md`: +~10 lines (drift checks including `ai-generated/` path validation).
- `docs/architecture.md`: +~15 lines documenting the two-axis namespace + AGENTS/CLAUDE contract.

**Net:** ~9 files touched/created, roughly +260 / −20 lines, plus the directory skeleton. No code-path changes to `lib/sciagent/` in this phase (the `new project` wiring is owned by the scaffolding phase; this phase only supplies the `docs/` template payload it copies).
