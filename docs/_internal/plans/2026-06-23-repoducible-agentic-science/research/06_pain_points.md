# Pain Points Catalog — Recurring Hand-Steering Requirements

**Date:** 2026-06-23  
**Scope:** Requirements analysis derived from owner's recurring hand-steering prompts and concrete example prompts A–E.  
**Method:** Cross-checked against canonical toolkit (`SciAgent-toolkit`), three live repos (`STING-cGAS-GSE329522`, `DC-nexus/DC_mouse_cancer`, `14839-DM-cGAS`), and internal plan artifacts.

---

## Theme 1 — Figure Legibility & Publication/Conference Dual-Context

| # | Req name | Latent durable rule | Where it should live | Current status | Acceptance criteria (no longer needs hand-steering) |
|---|---|---|---|---|---|
| 1.1 | **dual-context-legibility** | Every figure must be simultaneously legible (a) shrunk into a single-column journal layout and (b) projected to the back row of a conference room. This resolves to: base font ≥ 16 pt at configured canvas size; smallest rendered text ≥ 7 pt after journal reduction; no chartjunk; no crowding; thicker lines/larger points than default. | `AGENTS.md` template (project-level rule) + a `figure-standards` skill or system-prompt injected into the base role | **Partially present.** `DC_mouse_cancer/AGENTS.md §Figure standards` contains this rule verbatim with per-config-key guidance (`save_figure()`, `set_paper_style()`, `POINT_SIZE`, `CONFIG$figures.base_size`). `STING-cGAS-GSE329522/AGENTS.md` does NOT contain it. `14839-DM-cGAS/AGENTS.md` does NOT contain it directly (it lives inside the per-plan `CONVENTIONS_results-and-figures.md §3`). `SciAgent-toolkit/docs/guidelines/visualization.md` has `base_size=12` (too small) and no dual-context mandate. No base role skill encodes it. | Every new project AGENTS.md instantiated from the template already contains the dual-context rule, the `base_size ≥ 16` floor, and the `save_figure()` centralization mandate. A `doc-curator` compliance check (C4 rule) flags any viz script that has an inline `theme(text = element_size(...))` or `ggsave(width = <literal>...)`. |
| 1.2 | **centralized-figure-theme** | One theme entry point (`project_theme()` / `save_figure()` / equivalent) per project, reading sizes from config, never hardcoded per-script. All viz scripts call it; none duplicate theme blocks. | `AGENTS.md` template (Critical rule) + a `scrna-pipeline-conventions`-equivalent skill for bulk RNA-seq projects | **Partially present.** `STING-cGAS-GSE329522/AGENTS.md` mandates `save_figure()` and `set_paper_style()`. `DC_mouse_cancer/AGENTS.md §Figure standards` mandates `save_figure()` and `set_paper_style()`. `14839-DM-cGAS` plan `CONVENTIONS §3.2` mandates `project_theme()`/`save_figure()`/`style_running_sum()`. The toolkit `docs/guidelines/visualization.md` ships a `theme_publication(base_size=12)` — correct pattern, wrong default, not linked to a centralization mandate. No base role skill enforces the no-duplication rule. | AGENTS.md template contains explicit "centralize styling" rule with named entry points. `doc-curator` C4 compliance check (`grep -rn 'theme(text' / 'ggsave(width='` in viz scripts) is in the agent's default checklist. |
| 1.3 | **figure-series-overlay-ready** | When a series of figures will be compared side-by-side (across contrasts, databases, timepoints), panel geometry and axis scales must be invariant to per-figure content. Variable-width furniture (legends, long labels) must not float the panel edge. Shared scales come from config, not per-figure auto-scaling. | `AGENTS.md` template (Figure standards section) + per-plan `CONVENTIONS` document (already exists in `14839-DM-cGAS`) | **Ad-hoc / per-project.** Exists only in `14839-DM-cGAS` plan `CONVENTIONS §3.3` as a reusable pattern with concrete helpers (`style_running_sum()`, `running_sum_ylim`). Not in any other repo's AGENTS.md. Not in toolkit docs. | Template AGENTS.md has a rule: "When emitting a figure series (same analysis type, multiple contrasts/DBs/timepoints), panel geometry and axis limits must be fixed from config — no per-figure auto-scaling. Verify on `panel` cell bounding box, not page size." |
| 1.4 | **figure-critique-checklist** | A canonicalized figure quality checklist covering the recurring critique themes from prompt D: no axis-text truncation; no ambiguous glyph semantics; text never squished; column clustering for pattern recognition; cap to top-N; prefer residualized channel; README explains how to read the figure. | A `figure-audit` agent or a `code-reviewer`-extension checklist; or a `/figure-audit` slash command | **Absent.** `toolkit/agents/architect/graphic.md` does Tufte-principled graphic review but is for interactive dashboards (not static figures) and is in the `architect` role, not the `base` bioinformatics role. No checklist agent exists for static publication figures. The `doc-curator` C3 checks caption completeness, not figure rendering quality. | A `figure-audit` agent (base role) runs on a `03_results/<stage>/figures/` directory and checks the D-theme list: no truncated axis text, no ambiguous glyphs, column-clustering present when appropriate, top-N capping applied, README has reading instructions. Can run as a subagent pass after figure generation. |

---

## Theme 2 — Results Directory Routing: `03_results/{stage}/{figures,tables,README.md}`

| # | Req name | Latent durable rule | Where it should live | Current status | Acceptance criteria |
|---|---|---|---|---|---|
| 2.1 | **results-layout-canonical** | Every artifact goes to `03_results/<stage>/{tables,figures}/` (with optional `by_contrast/<contrast>/` and `_overview/` subdivisions). No artifact lands in an ad-hoc location or at the project root. Stage IDs are declared once in `analysis_config.yaml:stages`. | `AGENTS.md` template (Directory structure section, Critical rules) | **Present but inconsistent across repos.** `DC_mouse_cancer/AGENTS.md §Directory structure` and `14839-DM-cGAS/AGENTS.md §Directory structure` both specify the layout. `STING-cGAS-GSE329522/AGENTS.md` specifies it in the Normalize-then-visualize section. The template in `SciAgent-toolkit/templates/` (checked via `GEMINI.md.template`) does not contain a per-project directory structure skeleton. Stages declared in config is mandated but the `stages:` key itself is not scaffolded in the new-project template. | `sciagent new project` scaffold writes `AGENTS.md` with the full `03_results/` layout table and a starter `analysis_config.yaml` that includes a `stages:` key. Any agent writing to `03_results/` reads `stages:` first; writing to an undeclared stage ID is a hard error documented in AGENTS.md. |
| 2.2 | **compute-viz-strict-split** | Compute scripts (`NN_<name>.R` / `.py`) contain no `ggplot`/`ggsave`/`matplotlib` calls. Viz scripts (`NN_<name>_viz.R`) contain no statistical computation. They are always separate files with separate stage writes. | `AGENTS.md` template (Critical rules) + phase-plan template | **Present in two repos.** `STING-cGAS-GSE329522/AGENTS.md §Normalize then visualize` states this explicitly. `DC_mouse_cancer/AGENTS.md` Critical rule 2 states it. `14839-DM-cGAS/AGENTS.md` Critical rule 2 states it. The 2026-06-21-te/00_INDEX.md calls it "STRICT." But it is not in the project scaffold template, so new projects start without it. | AGENTS.md template Critical rule 2 is the compute/viz split. Phase-plan template enforces compute phases have "No figures" in their scope. A `code-reviewer` check (`grep ggsave NN_<name>.R` / `grep lmFit NN_<name>_viz.R`) catches violations. |
| 2.3 | **stage-id-discipline** | Stage IDs use two-digit numeric prefix + lowercase snake_case slug, declared in `analysis_config.yaml:stages`. New stages are appended there before any script writes to them. Pipelines operating in the same repo on disjoint biological questions (e.g., gene pipeline vs TE pipeline) use disjoint numeric prefixes. | `AGENTS.md` template + per-plan `00_INDEX.md` format | **Present only in per-plan artifacts.** `14839-DM-cGAS` TE plan explicitly uses `20_te_qc`–`24_te_viz` to be disjoint from gene `00`–`07`. But this decomposition principle is not in AGENTS.md or any AGENTS.md template. New plans have no enforced stage-ID naming convention. | AGENTS.md template specifies the `NN_<slug>` format and disjoint-prefix rule for parallel pipelines. Phase plan `00_INDEX.md` template requires a `stage_id` column; the INDEX is the source of truth checked before any script writes. |

---

## Theme 3 — README Caption Obligation After Figure Creation/Edit/Delete

| # | Req name | Latent durable rule | Where it should live | Current status | Acceptance criteria |
|---|---|---|---|---|---|
| 3.1 | **caption-after-figure** | Any agent that creates, modifies, or deletes a figure must update the sibling `README.md` in that stage directory before considering the task complete. Caption format: `## <relative/path>/<file>`, one-sentence scientific finding, `Script | Function | Config | Input` table. Never cite `docs/_internal/`. | `AGENTS.md` template (Artifact captions section, Critical rules) + `captions` agent description | **Partially present.** All three live repos have "Artifact captions" in AGENTS.md. `doc-curator` agent checks C3 (uncaptioned artifacts). `captions` agent exists and writes READMEs. BUT: (a) the `captions` agent is not automatically triggered by figure creation — it is dispatched by the human. (b) The agent's system prompt describes the README format correctly but does not reference the path-qualified `## <relative/path>/<file>` convention mandated in `14839-DM-cGAS` plans. (c) No hook or agent rule says "before returning from a task, verify adjacent README is updated." | The `captions` agent description states: "ALWAYS run as a mandatory cleanup subagent pass after any figure is created, modified, or deleted." The base role `AGENTS.md` contains: "A task is not complete until the sibling `README.md` has a caption entry for every file touched under `03_results/`." A post-task hook or the `doc-curator` C3 check catches violations. |
| 3.2 | **readme-reading-instructions** | Every `_overview` figure's README entry must explain how to read the figure: what each axis/glyph/channel means, what the sign convention is, what the claim tier is. Per-contrast figures must have a sidecar `README.md` with the direction legend. | `AGENTS.md` template (Artifact captions section) + `captions` agent | **Partially present.** The `CONVENTIONS_results-and-figures.md §2` in `14839-DM-cGAS` mandates per-contrast `README.md` direction legends via `write_contrast_readme()`. This is project-specific; no base toolkit rule. The `captions` agent's system prompt does not mention direction semantics or reading instructions. The figure critique theme from prompt D explicitly lists "READMEs must explain how to read the figure" as a recurring gap. | `captions` agent system prompt includes: "For each figure, after the scientific finding, add a 'How to read' section: what the sign convention is, what each glyph means, what claim tier the figure supports." |

---

## Theme 4 — Planning Artifacts: Location, Decomposition Rules, Model Assignment

| # | Req name | Latent durable rule | Where it should live | Current status | Acceptance criteria |
|---|---|---|---|---|---|
| 4.1 | **plan-artifact-location** | Planning artifacts go to `docs/_internal/plans/{date}-{slug}/` with an `00_INDEX.md` phase table plus numbered phase briefs `NN_<slug>.md`. Research notes from planning go to `docs/_internal/research/{date}-{slug}/`. Reasoning traces from plan reviews go to `docs/_internal/reasoning/`. None of these ever go to `docs/` (public) or `03_results/`. | `AGENTS.md` template (Documentation namespace table) | **Present in template AGENTS.md.** Both `DC_mouse_cancer/AGENTS.md §Documentation namespace` and `14839-DM-cGAS/AGENTS.md §Documentation namespace` have the routing table. The `00_INDEX.md` format is well-established by practice. BUT: the routing table does not mention `plans/` explicitly — it lists `sessions/`, `research/`, `reasoning/`, `scientific-context.md`, and `docs/plan/` (public), but not the `docs/_internal/plans/` internal planning directory. | AGENTS.md template documentation-namespace table includes: `plans/ → docs/_internal/plans/{date}-{slug}/` with format note "00_INDEX.md + NN_<slug>.md briefs." A planner that attempts to write plan files elsewhere fails the `doc-curator` one-way reference check. |
| 4.2 | **phase-bite-size-rule** | Each phase must fit inside ~30–40% of Sonnet's 200k context when read together with AGENTS.md and relevant skills. Rule of thumb: one phase = one downstream script; ≤3–5 files touched; one clearly testable output. Phases that exceed this must be split. | Phase-plan template (`00_INDEX.md`) + planning agent instructions | **Ad-hoc / per-plan.** `14839-DM-cGAS` 2026-06-21-te/00_INDEX.md states "one phase == one downstream script == one Sonnet brief" and "Decompose into bite-sized phases fitting ~30–40% of Sonnet's 200k context." This rule was composed by the owner in prompt B; it is not in any canonical AGENTS.md or planning tool definition. `SciAgent-toolkit/commands/architect/plan.md` specifies "≤3–5 files per phase" but does not mention context-budget framing or the Sonnet constraint. | A planning agent's system prompt or a phase-plan template explicitly states the Sonnet-context-budget rule. The `00_INDEX.md` template has a `context_budget` column or note. A reviewer step rejects any phase brief that references more than one script as the primary deliverable. |
| 4.3 | **model-role-assignment** | Planning and decomposition: done by larger model (Opus or equivalent). Implementation: one Sonnet per phase. Review every 2–3 phases: Opus. The planner never implements; the implementer never re-plans. Research exploration ("brainstorm / explore" waves): Opus fan-out. | Role definitions in the base planning workflow; could be an AGENTS.md Critical rule or a `/plan` command parameter | **Ad-hoc / absent from canonical docs.** The architect role agents (`architect.md`, `synth.md`, `meta-architect.md`, `divergent.md`) are already `model: opus`. The `captions` / `code-reviewer` / `insight-explorer` agents are `model: sonnet`. But the RULE "plan with Opus, implement with Sonnet, review with Opus every 2–3 phases" is not written anywhere in AGENTS.md or in any command definition. The `implement` command (`/implement`) is model-agnostic. | A `plan-and-implement` orchestration command (or the AGENTS.md "Planning" section) states: planner = Opus, implementer = Sonnet per phase, reviewer = Opus every N phases. The `00_INDEX.md` template has `model` column per phase. |
| 4.4 | **opus-review-cadence** | After every 2–3 implementation phases: an Opus reviewer checks plan adherence, consistency, code cleanliness, namespace separation, accidental complexity, drift/rot. Reviewer also verifies the implementation is runnable and produces the artifacts it claims. Lengthy phases may run in a monitored tmux session while later non-blocked phases proceed. | A `/review-phases` command or an orchestration command that enforces this cadence | **Absent as a shipped command.** `SciAgent-toolkit/commands/architect/review.md` exists but is scoped to the software-architect workflow (design docs + code, not scientific pipeline phases). No command encodes the "every 2–3 phases" cadence for science pipelines. The `code-reviewer` agent exists but is invoked manually. tmux pattern is used in `DC_Dictionary/CLAUDE.md` and `DC_hum_verse` but is not encoded in any agent. | A `/pipeline-review` command (or parameterized phase in the `/plan-pipeline` orchestration recipe) dispatches an Opus reviewer after phase N. The reviewer's checklist covers: plan adherence, runnable artifacts, namespace, drift. tmux-backgrounding of long phases is encoded in the reviewer's system prompt. |

---

## Theme 5 — Reproducibility: Persistent Reasoning Traces & Non-Ephemeral Scripts

| # | Req name | Latent durable rule | Where it should live | Current status | Acceptance criteria |
|---|---|---|---|---|---|
| 5.1 | **no-ephemeral-scripts** | All scripts written during analysis must be committed to the repo at their canonical path (`02_analysis/scripts/NN_<slug>.R`). No throwaway snippets run from `/tmp/`, the agent's scratch space, or a notebook cell that is not saved. Every invocation that produces a `03_results/` artifact must be reproducible from a committed script. | `AGENTS.md` template (Critical rules) | **Present by convention but not explicitly stated as a rule against ephemeral scripts.** The compute/viz split (Theme 2) implies this, but no AGENTS.md rule says "do not run throwaway scripts." The `_scratch/` directory in `03_results/` exists precisely as the ephemeral escape valve — but agents can and do drop files there without tracing them. `doc-curator` C3 does not catch scripts in `/tmp/`. | AGENTS.md Critical rule: "Never run a script that produces a `03_results/` artifact from a path outside `02_analysis/scripts/`. No `/tmp/` scripts, no in-session snippets. Every artifact must be reproducible by re-running the committed script." `doc-curator` checks that all artifacts under `03_results/` (excluding `_scratch/`) have a caption with a `Script` field pointing to a committed path. |
| 5.2 | **persistent-reasoning-traces** | Agent decision traces — why a particular method was chosen, why a phase boundary was drawn here, why a parameter was set — must be persisted in `docs/_internal/reasoning/YYYY-MM-DD_NN_topic.md`. These are not chat messages; they are files in the repo that survive session handoff. A downstream agent (or the owner) can read them to understand what was decided and why. | `AGENTS.md` template (Documentation namespace) + base agent descriptions | **Partially present.** `14839-DM-cGAS/docs/_internal/reasoning/` is actively populated (14+ files covering plan architecture, Opus reviews, cross-agent synthesis). The `README.md` in that directory defines the naming convention. But: (a) not all three live repos actively populate `reasoning/`; (b) the base AGENTS.md template for new projects does not instruct agents to write reasoning traces; (c) the `handoff` agent writes sessions but not general reasoning traces; (d) no agent has an explicit rule "before completing a non-trivial decision, write a reasoning trace." | Base AGENTS.md template has a rule: "For any non-trivial decision (method selection, phase boundary, parameter choice, unexpected finding), write a reasoning trace to `docs/_internal/reasoning/YYYY-MM-DD_NN_<topic>.md` before proceeding. A decision with no reasoning trace is non-reproducible." The `handoff` agent description includes writing a reasoning summary to `reasoning/` as part of its protocol. |
| 5.3 | **session-handoff-continuity** | Session handoff notes must record enough state that a new agent session can reproduce the same decisions: what scripts were run, what artifacts were produced, what decisions were deferred, what the next step is. References to `03_results/` paths and committed scripts, never ephemeral context. | `handoff` agent + AGENTS.md | **Present.** The `handoff` agent exists. `docs/_internal/sessions/` is defined. `AGENTS.md` in all live repos states "Session handoffs are written to `docs/_internal/sessions/YYYY-MM-DD_<slug>.md`." Gap: the handoff agent's system prompt (`handoff.md` — not read in this analysis) may not mandate reference to committed script paths as the reproducibility anchor. | The `handoff` agent description mandates: "Every handoff must list all scripts run since the last handoff (with committed paths), all artifacts produced (with `03_results/` paths), all open decisions with their `reasoning/` trace file." |
| 5.4 | **research-note-persistence** | Exploration outputs (web research, codebase reads, biology brainstorm results) must be saved as dated files in `docs/_internal/research/{date}-{slug}/` before being fed to planning agents. These are the injection corpus for planners; ephemeral in-session summaries are not. | AGENTS.md + research/exploration agent descriptions | **Present as convention in active use.** `14839-DM-cGAS/docs/_internal/research/` is heavily populated (22+ files across 4 dated slugs). The `bio-interpreter` agent writes to `web_notes.md`. The pattern from prompt B ("save interim explore artifacts to `docs/_internal/research/{date}-te-codebases` to serve as injection to planner agents") is followed in practice. Gap: no base agent description explicitly states "save to `docs/_internal/research/` before returning." | Base agent descriptions (`insight-explorer`, `bio-interpreter`) include: "Save research findings to `docs/_internal/research/{date}-{slug}/NN_<topic>.md` before returning. Do NOT summarize in chat only — a chat summary is ephemeral." |

---

## Reusable Orchestration Recipes (from Prompts A, B, C, E)

These are parameterized workflows the toolkit should ship as commands so the owner stops re-typing them. Each recipe is currently composed ad-hoc in each project.

### Recipe R1 — `/pipeline-plan` (from Prompt A, B)

**Purpose:** Plan a multi-phase analysis pipeline and orchestrate its implementation with model-appropriate agents.

**Parameters:**

| Parameter | Type | Description |
|---|---|---|
| `slug` | string | Plan slug, used as `docs/_internal/plans/{date}-{slug}/` |
| `scope_doc` | path | Research notes or context doc to feed the planner (e.g., `docs/_internal/research/{date}-{slug}/`) |
| `n_planners` | int (default 1) | Number of parallel Opus planner agents for decomposition |
| `context_budget_pct` | int (default 35) | Target % of Sonnet 200k context per phase brief |
| `review_cadence` | int (default 3) | Opus review after every N phases |
| `idempotency_contract` | bool (default true) | Require explicit idempotency notes per phase |

**Workflow:**

1. Opus planner reads `scope_doc` + existing AGENTS.md + relevant skills → writes `docs/_internal/plans/{date}-{slug}/00_INDEX.md` with phase table (seq, stage_id, slug, title, depends_on).
2. Opus planner decomposes into N phase briefs (`NN_<slug>.md`), each sized to fit `context_budget_pct` of Sonnet context (one script, ≤3–5 files, one testable output).
3. For each phase in dependency order: launch one Sonnet implementer. After every `review_cadence` phases: launch Opus reviewer (checks plan adherence, runnable artifacts, namespace, drift). If a phase is long-running and does not block later phases: reviewer starts it in a backgrounded tmux session.
4. All reasoning traces from review written to `docs/_internal/reasoning/`.

**Model assignments:** planner = Opus, implementer = Sonnet, reviewer = Opus.

---

### Recipe R2 — `/explore-and-plan` (from Prompt B, E)

**Purpose:** Multi-wave research → planning → implementation for a new analysis domain or dataset.

**Parameters:**

| Parameter | Type | Description |
|---|---|---|
| `slug` | string | Research slug, used as `docs/_internal/research/{date}-{slug}/` |
| `question` | string | The scientific question or analysis domain |
| `n_explorers` | int (default 3) | Parallel Opus exploration agents (web research, codebase reads, biology) |
| `n_planners` | int (default 1) | Opus agents for plan decomposition (usually 1, unless parallel feasibility required) |
| `idempotency_peer` | path | Path to a sibling pipeline whose idempotency must be preserved (e.g., gene pipeline path) |

**Workflow:**

1. **Wave 1 — Research fan-out:** Launch `n_explorers` Opus agents in parallel, each assigned a research lane (e.g., codebase architecture, biological literature, reference implementations). Each saves interim notes to `docs/_internal/research/{date}-{slug}/NN_<lane>.md`.
2. **Wave 2 — Synthesis:** One Opus synthesizer reads all Wave-1 notes → writes `docs/_internal/research/{date}-{slug}/_SYNTHESIS.md` with verified vs inferred claims tagged.
3. **Wave 3 — Plan decomposition:** Opus planner reads synthesis + AGENTS.md → writes `docs/_internal/plans/{date}-{slug}/00_INDEX.md` with idempotency contract (disjoint stage IDs, disjoint checkpoint names, disjoint master table namespaces). If `idempotency_peer` provided, planner explicitly verifies no overlap.
4. **Execution:** Use Recipe R1.

**Model assignments:** explorers = Opus, synthesizer = Opus, planner = Opus, implementers = Sonnet.

---

### Recipe R3 — `/add-figure-variant` (from Prompt C)

**Purpose:** Add a new figure family with its own namespace, fan-out implementation, and post-implementation review.

**Parameters:**

| Parameter | Type | Description |
|---|---|---|
| `slug` | string | Figure family slug (e.g., `te-genotype-panels`) |
| `stage_id` | string | Existing or new stage ID to write figures under |
| `namespace_token` | string | Grep-isolable token for all new identifiers (e.g., `gp`) |
| `n_implementers` | int (default 1) | Parallel Sonnet implementers (split compute/viz into separate agents) |
| `research_dir` | path | Where to save pre-implementation evidence notes |

**Workflow:**

1. **Evidence phase:** Launch Sonnet subagents (or one Opus) to gather evidence (existing stage layout, config keys, upstream checkpoints, reference implementations). Save to `docs/_internal/research/{date}-{slug}/`.
2. **Design phase:** Opus planner writes the plan INDEX with namespace declaration, idempotency contract, acceptance gate, and per-phase `script → stage_id → depends_on` table.
3. **Implementation fan-out:** If compute and viz are separate phases (strict compute→viz ordering), run compute Sonnet first, verify checkpoints on disk, then launch viz Sonnet.
4. **Review:** Opus reviewer checks: namespace isolation (no `grep <namespace_token>` hits in unrelated scripts), plan adherence, rendered figure content (open PDFs, confirm non-empty panels), caption completeness (README updated, no `_internal/` references).
5. **Mandatory cleanup:** `captions` agent pass on the affected stage directory.

**Model assignments:** evidence = Sonnet or Opus, planner = Opus, implementers = Sonnet, reviewer = Opus, captions = Sonnet.

---

### Recipe R4 — `/interpret-storm` (from Prompt E)

**Purpose:** Multi-wave biology interpretation and information-graphic design for a completed analysis.

**Parameters:**

| Parameter | Type | Description |
|---|---|---|
| `dataset_path` | path | Path to preprocessed dataset or results directory |
| `slug` | string | Interpretation slug |
| `n_interpreters` | int (default 3) | Parallel Opus interpretation agents (Wave 2) |
| `n_graphic_designers` | int (default 2) | Agents designing information-graphic figures or interactive viewers (Wave 3) |

**Workflow:**

1. **Wave 1 — Web research:** Parallel Sonnet/Opus agents gather biological literature context. Save to `docs/_internal/research/{date}-{slug}/web-*.md`.
2. **Wave 2 — Interpretation storm:** `n_interpreters` Opus agents independently interpret the dataset, each from a different angle (e.g., mechanism, pathway, clinical relevance). Each saves to `docs/_internal/research/{date}-{slug}/reason-*.md`. One Opus synthesizer writes `_SYNTHESIS.md` with verified vs inferred claims.
3. **Wave 3 — Information-graphic design:** `n_graphic_designers` agents each design figures that expose the interpretation layer (what the discovery is and how it was made). Each proposed figure has a namespace token, a data contract (what table it reads), and a claim-tier annotation. All design notes go to `docs/_internal/reasoning/{date}-<slug>-figures.md`.
4. **Implementation:** Use Recipe R3.
5. **All reasoning traces and interim notes saved to repo** before any implementation begins. No ephemeral summaries.

**Model assignments:** web research = Sonnet, interpreters = Opus, synthesizer = Opus, graphic designers = Opus or Sonnet, implementers = Sonnet.

---

## Gap Summary by Artifact Type

| Artifact | Gaps to close |
|---|---|
| `AGENTS.md` template | Missing: dual-context legibility rule (1.1), centralized-theme mandate (1.2), series-overlay-ready rule (1.3), no-ephemeral-scripts rule (5.1), reasoning-trace obligation (5.2), plans routing in namespace table (4.1) |
| Base role (project `AGENTS.md`) | Missing: figure-audit pass as mandatory post-figure step (3.1), reading-instructions in captions (3.2) |
| `captions` agent | Missing: path-qualified heading format, direction/reading-instructions section, "mandatory after figure creation" framing (3.1, 3.2) |
| `doc-curator` agent | Missing: C4 check (inline theme/ggsave duplication in viz scripts); C3 check could be extended to verify `Script` field points to a committed path (5.1) |
| `handoff` agent | Missing: reasoning-trace summary obligation; committed-script-path anchor in handoff (5.3) |
| `insight-explorer`, `bio-interpreter` agents | Missing: explicit "save to `docs/_internal/research/` before returning" rule (5.4) |
| New commands (absent) | `/pipeline-plan` (R1), `/explore-and-plan` (R2), `/add-figure-variant` (R3), `/interpret-storm` (R4) — none shipped as commands; all composed ad-hoc each session |
| Toolkit `docs/guidelines/visualization.md` | Has `base_size=12` (too small); no dual-context mandate; no centralization rule; no series-overlay pattern |
| `figure-audit` agent | Entirely absent; needed for D-theme critique checklist (1.4) |
