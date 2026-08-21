# Architecture — Reproducible Agentic Science: centralizing the owner's five standing conventions

**Date:** 2026-06-23 · **Author:** Opus (lead architect, synthesis of 8 research traces) · **Status:** design, for approval before `00_INDEX.md` implementation
**Scope (IN):** SciAgent-toolkit (`/data1/users/antonz/pipeline/scbio-docker/toolkits/SciAgent-toolkit`). A toolkit-owned mechanism so that a provider-agnostic agent dropped into an analysis repo follows the owner's five conventions without hand-steering.
**Scope (OUT):** Per-project content (no edits to STING/DC/14839 except an optional backfill phase). No new analysis methods. No replacement of the existing architect pipeline. No move off the bash/no-deps core.

---

## 1. Problem framing & first principles

The owner re-types the same five corrections into every repo: figure legibility, results placement, README adjacency, planning decomposition, and reproducibility/no-ephemeral-scripts. The drift analysis (`05`) is decisive about *why*: of the five, only README-adjacency is centralized (and even that excludes STING, which was never scaffolded from the template). The other four exist as **three-to-N incompatible ad-hoc copies** scattered across per-project `AGENTS.md` bodies, a never-loaded `docs/guidelines/` island, and — worst — buried inside per-project planning docs (`14839/.../CONVENTIONS_results-and-figures.md`). The planning-decomposition recipe "cannot drift because it was never recorded — it can only recur as ad-hoc steering" (`05 §4.4`).

The 2026-06-04 mind-palace plan already named the core insight and the fix, but stopped at proposal: **the environment is the interface.** An agent is tidy in DC not because it is told each turn but because the repo it wakes up in encodes the disposition. The corollary, sharpened by this synthesis:

> **The scaffold IS the interface. A convention the owner re-types is a convention that is not in the scaffold, not single-sourced, or not enforced.**

The toolkit already proves the right mechanism for centralized-and-propagating text: the `SCIAGENT:ROLES` managed block — single-sourced in the toolkit, stamped into every repo by `sciagent activate`, SHA1 drift-detected, re-rendered idempotently. The roles block content cannot drift because no one hand-edits it and `activate` re-renders it. We extend exactly this proven mechanism to craft.

### Design tenets (held throughout)

1. **One source of truth per convention.** Each rule has exactly one canonical home in the toolkit. Everything else references it. We delete the duplicate/contradictory homes (the `scrna-pipeline-conventions` flat layout; the `base_size=12` dead-end guideline).
2. **Enforcement by deterministic check over LLM goodwill.** A rule the agent *should* follow is weaker than a `validate` check that fails closed and a hook that blocks. "Hooks cannot hallucinate" (`08 Area 2`). Where deterministic enforcement is feasible and not chilling to exploration, prefer it.
3. **Progressive disclosure.** Always-on text (CRAFT block) stays terse (≤25 lines, respecting the ~150-line AGENTS.md budget and the "lost in the middle" 43% degradation finding in `08`). Depth lives in skills/helpers loaded on demand.
4. **Provider-agnostic via AGENTS.md.** AGENTS.md is the Linux-Foundation cross-tool standard (28+ tools, `08 Area 1`). The CRAFT block lands in AGENTS.md, not CLAUDE.md, so Codex/Cursor/Gemini/Pi all see it. Hooks are a Claude-Code enforcement *layer on top*, never the only line of defence — the rules degrade gracefully to "stated convention" on harnesses without hooks.
5. **Templatize the gold standard.** 14839 already solved all five. We promote its patterns (`project_theme`/`save_figure`/`save_overview`/`write_caption`/`append_master_table`/`style_running_sum`, the INDEX+phase plan format, research→synthesis→plan flow) into reusable toolkit artifacts rather than re-inventing.
6. **Do not duplicate existing machinery.** Per `03`, the toolkit already ships captions/doc-curator/handoff, a competent architect plan/implement/verify pipeline, the correct stage-based scaffold layout, and `graphic` (a Tufte reviewer). We *wire, extend, and reconcile* these — we do not rebuild them.

---

## 2. The three-layer enforcement model

Every convention is placed into exactly one of three mechanisms, by its nature:

- **(a) CENTRALIZED RULE** — always-on, terse, single-sourced text. Vehicle: a new toolkit-owned managed block `<!-- BEGIN SCIAGENT:CRAFT vN hash=... -->` stamped into AGENTS.md by `activate`/`inject`, alongside `SCIAGENT:ROLES`. Use for *disposition and standing rules* an agent must hold in context at all times (the mind-palace P1 recommendation). Mechanically this requires parameterizing `block.sh` (today its markers are the literal `SCIAGENT:ROLES`) to support a second block id; see §5.
- **(b) CAPABILITY** — a thing the agent reaches for: a skill (load-on-demand procedure), a helper library (code it calls), an agent (delegated worker), or a slash command (user-triggered workflow). Use for *procedures and reusable code* that should not occupy always-on context.
- **(c) GUARDRAIL** — a deterministic check that fails closed: a `sciagent validate` check (CI/pre-commit, harness-agnostic) and/or a Claude-Code hook (`PreToolUse` block via exit 2; `Stop`/`SubagentStop` for end-of-task sweeps). Use for *catching violations the LLM will make under pressure*.

**Placement principle:** a rule is *stated* once as (a), made *doable* by (b), and made *unskippable* by (c). The strongest conventions get all three layers; exploratory-friction-sensitive ones (see §7) get (a)+(b) only, with (c) as a warn-not-block.

Mapping each pain point to the model (justification per cell; explicit non-duplication noted):

| Pain point | (a) Centralized rule | (b) Capability | (c) Guardrail |
|---|---|---|---|
| **A. Figure legibility** | CRAFT block: dual-context floors (base ≥16pt screen / 5–7pt print floor), "centralize styling via the project theme entry point", "no inline theme/ggsave" | NEW `figure-style` skill (R+Py contract) + symlinked helper library `figure_helpers.{R,py}`; EXTEND `graphic` into base via a `figure-audit` capability | NEW `validate` check `figure-style` (grep inline `theme(`/`ggsave(width=`/hex in viz scripts; base_size floor in config); `PreToolUse` warn-hook on figure save |
| **B. Results placement** | CRAFT block: "every artifact under `03_results/<stage>/{figures,tables}/` with `by_contrast/`+`_overview/`; source table adjacent to its figure; compute never plots, viz never computes" | EXTEND helper library (`save_figure`/`save_overview`/`contrast_path`); RECONCILE `scrna-pipeline-conventions` to stage layout | NEW `validate` check `results-layout` (stage-id in config; figure has sibling table; no artifacts at results root) — already partly present as `_validate_docs_layout` |
| **C. README adjacency** | CRAFT block one line: "a task is not done until the sibling README captions every touched `03_results/` file (incl. how-to-read)" | EXTEND `captions` agent (path-qualified headings, how-to-read section, "mandatory after figure change"); `save_overview` writes the caption atomically | EXTEND `doc-curator` C3 + NEW `Stop`/`SubagentStop` hook: uncaptioned artifact under `03_results/` ⇒ flag; `validate` check |
| **D. Planning decomposition** | CRAFT block: "plans → `docs/_internal/plans/{date-slug}/`; Opus decomposes, Sonnet implements one phase, Opus reviews every N" + namespace-table row for `plans/` | NEW commands `/pipeline-plan`, `/explore-and-plan`, `/add-figure-variant`, `/interpret-storm`; NEW plan/research templates (INDEX + phase brief) | The commands' own runnable+artifacts acceptance gate (an Opus reviewer step that *runs the script and checks artifacts exist*), not just anchor-grep |
| **E. Reproducibility** | CRAFT block: "no ephemeral scripts — every `03_results/` artifact reproducible from a committed `02_analysis/scripts/NN_*`; log non-trivial decisions to `reasoning/` before proceeding" | NEW `reasoning-trace` skill + EXTEND `bio-interpreter`/`insight-explorer`/`handoff` to persist to `research/`+`reasoning/` before returning | NEW `PreToolUse` hook blocking writes to `/tmp`/scratch that produce results; `validate` check that captioned artifacts cite a committed `Script:` path |

**Non-duplication, explicit (`03`):** we do NOT add a new captioning agent (extend `captions`); do NOT add a new planning primitive (compose `/pipeline-plan` on top of `/plan`/`/implement`/`/verify`); do NOT add a new graphic reviewer for dashboards (`graphic` stays; we add a *static-figure* audit capability into the base role); do NOT re-scaffold the results layout (the template is already correct — we reconcile the one skill that contradicts it and add the missing `by_contrast/`+`_overview/` knowledge).

---

## 3. Per-pain-point design

### A — Figure legibility

**Centralized (CRAFT block, terse):** dual-context mandate ("figures are judged shrunk in a journal column AND projected to the back of a conference room — bigger, fewer, bolder"); the floors (screen base ≥16pt; smallest rendered text ≥5–7pt at print size; line ≥0.8 / 2pt projected); "style only through the project theme entry point — never inline `theme()`/`ggsave(width=…)`/hex"; top-N capping; no truncated axis labels; prefer the residualized channel. The numeric floors live in config; the block states the *rule*, config carries the *values*.

**Capability — the unified cross-language figure contract.** This is the central technical decision. 14839 has R helpers (`project_theme()`, `save_figure()`, `save_overview()`, `style_running_sum()`); DC has Python (`set_paper_style()`, `save_figure()`, `POINT_SIZE`, `purge_figures()`). They are the *same contract* in two languages. We promote a single named contract — the **figure-style contract** — implemented in both:

| Contract function | R (14839 origin) | Python (DC origin) | Reads from config |
|---|---|---|---|
| `project_theme()` / `set_paper_style()` | `project_theme(base_size, legend)` | `set_paper_style()` | `figures.{base_size,title_size,…}` |
| `save_figure(plot, stage, name, …)` | cairo_pdf + png | png@300 + pdf@600 | `figures.{width,height,dpi}` |
| `save_overview(plot, stage, name, table_df, finding, script, fn, config_kv, input)` | atomic fig+table+caption | same | layout + caption format |
| `contrast_path()` / `overview_path()` | mkdir by_contrast/_overview | same | `paths` |
| `style_series()` (was `style_running_sum`) | fixed axis/legend for comparability | same | `figures.{running_sum_ylim,…}` |
| `purge_figures(prefix)` | delete stale before write | DC's `purge_figures` | — |
| `write_caption()` / `append_master_table()` | idempotent caption / master append | same | schemas |

Delivery is the load-bearing propagation decision (§5): the contract ships as a **versioned helper library in the toolkit** (`lib/figure-style/figure_helpers.R` + `figure_helpers.py`), exposed to repos as a **symlink** (like skills/agents) rather than copied code — so a fix to the floors propagates to every repo on the next `activate`, exactly as agent content already does. A thin per-project shim sources it and overlays project config. The `figure-style` SKILL.md (concept) documents the contract, the dual print+screen variant pattern, and the named anti-patterns ("AdaW ggsave-inside-compute", "inline theme block", "truncated axis label"). New `viz`/`figure` tag attached (the `viz` tag currently has zero skills — `03`).

**Dual print+screen from one code path.** Per `08 Area 5`, Nature wants 5–7pt at 89mm; a projected slide wants ≥24pt. `save_figure()` gains a `variant` parameter (default both): it renders the same plot object twice from one config — a `print` variant (column width, print font tier `base_size_column`) and a `screen` variant (slide geometry, `base_size` tier) — into sibling outputs (`<name>.print.pdf` / `<name>.screen.png`). One code path, two artifacts; the floors are enforced per-variant. This directly retires the recurring "looks fine on my screen, illegible in print / vice-versa" steer.

**Guardrail.** `sciagent validate --check figure-style`: greps viz scripts for inline `theme(text`/`element_text(size=`/`ggsave(width=<number>`/raw hex `#[0-9a-f]{6}` and for absence of a `project_theme()`/`set_paper_style()` call; asserts `figures.base_size ≥ 16`. Soft-warn by default (exploration-safe), hard-fail under `--strict`. A `figure-audit` capability (subagent built from `graphic`'s Tufte checklist, but pointed at static `03_results/<stage>/figures/`) runs the D-theme checklist: truncation, ambiguous glyphs, top-N, column clustering, how-to-read present.

### B — Results placement

**Centralized (CRAFT block):** the canonical layout `03_results/<stage>/{figures,tables}/` with `by_contrast/<c>/` (per-contrast) and `_overview/` (cross-contrast); **source-table adjacency** ("a figure's source table lives at the same path stem — `_overview/<stem>.pdf` ⇒ `tables/_overview/<stem>.csv`"); strict compute→viz split. Stage IDs declared once in `analysis_config.yaml:stages`; disjoint numeric prefixes for parallel pipelines (gene `00–07`, TE `20–24`).

**Capability.** The same helper library: `save_overview()` *is* the adjacency mechanism (writes figure + sibling table + caption atomically — you cannot make the figure without the table). `contrast_path()`/`overview_path()` mkdir the sub-layout so the agent never hand-builds paths. The scaffold's `analysis_config.yaml.template` already has `stages:` and the `by_contrast`/`_overview` subdir keys are added.

**Resolving the `scrna-pipeline-conventions` contradiction.** `03` flags that this skill documents a *flat* `03_results/{tables,plots}/` layout contradicting the stage-based template. **Decision: the stage-based layout wins** (it is the owner's actual practice in all three repos and the template). We rewrite `scrna-pipeline-conventions` to the stage layout and make it *reference* the CRAFT block + figure-style contract rather than restate a competing one. The flat layout is deleted, not deprecated-in-place, to kill the SSOT split.

**Guardrail.** `sciagent validate --check results-layout`: every artifact under `03_results/` (excluding `objects/`,`master/`,`interactive/`,`_scratch/`) sits under a `<stage>/{figures,tables}/`; stage matches a `stages:` id; every figure has a same-stem table neighbor; no compute artifact at results root. Extends the existing `_validate_docs_layout` (which already warns on `.md` in results).

### C — README adjacency

**Centralized (CRAFT block, one line):** "A task is not complete until the sibling `README.md` captions every file you created/edited/deleted under `03_results/`, including *how to read* it (glyphs, sign convention, Δρ, claim tier)."

**Capability.** Extend `captions` (do not add a new agent — `03`): adopt 14839's path-qualified headings (`## figures/by_contrast/<c>/<file>`), add a mandatory **How to read** section (sign convention, glyph semantics, claim tier L0–L7), and add the framing "ALWAYS run as a mandatory cleanup pass after any figure is created/modified/deleted." `save_overview()` writes the caption atomically at figure-creation time, so the common path needs no separate agent call at all. Extend `handoff` Step 2 (already scans for uncaptioned artifacts) to *dispatch* captions when it finds gaps.

**Guardrail.** Extend `doc-curator` C3 (uncaptioned artifacts) to verify the caption's `Script:` field points to a committed path (ties C to E). Add a Claude-Code `SubagentStop`/`Stop` hook that, when the session touched `03_results/figures/`, checks the sibling README has a heading for each file and surfaces a blocking reminder. `validate --check captions` mirrors this harness-agnostically for CI.

### D — Planning decomposition

This is the convention that was *never written down* (`05 §4.4`), so the leverage is highest. The fix is two parts: (1) write the recipe into the CRAFT block + plan templates; (2) ship the orchestration commands (§4) so the owner stops re-typing the recipe.

**Centralized (CRAFT block):** plans live in `docs/_internal/plans/{date-slug}/` as `00_INDEX.md` + `NN_<slug>.md`; research → `docs/_internal/research/{date-slug}/`; reasoning → `docs/_internal/reasoning/`. The decomposition rule: "one phase == one bounded brief == one script == one Sonnet implementer (~30–40% of 200k ctx); Opus decomposes, Sonnet implements, Opus reviews every 2–3 phases checking plan adherence + that code *runs and produces its artifacts*." Add the missing `plans/` row to the namespace table.

**Capability — plan templates (promote 14839's gold-standard format).** Ship `templates/plan/00_INDEX.md.template` (title; Date·Scope IN/OUT; the "one phase==one brief" rule; phase table `seq|slug|title|tier|concern|depends_on`; `config_additions`/`toolkit_additions` block; Global notes: idempotency, compute→viz, claim ladder) and `templates/plan/NN_<slug>.md.template` (header with **figure-style contract** declaration VIZ/COMPUTE-only; §1 objective+out-of-scope-grep; §3 inputs cite research notes; §4 outputs; §5 implementation with exact helper names; §6 captions; §8 acceptance checks grep-verifiable + structural; §9 named gotchas). These templates are the durable encoding of the recipe; the commands fill them.

**Guardrail.** The acceptance gate is the Opus reviewer's runnable check, not a grep: the review step *executes* the phase's script (or confirms its logged run) and asserts the declared `03_results/` artifacts exist and are non-empty. This is the "code actually runs + produces its artifacts" requirement the owner re-types.

**Composition with existing pipeline.** `/pipeline-plan` does NOT replace `/plan`/`/implement`/`/verify` — it is a science-flavored orchestration *wrapper* that: uses the 14839 INDEX/phase templates (not the software data-models-first ordering); adds model-tiering and review cadence the software commands lack; and can delegate the per-phase mechanics to `/implement` and the drift check to `/verify`. The architect pipeline stays for software-design work; the new commands serve the analysis lane.

### E — Reproducibility / no ephemeral scripts

**Centralized (CRAFT block):** "No ephemeral scripts. Every `03_results/` artifact must be reproducible by re-running a committed `02_analysis/scripts/NN_<slug>.{R,py}`. No `/tmp` scripts, no unsaved snippets. For any non-trivial decision (method, phase boundary, parameter, surprise), write a reasoning trace to `docs/_internal/reasoning/YYYY-MM-DD_NN_<topic>.md` *before* proceeding — a decision with no trace is non-reproducible."

**Capability.** A `reasoning-trace` skill (concept, `provenance` tag) encoding the note format (title; Date·Role·Project; Scope; Sources; `## Decision N` evidence→decision→why-not) — promoting 14839's research-note→synthesis→plan flow and the mattpc "capture the answer, delete the shell" prototype discipline (`07 Pattern C`). Extend `bio-interpreter`/`insight-explorer` with "save to `research/` before returning — a chat summary is ephemeral", and `handoff` with "list every script run (committed path) + artifact produced + open decision with its reasoning trace." The multi-wave `interpret-storm` (E) is a command (§4) that persists every wave's notes to the repo before any implementation.

**Guardrail.** `PreToolUse` hook (Claude Code): block a `Bash` call that writes a `03_results/` artifact from a script outside `02_analysis/scripts/` (exit 2). `validate --check provenance`: every captioned artifact's `Script:` field resolves to an existing committed path; `_scratch/` is the only sanctioned ephemeral zone and is gitignored.

### Cross-cutting resolutions

- **R vs Python contract (14839 R helpers vs DC `plot_utils.py`):** resolved by the single named **figure-style contract** above — identical function names and semantics in both languages, both reading the same `analysis_config.yaml:figures` block, shipped as one symlinked toolkit library. Cross-language parity becomes a property of the toolkit, not a per-repo coincidence.
- **Flat vs stage layout contradiction:** resolved in favor of stage-based (§B); the contradicting skill is rewritten.

---

## 4. The orchestration-recipe system

Four parameterized commands replace the owner's repeated typing (catalog R1–R4 in `06`). Each is a slash command (user-invoked) that orchestrates Opus/Sonnet subagents with explicit model-tiering. All write durable artifacts to the repo (never chat-only). Defaults below.

### `/pipeline-plan` (R1) — plan + orchestrate a multi-phase pipeline

| Param | Default | Meaning |
|---|---|---|
| `slug` | required | `docs/_internal/plans/{date}-{slug}/` |
| `--scope-doc` | latest `research/{date}-{slug}/` | synthesis/context fed to the planner |
| `--n-planners` | 1 | parallel Opus decomposers |
| `--context-budget` | 35 | target % of Sonnet 200k per phase |
| `--review-every` | 3 | Opus review cadence (phases) |
| `--background` | auto | long non-blocking phases run in monitorable tmux |

**Flow:** Opus planner reads scope-doc + AGENTS.md (incl. CRAFT) + figure-style contract → writes `00_INDEX.md` + `NN_<slug>.md` from templates (§3D), each ≤ budget = one script. Then dependency-ordered execution: one Sonnet `/implement`-style worker per phase; after every `--review-every` phases an Opus reviewer (a) checks plan adherence/cleanliness/namespace/drift, (b) **runs the script and asserts artifacts exist** (the real gate), (c) writes the review to `reasoning/`. Long phases that don't block successors launch in tmux (`run_in_background`); later non-blocked phases proceed in parallel — dependency-aware overlap.

### `/explore-and-plan` (R2) — research fan-out → synthesis → plan

| Param | Default | Meaning |
|---|---|---|
| `question` / `slug` | required | scientific question / research slug |
| `--n-explorers` | 3 | parallel Opus explorers (codebase / literature / reference impl) |
| `--idempotency-peer` | none | sibling pipeline whose namespaces must stay disjoint |

**Flow:** Wave 1 explorers → `research/{date}-{slug}/NN_<lane>.md`. Wave 2 one Opus synthesizer → `_SYNTHESIS.md` (verified vs inferred tagged). Wave 3 → hands to `/pipeline-plan` with the synthesis as scope-doc; if `--idempotency-peer` given, planner verifies disjoint stage-ids/checkpoints/master namespaces. Encodes the 14839 "planner never re-explores; phase §3 cites the research note" discipline.

### `/add-figure-variant` (R3) — new figure family with its own namespace

| Param | Default | Meaning |
|---|---|---|
| `slug` / `stage-id` | required | family slug / target stage |
| `--namespace-token` | derived | grep-isolable token for new identifiers |
| `--n-implementers` | 1 | parallel Sonnets (split compute/viz) |

**Flow:** evidence → `research/`; Opus mini-plan (namespace decl, acceptance gate); compute Sonnet first → verify checkpoints on disk → viz Sonnet (figure-style contract, both print+screen variants); Opus review opens the PDFs / confirms non-empty panels + namespace isolation (`grep <token>` only in intended scripts) + caption completeness; mandatory `captions` cleanup pass.

### `/interpret-storm` (R4) — multi-wave interpretation + information-graphic design

| Param | Default | Meaning |
|---|---|---|
| `dataset-path` / `slug` | required | results to interpret |
| `--n-interpreters` | 3 | Opus interpreters (mechanism / pathway / clinical angles) |
| `--n-designers` | 2 | information-graphic / interactive-viewer designers |

**Flow:** Wave 1 web research → `research/{date}-{slug}/web-*.md`. Wave 2 interpreters → `reason-*.md` + Opus `_SYNTHESIS.md` (verified vs inferred). Wave 3 designers each propose a figure with namespace token + data contract + claim tier → `reasoning/{date}-{slug}-figures.md`. Wave 4 hands to `/add-figure-variant`. **Every wave persists to the repo before the next begins** — ephemerality is the failure mode this command exists to prevent.

### Role placement

- `/pipeline-plan`, `/explore-and-plan`, `/add-figure-variant`, `/interpret-storm` → a **new `science-architect` overlay role** (so analysis repos can `sciagent activate base science-architect` and get the orchestration suite without the full software-architect command set). They are NOT in `base` (they are deliberate, heavyweight, user-invoked) and NOT in `architect` (that lane is software-design-shaped).
- `reasoning-trace` skill, `figure-style` skill, the helper library, the extended `captions`/`doc-curator`/`handoff`/`bio-interpreter`/`insight-explorer` → `base` (always available in day-to-day analysis).
- The CRAFT block is rendered for **every** role (it is disposition, not a role feature).

---

## 5. Propagation & drift control

The diagnosis (`05`): symlinks already make agent/command/skill *content* drift-proof, but (i) AGENTS.md bodies are hand-frozen template copies, (ii) helper libraries are copied per-repo, (iii) repos pin stale toolkit commits (STING/DC 24 behind). The fixes:

1. **CRAFT as a managed block.** Promote craft to `<!-- BEGIN SCIAGENT:CRAFT v1 hash=… -->`, rendered by `activate`/`inject`/`eject` next to ROLES, SHA1 drift-detected, idempotently re-rendered. Mechanically: `block.sh` markers are currently the literal `SCIAGENT:ROLES`; parameterize the primitives to take a block-id (`block_read <file> <id>`, `block_write <file> <id> <body>`, etc.) and add a `craft.sh` renderer that emits the canonical craft text + the numeric floors pulled from a toolkit-owned `craft.yaml` SSOT. `activate` calls `block_render_and_write AGENTS.md ROLES` then `craft_render_and_write AGENTS.md`. Re-`activate` becomes the propagation verb: a floor change in `craft.yaml` reaches every repo on next activate. **This is the single highest-leverage change** — it moves the four un-centralized conventions into the one mechanism already proven to not drift.
2. **figure-style as a symlinked capability, not copied code.** The helper library lives in `lib/figure-style/` and is exposed via a symlink the same way skills are (`activate` adds `.sciagent`/helpers link or a `02_analysis/helpers/figure_helpers.{R,py}` symlink into the submodule). A fix to `save_figure()` propagates without per-repo edits. Per-project values stay in `analysis_config.yaml:figures`; the code is shared, the config is local — the correct seam.
3. **`sciagent update` / re-pin story.** Add a thin `sciagent update` verb: bump the repo's `01_modules/SciAgent-toolkit` submodule to a target commit, re-run `activate` (re-render both managed blocks + re-link helpers), and report what changed (new skills available, craft version bump). This addresses the 24-commit staleness directly and makes "get the latest conventions" a one-liner.
4. **Staleness validate check.** `validate --check freshness` compares the repo's pinned submodule commit and CRAFT block version against the toolkit HEAD and warns when behind (with the `sciagent update` hint). Also re-flags the DC archived-skill production link (`05 §1.6`).
5. **Idempotency.** `activate` re-render of both blocks must be a no-op when nothing changed (CRAFT hash matches) — same invariant the ROLES block already holds; covered by a new roundtrip test.

---

## 6. New tags / vocabulary

`tags.yaml` gains entries (closed-vocab, `validate`-enforced):

- **`figure`** — figure-style contract, legibility, dual-variant rendering, caption craft. (The `viz` tag exists but has zero skills and is plotting-generic; `figure` is the publication-figure-craft family. Attach `figure-style` skill to both `figure` and `viz`.)
- **`provenance`** — reasoning-trace persistence, no-ephemeral discipline, decision logs, session continuity. Home for `reasoning-trace`.
- **`planning`** — plan decomposition, orchestration recipes, phase templates. Home for the orchestration commands' companion skills (if any) and the plan templates' doc skill.

No tag is removed. `viz` keeps its definition; `figure` is the narrower craft family.

---

## 7. Risks, trade-offs, and explicit non-goals

**Risks / where this could get heavy or brittle:**
- *Block.sh parameterization risk.* Generalizing the single-block primitives to two blocks touches the most-tested, load-bearing file. Mitigation: keep ROLES bytes identical (regression tests are the contract); add CRAFT behind the same idempotency/drift invariants; ship it as its own phase with an Opus review checkpoint.
- *AGENTS.md budget.* CRAFT + ROLES must stay under the ~150-line/"lost-in-the-middle" budget. Mitigation: CRAFT ≤25 lines, floors in config, depth in skills; measure in a validate soft-warn.
- *Over-rigid guardrails chilling exploration.* Hard-failing on inline `ggsave` or "no `/tmp`" can block legitimate quick exploration. Mitigation: guardrails default to **soft-warn**; `--strict` (and CI) hard-fail. `_scratch/` and `$TMPDIR` remain sanctioned ephemeral zones (mattpc `07`). The no-ephemeral rule targets *result-producing* scripts, not throwaway probes whose answer is captured in a reasoning trace.
- *Symlinked helper library vs harness portability.* A symlinked `figure_helpers` assumes the submodule is present. Mitigation: thin per-project shim with a graceful fallback (`theme_minimal`) when the toolkit isn't linked, as 14839's `project_theme()` already does.
- *Orchestration commands assume Claude Code subagent/background semantics.* On other harnesses tmux-overlap degrades to sequential. Mitigation: the commands state the recipe as instructions (provider-agnostic), and only the background-execution optimization is Claude-specific.
- *Cost.* Opus-decompose + Opus-review-every-N is expensive. Mitigation: cadence is a parameter; Sonnet does the bulk; `--n-explorers`/`--n-planners` bounded.

**Explicit non-goals:**
- Not adopting Snakemake/Nextflow DAGs (`08` suggests it) — out of scope; the numbered-script + load_or_compute pattern is the owner's established, working idiom.
- Not auto-generating AGENTS.md prose with an LLM (the ETH finding: LLM-written instruction files *reduce* success — `08`). CRAFT is hand-authored once, toolkit-owned.
- Not adding a third stack slot or breaking the 2-role cap.
- Not editing existing repos beyond an optional, opt-in backfill phase.
- Not building a generic provenance/PROV-emitting subsystem — reasoning traces are markdown, by design lightweight (academic "minimum viable MLOps", `08 Area 4`).
- Not replacing `/plan`/`/implement`/`/verify` — compose, don't duplicate.

---

## 8. Acceptance criteria — "no longer needs hand-steering"

Per pain point, concrete and checkable:

- **A. Figure legibility.** `sciagent new project` → AGENTS.md already contains the dual-context CRAFT rule + `figures.base_size ≥ 16`; the figure-style helper library is linked; `save_figure()` emits print+screen variants from one call. `validate --check figure-style` passes on a conformant repo and flags an inline `theme()`/`ggsave(width=)`/raw-hex viz script. `figure-audit` returns the D-theme checklist. ⇒ The owner never re-types "make it legible projected / not truncated / centralize the theme."
- **B. Results placement.** Every artifact a fresh agent produces lands under `<stage>/{figures,tables}/` with `by_contrast/`+`_overview/` and a same-stem source table (because `save_overview` is the only sanctioned path). `scrna-pipeline-conventions` no longer contradicts the layout. `validate --check results-layout` passes. ⇒ No re-typing of the layout.
- **C. README adjacency.** A session that touches a figure leaves the sibling README captioned (incl. how-to-read) — enforced by `save_overview` on the common path, the `Stop` hook + `doc-curator` C3 as backstop. `validate --check captions` passes. ⇒ No re-typing of "update the README to explain the figure."
- **D. Planning decomposition.** `/pipeline-plan slug` produces `00_INDEX.md` + sized `NN_<slug>.md` phases from the gold-standard templates, with tier labels and Opus review checkpoints every N, and the review step actually runs scripts and checks artifacts. ⇒ The owner types `/pipeline-plan`, not the whole recipe.
- **E. Reproducibility.** No result is produced by a `/tmp`/unsaved script (PreToolUse hook + validate). Non-trivial decisions appear in `reasoning/` before implementation. `/interpret-storm` persists every wave. `validate --check provenance` confirms every captioned artifact cites a committed script. ⇒ No re-typing of "save the script / write down why / don't run throwaway code."

**Global meta-acceptance:** all five conventions have exactly one canonical home in the toolkit (`craft.yaml` + templates + helper library + skills), and a craft-rule change propagates to every repo via `sciagent update` + re-`activate` with no per-repo hand-edit. The drift table in `05 §5` flips every "Not centralized" row to "Centralized."
