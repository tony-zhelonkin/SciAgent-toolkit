# Plan INDEX — Reproducible Agentic Science: centralize, inject, enforce the five conventions

**Date:** 2026-06-23 · **Scope:** IN — SciAgent-toolkit changes that make the owner's five conventions (figure legibility, results placement, README adjacency, planning decomposition, reproducibility) centralized + injected + enforced so a dropped-in agent follows them without hand-steering. OUT — per-project edits (except opt-in backfill P15), new analysis methods, Snakemake/DAG adoption, LLM-authored AGENTS.md, replacing `/plan`/`/implement`/`/verify`.

**Rule:** one phase == one bounded brief == one Sonnet implementer (~30–40% of 200k ctx = one coherent change). Numbering is authoritative by **slug**, not by reading order. Opus REVIEW CHECKPOINT rows sit every 2–3 phases and each verifies: plan adherence, cleanliness/namespace, AND that `bash tests/run-all.sh` passes + the new artifacts exist and behave.

See `00_ARCHITECTURE.md` for design rationale. This INDEX is rich enough that a follow-up Opus can expand each `NN_<slug>.md` into a full brief without re-exploring.

---

## Phase table

| seq | slug | title | tier | concern | depends_on |
|----|------|-------|------|---------|------------|
| 01 | craft-block-primitives | Parameterize `block.sh` to support a second managed block id; keep ROLES bytes identical | Sonnet | A-E (mechanism) | — |
| 02 | craft-yaml-renderer | Add `craft.yaml` SSOT + `craft.sh` renderer; wire into `activate`/`inject`/`eject` | Opus | A-E (mechanism) | 01 |
| 03 | craft-content | Author the canonical CRAFT block text (all five conventions, terse, ≤25 lines) | Opus | A-E (rule) | 02 |
| R1 | review-craft-mechanism | **Opus REVIEW CHECKPOINT** — blocks render idempotently, no ROLES regression, tests green | Opus | review | 01,02,03 |
| 04 | config-figure-floors | Add `figures:` floors + `by_contrast`/`_overview` keys to `analysis_config.yaml.template` | Sonnet | A,B | 03 |
| 05 | figure-style-lib | Ship `lib/figure-style/figure_helpers.{R,py}` — unified cross-language contract | Opus | A,B | 04 |
| 06 | figure-style-link | Symlink the helper lib into repos via `activate`; thin per-project shim + fallback | Sonnet | A,B (propagation) | 05 |
| 07 | figure-style-skill | Add `figure-style` SKILL.md (concept) documenting contract + dual print/screen + anti-patterns | Sonnet | A,B | 05 |
| R2 | review-figure-stack | **Opus REVIEW CHECKPOINT** — helpers run in R+py, emit both variants, write adjacent table+caption | Opus | review | 04,05,06,07 |
| 08 | reconcile-scrna-layout | Rewrite `scrna-pipeline-conventions` to stage layout; delete flat-layout contradiction | Sonnet | B | 03 |
| 09 | tags-vocabulary | Add `figure`/`provenance`/`planning` tags to `tags.yaml` | Sonnet | A,D,E (vocab) | — |
| 10 | extend-captions | Extend `captions`: path-qualified headings, how-to-read section, mandatory-after-figure framing | Sonnet | C | 03 |
| 11 | extend-curator-handoff | Extend `doc-curator` C3 (Script-path committed) + `handoff` (scripts/artifacts/decisions) | Sonnet | C,E | 10 |
| 12 | figure-audit-agent | Add `figure-audit` agent (static-figure D-theme checklist from `graphic`'s lens) to base | Sonnet | A,C | 07,10 |
| R3 | review-cdocs-figaudit | **Opus REVIEW CHECKPOINT** — captions/curator/handoff/figure-audit run on a fixture stage dir | Opus | review | 08,10,11,12 |
| 13 | validate-checks | Add `validate` checks: figure-style, results-layout, captions, provenance, freshness | Opus | A,B,C,E (guardrail) | 04,08,10 |
| 14 | enforcement-hooks | Ship harness hooks: PreToolUse (no-ephemeral, figure-save warn) + Stop (caption sweep) | Sonnet | C,E (guardrail) | 13 |
| 15 | sciagent-update-verb | Add `sciagent update` (re-pin submodule + re-activate + report) | Sonnet | propagation | 02,06 |
| R4 | review-guardrails | **Opus REVIEW CHECKPOINT** — checks fail closed on violations, soft-warn by default, hooks fire | Opus | review | 13,14,15 |
| 16 | plan-templates | Ship `templates/plan/{00_INDEX,NN_slug}.md.template` (14839 gold-standard format) | Sonnet | D | 03 |
| 17 | reasoning-trace-skill | Add `reasoning-trace` SKILL.md + extend bio-interpreter/insight-explorer to persist | Sonnet | E | 09 |
| 18 | science-architect-role | New `science-architect` overlay role; wire orchestration commands into it | Sonnet | D | 16 |
| 19 | cmd-pipeline-plan | `/pipeline-plan` command (Opus-decompose → Sonnet-impl → Opus-review-every-N + runnable gate) | Opus | D | 16,18 |
| 20 | cmd-explore-and-plan | `/explore-and-plan` command (research fan-out → synthesis → /pipeline-plan) | Sonnet | D,E | 17,19 |
| R5 | review-planning-suite | **Opus REVIEW CHECKPOINT** — commands produce gold-standard plans; review step runs scripts | Opus | review | 16,17,18,19,20 |
| 21 | cmd-add-figure-variant | `/add-figure-variant` command (namespace + compute→viz fan-out + PDF-open review) | Sonnet | A,D | 07,19 |
| 22 | cmd-interpret-storm | `/interpret-storm` command (web → interpret → graphic-design waves, all persisted) | Opus | D,E | 17,20,21 |
| 23 | wire-craft-into-scaffold | Ensure `sciagent new project` + AGENTS.md template reference CRAFT/helpers; retire dead viz SSOT | Sonnet | A-E (integration) | 03,05,16 |
| R6 | review-orchestration | **Opus REVIEW CHECKPOINT** — full new-project smoke: scaffold→activate→plan→implement→figure→caption | Opus | review | 21,22,23 |
| 24 | backfill-repos (opt-in) | Backfill STING/DC/14839: re-activate (CRAFT), fix STING template, DC archived-skill, missing READMEs | Sonnet | A-E (rollout) | 23 |

**24 substantive phases + 6 Opus review checkpoints (R1–R6).** Checkpoints fall after every 2–4 phases at each subsystem boundary: mechanism (R1), figure stack (R2), captions/audit (R3), guardrails (R4), planning suite (R5), full integration (R6).

---

## toolkit_additions (new/changed files under the toolkit)

**lib/ (mechanism + guardrails)**
- `lib/sciagent/block.sh` — parameterize markers by block id (ROLES + CRAFT). [P01]
- `lib/sciagent/craft.sh` — `craft_render_and_write`, reads `craft.yaml`. [P02]
- `lib/sciagent/activate.sh`, `inject.sh`, `eject.sh` — call `craft_render_and_write` after ROLES. [P02]
- `lib/sciagent/validate.sh` — new `--check {figure-style,results-layout,captions,provenance,freshness}`. [P13]
- `lib/sciagent/update.sh` + bin/sciagent dispatch — `sciagent update`. [P15]
- `craft.yaml` (toolkit root) — SSOT for craft text params + numeric floors. [P02]

**lib/figure-style/ (symlinked helper library)**
- `lib/figure-style/figure_helpers.R`, `figure_helpers.py` — `project_theme`/`set_paper_style`, `save_figure` (print+screen variants), `save_overview`, `contrast_path`/`overview_path`, `style_series`, `purge_figures`, `write_caption`, `append_master_table`. [P05]

**skills/**
- `skills/figure-style/SKILL.md` — contract + dual-variant + anti-patterns; tags `[figure, viz]`. [P07]
- `skills/reasoning-trace/SKILL.md` — note format + no-ephemeral; tags `[provenance]`. [P17]
- `skills/scrna-pipeline-conventions/SKILL.md` — rewrite to stage layout. [P08]

**agents/analysis-base/**
- `figure-audit.md` — static-figure D-theme checklist reviewer (model sonnet). [P12]
- `captions.md`, `doc-curator.md`, `handoff.md`, `bio-interpreter.md`, `insight-explorer.md` — extensions. [P10,P11,P17]

**commands/**
- `commands/science/pipeline-plan.md`, `explore-and-plan.md`, `add-figure-variant.md`, `interpret-storm.md`. [P19,P20,P21,P22]

**roles/**
- `roles/science-architect.yaml` — overlay with the 4 orchestration commands. [P18]
- `roles/base.yaml` — add `figure-style`, `reasoning-trace` skills + `figure-audit` agent. [P07,P12,P17]

**templates/**
- `templates/plan/00_INDEX.md.template`, `NN_slug.md.template`. [P16]
- `templates/project/analysis/02_analysis/config/analysis_config.yaml.template` — `figures:` floors + subdir keys. [P04]
- `templates/project/analysis/AGENTS.md.template` — reference CRAFT/helpers; namespace-table `plans/` row. [P23]
- `templates/project/_common/.claude/hooks/` (or equivalent) — hook scripts + settings seed. [P14]

**tags.yaml** — add `figure`, `provenance`, `planning`. [P09]

**docs/** — retire/redirect `docs/guidelines/visualization.md` to the figure-style contract (kill the `base_size=12` dead-end). [P23]

**tests/** — `test_craft_block_roundtrip`, `test_craft_block_drift`, `test_activate_renders_craft`, `test_figure_helpers_contract`, `test_validate_figure_style`, `test_validate_results_layout`, `test_validate_provenance`, `test_update_verb`, `test_science_architect_role`, `test_tags_vocabulary` (extend). [each phase]

---

## Global notes

- **Enforcement philosophy: hooks/validate > goodwill.** Every convention is stated once (CRAFT), made doable (skill/helper/command), made unskippable (validate check + hook). Guardrails default soft-warn; `--strict`/CI hard-fail. `_scratch/` + `$TMPDIR` stay sanctioned ephemeral zones — the no-ephemeral rule targets *result-producing* scripts only.
- **Single source of truth.** `craft.yaml` is the one home for craft text + floors; `lib/figure-style/` is the one home for figure code; `templates/plan/` is the one home for the plan format; `tags.yaml` is the closed vocab. Deletions (flat scrna layout, dead viz SSOT) are part of the work — a second home is a drift source.
- **Propagation = re-activate.** Both managed blocks re-render on `activate`; the helper lib is symlinked not copied; `sciagent update` re-pins + re-activates so a floor change reaches every repo with no per-repo edit. Idempotent: re-activate with unchanged inputs is a no-op (CRAFT hash matches), same invariant ROLES already holds.
- **Dual print+screen contract.** `save_figure(..., variant="both")` renders one plot object twice from one config — `print` (column width, 5–7pt floor) and `screen` (slide geometry, ≥24pt floor) — enforcing floors per variant. One code path, two artifacts.
- **Backward-compat.** `/plan`/`/implement`/`/verify` and the architect role are untouched; `/pipeline-plan` composes on top (science-flavored INDEX/phase templates + model-tiering + runnable gate the software commands lack). ROLES block bytes stay identical (regression-tested) when CRAFT is added.
- **Idempotency + compute→viz (carried into every plan this toolkit generates).** Numbered scripts only; `load_or_compute` caching; compute scripts never plot, viz scripts never compute (named anti-pattern "ggsave-inside-compute"); idempotent `append_master_table`/`write_caption`.

---

## Inline phase sketches (6 most important; others expandable by a follow-up Opus)

### P01 — craft-block-primitives  *(Sonnet)*
- **Objective:** make `block.sh` handle two managed blocks (`ROLES`, `CRAFT`) without changing ROLES behavior.
- **Scope:** parameterize `block_read/write/remove/hash_check/stored_hash/line_range/render_and_write` to take a block id; keep a back-compat default of `ROLES`. OUT: writing CRAFT content (P03) or wiring activate (P02).
- **Key files:** `lib/sciagent/block.sh`, `lib/sciagent/stack.sh` (`block_render_and_write` signature), tests.
- **Acceptance:** `bash tests/run-all.sh` green; existing ROLES roundtrip/drift tests pass byte-identical; new `test_craft_block_roundtrip` writes+reads+drift-detects a CRAFT block independently of ROLES; both blocks coexist in one AGENTS.md.

### P05 — figure-style-lib  *(Opus)*
- **Objective:** one cross-language figure contract implemented in R and Python, both reading `analysis_config.yaml:figures`.
- **Scope:** `project_theme`/`set_paper_style`, `save_figure` (print+screen variants, cairo_pdf + png@300), `save_overview` (atomic figure+table+caption), `contrast_path`/`overview_path`, `style_series`, `purge_figures`, `write_caption`, `append_master_table` (idempotent, round_numeric_cols). OUT: the skill doc (P07), the symlink wiring (P06). Mine 14839 R helpers + DC `plot_utils.py` for exact semantics.
- **Key files:** `lib/figure-style/figure_helpers.R`, `figure_helpers.py`.
- **Acceptance:** in a fixture project, an R viz script and a py viz script each call `save_overview(...)` and produce `<stem>.print.pdf`+`<stem>.screen.png`+`tables/_overview/<stem>.csv`+README caption; floors enforced per variant; `test_figure_helpers_contract` asserts function parity across languages.

### P13 — validate-checks  *(Opus)*
- **Objective:** turn the four conventions into deterministic, fail-closed checks plus a freshness warning.
- **Scope:** `figure-style` (grep inline `theme(text`/`ggsave(width=<num>`/raw hex in `*_viz.*`; require a `project_theme`/`set_paper_style` call; assert `figures.base_size ≥ 16`); `results-layout` (artifact under `<stage>/{figures,tables}`; stage in `stages:`; figure has same-stem table); `captions` (every `03_results` artifact has a sibling-README heading); `provenance` (caption `Script:` resolves to committed path); `freshness` (pinned commit / CRAFT version vs HEAD). Soft-warn default; `--strict` hard-fail. Extends `_validate_docs_layout`.
- **Key files:** `lib/sciagent/validate.sh`, tests.
- **Acceptance:** each check passes on a conformant fixture and fails (or warns under default) on a crafted violation; `tests/run-all.sh` green; no false positives on the existing toolkit tree.

### P16 — plan-templates  *(Sonnet)*
- **Objective:** ship the 14839 gold-standard plan format as reusable templates.
- **Scope:** `00_INDEX.md.template` (title; Date·Scope IN/OUT; one-phase==one-brief rule; phase table `seq|slug|title|tier|concern|depends_on`; `config_additions`/`toolkit_additions`; Global notes idempotency/compute→viz/claim-ladder) and `NN_slug.md.template` (header w/ figure-style contract VIZ/COMPUTE declaration; objective+out-of-scope-grep; inputs cite research notes; outputs; implementation w/ exact helper names; captions; acceptance grep+structural; named gotchas). OUT: the commands that fill them (P19+).
- **Key files:** `templates/plan/00_INDEX.md.template`, `templates/plan/NN_slug.md.template`.
- **Acceptance:** templates render with `_subst`; a manual `/pipeline-plan` dry-run produces a structurally valid INDEX + phase brief matching the 14839 section schema.

### P19 — cmd-pipeline-plan  *(Opus)*
- **Objective:** replace the owner's hand-typed orchestration recipe with one command.
- **Scope:** params `slug` (req), `--scope-doc`, `--n-planners 1`, `--context-budget 35`, `--review-every 3`, `--background auto`. Flow: Opus planner → INDEX+phases from P16 templates (each ≤ budget = one script); dependency-ordered Sonnet implementers; Opus reviewer every N that checks adherence/namespace/drift AND **runs the script + asserts artifacts exist**; long non-blocking phases → tmux with dependency-aware overlap; reviews → `reasoning/`. Composes with `/implement`/`/verify`; does not replace them.
- **Key files:** `commands/science/pipeline-plan.md`, `roles/science-architect.yaml`.
- **Acceptance:** on a toy scope-doc, produces a date-slug plan dir with tiered phases and an Opus review row at phase N; the review step demonstrably runs a script and reports artifact presence; model-tiering explicit (planner Opus, impl Sonnet, review Opus).

### P14 — enforcement-hooks  *(Sonnet)*
- **Objective:** deterministic Claude-Code enforcement for the two conventions LLMs break under pressure (ephemerality, caption-skip).
- **Scope:** `PreToolUse` hook — block (exit 2) a Bash invocation that writes a `03_results/` artifact from a path outside `02_analysis/scripts/` (warn-not-block by default; block under strict); figure-save warn nudging the figure-style contract. `Stop`/`SubagentStop` hook — if the session touched `03_results/.../figures/`, check sibling README captions and surface a reminder. Ship as `templates/.../hooks/*.sh` + a `settings.local.json` seed registering them; degrade gracefully on non-Claude harnesses (validate covers CI).
- **Key files:** `templates/project/_common/.claude/hooks/`, settings seed, docs.
- **Acceptance:** hook scripts are self-contained bash (no deps), exit 2 on a crafted `/tmp` result-write under strict, exit 0 (warn) by default; `Stop` hook flags an uncaptioned figure in a fixture; not registered globally (project-scoped only).
