# Field inspection report

Counts use `find docs/_internal -type f`, including files inside nested `.git/` directories. “Tracked” is the exact parent-project command requested: `git -C <project> ls-files docs/_internal | wc -l`.

| Project | `docs/_internal` files | Parent tracked | Immediate subdirectories | Harness footprints present | Root `_scratch/` |
|---|---:|---:|---|---|---|
| Meta-Aging | 20 | 0 | `architecture`, `data`, `handoff`, `notes`, `overview`, `plans`, `reasoning`, `skills`, `web` | `.claude/`, `AGENTS.md`, `CLAUDE.md` | No |
| 14616-DM | 45 | 0 | `ai-generated`, `documentation`, `explainers`, `plans`, `reasoning` | `.claude/`, `.agents/`, `.gemini/`, `AGENTS.md`, `CLAUDE.md` | No |
| 14761-DM | 1,161 | 0 | `.git`, `.pytest_cache`, `handoffs`, `plans`, `reasoning`, `reports`, `research`, `sessions` | `.claude/`, `.agents/`, `AGENTS.md`, `CLAUDE.md` | No |
| 14782-DM | 2,881 | 0 | `.git`, `.pytest_cache`, `consensus-validation`, `handoffs`, `loop`, `plans`, `reasoning`, `reports`, `research`, `sessions` | `.claude/`, `.agents/`, `AGENTS.md`, `CLAUDE.md` | No |
| neuroimmune-receptor-atlas | 7 | 0 | `plans`, `reasoning`, `research`, `sessions` | `AGENTS.md`, `CLAUDE.md`; no harness directory | No |

All five projects have `AGENTS.md` and `CLAUDE.md`. None has `.codex/`.

## Per-project detail

### Meta-Aging

The 20 files occupy 424 KB. Three files sit directly in `docs/_internal/`; most named subdirectories are flat:

- `architecture` 4, `data` 1, `handoff` 1, `overview` 2, `reasoning` 1, `skills` 1.
- `notes` 2, `plans` 1, and `web` 4 place their files one level deeper in dated/topic folders.
- The web-research folder is literally named `2026-07-01 subject research`, including spaces.

Ignore provenance is hand-written: `.gitignore:65` has `_internal/` under “Internal / private notes.” There is no SCIO-managed marker around it.

`.claude/` contains two real files and no mounts: a large project-specific `settings.local.json` permission history and a live-looking `scheduled_tasks.lock`. These are clear harness-local accumulation.

Samples:

- `architecture/2026-07-01_meta-aging_structure_proposal.md` — durable architecture/decision record, although parts retain proposal status.
- `plans/2026-07-09-tf-motif-concordance-dive/README.md` — implementation plan.
- `reasoning/2026-08-19_integration-permouse-regulon-scoring-design.md` — durable decision/design record.
- `handoff/2026-08-20_permouse-regulon-arc.md` — current session handoff.
- `notes/.../README.md` — durable synthesis of an annotation ruling.
- `notes/.../transcript.md` — raw recovered Claude transcript; archival evidence with much lower information density than its companion summary.

No exact Claude scratchpad path was found. One document prescribes ordinary `/tmp/agy_prompt.txt`. No notebook contained a `/tmp/` reference. A non-root scratch area does exist at `integration/03_results/_scratch`.

### 14616-DM

The 45 files occupy 728 KB:

- `ai-generated` 19, nested as far as `ai-generated/plans/<date-slug>/...`.
- `documentation` 6, all under `documentation/docs.scarches.org/`.
- `explainers` 3 and `reasoning` 8 are flat.
- `plans` 9 uses `plans/<date-slug>/NN_*.md`.

Ignore handling has two rules. A hand-written `docs/_internal/` appears at line 19, while `git check-ignore -v` reports line 68 inside the current `SCIO:GITIGNORE` managed block because it is the last matching rule.

Harness footprint:

- `.claude/`: 6 regular files and 3 mounted symlinks. The real files include settings, local permissions, statusline, scheduled-task state, and two guardrail hooks.
- `.agents/`: 3 symlinks and no regular files.
- `.gemini/`: one real `settings.json` containing project MCP configuration.

Samples:

- `ai-generated/workflow-architecture/03-decisions.md` — durable architectural decisions.
- `ai-generated/plans/HANDOFFS.md` — append-only multi-agent execution journal; useful handoff history, but already 100 KB and mixes many completed phases.
- `plans/.../00_INDEX.md` — detailed annotation-transfer plan.
- `reasoning/2026-06-27_02f_consensus_decisions.md` — concise durable owner-decision record.
- `documentation/...treeArches-learning-celltype-hierarchy.md` — copied tutorial output, including notebook warnings and `/tmp/ipykernel_*` traces; reference artifact rather than project memory.
- `explainers/2026-06-27_foundations_faq.md` — durable project-specific explainer.

Seven files mention `/tmp/`, mainly logs, validation scripts, and copied notebook warnings. No exact `/tmp/claude-.../scratchpad/` path was found, and no notebook matched. Non-root scratch directories exist at `02_analysis/scripts/_scratch` and `03_results/_scratch`.

### 14761-DM

The headline 1,161-file count consists of:

- 960 nested `.git/` files.
- 4 `.pytest_cache` files.
- 197 other files: 3 immediate files, `plans` 122, `reasoning` 57, `reports` 3, `research` 2, `sessions` 9, and the empty-placeholder file under `handoffs`.

Content depth is concentrated in `plans/<date-slug>/`, with `_handoff/` and `_reviews/` beneath plan folders. There are 38 plan-local handoff files and 16 review files. Reasoning, reports, research, and sessions are otherwise flat.

The parent project ignores the tree through the legacy managed `SCIAGENT:GITIGNORE` block at `.gitignore:78`. The parent tracks zero files. However, `docs/_internal` is itself a clean nested Git repository:

- 195 files tracked in the nested repository.
- 2 current untracked reasoning files.
- 4 ignored cache files.
- Recent nested commits record reasoning and session handoffs.

Thus “zero tracked” is true only from the parent repository’s perspective.

Harness footprint:

- `.claude/`: 5 regular project guardrail/configuration files and 71 mounted symlinks.
- `.agents/`: 71 symlinks and no regular files.
- The regular Claude hook README is explicitly project-scoped; whether it was authored manually or materialized by tooling cannot be established from the project alone.

Samples:

- `sessions/2026-08-19_stage10-rederivation-and-the-first-supported-program.md` — strong session handoff with commits, findings, and next work.
- `reasoning/2026-08-19_stage10-rederivation-design-consult.md` — durable design record.
- `reports/2026-07-31_annotation-audit.md` — completed, evidence-heavy audit report.
- `plans/2026-08-06_annotation-first-refactor/STATE.md` — mutable plan position/handoff.
- `plans/2026-07-03-batch-preproc-annotate/00_INDEX.md` — historical plan; later state records show parts have been superseded, although this file is not clearly labeled as such.
- `plans/.../_handoff/2026-08-03_spec_leiden_t_subset_lens.md` — one-off scratch probe/specification that explicitly records an external Claude scratchpad.

Nine files mention temporary paths. The exact observed external pattern appears here:

`/tmp/claude-788715489/-workspaces-Meta-Aging/<uuid>/scratchpad/`

No notebook matched. `03_results/_scratch` also exists.

### 14782-DM

This is the largest and least document-like tree: 2,881 files and 231 MB.

- `.git` 2,161 files, 26 MB.
- `.pytest_cache` 4 files.
- `consensus-validation` 57 files, 196 MB.
- `loop` 173, `plans` 158, `reasoning` 89, `reports` 218, `research` 8, `sessions` 9.
- `handoffs` contains only `.gitkeep`.

Depth reaches six levels in `consensus-validation`, four under `plans`, and three under `loop`. The tree includes 212 JSON files, 25 CSVs, 5 logs, 3 Parquet files, Python source/bytecode, and a 196 MB model checkpoint. `reports/arbiter_cache_20260730` alone contains 210 hashed cache JSON files.

The parent ignores the tree through the legacy `SCIAGENT:GITIGNORE` block at `.gitignore:82` and tracks zero files. As in 14761-DM, `docs/_internal` is a nested Git repository:

- 708 files tracked in the nested repository.
- No untracked files.
- 12 ignored cache/generated files.

Harness footprint:

- `.claude/`: 6 regular files and 53 mounted symlinks.
- `.agents/`: 52 symlinks plus 8 regular files under the locally accumulated `marimo-pair` skill.
- `settings.local.json` is project-specific harness state. The provenance of the regular `marimo-pair` files cannot be determined from content alone.

Samples:

- `sessions/2026-08-03_palette-figures-widgets.md` — session handoff.
- `reasoning/annotation-architecture-rulings.md` — explicitly binding durable decision record.
- `plans/.../_superseded/01_architecture.md` — stale plan retained in an explicitly superseded directory.
- `loop/cards/73_alt9_eligibility/01_REPORT.md` — completed operational report built around a session-local `/tmp/73/` probe.
- `consensus-validation/.../out_annotate.log` — generated runtime log.
- One `reports/arbiter_cache_20260730/<hash>.json` — generated cache entry, not durable memory.

Fourteen files in the searched locations mention temporary paths. Exact Claude scratchpad paths occur in loop review reports, and several reports explicitly say scratch scripts lived outside the repository. No notebook matched. `03_results/_scratch` also exists.

A within-project naming conflict is documented in the binding rulings file: the top-level internal README specifies dateless ADR names, while `reasoning/README.md` specifies dated `YYYY-MM-DD_NN_topic.md` names.

### neuroimmune-receptor-atlas

This is a 36 KB scaffold rather than an accumulated memory tree:

- Three nonempty immediate files: `README.md`, `notes.md`, and `scientific-context.md`.
- `plans`, `reasoning`, `research`, and `sessions` each contain only `.gitkeep`.
- All content is flat.

Ignore provenance is the legacy managed `SCIAGENT:GITIGNORE` block at `.gitignore:75`. There are no harness-specific directories.

Samples:

- `scientific-context.md` — durable scientific framing.
- `README.md` — durable layout and naming convention scaffold.
- `notes.md` — short, typo-heavy scratch outline.

No temporary-path references were found in internal docs, analysis, or notebooks. A non-root `03_results/_scratch` directory exists.

## Consolidated naming variants

| Concept | Spellings/layouts observed | Projects |
|---|---|---|
| Session handoff | `handoff/`; root `handoff.md`; dated root `*_handoff.md`; `HANDOFFS.md`; `handoff_phase*.md`; `handoffs/`; `sessions/`; plan-local `_handoff/` | `handoff/` and root files: Meta-Aging; uppercase/file variants: 14616-DM; `handoffs/`: 14761-DM and 14782-DM, both placeholder-only; `sessions/`: 14761-DM, 14782-DM, neuroimmune; `_handoff/`: 14761-DM |
| Reasoning/decisions | `reasoning/`; `reports/`; `*-decisions.md`; `annotation-architecture-rulings.md`; `workflow-architecture/` | `reasoning/`: all five; `reports/`: 14761-DM and 14782-DM; decision/ruling filenames: 14616-DM and 14782-DM; workflow architecture: 14616-DM. No `decisions/` directory was observed |
| Plans | `plans/`; `ai-generated/plans/`; operational `loop/cards/` and `loop/proposed/` | `plans/`: all five; legacy AI-generated plans: 14616-DM; loop/card plans: 14782-DM |
| Dated plan/topic directory | `YYYY-MM-DD-slug`; `YYYY-MM-DD_slug`; `YYYY-MM-DD__slug`; `YYYY-MM-DD subject research` | Single hyphen: Meta-Aging, 14616-DM, 14761-DM, 14782-DM; underscore: Meta-Aging notes and newer 14761/14782 plans; double underscore: 14616-DM legacy plans; spaces: Meta-Aging web research |
| Reviews/audits | top-level `reports/`; plan-local `_reviews/`; ordinary `reviews/`; `*_audit.md` files | `reports/`: 14761-DM and 14782-DM; `_reviews/`: both; `reviews/`: 14782-DM loop/threads; audit filenames: both |
| Research/source material | `research/`; `web/`; `documentation/`; `notes/`; root `notes.md`; `explainers/` | `research/`: 14761-DM, 14782-DM, neuroimmune, plus an empty nested legacy directory in 14616-DM; `web/` and `notes/`: Meta-Aging; `documentation/` and `explainers/`: 14616-DM; `notes.md`: neuroimmune |
| Architecture/context | `architecture/`; `workflow-architecture/`; `overview/`; `scientific-context.md`; architecture/ruling files inside `reasoning/` | Meta-Aging; 14616-DM; Meta-Aging; 14761-DM, 14782-DM, neuroimmune; 14782-DM |
| Scratch | project-root `_scratch/`; `03_results/_scratch/`; `02_analysis/scripts/_scratch/`; `/tmp/...`; `/tmp/claude-.../scratchpad/` | Root variant absent everywhere; results variant exists in all five project scopes; script variant only 14616-DM; generic `/tmp` in Meta-Aging, 14616-DM, 14761-DM, 14782-DM; exact Claude pattern in 14761-DM and 14782-DM |

## Signal-to-noise estimates

These are purposive 3–6-file samples chosen to cover different categories, not random corpus estimates.

| Project | Sample assessment | Rough useful-signal proportion |
|---|---|---:|
| Meta-Aging | 3 durable decisions/syntheses, 1 handoff, 1 plan, 1 raw transcript artifact | about 5/6, 83% |
| 14616-DM | 3 durable decisions/explainers, 1 handoff journal, 1 plan, 1 copied tutorial artifact | about 5/6, 83% |
| 14761-DM | 2 durable records, 2 useful handoff/current-state records, 1 partly superseded plan, 1 scratch-derived probe/spec | about 4/6, 67% |
| 14782-DM | 1 durable decision record, 1 handoff, 1 completed operational report, 1 explicitly superseded plan, 2 generated artifacts | 2/6 lasting signal; 3/6 if the operational report counts |
| neuroimmune-receptor-atlas | 1 durable scientific context, 1 useful layout scaffold, 1 scratch note | 2/3 useful; only 1/3 scientific memory |

The structural evidence is harsher than the samples for 14782-DM: 210 cache JSONs, 54 non-Markdown validation artifacts, a model checkpoint, logs, Parquet data, Python bytecode, and Git plumbing dominate the tree.

## What could not be determined

- Exact authorship of regular files under harness directories. Symlinks are clearly mounts, while local settings, locks, hooks, Gemini configuration, and the regular `marimo-pair` skill are harness-specific accumulation; their human-versus-agent provenance is not recorded.
- A reliable whole-corpus semantic proportion. The trees were intentionally sampled rather than read exhaustively.
- Whether every old plan is superseded. I treated only explicitly labeled material as certainly stale and marked other cases as inferred.
- Actual historical `/tmp` usage beyond recorded textual traces. Absence of a reference does not prove an agent never wrote there.

## Observations for the designer

- The largest contradiction has already produced a field workaround: 14761-DM and 14782-DM keep `docs/_internal` ignored by the parent while tracking it in a nested repository. Parent-level audits report zero even though most internal content is independently versioned.
- Handoff information currently has at least seven structural/name variants, and both `handoffs/` directories are empty while real handoffs accumulate in `sessions/` and plan-local `_handoff/`.
- 14782-DM demonstrates that “internal memory” currently mixes decisions with executable code, caches, logs, binary data, a model checkpoint, and repository plumbing.
- Every project has a results-level scratch directory, yet exact external Claude scratchpad traces still occur in the two busiest projects.
