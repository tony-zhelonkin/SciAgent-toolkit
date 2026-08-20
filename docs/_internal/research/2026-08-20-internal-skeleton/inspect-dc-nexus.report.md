# Field survey: project-internal memory layout

Counts use `find … -type f`; symlinks are excluded. All six paths are independent Git worktrees.

| Project | `docs/_internal` files | Tracked | Immediate subdirectories | Harness dirs present | Root `_scratch/` |
|---|---:|---:|---|---|---|
| `DC-nexus` | 4,169 | 2 | `grant`, `overview`, `research` | `.claude` | No |
| `DC_Dictionary` | 0 — absent | 0 | — | `.claude`, `.gemini` | No |
| `DC_hum_verse` | 0 — absent | 0 | — | `.claude`, `.agents`, `.gemini` | No |
| `DC_mouse_cancer` | 188 | 0 | `.claude`, `.ref`, `handoffs`, `personal-notes`, `plans`, `reasoning`, `reports`, `research`, `sessions` | `.claude`, `.agents` | No |
| `14839-DM-cGAS` | 133 | 0 | `handoffs`, `pathway-explorer`, `plans`, `reasoning`, `reports`, `research`, `sessions` | `.claude` | No |
| `13036-DM_DMlab_summer_2025` | 0 — absent | 0 | — | `.claude`, `.gemini` | No |

Common pattern: every project has regular `AGENTS.md` and `CLAUDE.md`; none has `.codex/`. No project has the prescribed root `_scratch/`.

## Per-project findings

### `DC-nexus`

Layout:

| Location | Files | Depth |
|---|---:|---|
| Root of `docs/_internal` | 2 | `setup-notes.md`, `nvim-topology-setup.md` |
| `overview/` | 1 | Direct |
| `research/` | 2 | Direct |
| `grant/` | 4,164 | Nested as deep as 15 levels |

The directory is 568 MB. The exceptional depth and size come largely from `grant/local/cri/raw/.venv/`: 3,813 files, or 91.5% of the entire internal tree. The tree also contains 4 virtual-environment symlinks, compiled libraries, Python packages, PDFs, figures, proposal versions, and scripts. Only 356 files remain after excluding `.venv/`.

The two tracked files are:

- `research/Cell-to-Hum_similarity_REVIEW.md`
- `research/discussion.md`

`git check-ignore -v` identifies `.gitignore:48:_internal/` for an untracked probe. This is a handwritten generic rule under “Internal / private notes,” not a managed block.

Six sampled records:

- `overview/DC-NEXUS_overview.md` — durable project snapshot/seeding prompt, dated 2026-07-01.
- `setup-notes.md` — operational reference with unresolved Git/container state; potentially stale.
- `research/discussion.md` — raw Claude conversation/research scratch; tracked.
- `grant/.../2026-07-21-...-kickoff-plan.md` — extensive early grant plan later followed by multiple proposal revisions.
- `grant/local/cri-proposal/A3/NEXT_SESSION_PROMPT.md` — session handoff, subsequently overtaken by later B-series drafts.
- `grant/local/cri-proposal/B3/revision-notes.md` — durable decision/provenance record with explicit open verification items.

Harness footprint: `.claude/` contains three regular files (`settings.json`, `settings.local.json`, `statusline.sh`) and no mounted symlinks.

No `/tmp/` or Claude scratchpad reference was found in the searched internal, analysis, or notebook content.

### `DC_Dictionary`

There is no `docs/_internal/`, and `git check-ignore -v docs/_internal/_probe` returns no rule.

Agent memory instead accumulated at the repository root and in hidden directories:

- Root: tracked `plan.md` (1,417 lines), `done.md` (1,001 lines), and `handoff_20260220_025647.md`.
- `.handoff_archive/`: 54 files.
- `.plans_archive/`: 22 files.
- `.research_archive/`: 4 files.
- `.context/`: 6 files, including an untracked Claude compaction report.
- `.archive/`: 49 superseded scripts, documents, and plots.

Five samples:

- `plan.md` — durable scientific plan, but its title says v9.5 while metadata says v9.8.
- `handoff_20260220_025647.md` — current-at-time session handoff.
- `done.md` — durable but very large accumulated completion/history record.
- `.context/CLAUDE_compaction_summary.md` — stale harness-maintenance artifact.
- `.handoff_archive/CONSOLIDATION_RECORD.md` — historical provenance describing consolidation of nine earlier handoffs.

`.claude/` mixes 21 toolkit symlinks with 10 regular files. Regular authored/custom material includes four agent definitions, two command files, one output style, and two archived skills. Three command/style files are tracked; the agent definitions and archived skills are untracked. `.gemini/` contains only `settings.json`.

One archived handoff references `/tmp/test_dir/nested/path` as a generic filesystem test. No `/tmp/claude-…/scratchpad/` trace was found.

### `DC_hum_verse`

There is no `docs/_internal/`, and Git reports no ignore rule for it.

Memory is stored as tracked root files:

- `context.md`
- `plan.md`
- `tasks.md`
- `research_notes.md`
- `handoff_20260306_191330.md`

There are also six files under `.handoff_archive/`.

Five samples:

- `context.md` — durable scientific context, but its checklist still describes the project as largely unstarted.
- `plan.md` — detailed plan last updated 2026-02-18; later handoff evidence supersedes much of its status.
- `tasks.md` — stale/misapplied template: it identifies the project as `DC_cancer` and uses incompatible paths.
- `handoff_20260306_191330.md` — useful current-at-time session handoff with completed remediation and next steps.
- `research_notes.md` — substantive biological review/interpretation record.

Harness-specific accumulation is pronounced. `.claude/` has 49 toolkit symlinks plus 10 untracked regular files, including `.claude/plan.md` and six backup agent/command documents under `.claude/.bak/`. `.agents/` contains 48 symlinks and no regular files. `.gemini/` contains settings only.

No `/tmp/` or scratchpad reference was found.

### `DC_mouse_cancer`

This is the closest project to a populated common skeleton, but it contains several distinct payload types:

| Subdirectory | Files | Shape |
|---|---:|---|
| `.claude` | 1 | Direct backup settings file |
| `.ref` | 15 | All under `.ref/_scratch/{base,viz-llm}`; depth 3 |
| `handoffs` | 2 | One handoff plus `.gitkeep` |
| `personal-notes` | 1 | Direct |
| `plans` | 1 | `.gitkeep` only |
| `reasoning` | 153 | 11 direct; 142 under `2026-06-04-annotation/`, reaching depth 7 |
| `reports` | 1 | `.gitkeep` only |
| `research` | 3 | Direct |
| `sessions` | 9 | Direct |
| Internal root | 2 | `README.md`, `scientific-context.md` |

The 142-file annotation subtree contains only 8 Markdown files and 134 JSON/CSV/log/text artifacts. `.ref/_scratch/` contains 15 scripts or related scratch files.

`git check-ignore -v` reports `.gitignore:70:docs/_internal/`, inside the older `# BEGIN SCIAGENT:GITIGNORE` managed block. Nothing is tracked.

Six samples:

- `handoffs/2026-06-03_rerun-kickoff.md` — session handoff.
- `sessions/2026-06-05_annotation-viz-reconciliation-palette.md` — detailed session handoff.
- `reasoning/2026-06-05_annotation-label-normalization-architecture.md` — approved architectural decision record.
- `reasoning/2026-06-04-annotation/README.md` — durable run-provenance index.
- `personal-notes/2026-06-03_approach-to-annotating-myeloid-cells.md` — generic reference note with weak project specificity.
- `reasoning/2026-06-03_gpu-enablement.md` — useful diagnosis at creation time, now likely historical after subsequent GPU reruns.

The reasoning documents use frontmatter with `kind`, `status`, tags, and related-file links—a comparatively explicit convention.

At the root, `.agents/` contains only 48 toolkit symlinks. `.claude/` contains those mounts plus `settings.local.json`. A separate stale `.claude/settings.json.bak` sits inside `docs/_internal/`.

There is no root `_scratch/`, but `docs/_internal/.ref/_scratch/` exists. No `/tmp/` reference was found.

### `14839-DM-cGAS`

Layout:

| Subdirectory | Files | Shape |
|---|---:|---|
| `plans` | 83 | One `.gitkeep`, then `plans/<date-slug>/NN_*.md`; depth 2 |
| `research` | 28 | Two direct files plus five dated research bundles; depth 2 |
| `reasoning` | 16 | Direct |
| `pathway-explorer` | 1 | Direct |
| `handoffs` | 1 | `.gitkeep` only |
| `reports` | 1 | `.gitkeep` only |
| `sessions` | 1 | `README.md` only |
| Internal root | 2 | `README.md`, `scientific-context.md` |

Of 133 files, 127 are Markdown. `git check-ignore -v` reports `.gitignore:98:docs/_internal/` in a managed `SCIAGENT:GITIGNORE` block. Nothing is tracked.

Six samples:

- `plans/2026-06-17-init/00_INDEX.md` — detailed phase plan.
- `plans/2026-07-02-coupling-continuous-dependence/06_REVIEW_FINDINGS.md` — durable adversarial-review and correction record.
- `plans/2026-07-02-coupling-pipeline/00_INDEX.md` — consolidated execution plan that explicitly supersedes two retained source plans.
- `reasoning/2026-06-17_01_init-plan-architecture.md` — plan-architecture rationale.
- `research/.../2026-06-17_01_repo-conventions.md` — reusable codebase contract/research note.
- `pathway-explorer/pathway-explorer-input-spec.md` — durable artifact/schema contract.

The explicit supersession declaration is unusually good: 19 files in the two source plan bundles remain as provenance, while the consolidated index states which plan to execute.

`.claude/` contains only `settings.local.json`; no harness-authored documents were found there.

There is no root `_scratch/`, but both `03_results/_scratch/` and root `scratch_png/` exist. No `/tmp/` reference was found.

### `13036-DM_DMlab_summer_2025`

There is no `docs/_internal/`, and Git reports no ignore rule for it.

Tracked memory lives at the root as `context.md`, `plan.md`, `tasks.md`, `done.md`, and `handoff-2026-02-11.md`.

Five samples:

- `context.md` — untouched placeholder template; scratch/stale artifact.
- `plan.md` — substantive scientific and statistical plan.
- `tasks.md` — detailed execution plan, but it says high-priority work is unstarted while the later handoff says those analyses are complete.
- `handoff-2026-02-11.md` — consolidated session/status handoff.
- `done.md` — durable completion and troubleshooting record.

This repository uses capitalized stage names (`02_Analysis`, `03_Results`), another convention drift outside the internal skeleton.

`.claude/` mixes 18 toolkit symlinks with `settings.local.json`; `.gemini/` contains settings only. No harness-specific authored documents or `/tmp/` references were found.

## Consolidated naming variants

“Scaffold” means the directory contains only `.gitkeep` or a convention README.

| Concept | Spellings/locations observed | Projects |
|---|---|---|
| Session continuity | `docs/_internal/handoffs/`; `docs/_internal/sessions/`; root `handoff_<timestamp>.md`; root `handoff-YYYY-MM-DD.md`; `.handoff_archive/`; embedded `NEXT_SESSION_PROMPT.md` / `NEXT-SESSION.md` | `DC_mouse_cancer` uses both active `handoffs` and `sessions`; `14839` has scaffold versions of both; `DC_Dictionary`, `DC_hum_verse`, `13036`; embedded prompt form in `DC-nexus` |
| Plans | `docs/_internal/plans/<date-slug>/NN_*.md`; scaffold `plans/`; root `plan.md`; `.plans_archive/`; `*-plan.md` under `reasoning/`; grant-local kickoff plans | Active nested plans in `14839`; scaffold in `DC_mouse_cancer`; root plans in `DC_Dictionary`, `DC_hum_verse`, `13036`; archive in `DC_Dictionary`; reasoning-plan names in `DC_mouse_cancer` and `14839`; grant form in `DC-nexus` |
| Decisions/rationale | `reasoning/`; `reports/`; `personal-notes/`; review and revision notes embedded inside plan/grant bundles | Active `reasoning` in `DC_mouse_cancer`, `14839`; `reports` is scaffold-only in both; `personal-notes` in `DC_mouse_cancer`; embedded review/revision records in `14839`, `DC-nexus` |
| Research | `docs/_internal/research/`; nested `grant/.../research/lit/{raw,synthesis}`; root `research_notes.md`; `.research_archive/` | `DC-nexus`, `DC_mouse_cancer`, `14839`; root form in `DC_hum_verse`; archive form in `DC_Dictionary` |
| Project context | `overview/`; internal `scientific-context.md`; internal `README.md`; root `context.md`; `.context/`; root `done.md` | `DC-nexus`; `DC_mouse_cancer`, `14839`; `DC_hum_verse`, `13036`; `.context` and `done.md` in `DC_Dictionary`, with `done.md` also in `13036` |
| Scratch/transient work | Expected root `_scratch/` absent everywhere; `docs/_internal/.ref/_scratch/`; `03_results/_scratch/`; root `scratch_png/`; result `tmp/` directories; `.claude/.bak/` | Internal scratch in `DC_mouse_cancer`; results/root scratch variants in `14839`; result `tmp` in `DC_Dictionary`; harness backups in `DC_hum_verse` |
| Archives/supersession | `.archive/`; `.handoff_archive/`; `.plans_archive/`; `.research_archive/`; `.claude/.bak/`; retained date-slug source plans explicitly superseded by a consolidated plan | `DC_Dictionary`; `DC_hum_verse`; `14839` |

The survey found plural `handoffs/` and `sessions/`, but no singular `handoff/` or `session/` directory. The two plural forms coexist in the same projects, suggesting partially distinct intended meanings rather than simple spelling alternatives.

## Signal-to-noise estimates

These are rough estimates from the specified samples plus file-type/path counts, not exhaustive content judgments.

| Project | Estimate |
|---|---|
| `DC-nexus` | By file count, about 92% virtual-environment dependency material and 8% all other content. Among six sampled Markdown files, roughly one-third were durable records and two-thirds were chat, transient setup, plan, or superseded handoff material. |
| `DC_Dictionary` | No internal tree to score. In the five-file legacy-memory sample, about 40% was live/current-at-time planning or handoff material, 20% a durable historical summary, and 40% stale harness/consolidation artifacts. Hidden archive volume is high. |
| `DC_hum_verse` | About 40% useful/current-at-time signal and 60% stale or superseded state in the five-file sample; the wrong-project `tasks.md` is the clearest noise item. |
| `DC_mouse_cancer` | About 80% of files are machine run traces or `_scratch` content; about 20% are narrative documents/scaffolds. The sampled narrative subset itself is mostly high-signal, especially the frontmatter-bearing reasoning and handoff records. |
| `14839-DM-cGAS` | Nearly all files are human-readable Markdown. Roughly 60–70% appears durable/current planning, research, reasoning, or contracts; roughly 25–35% is historical or explicitly superseded plan provenance; the remainder is scaffolding. |
| `13036-DM_DMlab_summer_2025` | Three of five samples are substantive durable records; two are stale/template material: approximately 60% signal, 40% noise or superseded status. |

## Observations for the designer

- Location drift is larger than directory-name drift: three of six projects have no `docs/_internal/` and instead use tracked root files plus hidden archives.
- A shared directory name does not constrain payload. The surveyed trees contain decisions, handoffs, raw chat, machine traces, copied source, PDFs, figures, virtual environments, and scratch scripts.
- `handoffs/` and `sessions/` coexist without an evident cross-project semantic boundary.
- Two practices preserved useful provenance particularly well: status/frontmatter and related links in `DC_mouse_cancer`, and explicit “supersedes/retained as provenance” declarations in `14839`.
- The managed ignore marker observed is `SCIAGENT:GITIGNORE`, while `DC-nexus` uses a handwritten `_internal/` rule. Three projects with no internal directory have no matching ignore rule at all.
- Harness-neutrality is currently strongest in `14839` and `DC_mouse_cancer`, where harness directories are almost entirely mounts/config. `DC_Dictionary` and especially `DC_hum_verse` contain custom or backup documents inside `.claude/`.

What could not be determined:

- Filesystem inspection cannot reliably distinguish human authorship from agent authorship; “agent-authored” above is inferred from filenames, content, and placement.
- A missing `/tmp/claude-…/scratchpad/` reference does not prove such paths were never used; temporary external state may leave no repository trace.
- Currentness cannot be established for every plan without correlating all documents against implementation and Git history. Explicit status and supersession statements were used where present.
