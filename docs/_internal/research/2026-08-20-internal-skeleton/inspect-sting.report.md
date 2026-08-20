# `docs/_internal/` field survey

Counts use `find ... -type f`; tracked counts use `git ls-files docs/_internal`. Each assigned path is its own Git repository.

| Project | On-disk files | Tracked | Immediate subdirectories, recursive file counts | Harness dirs present | Root `_scratch/` |
|---|---:|---:|---|---|---|
| `STING-JR` | 61 | 0 | `_superceded` (3), `plans` (19), `reasoning` (23), `research` (1), `scratch` (7), `sessions` (4), `workflows` (1) | `.claude`, `.agents` | No |
| `human_pbmc_febrile` | 14 | 0 | `handoffs` (1), `plans` (3), `reasoning` (3), `reports` (2), `research` (2), `sessions` (1) | `.claude`, `.agents` | No |
| `human_ra_synovium` | 18 | 0 | `handoffs` (1), `plans` (2), `reasoning` (7), `reports` (2), `research` (3), `sessions` (1) | `.claude`, `.agents` | No |
| `human_treg_arthritis` | 28 | 0 | `handoffs` (1), `plans` (1), `reasoning` (17), `reports` (2), `research` (2), `sessions` (3) | `.claude`, `.agents` | No |
| `mouse_anchor` | 159 | 0 | `_scratch` (48), `handoffs` (1), `inbox` (17), `plans` (23), `reasoning` (46), `reply-package` (6), `reports` (2), `research` (10), `sessions` (4) | `.claude`, `.agents` | Yes, empty |
| `sting_positive_control` | 18 | 0 | `handoffs` (1), `plans` (1), `reasoning` (8), `reports` (3), `research` (2), `sessions` (1) | `.claude`, `.agents` | No |

All six also contain tracked `AGENTS.md` and `CLAUDE.md`; every `CLAUDE.md` is exactly `@AGENTS.md`. None contains `.codex/` or `.gemini/`.

## Layout, tracking, and project differences

Common child-project scaffold:

- The five child projects have `README.md` and `scientific-context.md` directly under `docs/_internal/`.
- Each has five `.gitkeep` files: one apiece in `handoffs`, `plans`, `reasoning`, `reports`, and `research`.
- Each also has `reasoning/README.md`, `research/README.md`, and `sessions/README.md`.
- Consequently, 10 files per child are scaffold/context before counting any substantive plan, decision, report, or handoff.
- `handoffs/` is empty except for `.gitkeep` in every child project. Where actual handoffs exist, agents wrote them under `sessions/`.

Depth and exceptions:

- `STING-JR`: three direct files—`STATUS.md`, `runbook.md`, and `reactive-notebooks-quickstart.md`. Most reasoning and session records sit directly in their subdirectory. Plans use dated directories such as `plans/2026-07-26-scoring-foundation-reset/NN_*.md`; the seven GMT scratch files likewise sit under one dated directory. Maximum file depth is `_internal/<subdir>/<dated-dir>/<file>`.
- `human_pbmc_febrile`: both substantive plans are nested under dated directories. There are no actual handoff or session records.
- `human_ra_synovium`: its single plan is nested under a dated directory. Reasoning, report, and research records are direct children. There are no actual handoff or session records.
- `human_treg_arthritis`: all files are direct children of their category directories. It has two actual records under `sessions/`; `plans/` contains only `.gitkeep`.
- `mouse_anchor`: the deepest and most heterogeneous tree. Plans, research, reasoning, inbox material, reply packages, and scratch previews all use another directory layer. Its 159 files include 49 PNGs, two PPTX files, one PDF, two Python scripts, and 99 Markdown files. The tree is about 15.2 MB. There is both an internal `_scratch/` containing 48 files and an empty project-root `_scratch/`.
- `sting_positive_control`: all categorized files are direct. `plans/`, `handoffs/`, and `sessions/` contain scaffolding only.

Git ignore provenance:

- Every project has zero tracked `docs/_internal` files.
- In the five child repositories, `git check-ignore -v` resolves to `docs/_internal/` inside a managed `# BEGIN SCIAGENT:GITIGNORE` block:
  - PBMC line 71
  - RA line 80
  - Treg line 124
  - Mouse line 107
  - SAVI line 77
- The umbrella differs: ignore provenance is a hand-written `.gitignore:69:_internal/` rule, preceded by “Internal / private notes.” It is outside a managed block and ignores any directory named `_internal`.

## Scratch behaviour

- `STING-JR/docs/_internal/scratch/` holds seven downloaded GMT gene-set payloads.
- `mouse_anchor/docs/_internal/_scratch/` holds `notes.md` plus 47 versioned preview PNGs in directories named `previews`, `previews_audit`, `previews_adv`, `previews_fix`, `previews_impl`, and `previews_tript`.
- Only `mouse_anchor` has project-root `_scratch/`, and it is empty.

References to external temporary storage:

- Four umbrella internal files mention scratch or `/tmp`:
  - a plan instructs the orchestrator to compose prompts “under the scratchpad”;
  - `reactive-notebooks-quickstart.md` and a session handoff route server logs to `/tmp`;
  - a reasoning record says diagnostics ran under `/tmp/aliasdiag/`.
- `mouse_anchor` has a forensic reasoning record whose claims depend on crops at `/tmp/forensics/seg_0..5.png`.
- The exact observed Claude scratchpad pattern appears in the umbrella’s `.claude/settings.local.json`: three permission-history entries reference `/tmp/claude-788715489/-workspaces-STING-JR/<uuid>/scratchpad/`, including a probe output and a `fix3/launch.sh`.
- No `/tmp`, `claude-<pid>`, or `scratchpad` reference was found in any surveyed project’s own `02_analysis/` tree, including its notebooks.
- I did not inspect the referenced `/tmp` artifacts because they are outside the assigned paths. Their survival and contents are therefore unknown.

## Harness footprints

The five child projects primarily contain managed mounts:

- PBMC, RA, and Treg: 51 symlinks in `.claude/` and 51 in `.agents/`.
- SAVI: 52 symlinks in each, reflecting one additional mounted skill.
- Mouse: 51 symlinks in each.
- `.agents/` contains no regular files in any child.
- PBMC, RA, Treg, and SAVI each have five regular `.claude/` files: two hooks, hook documentation, `settings.json`, and `statusline.sh`. Four have identical hashes across projects; the README differs by embedded project name. Their content identifies them as project-scoped SCIAGENT guardrails, so they appear provisioned rather than accumulated session memory.
- Mouse has no regular `.claude/` files.

The umbrella is the exception:

- `.claude/`: two regular files and one symlink. `settings.local.json` is harness-local permission history and contains the external scratchpad traces. `skills/delegate-cli/SKILL.md` identifies SciAgent-toolkit as its author.
- `.agents/`: eight regular files comprising a copied `marimo-pair` skill, with no symlinks.
- Attribution of the copied `marimo-pair` directory to a human, agent, or installer cannot be proven from its contents alone.

## Consolidated naming variants

| Concept | Every spelling/form observed | Projects |
|---|---|---|
| Session continuity / handoff | `sessions/` | All six |
|  | `handoffs/` | All five child projects; empty except `.gitkeep` |
| Durable reasoning / review / completion record | `reasoning/` | All six |
|  | `reports/` | All five child projects |
| Plans | `plans/` | All six |
| Research/reference synthesis | `research/` | All six |
| Scratch | `scratch/` | Umbrella |
|  | `_scratch/` under `docs/_internal` | Mouse |
|  | project-root `_scratch/` | Mouse, empty |
| Archive / lifecycle state | `_superceded/` | Umbrella |
|  | `*.SUPERSEDED.*` | Umbrella, inside the misspelled directory |
|  | `*.EXECUTED.md` | Umbrella |
| Project-wide status/context | `STATUS.md`, `runbook.md` | Umbrella |
|  | `README.md`, `scientific-context.md` | All five child projects |
| Collaborator intake/output | `inbox/`, `reply-package/` | Mouse |
| Workflow implementation | `workflows/` | Umbrella |

No singular `handoff/`, singular `session/`, or `decisions/` directory was observed. The strongest drift is semantic rather than merely grammatical: the scaffold creates both `handoffs/` and `sessions/`, but actual handoffs are written into `sessions/`.

The archive spelling also drifts within one project: the directory is `_superceded` while files inside it use the conventional spelling `SUPERSEDED`.

## Signal-to-noise estimates

These are bounded estimates from 3–6 content samples per project plus filename/category counts. “Signal” means a substantive durable or historically useful record; it does not guarantee that every claim remains current.

### `STING-JR` — sampled 6; estimated signal 65–75%

- `STATUS.md`: durable canonical status record, although later August records make freshness uncertain.
- scoring-foundation `00_INDEX.md`: detailed plan that explicitly supersedes earlier plans.
- `reasoning/2026-07-28_narrative-snapshot.md`: stale-by-design snapshot; it says it will be wrong within a week.
- `sessions/2026-07-28-hsr-decomp-handoff.md`: genuine session handoff.
- `_superceded/...EXECUTED.md`: explicitly stale executed plan.
- one GMT file under `scratch/`: scratch/reference payload rather than project memory.

Ten of 61 files are explicitly scratch or archived. Additional older plans remain in `plans/` after a newer index says they are superseded, lowering confidence that every non-archived file is live.

### `human_pbmc_febrile` — sampled 5; estimated signal 25–35%

- Two febrile plan indexes: plans, both subsequently superseded by the umbrella foundation reset.
- Kawasaki QC provenance audit: durable decision/provenance record.
- Febrile foundation report: durable completion/blocker report.
- `sessions/README.md`: scaffold only.

Only about four to five of 14 files are substantive project records, and at least one plan pair overlaps stale scope. There are no actual handoffs.

### `human_ra_synovium` — sampled 6; estimated signal 40–50%

- RA replication index: plan with a prominent partial-supersession banner.
- Ingest/QC note: durable method decision.
- Pseudobulk note: useful decision history, explicitly partially superseded.
- Seam review: durable defect/review record.
- T-cell annotation research: substantive research synthesis.
- stages 05–08 report: durable completion report.

Eight to nine of 18 files are substantive; the remainder is mostly scaffold. Several substantive records preserve corrected history rather than current state.

### `human_treg_arthritis` — sampled 6; estimated signal 55–65%

- Phase-1 scope/QC reasoning: durable historical decision record, with old engine-era numbers.
- Stale explorer-column note: durable unresolved hazard.
- Projection regeneration note: durable contract/provenance record.
- JIA sign-off report: durable completion record.
- Phase-1 go/no-go file: genuine handoff, now a dated numerical snapshot.
- CX2 coupling file: session/completion handoff.

This has the best substantive-to-scaffold ratio among the smaller child projects: roughly 18–19 of 28 files contain project-specific material, although several early records have become historical.

### `mouse_anchor` — sampled 6; estimated signal 50–60%

- `_scratch/notes.md`: raw conversational/probe artifact.
- reply synthesis: substantive research/collaborator audit.
- figure-sweep index: executed plan.
- CoReSh `RUN.md`: strong durable chain-of-custody record.
- claim-evidence memo: durable internal communication draft.
- DE handoff: genuine session handoff.

At least 48 of 159 files are explicitly scratch previews. Six more are duplicated dated reply-package builds, and several inbox files are source attachments rather than memory. The remaining reasoning, plans, research, reports, and handoffs are often high-signal, but the directory also functions as an artifact warehouse.

### `sting_positive_control` — sampled 6; estimated signal 40–50%

- Build/QC decisions: durable historical decision record.
- Seam migration review: useful defect audit, later superseded in outcome by repairs.
- share-export note: durable implementation decision.
- compute-defects report: durable repair/completion record.
- boundary-panel report: durable sign-off record.
- `sessions/README.md`: scaffold only.

Eight to nine of 18 files are substantive. There are no actual plans, handoffs, or session records despite all three scaffold locations existing.

## Observations for the designer

- Across these six repositories, 298 `docs/_internal` files exist on disk and zero are tracked.
- Empty scaffold dominates the smallest projects: 10 of 14 PBMC files, 10 of 18 RA files, and 10 of 18 SAVI files are scaffold/context rather than dated operational records.
- `handoffs/` and `sessions/` coexist, but agents consistently use `sessions/` for real handoffs.
- The umbrella and mouse projects demonstrate two different landfill modes: retained superseded plans in the umbrella, and binary previews/versioned reply packages in mouse.
- The most disciplined record sampled was mouse’s CoReSh `RUN.md`, which connects prompt, raw tool output, validation, derived artifact, and downstream consumer.
- The exact off-project Claude scratchpad pattern is confirmed by harness permission history even though no surviving scratchpad was inspected.
- Exact authorship, present-day correctness, and full-tree signal proportions cannot be established from a bounded sample. Explicit supersession banners helped; unbannered older records remain ambiguous.