# `docs/_internal/` field inspection

| Project | On-disk files | Tracked | Immediate subdirectories | Harness dirs present | `_scratch/` |
|---|---:|---:|---|---|---|
| `12868-EH` | 0 — absent | 0 | — | `.claude`, `.gemini` | Absent |
| `AdaW-eWAT-WL-bulkRNAseq` | 0 — absent | 0 | — | — | Absent |
| `Dakota_Black_NPC_DRP1` | 0 — absent | 0 | — | `.claude`, `.gemini` | Absent |
| `Gama_Vivian_DRP1_bulkRNAseq` | 0 — absent | 0 | — | `.claude`, `.gemini` | Absent |
| `13403-YD_Christina` | 100 | 0 | `app`, `contracts`, `handoffs`, `issues`, `notes`, `plans`, `reasoning`, `states`, `tools` | `.claude`, `.agents`, `.gemini` | Absent |
| `JBader_scHFD` | 0 — absent | 0 | — | `.claude`, `.gemini` | Absent |
| `scbio-docker` | 0 — absent | 0 | — | `.claude`, `.gemini` | Absent |

All seven paths are Git repository roots. `.codex/` is absent from every project.

## Project details

Six projects have no project-root `docs/_internal/` directory. For all six, `git check-ignore -v docs/_internal/probe` found no applicable ignore rule. Consequently, there were no internal files to sample and no internal naming vocabulary to extract from those projects.

### `13403-YD_Christina`

Layout:

| Subdirectory | Files | Placement |
|---|---:|---|
| `app/` | 15 | `app/<topic>/`; one topic also has `design/` beneath it |
| `contracts/` | 1 | Direct |
| `handoffs/` | 2 | Direct |
| `issues/` | 2 | Direct |
| `notes/` | 2 | Direct |
| `plans/` | 47 | All under `plans/yd-finalization/` |
| `reasoning/` | 22 | Direct |
| `states/` | 3 | All under `states/2026-06-18-cell-viewer/` |
| `tools/` | 6 | All under `tools/hdWGCNA/` |

There are no files directly in `docs/_internal/`.

All 100 files are untracked. `git check-ignore -v` resolves to the handwritten rule `.gitignore:5:docs/_internal/`; lines 4–5 contain duplicate forms:

```text
docs/_internal/*
docs/_internal/
```

There is no `SCIO:GITIGNORE` marker in that `.gitignore`.

Six sampled files:

- [handoff_20260203_032339.md](</data2/users/JCRLab/13403-YD-analysis/13403-YD_Christina/docs/_internal/handoffs/handoff_20260203_032339.md>) — session handoff with executable next steps. Later June records make it historical.
- [_campaign.md](</data2/users/JCRLab/13403-YD-analysis/13403-YD_Christina/docs/_internal/plans/yd-finalization/_campaign.md>) — durable orchestration ledger recording locked decisions, gates, dependencies, and completion status.
- [signature_reconciliation_SYNTHESIS_v3_LOCKED.md](</data2/users/JCRLab/13403-YD-analysis/13403-YD_Christina/docs/_internal/reasoning/signature_reconciliation_SYNTHESIS_v3_LOCKED.md>) — durable scientific decision record. It explicitly supersedes `reconciliation_SYNTHESIS_v2.md`, which remains present.
- [session-store-and-snapshots-state.md](</data2/users/JCRLab/13403-YD-analysis/13403-YD_Christina/docs/_internal/states/2026-06-18-cell-viewer/session-store-and-snapshots-state.md>) — detailed point-in-time artifact and implementation-state snapshot.
- [00_SYNTHESIS.md](</data2/users/JCRLab/13403-YD-analysis/13403-YD_Christina/docs/_internal/app/live-compute-contrast/00_SYNTHESIS.md>) — architectural synthesis containing user decisions and unresolved gates.
- [intro.md](</data2/users/JCRLab/13403-YD-analysis/13403-YD_Christina/docs/_internal/tools/hdWGCNA/intro.md>) — generic package/install reference material with relatively little project-specific memory.

An alternate memory convention exists beside this tree:

- `.handoff_archive/`: 5 files, all untracked.
- `docs/2026-05-25/`: 52 files, all tracked.

The local `.claude/commands/kickoff.md` routes agents to the tracked dated bundle for handoffs, decisions, and phase plans. This tracked bundle appears to be the project’s more reliable continuation record.

### Harness footprints

| Project | Root instruction files | Local regular harness content | Toolkit-style symlinks |
|---|---|---|---:|
| `12868-EH` | `CLAUDE.md` | Eight tracked `.claude` agent/skill documents plus local settings; `.gemini/settings.json` | 41 |
| `AdaW-eWAT-WL-bulkRNAseq` | `AGENTS.md`, `CLAUDE.md` | None | 0 |
| `Dakota_Black_NPC_DRP1` | `AGENTS.md`, `CLAUDE.md` | Claude/Gemini settings only | 30 |
| `Gama_Vivian_DRP1_bulkRNAseq` | `AGENTS.md`, `CLAUDE.md` | Claude settings, `scheduled_tasks.lock`, Gemini settings | 0 |
| `13403-YD_Christina` | `AGENTS.md`, `CLAUDE.md` | Untracked project-specific `kickoff.md` plus settings | 163 across `.claude` and `.agents` |
| `JBader_scHFD` | `AGENTS.md` | Six untracked deprecated Claude skill documents plus settings | 6 |
| `scbio-docker` | `AGENTS.md`, `CLAUDE.md` | Two untracked Claude agent stubs plus settings | 0 |

The regular content in `12868-EH`, `13403-YD_Christina`, `JBader_scHFD`, and `scbio-docker` is genuine harness-specific accumulation rather than mounted toolkit links. File provenance alone does not establish whether a human or an agent originally authored each item.

## Naming-variant table

Because only one project has `docs/_internal/`, this assignment cannot quantify cross-project spelling frequency. The following are all variants observed in the populated tree.

| Concept | Spellings or placements observed | Project |
|---|---|---|
| Handoff | `handoffs/`; filenames use both `HANDOFF_YYYY-MM-DD.md` and `handoff_YYYYMMDD_HHMMSS.md` | `13403-YD_Christina` |
| Session/state | `states/<date-slug>/*-state.md`; no `session/` or `sessions/` directory | `13403-YD_Christina` |
| Plan | `plans/<campaign>/`; also `app/<topic>/PLAN.md`; plan filenames include `phase-NN`, `impl-phase-NN`, `remediation-NN`, and underscore-prefixed campaign metadata | `13403-YD_Christina` |
| Reasoning/decision | `reasoning/`; decision-bearing records also occur in `plans/.../_campaign.md` and `app/.../00_SYNTHESIS.md`; no `decisions/` directory | `13403-YD_Christina` |
| Report/synthesis | `app/<topic>/00_SYNTHESIS.md` and `reasoning/*SYNTHESIS*.md`; no `reports/` directory | `13403-YD_Christina` |
| Application design | `app/`, with optional nested `design/` | `13403-YD_Christina` |
| Contract | `contracts/` | `13403-YD_Christina` |
| Issue/troubleshooting | `issues/` | `13403-YD_Christina` |
| Notes | `notes/` | `13403-YD_Christina` |
| Tool reference | `tools/<tool>/` | `13403-YD_Christina` |

Thus the observed drift is primarily within one project: inconsistent handoff filenames and overlapping placement of plans, decisions, designs, and syntheses. The requested cross-project contrasts such as `handoff/` versus `handoffs/` cannot be measured from this seven-project slice.

## Signal-to-noise estimates

| Project | Estimate |
|---|---|
| Six projects without `docs/_internal/` | Not applicable; no files |
| `13403-YD_Christina` | Sample: approximately 50% high-signal durable records, 33% historical handoff/state material, and 17% generic low-signal reference material |

The exact corpus composition by directory is 47% plans, 22% reasoning, 15% application/design reports, 5% handoffs and state snapshots, and 11% contracts/issues/notes/tool references.

The six-file sample was stratified rather than random, so the 50/33/17 estimate should not be extrapolated precisely. Definite stale residue exists: the locked v3 reasoning document says it supersedes v2 while v2 remains, and the February handoff says previous handoffs were archived while an older January handoff remains in the active directory. No sampled file resembled a one-off scratch probe.

## Scratch behavior

No project has a root `_scratch/` directory. No project-root `notebooks/` directory was found.

Two kinds of `/tmp/` use were observed:

- `AdaW-eWAT-WL-bulkRNAseq`: three files under `02_analysis/wl_subset/` reference `/tmp/`. Two preserve Claude task-output paths such as `/tmp/claude-788715489/-workspaces-AdaW-eWAT-WL-2025/tasks/<id>.output`; one references a conventional WGCNA log.
- `13403-YD_Christina`: two analysis scripts use `/tmp/starsolo_output` and `/tmp/starsolo_tmp_<pool>` as explicit compute storage.

No reference matched the exact `/tmp/claude-<pid>/<hashed-workspace>/<uuid>/scratchpad/` pattern. The `/tmp/` trees themselves were outside the assigned project paths and were not inspected.

## Observations for the designer

- Absence is the dominant field state: six of seven projects have no `docs/_internal/`.
- The one populated tree already spreads decisions and plans across several semantic containers and preserves explicitly superseded material.
- `13403-YD_Christina` demonstrates a parallel tracked convention—52 files under a dated `docs/` bundle—while its 100 internal files remain completely ignored.
- Harness-neutrality has a real migration surface: four projects contain local Claude-specific agent, command, or skill documents beyond settings and toolkit mounts.
- Precise authorship, currentness of all 100 internal files, and a corpus-wide signal ratio could not be determined within the required six-file sampling bound.