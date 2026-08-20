# Field findings — 24 provc-managed projects, 2026-08-20

Produced by four parallel read-only inspections (`inspect-*.report.md`). Read
those for per-project detail. This digest is what the design must answer to.

## 1. Absence is the dominant state

**10 of 24 projects have no `docs/_internal/` at all.** Six of seven in the
long-tail slice; two of six in dc-nexus; two more elsewhere. Any design that
assumes the directory exists is designing for a minority of the fleet.

## 2. The field already invented the nested repo — twice, independently

14761-DM and 14782-DM keep `docs/_internal/` ignored by the parent **and
tracked in a nested git repo inside it**:

| Project | nested commits | tracked files |
|---|---|---|
| 14761-DM | 126 | 195 |
| 14782-DM | 210 | 708 |

**Neither has a remote.** Both sit on branch `master` with no upstream. So the
memory is versioned and still cannot leave the machine — it dies with the disk,
not with the session. This corrects the framing in `BRIEF_consult.md` §1: "zero
tracked" is true only from the parent repo's perspective, and the gap is not
versioning but *reachability*.

That two projects converged on this without instruction is the strongest
available evidence for topology (c)/the nested-repo hypothesis. Weigh it.

## 3. The scaffold is mostly empty, and the empty parts mislead

The five STING child projects each carry ~10 scaffold files (`README.md`,
`scientific-context.md`, five `.gitkeep`, three per-directory `README.md`)
before any content. And in **every** child, `handoffs/` contains only
`.gitkeep` while the real handoffs were written to `sessions/`. Same in
14761-DM and 14782-DM: `handoffs/` holds a placeholder, `sessions/` holds the
content.

A directory agents never use is worse than absent: it is a signpost pointing
away from where the content is.

## 4. Naming drift, measured

| Concept | Spellings observed in the field |
|---|---|
| Handoff | `handoff/`, `handoffs/`, `sessions/`, `states/`, plan-local `_handoff/`, `HANDOFFS.md`, root `handoff.md`, `handoff_YYYYMMDD_HHMMSS.md`, `HANDOFF_YYYY-MM-DD.md`, `.handoff_archive/` |
| Reasoning | `reasoning/`, `reports/`, `*-decisions.md`, `*-rulings.md`, `workflow-architecture/`, `app/<topic>/00_SYNTHESIS.md` |
| Plans | `plans/`, `ai-generated/plans/`, `loop/cards/`, `app/<topic>/PLAN.md` |
| Dated dir | `YYYY-MM-DD-slug`, `YYYY-MM-DD_slug`, `YYYY-MM-DD__slug`, `YYYY-MM-DD subject research` (with spaces) |
| Research | `research/`, `web/`, `documentation/`, `notes/`, `explainers/`, `tools/<tool>/` |
| Scratch | `docs/_internal/scratch/`, `docs/_internal/_scratch/`, `03_results/_scratch/`, `02_analysis/scripts/_scratch/`, `/tmp/...` |
| Other one-offs | `contracts/`, `issues/`, `states/`, `inbox/`, `reply-package/`, `loop/`, `consensus-validation/`, `personal-notes/`, `_superceded/` (sic), `grant/` |

Also inconsistent: the ignore rule itself. Some projects carry
`docs/_internal/` inside a managed `SCIAGENT:GITIGNORE` block; others have a
hand-written generic `_internal/` (Meta-Aging:65, DC-nexus:48, STING-JR:69) that
ignores *any* directory of that name anywhere; 13403-YD_Christina has both
`docs/_internal/*` and `docs/_internal/` on adjacent lines.

## 5. The landfill is real and it is not documents

| Project | `docs/_internal/` size | What dominates |
|---|---|---|
| DC-nexus | **568 MB** | `grant/local/cri/raw/.venv/` — 3,813 files, **91.5% of the tree** |
| 14782-DM | **231 MB** | a 196 MB model checkpoint; 210 hashed cache JSONs; parquet, bytecode, logs |
| mouse_anchor | 15 MB | 49 PNGs, 2 PPTX, 1 PDF |

`docs/_internal/` has become the place things go when they have nowhere else —
including a Python virtualenv. Retention is not a refinement to add later; it is
the difference between a memory system and a dumping ground.

## 6. `_scratch/` at the project root does not exist

`craft.yaml` states `_scratch/` is "the only sanctioned throwaway zone."
**No project has one at the root.** One (mouse_anchor) has an empty one. What
exists instead: `03_results/_scratch/`, `02_analysis/scripts/_scratch/`,
`docs/_internal/scratch/`, `docs/_internal/_scratch/` — and `/tmp`.

The always-on instruction names a location the fleet does not have.

## 7. `/tmp` escape is confirmed, and it is harness-shaped

The exact pattern appears in 14761-DM and 14782-DM review reports:
`/tmp/claude-788715489/-workspaces-Meta-Aging/<uuid>/scratchpad/`. AdaW
preserves `/tmp/claude-788715489/-workspaces-.../tasks/<id>.output`. Several
14782-DM reports state outright that the scratch scripts lived outside the
repository.

Note the path is harness-chosen, keyed by PID and per-session UUID, outside
every project. No project-level convention can reach it after the fact.

## 8. Harness neutrality is not the current state

- `.codex/` exists in **zero** projects — Codex agents leave no project-level
  footprint at all, while Claude leaves `settings.local.json`,
  `scheduled_tasks.lock`, hooks, and status lines.
- `.gemini/settings.json` exists in 8+ projects.
- Locally-authored, non-mounted content has accumulated inside harness dirs:
  a `marimo-pair` skill under `.agents/` in 14782-DM, 8 tracked `.claude`
  agent/skill documents in 12868-EH, 6 deprecated skill docs in JBader_scHFD,
  an untracked `kickoff.md` command in 13403-YD_Christina.

Content authored *for one harness* is content the other harnesses cannot see.

## 9. Signal-to-noise, sampled

Purposive samples, not random: roughly 83% durable in the small trees
(Meta-Aging, 14616-DM), ~67% in 14761-DM, ~50% in 13403-YD_Christina, and
**2/6** in 14782-DM — where the structural evidence (caches, checkpoint,
bytecode) is harsher than the sample suggests. The pattern: the bigger the
tree, the lower the density.

Explicit stale residue is common and unlabelled — a `SYNTHESIS_v3_LOCKED.md`
that supersedes a v2 still sitting beside it; a February handoff saying prior
handoffs were archived while a January one remains active; a within-project
naming conflict where the top-level README mandates dateless ADR names and
`reasoning/README.md` mandates dated ones.

## 10. One project solved continuity a different way

13403-YD_Christina has 100 ignored files under `docs/_internal/` **and** 52
tracked files under a dated `docs/2026-05-25/` bundle. Its
`.claude/commands/kickoff.md` routes agents to the *tracked* bundle for
handoffs, decisions, and phase plans. The reliable record is the one outside
`_internal/`, because that is the one that survives.
