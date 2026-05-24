---
status: probationary
prune_review_at: 2026-11-20
---

# Audit-slice

Dispatch N concern-driven slicers in parallel over EXISTING code. Each slicer takes one cross-cutting concern and traces it through code + tests + docs + ADRs in execution-time order, tagging coupling by Page-Jones connascence type. Output is one slice file per concern under `docs/_meta/architecture-audit/{date}/slices/`.

This is **Shape γ — a sibling to `/map`, not a replacement.** `/map` is feature-driven, one agent, all surfaces, cheap. `/audit-slice` is concern-driven, N agents, depth along each concern, expensive. Pick the one that matches the work shape.

## Off-the-autopilot gate (read before Phase 0)

`/audit-slice` is deliberately OFF the navigation autopilot. There is no nav-contract auto-fire for it; no other command recommends it as an automatic next step; it **REFUSES auto-piped input** from another agent. The trigger is USER-AUTHORED: a hunch, an endpoint, a feature name, a smell description — typed by a human. The act of a human naming the concerns IS the judgment gate. If the concern list arrives synthesized from another command's output rather than typed by the user, refuse:

```
/audit-slice takes user-authored concerns only. The act of naming the concerns is
the judgment gate — it must not be auto-piped from another agent's output. Please
type the concern(s) you want traced.
```

## Invocation banner (always print first)

Before doing anything, print:

```
/audit-slice is a retrospective architecture-audit tool. A full audit (N slices +
/synthesize-audit) typically runs ~0.8-1M tokens — comparable to a 6-reviewer panel,
but it processes N CONCERNS, not N perspectives. It earns its keep only with an
empirical anchor (a smoke/demo pass, an in-vivo regression, a pre-pivot OSS-core
decision). If you are scoping a NEW feature, run /map instead — it is ~20× cheaper.

Proceed with the audit?
```

WAIT for confirmation before Phase 1. If the user is clearly mid-audit (an audit dir for today already exists), state that and proceed without re-prompting.

## Phase 0: Parse arguments

- `$ARGUMENTS` — one or more concern names / hunches. Multiple concerns may be space-separated bare slugs (`lasso-state theme-bundle-surface discoverability`) or quoted phrases (`'the /export route still calls the legacy writer'`).

If no concern is given, ask:
```
Name the concern(s) to slice. A concern is a cross-cutting thread — a hunch, an
endpoint, a piece of shared state, a smell, or an emergent concept — NOT a feature
and NOT a file. Examples:
  /audit-slice lasso-state theme-bundle-surface
  /audit-slice 'the export path is dead weight'
```

Normalise each concern to a kebab-case slug for filenames; keep the original phrasing to pass to the slicer.

## Phase 1: Verify / create the audit directory

The audit parent is `docs/_meta/architecture-audit/{date}/` where `{date}` is today (`YYYY-MM-DD`).

- If it does not exist, create it plus a `slices/` subdir.
- If it exists (a re-run today), state that and continue — slices are additive; do not clobber existing slice files. If a slice file for a given concern slug already exists, ask the user whether to overwrite or pick a new slug.

Announce:
```
Audit dir: docs/_meta/architecture-audit/{date}/
Slices dir: docs/_meta/architecture-audit/{date}/slices/
Concerns ({N}): {list}
```

## Phase 2: Dispatch slicers in parallel

**CRITICAL:** Use a SINGLE assistant message with MULTIPLE Agent tool calls (`subagent_type: "slicer"`), one per concern — this is what makes them run in parallel. Do not dispatch sequentially.

Default model per slicer is the agent's default (`opus`). For a concern the user flags as shallow (a single dead route, a narrow rename), note in that slicer's prompt that sonnet-class depth is acceptable — keep the traversal to the spine.

Number slices `01..N` in the order concerns were given.

Prompt template (one per concern):
```
Trace the concern `{original concern phrasing}` through the codebase.

Output path: docs/_meta/architecture-audit/{date}/slices/{NN}_{concern-slug}.md

Trace in EXECUTION-TIME order: input → transform → second-writer/reconciliation →
consumer. Re-explore the code; verify line numbers at HEAD. Cover code + tests +
docs + ADRs that lie ON this concern's spine; exclude everything off-spine.

Tag every smell with its Page-Jones connascence type (Name/Type/Convention/
Execution-Order/Location/Identity/Value/Algorithm) and state which edge type
(direct-call / shared-state / background-knowledge) it maps to for synthesis.

Produce: traversal, branching points, findings (E-NN with file:line + connascence
tag, including negative-space findings), cross-slice touchpoints (name the other
concern slugs), decoupling proposals (R-NN), open questions for synthesis.

Diagnose, do not redesign. Do not assign core/seam/removable — that is the
synthesis's judgment. One file only.
```

## Phase 3: Write the slice-index README

After all slicers complete, write `docs/_meta/architecture-audit/{date}/README.md` (overwrite if a re-run today added slices):

```markdown
---
date: YYYY-MM-DD
git_sha: {HEAD short SHA}
status: probationary
---

# Architecture audit — {date}

Concern-driven slices. Each traces one cross-cutting concern in execution-time order.
This audit is the JUDGMENT layer; pair it with `/components-extract` (the deterministic
substrate) and feed both to `/synthesize-audit`.

## Slices
| # | Concern | Findings | Decoupling proposals | Connascence flagged |
|---|---------|----------|----------------------|---------------------|
| 01 | [{concern}](slices/01_{slug}.md) | {E-count} | {R-count} | {tags} |
| ...|

## Cross-slice touchpoints
{Aggregate the touchpoints each slicer reported — where two concerns collide.}
```

## Phase 4: Route to synthesis

```
Audit slices complete ({N} files) under docs/_meta/architecture-audit/{date}/slices/.
Slice index: docs/_meta/architecture-audit/{date}/README.md

Per-slice headlines:
  01 {concern}: {first finding's headline}
  ...

Recommended next step:
  /synthesize-audit               — integrate these N slices + the /components-extract
                                    substrate into the full renderable components.json

If you have NOT yet run /components-extract on this repo, run it first — the synthesis
merges the deterministic substrate (physical files, static edges, metrics) with the
slices' judgment. Without the substrate the synthesis has no trust scaffold.
```

## Rules

1. **User-authored trigger only.** Refuse auto-piped concern lists. The human typing the prompt is the gate.
2. **Parallel dispatch is mandatory.** One message, multiple Agent calls. Sequential dispatch violates the design.
3. **Slicers re-explore; reviewers do not.** This is the distinguishing capability — do not constrain a slicer to read-only artifacts the way `/review` constrains reviewers.
4. **One concern per slicer, one file per slice.** If a slicer writes elsewhere or merges concerns, flag it.
5. **Additive, never clobbering.** A re-run today appends slices; never overwrite an existing slice without asking.
6. **This command is probationary** (`prune_review_at: 2026-11-20`). It is the right tool ~once per quarter per major project, not per-feature. Police the trigger; if you find yourself reaching for it to scope a new feature, you wanted `/map`.
