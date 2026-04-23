# Architect (standalone gate)

Re-run the architecture consistency/completeness gate on a feature's existing `design/*.md` without re-entering `/design`. Dispatches the `architect` agent; produces `docs/{feature}/design/review.md` with a Verdict (`READY FOR IMPLEMENTATION | NEEDS ITERATION | NEEDS DISCUSSION`).

Use when you:
- Hand-edited one or more `design/*.md` files and want a fresh Verdict without re-running all of `/design`.
- `/meta-design` was re-run (MADR revised) and want to check whether this feature's design still holds.
- `/meta-apply` was partial (MANUAL flag) and you want to re-gate after a subsequent `/design` revision.
- A peer feature's APPROVED design changed and you want a consistency re-check.

`/design` continues to dispatch architect internally as its Phase 5 — this command is the *standalone* entry point for the same gate. The `architect` agent's behaviour is identical in both cases; only the invocation path differs.

## Phase 0: Parse arguments

- `$ARGUMENTS[0]` — feature slug (required). Reject if missing: "Usage: /architect <feature>".
- `$ARGUMENTS[1+]` — optional redirect notes the user wants surfaced to architect (e.g., "focus on ADR-014's ship-gate language").

## Phase 1: Verify prerequisites

1. **`docs/{feature}/map.md` MUST exist.** If absent:

   ```
   docs/{feature}/map.md not found. Architect reads the map to check design against it.
   Run /map {feature} first.
   ```

   Stop.

2. **`docs/{feature}/design/` MUST exist and contain at least `README.md` + `01-architecture.md` + `02-behavior.md` + `03-decisions.md`.** If the directory is empty or only has `README.md`:

   ```
   docs/{feature}/design/ is incomplete ({found files}). Architect needs the full bundle.
   Run /design {feature} to draft the remaining docs first.
   ```

   Stop.

3. **Frontmatter sanity.** Every `docs/{feature}/design/*.md` should have `status: APPROVED` (if you're gating after a full `/design` approval) OR `status: DRAFT` (if you're gating a hand-edit in progress). Both are fine — architect can verdict on either. If `status:` is missing or something else, warn once and proceed.

## Phase 2: Announce the re-gate

```
Re-dispatching architect gate on docs/{feature}/design/*.md.

Bundle under review:
  - README.md         status: {status}   date: {date}
  - 01-architecture.md  status: {status}   date: {date}
  - 02-behavior.md    status: {status}   date: {date}
  - 03-decisions.md   status: {status}   date: {date}
  [+ conditional files if present]

Prior review ({if docs/{feature}/design/review.md exists}):
  Verdict: {prior Verdict line}
  Date: {prior date}

{If user supplied redirect notes in $ARGUMENTS[1+], quote them here.}

Dispatching. Wait for Verdict.
```

## Phase 3: Dispatch architect

Use the Agent tool with `subagent_type: "architect"`.

Prompt template (identical to `commands/design.md` Phase 5):

```
Review docs/{feature}/design/*.md against docs/{feature}/map.md.

Check internal consistency, completeness, pattern alignment with map.md, and design principles.

Also read `## User resolutions`, `## Binding chat decisions`, and `## Post-READY amendments` sections in every upstream artifact (map.md, synthesis.md, review/*.md, design/03-decisions.md). Treat them as equal-priority to ADRs; flag any design doc that silently contradicts a scribed decision.

{If _meta/design.md exists:}
Also verify every MADR in docs/_meta/design.md whose "Features affected" list contains `{feature}` is either honored locally or explicitly flagged in a "Meta overrides" block.

{If user supplied redirects in $ARGUMENTS[1+]:}
User redirects for this review:
{quote the redirects verbatim}

Write docs/{feature}/design/review.md with Verdict: READY FOR IMPLEMENTATION | NEEDS ITERATION | NEEDS DISCUSSION.
```

## Phase 4: Handle the verdict

Read `docs/{feature}/design/review.md` and parse the Verdict line.

**READY FOR IMPLEMENTATION:**

```
✅ Design approved and architect gate passed.

Design: docs/{feature}/design/
Review: docs/{feature}/design/review.md (Verdict: READY FOR IMPLEMENTATION)

If you hand-edited design docs after a prior APPROVED state:
  Frontmatter status may still read DRAFT — bump to APPROVED if the hand-edit was a deliberate
  revision rather than a typo fix. Main agent can do this directly; no /design re-run needed.

Next: /plan {feature}
```

If any `design/*.md` frontmatter has `status: DRAFT`, offer to flip it to APPROVED in one edit (main agent). Ask the user first:

```
design/*.md still have status: DRAFT. Flip to APPROVED now? (y/n)
```

Default behaviour: wait for the user's y/n before flipping. Do not auto-flip.

**NEEDS ITERATION:**

1. Surface architect's issues to the user (paste the Verdict's "Issues" section).
2. Ask: "Do you want me to apply the listed fixes, or will you hand-edit? After fixes land, re-run `/architect {feature}` to close the loop."
3. Wait. This command does NOT auto-loop — that's `/design`'s job. `/architect` is the standalone gate; looping belongs to the orchestration command.

**NEEDS DISCUSSION:**

1. Surface architect's open questions.
2. Wait for user input.
3. Once the user has answered, note which doc needs the update. Offer to edit or hand back for user edit. Re-run `/architect {feature}` after.

## Rules

1. **Standalone entry point, not a loop orchestrator.** `/architect` dispatches one architect pass and surfaces the Verdict. To loop iteratively on NEEDS_ITERATION, re-invoke `/architect {feature}` after each fix. The `/design` command does auto-loop (its Phase 6) because the drafting and gating are coupled there; standalone gating is deliberately one-shot so the user stays in control between iterations.
2. **Prerequisites are lighter than `/design`.** No upstream-artifact re-read (no map, synthesis, reviews). Architect already reads those from inside its own agent prompt.
3. **No frontmatter auto-flips.** `DRAFT → APPROVED` happens on explicit user consent (Phase 4), never implicitly. Users hand-edit design docs often; their intent between "tweak" and "ship" must remain explicit.
4. **Never dispatch any other agent.** Only `architect`. If the verdict says "re-run synthesis" or "re-map", surface the recommendation; do not act on it from this command.
5. **Respect the meta-layer when present.** If `docs/_meta/design.md` exists, architect must treat its MADRs as binding on this feature per the MADRs' "Features affected" lists. The prompt template in Phase 3 already wires this in.
