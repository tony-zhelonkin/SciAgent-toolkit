# Meta-apply

Propagate every accepted MADR in `docs/_meta/design.md` to every in-flight per-feature design in a single parallel batch. Dispatches N `feature-reviser` subagents in parallel, surfaces one compound diff gate, then dispatches N `architect` subagents in parallel for the architecture-gate step.

Use after `/meta-design` Gate 1 approval, in place of N serial `/design` re-runs, whenever the propagation is mechanical (apply MADR consequences to existing per-feature design/*.md). Features flagged `MANUAL_REDESIGN_NEEDED` by a sub-agent are skipped cleanly — the user runs `/design <slug>` for those.

## Phase 0: Parse arguments

- `$ARGUMENTS` — optional comma-separated list of feature slugs.
  - Empty → auto-select every feature in `docs/_meta/design.md`'s frontmatter `features:` list that ALSO has a `docs/{slug}/design/` directory.
  - Explicit list → those slugs. Every slug MUST appear in `_meta/design.md`'s `features:` frontmatter (else reject: "slug {x} not in `_meta/design.md` scope; re-run /meta-design with updated scope or fix the slug").

Announce resolution:

```
Resolved meta-apply scope:
  features with design/       = [...]       ← feature-revisers will be dispatched for these
  features without design/    = [...]       ← skipped; run /design <slug> for first-draft
  output                      = edits to docs/{slug}/design/*.md (each in-scope feature)
                              + append to docs/_meta/deferred.md
```

If zero features have `design/`, reject: "No features in scope have a `design/` directory. Use `/design <slug>` for each — first-drafts are outside `/meta-apply`'s scope."

## Phase 1: Verify prerequisites

1. **`docs/_meta/design.md` MUST exist with `status: APPROVED`** in frontmatter. If missing or `status: DRAFT`:

   ```
   /meta-apply requires docs/_meta/design.md with status: APPROVED.

   Found: {absent | status: DRAFT | status: OTHER}

   Run /meta-design and approve at Gate 1 first.
   ```

   Stop.

2. **At least one feature in scope must have a `design/` directory** (enforced in Phase 0 already; restate if list is empty after filtering).

3. **No concurrent run in progress.** If any `docs/{slug}/design/*.md` in scope has been modified in the last 60 seconds (mtime-based sanity check), warn:

   ```
   Warning: docs/{slug}/design/*.md was modified recently. Another /meta-apply or /design may be running. Continue anyway? (y/n)
   ```

   WAIT. Otherwise continue.

## Phase 2: Allocate deferred-ID blocks (serialised, pre-dispatch)

**This phase runs in the main agent, BEFORE parallel dispatch.** Parallel appends to `_meta/deferred.md` would race; pre-allocation prevents that.

1. Read `docs/_meta/deferred.md` and find the current maximum ID `D-NNN`. If the file does not exist yet, the watermark is `D-000` (next available is `D-001`).
2. Assign each feature-reviser a contiguous 20-ID block starting above the watermark:
   - feature A (first in scope) → `D-{watermark+1}..D-{watermark+20}`
   - feature B → `D-{watermark+21}..D-{watermark+40}`
   - feature C → `D-{watermark+41}..D-{watermark+60}`
   - ...
3. Announce allocations before dispatch:

   ```
   Deferred-ID block allocations (exclusive per sub-agent):
     normalisation        → D-013..D-032
     umap                 → D-033..D-052
     color-encoding       → D-053..D-072
     theme-bundles        → D-073..D-092
     biological-workflow  → D-093..D-112

   Watermark after this run (if all blocks fill): D-112. Actual new rows will be ≤100.
   ```

The `20` block size is calibrated so a well-scoped sub-agent never exhausts it. A sub-agent reporting its block exhausted is a signal of overreach — the sub-agent will emit `MANUAL_REDESIGN_NEEDED` rather than steal from a sibling's block.

## Phase 3: Parallel dispatch of `feature-reviser` sub-agents

**CRITICAL:** Use a SINGLE assistant message with MULTIPLE `Agent` tool calls — this is what makes them run in parallel. Same pattern as `/review --as all` (see `commands/review.md § Phase 3`). Sequential dispatch violates the design and defeats the cost savings.

For each feature in scope, derive these fields from `_meta/design.md`:

- **`MADRs to apply`**: every MADR section whose "Consequences per feature" block contains this slug on any bullet (grep-style scan; cite the MADR number + one-line consequence quote).
- **`Reader-map sections`**: entries under `_meta/design.md § Reader map` keyed by this slug.

Then, for each feature, issue one Agent call with `subagent_type: "feature-reviser"` and this exact prompt body (parameterised per feature):

```
Apply docs/_meta/design.md MADRs to feature `{slug}`.

Your assigned deferred-ID block: D-{first}..D-{last} (use ≤20 IDs; do not exceed).

MADRs to apply (and the `Consequences per feature: {slug}` lines that name this slug):
  - MADR-{N1} (§{section-header-1}) — consequence for {slug}: "{one-line quote}"
  - MADR-{N2} (§{section-header-2}) — consequence for {slug}: "{one-line quote}"
  - ...

Reader-map sections to re-read (from _meta/design.md § Reader map):
  - {path § ADR-N}
  - ...

Read in this order:
  1. docs/_meta/design.md (binding; all listed MADRs)
  2. docs/_meta/map.md (inter-feature dependencies — cite, never modify)
  3. docs/_meta/deferred.md (current state; your assigned ID block is exclusive to you)
  4. docs/{slug}/map.md, synthesis.md (if present), review/*.md, design/*.md

Produce:
  - Edits under docs/{slug}/design/*.md for every MADR consequence you can apply mechanically
  - Row appends to docs/_meta/deferred.md using your ID block
  - Diff summary as your final message (exact format in your agent definition)

If ANY listed MADR cannot be applied mechanically, emit `status: MANUAL_REDESIGN_NEEDED` in your
diff summary with a one-paragraph explanation per conflict. Do NOT write speculative revisions.
```

Dispatch all N calls in one message. Do not wait for one before sending the next.

## Phase 4: Collect diffs + surface compound report

After all parallel feature-revisers return, read each agent's final diff summary. Then, for each feature, read the post-dispatch state of `docs/{slug}/design/*.md` to confirm the reported edits actually landed (trust but verify — an agent's self-report describes intent, not necessarily effect).

Print the compound report. Table first, per-feature highlights second, MANUAL flags quoted prominently at the end:

```
Meta-apply dispatched {N} feature-revisers. Results:

┌──────────────────────┬─────────┬────────────────┬──────────────┬──────────────┬──────────────┬──────────┐
│ Feature              │ Status  │ Files touched  │ ±Lines       │ Revisions    │ Deferrals    │ Conflict │
├──────────────────────┼─────────┼────────────────┼──────────────┼──────────────┼──────────────┼──────────┤
│ normalisation        │ OK      │ 03, 05         │ +89 / -12    │ 2 ADRs       │ D-013..014   │ —        │
│ umap                 │ OK      │ README, 01, 03 │ +45 / -31    │ 3 ADRs       │ —            │ —        │
│ color-encoding       │ MANUAL  │ (none written) │ —            │ —            │ —            │ ADR-005  │
│ theme-bundles        │ OK      │ README, 03     │ +110 / -46   │ first-draft  │ D-073..075   │ —        │
│ biological-workflow  │ OK      │ 03, 04         │ +120 / -18   │ 3 ADRs       │ D-093..099   │ —        │
└──────────────────────┴─────────┴────────────────┴──────────────┴──────────────┴──────────────┴──────────┘

Per-feature highlights:
  - normalisation — ADR-006 Status OPEN → ACCEPTED (Option B), authorised by MADR-006;
                    ADR-013 opacity-for-degenerate re-homed to marker overlay per MADR-002.
  - umap — ADR-001/ADR-008 gene-set fields renamed per MADR-004; MADR-003 cited in ADR-003.
  - color-encoding — MANUAL_REDESIGN_NEEDED: feature has only design/README.md (no 03-decisions.md
                     yet); planned ADR-005 (opacity = padj-percentile) conflicts with MADR-002 at
                     portfolio level. First-draft + structural rework needed — not mechanical
                     edit. Run /design color-encoding.
  - theme-bundles — runtime column resolver dropped per MADR-001; AC revised for MADR-002 priority
                    chain; bundle-dim composition rule added.
  - biological-workflow — opacity priority chain added; ADR-010 hue revised for MADR-006;
                          ENTITY_PROFILES consumption captured in ADR-011.

MANUAL_REDESIGN_NEEDED flags (run /design <slug> for each — /meta-apply will not touch these):
  - color-encoding: {one paragraph explanation from sub-agent, verbatim}
```

## Phase 5: Compound approval gate (Gate 1 of /meta-apply)

Single gate. Five options:

```
Please review the compound diff above and pick:

  1. Approve all OK-status features
     → I'll flip DRAFT → APPROVED on each approved feature's design/*.md frontmatter
     → Dispatch architect gate in parallel on every approved feature's design/*.md
     → MANUAL-flagged features are left untouched; run /design <slug> manually for each

  2. Approve some, reject others
     → Name the features to reject; I'll revert those edits (git restore on
        docs/{slug}/design/*.md and strip the feature's deferred rows from _meta/deferred.md)
     → Remaining OK features proceed to architect gate

  3. Redo a specific feature
     → Name the feature + a redirect; I'll re-dispatch its feature-reviser
       (re-uses the same deferred-ID block — no re-allocation)

  4. Review a specific feature's diff
     → Name the feature; I'll show file-by-file diff

  5. Questions
```

WAIT for the user.

## Phase 6: On approval — parallel architect dispatch

For every feature that survived the compound gate with `Status: OK`:

1. Bump `status: DRAFT → APPROVED` in the frontmatter of every `docs/{slug}/design/*.md` file (main agent does this directly — trivial, not worth re-dispatching).
2. Dispatch `architect` subagent via Agent tool. **Use a SINGLE assistant message with MULTIPLE Agent tool calls** (same parallel-dispatch pattern as Phase 3 and `/review`). One architect per feature.

Architect prompt template (identical to `commands/design.md § Phase 5`):

```
Review docs/{slug}/design/*.md against docs/{slug}/map.md.

Check internal consistency, completeness, pattern alignment with map.md, and design principles.
Pay special attention to the "Inherited meta-design constraints" block — the feature-reviser
just added or extended it; verify every MADR citation maps to a coherent ADR revision.

Write docs/{slug}/design/review.md with Verdict: READY FOR IMPLEMENTATION | NEEDS ITERATION | NEEDS DISCUSSION.
```

Dispatch all N calls in one message.

## Phase 7: Collect verdicts + emit next-step recommendation

After all parallel architects return, read each `docs/{slug}/design/review.md` and parse the Verdict line. Surface per-feature outcomes:

```
Meta-apply complete — {N} features revised, architect gate run in parallel:

  ✅ normalisation            READY FOR IMPLEMENTATION
  ✅ umap                     READY FOR IMPLEMENTATION
  ⚠️ biological-workflow      NEEDS ITERATION (see docs/biological-workflow/design/review.md)
  ✅ theme-bundles            READY FOR IMPLEMENTATION
  ⏭️ color-encoding           (skipped — MANUAL_REDESIGN_NEEDED at Phase 4)

Deferred rows appended: D-013..D-099 (N new rows total; _meta/deferred.md watermark now D-099).

Recommended next steps:
  /design color-encoding              — for the MANUAL-flagged feature (first-draft + re-architecting)
  /design biological-workflow --iterate — architect flagged NEEDS_ITERATION; loop to READY
  /meta-plan                          — once every feature has Verdict READY, sequence the portfolio
```

For each feature with a NEEDS ITERATION or NEEDS DISCUSSION verdict, **roll back the status bump**: flip `status: APPROVED → DRAFT` on that feature's `design/*.md` frontmatter. The verdict means the architect isn't satisfied; `APPROVED` would be a lie. A subsequent `/design <slug> --iterate` run will re-advance it when Verdict = READY.

For READY features, leave `status: APPROVED` in place — `/plan` requires it.

## Phase 8: Status-field housekeeping (trivial, main-agent)

After Phase 7, the state of every feature's frontmatter is:

| Verdict | design/*.md `status:` |
|---|---|
| READY FOR IMPLEMENTATION | `APPROVED` |
| NEEDS ITERATION | `DRAFT` |
| NEEDS DISCUSSION | `DRAFT` |
| (skipped, MANUAL) | unchanged from pre-run (whatever it was) |

Confirm in one line:

```
Frontmatter status set: {READY count} APPROVED, {NEEDS count} DRAFT, {MANUAL count} unchanged.
```

## Rules

1. **Parallel dispatch is mandatory.** Phase 3 (feature-revisers) and Phase 6 (architects) each use a single assistant message with multiple Agent tool calls. Sequential dispatch defeats the entire point of this command.
2. **Compound gate ONCE, architect gate PER FEATURE.** Gate 1 approves the batch of edits. Then architect runs in parallel on each approved feature and returns per-feature verdicts — no second compound gate on architect output.
3. **MANUAL features are skipped, not failed.** A sub-agent emitting `MANUAL_REDESIGN_NEEDED` is behaving correctly — it refused to speculate. The user's cue is to run `/design <slug>` for that specific feature.
4. **Deferred-ID blocks are pre-allocated and exclusive.** Sub-agents never allocate their own IDs. A sub-agent reporting its block exhausted is overreach → that feature is MANUAL.
5. **Sub-agents never modify `_meta/design.md` or `_meta/map.md`.** Those belong to `meta-architect`. `_meta/deferred.md` is the only shared meta-layer file sub-agents may append to.
6. **Do not auto-progress to `/meta-plan`.** Surface the recommendation; the user starts when ready.
7. **Do not self-loop on NEEDS_ITERATION.** Surface the verdict; the user runs `/design <slug> --iterate` (or re-runs `/meta-apply <slug>` with a redirect) manually. Self-looping would hide architect signal behind another compound turn.
8. **Large diffs (>200 lines touched in a single feature) are a smell, not a defect.** Surface them; suggest the user drill in via option 4 or, if the edits feel aggressive, re-run `/design <slug>` manually for a careful hand pass.
