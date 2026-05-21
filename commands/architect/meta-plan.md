# Meta-plan

Sequence per-feature implementation phases across the portfolio so dependencies land before their consumers, and so simultaneously active phases don't collide on the same files.

## Phase 0: Parse arguments

- `$ARGUMENTS` — optional comma-separated list of feature slugs
  - Empty → use the same scope as `_meta/design.md` (read its `features:` frontmatter)
  - Explicit list → must be a subset of `_meta/design.md`'s scope

Announce:

```
Resolved meta-plan scope:
  features = [{slug}, {slug}, ...]
  output    = docs/_meta/plan.md
```

## Phase 1: Verify prerequisites

1. **`docs/_meta/design.md` MUST exist with `status: APPROVED`.** If absent or DRAFT:

   ```
   /meta-plan requires docs/_meta/design.md with status APPROVED.

   Current: {missing | DRAFT}
   Run /meta-design (and approve at Gate 1) first.
   ```

   Stop.

2. **Per-feature plan input.** For each feature in scope, check whether `docs/{slug}/plan/phase-NN.md` files exist. If a feature has no plan files:

   ```
   Warning: feature {slug} has no plan/ directory — its phases cannot be sequenced.

   Choose:
     1. Continue — meta-plan will note "{slug}: no phases yet" and skip sequencing for it
     2. Stop — run /plan {slug} first, then re-invoke /meta-plan
   ```

   WAIT.

3. **At least one feature must have plan files.** If every feature in scope lacks plans, reject — there is nothing to sequence.

## Phase 2: Check existing state

- **`docs/_meta/plan.md` exists** → ask:

  ```
  docs/_meta/plan.md already exists ({date}, covers {N} portfolio phases).

  Choose:
    1. Use existing
    2. Regenerate
    3. Iterate with redirects
  ```

  WAIT.

- **Does not exist** → state `No existing _meta/plan.md — proceeding to meta-architect dispatch.` and continue.

## Phase 3: Dispatch meta-architect

Use the Agent tool with `subagent_type: "meta-architect"`.

Prompt template:

```
Produce docs/_meta/plan.md.

Features in scope: {comma-separated slugs}

Read in this order:
  1. docs/_meta/design.md (REQUIRED — MADR consequences drive ordering)
  2. docs/_meta/map.md (for inter-feature dependencies)
  3. For each feature with a plan/:
     - docs/{slug}/plan/README.md (phase summary table)
     - docs/{slug}/plan/phase-NN.md (per-phase Files-to-Create / Files-to-Modify)
     - docs/{slug}/design/03-decisions.md (for ADR cross-references)
  4. For each feature WITHOUT a plan/:
     - Note "{slug}: no phases yet" in the output; do not invent phases

Produce docs/_meta/plan.md following the _meta/plan.md template in your system prompt. It MUST contain:
  - Dependency graph (Mermaid) — feature-phases as nodes, blocking relationships as edges
  - Portfolio-phase order — table assigning every feature-phase to a portfolio phase P0..PN
  - Shared-file collision check — for each portfolio phase, list files modified by every active phase and flag overlaps
  - Parking items — items from _meta/deferred.md flagged for a future wave
  - Reader map

A phase that depends on another (per MADR consequence or _meta/map.md inter-feature dependency) MUST land in a later portfolio phase than its dependency.

If two simultaneously active feature-phases touch overlapping line ranges in the same file, flag it as "overlap — split or sequence" rather than silently sequencing.

Read but do not modify per-feature plan.md or design.md files.
```

## Phase 4: Gate — present sequencing for approval

Read `docs/_meta/plan.md` and present:

```
Meta-plan drafted at docs/_meta/plan.md.

Portfolio phase order:
  P0 — {features and their phases}    rationale: {one line}
  P1 — ...
  ...

Dependency graph:
  {render the Mermaid graph as ASCII or quote it raw}

Collision check:
  - {N} clean phases
  - {M} overlap warnings (these are blockers):
      P{N}: {file} — touched by {fA phase-X, fB phase-Y}

Parking items: {count from deferred.md activated for future waves}

Please review and:
  1. Approve — sequencing is correct
  2. Request resequencing — name the swap (e.g., "move biological-workflow phase-01 to P1")
  3. Resolve a collision — name the file and which feature should land first
  4. Questions
```

WAIT.

## Phase 5: Iteration

If the user requests changes, capture redirects and re-dispatch (Phase 3) with redirects appended. Loop until approved.

## Phase 6: Surface downstream

```
Meta-plan approved: docs/_meta/plan.md

Portfolio-phase assignments:
  - {feature}: phase-01 → P0, phase-02 → P2
  - {feature}: phase-01 → P2
  - {feature}: no plan/ — run /plan {slug} before participating

Recommended next steps:
  /plan {feature}                     — for any feature whose plan does not yet exist
  /implement {feature} 1              — start P0 work (per portfolio sequence)
```

## Rules

1. **Meta-plan is sequencing, not phase decomposition.** Per-phase decomposition is `/plan`'s job.
2. **No collision is silently sequenced.** A flagged overlap requires user resolution before any P{N} containing it goes live.
3. **Dependency edges come from `_meta/design.md` (MADR consequences) and `_meta/map.md` (inter-feature dependencies).** Do not invent dependencies; cite the source.
4. **A feature without a plan is noted, not faked.** Do not generate placeholder phases.
5. **Do not auto-progress to `/implement`.** Surface the recommendation; let the user start when ready.
