# Plan

Decompose an approved design into ordered implementation phases.

## Phase 0: Parse arguments

- `$ARGUMENTS[0]` — feature slug (required)

## Phase 1: Verify prerequisites

Read:
1. `docs/{feature}/design/README.md` — must have `**Status**: APPROVED`
2. `docs/{feature}/design/01-architecture.md`
3. `docs/{feature}/design/02-behavior.md`
4. `docs/{feature}/design/03-decisions.md`
5. `docs/{feature}/design/review.md` — must have Verdict `READY FOR IMPLEMENTATION`
6. `docs/{feature}/map.md` — for existing patterns to follow
7. `docs/_meta/plan.md` — if present, READ. Identify the portfolio-phase slot assigned to this feature; decompose within that slot. Do NOT produce phases that violate the portfolio sequencing (e.g., do not assume a primitive that the meta-plan says another feature ships first).

**Scribe-on-latest read rule.** For every upstream artifact above, also read any `## User resolutions`, `## Binding chat decisions`, or `## Post-READY amendments` section present. Treat those as equal-priority to design ADRs and cite them explicitly when they shape a phase boundary.

If design status or architect verdict is not ready:
```
Design is not ready for planning. Current state:
  - README status: {status}
  - Architect verdict: {verdict}

Run /design {feature} to iterate until both are ready.
```

Stop.

## Phase 2: Decompose into phases

Break the implementation into ordered phases where:
- Each phase produces something testable
- Each phase builds on the previous one
- No phase touches more than ~3-5 files
- Data models / types come first; adapters / UI come last

Typical ordering:
1. Data models / types / schemas (the nouns)
2. Core business logic (the verbs)
3. Storage / persistence
4. API / interface layer
5. Integration and wiring
6. Error handling / edge cases

## Phase 3: Write plan files

Output directory: `docs/{feature}/plan/`

### README.md
```markdown
---
feature: {slug}
status: APPROVED
date: YYYY-MM-DD
---

# Implementation Plan: {feature}

## Overview
1-2 sentences describing what's being built.

## Phase Summary
| Phase | Focus | Files | Depends on |
|-------|-------|-------|------------|
| 1 | Data models | ... | — |
| 2 | Core logic | ... | Phase 1 |
| ... |

## Success Criteria (from design README)
- [ ] ...
```

### phase-NN.md (one per phase)
```markdown
# Phase {N}: {Focus}

## Context
What previous phases produced — types, interfaces available now.

## Files to Create
### `path/to/new_file.ext`
- Purpose
- Key implementation details
- Patterns to follow (reference map.md or existing design docs)
- Interfaces to implement

## Files to Modify
### `path/to/existing.ext`
- What changes
- Lines affected (from map.md references)

## Verification
- [ ] Phase-specific check 1
- [ ] Phase-specific check 2

## References
- design/01-architecture.md §...
- design/02-behavior.md §...
- map.md for patterns to match
```

## Phase 4: Surface plan with a glance-check

The plan is **already stamped `status: APPROVED`** in Phase 3 — this is not a gate on the plan's status. It's a one-shot ack so the meta-architect can catch a wrong phase boundary before implementation momentum starts. Silence or "proceed" moves on; "changes" loops back to revise.

```
Implementation plan drafted at docs/{feature}/plan/ and stamped APPROVED (architect verdict READY on design implies plan is derivative):
  - README.md
  - phase-01.md  ({focus})
  - phase-02.md  ({focus})
  ...

Phase summary:
  1. {focus}
  2. {focus}
  ...

Glance-check — is a phase boundary obviously wrong?
  - "proceed" / silence  — I'll recommend /implement
  - "change phase N: {what}"  — I'll revise that phase and re-stamp
  - "questions"           — ask about specific phases

(Plan is already APPROVED; this ack does not re-gate.)

Next: /implement {feature} 1   (single phase, hand-paced)
      /implement {feature} --auto  (end-to-end, stops only at failure/uncertainty/completion)
```

Do NOT wait indefinitely for approval — surface the plan, recommend the two `/implement` invocations, and stop. If the user says "change phase N" within this turn, revise in place; otherwise the command is complete.

## Rules

1. **Plan is NOT design.** Do not re-decide architecture. If you find yourself redesigning, stop and run `/design` again.
2. **Each phase is testable.** If a phase produces untestable output, split it or fold it into a neighbour.
3. **Reference, don't duplicate.** Phase files reference design docs and map.md — they don't restate them.
4. **Plan is stamped APPROVED at draft time.** Phase 4's ack is a glance-check, not a status gate. The real gate was `/design`'s Phase 4 + architect `READY` verdict; re-gating here would re-litigate an already-made decision. If the user spots a wrong phase boundary, revise in place — do not flip status back to DRAFT.
