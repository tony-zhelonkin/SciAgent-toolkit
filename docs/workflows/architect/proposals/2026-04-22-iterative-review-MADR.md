---
doc: iterative-review-MADR
date: 2026-04-22
status: ACCEPTED
author: user + assistant (shadartist@gmail.com session)
scope: SciAgent-toolkit/commands/{review,synthesize}.md + docs/workflows/architect/*.md + skills/architecture-first-dev/SKILL.md
motivation: During a `leiden-stability` review the user organically invented a pattern the harness does not name — a second review round where each reviewer saw the other round-1 outputs + the synthesis + the user's locked design direction, to reconcile after round-1 divergence. Today `/review` silently overwrites prior reviews and has no concept of "round" or "cross-informed." Round-1 reviews are deliberately parallel-independent to surface divergence; a second round serves a different purpose (reconcile) and needs a different regime.
---

# MADR-H5 — Iterative review (`--iterate`) with `.history/` archive and cross-informed regime

**Status.** ACCEPTED 2026-04-22.

**Context.** `commands/review.md` today writes `docs/{feature}/review/<reviewer>.md` once per reviewer. Re-running `/review` overwrites those files — the round-1 positions are lost. In practice users iterate: after reading `/synthesize` output, they lock a design direction, re-map to see the blast radius, and want another reviewer pass that reconciles the round-1 divergences in light of the locked direction.

Two distinct primitives emerge from this pattern:

1. **Review rounds** — multiple reviewer passes over a feature, each preserved for audit (why `bioinf` shifted position between round 1 and round 2 should be citeable).
2. **Cross-informed vs. independent regime** — round-1 reviewers are deliberately blind to each other (surfacing divergence is the goal); round-2 reviewers read each other's round-1 outputs + the synthesis + the user's locked direction (reconciliation is the goal). These are qualitatively different regimes that should not silently leak into each other.

The harness has two adjacent mechanisms it composes with:

- **MADR-H3 scribe-on-latest** (just shipped) puts the user's binding chat decisions into `synthesis.md § User resolutions`. Round-2 reviewers reading `synthesis.md` see the locked direction automatically — no new mechanism needed for direction-injection.
- **ADR-018 gate-placement discipline** — judgment gates are preserved; this MADR adds a judgment-gate-respecting iteration primitive, not a bypass.

The existing OQ 2 in `docs/workflows/architect/03-decisions.md` anticipated this exact question; this MADR resolves it.

**Options.**

- (a) Suffix in filename (`bioinf_v1.md`, `bioinf_v2.md`) — rejected. Downstream readers (design, architect, plan) would need a max-version computation, and every citation in a design doc must pick a version. Worst of the three.
- (b) Subdir by round (`review/round-1/bioinf.md`, `review/round-2/bioinf.md`) — transparent, but forces every downstream reader to learn `glob round-{max}/*.md` instead of `review/*.md`.
- (c) **`.history/` archive** — `review/bioinf.md` is always the latest; on iterate, prior round is moved to `review/.history/round-N/bioinf.md`. Downstream readers need zero changes; MADR-H3's "latest upstream artifact" rule composes without edits; git already handles cross-feature audit.
- (d) Always cross-inform — rejected. Round-1's job is divergence-surfacing; cross-informing would reduce it to a premature-consensus machine.

**Decision.** (c) + distinct regimes.

### File layout (`.history/` archive)

```
docs/{feature}/
├── review/
│   ├── bioinf.md                     ← always latest round
│   ├── wetlab.md
│   ├── ...
│   └── .history/
│       ├── round-1/
│       │   ├── bioinf.md             ← archived on iterate
│       │   └── ...
│       └── round-2/
│           └── ...
├── synthesis.md                      ← always latest round
└── synthesis/
    └── .history/
        ├── round-1.md
        └── round-2.md
```

Round-1 is the common case; `/review` without `--iterate` writes `review/<name>.md` directly (no history yet — no prior round to archive). First `--iterate` creates `.history/round-1/` and writes the new round to `review/`. Subsequent iterates increment `round-N`.

### `/review --iterate` — cross-informed regime

New flag on `commands/review.md`:

- `/review <feature> --iterate [--as <spec>] [--but <csv>]` — dispatches a cross-informed reviewer round.
- Before dispatch: move every existing `review/*.md` to `review/.history/round-{max+1-1}/<name>.md` (i.e., archive current round as round-N; new round is N+1).
- Dispatch prompt is *different* from round-1: each reviewer receives a composite context — their own prior-round output, every other reviewer's prior-round output, the current `synthesis.md` (including `§ User resolutions` where the locked direction lives), the current `map.md`, and any `design/` drafts that exist.
- The dispatch prompt explicitly tells each reviewer: *"You are refining, not re-doing. Update your position in light of the sibling reviews and the user's locked direction. For each position you update, note it explicitly as an update; for each position you stand by, say you stood by it — do not silently re-write your round-1 as if you'd always held the round-2 view."*
- Reviewer agents are unchanged; the iterate discipline lives at the call site, not in the agent system prompts. This keeps round-1 behaviour clean and makes the iterate-regime explicit-at-call.

### `/synthesize` mirrors the pattern

- First-run: writes `synthesis.md` (unchanged from today).
- On re-run (when `review/*.md` files are newer than `synthesis.md` — i.e., an iterate round just landed): move existing `synthesis.md` to `synthesis/.history/round-{max}.md`, then write the new `synthesis.md`.
- Downstream readers always read `synthesis.md` = latest. MADR-H3's `§ User resolutions` section is carried into the new synthesis by the `synth` agent (it reads the prior synthesis as an input on iterate runs).

### Downstream readers — no changes required

`commands/{design,architect,plan,implement,verify,diagram}.md` all continue to read `review/*.md` and `synthesis.md`. The "primary file is always latest" convention makes the round abstraction transparent to them. The MADR-H3 scribe-on-latest rule composes without modification.

**Trade-offs.**

- **Gains.** Named, non-destructive iteration primitive. Round-1 divergence and round-2 reconciliation are distinct regimes with distinct purposes — the harness reflects the epistemic structure rather than conflating. Audit trail is preserved cheaply (within-feature history file, not just git). MADR-H3 + ADR-018 compose cleanly; no breaking changes downstream.
- **Costs.** `.history/` adds one directory per iterated feature. Reviewer dispatch prompt for `--iterate` is longer (sibling reviews injected) — token cost per reviewer on iterate rounds is ~2-3× a round-1 reviewer. This is acceptable: iterate rounds are rare (you don't iterate a feature you got right the first time) and each reviewer is deliberately doing harder work (reconciliation, not fresh assessment).
- **Risk.** Round-2 reviewers reading sibling round-1s could drift toward consensus prematurely. Mitigated by the explicit dispatch-prompt discipline ("refine, don't re-do; update *or* stand by, explicitly") and by preserving round-1 in `.history/` so the degree of position-shift is visible.

**Harness edits.**

- `commands/review.md` — add `--iterate` flag, Phase 2.5 (archive prior round on iterate), Phase 3b (cross-informed dispatch prompt for iterate rounds); Rule 6 ("iterate is always cross-informed; `/review` without `--iterate` is always parallel-independent"); Rule 7 ("never overwrite a review file; always archive to `.history/round-N/` before replacement").
- `commands/synthesize.md` — Phase 1.5 (archive prior synthesis on re-run if `review/*.md` mtime > `synthesis.md` mtime); Rule 4 ("synthesize never destroys; prior synthesis archived to `synthesis/.history/round-N.md`"). The `synth` agent receives both the new reviews and the prior synthesis as inputs on iterate runs.
- `docs/workflows/architect/00-quickstart.md` — add `/review --iterate` row in commands-at-a-glance; add per-stage tip explaining the cross-informed regime and when to reach for it.
- `docs/workflows/architect/03-decisions.md` — close OQ 2 (cite this MADR); add ADR-019 documenting the `.history/` + cross-informed decision.
- `skills/architecture-first-dev/SKILL.md` — cheat-sheet row for "synthesis.md exists, user locked a design direction, want to reconcile round-1 divergence" → `/review --iterate`; framing sentence on the round-1 vs round-2 regime distinction.
- No agent YAML changes. No role YAML changes.

**Adoption order.** Single MADR; single landing. The edits are scoped to two commands + three docs + the skill. No migration of existing feature directories is required — `.history/` is only created lazily on first iterate.

## Approval record

- [x] User approval: 2026-04-22
- [ ] Harness edits committed: commit SHA
- [ ] Tested on first iterate-capable feature: feature slug + date (leiden-stability is the natural first target)
