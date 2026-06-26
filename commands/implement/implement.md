# Implement

Execute a **multi-agent, multi-pass implementation campaign** over a solidified
phase-plan produced by `/decompose` — dispatching dedicated Opus implementer +
reviewer pairs per phase, sliding-window architecture reviewers across seams,
and maintaining a durable ledger that makes the campaign resumable after any
interruption.

This command is the **continuation of `/decompose`**: it consumes the
`DECOMPOSE-COMPLETE` handoff block and drives every phase from stub to verified
artifact. Implementation is explicitly its only scope — do **not** use it for
planning (use `/decompose`) or mechanical drift-checks (use `/verify`).

Every subagent dispatched here is an **Opus** subagent (`model: opus`), except
purely mechanical fixes surfaced by reviewers, which may be delegated to a
Sonnet implementer to avoid burning Opus budget on line-edits. The main agent
is the conductor: it dispatches, gates, and records — it does not itself
implement phases.

## Usage

```
/implement <slug> [--phase N] [--auto] [--resume <slug>]
```

- **`<slug>`** (required) — the kebab-case plan identifier matching
  `docs/_internal/plans/<slug>/` produced by `/decompose`. Also accepted as a
  bare path to the plan directory.
- **`--phase N`** — run exactly phase N in isolation (hand-paced, default mode
  when no `--auto`). Stops after that phase and waits.
- **`--auto`** — autonomous end-to-end mode: drive all remaining phases to
  completion, stopping only at a hard FAIL, an escalation the conductor cannot
  resolve, or campaign completion. Implied whenever `/implement` runs inside a
  goal-driven loop that drives to completion without checking back.
- **`--resume <slug>`** — re-read `docs/_internal/plans/<slug>/_implementation.md`
  and pick up at the first incomplete phase. If `<slug>` was already given as
  the first positional argument, the flag is redundant but accepted.

The simplest valid call after a `DECOMPOSE-COMPLETE` handoff:

```
/implement figma-relay-recs --auto
```

Or hand-paced, one phase at a time:

```
/implement figma-relay-recs --phase 1
```

### Consuming the DECOMPOSE-COMPLETE handoff

When `/decompose` emits:

```
DECOMPOSE-COMPLETE
  status:      SOLIDIFIED
  plan_dir:    docs/_internal/plans/{slug}/
  phases:      {N}
  entrypoint:  docs/_internal/plans/{slug}/phase-01.md
  next:        implement phases 1..{N} in dependency order (e.g. /implement {slug} 1 --auto)
```

`/implement <slug> --auto` or `/implement <slug> --phase 1` is the correct
continuation. The slug, phase count, and plan directory are read directly from
`docs/_internal/plans/<slug>/README.md` — the handoff block need not be
reproduced; the plan directory is the authoritative source.

## Autonomous mode (`--auto`, or running under a goal)

When `--auto` is set — **or when this command runs inside a goal-driven loop
that drives to completion without checking back** — every human-in-the-loop
gate degrades to *"choose the sensible default, record the choice in
`_implementation.md`, and proceed."* Specifically:

| Interactive gate (default mode) | Autonomous behaviour |
|---|---|
| Echo plan summary, await confirmation before phase 1 | Adopt the plan as-is, log it, begin phase 1 |
| Escalation: spec contradicts ground truth | Conductor applies the conservative reading, logs it, flags for human review post-campaign; proceeds only if the risk is bounded — stops if the discrepancy could corrupt downstream phases |
| Phase reviewer returns NEEDS-WORK | Conductor dispatches a Sonnet fixer, re-runs the reviewer once; if still NEEDS-WORK, stop and surface |
| Architecture reviewer flags a seam defect | Conductor escalates to the human (hard stop — architecture calls cannot be automated away) |
| Window consolidation reveals a cross-phase break | Conductor patches both phases, re-runs the window reviewer once; if unresolved, stop |
| Existing `_implementation.md` present | **Resume** at the first incomplete phase |
| All phases SOUND | Emit `IMPLEMENT-COMPLETE` handoff and yield |

What **does not** change in autonomous mode: the scope boundary and the
escalation contract. The conductor never silently deviates from a frozen spec —
it escalates or stops. The disciplines in the **Rules** section are hard in
both modes; `--auto` removes human *gates*, never the *safety checks*.

In autonomous mode, the conductor emits a machine-readable handoff at
campaign completion so the outer goal loop can continue:

```
IMPLEMENT-COMPLETE
  status:      SHIPPED
  plan_dir:    docs/_internal/plans/{slug}/
  phases:      {N}
  ledger:      docs/_internal/plans/{slug}/_implementation.md
  next:        /verify {slug}   (mechanical drift check against plan + design)
```

## Phase 0: Parse arguments & resolve plan

Separate `<slug>`, `--phase`, `--auto`, `--resume` from `$ARGUMENTS`. If
`<slug>` is missing, reject: `Usage: /implement <slug> [--phase N] [--auto]`.

**Resolve the plan directory.** Read
`docs/_internal/plans/<slug>/README.md` — verify `status: SOLIDIFIED` is
stamped. If the status is not SOLIDIFIED, stop:
```
Plan <slug> is not SOLIDIFIED — run /decompose first.
```

**Resolve resume state.** If `_implementation.md` already exists, report the
current phase state (which phases are `shipped`, `in-progress`, `blocked`) and
determine the first incomplete phase. In `--auto` mode: resume silently, log
the decision. In default mode: surface current state and ask whether to resume
or restart.

**Echo the dispatch plan.** Before any subagent fires:
```
Implementing <slug>: {N} phases, starting at phase {start}.
  Plan:    docs/_internal/plans/<slug>/
  Ledger:  docs/_internal/plans/<slug>/_implementation.md
  Mode:    [single-phase N | end-to-end from N | end-to-end all]
```

In `--auto` mode: log and proceed. In default mode: await confirmation.

**Seed the ledger.** If `_implementation.md` does not exist, create it now
with the plan summary, agent roster, start timestamp, and an empty phase log.
The ledger is the source of truth for resume: every dispatch, gate, verdict,
escalation, and recovery must be recorded before the conductor moves on.

## Stage A — Dependency-ordered dispatch of implementer + reviewer pairs

For each phase in dependency order (phases whose dependencies are already
`shipped` or `deferred` are eligible; others wait):

### A1 — Opus Implementer

Dispatch one Opus `implementer` subagent per phase. Independent phases
(no unshipped dependencies on each other) may be dispatched concurrently.

Implementer dispatch prompt (per phase N):
> You are implementing **Phase N: {charter}** of {slug}.
> Read:
>   1. `docs/_internal/plans/{slug}/phase-NN.md` — your specification.
>   2. `docs/_internal/plans/{slug}/README.md` — global context and seams to
>      honour.
>   3. The phase docs for every phase you depend on — what they expose to you.
>   4. The actual codebase/artifacts at your slice: {slice}.
>
> Execute the specification **exactly as written**. If you find a contradiction
> between the spec and ground truth (a gate direction inverted, a routing key
> wrong, a file that does not exist at the stated path), **STOP immediately**
> and return an ESCALATION report — do not guess or patch around it. Otherwise:
> create/modify files as specified, follow existing code patterns, keep surface
> area minimal (only what the phase specifies), and run the project's tests if
> test infrastructure exists.
>
> When done, stamp `phase-NN.md` frontmatter:
> ```yaml
> ---
> phase: <N>
> feature: <slug>
> status: shipped              # in-progress if partial; blocked if hard-stopped
> completed: YYYY-MM-DD
> commit: <SHA or empty>
> files_touched:
>   - <path1>
>   - <path2>
> verification: pass | partial | deferred
> notes: <optional one-liner>
> ---
> ```
>
> Move any component superseded by this phase to a `_superseded/` subfolder
> within its directory and leave a one-line trace comment (filename + reason).
>
> Return: DONE with the list of files touched, or ESCALATION with the specific
> contradiction found. Do not return a partial DONE — if any step is incomplete,
> return ESCALATION or stamp `status: in-progress` and explain in notes.

**If the implementer returns ESCALATION:** the conductor reads the report,
verifies the contradiction on disk, applies the conservative reading or
requests human input (see Rules §1), logs the resolution in `_implementation.md`,
and re-dispatches the implementer with the clarified spec. This loop runs at
most twice without human input; on the third ESCALATION the conductor surfaces
to the human and halts.

### A2 — Opus Reviewer (paired, per phase)

Once the implementer returns DONE, dispatch one Opus `reviewer` subagent.

Reviewer dispatch prompt (per phase N):
> Review the **on-disk implementation of Phase N: {charter}** of {slug}.
> Do not rely on the implementer's self-report — **verify on disk directly**:
> read every file listed in `phase-NN.md` `files_touched`, re-run the phase's
> verification checklist yourself, re-derive any quoted counts or metrics.
>
> Check against:
>   - `phase-NN.md` — the specification (are all items delivered?).
>   - `docs/_internal/plans/{slug}/README.md` — seams to honour (does this
>     phase respect its declared interfaces to neighbours?).
>   - The actual codebase — patterns followed? Minimal surface? No accidental
>     complexity? No ornamental code? `_superseded/` used where appropriate?
>
> Adversarial verification where applicable: if a guard or gate was
> implemented, inject the failure condition and confirm it fires — do not just
> confirm it passes normally.
>
> Return: SOUND (phase is correct and complete), or NEEDS-WORK with **concrete
> edits** (not vibes) — specific file, line-range, and the exact change needed.
> The conductor applies NEEDS-WORK edits directly or dispatches a Sonnet fixer.

**If reviewer returns NEEDS-WORK:** the conductor applies the listed edits (or
dispatches a Sonnet `fixer` for mechanical changes) and re-runs the reviewer
**once**. If still NEEDS-WORK after one fix pass, the conductor stops and
surfaces to the human.

**A phase is `shipped` only after its reviewer returns SOUND.** Record the
pair verdict in `_implementation.md`.

### A3 — Architecture reviewer (for complex or seam-heavy phases)

For phases flagged in `README.md` as high-complexity or seam-heavy —
specifically: phases touching devops infrastructure, stateful session handling,
cross-phase integration contracts, or any phase the conductor judges as
carrying load-bearing architectural decisions — dispatch a **+1 Opus
architecture reviewer** operating in **extended thinking** (high reasoning).
This is a distinct agent from the line reviewer in A2.

Architecture reviewer dispatch prompt:
> You are doing a deep architectural review of **Phase N: {charter}** of
> {slug} — not a line review, an architecture review. Think carefully about:
>
>   1. **Contract correctness.** Does this phase's exported interface (the one
>      downstream phases depend on) match what `README.md` and the adjacent
>      `phase-NN.md` files declare? Name any discrepancy exactly.
>   2. **Load-bearing decisions.** Identify the 1–3 decisions in this phase
>      that, if wrong, would force rework in two or more later phases. Are they
>      correct? Are they reversible if not?
>   3. **Accidental coupling.** Does this phase introduce hidden state, side
>      channels, or leaked abstractions that the line reviewer might miss?
>   4. **Autonomy hazard.** Is there anything here the implementer may have
>      patched around (rather than escalating) that introduces a latent bug?
>
> Return: SOUND (architecture is clean), or ARCH-CONCERN with the load-bearing
> decision at risk, the specific file/line it lives in, and a proposed
> resolution. Do not return vague concerns — name exactly what breaks and where.

**If architecture reviewer returns ARCH-CONCERN:** this is a **hard stop** —
the conductor surfaces to the human regardless of `--auto` mode. Architecture
calls cannot be automated away.

## Stage B — Sliding-window seam review

Once a logical block of consecutive phases is all `shipped`, run Opus `seam-
reviewer` agents over **overlapping sliding windows** of `--window` consecutive
phases (default 3, mirroring `/decompose`'s Stage D): phases {1,2,3}, {3,4,5},
{5,6,7}, … — windows overlap by one so every boundary is double-covered.

Window reviewer dispatch prompt (per window):
> You are auditing the **on-disk implementation of phases {window}** of {slug}
> together. You are not re-reviewing a single phase — you are hardening the
> **seams between them** in the actual delivered code.
>
> Verify on disk (read the files, do not trust reports):
>   - Things that should be cohesive across these phases actually are — no
>     silent fragmentation of one concern across three files or modules.
>   - Things that should be decoupled actually are — no hidden coupling, no
>     leaked abstraction, no re-entrant state that was not declared.
>   - The interfaces between phases are exactly as `README.md` specified — no
>     undeclared additions, no quietly dropped contracts.
>   - No slop accumulated across phase boundaries: duplicate helpers,
>     inconsistent naming, copy-paste logic that should be shared.
>
> Return: SEAM-SOUND (all boundaries clean), or SEAM-DEFECT with the specific
> files and seam that is broken, and the concrete cross-phase fix needed.

**If SEAM-DEFECT:** the conductor patches the affected phases, updates their
`phase-NN.md` `files_touched`, and re-runs the window reviewer once. If still
SEAM-DEFECT, the conductor escalates to the human (hard stop in `--auto`).

Log every cross-phase decision in `_implementation.md`'s seam log section.

Optionally run a single final window reviewer over the **whole** campaign when
the phase count exceeds 2×window, to catch end-to-end coherence the rolling
windows cannot see.

## Stage C — Close out

When all phases are `shipped` and all window reviews are `SEAM-SOUND`, close
the campaign:

1. Stamp `docs/_internal/plans/{slug}/README.md` `status: SHIPPED`.
2. Record final summary in `_implementation.md`: phase count, dates, all agent
   verdicts, any escalations and how they were resolved.
3. Surface the campaign summary:

```
Implementation campaign complete for {slug} — SHIPPED.

  docs/_internal/plans/{slug}/
    README.md            plan + status (SHIPPED)
    phase-01.md … NN    each stamped shipped/deferred/blocked
    _implementation.md  orchestration ledger (all dispatches, verdicts, escalations)

Passes applied:
  A implementer+reviewer pairs   ({n} pairs, all SOUND)
  A+ architecture reviews        ({n} arch reviews, all SOUND)
  B sliding-window seam reviews  ({windows}, all SEAM-SOUND)

Next: /verify {slug}   (mechanical drift check against plan + design)
```

In **autonomous mode** (`--auto` / goal-driven), append the `IMPLEMENT-COMPLETE`
handoff block (see *Autonomous mode*) so the outer loop can continue — then
yield. In **default mode**, stop here. **Either way, do not run `/verify`
automatically** — verification is a separate step.

## Rules

1. **Escalate, don't improvise on a frozen contract.** If an implementer finds
   the spec contradicts ground truth (a gate direction inverted, a routing key
   wrong, a file at the wrong path), it STOPS and returns an ESCALATION — it
   does not guess or silently patch around it. The conductor verifies the
   discrepancy on disk, applies the conservative reading, and logs the
   resolution. This discipline caught an inverted barcode-subset gate and a
   routing-seam bug in the campaign that produced this methodology.

2. **Verify on-disk ground truth; never trust a report.** Reviewers re-run
   tests themselves, re-derive counts, re-parse artifacts from the file system.
   The conductor spot-checks claimed outputs before recording SOUND verdicts.
   This discipline recovered two session-limit truncations — one hiding a
   half-finished run and one masking a live zombie process.

3. **Adversarial verification.** To confirm a guard or gate works correctly,
   INJECT the failure condition (e.g. malformed input, missing file, inverted
   flag) and prove the guard fires — do not just confirm it passes on the happy
   path. Independent re-derivation over any quoted metrics or counts.

4. **The ledger is the source of truth.** `_implementation.md` records every
   dispatch, gate verdict, escalation, recovery, and architectural decision
   before the conductor moves on. A campaign interrupted at any point must be
   resumable from this file alone, without relying on conversation history.

5. **Surface to the human only at critical inflection points.** These are:
   (a) an ESCALATION the conductor cannot resolve by conservative reading;
   (b) an ARCH-CONCERN from the architecture reviewer;
   (c) a NEEDS-WORK or SEAM-DEFECT that survives a single fix pass.
   All other decisions are made by the conductor, logged, and proceeded.
   Surfacing on every phase transition is noise — reserve it for calls only
   the human meta-architect can make.

6. **Move superseded components to `_superseded/` subfolders.** When a phase
   replaces a prior implementation (a script, a config, a handler), the old
   version moves to `_superseded/<original-name>` within its directory. A
   one-line trace comment (filename, the phase that superseded it, the reason)
   is left in the `_superseded/` folder as `README.md` or inline. Never delete
   silently.

7. **Opus throughout for judgment; Sonnet for mechanical fixes only.** Every
   implementer, reviewer, architecture reviewer, and seam reviewer is `model:
   opus`. A Sonnet `fixer` may be dispatched by the conductor to apply
   concrete mechanical edits listed by a reviewer — but only when the edits are
   fully specified and require no judgment (exact file, line-range, replacement
   text). If the fix requires any architectural reasoning, re-dispatch Opus.

8. **One implementer owns one phase.** Dedicated per-phase implementers are
   the point — do not let one agent implement several phases; context is lost
   and seam errors multiply. The only exception is a pair of trivially
   coupled micro-phases (e.g. a one-line config change and its test) — these
   may be bundled at the conductor's discretion, recorded in `_implementation.md`.

9. **Minimal surface area.** Implementers do not add "while we're here"
   refactors, speculative helpers, or defensive code for scenarios outside the
   phase spec. If something looks like it should change, they escalate; they do
   not silently improve. The north-star from the plan's README is a gate, not a
   garnish.

10. **The scope boundary is hard in both modes.** `/implement` never plans or
    re-designs phases. If a phase doc is architecturally broken, the implementer
    escalates; the conductor surfaces to the human. Under a goal, `/implement`
    produces shipped artifacts and yields an `IMPLEMENT-COMPLETE` handoff; it
    does not re-enter `/decompose` territory, and it does not idle waiting for
    a human who isn't coming.
