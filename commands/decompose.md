# Decompose

Orchestrate a **multi-agent, multi-pass planning campaign** over a large,
multi-artifact problem — several codebases, long documents, mixed scripts and
specs — and emit a solidified, architecturally-validated **phase-plan document
set**. The output is a plan ready to hand to implementers; **implementation is
explicitly out of scope for this command.**

Use this instead of the architect skill's single-pass plan route when the problem
is too big for one agent to hold:
when decomposition itself needs review, when each phase deserves its own
dedicated planner working against its own slice of the artifacts, and when the
seams between phases must be audited before a single line is written.

Every subagent dispatched here is an **Opus** subagent (`model: opus`). The main
agent is the conductor: it dispatches, gates, and records — it does not itself
plan the phases.

## Usage

```
/decompose <free-form task statement>                  — just describe the job
/decompose <task...> --north-star '<rubric or path>'   — supply a design rubric
/decompose <task...> --artifacts <pathA,pathB,...>     — pin the in-scope codebases/docs
/decompose <task...> --slug <kebab-id>                  — name the output folder yourself
/decompose <task...> --window 3                        — sliding-window size (default 3)
/decompose <task...> --auto                             — autonomous: no human gates, drive to SOLIDIFIED
/decompose --resume <slug>                             — continue an in-flight campaign
```

Everything you type after `/decompose` is the **task statement** — write it in plain
prose, no quoting required. Flags (`--north-star`, `--artifacts`, `--slug`,
`--window`, `--resume`) may appear anywhere; the prose around them is the task.
You do **not** invent a slug — it's derived from the task and echoed back for
confirmation. So the simplest valid call is just one sentence:

```
/decompose Address 100% of the recommendations in docs/testing/claude-code-interaction
           against the figma-relay codebase --north-star 'laconic, agent-first Figma
           controls for a blind operator'
```

- **task statement** (required) — all the non-flag prose in `$ARGUMENTS`. Describe
  what to plan; plain text, no quotes needed.
- `--slug <kebab-id>` — optional explicit name for the `docs/_internal/plans/<slug>/` output
  folder. **If omitted, derive a kebab-case slug from the task statement and echo
  it back** ("planning into `docs/_internal/plans/figma-relay-recs/` — ok?") before any work.
- `--artifacts <csv>` — the large artifacts in scope (repos, dirs, design docs).
  If omitted, the pioneer infers them and **lists what it chose for confirmation.**
- `--north-star <path|quoted>` — a domain-specific design rubric that the
  populate and consolidate passes optimise against. Defaults to the **House
  rubric** below.
- `--window N` — consecutive-phase window for the consolidation pass (default 3).
- `--auto` — **autonomous mode** (see below). Implied whenever `/decompose` is
  invoked under a goal-driven / drive-to-completion run.
- `--resume <slug>` — re-read `docs/_internal/plans/<slug>/` and pick up at the first
  incomplete stage (slug required here so the campaign is unambiguous).

## Autonomous mode (`--auto`, or running under a goal)

When `--auto` is set — **or when this command runs inside a goal-driven loop that
drives to completion without checking back** — every human-in-the-loop gate
degrades to *"choose the sensible default, record the choice in `_campaign.md`,
and proceed."* Specifically:

| Interactive gate (default mode) | Autonomous behaviour |
|---|---|
| Echo derived slug, await confirmation | Adopt the derived slug, log it, continue |
| Pioneer lists inferred artifacts for confirmation | Adopt the inferred set, log it; if truly nothing is in scope, that's a hard FAIL — stop |
| Stage B resplit "glance-check" before proceeding | Apply the resplit, log it, continue (still bounded to **two** Stage-B passes) |
| Existing `docs/.../plan/` → ask resume-or-restart | **Resume** at the first incomplete stage |
| Stage E "surface and stop" | Surface, then emit the handoff signal below — do not idle |

What **does not** change in autonomous mode: the hard scope boundary. `/decompose`
still ends at a `SOLIDIFIED` plan and **never crosses into implementation**, even
under a goal. Instead of silently halting it yields control back with a
machine-readable handoff so the outer goal loop can pick up the implementation:

```
DECOMPOSE-COMPLETE
  status:      SOLIDIFIED
  plan_dir:    docs/_internal/plans/{slug}/
  phases:      {N}
  entrypoint:  docs/_internal/plans/{slug}/phase-01.md
  next:        implement phases 1..{N} in dependency order (e.g. /implement {slug} 1 --auto)
```

The two-pass Stage-B cap and the "is anything in scope?" check are retained as
**hard stops** even in autonomous mode: a goal-driven loop should fail loudly on a
malformed decomposition rather than burn the whole budget on a bad skeleton.

All artifacts land under `docs/_internal/plans/{slug}/`.

## The House rubric (default north-star)

Every populate and consolidate agent optimises the plan not only for **function**
but for **structure, architecture, maintainability, and the beauty of minimally
necessary complexity.** The tool being planned is, by default, held by a *blind*
operator — an agent. So the interface the plan carves out must read as **laconic,
ergonomic, and effortless to wield**: maximum control with minimum surface,
nothing accidental, nothing ornamental. Couple what is cohesive; decouple what is
separable; spend complexity only where the problem genuinely demands it. When the
user supplies `--north-star`, that rubric **augments** (never silently replaces)
this one — cite both.

## Phase 0: Parse arguments & resolve scope

Separate the **task statement** (all non-flag prose) from the flags
(`--slug`, `--artifacts`, `--north-star`, `--window`, `--resume`). If the task
statement is empty and `--resume` was not passed, ask what to plan.

**Resolve the slug.** If `--slug` was given, use it. Otherwise derive a kebab-case
slug from the task statement (e.g. "Address the figma-relay recommendations" →
`figma-relay-recs`) and **echo it for confirmation** before doing any work — the
user may rename it in their next message. On `--resume <slug>`, the slug is given.

Create `docs/_internal/plans/{slug}/` if absent. If it already exists and `--resume` was not
passed, report current stage state (which `phase-NN.md` exist, which gates passed
per `_campaign.md`) and ask whether to **resume** or **restart**.

## Stage A — Pioneer decomposition (1–2 Opus, parallel-independent)

Dispatch **one or two** Opus `pioneer` subagents. Two only when the task statement
visibly contains **multiple independent sub-problems** — in which case each pioneer
owns one sub-problem and decomposes it **recursively** (a sub-problem with its own
internal structure gets its own nested phase list). One pioneer otherwise.

Pioneer dispatch prompt (per pioneer):
> You are decomposing a large multi-artifact problem into an **ordered list of
> implementation phases** — you are NOT writing the plans, only the skeleton.
> Read the task statement and reconnoitre the in-scope artifacts: {artifacts}.
> Produce an ordered phase list where **each phase is a unit a single Opus agent
> can fully plan within one context window** — neither a one-liner nor a
> subsystem. For each phase give: a title, a one-sentence charter, the specific
> artifact slices it touches, and its hard dependencies on earlier phases. If the
> problem is itself several sub-problems, decompose each recursively and label the
> nesting. Surface the **seams** you expect between phases (what must stay
> cohesive across a boundary; what must stay decoupled). Return structured phase
> stubs only — no populated plans. Optimise the decomposition against this
> north-star: {north-star}.

The main agent merges pioneer output (if two ran, reconcile overlapping phases and
renumber into one global order) and writes:

- `docs/_internal/plans/{slug}/README.md` — overview, the **Phase Summary** table (Phase |
  Charter | Artifacts | Depends-on), the resolved artifact set, and the active
  north-star (House + any `--north-star`).
- `docs/_internal/plans/{slug}/phase-NN.md` — one **stub** per phase: charter, artifact
  slice, dependencies, expected seams. Body left for Stage C to populate.
- `docs/_internal/plans/{slug}/_campaign.md` — the orchestration ledger: stage gates, agent
  roster, and the decomposition rationale. Seed it with Stage A's verdict.

## Stage B — Feasibility gate (1 Opus reviewer)

Dispatch one Opus `feasibility` reviewer over the full stub set + README.

> Review this phase decomposition for **single-agent feasibility**. For each
> phase, judge whether one Opus agent could populate a complete, high-quality
> plan for it within one context window given its artifact slice. Flag phases
> that are **too large** (recommend a split, with the split boundary) or **too
> thin** (recommend a merge with a named neighbour). Check the dependency graph is
> acyclic and the ordering is buildable. Check the seams the pioneer named are the
> real ones. Return a verdict per phase plus an overall READY / NEEDS-RESPLIT.

If `NEEDS-RESPLIT`: apply the recommended splits/merges, renumber, update README
+ stubs + `_campaign.md`, and re-run Stage B **once**. Surface the resplit to the
user as a glance-check before proceeding. Do not loop more than twice without user
input.

## Stage C — Populate, in planner+reviewer pairs (Opus, pipelined)

Each phase is planned by its **own** Opus agent — this is deliberate: a dedicated
planner per phase keeps maximum control over its stage. Then each populated phase
is reviewed by a paired Opus reviewer. Planner and reviewer move as a **pair, phase
by phase**; pairs for independent phases may run concurrently, but a phase whose
dependency is still unpopulated waits for it.

**Planner** dispatch prompt (per phase N):
> Populate the full implementation plan for **Phase N: {charter}** of {slug}.
> Read the phase stub, the README (for global context and the seams you must
> honour), the plans of phases you depend on, and **your phase's artifact slice
> against the actual codebase/docs**: {slice}. Write `phase-NN.md` with: Context
> (what upstream phases make available), Files to Create / Modify (with the real
> paths and the patterns to follow, grounded in the artifacts), the interfaces
> this phase exposes downstream, Verification checks, and References to the
> artifacts. Plan ONLY this phase — do not redesign neighbours; honour their
> seams. Optimise against the north-star: {north-star}. Emit a plan, not code.

**Reviewer** dispatch prompt (paired, per phase N):
> Review the populated `phase-NN.md` against the **actual artifacts** ({slice})
> and the README seams. Trim nicks: paths that don't exist, patterns that don't
> match the codebase, interfaces that contradict a neighbour, accidental
> complexity, anything ornamental. Confirm the phase is architecturally sound in
> isolation and at its declared boundaries. Return concrete edits (not vibes); the
> main agent applies them. Verdict: SOUND / NEEDS-WORK with specifics.

Apply each reviewer's edits to its `phase-NN.md`. Record per-pair verdicts in
`_campaign.md`. A phase is "populated" only after its reviewer returns SOUND.

## Stage D — Sliding-window consolidation (Opus, overlapping windows)

Once all phases are populated, run Opus `consolidator` agents over **overlapping
sliding windows** of `--window` consecutive phases (default 3): phases {1,2,3},
{3,4,5}, {5,6,7}, … — windows overlap by one so every boundary is double-covered.
These are the meta-agents that make the plan *robust*.

Consolidator dispatch prompt (per window):
> You are auditing phases {window} **together** against the artifacts and the
> overall architecture. You are not re-planning a single phase — you are
> hardening the **seams between them**. Verify: things that should be cohesive
> across these phases actually are (no silent fragmentation of one concern across
> three phases); things that should be decoupled actually are (no hidden coupling,
> no leaked abstraction). Deepen the architectural thinking — name the load-bearing
> decisions and check each phase honours them. Hunt accidental complexity across
> the window and propose the simpler shape where one exists. Judge the window
> against the north-star: {north-star} — is the interface this set carves out
> laconic, ergonomic, beautiful to wield for a blind agent operator? Return
> cross-phase edits (which `phase-NN.md`, what change, why) plus a window verdict.

Apply window edits. Where two overlapping windows touch the same shared phase,
reconcile their edits (the shared phase is the overlap by design — prefer the edit
that strengthens the seam). Log every cross-phase decision in
`docs/_internal/plans/{slug}/_consolidation.md` (the seam ledger: coupling/decoupling calls,
complexity-budget verdicts, load-bearing decisions). Update `_campaign.md` gate.

Optionally run a single final consolidator over the **whole** plan when the phase
count exceeds 2×window, to catch end-to-end coherence the windows can't see.

## Stage E — Surface the solidified plan

The plan is now multi-pass validated. Stamp `docs/_internal/plans/{slug}/README.md` with
`status: SOLIDIFIED` and surface:

```
Phase-plan campaign complete for {slug} — SOLIDIFIED.

  docs/_internal/plans/{slug}/
    README.md          overview + phase summary + north-star
    phase-01.md … NN   populated, pair-reviewed, window-consolidated
    _decomposition…    pioneer rationale (in _campaign.md)
    _consolidation.md  seam ledger (coupling/decoupling, complexity calls)
    _campaign.md       orchestration ledger (every gate + verdict)

Passes applied:
  A pioneer decompose   ({n} pioneers)
  B feasibility gate    ({READY|resplit×k})
  C populate pairs      ({n} planner+reviewer pairs, all SOUND)
  D sliding windows     ({windows}, all verdicts)

This plan is ready to hand to implementers.
Implementation is out of scope for /decompose — kick it off separately
(e.g. /implement {slug} <phase> under the architect role).
```

In **autonomous mode** (`--auto` / goal-driven), append the `DECOMPOSE-COMPLETE`
handoff block (see *Autonomous mode*) so the outer loop can continue into
implementation — then yield. In **default mode**, stop here. **Either way, do not
begin implementing** — the scope boundary is hard in both modes; only the *exit
behaviour* differs (idle-and-wait vs. yield-with-handoff).

## Rules

1. **Planning only.** This command never writes product code. If a subagent
   starts implementing, stop it — its job is the plan.
2. **Opus throughout.** Every dispatched subagent is `model: opus`. The main agent
   conducts; it does not author phase plans itself.
3. **One planner owns one phase.** Dedicated per-phase planners are the point —
   do not let one agent populate several phases; control is lost.
4. **Pairs, then windows.** A phase is reviewed in isolation (Stage C) *and* at its
   seams (Stage D). Both gates are mandatory before SOLIDIFIED.
5. **Windows overlap.** Never tile windows edge-to-edge; the one-phase overlap is
   what double-covers every boundary.
6. **Ground every claim in artifacts.** Planners and reviewers cite real paths and
   real patterns. A plan that references files that don't exist is NEEDS-WORK.
7. **The ledger is the source of truth.** `_campaign.md` records every gate and
   verdict so `--resume` can pick up deterministically and the user can audit the
   campaign after the fact.
8. **Minimally necessary complexity.** Spend complexity only where the problem
   demands it; the north-star is a gate, not a garnish.
9. **The scope boundary is hard in both modes.** Autonomous mode removes the human
   *gates*, never the *planning-only* contract. Under a goal, `/decompose` produces
   the plan and yields a `DECOMPOSE-COMPLETE` handoff; it does not implement, and it
   does not idle waiting for a human who isn't coming.
