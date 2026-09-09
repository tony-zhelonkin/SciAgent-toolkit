# Pipeline Plan

## What this is for

Turn a synthesized scope into a phased analysis pipeline, then run it with a
gate that **executes each stage and asserts its artifacts exist**.

That gate is the whole reason the command exists. A plan can be followed, the
diff can look right, and the stage can still not run — or run and write nothing.
Checking that by hand, phase after phase, is the work this automates. Everything
else here is scaffolding around it.

## When to use

- A question has been researched and synthesized, and the work ahead is several
  stages long with real dependencies between them.
- You want the phase boundaries written down before implementation starts.

## When not to use

- The scope is not yet synthesized → run `/explore-and-plan` first. This command
  refuses to invent context.
- The work is one stage → write the stage. A plan directory for a single phase is
  overhead.
- The change is to software rather than an analysis pipeline → the
  `architecture-first-dev` skill routes design, planning and verification for
  that track.

## What it writes

```
docs/_internal/_project/plans/{date}-{slug}/
├── 00_INDEX.md          # the phase table: seq | slug | title | tier | concern | depends_on
└── NN_<slug>.md         # one bounded brief per phase
docs/_internal/_project/{checkpoint-slug}-review.md   # one per review checkpoint
```

Plus the stages and artifacts the phases declare. Nothing outside the declared
paths.

## When arguments are missing

`slug` is the only required argument. If it is absent, **ask what the pipeline is
for** and propose a slug from the answer. Do not print a usage line — the caller
is in a conversation, and the thing actually missing is intent, not syntax.

The same applies to a missing scope document, with one difference: that one is a
stop, not a question. If neither `--scope-doc` nor
`docs/_internal/_project/{slug}-synthesis.md` exists, say so and recommend
`/explore-and-plan "<question>" {slug}`. A plan built on invented context wastes
every phase downstream of it.

## Arguments

| Param | Shape | Default | Meaning |
|---|---|---|---|
| `slug` | positional | ask | Plan slug. Resolves the plan dir under `docs/_internal/_project/plans/{today}-{slug}/`. |
| `--scope-doc <path>` | flag | `docs/_internal/_project/{slug}-synthesis.md` | The synthesis the planner decomposes. The source of truth for *what* to build. |
| `--n-planners <int>` | flag | `1` | Parallel decompositions. Above 1, reconcile into one authoritative INDEX; competing plans must not be left on disk. |
| `--context-budget <int>` | flag | `35` | Target share of one implementer's context per phase. Drives how finely phases split. |
| `--review-every <int>` | flag | `3` | Review checkpoint cadence, in phases. |
| `--background <auto\|on\|off>` | flag | `auto` | Overlap long phases whose successors do not wait on them. `auto` decides per phase from `depends_on`. Degrades to sequential where the harness has no background semantics. |

Flag order does not matter. Announce the resolved configuration before starting.

## Phase 1 — decompose

Dispatch the planner(s) at the open-judgement tier. Each one reads only:

- the resolved scope document;
- the repository's `AGENTS.md`, including the `SCIO:CRAFT` managed block;
- the `figure-style` skill and its helper library, so each phase declares
  COMPUTE-ONLY / VIZ-ONLY / MIXED correctly;
- `02_analysis/config/analysis_config.yaml` (`stages:`, `figures:`), to ground
  stage ids and floors.

It writes the plan from `templates/plan/00_INDEX.md.template` and
`templates/plan/NN_slug.md.template` — same section order, same table columns,
because reviewers grep the headings.

Sizing rule: **one phase == one stage == one implementer**, within
`--context-budget`. A phase that needs more is too big; split it. Compute and viz
are always separate phases.

Insert a review checkpoint row every `--review-every` substantive phases, with
`depends_on` naming the phases it gates and a title stating the gate verbatim.

Surface the phase table for a glance-check — a wrong phase boundary is cheapest
to catch here — then proceed. Do not block indefinitely on approval.

## Phase 2 — execute

Walk the phase table in dependency order. For each substantive phase, dispatch
one implementer with the phase brief as its contract, delegating per-phase
mechanics to `/implement` where that fits. A phase starts only once every phase
in its `depends_on` has produced its declared artifacts. Stamp each brief's
`status:` frontmatter as its implementer finishes.

When the next checkpoint's dependencies are all complete, pause new dispatch and
run Phase 3 before continuing.

## Phase 3 — review, and the gate

At each checkpoint, one reviewer at the open-judgement tier does three things:

**(a) Adherence.** Each gated phase did what its brief specified: no out-of-scope
edits, no second homes, no dead code, and identifiers grep-isolable from peers.

**(b) Run the stage.** Execute each gated phase's `02_analysis/stages/NN_*`, or
confirm a committed and logged run, then assert every artifact its §4 Outputs
declares exists and is non-empty. **A phase whose stage does not run, or whose
declared artifacts are absent or empty, fails — however good the diff looks.**

**(c) Persist the verdict.** Write it with its evidence — commands run, artifact
listings and sizes, failures — to
`docs/_internal/_project/{checkpoint-slug}-review.md`. A review that lives only
in chat is not reproducible.

A failed checkpoint stops downstream dispatch. Report the phase, the missing or
empty artifact, and the command that failed. Do not advance past it.

## Phase 4 — report

Per phase: ran / artifacts present. Per checkpoint: the verdict and where it was
written. Then the mechanical drift check between the plan, the scope document and
the code as it now stands.

If stopped at a failed checkpoint, report that instead, and the fix-then-re-review
loop.

## Rules

1. **The gate runs code.** Not a grep, not a clean diff. The stage runs and its
   declared artifacts exist and are non-empty.
2. **One phase, one stage, one implementer.** Over budget means split. Compute
   never plots; viz never computes.
3. **Fill the templates.** The plan *is* those templates filled, in their section
   order.
4. **Ground in the scope document.** No scope, no plan. Briefs cite it precisely
   enough that an implementer does not re-explore.
5. **Persist every review.** A verdict with no trace is not reproducible.
6. **Overlap only where dependencies allow it**, and treat sequential fallback as
   correct rather than degraded.
7. **State the tier at every dispatch.** Never silently drop the reviewer to a
   cheaper tier — the review is the product.
