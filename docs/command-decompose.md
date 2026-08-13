# Decompose command

A single heavy-orchestration command for the **planning phase** of large,
multi-artifact problems — several codebases, long documents, mixed scripts and
specs — where decomposition itself needs review and each phase deserves its own
dedicated planner.

| Command | What it does |
|---|---|
| `/decompose <task...>` | Run a multi-agent, multi-pass planning campaign → `docs/_internal/plans/{slug}/phase-NN.md` (SOLIDIFIED). Slug auto-derived from the task; `--north-star`, `--artifacts`, `--slug`, `--window`, `--resume` optional. |

## How it differs from `/plan`

`/plan` (architect pipeline) is a **single-agent** decomposition of an
already-approved design into phase files. `/decompose` is a **multi-agent
campaign** for problems too big for one agent to hold:

```
A  pioneer decompose      1–2 Opus, recursive on multi-problem inputs
B  feasibility gate       1 Opus — is each phase single-agent-sized?
C  populate (pairs)       Opus planner + Opus reviewer, per phase
D  sliding-window         Opus consolidators over overlapping 3-phase windows
E  surface                SOLIDIFIED plan, ready to hand to implementers
```

**Implementation is out of scope.** `/decompose` ends at a validated plan; run
implementation separately (e.g. `/implement` under the architect role).

**Goal-composable.** Under a goal-driven / drive-to-completion run (or `--auto`),
every human gate degrades to a logged default and the command ends with a
`DECOMPOSE-COMPLETE` handoff (plan dir + entrypoint) instead of idling — so the
outer goal loop can pick up implementation from the durable plan. The
planning-only scope boundary stays hard in both modes.

Every dispatched subagent is an Opus subagent; the main agent conducts.
