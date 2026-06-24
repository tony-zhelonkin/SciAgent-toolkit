# Pipeline Plan

Plan AND orchestrate a multi-phase analysis pipeline with explicit model-tiering and a runnable acceptance gate. Opus decomposes a scope-doc into bite-sized phases (one phase == one script == one implementer); Sonnet implementers execute them dependency-ordered; Opus reviews every N phases by **actually running the script and asserting its artifacts exist**.

This command **composes on top of** `/plan`/`/implement`/`/verify` — it does NOT replace them. It adds science-flavored decomposition (the 14839 gold-standard INDEX/phase format), model-tiering (Opus plan, Sonnet implement, Opus review), and a runnable+artifacts gate the software commands lack. Where per-phase mechanics fit the software flow, delegate to `/implement`; the drift check to `/verify`.

## Phase 0: Parse arguments

| Param | Shape | Default | Meaning |
|---|---|---|---|
| `slug` | `$ARGUMENTS[0]` (positional, **required**) | — | Plan slug. Resolves the plan dir to `docs/_internal/plans/{today}-{slug}/` (date = `date +%F`). |
| `--scope-doc <path>` | flag | latest `docs/_internal/research/{date}-{slug}/` synthesis (prefer `_SYNTHESIS.md`, else newest `*.md`) | The synthesis/context the Opus planner reads to decompose. |
| `--n-planners <int>` | flag | `1` | Number of parallel Opus decomposers. `>1` fan-outs alternative decompositions; you then pick/merge into one INDEX. |
| `--context-budget <int>` | flag | `35` | Target % of a 200k Sonnet context per phase == one script == one bounded brief. Drives how finely the planner splits phases. |
| `--review-every <int>` | flag | `3` | Opus REVIEW CHECKPOINT cadence, in phases. A checkpoint row is inserted after every this-many substantive phases. |
| `--background <auto\|on\|off>` | flag | `auto` | Long non-blocking phases run in monitorable tmux (`run_in_background`). `auto` = decide per phase from `depends_on`; `on` = always background long phases; `off` = strictly sequential. |

Flag parsing is order-independent. If `slug` is missing, reject:
```
Usage: /pipeline-plan <slug> [--scope-doc <path>] [--n-planners <int>] [--context-budget <int>] [--review-every <int>] [--background auto|on|off]
```

Resolve and announce the run config:
```
/pipeline-plan {slug}
  Plan dir:        docs/_internal/plans/{date}-{slug}/
  Scope-doc:       {resolved path}
  Planners:        {n-planners} (Opus)
  Context budget:  {context-budget}% per phase
  Review cadence:  every {review-every} phases (Opus)
  Background:      {background}
```

**STOP condition — no scope-doc.** If no `--scope-doc` is given and no `docs/_internal/research/{date}-{slug}/` synthesis exists, do NOT invent context. Stop:
```
No scope-doc found at docs/_internal/research/{date}-{slug}/ and none passed via --scope-doc.
A pipeline plan must be grounded in a synthesized scope. Run:
  /explore-and-plan "<question>" {slug}
to produce the research fan-out + synthesis, which hands off to /pipeline-plan automatically.
(Or pass --scope-doc <path> to an existing synthesis.)
```
Stop.

## Phase 1: Plan — Opus decomposition (model tier: **Opus**)

Dispatch `--n-planners` Opus planner(s). Each planner:

1. **Reads** (do not re-explore beyond these):
   - the resolved `--scope-doc` (the verified-vs-inferred synthesis is the source of truth for *what* to build),
   - the repo's `AGENTS.md` — **including the `SCIAGENT:CRAFT` managed block** (the five standing conventions: figures, results placement, README adjacency, planning decomposition, reproducibility),
   - the **figure-style contract** (`skills/figure-style/SKILL.md` + `lib/figure-style/figure_helpers.{R,py}`) — so each phase declares COMPUTE-ONLY / VIZ-ONLY / MIXED correctly,
   - `02_analysis/config/analysis_config.yaml` (`stages:`, `figures:`) — to ground stage-ids and floors.

2. **Writes the plan FROM THE P16 TEMPLATES** into `docs/_internal/plans/{date}-{slug}/`:
   - `00_INDEX.md` from `templates/plan/00_INDEX.md.template`,
   - one `NN_<slug>.md` per phase from `templates/plan/NN_slug.md.template`.

3. **Sizes each phase** so a single Sonnet implementer stays **≤ `--context-budget`%** of a 200k context == **one script == one bounded brief**. A phase that needs more than that is too big — split it. Compute and viz are separate phases (compute never plots, viz never computes).

4. **Fills the INDEX phase table** with columns `seq | slug | title | tier | concern | depends_on`. Substantive phases are `tier: Sonnet`; open-judgement/design phases are `tier: Opus`. `depends_on` lists prior `seq` ids (comma-separated, `—` for none).

5. **Inserts an Opus REVIEW CHECKPOINT row every `--review-every` phases** (e.g. `R1`, `R2`, …), `tier: Opus`, `concern: review`, `depends_on` = the seq ids it gates. Each checkpoint row's title states the gate verbatim: *plan adherence + cleanliness/namespace + tests green + the phase's `03_results/` artifacts exist and are non-empty*.

6. Each `NN_<slug>.md` fills all fixed sections of the template: figure-style contract declaration; §1 Objective; §2 Scope (with out-of-scope greps); §3 Inputs (citing the scope-doc precisely — the implementer must act without re-exploring); §4 Outputs (exact `03_results/<stage>/...` artifact paths); §5 Implementation (exact helper names — `save_overview()`, `contrast_path()`, …); §6 Captions; §7 Acceptance checks (grep + structural, runnable verbatim); §8 Gotchas.

If `--n-planners > 1`: collect the alternative decompositions, then (as Opus) reconcile into **one** authoritative INDEX + phase set — do not leave competing plans on disk.

After the INDEX is written, **surface the phase table for a glance-check** (wrong phase boundary catch), then proceed to execution — do NOT block indefinitely on approval.

## Phase 2: Execute — dependency-ordered Sonnet implementers (model tier: **Sonnet**)

Walk the INDEX phase table in dependency order (topological by `depends_on`). For each substantive (non-review) phase:

1. Dispatch **one Sonnet implementer** for that phase. Delegate the per-phase mechanics to `/implement` where natural (the phase brief is the contract; the implementer creates exactly the script + artifacts the brief declares, following `map.md`/existing patterns, minimal surface area).
2. **Respect `depends_on`.** A phase does not start until every phase in its `depends_on` is complete (its declared artifacts exist).
3. **Background overlap (`--background`).** A phase that is long-running AND does not block its successors launches in tmux via `run_in_background` (monitorable). Later phases whose `depends_on` is already satisfied proceed in parallel — dependency-aware overlap. Under `--background off`, run strictly sequentially. Under `auto`, background only phases with no un-started dependents.
4. Stamp each phase's `NN_<slug>.md` frontmatter `status:` as the implementer finishes (same canonical schema `/implement` uses), so the reviewer and `/verify` can read completion state.

When the running phase count reaches a multiple of `--review-every` (i.e. the next `RN` checkpoint's dependencies are all complete), pause new dispatch and run Phase 3 for that checkpoint before continuing.

## Phase 3: Review — Opus, every `--review-every` phases — THE REAL GATE (model tier: **Opus**)

At each `RN` checkpoint, dispatch **one Opus reviewer** over the phases the checkpoint gates. The reviewer MUST do all three:

> **(a) Adherence + hygiene.** Confirm each gated phase did what its `NN_<slug>.md` brief specified — plan adherence, code cleanliness, namespace separation (the phase's identifiers/stage-ids are grep-isolable and disjoint from peers), and no drift/rot (no out-of-scope edits, no second homes, no dead code left behind).
>
> **(b) Runnable + artifacts gate — RUN THE SCRIPT.** For each gated phase, **actually execute the phase's `02_analysis/scripts/NN_*` script (or confirm its committed/logged run), then assert every artifact the phase's §4 Outputs declares under `03_results/` exists AND is non-empty.** A phase whose script does not run, or whose declared artifacts are absent or empty, FAILS the checkpoint — regardless of how good the diff looks. This runnable+artifacts assertion is the thing the owner otherwise re-types by hand; it is the point of the review.
>
> **(c) Persist the review.** Write the review verdict + evidence (commands run, artifact `ls`/size output, any failures) to `docs/_internal/reasoning/{date}_RN_{checkpoint-slug}.md`. A review that lives only in chat is non-reproducible.

If a checkpoint fails, STOP dispatch of downstream phases and surface the failure with the exact phase, the missing/empty artifact, and the script command that failed. Do not paper over it by advancing.

## Phase 4: Output summary

When all phases and checkpoints are green:
```
/pipeline-plan {slug} complete.

Plan dir: docs/_internal/plans/{date}-{slug}/
  00_INDEX.md
  NN_<slug>.md  × {phase_count} substantive phases
  Review checkpoints: {checkpoint_count} (every {review-every} phases)

Model tiering:  planner = Opus · implementers = Sonnet · reviewers = Opus
Reviews persisted to: docs/_internal/reasoning/{date}_RN_*.md

Per-phase:
  01 {slug} ✅ — {script} ran, artifacts present
  02 {slug} ✅ — ...
  R1 ✅ — adherence + tests green + artifacts non-empty
  ...

Next: /verify {slug}   (mechanical drift check: plan + scope-doc vs current code)
```

If stopped at a failed checkpoint, report the failing phase, the missing/empty artifact path, and the failed script command instead, and recommend the fix-then-re-review loop.

## Rules

1. **Model tiering is explicit, throughout.** Planner = **Opus** (open-judgement decomposition). Implementers = **Sonnet** (bounded per-phase edits). Reviewers = **Opus** (the gate). State the tier at every dispatch — never silently downgrade the reviewer.
2. **The review gate RUNS code.** A checkpoint is not satisfied by a grep or a clean-looking diff. It is satisfied only when the phase's script runs and its declared `03_results/` artifacts exist and are non-empty. This is non-negotiable — it is the convention the owner re-types.
3. **One phase == one script == one implementer (~`--context-budget`%).** If a phase exceeds the budget, split it. Compute and viz are always separate phases (compute never plots; viz never computes).
4. **Composition, not replacement.** `/pipeline-plan` wraps `/plan`/`/implement`/`/verify` — it adds the science INDEX/phase templates, model-tiering, and the runnable+artifacts gate. It delegates per-phase mechanics to `/implement` and the final drift check to `/verify`. The software architect pipeline is untouched.
5. **Fill the templates; don't free-form.** The plan IS the P16 templates (`templates/plan/00_INDEX.md.template` + `NN_slug.md.template`) filled — same fixed section order, same `seq|slug|title|tier|concern|depends_on` table. Reviewers grep section headers.
6. **Ground in the scope-doc; never invent context.** Missing scope-doc → stop and recommend `/explore-and-plan`. Phase briefs cite the scope-doc precisely so implementers don't re-explore.
7. **Persist every review.** Reviews go to `docs/_internal/reasoning/`. A decision or verdict with no trace is non-reproducible (CRAFT reproducibility rule).
8. **Respect dependencies; overlap only when safe.** Background/tmux overlap is an optimization gated on `depends_on` and `--background`; on harnesses without background semantics it degrades to sequential, and that is fine.
