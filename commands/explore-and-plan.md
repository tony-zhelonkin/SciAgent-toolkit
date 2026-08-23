# Explore and Plan

## What this is for

Research a question across several lanes at once, reconcile the lanes into one
synthesis, and hand that synthesis to `/pipeline-plan` so the planner never
re-explores.

The reason it exists is the persistence rule: **every wave writes to disk before
the next begins.** A finding that lives only in chat cannot be cited by a phase
brief, cannot be checked by a reviewer, and is gone at the end of the session.
Everything else here is scaffolding around that.

## When to use

- A question is open enough that one agent reading one way would miss most of
  the answer — the codebase, the literature, and someone else's implementation
  each hold a different part of it.
- The work downstream is a multi-stage pipeline, so the research has to survive
  into planning as a citable document.

## When not to use

- The question is already answered and you know the shape of the work → go
  straight to `/pipeline-plan`.
- One lane would do → read it and write a topic note. Three explorers for a
  question one grep answers is a fan-out that costs more than it finds.

## What it writes

```
docs/_internal/_project/{slug}-research-<lane>.md   # one per explorer
docs/_internal/_project/{slug}-synthesis.md         # what the planner reads
```

Flat `_project/` topic notes, because the research spans stages and precedes
them. `{slug}-synthesis.md` is the path `/pipeline-plan` already resolves
`--scope-doc` to by default, so the two commands name one file.

Re-running the same slug updates these in place. The memory is its own git
repository, which is what makes that safe.

## When arguments are missing

Both positionals are required and neither has a sensible default. If either is
absent, **ask what should be researched** and propose a slug from the answer —
what is missing is the question, not the syntax.

A question under five words cannot anchor a fan-out; three explorers would each
invent their own reading of it. Say so and ask for a fuller one:

```
"{question}" is too short to anchor a research fan-out. What is the
question in full? For example:
  /explore-and-plan "How should we detect TE-derived transcripts in bulk RNA-seq?" te-detection
```

## Arguments

| Param | Shape | Default | Meaning |
|---|---|---|---|
| `question` | positional | ask | The scientific question. Sets the frame for every explorer and the synthesizer. |
| `slug` | positional | ask | Names both output files, and is forwarded to `/pipeline-plan`. |
| `--n-explorers <int>` | flag | `3` | Parallel explorers, one lane each. |
| `--idempotency-peer <slug>` | flag | none | A sibling pipeline whose stage-ids, checkpoint names and column prefixes the new plan must not clobber. Turns disjointness into a planning gate. |

Flag order does not matter. Announce the resolved configuration before starting:

```
/explore-and-plan "{question}"  slug={slug}
  Notes:            docs/_internal/_project/{slug}-*.md
  Explorers:        {n-explorers} (Opus, parallel)
  Lanes:            codebase · literature · reference-implementation
  Idempotency peer: {peer slug | none}
```

## Wave 1 — explore (Opus, fan-out)

Dispatch `--n-explorers` Opus explorers **in parallel**, one lane each:

| Lane | Focus |
|---|---|
| `codebase` | Repo layout, existing `02_analysis/stages/`, `analysis_config.yaml` (stages, figures, paths), helper libraries, prior plans under `docs/_internal/_project/plans/`, existing artifacts under `03_results/` — the internal ground truth. |
| `literature` | Published methods, benchmarks, canonical approaches. Surface only what is verifiable by search or an in-repo reference; do not invent citations. |
| `reference-impl` | External repos or toolkits solving an analogous problem: patterns, parameter choices, conventions worth importing or avoiding. |

At a count other than 3, distribute lanes evenly and give each explorer a
non-overlapping focus described in its prompt.

**Each explorer writes its own trace** to
`docs/_internal/_project/{slug}-research-<lane>.md`, and the file exists on disk
before the explorer returns. If an explorer cannot write, the orchestrator
writes the file on receiving its output — before any later wave starts.

Trace format:

```markdown
# Research lane: <lane> — {slug}

**Date:** {date} · **Explorer:** Opus, {lane} lane

## Question
{the question as posed}

## Sources consulted
- {file path / URL / DOI for each source read or searched}

## Findings
### <short label>
**Evidence:** {what was observed or read}
**Status:** verified | inferred
**Relevance:** {why this matters for the question}

## Open questions for the synthesizer
- {what this lane could not resolve}
```

**Stop if a trace is absent or empty** after its explorer returns. Name the lane
and its file path; a chat-only lane cannot be synthesized.

## Wave 2 — synthesize (Opus, one agent)

One Opus synthesizer reads every `{slug}-research-*.md` and writes
`docs/_internal/_project/{slug}-synthesis.md`. It does not re-explore: the lane
traces are this wave's ground truth. If a trace is missing, stop rather than
silently synthesizing without it.

Every claim carries `[verified]` — directly observed in a cited source — or
`[inferred]` — a conclusion drawn from evidence, plausible but not confirmed.
**No synthesis-original claims:** each finding traces to at least one lane.

```markdown
# Synthesis: {slug}

**Date:** {date} · **Lanes read:** {list}

## Summary
1–3 sentences: the most important insight for planning.

## Verified findings
- [verified] {finding} — {source lane}

## Inferred findings
- [inferred] {finding} — {source lane} · {basis for inference}

## Open questions
- {unresolved items the planner should know about}

## Peer namespace
<!-- present only when --idempotency-peer is set -->
Peer pipeline: {peer-slug}
Known stage-ids / checkpoint names / column prefixes: {lists}
The new plan must not reuse any of these identifiers.

## Recommended approach
{1–2 paragraphs: the most grounded implementation path the three lanes together
suggest, citing the source lane for each recommendation.}
```

**Stop if the synthesis is not on disk** when the synthesizer returns. The file
is the handoff artifact; synthesis content is never passed inline.

## Wave 3 — hand off

```
/pipeline-plan {slug} --scope-doc docs/_internal/_project/{slug}-synthesis.md
```

The flag is passed explicitly even though it matches `/pipeline-plan`'s default,
so the contract between the two commands is visible in the invocation.

With `--idempotency-peer` set, add to the handoff context — `/pipeline-plan`
takes no such flag:

> The scope-doc carries a `## Peer namespace` section. Verify that every
> stage-id, checkpoint name and column prefix in the new `00_INDEX.md` and
> `NN_<slug>.md` phase files is disjoint from it. A collision with the peer
> pipeline is a planning failure — resolve it before the execute phase.

`/pipeline-plan` reads only the synthesis. The lane traces are the explorers'
working notes; the synthesis is what phase briefs cite.

## Report

```
/explore-and-plan {slug} complete.

  docs/_internal/_project/
    {slug}-research-codebase.md        (Wave 1, codebase lane)
    {slug}-research-literature.md      (Wave 1, literature lane)
    {slug}-research-reference-impl.md  (Wave 1, reference-impl lane)
    {slug}-synthesis.md                (Wave 2, verified vs inferred)

  Verified: {count} · Inferred: {count} · Open: {count}
  Idempotency peer: {peer slug | none}

Handed off to /pipeline-plan — see its output for the phase table.
```

Interrupted at a wave boundary, report which wave failed and the missing file,
and recommend re-running from there.

## Rules

1. **Persist every wave before the next begins.** Lane traces exist before the
   synthesizer starts; the synthesis exists before `/pipeline-plan` is invoked.
   This is the point of the command, not an optimization.
2. **Verified or inferred, tagged on every claim.** The distinction is what
   keeps the planner from treating a guess as a fact.
3. **The planner never re-explores.** It reads the synthesis and nothing else.
   A phase brief cites the synthesis precisely enough that its implementer does
   not go looking either.
4. **Distinct lanes, no overlap.** Explorers do not read each other's traces —
   reconciling them is the synthesizer's job.
5. **Peer disjointness is a gate, not a warning.** With `--idempotency-peer`
   set, a stage-id or column-prefix collision stops the plan.
6. **State the tier at every dispatch.** Explorers and synthesizer are Opus:
   open-ended research and cross-lane reconciliation are judgement.
7. **Partial runs are resumable.** Traces already written are valid. On a
   re-run, detect them and confirm before re-exploring a completed lane.
