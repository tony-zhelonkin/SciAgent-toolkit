# Explore and Plan

Research fan-out → synthesis → pipeline plan. Runs `--n-explorers` parallel Opus explorers across distinct lanes (codebase / literature / reference implementations), persists every lane trace to disk, synthesizes the traces into a single verified-vs-inferred `_SYNTHESIS.md`, and hands off to `/pipeline-plan` so the planner never re-explores. Every wave writes to the repo before the next begins — a chat-only finding is a failure.

## Phase 0: Parse arguments

| Param | Shape | Default | Meaning |
|---|---|---|---|
| `question` | `$ARGUMENTS[0]` (positional, **required**) | — | The scientific question driving the research fan-out. Sets the frame for all explorers and the synthesizer. |
| `slug` | `$ARGUMENTS[1]` (positional, **required**) | — | Research slug. Resolves the research dir to `docs/_internal/research/{today}-{slug}/` (date = `date +%F`). The same slug is forwarded to `/pipeline-plan`. |
| `--n-explorers <int>` | flag | `3` | Number of parallel Opus explorers. Each explorer owns a distinct lane; the default three lanes are **codebase** (repo layout, existing stages, stage-ids, config), **literature** (published methods, known approaches, relevant papers), and **reference-implementation** (external repos or toolkits that solve an analogous problem). |
| `--idempotency-peer <slug>` | flag | none | A sibling pipeline slug whose namespaces (stage-ids, checkpoint names, master-table column prefixes) the new plan must not clobber. When set, the synthesizer annotates the `_SYNTHESIS.md` with peer namespace facts and the `/pipeline-plan` handoff includes an explicit instruction to verify disjointness. |

Flag parsing is order-independent. If either positional argument is missing, reject:
```
Usage: /explore-and-plan "<question>" <slug> [--n-explorers <int>] [--idempotency-peer <slug>]
```

Resolve and announce the run config:
```
/explore-and-plan "{question}"  slug={slug}
  Research dir:    docs/_internal/research/{date}-{slug}/
  Explorers:       {n-explorers} (Opus, parallel)
  Lanes:           codebase · literature · reference-implementation
  Idempotency peer: {peer slug | none}
```

**STOP condition — ambiguous question.** If `question` is empty or fewer than five words, reject:
```
The question "{question}" is too short to anchor a research fan-out.
Provide a full scientific question (≥5 words), e.g.:
  /explore-and-plan "How should we detect TE-derived transcripts in bulk RNA-seq?" te-detection
```
Stop.

## Phase 1: Initialize research directory

Create `docs/_internal/research/{date}-{slug}/` if it does not already exist. Do NOT proceed if the directory cannot be created (permissions, path collision). Announce:
```
Research dir ready: docs/_internal/research/{date}-{slug}/
```

This directory must exist on disk before any explorer writes — explorers MUST NOT buffer output in memory and write at the end; each explorer writes its trace as it completes so a partial run is recoverable.

## Phase 2: Wave 1 — explore (model tier: **Opus**, fan-out)

Dispatch `--n-explorers` Opus explorers **in parallel**. Each explorer owns one distinct lane. With the default of 3 explorers, the lane assignments are:

| Explorer | Lane | Focus |
|---|---|---|
| 01 | `codebase` | Repo layout, existing `02_analysis/stages/`, `analysis_config.yaml` (stages, figures, paths), helper libraries, prior plan dirs under `docs/_internal/plans/`, any existing results artifacts under `03_results/` — the internal ground truth. |
| 02 | `literature` | Published methods, benchmarks, canonical approaches, key papers for the scientific question. Do NOT invent citations; surface only what is verifiable via web search or in-repo references. |
| 03 | `reference-impl` | External repos or toolkits that solve an analogous problem: code patterns, parameter choices, conventions worth importing or avoiding. |

If `--n-explorers` differs from 3, distribute lanes evenly across explorers; assign each explorer a non-overlapping focus; describe its lane in its prompt.

**Persistence rule — each explorer writes its own trace.** Each explorer MUST write its output to:
```
docs/_internal/research/{date}-{slug}/NN_<lane>.md
```
where `NN` is the zero-padded explorer index (`01`, `02`, `03`, …) and `<lane>` is the lane name (`codebase`, `literature`, `reference-impl`). **The file must exist on disk before the explorer returns.** If a read-only-mode explorer cannot write files directly, the orchestrator writes the file immediately upon receiving the explorer's output — before starting any later wave. A lane trace that exists only in chat is a failure of this command's core purpose.

Lane trace format (each `NN_<lane>.md`):
```markdown
# Research Lane: <lane>

**Date:** {date}  ·  **Explorer:** Opus (lane {NN})  ·  **Slug:** {slug}

## Question
{the question as posed}

## Sources consulted
- {file path / URL / paper DOI for each source read or searched}

## Findings
### Finding 1 — <short label>
**Evidence:** {what was observed or read}
**Status:** verified | inferred
**Relevance:** {why this matters for the question}

### Finding 2 — <short label>
[repeat the evidence → status → relevance block]

## Open questions for the synthesizer
- {anything this explorer could not resolve}
```

**STOP condition — explorer failure.** If an explorer fails to produce a trace file (the file is absent or empty after the explorer returns), stop:
```
Explorer {NN} ({lane}) did not write docs/_internal/research/{date}-{slug}/NN_<lane>.md.
A chat-only finding cannot be synthesized. Re-run or investigate the explorer.
```
Stop.

After all `--n-explorers` lane traces are on disk, announce:
```
Wave 1 complete — {n-explorers} lane traces written:
  docs/_internal/research/{date}-{slug}/01_codebase.md
  docs/_internal/research/{date}-{slug}/02_literature.md
  docs/_internal/research/{date}-{slug}/03_reference-impl.md
  [additional lanes if n-explorers > 3]
Proceeding to Wave 2 — synthesis.
```

## Phase 3: Wave 2 — synthesize (model tier: **Opus**, 1 synthesizer)

Dispatch **one Opus synthesizer**. The synthesizer:

1. **Reads** all `NN_<lane>.md` files from `docs/_internal/research/{date}-{slug}/`. It does NOT re-explore the codebase or the literature — the lane traces are the source of truth for this wave. If a trace is absent, stop (do not silently omit it).

2. **Produces** `docs/_internal/research/{date}-{slug}/_SYNTHESIS.md` — the single document the planner will read. The planner trusts the synthesis; it never re-explores downstream.

3. **Tags every claim** as one of:
   - `[verified]` — directly observed in a source the explorer cited (file content, committed artifact, canonical paper)
   - `[inferred]` — a conclusion the explorer drew from evidence; plausible but not directly confirmed

   No synthesis-original claims are permitted: every finding traces to at least one lane trace.

4. **If `--idempotency-peer` is set**, the synthesizer adds a `## Peer namespace` section to `_SYNTHESIS.md` listing the peer pipeline's known stage-ids, checkpoint names, and master-table column prefixes (sourced from the codebase lane's findings). The planner will use this section to verify disjointness.

`_SYNTHESIS.md` structure:
```markdown
# Synthesis: {slug}

**Date:** {date}  ·  **Synthesizer:** Opus  ·  **Lanes read:** {list}

## Summary
1–3 sentences: the single most important insight for planning.

## Verified findings
- [verified] {finding} — {source lane}
- ...

## Inferred findings
- [inferred] {finding} — {source lane} · {basis for inference}
- ...

## Open questions
- {unresolved items the planner should be aware of, sourced from explorer open-questions}

## Peer namespace
<!-- present only when --idempotency-peer is set -->
Peer pipeline: {peer-slug}
Known stage-ids: {list}
Known checkpoint names: {list}
Known master-table column prefixes: {list}
The new plan MUST NOT reuse any of these identifiers.

## Recommended approach
{1–2 paragraphs: what the codebase, literature, and reference-impl lanes together suggest as the most grounded implementation path. Cite the source lane for each recommendation.}
```

**Persistence rule — synthesis must be written before handoff.** `_SYNTHESIS.md` MUST exist on disk before Phase 4 begins. Do not pass synthesis content inline to `/pipeline-plan`; the file is the handoff artifact. Announce:
```
Wave 2 complete — synthesis written:
  docs/_internal/research/{date}-{slug}/_SYNTHESIS.md
  Verified findings: {count}
  Inferred findings: {count}
  Open questions:    {count}
Proceeding to Wave 3 — /pipeline-plan handoff.
```

**STOP condition — synthesis absent.** If `_SYNTHESIS.md` is not on disk after the synthesizer returns, stop:
```
Synthesizer did not write docs/_internal/research/{date}-{slug}/_SYNTHESIS.md.
The planner cannot proceed without a persisted synthesis.
```
Stop.

## Phase 4: Wave 3 — plan handoff

Hand off to `/pipeline-plan` by invoking:
```
/pipeline-plan {slug} --scope-doc docs/_internal/research/{date}-{slug}/_SYNTHESIS.md
```

If `--idempotency-peer` is set, append the instruction to the `/pipeline-plan` invocation context (not as a flag — `/pipeline-plan` does not accept `--idempotency-peer` directly):

> The `_SYNTHESIS.md` scope-doc contains a `## Peer namespace` section. The Opus planner MUST verify that every stage-id, checkpoint name, and master-table column prefix in the new `00_INDEX.md` and `NN_<slug>.md` phase files is disjoint from the peer namespace listed there. A collision with the peer pipeline is a planning failure — resolve it before proceeding to the execute phase.

`/pipeline-plan` reads only the `_SYNTHESIS.md`; it does NOT re-read the lane traces. The lane traces are the explorers' working notes; the synthesis is the authoritative source of truth. This is the 14839 rule: **the planner never re-explores; each phase §3 cites the research note.**

## Phase 5: Output summary

When all three waves complete and `/pipeline-plan` has been handed off:
```
/explore-and-plan {slug} complete.

Research dir: docs/_internal/research/{date}-{slug}/
  01_codebase.md       (Wave 1 — Opus explorer, codebase lane)
  02_literature.md     (Wave 1 — Opus explorer, literature lane)
  03_reference-impl.md (Wave 1 — Opus explorer, reference-impl lane)
  _SYNTHESIS.md        (Wave 2 — Opus synthesizer, verified-vs-inferred)

Synthesis stats:
  Verified findings: {count}
  Inferred findings: {count}
  Open questions:    {count}
  Idempotency peer:  {peer slug | none}

Handoff:
  /pipeline-plan {slug} --scope-doc docs/_internal/research/{date}-{slug}/_SYNTHESIS.md

/pipeline-plan is now running. See its output for plan dir, phase table, and execution progress.
```

If the run was interrupted at a wave boundary (a lane trace or the synthesis is absent), report which wave failed and the missing file path, and recommend re-running from the failed wave.

## Rules

1. **Persist every wave before the next begins.** Wave 1 lane traces exist on disk before the synthesizer starts. `_SYNTHESIS.md` exists on disk before `/pipeline-plan` is invoked. This is not an optimization — it is the point of the command. A chat-only finding is a failure.
2. **Verified vs inferred, explicitly tagged.** Every claim in `_SYNTHESIS.md` carries `[verified]` or `[inferred]`. The distinction protects the planner from treating guesses as facts. No synthesis-original claims.
3. **The planner never re-explores.** `/pipeline-plan` reads only `_SYNTHESIS.md`. The 14839 discipline: phase §3 cites the research note precisely so implementers don't re-explore. The lane traces are the explorers' working notes, not a second scope-doc.
4. **Distinct lanes, no overlap.** Each explorer covers a non-overlapping part of the question space. Explorers do not read each other's traces — that is the synthesizer's job.
5. **Idempotency-peer disjointness is a planning gate, not a warning.** When `--idempotency-peer` is set, the planner must verify namespace disjointness before proceeding to Phase 2 (execute) of `/pipeline-plan`. A stage-id or column-prefix collision with the peer pipeline is a hard stop.
6. **Model tiering is explicit.** Explorers = **Opus** (open-ended research requires judgment). Synthesizer = **Opus** (cross-lane reconciliation and claim-tier assignment require judgment). `/pipeline-plan` inherits its own model tiering (Opus plan, Sonnet implement, Opus review). State the tier at every dispatch.
7. **Recoverable on partial run.** If the command is interrupted mid-wave, the already-written lane traces are valid. Re-running should detect existing traces and skip re-exploration for lanes already completed, or prompt the user to confirm re-exploration.
