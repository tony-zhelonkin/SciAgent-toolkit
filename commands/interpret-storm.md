# Interpret Storm

Multi-wave interpretation of a result set → information-graphic / interactive-viewer design. Runs four waves — **web research** (fan-out), **interpret** (`--n-interpreters` Opus interpreters across distinct angles + one Opus synthesizer), **design** (`--n-designers` Opus designers proposing information graphics or interactive viewers), and **build** (hand off each accepted design to `/add-figure-variant`). **Every wave writes to the repo before the next begins.** Ephemerality is the exact failure mode this command exists to prevent: a wave that lives only in chat is a failure, not a step. See `skills/reasoning-trace/SKILL.md` — capture the answer, delete the shell; a chat summary is ephemeral.

## Phase 0: Parse arguments

| Param | Shape | Default | Meaning |
|---|---|---|---|
| `dataset-path` | `$ARGUMENTS[0]` (positional, **required**) | — | The results to interpret — a checkpoint, master table, or `03_results/<stage-id>/` artifact directory. This is the internal ground truth the interpreters read alongside the Wave-1 web notes. |
| `slug` | `$ARGUMENTS[1]` (positional, **required**) | — | Interpretation slug. Resolves the research dir to `docs/_internal/research/{today}-{slug}/` (date = `date +%F`) and the design doc to `docs/_internal/reasoning/{today}-{slug}-figures.md`. The slug is forwarded (prefixed) to each `/add-figure-variant` hand-off in Wave 4. |
| `--n-interpreters <int>` | flag | `3` | Number of Opus interpreters in Wave 2. Each takes a distinct interpretive angle; the default three angles are **mechanism**, **pathway**, and **clinical**. |
| `--n-designers <int>` | flag | `2` | Number of Opus designers in Wave 3. Each proposes one information graphic or interactive viewer that answers a synthesis gap. |

Flag parsing is order-independent. If either positional argument is missing, reject:

```
Usage: /interpret-storm <dataset-path> <slug> [--n-interpreters <int>] [--n-designers <int>]
```

Resolve and announce the run config:

```
/interpret-storm {dataset-path}  slug={slug}
  Dataset:       {dataset-path}
  Research dir:  docs/_internal/research/{date}-{slug}/
  Design doc:    docs/_internal/reasoning/{date}-{slug}-figures.md
  Interpreters:  {n-interpreters} (Opus, parallel) — angles: mechanism · pathway · clinical
  Designers:     {n-designers} (Opus, parallel)
  Web tools:     {available | unavailable — degrades to "cite what you can reach"}
```

**STOP condition — dataset absent.** If `dataset-path` does not exist on disk (file or directory), do NOT fabricate results to interpret. Stop:

```
Dataset {dataset-path} does not exist. There is nothing to interpret.
Pass an existing checkpoint, master table, or 03_results/<stage-id>/ artifact dir.
```

Stop.

## Phase 1: Initialize research directory

Create `docs/_internal/research/{date}-{slug}/` if it does not already exist. Do NOT proceed if the directory cannot be created (permissions, path collision). Announce:

```
Research dir ready: docs/_internal/research/{date}-{slug}/
```

This directory must exist on disk before any wave writes. No wave buffers its output in memory and flushes at the end — each strand, each note, each synthesis lands on disk as it completes, so a partial run is recoverable from exactly the wave that was interrupted.

## Phase 2: Wave 1 — web research (model tier: **Opus**, fan-out)

Gather external context that the interpreters will read alongside the dataset: published literature, pathway/interaction databases, prior art, canonical methods. Dispatch a fan-out of Opus research strands (one per topic — e.g. the genes/pathways/cell-types that dominate the dataset, the disease or perturbation context, the methods used to generate the results).

**Web tools.** Use the harness's web tools (`WebSearch` / `WebFetch`) where available, plus any in-harness literature tools (e.g. PubMed). **When web is unavailable, the command degrades gracefully: cite what you can reach** — in-repo references, the dataset's own metadata, prior research notes under `docs/_internal/research/` — and say so explicitly in each note. The command ALWAYS persists what it found, even when "what it found" is "web unavailable; grounded only in {in-repo sources}." A degraded Wave 1 is still a persisted Wave 1.

**Persistence rule — each strand writes its own note.** Each research strand MUST write to:

```
docs/_internal/research/{date}-{slug}/web-<topic>.md
```

where `<topic>` is a short kebab-case label (`web-tnf-signaling.md`, `web-cd8-exhaustion.md`, `web-prior-art.md`). **The file must exist on disk before Wave 2 begins.** If a read-only strand cannot write directly, the orchestrator writes the note immediately on receiving the strand's output — before any interpreter is dispatched. A web finding that lives only in chat is a failure of this command's core purpose.

Each `web-<topic>.md` follows the `reasoning-trace` note format:

```markdown
# Web research: <topic>

**Date:** {date}  ·  **Wave:** 1 (web research, Opus)  ·  **Slug:** {slug}

## Scope
One sentence: what external question this strand covers and which part of the dataset it serves.

## Sources
- {URL / DOI / PMID / in-repo path for each source reached}
- Web status: {available | unavailable — grounded only in {in-repo sources}}

## Findings
### Finding 1 — <short label>
**Evidence:** {what the source said}
**Citation:** {URL / DOI / PMID — exact, never invented}
**Relevance:** {why this matters for interpreting the dataset}

### Finding 2 — <short label>
[repeat the evidence → citation → relevance block]

## Open questions for the interpreters
- {anything this strand could not resolve from external context}
```

**Never invent a citation.** If a claim cannot be tied to a reachable source, mark it `[uncited — interpreter to verify against dataset]` rather than attaching a fabricated DOI.

**STOP condition — web wave not persisted.** If no `web-*.md` note exists on disk after Wave 1 returns, stop:

```
Wave 1 wrote no docs/_internal/research/{date}-{slug}/web-*.md notes.
Even a web-unavailable run must persist what it could reach. A chat-only web wave
cannot be read by the interpreters. Re-run or write the notes manually.
```

Stop. After the web notes are on disk, announce:

```
Wave 1 complete — web research persisted:
  docs/_internal/research/{date}-{slug}/web-<topic>.md  (×N)
  Web status: {available | unavailable}
Proceeding to Wave 2 — interpret.
```

## Phase 3: Wave 2 — interpret (model tier: **Opus**, `--n-interpreters` + 1 synthesizer)

Dispatch `--n-interpreters` Opus interpreters **in parallel**. Each interpreter takes a distinct angle. With the default of 3 interpreters, the angle assignments are:

| Interpreter | Angle | Focus |
|---|---|---|
| mechanism | molecular/cellular mechanism | What biological mechanism do these results imply? Which genes/regulators/cell states drive the observed effect, and by what proposed causal chain? |
| pathway | pathway / network | Which pathways, gene sets, or interaction networks are coherent across the results? Where do the enrichments converge or contradict? |
| clinical | clinical / translational | What is the disease, biomarker, or therapeutic relevance? What would a clinician or translational scientist take from this, and at what confidence? |

If `--n-interpreters` differs from 3, distribute angles so each interpreter owns a distinct, non-overlapping interpretive frame; describe its angle in its prompt.

Each interpreter reads (1) the dataset at `dataset-path` and (2) all Wave-1 `web-*.md` notes — it does NOT re-do the web research; the persisted notes are the external source of truth for this wave. Each interpreter writes:

```
docs/_internal/research/{date}-{slug}/reason-<angle>.md
```

(`reason-mechanism.md`, `reason-pathway.md`, `reason-clinical.md`). **The file must exist on disk before the synthesizer is dispatched.** Format:

```markdown
# Interpretation: <angle>

**Date:** {date}  ·  **Wave:** 2 (interpret, Opus)  ·  **Angle:** <angle>  ·  **Slug:** {slug}

## Scope
One sentence: what this angle interprets about the dataset.

## Sources
- Dataset: {dataset-path} — {what was read: columns, objects, tables}
- Web notes read: {list of web-<topic>.md}

## Interpretation
### Claim 1 — <short label>
**Evidence:** {dataset observation + supporting web finding, with citation}
**Status:** verified (directly in the dataset/cited source) | inferred (a conclusion drawn)
**Reading:** {what this means under this angle}

### Claim 2 — <short label>
[repeat the evidence → status → reading block]

## GAPS / open questions
- {what the data cannot yet answer; what would resolve it; which would most change the conclusions}
```

**STOP condition — interpreter trace absent.** If any `reason-<angle>.md` is absent or empty after its interpreter returns, stop:

```
Interpreter <angle> did not write docs/_internal/research/{date}-{slug}/reason-<angle>.md.
A chat-only interpretation cannot be synthesized. Re-run or investigate the interpreter.
```

Stop.

After all `--n-interpreters` traces are on disk, dispatch **one Opus synthesizer**. The synthesizer:

1. **Reads** all `reason-<angle>.md` traces (and may consult the `web-*.md` notes for citations). It does NOT re-interpret the dataset from scratch — the angle traces are the source of truth for this wave.
2. **Produces** `docs/_internal/research/{date}-{slug}/_SYNTHESIS.md` — the single document Wave 3 reads.
3. **Tags every claim** as one of:
   - `[verified]` — directly observed in the dataset or a cited source the interpreter named
   - `[inferred]` — a conclusion an interpreter drew; plausible but not directly confirmed

   No synthesis-original claims: every finding traces to at least one angle trace.
4. **Lists the highest-value gaps** — ranked. These are the gaps the Wave-3 designers will target: each accepted figure must answer one of them.

`_SYNTHESIS.md` structure:

```markdown
# Synthesis: {slug}

**Date:** {date}  ·  **Wave:** 2 (synthesize, Opus)  ·  **Angles read:** {list}

## Summary
1–3 sentences: the single most important interpretation across all angles.

## Verified claims
- [verified] {claim} — {angle} · {dataset/source}
- ...

## Inferred claims
- [inferred] {claim} — {angle} · {basis for inference}
- ...

## Highest-value gaps (ranked — Wave-3 design targets)
1. **<gap label>** — {what is unknown; why resolving it matters most; what figure or viewer would surface it}
2. **<gap label>** — ...
```

**Persistence rule — synthesis before design.** `_SYNTHESIS.md` MUST exist on disk before Wave 3 begins. Do not pass synthesis content inline to the designers; the file is the hand-off artifact. Announce:

```
Wave 2 complete — interpretation persisted:
  reason-<angle>.md  (×{n-interpreters})
  _SYNTHESIS.md      (verified: {count} · inferred: {count} · gaps: {count})
Proceeding to Wave 3 — design.
```

**STOP condition — synthesis absent.** If `_SYNTHESIS.md` is not on disk after the synthesizer returns, stop:

```
Synthesizer did not write docs/_internal/research/{date}-{slug}/_SYNTHESIS.md.
The designers cannot proceed without a persisted synthesis of the gaps.
```

Stop.

## Phase 4: Wave 3 — design (model tier: **Opus**, `--n-designers`)

Dispatch `--n-designers` Opus designers **in parallel**. Each designer reads `_SYNTHESIS.md` (specifically the ranked gaps) and proposes **one** information graphic or interactive viewer that answers a distinct synthesis gap. Two designers, two distinct gaps — no two proposals target the same gap.

Each designer consults `skills/figure-style/SKILL.md` (the single-tier, dual-format PDF + PNG contract, the legibility floors, the `save_overview()` figure+table+caption discipline) so the proposal is buildable by `/add-figure-variant` without redesign.

Every proposal MUST declare, concretely:

- **Namespace token** — a grep-isolable token (derived from `{slug}` + a short suffix, e.g. `{slug}_mechflow`) that will prefix every new identifier, filename, and artifact stem. This is the token Wave 4 passes to `/add-figure-variant --namespace-token`.
- **Data contract** — exactly what the figure needs: which columns / objects / tables from `dataset-path`, the checkpoint the compute step would write, and the viz inputs. (This is the data contract `/add-figure-variant`'s mini-plan will inherit.)
- **Figure-style plan** — the figure plan: panel layout, what `<stem>.pdf` (vector) and `<stem>.png` (raster) — same geometry, one plot object — each carry, the sub-layout (`_overview/` or `by_contrast/<c>/`), and (for an interactive viewer) the static fallback panel that still satisfies the figure-style contract.
- **Claim tier (L0–L7)** — the epistemic level the figure supports, on the project's interpretation ladder: **L0** raw data · **L1** QC metric · **L2** normalized counts · **L3** statistical test result · **L4** pathway/gene-set enrichment · **L5** comparative claim · **L6** mechanistic inference · **L7** proposed mechanism. A figure that targets an `[inferred]` synthesis gap MUST NOT claim a tier above L6 — it is a proposed reading, not a verified fact.
- **Which gap it answers** — the exact ranked gap from `_SYNTHESIS.md` this design resolves.

**Persistence rule — designs to disk before build.** All proposals are persisted to a single design doc:

```
docs/_internal/reasoning/{date}-{slug}-figures.md
```

**This file must exist on disk before Wave 4 begins.** Format:

```markdown
# Figure designs: {slug}

**Date:** {date}  ·  **Wave:** 3 (design, Opus)  ·  **Synthesis:** docs/_internal/research/{date}-{slug}/_SYNTHESIS.md

## Design 1 — <token suffix>
- **Answers gap:** {ranked gap N from _SYNTHESIS.md}
- **Kind:** information graphic | interactive viewer
- **Namespace token:** `{slug}_<suffix>`
- **Data contract:** compute inputs {columns/objects} → checkpoint {path/schema} → viz inputs {checkpoint}
- **Variant plan (print+screen):** {print panel/geometry} · {screen panel/geometry} · sub-layout {_overview | by_contrast/<c>}
- **Claim tier:** L<n> — {one line: why this tier; never above L6 for an inferred gap}
- **Status of the underlying claim:** [verified] | [inferred] (from _SYNTHESIS.md)

## Design 2 — <token suffix>
[repeat the block — must answer a different gap]
```

**STOP condition — design doc absent.** If `docs/_internal/reasoning/{date}-{slug}-figures.md` is not on disk after the designers return, stop:

```
Wave 3 wrote no docs/_internal/reasoning/{date}-{slug}-figures.md.
The /add-figure-variant hand-offs cannot proceed without persisted, namespaced designs.
```

Stop. After the design doc is on disk, announce:

```
Wave 3 complete — designs persisted:
  docs/_internal/reasoning/{date}-{slug}-figures.md  ({n-designers} designs)
Proceeding to Wave 4 — build (/add-figure-variant hand-offs).
```

## Phase 5: Wave 4 — build (hand off to `/add-figure-variant`)

For **each accepted design** in `docs/_internal/reasoning/{date}-{slug}-figures.md`, hand off to `/add-figure-variant`. Each hand-off runs that command's full compute → viz → Opus-review pipeline; `/interpret-storm` does NOT re-implement figures itself.

```
/add-figure-variant {slug}-<token-suffix> <stage-id> --namespace-token <token>
```

where `<token-suffix>` and `<token>` come from the design's namespace declaration and `<stage-id>` is the results stage the figure belongs to (the stage that owns `dataset-path`). The design doc IS the rationale `/add-figure-variant` reads — its data contract and variant plan seed the mini-plan, so the designs are not re-derived.

Hand off the designs in order. If a design's `<stage-id>` is not yet declared, `/add-figure-variant`'s own Phase-0 STOP fires — declare the stage first, then re-hand-off that design. The other designs are unaffected.

Each `/add-figure-variant` run carries its own model tiering (Opus mini-plan → Sonnet compute → Sonnet viz → Opus review → `captions` cleanup). Do not collapse it — it is the build gate that opens the produced PDFs and verifies namespace isolation.

## Phase 6: Output summary

When all four waves complete and every design has been handed to `/add-figure-variant`:

```
/interpret-storm {slug} complete.

Research dir: docs/_internal/research/{date}-{slug}/
  web-<topic>.md     (×N)              — Wave 1 (Opus web research, citations)
  reason-<angle>.md  (×{n-interpreters}) — Wave 2 (Opus interpreters: mechanism/pathway/clinical)
  _SYNTHESIS.md                          — Wave 2 (Opus synthesizer, verified-vs-inferred + ranked gaps)

Design doc:
  docs/_internal/reasoning/{date}-{slug}-figures.md  — Wave 3 ({n-designers} Opus designs: token + data contract + variant plan + claim tier)

Synthesis stats:
  Verified claims:  {count}
  Inferred claims:  {count}
  Ranked gaps:      {count}
  Web status:       {available | unavailable}

Build hand-offs (Wave 4):
  /add-figure-variant {slug}-<suffix-1> <stage-id> --namespace-token <token-1>
  /add-figure-variant {slug}-<suffix-2> <stage-id> --namespace-token <token-2>
  [one per accepted design]

/add-figure-variant is now running per design. See its output for scripts, artifacts, and the review verdict.
```

If the run was interrupted at a wave boundary (a `web-*.md`, a `reason-<angle>.md`, the `_SYNTHESIS.md`, or the design doc is absent), report which wave failed and the missing file path, and recommend re-running from the failed wave. The already-persisted waves are valid — that is the recoverability the persist-every-wave discipline buys.

## Rules

1. **Persist every wave before the next begins — this is the POINT.** Wave 1 `web-*.md` notes exist on disk before any interpreter is dispatched. Wave 2 `reason-<angle>.md` traces exist before the synthesizer; `_SYNTHESIS.md` exists before any designer. Wave 3 `{date}-{slug}-figures.md` exists before any `/add-figure-variant` hand-off. This is not an optimization — this command exists *because* agents otherwise lose the reasoning in chat. A chat-only wave is a failure, not a step. (`reasoning-trace`: capture the answer, delete the shell.)
2. **Web research degrades, never disappears.** Use `WebSearch`/`WebFetch` (and PubMed-style tools) where available; when web is unavailable, cite what you can reach (in-repo references, dataset metadata, prior notes) and say so. ALWAYS persist what was found. Never invent a citation — mark unsupported claims `[uncited]` for the interpreters to verify against the dataset.
3. **Verified vs inferred, explicitly tagged.** Every claim in `_SYNTHESIS.md` carries `[verified]` or `[inferred]`; the distinction flows into the design doc (claim tier) and protects downstream figures from presenting a guess as a fact. No synthesis-original claims — every claim traces to an angle trace.
4. **Distinct angles, distinct gaps.** Each interpreter owns a non-overlapping interpretive frame; each designer answers a distinct ranked gap. Interpreters do not read each other's traces — that is the synthesizer's job.
5. **Every design is buildable by `/add-figure-variant` without redesign.** A proposal MUST declare a namespace token, a data contract, a figure-style plan (per the `figure-style` contract), and a claim tier (L0–L7). A design missing any of these cannot be handed off.
6. **The build gate is `/add-figure-variant` — do not collapse it.** Wave 4 hands each design to `/add-figure-variant`, which opens the produced PDFs, runs `figure-audit`, verifies namespace isolation, and runs the mandatory `captions` pass. `/interpret-storm` never writes figure scripts itself.
7. **Model tiering is explicit.** Web research = **Opus**. Interpreters = **Opus** (interpretation requires judgment). Synthesizer = **Opus** (cross-angle reconciliation + claim-tier assignment). Designers = **Opus** (design judgment + data-contract reasoning). Each `/add-figure-variant` hand-off carries its own Opus/Sonnet tiering. State the tier at every dispatch.
8. **Recoverable on partial run.** If interrupted mid-wave, the already-written notes are valid. Re-running detects existing notes and skips re-doing completed waves, or prompts to confirm re-running them.
