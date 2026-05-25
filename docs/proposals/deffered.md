# Deferred design topics

Topics the grilling has explicitly deferred. Not killed — deferred until a concrete consumer or pain point arrives. Each entry records *what was proposed*, *why it was deferred*, *what (if anything) replaces it for now*, and *the revisit trigger*.

Decisions land here when "no current concrete use case" is the conclusion of grilling against the **minimal but no less** preference. Items here can be promoted back into the active proposal queue any time the revisit trigger fires.

---

## ADR-005 — trace recording + workflow-skill-creator infrastructure

**Deferred**: 2026-05-24.

**What was proposed**
A JSONL conversation-trace recorder spanning Claude Code and Pi sessions, plus a `workflow_skill_creator` distiller modelled on Google Antigravity's interview-driven four-phase gate model. The recorder writes session events to disk; the distiller reads those traces (or, in the Antigravity model, interview transcripts) and scaffolds new skills.

**Deferral scope**
The entire ADR — both the recorder (proposed 005b) and the interview-driven distiller (proposed 005a). The earlier recommendation in `ai-research/11-architectural-backbone-and-trace-necessity.md` §3 was to ship 005a and defer 005b. Anton chose to defer both.

**What replaces it (in scope now)**
A single command/agent that ingests **documentation and/or a codebase** as input and emits a skill or family of skills covering the capabilities therein. Matches Anton's existing personal workflow for skill creation. No conversation tracing, no interview gates, no four-phase distillation pipeline — just a generator with a clear input/output contract.

Treat this as a **future ADR placeholder**, not an ADR-005 reformulation. When it is formalised it will be a fresh ADR with its own grilling pass, not a rebrand of the deferred design.

**Why deferred**
- **No second consumer for traces.** ADR-004 is not committed (still grilling deployment shape). ADR-006 benchmarking is being scoped down to MVP-or-none. Building trace infrastructure for hypothetical consumers violates the **minimal but no less** preference.
- **The interview-driven distiller is heavier than the current workflow needs.** Anton's go-to method — point an agent at docs/code and ask for skills — produces working skills without a four-phase gate model. The gates would be ceremony, not value, for solo iteration.
- **Recorder is harness-coupled by construction.** Adding it now requires Claude Code hook instrumentation *and* Pi event handler instrumentation. The engineering cost compounds against the **harness-agnostic preference** without an offsetting concrete benefit.

**Revisit triggers**
- A second consumer of conversation traces materialises. Specifically: (a) the ERA optimizer overlay (ADR-004) ships and needs replay capability, *or* (b) the benchmark harness (ADR-006) needs trace artefacts as input, *or* (c) a multi-user setup needs trace sharing for collaboration.
- The docs/code → skill generator proves insufficient and the pain points specifically point at "I wish I had a recording of that session."
- A regulatory or audit obligation makes trace recording table stakes (unlikely in solo research; possible if sciagent is adopted in a regulated org).

**Cross-references**
- Background research: `ai-research/01-trace-recording-claude-code.md`, `ai-research/02-trace-recording-pi-coding-agent.md`, `ai-research/03-cross-harness-trace-standards.md`, `ai-research/04-anthropic-skills-and-workflow-skill-creator.md`, `ai-research/11-architectural-backbone-and-trace-necessity.md` §3.
- Spec section deferred: ADR-005 in `sciagent-extension-design-spec.md`.

---

## ADR-004 — ERA-style optimizer overlay (and ADR-006 — benchmark harness)

**Deferred**: 2026-05-24.

**What was proposed**
ADR-004: a `scorable` overlay carrying mutator/scorer/selector sub-agents implementing Aygün et al.'s Flat UCB Tree Search to drive automated optimisation of skill/role configurations. ADR-006: a two-tier benchmark harness (tier-1 hand-curated tasks, tier-2 external bio-benchmark suite) providing the fitness signal the optimizer consumes.

**Deferral scope**
Both ADRs deferred together. They form a single dependency chain — ADR-006 exists only to feed ADR-004; deferring 004 without 006 leaves a benchmark with no consumer, and deferring 006 without 004 leaves an optimizer with no signal. The two only cohere as a pair, so they leave or stay as a pair.

**What replaces it (in scope now)**
Nothing. Sciagent ships as a fully opinionated personal context manager. Role and skill quality is judged by the felt friction of daily work on real projects (DC_hum_verse, pathway-explorer, coresh-shared-cache). The fitness signal is "does Anton ship faster," delivered every working day at zero marginal cost.

**Why deferred**
- **No second consumer for an optimizer.** Solo researcher, single primary harness. The "shareable value" of objective scoring is hypothetical; the cost (benchmark authoring, grader maintenance, CI wiring) is real today.
- **Anton's own stated rule applied honestly.** From auto-memory: "each addition must justify against a *current concrete* use case." ADR-004 fails that test first; ADR-006 fails it as a transitive consequence. Honoring the rule means deferring both.
- **No clear benchmark tasks exist yet.** Anton: *"I'm not really sure how to benchmark which tasks to benchmark on. So for now we keep this opinionated."* When the right tasks aren't obvious, authoring them prematurely encodes guesses as gold and then Goodharting takes over.
- **N=1 self-evaluation is the wrong measurement protocol.** Anton would be task author, grader, and consumer. Three layers of the same bias compound; no information beyond felt friction is extracted.
- **Zero lock-in cost from waiting.** ERA's Flat UCB Tree Search is 154 LOC under Apache-2.0 (`docs/.ref/era/implementation/futs.py`). Borrowing it later costs the same as borrowing it now; the algorithm does not get harder to port.
- **The MVP-minimal benchmark middle path is the worst of three options.** It pays setup cost (5–10 tasks authored + graded by Anton) with no optimizer to consume the signal and no second voice to disagree. Pure overhead.

**Revisit triggers** (any one of):
- A second human user adopts sciagent and requests configurable optimisation.
- Pi-harness and Claude-Code-harness produce divergent role rankings Anton cannot reconcile by inspection (concrete A/B disagreement, not abstract concern).
- Sciagent grows ≥3 role overlays whose pairwise differences are too subtle to eyeball-rank.
- A funded collaboration or paper submission requires objective metric reporting on toolkit choices.

**Reactivation surface (when triggered)**
- Copy `implementation/futs.py` from `github.com/google-research/era` (Apache-2.0, ~150 LOC).
- Preserve the Apache-2.0 header; cite arXiv:2509.06503 in any derived doc.
- Ignore the rest of the ERA repo — the published "applications" are PDFs only, no source. The "evaluation benchmark" mentioned in ERA's README ships as two example Jupyter notebooks, not a runnable suite.
- Wire `generate_fn` / `execute_fn` (the two callables the algorithm needs) to whatever benchmark tasks exist at reactivation time.
- Reopen ADR-006 first; the optimizer needs a fitness signal before it makes sense.

**Cross-references**
- ERA inspection report: Opus agent 2026-05-24 (license: Apache-2.0 maximally permissive; algorithm is 154 LOC; no orchestration shell to borrow).
- Background research: `ai-research/06-era-tree-search-paper.md`, `ai-research/08-benchmark-harnesses-bio.md`, `ai-research/11-architectural-backbone-and-trace-necessity.md` §4, `ai-research/11-architectural-backbone-and-trace-necessity.md` §4a.
- Spec sections deferred: ADR-004 and ADR-006 in `sciagent-extension-design-spec.md`.
- Worked example using ADR-004 (kickoff §4 sub-agent bisection): preserved as illustrative of the [Arch]/[Eng]/[Taste] taxonomy; the decisions themselves are moot under this deferral.

---

## ADR-008 — `lab-loop` overlay template

**Dropped**: 2026-05-24.

**What was proposed**
`lab-loop` as one of five "swap-target stub" overlays the ADR-004 `scorable` optimizer would sample over — a stand-in workflow shape (data-exploration loops) the optimizer could mutate skills/agents in and around.

**Status**: dropped from the active spec. Not killed forever — reactivation is conditional on ADR-004 reactivating, since `lab-loop`'s justification was *as a swap target for the optimizer*. With no optimizer, no swap target.

**What replaces it (in scope now)**: nothing. If a generic "lab-loop" role turns out to be useful on its own merits (separate from the optimizer story), it can be authored later as a standalone role with its own justification — but no such authoring is in scope from this grilling round.

**Why dropped**
- Original motivation was downstream of ADR-004; that ADR is deferred. With the consumer gone, the producer is unjustified.
- Reauthoring `lab-loop` as a "standalone useful role" would require fresh justification (concrete workflow it supports). No such case was named during grilling.

**Revisit trigger**: ADR-004 reactivates and the optimizer needs swap-target stubs. At that point `lab-loop` (and the four other stubs originally proposed) get re-grilled together as part of ADR-004's reactivation surface.

---

## ADR-009 — `uv` pinning strategy for skill Python deps

**Deferred**: 2026-05-24.

**What was proposed**
A toolkit-level decision on whether skills should declare Python dependencies via `uv` (mandatory, optional, or container-fast-path variants). The original spec considered making `uv` a first-class part of how skills bring their own deps into a project.

**Status**: deferred. Not killed — venv management is a real concern in non-containerized workflows.

**What replaces it (in scope now)**: nothing. Skills currently in the toolkit declare deps in whatever ad-hoc form their author chose; sciagent does not enforce or normalize.

**Why deferred** (Anton verbatim)
*"UV seems to be optional at the moment. I personally work mostly within Docker containers. I have everything pre-installed there and I'm not that concerned with managing virtual environments."*

Workflow assumption baked into the deferral: Anton's current setup is Docker-container-with-pre-installed-envs. Skills run in that container; deps are managed by the container build, not by per-skill venv management. `uv` would add ceremony without solving a current pain.

**Revisit triggers**
- Sciagent adopted by a workflow that isn't containerized (a collaborator with bare-metal Python, a CI environment without the project's container).
- A skill needs a Python dep that isn't in the container and adding it to the container is more friction than per-skill pinning would be.
- Reproducibility-for-publication requirements force per-skill pin manifests.

---

*Format note: future deferrals append new top-level sections (`## ADR-XXX — <topic>`) to this file. Do not re-order; chronology is itself useful context.*
