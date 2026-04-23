---
parent: ./README.md
view: decision
---

# Decisions — architect role

The design choices behind the harness, with rationale and trade-offs.

---

## ADR-001: Single cartography artifact; reviewers never re-explore

### Context
2026 agentic-coding research shows 60–80% of agent tokens go to orientation (figuring out where code lives, what imports what). A naive multi-agent pipeline where each agent re-explores the codebase compounds this tax: N reviewers × full-codebase context.

### Decision
One agent (`mapper`) produces `docs/{feature}/map.md` once per feature. All downstream agents (`bioinf`, `wetlab`, `graphic`, `stat`, `divergent`, `ml`, `synth`, `architect`, plus the main agent during `/design`/`/plan`/`/implement`) read this artifact instead of re-grepping the code.

### Trade-offs
- **Gains:** 1 cartography pass amortized across 1–6 reviewers + design + plan + implement. Agents have smaller, focused contexts. Map can be hand-edited to correct mapper errors before reviewers see it.
- **Costs:** If map.md is stale, all downstream work inherits the staleness. If map.md is incomplete, reviewers can miss things. Mitigation: reviewers add "Open questions for mapper" sections; re-running `/map` is cheap.

### Alternatives considered
- Separate `/scope` and `/research` commands — rejected; scope is a subset of map's blast-radius section, separate commands would double-pay the orientation tax
- Let every agent grep freely — rejected; this is the 60–80% waste pattern documented in the literature
- Share a vector index / semantic map instead of a markdown file — rejected for now; markdown is human-reviewable at a gate, which a binary index is not

### Sources
- RDEL #138 "Where do all the tokens go in agentic software engineering" (2026)
- Nesler, "Your AI Coding Agent Wastes 80% of Its Tokens Just Finding Things" (Medium, Feb 2026)
- Anthropic, "How and when to use subagents in Claude Code"

---

## ADR-002: Composable reviewer panel, not all-or-nothing

### Context
Not every feature needs all lenses. A small typography change doesn't need a statistician. A normalisation audit doesn't need an information graphic designer. An all-or-nothing review would either be overkill for small features or underpowered for cross-cutting ones.

### Decision
`/review` accepts three invocation shapes:
- `--as all` — full 6-reviewer panel
- `--as all --but <csv>` — exclusion syntax
- `--as <csv>` — inclusion (subset)

Synthesis is conditional: `/synthesize` refuses to run with <2 reviewers (reading a single review is cheaper than reading a synthesis of it).

### Trade-offs
- **Gains:** Small features pay only for the lenses they need. The all-with-exclusions syntax is ergonomic for common cases ("everything except divergent" / "everything except wetlab when the change is purely UI")
- **Costs:** User has to choose. Interactive fallback mitigates this.

### Rationale
Naming the choice forces the user to think about what lenses the feature actually needs — this is intentional friction. An auto-chosen set would hide the decision.

---

## ADR-003: Short agent names

### Context
Long agent names (`computational-biologist`, `machine-learning-manifold-scientist`) clutter the `/review --as` command and add typing friction. Claude Code agent names appear throughout the UI and logs.

### Decision
All agents use short names: `bioinf`, `wetlab`, `graphic`, `stat`, `divergent`, `ml`, `mapper`, `architect`, `synth`, `meta-architect`.

### Trade-offs
- **Gains:** `/review umap --as all --but stat,divergent` is readable. Agents fit in tight UIs. Names compose naturally in prompts.
- **Costs:** Less self-descriptive to a newcomer. Mitigation: the agent file's `description` field carries the long form and usage examples.

### Naming conventions
- Reviewers: descriptor of lens (`bioinf`, `wetlab`, `graphic`, `stat`, `ml`, `divergent`)
- Pipeline agents: verb/role (`mapper`, `architect`, `synth`, `meta-architect`)

---

## ADR-004: Human gates at every stage boundary

### Context
Autonomous multi-agent pipelines are tempting but risky: errors compound silently, architectural drift happens without human awareness, and debugging is painful when three stages have elapsed since the bad decision.

### Decision
Every stage ends with a hand-off to the user. The assistant waits. There is no "run the whole pipeline" macro.

### Trade-offs
- **Gains:** User intervenes early; small corrections at map stage are cheap, corrections at implement stage are expensive. User sees each artifact before the next stage consumes it.
- **Costs:** Slower wall-clock time for a perfect run. User must be present.

### Rationale
The purpose of the harness is to slow the assistant down to the user's thinking pace at architectural moments, then speed it back up during implementation. Gates are the feature.

---

## ADR-005: Markdown artifacts as the human interface

### Context
The design review process could be interactive chat, structured JSON, a diff view, or markdown files.

### Decision
Every stage produces a human-readable markdown file under `docs/{feature}/`.

### Trade-offs
- **Gains:** Version-controllable. Diffable. Hand-editable. Readable outside Claude Code (any markdown viewer / GitHub). Durable across sessions. Can be read by humans, by future agents, by LLMs of any vendor.
- **Costs:** Less structured than JSON. Mapper/architect agents have to follow formatting conventions (not type-checked).

### Rationale
The markdown artifact *is* the architecture. If the architect harness disappeared tomorrow, `docs/{feature}/` would still be useful engineering documentation.

---

## ADR-006: Commands as plain markdown templates

### Context
Claude Code supports slash commands at `.claude/commands/<name>.md`. The file content is injected as the user prompt when the slash command fires.

### Decision
Each command is a plain markdown file (no YAML frontmatter needed) that reads as a detailed set of phased instructions. Claude Code injects these verbatim; the main agent then follows them.

### Trade-offs
- **Gains:** No custom DSL. Command authors write in natural English. Commands are self-documenting (reading `commands/review.md` tells you exactly what `/review` does).
- **Costs:** Verbose. No compile-time validation of references (e.g., if a command names an agent that doesn't exist, it fails at runtime).

### Anti-pattern avoided
Earlier prototypes put the command logic inside the main agent's prompt. This made commands opaque — you couldn't read what a command did without asking Claude. Putting dispatch logic in a per-command markdown file makes the behavior inspectable.

---

## ADR-007: Parallel reviewer dispatch is mandatory

### Context
Six reviewers dispatched sequentially is 6× the wall-clock time. Sequential dispatch also doubles the communication cost (each reviewer's summary gets re-read by the main agent before the next dispatch).

### Decision
`/review --as all` dispatches reviewers via a SINGLE assistant message containing multiple Agent tool calls. This is the Claude Code parallel-dispatch pattern.

### Trade-offs
- **Gains:** Near-constant wall-clock time regardless of N reviewers. Each reviewer sees the same snapshot of upstream artifacts.
- **Costs:** No reviewer can build on another reviewer's output within a single run. Mitigation: `/synthesize` is the integration point, and re-running `/review --as <reviewer>` after synthesis lets a specific reviewer react to it.

---

## ADR-008: Role YAML as the only manifest

### Context
A role could be a Python class, a JSON manifest, a set of environment variables, or a declarative YAML file.

### Decision
`roles/<name>.yaml` with three lists: `agents`, `skills`, `commands`. Activation script symlinks by name.

### Trade-offs
- **Gains:** Declarative, inspectable, git-diffable. Easy to write by hand. No language runtime needed. Consistent with existing roles in the toolkit (`base.yaml`, `pathway-signature-agent.yaml`).
- **Costs:** No schema validation beyond what the activation script does (warn on missing files). No cross-role inheritance.

### Backward compatibility
The `commands:` key was added for this role; existing roles that omit it work unchanged (the script silently no-ops if the key is absent).

---

## ADR-009: `divergent` and `ml` use opus; others use sonnet

### Context
Opus is slower and more expensive than sonnet, but better at math-heavy reasoning and multi-hop logical deduction.

### Decision
- `divergent` — opus (needs to reason about hidden assumptions, multi-step failure modes)
- `ml` — opus (manifold math, distance-metric theory)
- `architect` — opus (design consistency checking across multiple documents)
- `synth` — opus (integrates multiple reviewer positions without flattening)
- `bioinf`, `wetlab`, `graphic`, `stat`, `mapper` — sonnet

### Trade-offs
- **Gains:** Spending compute where it matters. Sonnet is plenty for fact-gathering and lens-specific critique; opus pays off on math and contrarian reasoning.
- **Costs:** Slower runs for features that invoke the opus agents.

### Revisit trigger
If Haiku 4.5+ approaches Sonnet quality, consider downgrading `mapper` (pure fact-gathering) to Haiku.

---

## ADR-010: Pinned frontmatter vocabulary for design docs

### Context
First two `/design` runs (umap, normalisation) showed frontmatter drift: `status:` present in one feature's frontmatter but missing in the other, `view: behavioural` (UK) vs `view: behavioral` (US), and ad-hoc values like `view: data` vs `view: pipeline` for the same document role. A frontmatter field that drifts is worse than no frontmatter at all — consumers can't rely on it, and the indexing benefit evaporates.

### Decision
`commands/design.md` Phase 2.5 pins the frontmatter exactly:
- `feature: {slug}`
- `view:` from a fixed vocabulary — exactly `structural | behavioral | decision | data-flow | api-contract` (US spelling; README.md omits `view:` because it's the index, not a view)
- `status:` lives ONLY in frontmatter — never duplicated as `**Status**:` in body
- `date:` ISO form `YYYY-MM-DD`

### Trade-offs
- **Gains:** Frontmatter becomes reliable enough for scripted consumers (grep-by-view, status dashboards). No UK/US drift. Status single-source-of-truth prevents frontmatter and body from disagreeing.
- **Costs:** Slightly more rigid. If a novel view emerges, the vocabulary must be extended in the command spec before it can be used.

### Revisit trigger
If two features drift on a new view dimension (e.g., need a `migration` or `testing` view), extend the vocabulary in `commands/design.md` rather than letting the main agent free-form.

---

## ADR-011: Staged drafting of design docs with scope-approval gate

### Context
The normalisation `/design` run drafted ~1600 lines across 6 documents over ~11 minutes before reaching the first human gate. If the user rejects the scope interpretation or a foundational ADR at that point, the downstream docs (01/02/03/04/05) all need re-drafting — an expensive invalidation.

### Decision
`/design` splits into two drafting phases with a gate between them:

- **Phase 3a** — draft README.md only, including scope, acceptance criteria, conditional-doc decisions, and a one-sentence "approach summary" per planned ADR. User approves scope (Gate 1).
- **Phase 3b** — only after Gate 1, draft 01-architecture / 02-behavior / 03-decisions + approved conditional files. User approves full bundle (Gate 2). Then dispatch architect.

### Trade-offs
- **Gains:** Cheap correction point (~150 lines of README) before expensive drafting (~1500 lines of ADRs). The user sees the approach headline for every ADR before any are written out in full. Small features still finish in one turn — user can approve Gate 1 and Gate 2 back-to-back.
- **Costs:** One extra gate per `/design` run. For tiny features this is friction.

### Alternatives considered
- Single gate at end (original design) — rejected; the cost of rework dominates the cost of an extra gate
- Auto-approve Gate 1 if the command detects a "small" feature — rejected; defining "small" is fraught, and the user can always just type "approve" twice

### Revisit trigger
If Gate 1 feedback is substantive <20% of the time in practice, consider merging Gates 1 and 2 for low-complexity features (flag: `--single-gate`).

---

## ADR-012: ADRs and Open Questions are distinct numbered lists

### Context
The umap `/design` run used `O-1 / O-2 / O-3` for unresolved items in a dedicated section. The normalisation run surfaced unresolved items *as ADRs* (ADR-006, ADR-010, etc.) whose Decision lines said "blocking — needs user input". A downstream reader sees two conventions and can't quickly count "how many things still need a decision?".

### Decision
`03-decisions.md` has two distinct numbered lists:
- **ADRs (ADR-NNN)** — decisions already made. Status = ACCEPTED. Decision line is explicit.
- **Open Questions (OQ-N)** — unresolved items. Explicitly state what blocks a decision and who should decide.

An ADR with a decision still pending is a defect — it should have been written as an OQ.

### Trade-offs
- **Gains:** Reader can count blocking items at a glance. Architect can flag any ADR-without-Decision as a review-blocker.
- **Costs:** Slight rigidity during drafting — the author must classify each item up front. Mitigation: promoting an OQ to an ADR once decided is a one-line edit.

---

## ADR-013: "Deliberate departures" is a required section when synthesis or audit exists

### Context
The normalisation `/design` run surfaced a section called "Deliberate departures from synthesis / audit" listing every upstream recommendation the design did NOT adopt, each with a one-sentence reason (Qn deferred, no deprecation alias, ENTITY_PROFILES not introduced, etc.). This is one of the most valuable artifacts for a reviewer — it names the design's conscious disagreements with the reviewer panel. The umap run didn't include it. Inconsistent.

### Decision
If `docs/{feature}/synthesis.md` OR `docs/{feature}/audit.md` exists upstream, `03-decisions.md` must end with a "Deliberate departures from synthesis/audit" section. Every upstream recommendation the design chose not to adopt gets one line with its rationale. Architect treats a missing section (when upstream synthesis/audit exists) as a NEEDS ITERATION trigger.

### Trade-offs
- **Gains:** Forces the design author to consciously address every upstream recommendation, making silent disagreement visible. Reviewer effort pays off — the panel sees how its recommendations landed.
- **Costs:** Extra drafting step. If the design adopts every upstream recommendation, the section can be a single line ("All upstream recommendations adopted").

### Revisit trigger
If users consistently find the section adds no signal (all designs adopt everything), consider making it optional. Unlikely in practice — divergent reviewers especially tend to raise recommendations the design won't fully adopt.

---

## ADR-014: Meta-layer is single-agent, single-directory

### Context
When two or more features are in flight at the same time and touch overlapping surfaces, each feature's `/design` makes locally-optimal choices that can quietly contradict its siblings. The architect ends up doing the cross-reference manually by reading N independent synthesis docs and refraction-inverting in their head ("the swirl"). Empirically observed during a wave of 5 pathway-explorer features: the same defect surfaced through ≥3 different lenses each time.

The meta-layer formalises the cross-reference. Three structural choices were available:

1. **Many meta agents in parallel** (mirror the per-feature reviewer panel at the portfolio level — meta-bioinf, meta-ml, etc.)
2. **One meta agent that produces all three meta artifacts in one invocation**
3. **One meta agent that produces exactly one meta artifact per command invocation**

And for the artifact directory:

A. **Each per-feature dir gets its own `meta/` subdir** (e.g., `docs/normalisation/meta-relations.md`)
B. **A single sibling directory `docs/_meta/` holds all portfolio artifacts**

### Decision

- **One agent, `meta-architect`** (option 3) — pure orchestrator, one file per invocation.
- **Single sibling directory `docs/_meta/`** (option B) — holds `map.md`, `design.md`, `plan.md`, `deferred.md`. Underscore prefix signals "framework-level, not a feature."

Three commands (`/meta-map`, `/meta-design`, `/meta-plan`) each dispatch the same agent with a different target file.

### Trade-offs

- **Gains:** Mirrors the per-feature command structure (`/map → /design → /plan`), so users learn the meta-layer for free. One agent means one system prompt to maintain. One directory means meta artifacts are inspectable in a single place; consumers (a per-feature `/design`) check exactly one path. The agent's per-invocation file boundary keeps each output focused.
- **Costs:** No specialist meta lenses — meta-architect is a generalist, so meta-level statistical or wet-lab critique is not surfaced separately. (If this matters, a user can always invoke `/review _meta` after pretending `_meta/` is a feature — escape hatch, not a v1 feature.)

### Alternatives considered

- **Option 1 (parallel meta-reviewers)** rejected: doubles every reviewer, and the reviewers' lenses are already applied at the per-feature level. Re-running them at portfolio level surfaces no new information cheaply enough to justify the cost.
- **Option 2 (one invocation produces all three meta files)** rejected: violates the per-command human-gate principle. The user wants to approve `_meta/map.md` before `_meta/design.md` is drafted, just as they approve `map.md` before `synthesis.md`.
- **Option A (per-feature `meta/` subdirs)** rejected: scatters portfolio artifacts; the architect would have to enumerate N directories to read the portfolio state.

### Revisit trigger
If the meta-architect's outputs feel under-critiqued (e.g., MADRs that should have been flagged by a statistician slip through), introduce `/review _meta --as <list>` as an escape hatch and let users pretend `_meta/` is a feature for review purposes.

---

## ADR-015: Deferred items live in `_meta/deferred.md`, not per-feature

### Context
Per-feature `/design` Phase 3b drafts a "Deliberate departures from synthesis/audit" section listing every upstream recommendation the design did NOT adopt (with one-sentence reason). Some of these are permanent rejections; some are time-boxed deferrals ("we'll revisit in round 2"). The deferrals are roadmap items; the rejections are not.

If both classes stay only in the per-feature `03-decisions.md`, an architect re-prioritising work has to read N feature design archives to find what was parked. This is the same swirl problem at the deferral level.

### Decision

- **Permanent rejections** stay in their source doc (the feature's `03-decisions.md § Deliberate departures` for per-feature; `_meta/design.md § Deliberate meta-level rejections` for portfolio) and only there.
- **Time-boxed deferrals** ALSO get appended to `docs/_meta/deferred.md` as a row with full provenance: ID, source feature, source doc, description, reason deferred, cost (S/M/L), depends-on, status (`DEFERRED | IN-CONSIDERATION | REJECTED | ACTIVATED`).
- The file is **append-only** with monotonically increasing IDs (`D-001`, `D-002`, …). The dispatching command reads the current max ID before appending.
- **Two append paths**, sharing one schema:
  - `commands/design.md` Phase 3b — per-feature `/design` appends with `Source feature: <slug>`.
  - `commands/meta-design.md` Phase 6.5 — `/meta-design` appends with `Source feature: _meta` for items pulled from the meta-design's `Deliberate meta-level rejections § Time-boxed sub-items` subsection, deferred `OQ-M-N` items, and MADR Consequences referencing future rounds.
- `/meta-plan` consumes `deferred.md` for its "Parking items" section regardless of which command appended each row.

### Trade-offs

- **Gains:** A single grep across `_meta/deferred.md` answers "what have we parked across the portfolio?". Each row carries enough provenance that an architect can re-prioritise without re-reading the source feature's design. Schema is pinned — table columns don't drift across features.
- **Costs:** Both `/design` and `/meta-design` do extra work (read max ID, append row). The classification (permanent vs time-boxed) is a judgment call the design author has to make; the user can re-classify. Two writers to one file means coordination is convention, not enforcement — but append-only + ID-monotonic protects against silent overwrite.

### Alternatives considered

- **Single `roadmap.md` instead of `deferred.md`** rejected: the file records what's parked, not an active roadmap. "Roadmap" implies a forward plan; "deferred" is honest.
- **Auto-extract from each feature's "Deliberate departures" via a script** rejected: introduces a hidden parser. Plain-text append is inspectable.
- **Classify everything as roadmap (no separate "permanent rejection" category)** rejected: feature designs need a permanent record of "we considered X and chose not to do it" without polluting the portfolio roadmap with non-actionable noise.

### Revisit trigger
If `_meta/deferred.md` grows past ~100 rows and becomes hard to scan, introduce status-filtered views (e.g., `/meta-roadmap` command that reads `deferred.md` and groups by status). Premature today.

---

## ADR-016: No meta-architect gate in v1

### Context
The per-feature pipeline has an automated `architect` gate after `/design`: the design docs are reviewed by a separate agent, which emits a Verdict (`READY FOR IMPLEMENTATION | NEEDS ITERATION | NEEDS DISCUSSION`) parsed by `/plan`. This gate catches consistency/completeness issues without requiring the user to read every line.

The meta-layer could mirror this with a `meta-architect-gate` agent that reviews `_meta/design.md` for internal consistency before MADRs are flipped from PROPOSED to ACCEPTED.

### Decision
**No meta-architect gate in v1.** `/meta-design` is human-gated only — the user reviews MADR headlines at Gate 1 and approves/redirects/rejects directly.

### Trade-offs

- **Gains:** Simpler harness (one fewer agent, one fewer prompt to maintain). The user is already paying close attention at Gate 1 because the MADRs constrain every per-feature design downstream — the human review is more thorough than at per-feature gates.
- **Costs:** A consistency defect inside `_meta/design.md` (e.g., MADR-002 contradicts MADR-005) might slip past the user. If it does, the next per-feature `/design` will surface the contradiction as a meta-conflict and the user re-runs `/meta-design` to resolve.

### Alternatives considered

- **Add a `meta-architect-gate` agent** rejected for now — the cost of failure is bounded (per-feature `/design` catches it on first inheritance), and the user is the most reliable consistency-checker for portfolio-level decisions today.
- **Have `architect` itself do double duty** (reviewing both per-feature designs and meta designs) rejected: their inputs are different shapes; conflating them weakens the per-feature gate.

### Revisit trigger
If across 2–3 portfolio waves, ≥2 MADR-vs-MADR contradictions land in `_meta/design.md` and are only caught after they cascade into per-feature designs, introduce a `meta-architect-gate` agent. Until then, the user is the gate.

---

## ADR-017: Parallel propagation of accepted MADRs via `/meta-apply`

### Context
Empirical observation from the pathway-explorer 5-feature wave: after `/meta-design` Gate 1 approval, propagation of accepted MADRs into per-feature designs was left to N serial `/design` re-runs. Measured cost for N = 5:

- ~300–400K tokens of redundant re-reading. Each `/design` re-reads the same `_meta/design.md`, `_meta/map.md`, and its own feature's full upstream chain (map.md, synthesis.md, review/*.md).
- ~10–15 user gate turns (scope gate + full-bundle gate per `/design`, ~2–3 turns per feature).
- Wall-clock dominated by serial main-agent dispatch — one feature finishes before the next begins.
- Temporal incoherence: features sit in "MADRs accepted but not yet propagated" state for the duration of the serial pass.

The propagation is highly parallel: each feature's revision is independent once MADRs are ACCEPTED. The work is narrow — apply MADR consequences to one feature's existing design. This is structurally identical to `/review --as all` (see ADR-007) at the propagation layer.

### Decision
Add `/meta-apply` as a command that dispatches N `feature-reviser` subagents in parallel (one per in-scope feature with an existing `design/` directory), surfaces a compound diff gate, then dispatches N `architect` subagents in parallel for the architecture-gate step. Keep `/design` unchanged as the serial, careful path for first-drafts, deep re-architecting, and `NEEDS_ITERATION` loops.

### Trade-offs

- **Gains:**
  - ~3× token reduction: sub-agents read only `_meta/design.md` + their own feature's docs, not the main-agent's full upstream chain × N.
  - ~5× wall-clock reduction: Phase 3 and Phase 6 are parallel fan-outs, not serial.
  - 1 compound gate replaces ~2N per-feature gates (scope + full-bundle per serial `/design`).
  - Temporal coherence: every MADR lands across the portfolio in the same turn.
- **Costs:**
  - Compound gate is harder to scan than N focused gates. Mitigated by per-feature summary table + drill-down on request (Phase 5 option 4) + prominent MANUAL flag call-out.
  - Sub-agent scope is narrower than `/design`'s main agent. The `feature-reviser` cannot re-architect; it can only edit-and-cite. If a MADR actually requires re-architecting a feature, the sub-agent flags `MANUAL_REDESIGN_NEEDED` and the user falls back to `/design` for that specific feature. Net: no lost capability, cleanly partitioned by work type.
  - Parallel writes to `_meta/deferred.md` require pre-allocation of per-sub-agent ID blocks. Main agent handles allocation in Phase 2 before dispatch; no cross-agent coordination needed.

### Alternatives considered

- **Extend `/meta-design` with a propagation Phase 8** — rejected. Makes `/meta-design` large and mixes authoring with propagation concerns. Propagation is sometimes wanted independently (e.g., after a partial MADR revision without re-running the full meta-design draft).
- **Auto-trigger `/meta-apply` immediately after `/meta-design` Gate 1** — rejected. Violates the human-gate discipline; the user may want to iterate on MADRs before applying them. Explicit separation preserves the gate.
- **Main-agent sequential propagation (keep one command but run serially)** — rejected. Zero parallelism, zero token savings, defeats the purpose.
- **Per-feature commands `/meta-apply-{slug}`** — rejected. N times the command-catalog size for one semantic operation. A single variadic command is cleaner.
- **Give `feature-reviser` the full `/design` scope (allow re-architecting)** — rejected. Widens the sub-agent's blast radius; speculative writes become likely; compound gate becomes harder to trust. Keep sub-agent narrow; flag MANUAL and let `/design` handle judgment-heavy cases.

### Revisit trigger
If ≥30% of `/meta-apply` invocations return majority-`MANUAL_REDESIGN_NEEDED`, the sub-agent's mechanical-application heuristic is failing — either MADRs are systematically too aggressive (meta-design is doing work that per-feature design should own) or `feature-reviser` is too conservative. Introduce either an escape-hatch for sub-agents to ask `architect` for judgment mid-run, or scale back `/meta-apply` to be opt-in per feature. Until then, MANUAL is a correctness signal, not a defect.

### Implementation pointer
See `commands/meta-apply.md` for the 8-phase command logic, `agents/feature-reviser.md` for the sub-agent's scope contract and output format, and [02-behavior.md § Stage 4c](./02-behavior.md) for the flow diagram.

---

## ADR-018: Gates live where judgment is required, not where execution is merely serialised

### Context
Through the first months of live use, a pattern emerged: users self-describing as meta-architects (concerned with MADRs, ADR content, cross-feature sequencing, architectural continuities) found themselves repeatedly asked to approve *mechanical* artifacts — plan decomposition and per-phase implementation cadence. Real `/status` runs on multi-feature portfolios showed 4+ features stalled at "plan gate pending," awaiting a rubber-stamp that the meta-architect had no useful judgment to apply: the architect verdict on the design was already `READY`, and the plan is derivative of the design.

The original harness treated every stage ending as a gate. Empirically, that over-counts judgment opportunities. Most stages have a *judgment moment* (where the human's input genuinely shapes the outcome) and possibly a *serialisation moment* (where the outcome is already determined and the human is merely the scheduler). Conflating those two produces friction without signal.

### Decision
Gates live where human judgment is required, not where execution is merely serialised:

- **Judgment gates (preserved):** `/map` review, `/review` lens selection, `/synthesize` output, `/design` Phase 3a (scope) and Phase 4 (ADR content + OQ defaults), architect `NEEDS_ITERATION`/`NEEDS_DISCUSSION` verdicts, `/verify` `INCOMPLETE`/`DRIFT`/`NEEDS REVIEW` verdicts, `/meta-design` Gate 1 (MADRs), `/meta-apply` compound edit gate, `/meta-plan` portfolio sequencing, `/implement` stop-at-uncertainty.
- **Serialisation moments (relaxed):** `/plan` status transition (auto-stamped `APPROVED` at draft time; Phase 4 is a glance-check, not a status gate); per-phase cadence in `/implement` (default stays hand-paced, but `--auto` opts into end-to-end execution that halts only on phase failure, stop-at-uncertainty, or completion).

### Trade-offs
- **Gains.** Meta-architects stop seeing mechanical approval prompts at rhythms mismatched with their role. Four-feature portfolios no longer stall at "plan gate pending" for no informational reason. The `/implement ... --auto` mode reduces 12-phase features from 12 manual dispatches to 1. Hand-pacers who *want* fine-grained control keep today's experience unchanged — per-phase is still the default.
- **Costs.** A silent `--auto` run could hide a bad-but-not-failing phase (the test passes, but the code is subtly wrong). Mitigated by `/verify` as a post-run mechanical check and `/review` as a post-run quality check. Neither are new; they become the natural follow-ups to `--auto`.
- **Risk.** The harness has historically framed "every gate is required" as a feature. Relaxing two gates could seed cargo-cult relaxation of other gates that genuinely need judgment. Mitigated by the explicit preserved-gates list above — this ADR is the evidence that the architect considered which gates matter, and this is the boundary.

### Alternatives considered
- **Keep all gates.** Rejected: empirical evidence (portfolio stalls on rubber-stamp gates) overrides the original uniform-gate design. Uniform is simpler; non-uniform with a principled split is more accurate.
- **Default `/implement` to `--auto`.** Rejected: hand-pacing is legitimately the right mode for unfamiliar territory, debugging, or regression-risk phases. Flipping the default silently would surprise users who built muscle memory around per-phase. `--auto` as opt-in preserves the existing contract while opening the fast lane.
- **Remove the `/plan` Phase 4 ack entirely.** Rejected: even auto-stamped, a one-shot glance-check is cheap and catches the occasional wrong phase boundary before implementation momentum starts. Silence-means-proceed keeps it frictionless when there's nothing to flag.
- **Add a config setting** (`gate_mode: strict|relaxed`). Rejected: adds surface area and session-level mutable state. Positional arguments on `/implement` (no flag vs `--auto`) carry the mode explicitly and without persistence.

### Consequences
- `commands/plan.md` Phase 4 reworded; plan README template now defaults to `status: APPROVED`.
- `commands/implement.md` gains `--auto` flag parsing and an end-to-end execution loop with one-liner per-phase reporting; single-phase default unchanged.
- `skills/architecture-first-dev/SKILL.md` gains a gate-placement framing sentence and an anti-pattern entry ("Asking the user to approve plan decomposition — the `/design` architect verdict is the real gate").
- This ADR is the citeable reference for why those edits are not a relaxation of the harness's human-gate discipline, but a sharpening of where that discipline applies.

### Revisit trigger
If users start re-invoking `/implement` per phase anyway after getting familiar with `--auto` (i.e., they prefer the hand-paced rhythm even when the feature is familiar), the default choice is wrong — flip `--auto` to default and add `--paced` as the opt-in. Until then, per-phase-default is the right safe shape.

---

## ADR-019: Iterative review with `.history/` archive and cross-informed regime

### Context
`/review` was designed for a single pass: each reviewer provides an independent lens, the user synthesizes and moves on. Live use (leiden-stability, 2026-04-22) surfaced that users iterate: read the synthesis, lock a design direction, re-map, and want reviewers to reconcile the round-1 divergence under the locked constraint. Today's `/review` silently overwrites prior files, losing the round-1 positions that justified the direction-lock. OQ-2 anticipated this exact question; this ADR resolves it along with a second concern not in OQ-2: round-1 (divergence-surfacing) and iterate (reconciliation) are qualitatively different regimes and should not be silently conflated.

### Decision
Two primitives land together:

1. **`.history/` archive** — `review/<name>.md` and `synthesis.md` always name the latest round. On iterate, the prior round moves to `review/.history/round-N/<name>.md` and `synthesis/.history/round-N.md`. Downstream readers (`/design`, `/architect`, `/plan`, `/implement`, `/verify`, `/diagram`) read the primary files and need zero changes — MADR-H3's "latest upstream artifact" rule composes without modification.
2. **`--iterate` regime** — new flag on `/review`. When set: every reviewer's prior-round output + every sibling reviewer's prior-round output + the current synthesis (incl. `§ User resolutions`) + the map are injected into the dispatch prompt. Each reviewer is required to tag every prior position as UPDATED / STOOD BY / RETIRED in a `## Round N iterate summary` section at the top of the new review. Round-1 without `--iterate` is unchanged (parallel-independent); the two regimes never mix.

### Trade-offs
- **Gains.** Named, non-destructive iteration. Divergence-surfacing (round-1) and reconciliation (iterate) are distinct regimes with distinct dispatch prompts; the harness reflects the epistemic structure instead of conflating it. Position-shift is auditable via the required UPDATED/STOOD-BY/RETIRED tags — a round-2 reviewer cannot silently pretend they always held the round-2 view. Audit trail cheap: `.history/` is lazy (only created on first iterate).
- **Costs.** Iterate reviewer dispatch prompts are 2-3× a round-1 prompt (sibling reviews injected). This is acceptable — iterate rounds are rare (you don't iterate a feature you got right) and each reviewer is doing harder work (reconciliation > fresh assessment). One extra directory per iterated feature.
- **Risk.** Cross-informed round-2 could drift toward premature consensus. Mitigated by explicit dispatch-prompt discipline (refine-or-stand-by; never silently rewrite) and by preserving round-1 in `.history/` so the degree of position-shift is visible after the fact.

### Alternatives considered
- **`_v1`/`_v2` filename suffix** — rejected. Every downstream reader would need max-version logic, and every citation in a design doc would need to pick a version. Worst of the three layouts.
- **Subdir by round, no primary-latest file** (`review/round-1/`, `review/round-2/`) — rejected. Forces downstream readers to learn `glob round-{max}/*.md`, breaks backwards compatibility with `/design` and `/architect`. The primary-file-is-always-latest convention is load-bearing.
- **Always cross-inform, even in round-1** — rejected. Round-1's job is divergence-surfacing; cross-informing on round-1 produces premature-consensus and defeats the parallel-independent design (ADR-002).
- **Make iterate implicit (re-running `/review` auto-iterates if prior exists)** — rejected. Users who hand-edit one reviewer file and want to re-run the others would find their edits silently archived. Explicit `--iterate` preserves intent.

### Consequences
- `commands/review.md` gains `--iterate` flag, Phase 2.5 (archive prior round), and a distinct iterate-dispatch prompt with required position-shift tagging. Rules 6–9 added documenting regime discipline.
- `commands/synthesize.md` gains Phase 1.5 (archive prior synthesis on re-run) and updates the `synth` dispatch prompt to consume the prior synthesis + prior reviewer round when producing an iterate synthesis. Rules 4–5 added.
- `docs/workflows/architect/00-quickstart.md` documents both regimes and their distinct purposes. `skills/architecture-first-dev/SKILL.md` gains a cheat-sheet row for the "locked direction → iterate" flow.
- Reviewer agents (`bioinf`, `wetlab`, `graphic`, `stat`, `divergent`, `ml`) gain a single `## Iterate rounds` subsection that conditionally asks for a `## Round N iterate summary` prepend with UPDATED/STOOD BY/RETIRED tags. Round-1 behaviour is unchanged (the subsection triggers only when the dispatch prompt names the round ≥ 2). This is a minor revision of MADR-H5's original "no agent YAML changes" clause: empirically, reviewer output templates are prescriptive enough that a hardcoded `## Summary`-first template fights a dispatch-prompt "prepend a new section" instruction, especially for Sonnet-backed reviewers. The conditional prepend is narrow — regime-gated on the dispatch prompt, no round-1 leakage — and makes the iterate contract robust rather than fragile. No role YAML changes.
- `agents/synth.md` gains permission to read `.history/` archives when the dispatch prompt names them, a `## Position shifts` section (iterate-only) and a `## User resolutions` section (always-present) in its output template, and a carve-out in the "No invention" hard rule so scribing user chat decisions verbatim is not an "invention." This is a permanent capability change, not regime-leakage: `synth` genuinely has new inputs (history files) and new output sections (User resolutions) that outlast any single iterate run.
- This ADR closes OQ-2; the `synthesis.md` primary-file convention keeps ADR-018 (gate-placement) and MADR-H3 (scribe-on-latest) composable without edit.

### Revisit trigger
If iterate rounds start producing position-shifts where ≥50% of reviewers mark everything UPDATED (i.e., round-2 universally overrides round-1), the locked direction is doing too much work and round-1 is vestigial. Consider whether round-1 should be shorter or lens-filtered. Conversely, if ≥80% are STOOD BY, the iterate round isn't earning its token cost — either round-1 was already good enough or the locked direction is too vague to bind.

---

## Preserved disagreements

None at present. Future decisions may introduce tensions (e.g., if a new architecture baseline bootstrap flow is added); record them here.

---

## Open questions

1. **When to retire stale maps?** Currently no TTL on `map.md`. A feature revisited 6 months later will have an outdated map. Consider a `--check-age` flag on `/map` that warns on old files.
2. ~~**Should `/review` support staged re-review?**~~ **RESOLVED 2026-04-22 by MADR-H5 (see proposal `2026-04-22-iterative-review-MADR.md`) and ADR-019 below.** Adopted: `.history/` archive pattern + cross-informed `--iterate` regime. Rejected: `_v1`/`_v2` filename suffix (pushes max-version computation onto every downstream reader). Prior rounds are preserved under `review/.history/round-N/` and `synthesis/.history/round-N.md`; `review/<name>.md` and `synthesis.md` always name the latest.
3. **Cross-feature architecture baseline?** The `.bak` harness had a `/design-architecture-overall` command that synthesized multiple features' designs into a system-wide future-state doc. Useful for large projects with many concurrent features; deferred until needed.

---

## Constraints — things this harness assumes

- Claude Code as the execution environment (subagent dispatch, slash commands, Agent tool)
- A git-tracked project where `docs/` is an acceptable place for design artifacts
- Users comfortable reading markdown (not just interacting via chat)
- Features small enough that a single map.md captures their scope (very large cross-cutting initiatives may need a project-level map too)
