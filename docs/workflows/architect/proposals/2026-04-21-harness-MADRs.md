---
doc: harness-MADRs
date: 2026-04-21
status: ACCEPTED
author: user + assistant (shadartist@gmail.com session)
scope: SciAgent-toolkit/docs/workflows/architect/ + commands/ + skills/architecture-first-dev/
motivation: Gaps surfaced by the contrast-adapter session under pathway-explorer — non-canonical files appeared (00-chat-decisions.md, proposal.md); scope-approval happened without diagrams; no current-vs-proposed architectural view; verbal binding decisions had no defined home.
---

# Harness MADRs — architecture-first-dev improvements

Four MADRs against the architect harness itself. Each is self-contained; they compose
but can be accepted independently. All four are PROPOSED; status flips to ACCEPTED
on user approval, at which point the referenced harness files are edited.

---

## MADR-H1 — Mermaid diagrams move from Gate 2 to Gate 1

**Status.** ACCEPTED 2026-04-22.

**Context.** `commands/design.md:147,152` require `01-architecture.md` to carry a
module Mermaid and `02-behavior.md` to carry a data-flow Mermaid. Both are drafted
in Phase 3b. Gate 1 (scope approval, Phase 3a) is README-only — prose, ADR
headlines, acceptance criteria. The user is asked to approve architectural scope
for a multi-surface feature (Protocol + adapter + HTML + JS, 16 ADRs) with no
picture. Empirically (contrast-adapter session, 2026-04-21) this produced "I
can't reason about architecture without a diagram" friction at exactly the gate
where corrections are cheapest.

**Options.**
- (a) Leave as-is; diagrams are a Phase 3b concern.
- (b) Add a lightweight module sketch + one-arrow data-flow to the Gate-1 README.
- (c) Move all diagram work to Phase 3a.

**Decision.** (b). README at Gate 1 MUST include two Mermaid blocks, both
rough/scope-level:

1. **Module sketch** — one node per module or surface being added/modified,
   edges showing "calls" / "produces" / "consumes". Coarse: actor-level, not
   function-level. Purpose: lets the user see the *shape* of the addition.
2. **Data-flow sketch** — one arrow per pipeline stage, from raw input to
   rendered output. Purpose: lets the user see *how many hops* the design
   introduces.

Detailed diagrams (full module dependency graph, sequence diagrams, state
transitions) still land in Phase 3b as today; the Gate-1 versions are
intentionally coarse and ~5–10 lines of Mermaid each.

**Trade-offs.**
- Gain: scope review becomes architect-grade, not prose-grade. The highest-ROI
  correction window (~150 lines of README vs ~1500 of full bundle) gains its
  missing affordance.
- Cost: ~5–15 minutes of Mermaid authoring per feature, inside the same
  `/design` Phase 3a token budget. Negligible.
- Risk: premature commitment to a shape that Phase 3b would refine. Mitigated
  by declaring the Gate-1 diagrams as "scope-level, non-binding" in the README
  template — architect does not verify them against 01-architecture.md.

**Harness edits.**
- `commands/design.md` Phase 3a template: add "Module sketch (Mermaid)" and
  "Data-flow sketch (Mermaid)" as required README sections before "Acceptance
  criteria".
- `00-quickstart.md:119,124` tip for `/design` Phase 3a: mention that the
  Gate-1 diagrams are scope-level and will be refined in Phase 3b.
- `skills/architecture-first-dev/SKILL.md`: add one line in the `/design`
  section — "Gate-1 README now carries coarse diagrams; Phase 3b detailed
  diagrams build on them."

---

## MADR-H2 — `01-architecture.md` requires a current-vs-proposed delta view

**Status.** ACCEPTED 2026-04-22.

**Context.** `01-architecture.md` today shows the *proposed* module structure.
It does not show the *current* structure or the *delta*. For a feature that
extends an existing system (the common case), the architect-director needs to
see what changes: which edges are new, which are removed, which modules are
touched. Reading the proposed diagram against `map.md` prose is a mental-diff
exercise that the harness should do for the user.

**Options.**
- (a) Keep as-is; user cross-references `map.md` and `01-architecture.md` by hand.
- (b) Require a side-by-side "Current" / "Proposed" pair of diagrams.
- (c) Require a single diagram with delta annotations (added edges bold,
  removed edges dashed + red, new nodes highlighted via Mermaid `classDef`).

**Decision.** (c). One diagram, annotated. Rationale: two diagrams double the
visual-parsing cost and force the reader to diff manually. A single annotated
diagram keeps the cognitive load low and makes the delta the figure-ground.

**Required subsection** in `01-architecture.md`:

```markdown
## Delta view — what changes

<mermaid graph with:
  classDef new fill:#9f6,stroke:#333;
  classDef removed fill:#f99,stroke:#333,stroke-dasharray:5 5;
  classDef touched fill:#ff9,stroke:#333;
  new edges: style N stroke:#3a3,stroke-width:2px;
  removed edges: style N stroke:#a33,stroke-dasharray:5 5;
>

### Legend
- Green fill / solid green edge: added by this feature
- Yellow fill: existing module, modified by this feature
- Red dashed: removed by this feature
- Unstyled: existing, unchanged
```

**Trade-offs.**
- Gain: "what am I actually approving?" becomes answerable from one figure.
  Reviewers (divergent, architect) can detect over-scoping by counting green
  nodes/edges.
- Cost: 10–20 lines of Mermaid per feature. The raw material already exists
  in `map.md` ("Blast radius") and the proposed-architecture section.
- Risk: delta diagrams drift from reality during implementation. Mitigated by
  `/verify` which can optionally re-render the delta from code (deferred; not
  blocking this MADR).

**Harness edits.**
- `commands/design.md` Phase 3b, 01-architecture.md section: add "Delta view"
  as a required subsection with the Mermaid class-def legend above.
- `docs/workflows/architect/01-architecture.md` (the harness self-doc): add
  a worked example of a delta diagram for one of the existing features.

---

## MADR-H3 — Binding chat decisions scribed to the latest upstream artifact ("scribe-on-latest")

**Status.** ACCEPTED 2026-04-22.

**Context.** Between `/synthesize` and `/design`, users routinely make binding
architectural decisions in chat (reviewing `synthesis.md` and asking "actually
let's retire the tier badge entirely"). The harness has no defined home for
these decisions, so the assistant invents siblings (`00-chat-decisions.md` in
the contrast-adapter session). This pollutes the canonical file set and means
downstream commands (architect, plan, verify) don't know to read them.

"Write into `03-decisions.md`" doesn't work — that file doesn't exist yet when
the chat happens. And chat decisions can happen at *any* pipeline stage, not
just post-synthesis:

| Stage when chat happens | Only upstream artifact(s) available |
|---|---|
| Post-`/map`, pre-`/review` | `map.md` |
| Post-single-reviewer (synth refused) | `review/*.md` |
| Post-`/synthesize` | `synthesis.md` (+ upstream) |
| Mid-`/design` iteration (NEEDS_ITERATION loop) | `design/03-decisions.md` |
| Post-READY, pre-`/plan` | `design/03-decisions.md` (frozen) |

**Options.**
- (a) Ban chat decisions; force the user to re-run commands.
  Rejected: hostile to the "ask, you answer" navigation contract in
  `skills/architecture-first-dev/SKILL.md:12`.
- (b) Single sibling file (`chat-decisions.md`).
  Rejected: adds a canonical file; downstream commands need new read rules;
  sequence information (when was this decided?) is lost.
- (c) New journal file (`decisions-journal.md`) across all stages.
  Rejected: yet another canonical file; duplicates content that belongs in-situ.
- (d) **Scribe-on-latest.** Every canonical doc gets an optional
  `## User resolutions` (or `## Binding chat decisions` for design docs) section.
  The assistant appends binding verbal decisions to the **most recent upstream
  artifact the next command will consume**. Every downstream command reads the
  section and treats it as authoritative.

**Decision.** (d).

**Target-selection rule** (assistant-internal):

| State at time of chat | Append to | Section heading |
|---|---|---|
| `map.md` exists; no reviews | `map.md` | `## User resolutions` |
| Reviews exist; no synthesis | latest `review/*.md` or a new `review/user.md` | `## User resolutions` |
| `synthesis.md` exists | `synthesis.md` | `## User resolutions` |
| `design/03-decisions.md` exists (DRAFT) | `design/03-decisions.md` | `## Binding chat decisions` |
| `design/` READY, pre-`/plan` | `design/03-decisions.md` | `## Post-READY amendments` + re-dispatch architect |

**Cross-stage invalidation rule.** If a chat decision during a *later* stage
logically invalidates an *earlier* frozen artifact, the assistant STOPs and
surfaces the conflict to the user, same pattern as the existing MADR-conflict
stop-gate in `/design` Phase 1. Example: a decision during `/design` that
changes reviewer consensus — assistant must say "this changes synthesis-level
consensus; re-run `/synthesize` or override locally?" before writing.

**Downstream-read rule.** Every command's Phase 1 gains one line:
> "Also read `§ User resolutions` / `§ Binding chat decisions` /
> `§ Post-READY amendments` from each upstream artifact; treat them as
> equal-priority to reviewer consensus and cite them explicitly when they
> shape a decision."

**Trade-offs.**
- Gain: closed canonical file set stays closed. Decisions live next to the
  context that motivated them (audit-honest, grep-friendly). Token-efficient —
  next command already reads those files; no extra I/O or duplication.
- Cost: one-time template update to `map.md`, `review/*.md` (optional section),
  `synthesis.md`, `03-decisions.md`; one-line Phase 1 rule per downstream
  command; one selection-rule block in the skill.
- Risk: the assistant picks the wrong target (e.g., appends to `synthesis.md`
  when it should have appended to `03-decisions.md` mid-design).
  Mitigation: the selection rule is deterministic (keyed on file existence,
  not judgment); the skill documents it explicitly and includes failure-mode
  examples.

**Harness edits.**
- `commands/synthesize.md`: output template gains an always-present but
  initially-empty `## User resolutions` section.
- `commands/map.md`: output template gains the same section.
- `commands/design.md` 03-decisions.md template: add
  `## Binding chat decisions` and `## Post-READY amendments` (both initially
  empty, appended to as the pipeline progresses).
- All downstream commands (`design`, `plan`, `implement`, `verify`, architect
  agent): add Phase-1 read rule.
- `skills/architecture-first-dev/SKILL.md`: add "Scribe-on-latest" subsection
  with the selection table and the cross-stage invalidation rule.

---

## MADR-H4 — `/diagram <slug>` utility command

**Status.** ACCEPTED 2026-04-22.

**Context.** During long design conversations the user wants to see "the
current state of the architecture, as a picture" without triggering `/design`
re-entry, Phase 3b drafting, or an architect gate. Today the only way is to
open `01-architecture.md` and `02-behavior.md` in a Mermaid-aware viewer and
cross-reference. There is no single-call "render me the current picture"
affordance.

**Options.**
- (a) Do nothing; the user opens `.md` files.
- (b) Add `/diagram <slug>` as a read-only utility command, parallel to
  `/status` — read existing artifacts, produce a single consolidated
  `docs/{slug}/diagrams.md` with module sketch, sequence, data-flow, delta
  view. Writes one file; no gate; no sub-agents.
- (c) Extend `/status` to optionally include diagrams.

**Decision.** (b). Rationale: `/status` is portfolio-level (<30s, no file
writes); `/diagram` is feature-level and writes a single artifact. Separate
concerns, separate commands.

**Command contract.**
- Input: `<slug>`, reads `docs/{slug}/map.md`, `docs/{slug}/design/*.md`
  (whichever exist).
- Output: `docs/{slug}/diagrams.md` — consolidated Mermaid blocks (module
  sketch from `01-architecture.md`, data-flow from `02-behavior.md`, delta
  view from `01-architecture.md § Delta view`, sequence diagrams where
  present). Plus a short "missing" list if any expected diagram source is
  absent.
- No ADRs, no prose analysis, no gate. Pure render-and-collect.
- Main-agent only; no sub-agent dispatch. Target wall-clock: <10s.

**Trade-offs.**
- Gain: cheap "show me the shape" call at any point mid-flow. Supports the
  director-level reasoning mode the user asked for.
- Cost: one new command file (`commands/diagram.md`); one line in
  `00-quickstart.md` commands-at-a-glance. No new gate, no new agent.
- Risk: diagrams drift from `.md` sources if `diagrams.md` is hand-edited.
  Mitigation: file carries a `regenerated-from:` frontmatter with source
  file hashes; next `/diagram` run overwrites without warning if hashes match.

**Harness edits.**
- New: `commands/diagram.md` with the contract above.
- `00-quickstart.md` commands-at-a-glance: add `/diagram <slug>` row in the
  per-feature section.
- `skills/architecture-first-dev/SKILL.md`: one-line entry under "What to
  ask the assistant for" — *"Show me the architecture as it stands" → `/diagram <slug>`*.

---

## Cross-MADR consistency

- H1 and H2 compose: Gate-1 coarse diagrams (H1) and Phase-3b delta view (H2)
  share Mermaid as the single diagram format. No separate tooling.
- H3 is orthogonal to H1/H2/H4.
- H4 consumes the output of H1 and H2: richer Gate-1 diagrams + mandatory
  delta view make `/diagram` more useful.

## Adoption order

1. **H3 first** (highest recurrence risk; the contrast-adapter session already
   produced a non-canonical sibling). Ship before the next `/design` invocation.
2. **H1** next (cheap; unblocks director-level Gate-1 review immediately).
3. **H2** (requires a worked example; slightly more authoring).
4. **H4** (nice-to-have; depends on H1+H2 to be maximally useful).

## Approval record

- [x] User approval: 2026-04-22
- [x] Harness edits applied: 2026-04-22
  - H3: `commands/map.md`, `commands/synthesize.md`, `commands/design.md` (03-decisions.md template + Phase-1 read rule + Phase-5 architect-dispatch rule), `commands/plan.md`, `commands/implement.md`, `commands/verify.md`, `commands/architect.md`, `skills/architecture-first-dev/SKILL.md` (selection table + cross-stage invalidation rule)
  - H1: `commands/design.md` Phase 3a (Module sketch + Data-flow sketch), `00-quickstart.md` tips, `SKILL.md`
  - H2: `commands/design.md` Phase 3b (`01-architecture.md § Delta view` required), `docs/workflows/architect/01-architecture.md` worked example, `SKILL.md`, `00-quickstart.md` tips
  - H4: `commands/diagram.md` (new), `roles/architect.yaml` commands list, `00-quickstart.md` (commands-at-a-glance + tips + "what to ask"), `SKILL.md` navigation table, `docs/workflows/architect/01-architecture.md` commands catalog
- [ ] Tested on first real feature post-adoption: feature slug + date
