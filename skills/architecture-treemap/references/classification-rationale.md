# Classification rationale

This file explains the `core / seam / removable` classification system used by
the architecture-treemap renderer and the architecture-first audit cadence.

---

## The three-class model

Every component in `components.json` carries exactly one of three
classifications. The classification is a **judgment call** — it is not
derivable from the code alone, and it is not produced by the extractor. It is
authored by a human architect (or an LLM acting as one) during the
`/audit-slice` + `/synthesize-audit` phases.

### core

A *core* component is load-bearing: the application cannot function without it.
Core components implement the essential data transformations, domain logic, or
protocol handling that defines what the application is. They may have high
fan-in (many callers depend on them) or high cyclomatic complexity (they handle
real domain variation).

**Colour in the treemap:** green family (`#2a7a2a` base, lightened by LOC ratio).

**Audit implication:** core components are the last things to refactor. Seam
them before pruning anything that depends on them. If a core component has a
background-knowledge edge, encoding that invariant explicitly is high-value
work.

### seam

A *seam* component is a deliberately thin boundary: it connects two larger
parts of the system through a well-defined interface (a Protocol, a typed
contract, a declarative configuration table). The seam component itself is
small; its value is in the coupling it prevents.

**Colour in the treemap:** amber family (`#c87800` base).

**Audit implication:** seams are the target state for the *modularize* phase of
a prune-modularize-extend frame. Turning a tangled direct-call or shared-state
coupling into a seam weakens the connascence and makes the two sides
independently changeable. Seams should stay small; if a seam grows in LOC it
may be accumulating logic that belongs in core.

### removable

A *removable* component is a STRIP candidate: it is no longer used by any live
caller, its feature has been superseded, or it is dead-on-arrival (the payload
slot it writes is never read). Removing it simplifies the codebase without
changing observable behavior.

**Colour in the treemap:** red family (`#a02020` base).

**Audit implication:** prune removable components first. The `prune_candidates`
array in `components.json` gives the LOC count, the downstream simplifications
each removal unblocks, and any latent bug the removal fixes. The treemap surface
makes the "architectural prize" visible: a large red tile with background-
knowledge edges pointing into core means pruning it collapses multiple implicit
invariants at once.

---

## How the renderer derives colour

The renderer does not use any `color_hint` field (that field is not in the
schema). Colour is derived from two data points:

1. **Classification** → base hue (`core` = green, `seam` = amber, `removable` = red).
2. **LOC relative to the maximum LOC in the same view** → brightness. A
   component at maximum LOC gets the full base colour; a component near zero
   gets a near-neutral dark tone. This keeps small tiles legible without
   washing out the classification signal.

Physical tiles use the **highest-priority classification among their logical
owners**: `removable > seam > core`. This means a physical file shared between
a core and a removable component shows in the removable red, drawing attention
to the coupling.

---

## The "judgement only" pill

If a logical component has no `smoke_findings`, no `audit_slices`, and no entry
in `prune_candidates`, the side panel displays a small grey "judgement only"
pill. This signals that the component's classification is not backed by any
audited evidence in the current snapshot — the verdict may still be correct, but
the reader should treat it as an unexamined assumption rather than a confirmed
finding.
