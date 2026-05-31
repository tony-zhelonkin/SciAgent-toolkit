# Metrics architecture — single source of truth (ADR)

Status: accepted (2026-05-24)

## Context

The treemap stack annotates components with per-component metrics. The first
cut hardcoded the same metric list in three layers:

- the extractor's compute-list (`extract_components.py`),
- the schema's closed `metrics` property list (`components.schema.json`),
- the renderer's badge list and side-panel row list (`treemap.bundle.js`).

Adding one metric meant editing all three in lockstep. That is Connascence of
Value/Algorithm smeared across three files — the exact smell this skill exists
to detect. The metric set is also expected to grow (more base metrics, and
DERIVED signals such as a refactor-pressure score), so the seam needed fixing
before it calcified.

## Decision

Introduce a **metric descriptor registry** as the single source of truth:
`scripts/metric_registry.py`. Each metric is declared once as a
`MetricDescriptor` (key, label, short badge label, unit, `higher_is_better`,
`kind` = base|derived, `render` = badge|hidden, badge threshold/colour, and a
`formula` for derived metrics). The three layers derive from it:

- **Extractor** iterates the registry. It computes base metrics with the
  AST/git handles it owns, then calls `compute_derived()` to layer on derived
  metrics via their declared formulas. It emits the registry into the manifest
  as a top-level `metric_descriptors` block.
- **Schema** stops enumerating a closed property list. `metrics` becomes an
  **open numeric map** (`additionalProperties: {type: number, minimum: 0}`):
  any registered metric validates automatically. The named properties remain in
  the schema as documentation only, not as a gate. A `metric_descriptors`
  block `$def` accepts the serialised registry.
- **Renderer** reads `COMPONENTS_DATA.metric_descriptors` to decide which
  badges to draw, what to label panel rows, and the good/bad badge semantics
  (`higher_is_better`). The rendered HTML is therefore a **pure function of the
  manifest** — the descriptors travel inside the file.

### Why open-map over closed-properties

A closed property list is a second copy of the registry that must be kept in
sync by hand. The whole point is to delete that copy. The open map delegates
"which keys are legal" to the registry (the single source of truth), and the
schema only enforces the invariant that always holds: every metric value is a
non-negative number. Validation rigour is preserved (type + non-negativity);
the brittle coupling is removed.

### Why descriptors in the manifest (not vendored into the renderer)

The manifest is the wire contract and the thing people email around. Emitting
descriptors into it keeps the render self-describing: a manifest produced by a
newer extractor (with a metric the renderer has never heard of) still renders
that metric's badge and panel row correctly, because the semantics rode along
in the file. The renderer keeps a built-in default descriptor table ONLY as a
backward-compatibility fallback for pre-registry manifests that omit the block.

### Base vs derived

- **Base** metrics are computed directly from the substrate (loc, fan_in,
  fan_out, cyclomatic, churn_90d, test_ratio).
- **Derived** metrics are pure functions of base metrics, declared with a
  `formula`. They are computed last, after all base metrics are present, and
  are inserted only when the formula returns non-`None` — so a derived metric
  is never fabricated from partial inputs. `refactor_pressure` landed as the
  first derived metric: `churn_90d * cyclomatic * fan_in / test_ratio`
  (test_ratio floored at a small epsilon so "no tests" amplifies pressure
  rather than dividing by zero). It surfaces the highest-leverage refactor
  targets — hot, complex, depended-upon, under-tested code.

## The "add a metric in one place" property

To add a metric X today, you touch **only `scripts/metric_registry.py`**:
append one `MetricDescriptor`. After that:

- if `kind="base"` with a `compute` hook bound in the extractor, the extractor
  computes it; if `kind="derived"` with a `formula`, it is computed from base
  metrics automatically;
- the schema accepts it (open numeric map) — no property to add;
- the renderer shows its badge and panel row (descriptors are in the manifest)
  — no badge-list or row-list to edit.

The three hardcoded lists that used to drift are gone. The registry is a plain
declarative table — a monolith, deliberately: the goal is a clean internal
seam, not a plugin loader or microservice. Adding a metric is registration, not
infrastructure.

## Consequences

- Backward compatible: manifests using the old named metric keys (and omitting
  `metric_descriptors`) still validate and render — proven by the golden
  fixture's pre-registry sibling.
- The golden fixture now carries a `metric_descriptors` block and one derived
  `refactor_pressure` value to exercise the new path.
- A derived metric whose base inputs are missing (e.g. no radon → no
  cyclomatic → no refactor_pressure) is simply absent. Honest by construction.
