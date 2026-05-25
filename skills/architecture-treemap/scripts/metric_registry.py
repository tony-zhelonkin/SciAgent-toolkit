#!/usr/bin/env python3
"""Single source of truth for the per-component metric model.

Every metric the treemap stack knows about is declared ONCE here, as a
``MetricDescriptor``. The three layers that touch metrics each derive from
this table instead of hardcoding their own list:

  * the extractor (``extract_components.py``) iterates the BASE descriptors to
    decide what to compute, and the DERIVED descriptors to compute formulas;
  * the schema accepts any numeric metric key (open map) — the descriptor's
    ``label``/``unit`` no longer need a closed property list to validate;
  * the renderer reads the descriptors emitted into the manifest's
    ``metric_descriptors`` block, so the badge set and the side-panel rows are
    a pure function of the file, not of hardcoded JavaScript.

Adding a metric is a ONE-PLACE change: append a ``MetricDescriptor`` here.
If it is ``kind="base"`` the extractor computes it (given a registered
compute hook); if it is ``kind="derived"`` it is computed from base metrics
via ``formula``. The schema already accepts it (open numeric map), and the
renderer shows its badge automatically (descriptors are emitted into the
manifest). No edits to the schema property list, the renderer badge list, or
the extractor compute-list are required — those lists no longer exist.

The descriptor table is intentionally a plain declarative list (a monolith,
not a plugin loader): the goal is a clean internal seam, not infrastructure.
"""

from __future__ import annotations

import dataclasses
from typing import Callable, Optional


@dataclasses.dataclass(frozen=True)
class MetricDescriptor:
    """One metric's complete declaration — the single source of truth.

    Fields:
      key              wire-format key under ``metrics`` (e.g. "fan_in").
      label            full human label for the side-panel breakdown.
      short_badge_label compact label for the on-tile badge (kept short so it
                       fits a ~40px tile); None -> no badge ever drawn.
      unit             display unit suffix ("", "loc", "commits", ...).
      value_kind       "int" or "number" — drives schema-side numeric typing.
      higher_is_better whether a high value is GOOD (True) or a smell (False).
                       Drives badge semantics: a badge lights when a value is
                       on the WORSE side of ``badge_threshold``.
      kind             "base" (extractor computes it directly) or "derived"
                       (computed from other metrics via ``formula``).
      render           "badge" (eligible for an on-tile badge) or "hidden"
                       (carried in the manifest + side panel, but no badge).
      badge_threshold  numeric threshold at which the badge lights; None -> no
                       badge regardless of ``render``.
      badge_color      hex colour for the badge chip.
      compute          OPTIONAL callable the extractor calls for a base metric.
                       Bound by the extractor at runtime (it owns AST/git
                       handles), so it is left None here and injected. Derived
                       metrics ignore it and use ``formula`` instead.
      formula          OPTIONAL callable(metrics_dict)->Optional[number] for a
                       derived metric. Returns None when inputs are missing so
                       the metric is simply absent (never fabricated).
      depends_on       base metric keys a derived metric reads (documentation +
                       a cheap presence check before calling ``formula``).
    """

    key: str
    label: str
    short_badge_label: Optional[str]
    unit: str
    value_kind: str  # "int" | "number"
    higher_is_better: bool
    kind: str  # "base" | "derived"
    render: str  # "badge" | "hidden"
    badge_threshold: Optional[float] = None
    badge_color: str = "#666"
    compute: Optional[Callable] = None
    formula: Optional[Callable] = None
    depends_on: tuple = ()

    def badge_lights(self, value) -> bool:
        """True if ``value`` is on the worse side of the badge threshold.

        For higher_is_better metrics (e.g. test_ratio) the badge is a WARNING
        that lights when the value falls BELOW the threshold. For
        lower-is-better metrics (fan_out, cyclomatic, churn, refactor_pressure)
        it lights when the value rises AT OR ABOVE the threshold.
        """
        if self.render != "badge" or self.badge_threshold is None:
            return False
        if value is None:
            return False
        if self.higher_is_better:
            return value < self.badge_threshold
        return value >= self.badge_threshold


# ---------------------------------------------------------------------------
# Derived-metric formulas
# ---------------------------------------------------------------------------

def _instability(metrics: dict):
    """Instability I = fan_out / (fan_in + fan_out)  (Robert C. Martin).

    A deterministic coupling metric in [0, 1]:
      * I = 0  — maximally STABLE: many modules depend on this one, it depends
                 on nothing (a pure sink). Hard/expensive to change.
      * I = 1  — maximally UNSTABLE: this depends on many, nothing depends on
                 it (a pure source/leaf). Cheap to change.

    Teaches the Stable-Dependencies Principle: dependencies should point toward
    MORE-stable (lower-I) modules. Both inputs come from the same AST import
    graph that produces fan_in / fan_out, so this is pure substrate (no new
    tool, no judgment).

    Returns None when both fan_in and fan_out are absent OR their sum is zero
    (an isolated node has no defined instability — guard the div-by-zero rather
    than fabricate a 0 or 1).
    """
    fan_in = metrics.get("fan_in")
    fan_out = metrics.get("fan_out")
    if fan_in is None and fan_out is None:
        return None
    fan_in = fan_in or 0
    fan_out = fan_out or 0
    total = fan_in + fan_out
    if total == 0:
        return None  # isolated node: instability is undefined, not 0.
    return round(fan_out / total, 3)


def _refactor_pressure(metrics: dict):
    """refactor-pressure ~= churn * complexity * fan_in / test_ratio.

    A single scalar that surfaces "look here first": code that changes a lot
    (churn), is hard to reason about (cyclomatic), is depended on by many
    (fan_in), and is under-tested (low test_ratio) is the highest-pressure
    refactor target. Returns None if any required base metric is absent so the
    derived value is never fabricated from partial data.

    test_ratio is floored at a small epsilon so an untested-but-otherwise-hot
    component does not divide-by-zero; the floor makes "no tests" amplify
    pressure rather than crash it.
    """
    churn = metrics.get("churn_90d")
    cyclo = metrics.get("cyclomatic")
    fan_in = metrics.get("fan_in")
    if churn is None or cyclo is None or fan_in is None:
        return None
    test_ratio = metrics.get("test_ratio")
    denom = test_ratio if (test_ratio is not None and test_ratio > 0.05) else 0.05
    pressure = (churn * cyclo * max(fan_in, 1)) / denom
    return round(pressure, 2)


# ---------------------------------------------------------------------------
# THE REGISTRY — add a metric by appending one descriptor here.
# ---------------------------------------------------------------------------

METRIC_REGISTRY = [
    MetricDescriptor(
        key="loc",
        label="LOC",
        short_badge_label=None,  # size is encoded as tile area, not a badge
        unit="loc",
        value_kind="int",
        higher_is_better=False,
        kind="base",
        render="hidden",
    ),
    MetricDescriptor(
        key="fan_in",
        label="Fan-in",
        short_badge_label=None,  # shown in panel; not a standalone badge
        unit="",
        value_kind="int",
        higher_is_better=False,
        kind="base",
        render="hidden",
    ),
    MetricDescriptor(
        key="fan_out",
        label="Fan-out",
        short_badge_label="hi fan-out",
        unit="",
        value_kind="int",
        higher_is_better=False,
        kind="base",
        render="badge",
        badge_threshold=5,
        badge_color="#6060a0",
    ),
    MetricDescriptor(
        key="cyclomatic",
        label="Cyclomatic",
        short_badge_label="hi CC",
        unit="",
        value_kind="number",
        higher_is_better=False,
        kind="base",
        render="badge",
        badge_threshold=10,
        badge_color="#a06020",
    ),
    MetricDescriptor(
        key="churn_90d",
        label="Churn 90d",
        short_badge_label="hi churn",
        unit="commits",
        value_kind="int",
        higher_is_better=False,
        kind="base",
        render="badge",
        badge_threshold=10,
        badge_color="#405080",
    ),
    MetricDescriptor(
        key="test_ratio",
        label="Test ratio",
        short_badge_label="low tests",
        unit="",
        value_kind="number",
        higher_is_better=True,
        kind="base",
        render="badge",
        badge_threshold=0.3,
        badge_color="#a04040",
    ),
    # ---- derived ----
    MetricDescriptor(
        key="instability",
        label="Instability (I)",
        short_badge_label=None,  # no badge: neither high nor low I is a smell
        unit="",
        value_kind="number",
        higher_is_better=False,  # nominal; no threshold so semantics are inert
        kind="derived",
        render="hidden",
        badge_threshold=None,
        formula=_instability,
        depends_on=("fan_in", "fan_out"),
    ),
    MetricDescriptor(
        key="refactor_pressure",
        label="Refactor pressure",
        short_badge_label="hi pressure",
        unit="",
        value_kind="number",
        higher_is_better=False,
        kind="derived",
        render="badge",
        badge_threshold=400,
        badge_color="#902020",
        formula=_refactor_pressure,
        depends_on=("churn_90d", "cyclomatic", "fan_in", "test_ratio"),
    ),
    MetricDescriptor(
        key="refactor_pressure_loc_proxy",
        label="Refactor pressure (LOC proxy — radon absent)",
        short_badge_label=None,  # degraded proxy: panel/lens only, never a tile badge
        unit="",
        value_kind="number",
        higher_is_better=False,
        kind="derived",
        render="hidden",
        badge_threshold=None,
        # No formula here: this DEGRADED proxy is computed by the extractor only
        # when radon is absent (so real cyclomatic is missing), substituting LOC
        # for complexity. It is emitted with a distinct key and label so the
        # renderer can show it CLEARLY MARKED as degraded, never conflated with
        # the real refactor_pressure. compute_derived skips it (no formula).
        depends_on=("churn_90d", "loc", "fan_in", "test_ratio"),
    ),
]

# Fast lookups
REGISTRY_BY_KEY = {d.key: d for d in METRIC_REGISTRY}
BASE_METRICS = [d for d in METRIC_REGISTRY if d.kind == "base"]
DERIVED_METRICS = [d for d in METRIC_REGISTRY if d.kind == "derived"]


def refactor_pressure_loc_proxy(metrics: dict):
    """DEGRADED refactor-pressure when radon is absent (no cyclomatic).

    Substitutes loc/100 for the missing cyclomatic factor so the "look here
    first" lens still ranks something instead of going blank. Returns None if
    the genuine cyclomatic-based pressure is computable (radon present) or if
    churn is missing. The caller is responsible for emitting this under the
    distinct `refactor_pressure_loc_proxy` key and labelling it as degraded —
    it is NOT computed by compute_derived (it has no formula in the registry),
    so it can never be silently mistaken for the real metric.
    """
    if metrics.get("cyclomatic") is not None:
        return None  # real pressure is computable; do not emit the proxy.
    churn = metrics.get("churn_90d")
    fan_in = metrics.get("fan_in")
    loc = metrics.get("loc")
    if churn is None or loc is None:
        return None
    test_ratio = metrics.get("test_ratio")
    denom = test_ratio if (test_ratio is not None and test_ratio > 0.05) else 0.05
    complexity_proxy = loc / 100.0
    pressure = (churn * complexity_proxy * max(fan_in or 1, 1)) / denom
    return round(pressure, 2)


def compute_derived(metrics: dict) -> dict:
    """Return a NEW dict with all derivable derived metrics added.

    Iterates DERIVED descriptors, calls each formula against the base metrics
    present, and inserts the result only when the formula returns non-None.
    Leaves base metrics untouched. Idempotent.
    """
    enriched = dict(metrics)
    for desc in DERIVED_METRICS:
        if desc.formula is None:
            continue
        value = desc.formula(enriched)
        if value is not None:
            enriched[desc.key] = value
    return enriched


def descriptor_manifest_block() -> dict:
    """Serialise the registry into the manifest ``metric_descriptors`` block.

    Only the fields the renderer needs are emitted (no Python callables). This
    makes the rendered HTML a pure function of the manifest: the badge set, the
    side-panel rows, and the good/bad semantics all travel inside the file.
    """
    block = {}
    for desc in METRIC_REGISTRY:
        block[desc.key] = {
            "label": desc.label,
            "short_badge_label": desc.short_badge_label,
            "unit": desc.unit,
            "value_kind": desc.value_kind,
            "higher_is_better": desc.higher_is_better,
            "kind": desc.kind,
            "render": desc.render,
            "badge_threshold": desc.badge_threshold,
            "badge_color": desc.badge_color,
        }
    return block
