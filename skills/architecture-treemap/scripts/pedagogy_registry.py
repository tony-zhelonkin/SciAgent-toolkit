#!/usr/bin/env python3
"""Single source of truth for the treemap's explainer + glossary pedagogy.

This is the sibling of ``metric_registry.py`` for everything that is NOT a
per-component metric: the core/seam/removable CLASSIFICATIONS, the three EDGE
TYPES, the two EVIDENCE CLASSES, the direction encoding, and a concept GLOSSARY.

Why a registry (and why ONE file, not three):
    The first cut hardcoded edge-type help in JavaScript (``EDGE_TYPE_HELP``) and
    classification/evidence help in the template HTML. That is the same
    Connascence-of-Value-smeared-across-files smell the metrics ADR
    (references/metrics-architecture.md) was written to kill: the same prose
    lived in two places and drifted. This registry collapses it to one
    declarative table, serialised into the manifest's ``pedagogy`` block, so the
    renderer draws every ``(?)`` explainer and the glossary as a PURE FUNCTION of
    the file — exactly the invariant the metric descriptors already hold.

    One file rather than three (classification/edge/evidence) is the
    minimal-but-no-less call: the three groups are small, share the same
    explainer contract, and feed one manifest block. Three files would mirror
    three wire-blocks, but there is only one wire-block here, so one file wins.

The explainer contract — every entry carries the same four fields:
    what            definition (one short sentence).
    how             how the claim is established / where it comes from. For
                    metric-anchored entries this names the ACTUAL formula + tool
                    (pulled from metric_registry where a metric backs it).
    teaches         the architectural lesson AND what would bias or fool the
                    signal. The honesty about what fools it IS the pedagogy.
    epistemic_source  one of:
        "measured"             deterministic; an AST/git/tool proved it.
        "metric-anchored"      judgment biased hard by a metric, but still
                               judgment (the metric does not deliver the verdict).
        "requires-your-intent" extrinsic; only a human stating a concern can
                               supply it (a structural pass can NEVER assert it).

    epistemic_source is the load-bearing channel. It tells the reader how much
    of each claim lives in the code vs in someone's head — the spine of the
    whole pedagogy.

Adding/Changing pedagogy is a ONE-PLACE change here; the schema accepts the
``pedagogy`` block as an open object (backward compatible — manifests without it
fall back to the renderer's built-in defaults), and the renderer reads it.
"""

from __future__ import annotations

import metric_registry


# Epistemic-source enum, surfaced as a chip on every classification and edge.
EPISTEMIC_MEASURED = "measured"
EPISTEMIC_METRIC_ANCHORED = "metric-anchored"
EPISTEMIC_REQUIRES_INTENT = "requires-your-intent"

EPISTEMIC_SOURCES = {
    EPISTEMIC_MEASURED: {
        "label": "measured",
        "what": "Deterministic: an AST, git, or metric tool proved this directly.",
        "color": "#3a7a3a",
    },
    EPISTEMIC_METRIC_ANCHORED: {
        "label": "metric-anchored judgment",
        "what": (
            "Judgment biased hard by a metric (import graph, fan-in, churn), "
            "but still judgment — the metric points; it does not decide."
        ),
        "color": "#a07020",
    },
    EPISTEMIC_REQUIRES_INTENT: {
        "label": "requires your intent",
        "what": (
            "Extrinsic: not in the code. Only a human naming a concern (via "
            "/audit-slice) can ground this; a structural pass cannot assert it."
        ),
        "color": "#8a4a8a",
    },
}


# ---------------------------------------------------------------------------
# Classifications — core / seam / removable.
# The structural pass may PROPOSE core/seam as provisional, metric-anchored
# hypotheses (showing the metrics that drove them); it must NEVER assert
# removable (that needs intent). The renderer reads provisional_label to render
# a CANDIDATE marker rather than a hard verdict.
# ---------------------------------------------------------------------------

CLASSIFICATIONS = {
    "core": {
        "what": "Load-bearing: the application cannot function without it.",
        "how": (
            "Proposed from the import graph + metrics: high fan-in, central in "
            "the dependency graph, often high churn or complexity. A structural "
            "pass PROPOSES this; a slice confirms it."
        ),
        "teaches": (
            "Core is the last thing to refactor; seam around it before pruning "
            "its dependents. Fooled by re-export hubs (high fan-in but trivial) "
            "and by product intent — kill a feature and a former-core module "
            "can become removable. 'Essential' is partly a product claim, not a "
            "pure code fact."
        ),
        "epistemic_source": EPISTEMIC_METRIC_ANCHORED,
        "color": "#2a7a2a",
        "provisional_label": "core (structural hypothesis)",
    },
    "seam": {
        "what": "A deliberately thin boundary that connects two larger parts.",
        "how": (
            "Weakly detectable: small LOC + meaningful fan-in + a Protocol/typed "
            "contract HINT at a seam, but DELIBERATENESS is not measurable. A "
            "structural pass can only flag a CANDIDATE; intent confirms it."
        ),
        "teaches": (
            "Seams weaken connascence and make two sides independently "
            "changeable — the target state of a modularize phase. Thinness is "
            "measurable; deliberateness and boundary-ness are not, so seam is "
            "the most judgment-dependent class to DETECT. Always a candidate "
            "('possible extension point — confirm intent'), never an asserted "
            "detection."
        ),
        "epistemic_source": EPISTEMIC_REQUIRES_INTENT,
        "color": "#c87800",
        "provisional_label": "seam? (confirm intent)",
    },
    "removable": {
        "what": "A strip candidate: dead, superseded, or dead-on-arrival.",
        "how": (
            "NOT in the code. Requires a user-stated concern (/audit-slice) that "
            "establishes a feature is gone or a path is unused. A structural or "
            "LLM pass must NEVER auto-classify removable — it would have to "
            "fabricate the intent it was never given."
        ),
        "teaches": (
            "Removable is meaningless without intent: 'unused' is relative to a "
            "product surface only a human defines. Until a slice grounds it, a "
            "removable verdict renders BLANK ('requires your judgment — run "
            "/audit-slice'), never as an LLM guess dressed as a finding."
        ),
        "epistemic_source": EPISTEMIC_REQUIRES_INTENT,
        "color": "#a02020",
        "provisional_label": "requires your judgment — run /audit-slice",
    },
}


# ---------------------------------------------------------------------------
# Edge types — the three coupling categories (= Page-Jones connascence forms).
# Content condensed from references/edge-types.md.
# ---------------------------------------------------------------------------

EDGE_TYPES = {
    "direct-call": {
        "what": "One component names another's identifier (a call/import).",
        "how": (
            "Detected statically by the AST import graph (build_import_graph). "
            "The strongest evidence the tool has: the caller's source literally "
            "contains the callee's name."
        ),
        "teaches": (
            "Connascence of Name — the weakest, normal form of coupling. Not a "
            "smell. Only a problem when the callee is itself removable. "
            "Edge direction is the IMPORT direction (A imports B ⇒ A→B), not "
            "data-flow: a shared config that feeds constants to many callers will "
            "show many arrows IN (callers depend on it) and none out — it is a "
            "dependency SINK even though data flows out of it. "
            "Fooled by re-export hubs and __init__ facades that inflate the count."
        ),
        "epistemic_source": EPISTEMIC_MEASURED,
        "color": "#888",
    },
    "shared-state": {
        "what": "Two components both read/write the same data structure.",
        "how": (
            "Audit-asserted: NOT statically detectable from imports. A human "
            "slice must observe both writers touching one global/singleton/DOM "
            "node without a single-writer contract."
        ),
        "teaches": (
            "Connascence of Identity or Value. Tolerable at small scale; "
            "corrosive as the structure grows — any shape change forces "
            "coordinated edits and out-of-order writes cause silent bugs. "
            "Fooled by assuming co-location implies coordination."
        ),
        "epistemic_source": EPISTEMIC_METRIC_ANCHORED,
        "color": "#3a86c0",
    },
    "background-knowledge": {
        "what": "An unenforced invariant one component assumes about another.",
        "how": (
            "Requires intent: the invariant lives only in a human's head, an "
            "ADR, or a drifting comment. No code, type, or test enforces it, so "
            "nothing static can detect it — only a human who knows the contract."
        ),
        "teaches": (
            "Connascence of Convention/Algorithm — the most expensive form, "
            "always a smell: the compiler, type checker, and test suite cannot "
            "catch a violation. Pruning a removable component that carries one "
            "deletes the invariant entirely — the architectural prize."
        ),
        "epistemic_source": EPISTEMIC_REQUIRES_INTENT,
        "color": "#8a5cb0",
    },
}


# ---------------------------------------------------------------------------
# Evidence classes — how an edge was established (orthogonal to edge type).
# ---------------------------------------------------------------------------

EVIDENCE_CLASSES = {
    "static": {
        "what": "The AST/import graph proved this edge.",
        "how": (
            "Emitted by build_import_graph when an import resolves to a repo "
            "file. Renders as a SOLID stroke."
        ),
        "teaches": (
            "Solid = the tool can prove it. Cannot prove runtime/dynamic "
            "coupling (reflection, dependency injection, string-keyed dispatch)."
        ),
        "epistemic_source": EPISTEMIC_MEASURED,
    },
    "audit-asserted": {
        "what": "A human asserted this edge from a slice.",
        "how": (
            "Authored during /audit-slice + /synthesize-audit. Renders as a "
            "DASHED stroke. As trustworthy as the slice behind it."
        ),
        "teaches": (
            "Dashed = judgment, not proof. Shared-state and background-knowledge "
            "edges are always audit-asserted because they are not statically "
            "detectable; that is precisely why they are the expensive forms."
        ),
        "epistemic_source": EPISTEMIC_METRIC_ANCHORED,
    },
}


# Direction encoding — a display convention, not a claim about the code.
DIRECTION_ENCODING = {
    "what": "Outgoing edges draw full-colour; incoming edges draw lighter.",
    "how": "A renderer convention (opacity + colour-lightness), not data.",
    "teaches": (
        "Direction encodes the DEPENDENCY direction (A imports B ⇒ A→B), NOT "
        "the data-flow direction (where values travel) — the two point opposite "
        "ways. A single-source-of-truth config is a high-fan-in, zero-fan-out "
        "SINK: everyone depends on it (arrows IN), it depends on nothing (no "
        "arrows out). The data flows outward to callers, but the dependency "
        "arrows point inward toward the config. Arrows OUT of a config — i.e., "
        "the config importing its consumers — would be the architectural smell. "
        "Use this to read stability: depend TOWARD more-stable (lower-instability) "
        "modules (Stable-Dependencies Principle)."
    ),
    "epistemic_source": EPISTEMIC_MEASURED,
}


# ---------------------------------------------------------------------------
# Glossary — abstract concepts, each TETHERED to the proxy metric that surfaces
# it in THIS tool. The honesty rule: where a metric is a weak proxy, the entry
# SAYS SO. No concept enters the glossary unless it maps to a metric the tool
# emits (or a classification/edge it renders). Sourced from Page-Jones
# (connascence), Martin (instability / SDP / ADP), Ousterhout (deep/shallow).
# `anchor` names the metric key / classification / edge it points at so the
# renderer can make it a clickable tether.
# ---------------------------------------------------------------------------

GLOSSARY = [
    {
        "term": "Coupling",
        "anchor": "fan_in,fan_out",
        "definition": (
            "How much one component must know about another. Surfaced directly "
            "as fan-in (importers) + fan-out (imports) — counts of import edges."
        ),
        "proxy_honesty": "Direct for static imports; blind to runtime coupling.",
    },
    {
        "term": "Cohesion",
        "anchor": None,
        "definition": (
            "How single-purpose a module is (one reason to change). The tool "
            "does NOT ship a measured cohesion metric: a clean code-intrinsic "
            "proxy (within-file vs cross-file reference ratio) is not available "
            "from the file-level import graph, and len(logical_owners) is "
            "judgment-derived (authored by synthesis), not substrate — presenting "
            "it as 'measured' would dress judgment as measurement."
        ),
        "proxy_honesty": (
            "No measured proxy here by design. Read it qualitatively: a file "
            "owned by many logical components (a synthesis judgment) hints at low "
            "cohesion, but that is a verdict, not a number."
        ),
    },
    {
        "term": "Decoupling",
        "anchor": "seam",
        "definition": (
            "Reducing what two components must know about each other, usually by "
            "inserting a seam (a typed boundary). The seam classification marks "
            "where this has been done deliberately."
        ),
        "proxy_honesty": "Seam-ness is judgment; only thinness (LOC) is measured.",
    },
    {
        "term": "Instability (I)",
        "anchor": "instability",
        "definition": (
            "I = fan_out / (fan_in + fan_out), in [0,1]. I=0 is maximally stable "
            "(a depended-on sink); I=1 is maximally unstable (a leaf)."
        ),
        "proxy_honesty": (
            "Exact for static imports. Abstractness A (the other axis of the "
            "main sequence) is not reliably computable in Python, so the tool "
            "surfaces I but NOT main-sequence distance — and says so."
        ),
    },
    {
        "term": "Stability",
        "anchor": "fan_in",
        "definition": (
            "Resistance to change. Proxied by low instability + high fan-in: "
            "many dependents make a module expensive to change."
        ),
        "proxy_honesty": (
            "Partial — stability is also about change RATE; pair with churn_90d. "
            "A low-I module that churns constantly is not really stable."
        ),
    },
    {
        "term": "SRP / one-reason-to-change",
        "anchor": "churn_90d",
        "definition": (
            "A module should have one reason to change. Proxied by many logical "
            "owners (judgment) plus high churn — many owners changing often "
            "approximates many reasons to change."
        ),
        "proxy_honesty": (
            "Weak proxy: churn counts commits, not distinct REASONS; a mass "
            "reformat inflates it. The owner count is a judgment, not a metric."
        ),
    },
    {
        "term": "Acyclic Dependencies Principle (ADP)",
        "anchor": "cycle",
        "definition": (
            "The dependency graph should have no cycles. Detected exactly by "
            "Tarjan strongly-connected-components over the static import graph; "
            "any component in a cycle of size >1 is flagged."
        ),
        "proxy_honesty": (
            "Strong — SCC is exact for STATIC imports. A cycle formed only "
            "through runtime/dynamic dispatch is invisible to it."
        ),
    },
    {
        "term": "Stable-Dependencies Principle (SDP)",
        "anchor": "instability",
        "definition": (
            "Dependencies should point toward more-stable (lower-I) components. "
            "Read the instability gradient along an edge: depending on a more "
            "unstable module than yourself is the violation."
        ),
        "proxy_honesty": "Proxy via the I gradient; ignores dynamic coupling.",
    },
    {
        "term": "Connascence (Page-Jones)",
        "anchor": "direct-call,shared-state,background-knowledge",
        "definition": (
            "Two components are connascent if a change in one forces a change in "
            "the other. Strength taxonomy (weak→strong): Name → Type → Meaning/"
            "Convention → Position → Algorithm → Execution-Order → Timing → "
            "Identity. The three edge types map onto it: direct-call = Name/Type "
            "(weakest), shared-state = Identity/Value, background-knowledge = "
            "Convention/Algorithm (strongest)."
        ),
        "proxy_honesty": (
            "Name/Type are detected exactly (AST imports). Convention/Algorithm "
            "cannot be measured at all — which is WHY they are the most "
            "expensive forms. Locality and other forms have no proxy here."
        ),
    },
    {
        "term": "Deep vs shallow modules (Ousterhout)",
        "anchor": "loc,fan_in",
        "definition": (
            "A deep module hides a lot of implementation behind a narrow "
            "interface (high LOC, low fan-in is a hint). A shallow module adds "
            "little behind a wide interface (low LOC, high fan-out is a hint)."
        ),
        "proxy_honesty": (
            "Weak proxy: interface NARROWNESS is not directly measured. LOC and "
            "fan counts are only hints — a file can be large and shallow."
        ),
    },
]


# ---------------------------------------------------------------------------
# Manifest serialisation — emit one self-describing `pedagogy` block. The
# renderer reads it the same way it reads metric_descriptors, so the explainers
# and glossary travel inside the file (newer manifest → newer pedagogy renders).
# ---------------------------------------------------------------------------

def _metric_explainers() -> dict:
    """Build explainer entries for each metric, pulling the formula/tool from
    metric_registry so the metric's "how" is never duplicated here.

    Static descriptions of how each base metric is computed (the tool/source);
    derived metrics state their formula. Keyed by metric key so the renderer can
    attach a (?) to each metric badge and panel row.
    """
    # How each metric's value is established (the source / tool / formula).
    how_by_key = {
        "loc": ("Non-blank, non-comment physical lines (stdlib counter).", EPISTEMIC_MEASURED),
        "fan_in": ("Distinct repo modules that import this one (AST import graph).", EPISTEMIC_MEASURED),
        "fan_out": ("Distinct repo modules this one imports (AST import graph).", EPISTEMIC_MEASURED),
        "cyclomatic": ("Mean radon cc_visit complexity over blocks; ABSENT if radon is not installed.", EPISTEMIC_MEASURED),
        "churn_90d": ("Commits touching the file in the window (git log --since).", EPISTEMIC_MEASURED),
        "test_ratio": ("Repo test-LOC / source-LOC, attached to modules a test imports (a rollup proxy).", EPISTEMIC_MEASURED),
        "instability": ("I = fan_out / (fan_in + fan_out); undefined for isolated nodes.", EPISTEMIC_MEASURED),
        "refactor_pressure": ("churn * cyclomatic * max(fan_in,1) / max(test_ratio, 0.05); ABSENT if any input (esp. cyclomatic) is missing.", EPISTEMIC_METRIC_ANCHORED),
    }
    teaches_by_key = {
        "loc": "Size != complexity. A huge file is a prompt to LOOK, not a verdict. Fooled by generated code and long string literals.",
        "fan_in": "High fan-in = widely depended-on, so change is expensive (stability). Fooled by re-export hubs / __init__.",
        "fan_out": "High fan-out = this knows about many things (instability). Fooled by a facade that legitimately wires many parts.",
        "cyclomatic": "Branching complexity: hard to test, hard to reason about. When radon is absent this is UNAVAILABLE, not 0 — do not read a missing value as 'simple'.",
        "churn_90d": "A change magnet: still-evolving or unstable code. Fooled by mass reformatting and renames.",
        "test_ratio": "A repo-level rollup, not per-file truth: import != assertion. Admit the proxy.",
        "instability": "Teaches Stable-Dependencies: depend toward lower-I modules. Ignores dynamic coupling; A (abstractness) is not measured, so no main-sequence distance.",
        "refactor_pressure": "Biases attention to hot/complex/depended/untested code — 'look here first', never 'refactor this'. None if any input is missing (fully absent when radon is).",
    }
    explainers = {}
    for desc in metric_registry.METRIC_REGISTRY:
        how, source = how_by_key.get(desc.key, ("Computed by the extractor.", EPISTEMIC_MEASURED))
        explainers[desc.key] = {
            "what": desc.label,
            "how": how,
            "teaches": teaches_by_key.get(desc.key, ""),
            "epistemic_source": source,
        }
    return explainers


def pedagogy_manifest_block() -> dict:
    """Serialise the full pedagogy registry into the manifest `pedagogy` block.

    The renderer reads this to draw the (?) explainers on classifications, edge
    types, evidence classes, metric badges, the direction legend, and the
    glossary panel. Self-describing: a newer manifest carrying a new edge type
    or glossary entry renders its explainer because the content rode along.
    """
    return {
        "epistemic_sources": EPISTEMIC_SOURCES,
        "classifications": CLASSIFICATIONS,
        "edge_types": EDGE_TYPES,
        "evidence_classes": EVIDENCE_CLASSES,
        "direction_encoding": DIRECTION_ENCODING,
        "metrics": _metric_explainers(),
        "glossary": GLOSSARY,
    }
