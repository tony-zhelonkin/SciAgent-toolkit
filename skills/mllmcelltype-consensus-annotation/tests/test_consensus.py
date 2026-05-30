"""#1 — Python-recomputed consensus metrics: invariants, not float literals."""

import math

from mllmct.core.consensus import compute_consensus_metrics, recompute_lens_metrics


def test_unanimous():
    lbl, cp, ent = compute_consensus_metrics({"a": "X", "b": "X", "c": "X"})
    assert lbl == "X"
    assert cp == 1.0
    assert ent == 0.0


def test_all_different():
    lbl, cp, ent = compute_consensus_metrics({"a": "X", "b": "Y", "c": "Z"})
    assert cp == 1 / 3
    assert math.isclose(ent, math.log2(3))


def test_tie_breaks_lexicographically():
    lbl, cp, ent = compute_consensus_metrics({"a": "Y", "b": "X"})
    assert lbl == "X"  # tie → smallest label
    assert cp == 0.5


def test_ranges():
    lbl, cp, ent = compute_consensus_metrics({"a": "X", "b": "X", "c": "Y"})
    assert 0.0 <= cp <= 1.0
    assert ent >= 0.0


def test_empty():
    lbl, cp, ent = compute_consensus_metrics({})
    assert lbl == "" and math.isnan(cp) and math.isnan(ent)


def test_recompute_over_clusters():
    cp, ent, maj = recompute_lens_metrics({
        "m1": {"0": "X", "1": "Y"},
        "m2": {"0": "X", "1": "Z"},
    })
    assert cp["0"] == 1.0 and ent["0"] == 0.0
    assert cp["1"] == 0.5 and ent["1"] > 0.0
