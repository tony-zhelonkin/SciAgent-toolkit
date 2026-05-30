"""Two-axis join: every branch + canonical ordering invariance."""

import pandas as pd

from mllmct.reconcile import JoinConfig, canonical_two_axis, join_cell, reconcile_table

CFG = JoinConfig(
    distinct_labels=["Stressed"],
    order=["Stressed", "Apoptotic", "Activated"],
    uninformative_labels=["Ambiguous_LowSignal", "NA"],
    conf_floor=0.5,
    fallback_label="Ambiguous_LowSignal",
)


def test_both_uninformative():
    assert join_cell("Ambiguous_LowSignal", 1.0, False,
                     "Ambiguous_LowSignal", 1.0, False, CFG) == ("Ambiguous_LowSignal", "both_uninformative")


def test_fallthrough_b():
    assert join_cell("Ambiguous_LowSignal", 1.0, False, "Activated", 1.0, False, CFG) == ("Activated", "fallthrough_b")


def test_fallthrough_a():
    assert join_cell("Activated", 1.0, False, "NA", 1.0, False, CFG) == ("Activated", "fallthrough_a")


def test_agree():
    assert join_cell("Activated", 1.0, False, "Activated", 1.0, False, CFG) == ("Activated", "agree")


def test_two_axis():
    lab, method = join_cell("Stressed", 1.0, False, "Activated", 1.0, False, CFG)
    assert method == "two_axis"
    assert lab == "Stressed · Activated"


def test_defer_b():
    # A informative but NOT distinct → defer to the finer B
    assert join_cell("Activated", 1.0, False, "Apoptotic", 1.0, False, CFG) == ("Apoptotic", "defer_b")


def test_low_conf_is_uninformative():
    assert join_cell("Stressed", 0.4, False, "Activated", 1.0, False, CFG) == ("Activated", "fallthrough_b")


def test_canonical_order_invariant():
    assert canonical_two_axis("Activated", "Stressed", CFG) == canonical_two_axis("Stressed", "Activated", CFG)
    assert canonical_two_axis("Activated", "Stressed", CFG) == "Stressed · Activated"


def test_reconcile_table_branches(fixtures_dir):
    df = pd.read_csv(fixtures_dir / "axes_synth.csv")
    out = reconcile_table(df, CFG)
    methods = list(out["joint_method"])
    assert methods == ["both_uninformative", "fallthrough_b", "agree",
                       "two_axis", "defer_b", "fallthrough_b"]
    # compound 'Stressed / Apoptotic' (row 6) is resolved + treated as uninformative
    assert out.iloc[5]["joint_label"] == "Activated"
