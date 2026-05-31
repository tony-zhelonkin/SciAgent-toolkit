"""#6/#10 — trace tokens: live capture, cached-re-run non-clobber, cost_summary aggregate."""

import json

from mllmct.core.cost import write_cost_summary
from mllmct.core.trace import TraceWriter


def test_live_capture_then_cached_preserves(tmp_path):
    base = tmp_path / "trace"
    tw = TraceWriter("lensA", base_dir=base / "lensA")

    # live capture → real totals
    tw.write_tokens({"gemini-2.5-flash": {"prompt": 4000, "output": 200, "total": 4200, "calls": 1}})
    rec = json.loads((base / "lensA" / "tokens.json").read_text())
    assert rec["from_live_calls"] is True
    assert rec["totals"]["total"] == 4200

    # fully-cached re-run (empty) → must NOT clobber; prior preserved under last_known
    tw.write_tokens({})
    rec2 = json.loads((base / "lensA" / "tokens.json").read_text())
    assert rec2["cache_hit_no_new_tokens"] is True
    assert rec2["last_known"]["totals"]["total"] == 4200


def test_unavailable_when_no_prior(tmp_path):
    base = tmp_path / "trace"
    tw = TraceWriter("lensB", base_dir=base / "lensB")
    tw.write_tokens({})
    rec = json.loads((base / "lensB" / "tokens.json").read_text())
    assert rec["status"] == "unavailable"


def test_cost_summary_aggregates_two_lenses(tmp_path):
    base = tmp_path / "trace"
    TraceWriter("lensA", base_dir=base / "lensA").write_tokens(
        {"gemini-2.5-flash": {"prompt": 4000, "output": 200, "total": 4200, "calls": 1}})
    TraceWriter("lensB", base_dir=base / "lensB").write_tokens(
        {"gemini-2.5-flash": {"prompt": 1000, "output": 50, "total": 1050, "calls": 1}})
    summary = write_cost_summary(["lensA", "lensB"], ["gemini-2.5-flash"], base_dir=base)
    assert summary["grand_total_tokens"]["total"] == 5250
    assert summary["n_lenses_with_cost"] == 2
    assert summary["grand_total_usd"] > 0


def test_meta_and_prompt(tmp_path):
    base = tmp_path / "trace" / "lensC"
    tw = TraceWriter("lensC", base_dir=base)
    tw.write_prompt("hello prompt")
    tw.write_meta(models=["m"], sources={"markers": "x.csv"}, extra={"profile": "p"})
    assert (base / "prompt.txt").read_text() == "hello prompt"
    meta = json.loads((base / "meta.json").read_text())
    assert meta["lens"] == "lensC" and meta["profile"] == "p"
