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


def test_meta_records_the_arbiter_that_resolved_the_disagreements(tmp_path):
    """The arbiter answers every non-unanimous cluster and its answers ship in the delivered
    labels, so a meta.json naming only the annotators describes a smaller panel than the one that
    ran. In the 14782-DM 2026-07-30 run the arbiter was 49.5% of total spend and appears nowhere
    in its own provenance file."""
    tw = TraceWriter("celltype", base_dir=tmp_path)
    meta = json.loads(tw.write_meta(models=["a/x", "b/y"],
                                    consensus_model="c/arbiter").read_text())

    assert meta["consensus_model"] == "c/arbiter"
    assert meta["models"] == ["a/x", "b/y"]


def test_meta_states_an_unset_arbiter_as_a_fact_rather_than_omitting_it(tmp_path):
    """Unset, the library resolves an arbiter of its own and answers anyway. A missing key reads
    as "no arbiter ran"; the recorded string says a model ran and its identity was not captured."""
    tw = TraceWriter("celltype", base_dir=tmp_path)
    meta = json.loads(tw.write_meta(models=["a/x"]).read_text())

    assert "consensus_model" in meta
    assert "unset" in meta["consensus_model"]
