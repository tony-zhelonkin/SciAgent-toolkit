"""#7 — cost reporting: native-vs-rate preference, fallback rate, fully-cached path."""

from mllmct.core.cost import report_cost


def test_rate_estimate_for_known_model():
    out = report_cost({"gemini-2.5-flash": {"prompt": 1_000_000, "output": 0, "total": 1_000_000}},
                      ["gemini-2.5-flash"])
    assert out["per_model"]["gemini-2.5-flash"]["cost_source"] == "rate_estimate"
    assert out["total_usd"] > 0
    assert out["fully_cached"] is False


def test_native_cost_preferred():
    out = report_cost({"openai/gpt-5": {"prompt": 10, "output": 5, "total": 15, "cost_usd": 0.42}},
                      ["openai/gpt-5"])
    pm = out["per_model"]["openai/gpt-5"]
    assert pm["cost_source"] == "native_provider"
    assert pm["usd"] == 0.42


def test_unknown_model_fallback_rate():
    out = report_cost({"mystery-model": {"prompt": 100, "output": 100, "total": 200}},
                      ["mystery-model"])
    assert out["per_model"]["mystery-model"]["rate_fallback"] is True


def test_fully_cached_never_aborts():
    out = report_cost({}, ["gemini-2.5-flash"])
    assert out["fully_cached"] is True
    assert out["total_usd"] == 0.0
