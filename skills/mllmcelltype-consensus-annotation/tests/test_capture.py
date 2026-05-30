"""#2-determinism / #7 — config rewrite + OpenRouter temp-inject + usage accumulate."""

import json

from mllmct.core.capture import (
    GeminiTokenCapture,
    OpenRouterCapture,
    force_deterministic_config,
)


def test_force_deterministic_dict_preserves_other_fields():
    out = force_deterministic_config({"temperature": 0.7, "max_output_tokens": 4096})
    assert out["temperature"] == 0.0
    assert out["max_output_tokens"] == 4096   # MERGED, not clobbered
    assert out["seed"] == 0


def test_force_deterministic_pydantic_config():
    from google.genai import types
    cfg = types.GenerateContentConfig(temperature=0.7, max_output_tokens=4096)
    out = force_deterministic_config(cfg)
    assert out.temperature == 0.0
    assert out.max_output_tokens == 4096
    assert out.seed == 0


def test_force_deterministic_none_passthrough():
    assert force_deterministic_config(None) is None


def test_gemini_accumulate():
    from types import SimpleNamespace
    cap = GeminiTokenCapture()
    cap._accumulate("gemini-2.5-flash",
                    SimpleNamespace(prompt_token_count=100, candidates_token_count=20,
                                    total_token_count=120))
    cap._accumulate("gemini-2.5-flash",
                    SimpleNamespace(prompt_token_count=50, candidates_token_count=10,
                                    total_token_count=60))
    u = cap.usage["gemini-2.5-flash"]
    assert u["prompt"] == 150 and u["output"] == 30 and u["total"] == 180 and u["calls"] == 2


def test_openrouter_injects_temp_and_reads_usage():
    import mllmcelltype.providers.openrouter as orm

    captured = {}

    class FakeResp:
        def json(self):
            return {"choices": [{"message": {"content": "x"}}],
                    "usage": {"prompt_tokens": 10, "completion_tokens": 5,
                              "total_tokens": 15, "cost": 0.001}}

    def fake_post(*a, **k):
        captured["data"] = k.get("data")
        return FakeResp()

    orig = orm.requests.post
    orm.requests.post = fake_post
    try:
        with OpenRouterCapture() as cap:
            orm.requests.post(
                "https://example/api",
                data=json.dumps({"model": "openai/gpt-5",
                                 "messages": [{"role": "user", "content": "hi"}]}),
            )
        body = json.loads(captured["data"])
        assert body["temperature"] == 0
        assert body["seed"] == 0
        u = cap.usage["openai/gpt-5"]
        assert u["total"] == 15 and u["cost_usd"] == 0.001 and u["calls"] == 1
    finally:
        orm.requests.post = orig
