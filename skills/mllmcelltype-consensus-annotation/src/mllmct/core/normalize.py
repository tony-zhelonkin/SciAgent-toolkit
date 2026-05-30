"""Shared token-usage shape normalization.

Both capture wrappers (Gemini SDK, OpenRouter HTTP) emit per-model usage dicts. This
module collapses either shape into one self-describing record so trace/cost code is
provider-agnostic. Ported verbatim from the original ``cellstate_obs._normalize_token_usage``.
"""

from __future__ import annotations

from typing import Any


def normalize_token_usage(token_usage: dict[str, Any] | None) -> dict[str, Any]:
    """Normalize a capture's ``usage`` dict into a self-describing record.

    Accepts ``{model: {prompt, output, total, calls}}`` (Gemini capture shape) or
    ``{model: {prompt, output, total, calls, cost_usd}}`` (OpenRouter capture shape,
    which adds a native USD cost) and returns::

        {
          "per_model": {model: {prompt, output, total, n_live_calls, cost_usd}},
          "totals":    {prompt, output, total, n_live_calls, cost_usd},
        }

    Missing/extra keys are tolerated; a falsy input yields all-zero totals. The native
    ``cost_usd`` is carried through so ``report_cost`` can PREFER it over rate-estimated
    cost when present (no pricing assumptions needed for a provider that reports cost).
    """
    per_model: dict[str, dict[str, Any]] = {}
    tot = {"prompt": 0, "output": 0, "total": 0, "n_live_calls": 0, "cost_usd": 0.0}
    for model, usage in (token_usage or {}).items():
        if not isinstance(usage, dict):
            continue
        prompt = int(usage.get("prompt", usage.get("input", 0)) or 0)
        output = int(usage.get("output", usage.get("candidates", 0)) or 0)
        total = int(usage.get("total", 0) or 0)
        if total == 0:
            total = prompt + output
        calls = int(usage.get("calls", usage.get("n_live_calls", 0)) or 0)
        try:
            cost_usd = float(usage.get("cost_usd", 0.0) or 0.0)
        except (TypeError, ValueError):
            cost_usd = 0.0
        per_model[str(model)] = {
            "prompt": prompt, "output": output,
            "total": total, "n_live_calls": calls,
            "cost_usd": cost_usd,
        }
        tot["prompt"] += prompt
        tot["output"] += output
        tot["total"] += total
        tot["n_live_calls"] += calls
        tot["cost_usd"] += cost_usd
    return {"per_model": per_model, "totals": tot}
