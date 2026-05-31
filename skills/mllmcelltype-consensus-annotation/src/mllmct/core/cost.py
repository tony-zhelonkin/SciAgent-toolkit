"""Logging-only USD cost estimate + aggregate summary (remediation fix #7).

Replaces a dead cost guard that depended on token totals no provider ever logged. These
functions NEVER abort, threshold, or call ``sys.exit`` — cost is reported, not gated. When
a provider reports native USD (OpenRouter ``usage.cost``), it is PREFERRED over the
rate-estimate so no pricing assumption is needed. Decoupled from project paths
(``base_dir`` is an explicit arg). Ported from ``cellstate_obs``.
"""

from __future__ import annotations

import logging
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from .normalize import normalize_token_usage
from .trace import atomic_write_json, extract_token_totals, git_head, read_prior_tokens

log = logging.getLogger(__name__)

# USD per 1,000,000 tokens, split input (prompt) vs output (candidates). DOCUMENTED
# ASSUMPTIONS for cost *reporting only* — they never gate a run. Unknown models fall back
# to DEFAULT_RATE. Update against current provider pricing; OpenRouter native cost is
# preferred when present so these mostly matter for the direct-Gemini path.
MODEL_RATES_USD_PER_1M: dict[str, dict[str, float]] = {
    "gemini-3.1-pro-preview": {"input": 1.25, "output": 10.00},
    "gemini-2.5-pro": {"input": 1.25, "output": 10.00},
    "gemini-2.5-flash": {"input": 0.30, "output": 2.50},
    "gemini-3-flash-preview": {"input": 0.30, "output": 2.50},
    "gemini-3.1-flash-lite": {"input": 0.10, "output": 0.40},
}
DEFAULT_RATE_USD_PER_1M: dict[str, float] = {"input": 0.50, "output": 1.50}


def report_cost(token_usage: dict[str, Any] | None, models: list[str]) -> dict:
    """Best-effort USD cost estimate from token counts. LOGGING ONLY (no abort).

    Returns ``{"total_usd", "per_model", "rates_used", "fully_cached", "n_live_calls"}``,
    or a ``fully_cached`` record when this run made zero live calls (cost of NEW work = $0;
    the prior tokens.json holds the original cost).
    """
    norm = normalize_token_usage(token_usage)
    n_live_calls = norm["totals"]["n_live_calls"]
    has_tokens = norm["totals"]["total"] > 0

    if not has_tokens and n_live_calls == 0:
        log.info("LLM cost this run: $0.00 — FULLY CACHED (0 live calls). New token cost "
                 "is $0; the original cost is preserved per-run in tokens.json (last_known) "
                 "and aggregated in cost_summary.json.")
        return {"fully_cached": True, "n_live_calls": 0, "total_usd": 0.0,
                "per_model": {}, "rates_used": {}}

    per_model: dict[str, dict[str, Any]] = {}
    total_usd = 0.0
    rates_used: dict[str, dict[str, float]] = {}

    for model, usage in norm["per_model"].items():
        in_tok = int(usage.get("prompt", 0) or 0)
        out_tok = int(usage.get("output", 0) or 0)
        native_cost = float(usage.get("cost_usd", 0.0) or 0.0)
        rate = MODEL_RATES_USD_PER_1M.get(model, DEFAULT_RATE_USD_PER_1M)
        rates_used[model] = rate
        est_cost = in_tok * rate["input"] / 1_000_000 + out_tok * rate["output"] / 1_000_000
        use_native = native_cost > 0
        cost = native_cost if use_native else est_cost
        per_model[model] = {
            "input_tokens": in_tok,
            "output_tokens": out_tok,
            "usd": round(cost, 6),
            "cost_source": "native_provider" if use_native else "rate_estimate",
            "rate_fallback": (not use_native) and (model not in MODEL_RATES_USD_PER_1M),
        }
        total_usd += cost

    breakdown = {
        "total_usd": round(total_usd, 6),
        "per_model": per_model,
        "rates_used": rates_used,
        "fully_cached": False,
        "n_live_calls": n_live_calls,
    }

    log.info("LLM cost estimate (logging only, no cutoff): $%.4f total across %d model(s) "
             "from %d live call(s)", total_usd, len(per_model), n_live_calls)
    for model, info in per_model.items():
        note = " [fallback rate]" if info["rate_fallback"] else ""
        log.info("  %s: in=%d out=%d -> $%.4f%s",
                 model, info["input_tokens"], info["output_tokens"], info["usd"], note)
    # Only warn when a model ACTUALLY fell back to a rate estimate (no native cost AND not
    # in the rate table) — a model billed via native provider cost is exact, no caveat needed.
    if any(info["rate_fallback"] for info in per_model.values()):
        log.info("  NOTE: one or more models used the documented fallback rate; verify "
                 "MODEL_RATES_USD_PER_1M against current pricing.")
    return breakdown


def write_cost_summary(
    lenses: list[str],
    models: list[str],
    base_dir: Path | str,
    git_cwd: str | Path | None = None,
) -> dict:
    """Aggregate every lens's LAST-KNOWN token totals → one ``cost_summary.json``.

    For each lens reads ``<base_dir>/<lens>/tokens.json`` and pulls the real per-model
    totals — live OR preserved ``last_known`` if cached this run — then sums tokens + USD
    and writes ``<base_dir>/cost_summary.json``. Makes the cost of a (possibly fully-cached)
    pipeline knowable from disk. Returns the summary dict.
    """
    base_dir = Path(base_dir)
    base_dir.mkdir(parents=True, exist_ok=True)

    per_lens: dict[str, Any] = {}
    grand_per_model: dict[str, dict[str, int]] = {}
    lenses_with_cost = 0
    lenses_missing_cost: list[str] = []

    for lens in lenses:
        tok_path = base_dir / lens / "tokens.json"
        prior = read_prior_tokens(tok_path)
        pm = extract_token_totals(prior)
        cached = bool(isinstance(prior, dict) and prior.get("cache_hit_no_new_tokens"))
        if not pm:
            lenses_missing_cost.append(lens)
            per_lens[lens] = {"status": "no_real_tokens_on_disk",
                              "cache_hit_no_new_tokens": cached}
            continue
        lenses_with_cost += 1
        lens_cost = report_cost(pm, models)
        per_lens[lens] = {
            "per_model_tokens": pm,
            "total_usd": lens_cost.get("total_usd", 0.0),
            "source": "last_known (cached this run)" if cached else "live this run",
        }
        for model, u in pm.items():
            g = grand_per_model.setdefault(
                model, {"prompt": 0, "output": 0, "total": 0, "n_live_calls": 0, "cost_usd": 0.0})
            for k in ("prompt", "output", "total", "n_live_calls"):
                g[k] += int(u.get(k, 0) or 0)
            try:
                g["cost_usd"] += float(u.get("cost_usd", 0.0) or 0.0)
            except (TypeError, ValueError):
                pass

    grand_cost = report_cost(grand_per_model, models) if grand_per_model else {}

    summary = {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "git_head": git_head(git_cwd),
        "models": models,
        "lenses": lenses,
        "n_lenses_with_cost": lenses_with_cost,
        "lenses_missing_cost": lenses_missing_cost,
        "grand_total_tokens": {
            "prompt": sum(u["prompt"] for u in grand_per_model.values()),
            "output": sum(u["output"] for u in grand_per_model.values()),
            "total": sum(u["total"] for u in grand_per_model.values()),
        } if grand_per_model else {"prompt": 0, "output": 0, "total": 0},
        "grand_total_usd": grand_cost.get("total_usd", 0.0),
        "per_model_totals": grand_per_model,
        "per_lens": per_lens,
        "note": ("Aggregates each lens's LAST-KNOWN real token capture (live this run, or "
                 "preserved from the first cache-miss run). USD is a best-effort estimate at "
                 "MODEL_RATES_USD_PER_1M (or native provider cost where reported)."),
    }
    out_path = base_dir / "cost_summary.json"
    atomic_write_json(out_path, summary)
    log.info("Cost summary written: %s — grand total $%.4f across %d/%d lens(es) with cost data",
             out_path, summary["grand_total_usd"], lenses_with_cost, len(lenses))
    return summary
