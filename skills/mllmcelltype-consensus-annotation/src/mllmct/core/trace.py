"""Self-describing per-run trace tree + idempotent token non-clobber (fixes #6, #10).

For any produced label a reviewer can open, under ``<base_dir>/``:

    prompt.txt            full rendered prompt we sent
    model_responses.json  {model: {"raw": ..., "parsed": {...}}}
    discussion.json       cross-model discussion rounds, if any
    tokens.json           per-model usage + totals (preserved across cached re-runs)
    meta.json             timestamp, lens, models, git HEAD, source files
    llmcelltype_debug.log  raw library DEBUG capture (written by core.debug_capture)

Every method is idempotent (overwrites its one target) and defensive (missing input is
recorded as "unavailable", never raised). The token writer NEVER clobbers a prior real
capture with "unavailable" — a fully-cached re-run preserves the original under
``last_known`` (fix #10). Decoupled from project paths: ``base_dir`` is an explicit arg.
Ported from ``cellstate_obs``.
"""

from __future__ import annotations

import json
import logging
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from .normalize import normalize_token_usage

log = logging.getLogger(__name__)


def git_head(cwd: str | Path | None = None) -> str:
    """Best-effort current git HEAD; 'unavailable' if not resolvable."""
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "HEAD"],
            cwd=str(cwd) if cwd else None,
            stderr=subprocess.DEVNULL,
        ).decode().strip()
    except Exception:
        return "unavailable"


def atomic_write_json(path: Path, payload: Any) -> None:
    """Write JSON defensively; never raise on serialization issues."""
    try:
        path.write_text(json.dumps(payload, indent=2, default=str))
    except Exception as exc:  # pragma: no cover - defensive
        log.warning("trace: failed to write %s: %s", path.name, exc)


def read_prior_tokens(path: Path) -> dict | None:
    """Read an existing tokens.json (or None). Never raises."""
    if not path.exists():
        return None
    try:
        return json.loads(path.read_text())
    except Exception:
        return None


def prior_has_real_tokens(prior: dict | None) -> bool:
    """True iff a prior tokens.json holds real (>0) captured token totals.

    Handles both the normalized shape (``totals.total``) and a preserved record
    (``last_known.totals.total``), so a chain of cached re-runs keeps protecting the
    original numbers.
    """
    if not isinstance(prior, dict):
        return False
    for node in (prior, prior.get("last_known")):
        if isinstance(node, dict):
            totals = node.get("totals")
            if isinstance(totals, dict) and int(totals.get("total", 0) or 0) > 0:
                return True
    return False


def extract_token_totals(prior: dict | None) -> dict[str, dict[str, int]] | None:
    """Pull per-model token totals out of a tokens.json record for aggregation.

    Returns ``{model: {prompt, output, total, n_live_calls}}`` from the live record or,
    if this run was cached, from the preserved ``last_known`` block. None when no real
    per-model numbers are present.
    """
    if not isinstance(prior, dict):
        return None
    for node in (prior, prior.get("last_known")):
        if isinstance(node, dict):
            pm = node.get("per_model")
            if isinstance(pm, dict) and pm:
                return pm
    return None


class TraceWriter:
    """Writes a self-describing LLM trace tree for one lens/run under ``base_dir``."""

    def __init__(self, lens: str, base_dir: Path | str, git_cwd: str | Path | None = None) -> None:
        self.lens = lens
        self.dir = Path(base_dir)
        self.dir.mkdir(parents=True, exist_ok=True)
        self._git_cwd = git_cwd

    # -- prompt -----------------------------------------------------------
    def write_prompt(self, prompt: str | list[str] | None) -> Path:
        """Persist the full rendered prompt(s) we sent (we build it, so we log it)."""
        path = self.dir / "prompt.txt"
        if prompt is None:
            text = "unavailable — prompt not supplied to TraceWriter\n"
        elif isinstance(prompt, (list, tuple)):
            text = "\n\n===== PROMPT SEGMENT =====\n\n".join(str(p) for p in prompt)
        else:
            text = str(prompt)
        try:
            path.write_text(text)
        except Exception as exc:  # pragma: no cover - defensive
            log.warning("TraceWriter: failed to write prompt.txt: %s", exc)
        return path

    # -- per-model responses ---------------------------------------------
    def write_model_responses(
        self,
        model_annotations: dict[str, Any] | None,
        raw_responses: dict[str, Any] | None = None,
    ) -> Path:
        """Persist each model's response (parsed labels + optional raw text)."""
        path = self.dir / "model_responses.json"
        model_annotations = model_annotations or {}
        raw_responses = raw_responses or {}
        if not model_annotations and not raw_responses:
            atomic_write_json(path, {"status": "unavailable",
                                     "note": "no model annotations or raw responses provided"})
            return path
        merged: dict[str, dict[str, Any]] = {}
        keys = set(model_annotations) | set(raw_responses)
        for k in keys:
            merged[k] = {
                "raw": raw_responses.get(k, "unavailable — raw response not captured "
                                            "(see llmcelltype_debug.log)"),
                "parsed": model_annotations.get(k, "unavailable"),
            }
        atomic_write_json(path, merged)
        return path

    # -- discussion -------------------------------------------------------
    def write_discussion(self, discussion_logs: dict[str, Any] | None) -> Path:
        """Persist the discussion rounds (``result['discussion_logs']``)."""
        path = self.dir / "discussion.json"
        if not discussion_logs:
            atomic_write_json(path, {"status": "unavailable",
                                     "note": "no discussion rounds (no controversial "
                                             "clusters, or discussion_logs not returned)"})
            return path
        atomic_write_json(path, discussion_logs)
        return path

    # -- tokens -----------------------------------------------------------
    def write_tokens(self, token_usage: dict[str, Any] | None) -> Path:
        """Persist token usage per model + totals — PRESERVING prior real captures (#10).

        Captures only happen on cache-MISS provider calls; the mllmcelltype on-disk cache
        short-circuits before the provider runs, so a re-run dominated by cache hits
        captures ZERO tokens. We must NOT overwrite an earlier cost-bearing tokens.json
        with "unavailable":
          * real capture this run → write, flagged ``from_live_calls: true``;
          * nothing captured but prior real numbers on disk → preserve them under
            ``last_known``, flag ``cache_hit_no_new_tokens: true``;
          * nothing captured and no prior data → honest "unavailable".
        """
        path = self.dir / "tokens.json"
        captured = normalize_token_usage(token_usage)
        prior = read_prior_tokens(path)

        if captured["totals"]["total"] > 0 or captured["totals"]["n_live_calls"] > 0:
            payload = dict(captured)
            payload["from_live_calls"] = True
            payload["cache_hit_no_new_tokens"] = False
            payload["note"] = (
                f"captured from {captured['totals']['n_live_calls']} live "
                f"(cache-miss) provider call(s) this run"
            )
            atomic_write_json(path, payload)
            return path

        if prior_has_real_tokens(prior):
            payload = {
                "status": "cached_no_new_tokens",
                "cache_hit_no_new_tokens": True,
                "from_live_calls": False,
                "note": ("this run made 0 live (cache-miss) calls — all clusters served "
                         "from the on-disk cache, so no new usage was produced. The "
                         "original cost-bearing capture is preserved below under 'last_known'."),
                "last_known": prior,
            }
            atomic_write_json(path, payload)
            return path

        atomic_write_json(path, {
            "status": "unavailable",
            "cache_hit_no_new_tokens": True,
            "from_live_calls": False,
            "note": ("token usage not captured this run (0 live calls — fully cached or "
                     "capture disabled) and no prior real capture on disk; cost is not "
                     "recoverable. Re-run with a cleared cache to force cache-miss calls."),
        })
        return path

    # -- meta -------------------------------------------------------------
    def write_meta(
        self,
        models: list[str] | None = None,
        sources: dict[str, str] | None = None,
        extra: dict[str, Any] | None = None,
    ) -> Path:
        """Persist run metadata (timestamp, lens, models, git HEAD, source files)."""
        path = self.dir / "meta.json"
        payload: dict[str, Any] = {
            "lens": self.lens,
            "timestamp_utc": datetime.now(timezone.utc).isoformat(),
            "git_head": git_head(self._git_cwd),
            "models": models or "unavailable",
            "sources": sources or "unavailable",
        }
        if extra:
            payload.update(extra)
        atomic_write_json(path, payload)
        return path
