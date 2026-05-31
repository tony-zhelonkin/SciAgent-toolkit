"""Token-usage capture + forced determinism (remediation fixes #2-determinism, #7).

mLLMCelltype talks to providers but DISCARDS token usage, and it hardcodes
``temperature=0.7`` with no seed — so by default a run is both un-costable and
non-reproducible. We recover usage and force ``temperature=0, seed=0`` by wrapping the
provider call seam WITHOUT editing the vendored library:

  - Gemini (bare model ids) → monkeypatch ``google.genai.models.Models.generate_content``.
  - OpenRouter (slugs with '/') → monkeypatch ``requests.post`` as seen by
    ``mllmcelltype.providers.openrouter`` (and read native USD cost from the response).

Both context managers expose ``self.usage`` in the shared shape consumed by
``core.normalize`` / ``core.trace`` / ``core.cost``. Both are safe to construct offline
(patch-target imports are guarded). Ported from the original
``cellstate_llm.GeminiTokenCapture`` / ``cellstate_obs.OpenRouterCapture``.
"""

from __future__ import annotations

import json
import logging
from typing import Any

log = logging.getLogger(__name__)

# Deterministic sampling applied to EVERY intercepted provider call. The vendored
# providers hardcode temperature=0.7 with no seed, so re-runs drift. seed is best-effort
# per Google's own field docstring (not an absolute guarantee).
DETERMINISTIC_TEMPERATURE: float = 0.0
DETERMINISTIC_SEED: int = 0


# ---------------------------------------------------------------------------
# Gemini config rewrite (merge, do not clobber)
# ---------------------------------------------------------------------------
def force_deterministic_config(config: Any) -> Any:
    """Return ``config`` with ``temperature=0`` (+ ``seed``), preserving all other fields.

    Handles three shapes the library might pass:
      * a ``GenerateContentConfig`` pydantic BaseModel → ``model_copy(update=...)``
        (or ``copy(update=...)`` on pydantic v1), preserving e.g. ``max_output_tokens``;
      * a plain ``dict`` (``GenerateContentConfigDict``) → copy + key override;
      * ``None`` / anything unrecognised → returned unchanged (never break the call).

    ``seed`` is injected only when the config actually exposes that field.
    """
    if config is None:
        return config

    overrides: dict[str, Any] = {"temperature": DETERMINISTIC_TEMPERATURE}

    # Plain dict config (checked FIRST: dict also has ``.copy``, so it would otherwise
    # be mistaken for the pydantic-BaseModel branch below).
    if isinstance(config, dict):
        merged = dict(config)
        merged["temperature"] = DETERMINISTIC_TEMPERATURE
        merged["seed"] = DETERMINISTIC_SEED
        return merged

    if hasattr(config, "model_copy") or hasattr(config, "copy"):
        field_names: set[str] = set()
        cls = type(config)  # access on the CLASS (instance access is deprecated in pydantic v2.11+)
        if hasattr(cls, "model_fields"):  # pydantic v2
            field_names = set(cls.model_fields)
        elif hasattr(cls, "__fields__"):  # pydantic v1
            field_names = set(cls.__fields__)
        if "seed" in field_names:
            overrides["seed"] = DETERMINISTIC_SEED
        try:
            if hasattr(config, "model_copy"):
                return config.model_copy(update=overrides)
            return config.copy(update=overrides)  # pydantic v1 fallback
        except Exception as exc:  # pragma: no cover - defensive
            log.debug("force_deterministic_config: model_copy failed (%s)", exc)
            return config

    log.debug("force_deterministic_config: unrecognised config type %s — left as-is", type(config))
    return config


class GeminiTokenCapture:
    """Capture token usage AND force deterministic sampling on the Gemini SDK call.

    SOURCE OF TRUTH: ``google.genai`` ``GenerateContentResponse.usage_metadata``. The
    vendored provider reads only ``response.text`` and throws usage away, so wrapping the
    SDK call is the only non-invasive way to recover token counts. The same wrapper
    rewrites the incoming ``config`` to ``temperature=0`` + a fixed ``seed`` (merging over
    the library's config so ``max_output_tokens`` survives).

    Usage::

        with GeminiTokenCapture() as cap:
            result = interactive_consensus_annotation(...)
        token_usage = cap.usage   # {model: {prompt, output, total, calls}}
    """

    def __init__(self) -> None:
        self.usage: dict[str, dict[str, int]] = {}
        self._orig = None
        self._models_cls = None

    def _accumulate(self, model: str, um: Any) -> None:
        if um is None:
            return
        bucket = self.usage.setdefault(
            str(model), {"prompt": 0, "output": 0, "total": 0, "calls": 0}
        )
        bucket["prompt"] += int(getattr(um, "prompt_token_count", 0) or 0)
        bucket["output"] += int(getattr(um, "candidates_token_count", 0) or 0)
        bucket["total"] += int(getattr(um, "total_token_count", 0) or 0)
        bucket["calls"] += 1

    def __enter__(self) -> "GeminiTokenCapture":
        try:
            from google.genai import models as genai_models
        except Exception as exc:  # pragma: no cover - SDK absent
            log.info("GeminiTokenCapture: google-genai not importable (%s) — token capture disabled", exc)
            return self

        self._models_cls = genai_models.Models
        self._orig = genai_models.Models.generate_content
        capture = self

        def _wrapped(self_models, *args, **kwargs):  # noqa: ANN001
            # Force deterministic sampling, merging over the library's config. The
            # vendored provider passes config as a KWARG, but handle a positional config
            # too (signature: generate_content(model, contents, config)).
            if "config" in kwargs:
                kwargs["config"] = force_deterministic_config(kwargs["config"])
            elif len(args) >= 3:
                args = list(args)
                args[2] = force_deterministic_config(args[2])
                args = tuple(args)

            response = capture._orig(self_models, *args, **kwargs)
            try:
                model = kwargs.get("model") or (args[0] if args else "unknown")
                capture._accumulate(model, getattr(response, "usage_metadata", None))
            except Exception as exc:  # pragma: no cover - defensive
                log.debug("GeminiTokenCapture: usage read failed: %s", exc)
            return response

        genai_models.Models.generate_content = _wrapped
        return self

    def __exit__(self, *exc) -> None:
        if self._orig is not None and self._models_cls is not None:
            self._models_cls.generate_content = self._orig
        self._orig = None
        self._models_cls = None


class OpenRouterCapture:
    """Capture per-vendor tokens + native USD cost AND force temp=0 on the OpenRouter
    HTTP path — without editing the vendored library.

    The OpenRouter provider posts via ``requests.post`` and the response parser discards
    the sibling ``usage`` block (tokens + OpenRouter's native USD ``cost``); it also omits
    ``temperature``. We patch ``mllmcelltype.providers.openrouter.requests.post`` (which IS
    the global ``requests.post``): inject ``temperature=0``/``seed=0``/``usage.include`` into
    the outgoing body, then read ``usage`` + ``cost`` off the returned response (re-parsing
    cached bytes, so the library's own ``.json()`` still works).

    ``self.usage`` is ``{model: {prompt, output, total, calls, cost_usd}}`` — the Gemini
    shape plus a native ``cost_usd`` that ``report_cost`` prefers when > 0.
    """

    def __init__(self) -> None:
        self.usage: dict[str, dict[str, Any]] = {}
        self._orig_post = None
        self._requests_mod = None

    def _accumulate(self, model: str, response: Any) -> None:
        """Read usage + native USD cost off a returned requests.Response. Never raises."""
        if not model:
            model = "unknown"
        try:
            payload = response.json()
        except Exception as exc:  # pragma: no cover - defensive
            log.debug("OpenRouterCapture: response.json() failed: %s", exc)
            return
        if not isinstance(payload, dict):
            return

        usage = payload.get("usage") or {}
        if not isinstance(usage, dict):
            usage = {}

        prompt = int(usage.get("prompt_tokens", 0) or 0)
        output = int(usage.get("completion_tokens", 0) or 0)
        total = int(usage.get("total_tokens", 0) or 0)
        if total == 0:
            total = prompt + output

        cost = usage.get("cost")
        if cost is None:
            cost = payload.get("cost")
        try:
            cost_usd = float(cost) if cost is not None else 0.0
        except (TypeError, ValueError):
            cost_usd = 0.0

        bucket = self.usage.setdefault(
            str(model),
            {"prompt": 0, "output": 0, "total": 0, "calls": 0, "cost_usd": 0.0},
        )
        bucket["prompt"] += prompt
        bucket["output"] += output
        bucket["total"] += total
        bucket["calls"] += 1
        bucket["cost_usd"] += cost_usd

    def __enter__(self) -> "OpenRouterCapture":
        try:
            import mllmcelltype.providers.openrouter as _or_mod
        except Exception as exc:  # pragma: no cover - library absent
            log.info("OpenRouterCapture: mllmcelltype OpenRouter provider not importable "
                     "(%s) — usage/temp capture disabled", exc)
            return self

        self._requests_mod = _or_mod.requests
        self._orig_post = self._requests_mod.post
        capture = self
        orig_post = self._orig_post

        def _wrapped_post(*args, **kwargs):  # noqa: ANN002, ANN003
            # (a) Inject temperature=0 (+ seed, + usage.include) into the chat-completions
            #     body BEFORE posting. Body travels as data=json.dumps(body).
            model_slug = ""
            try:
                data = kwargs.get("data")
                if isinstance(data, (str, bytes)):
                    body = json.loads(data)
                    if isinstance(body, dict) and "messages" in body:
                        body["temperature"] = 0
                        body.setdefault("seed", 0)
                        body.setdefault("usage", {"include": True})
                        model_slug = str(body.get("model", "") or "")
                        kwargs["data"] = json.dumps(body)
                elif isinstance(kwargs.get("json"), dict) and "messages" in kwargs["json"]:
                    body = kwargs["json"]
                    body["temperature"] = 0
                    body.setdefault("seed", 0)
                    body.setdefault("usage", {"include": True})
                    model_slug = str(body.get("model", "") or "")
            except Exception as exc:  # pragma: no cover - never break the call
                log.debug("OpenRouterCapture: temp-injection skipped: %s", exc)

            response = orig_post(*args, **kwargs)

            # (b) Read usage + native cost off the RETURNED response (does not consume it).
            try:
                if model_slug:
                    capture._accumulate(model_slug, response)
            except Exception as exc:  # pragma: no cover - never break the call
                log.debug("OpenRouterCapture: usage read failed: %s", exc)
            return response

        self._requests_mod.post = _wrapped_post
        return self

    def __exit__(self, *exc) -> None:
        if self._orig_post is not None and self._requests_mod is not None:
            self._requests_mod.post = self._orig_post
        self._orig_post = None
        self._requests_mod = None
