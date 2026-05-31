"""Open/closed vocabulary harmonization with a pluggable guard (fixes #8, #10).

Generalized from the project-specific ``cellstate_llm.harmonize_label`` so ONE function
serves both modes:

  - cell-STATE: a curated ``vocab`` + ``synonyms``, ``mode="open"`` (unmapped kept as
    ``Novel:<normalized>`` for discovery), a domain ``guard_fn`` (e.g. lineage reject),
    ``fallback_label="Ambiguous_LowSignal"``.
  - cell-TYPE: ``vocab=[]`` / open ontology, ``novel_prefix=""`` + ``clean_novel=False`` so
    the LLM's free-text label passes through verbatim; no guard.

All policy is arguments; nothing here names a biology. Sourced by the engine from the
active profile.
"""

from __future__ import annotations

import logging
import re
from typing import Callable

log = logging.getLogger(__name__)


def normalize_novel_token(raw: str) -> str:
    """Clean a kept-novel label: drop a trailing ' cells', collapse whitespace/'/'→'_'.

    e.g. ``'DNA Damage Response cells' -> 'DNA_Damage_Response'``,
    ``'Stress / Recovery cells' -> 'Stress_Recovery'``. The RAW form should be preserved
    by the caller (e.g. in a trace JSON) for audit; this only shapes the stored token.
    """
    s = str(raw).strip()
    s = re.sub(r"\s+[Cc]ells?$", "", s).strip()
    s = re.sub(r"\s*/\s*", "_", s)
    s = re.sub(r"\s+", "_", s)
    s = re.sub(r"_+", "_", s).strip("_")
    return s


def is_novel_label(label: str, novel_prefix: str = "Novel:") -> bool:
    """True iff ``label`` is an open-vocab preserved novel state (``<novel_prefix>*``)."""
    return bool(novel_prefix) and isinstance(label, str) and label.startswith(novel_prefix)


def harmonize_label(
    raw: str,
    vocab: list[str],
    synonyms: dict[str, str],
    *,
    guard_fn: Callable[[str], bool] | None = None,
    mode: str = "open",
    fallback_label: str = "Ambiguous_LowSignal",
    novel_prefix: str = "Novel:",
    clean_novel: bool = True,
) -> str:
    """Map a raw LLM label → canonical term, honoring ``mode`` and the optional ``guard_fn``.

    Algorithm:
      1. blank → ``fallback_label``.
      2. exact case-insensitive match against ``vocab`` → the canonical-cased vocab term.
      3. longest-key substring match against ``synonyms`` → its target.
      4. GUARD (both modes): if ``guard_fn(lowered_raw)`` is True, return ``fallback_label``
         (e.g. a lineage restatement that must never become a state). Runs AFTER vocab/
         synonym matching (those are real terms) but BEFORE the open-mode novel fallback.
      5. unmapped, un-guarded:
           - ``open``   → ``f"{novel_prefix}{token}"`` where ``token`` is the normalized
             (``clean_novel=True``) or stripped (``False``) raw — discovery-safe; with
             ``novel_prefix=""`` + ``clean_novel=False`` this is verbatim passthrough (cell-type).
           - ``closed`` → ``fallback_label``.
    """
    if not raw or not str(raw).strip():
        return fallback_label

    lowered = str(raw).strip().lower()

    # 1. exact vocab match (case-insensitive) → canonical casing
    for v in vocab:
        if v.lower() == lowered:
            return v

    # 2. longest-key substring synonym match
    best_key: str | None = None
    best_len = 0
    for key, target in synonyms.items():
        if key in lowered and len(key) > best_len:
            best_key = key
            best_len = len(key)
    if best_key is not None:
        return synonyms[best_key]

    # 3. guard (reject in BOTH modes)
    if guard_fn is not None and guard_fn(lowered):
        log.warning("GUARD-REJECTED (→ %s): %r", fallback_label, raw)
        return fallback_label

    # 4. unmapped, un-guarded → mode-dependent
    if mode == "open":
        token = normalize_novel_token(raw) if clean_novel else str(raw).strip()
        label = f"{novel_prefix}{token}"
        if novel_prefix:
            log.info("OPEN-VOCAB: preserving unmapped label as %r (raw=%r)", label, str(raw).strip())
        return label

    log.warning("UNMAPPED label (closed mode → %s): %r", fallback_label, raw)
    return fallback_label
