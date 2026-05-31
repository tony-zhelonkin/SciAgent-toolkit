"""Build label-guard functions + post-hoc label guardrails (generic).

The project encoded one specific guard — "an identity that is FIXED for the dataset (a
lineage) must never be emitted as a cell STATE". Generalized here to:

  - ``make_reject_guard(patterns)`` → a ``guard_fn(lowered_label) -> bool`` that returns
    True when the label matches any reject pattern (whole-word regexes). Threaded into
    ``harmonize_label(..., guard_fn=...)``. Empty patterns → a guard that never fires.
  - ``check_restricted_states(label, group, restricted)`` → warnings when a state is
    assigned to a ``group`` (e.g. lineage) outside its allow-list.

Patterns and the restricted map come from the active profile's ``domain_guards`` — no
biology is hardcoded here.
"""

from __future__ import annotations

import re
from typing import Callable


def make_reject_guard(patterns: list[str]) -> Callable[[str], bool]:
    """Compile ``patterns`` (regex strings, matched against the LOWERCASED label) into a
    ``guard_fn(lowered) -> bool``. Returns True iff any pattern matches. ``[]`` → never fires.
    """
    compiled = [re.compile(p) for p in (patterns or [])]

    def guard(lowered: str) -> bool:
        if not compiled:
            return False
        return any(p.search(lowered) for p in compiled)

    return guard


def check_restricted_states(
    label: str,
    group: str,
    restricted: dict[str, list[str]],
) -> list[str]:
    """Warn if ``label`` is assigned to a ``group`` not in its allow-list.

    ``restricted`` is ``{state_label: [allowed_group, ...]}``. A state with no entry is
    unrestricted (no warning). Returns a list of human-readable warning strings (empty ==
    clean); warnings DO NOT mutate the label — the caller collects them for review.
    """
    warnings: list[str] = []
    if not isinstance(label, str):
        return warnings
    allowed = restricted.get(label)
    if allowed is not None and group not in allowed:
        warnings.append(
            f"RESTRICTED-STATE VIOLATION: {label!r} assigned to group {group!r} "
            f"but is restricted to {allowed} — flagged for human review."
        )
    return warnings
