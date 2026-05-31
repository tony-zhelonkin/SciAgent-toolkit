"""Two-axis semantic join (generalized 15b) — deterministic, no LLM.

When each cell carries TWO orthogonal annotations (axis A and axis B, each from its own
lens), collapsing them to one consensus throws away biology. This module joins them per
cell into a single label + an audit trail, using only config from ``profile.join``:

  - both axes uninformative        → fallback label                 (both_uninformative)
  - exactly one informative        → the informative axis label     (fallthrough_a / _b)
  - agree                          → that shared label              (agree)
  - disagree, A is a DISTINCT state → canonically-ordered "A · B"   (two_axis)
  - disagree, A generic            → the (finer) B label            (defer_b)

"Distinct" / "uninformative" / the canonical order are all profile lists, so nothing here
names a biology. Axis A is the one whose ``distinct_labels`` trigger a two-axis keep; axis B
is the finer fallback. Canonical ordering guarantees the unordered pair {A,B} renders ONE
string (fixes the v3 'A · B'/'B · A' duplication). The engine that builds the per-cell axes
table from your data is project-specific and lives in your analysis repo.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import pandas as pd

DEFAULT_SEP = " · "  # U+00B7 middle dot


@dataclass
class JoinConfig:
    distinct_labels: list[str] = field(default_factory=list)
    order: list[str] = field(default_factory=list)            # canonical priority for two-axis naming
    uninformative_labels: list[str] = field(default_factory=list)
    conf_floor: float = 0.5
    fallback_label: str = "Ambiguous_LowSignal"
    sep: str = DEFAULT_SEP

    @classmethod
    def from_profile_join(cls, join: dict[str, Any], fallback_label: str) -> "JoinConfig":
        return cls(
            distinct_labels=list(join.get("distinct_labels", [])),
            order=list(join.get("order", [])),
            uninformative_labels=list(join.get("uninformative_labels", [])),
            conf_floor=float(join.get("conf_floor", 0.5)),
            fallback_label=str(join.get("fallback_label", fallback_label)),
            sep=str(join.get("sep", DEFAULT_SEP)),
        )


def resolve_compound(label: str) -> tuple[str, bool]:
    """A lens may emit a compound consensus like 'A / B' (a split decision). Take the FIRST
    token as representative and flag it so the join treats it as not-fully-trusted."""
    if isinstance(label, str) and "/" in label:
        return label.split("/")[0].strip(), True
    return label, False


def is_uninformative(label: str, conf: float, was_compound: bool, cfg: JoinConfig) -> bool:
    """Uninformative if the label carries no state info, conf < floor, or it was compound."""
    if label in cfg.uninformative_labels:
        return True
    try:
        if float(conf) < cfg.conf_floor:
            return True
    except (TypeError, ValueError):
        return True
    return bool(was_compound)


def canonical_two_axis(a_lbl: str, b_lbl: str, cfg: JoinConfig) -> str:
    """Deterministically ordered two-axis name; {A,B} renders exactly one string.

    Lower index in ``cfg.order`` = left slot; labels not in the order list sort AFTER all
    ordered ones (so a distinct label always precedes a generic one).
    """
    def rank(x: str) -> int:
        return cfg.order.index(x) if x in cfg.order else len(cfg.order)

    a, b = sorted([a_lbl, b_lbl], key=rank)
    return f"{a}{cfg.sep}{b}"


def join_cell(
    a_lbl: str, a_conf: float, a_cmp: bool,
    b_lbl: str, b_conf: float, b_cmp: bool,
    cfg: JoinConfig,
) -> tuple[str, str]:
    """Apply the join rule to one cell. Returns ``(joint_label, method)``."""
    a_bad = is_uninformative(a_lbl, a_conf, a_cmp, cfg)
    b_bad = is_uninformative(b_lbl, b_conf, b_cmp, cfg)

    if a_bad and b_bad:
        return cfg.fallback_label, "both_uninformative"
    if a_bad and not b_bad:
        return b_lbl, "fallthrough_b"
    if b_bad and not a_bad:
        return a_lbl, "fallthrough_a"
    if a_lbl == b_lbl:
        return a_lbl, "agree"
    if a_lbl in cfg.distinct_labels:
        return canonical_two_axis(a_lbl, b_lbl, cfg), "two_axis"
    return b_lbl, "defer_b"


def reconcile_table(
    df: pd.DataFrame,
    cfg: JoinConfig,
    *,
    a_label_col: str = "axis_a_label",
    a_conf_col: str = "axis_a_conf",
    b_label_col: str = "axis_b_label",
    b_conf_col: str = "axis_b_conf",
) -> pd.DataFrame:
    """Add ``joint_label`` / ``joint_method`` / ``joint_conf`` columns to a per-cell axes table.

    Each row is one cell with its two axis labels + confidences. ``joint_conf`` is the weaker
    (min) of the two contributing confidences (conservative). Compound 'A / B' labels are
    resolved per cell and flagged low-confidence.
    """
    joint, method, jconf = [], [], []
    for _, row in df.iterrows():
        a_lbl, a_cmp = resolve_compound(str(row[a_label_col]))
        b_lbl, b_cmp = resolve_compound(str(row[b_label_col]))
        try:
            a_conf = float(row.get(a_conf_col, 0.0))
        except (TypeError, ValueError):
            a_conf = 0.0
        try:
            b_conf = float(row.get(b_conf_col, 0.0))
        except (TypeError, ValueError):
            b_conf = 0.0
        j, m = join_cell(a_lbl, a_conf, a_cmp, b_lbl, b_conf, b_cmp, cfg)
        joint.append(j)
        method.append(m)
        jconf.append(min(a_conf, b_conf))

    out = df.copy()
    out["joint_label"] = joint
    out["joint_method"] = method
    out["joint_conf"] = jconf
    return out
