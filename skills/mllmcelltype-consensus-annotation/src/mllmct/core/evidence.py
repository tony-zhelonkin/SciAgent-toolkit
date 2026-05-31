"""EvidenceProvider ABC — the "inject extra metadata" extension point (fix #7 knob).

Cell-TYPE annotation feeds the LLM marker genes only. Cell-STATE annotation additionally
injects per-cluster EVIDENCE (decoded programs, de-saturated signatures, binned axis
panels, …) so the model reasons over discriminative biology, not just markers. That
evidence assembly is the ONLY genuinely domain-specific step, so it lives behind this
interface: a profile names an ``EvidenceProvider`` and its ``options``; the engine calls
``assemble`` only when ``inject_evidence`` is true.

``core`` ships the ABC and a generic markers loader; a worked reference implementation is
``mllmct.plugins.panel_evidence`` (synthetic/neutral). Other domains write their own provider
and point a profile at it — no edits to ``core`` or the engine.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from importlib import import_module
from pathlib import Path

import pandas as pd

# Column names the generic markers loader will try, in order.
_MARKER_COLUMNS = ("markers", "markers_top", "markers_top15", "marker_genes")


class EvidenceProvider(ABC):
    """Turn a per-cluster evidence table into ``{cluster_id: evidence_string}``.

    Implementations read only files/options handed to them (no global project config), so
    they are unit-testable from fixtures.
    """

    def __init__(self, options: dict | None = None, aux_path: str | None = None) -> None:
        """``options``: provider-specific knobs (from ``profile.evidence_options``).
        ``aux_path``: optional secondary table path (e.g. a program-decode map)."""
        self.options = dict(options or {})
        self.aux_path = aux_path

    @abstractmethod
    def assemble(self, evidence_csv: str) -> dict[str, str]:
        """Return ``{cluster_id: evidence_string}`` (string cluster ids). Empty dict ⇒
        markers-only (the engine then renders no evidence block)."""
        raise NotImplementedError


def load_provider(
    dotted: str,
    options: dict | None = None,
    aux_path: str | None = None,
) -> EvidenceProvider:
    """Instantiate an ``EvidenceProvider`` from a ``"module.path:ClassName"`` string."""
    if ":" not in dotted:
        raise ValueError(f"evidence_provider must be 'module:Class', got {dotted!r}")
    mod_name, cls_name = dotted.split(":", 1)
    cls = getattr(import_module(mod_name), cls_name)
    if not (isinstance(cls, type) and issubclass(cls, EvidenceProvider)):
        raise TypeError(f"{dotted} is not an EvidenceProvider subclass")
    return cls(options=options, aux_path=aux_path)


def load_markers_csv(
    path: str | Path,
    n_markers: int,
    columns: tuple[str, ...] = _MARKER_COLUMNS,
) -> dict[str, list[str]]:
    """Load ``{cluster_id: [top marker symbols]}`` from a CSV.

    Requires a ``cluster`` column and one ';'-joined markers column (the first of
    ``columns`` present). Cluster ids are stringified to match ``adata.obs`` categories.
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"markers table missing: {path}")
    df = pd.read_csv(path)
    if "cluster" not in df.columns:
        raise ValueError(f"{path} has no 'cluster' column")
    marker_col = next((c for c in columns if c in df.columns), None)
    if marker_col is None:
        raise ValueError(f"{path} has none of the markers columns {columns}")
    out: dict[str, list[str]] = {}
    for _, row in df.iterrows():
        cid = str(row["cluster"])
        genes = [g.strip() for g in str(row.get(marker_col, "")).split(";") if g.strip()]
        out[cid] = genes[:n_markers]
    return out
