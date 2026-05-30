"""PanelEvidenceProvider — a domain-NEUTRAL reference EvidenceProvider.

This is the worked example of ``mllmct.core.evidence.EvidenceProvider``. It demonstrates,
on generic/abstract concepts (no real biology), the four evidence-assembly mechanics that
made the original cell-state run trustworthy:

  1. PROGRAM DECODE — opaque program ids (``P1``, ``P7``) are decoded to a human-readable
     ``program[category:genes]`` token via a program-map table, so the LLM never sees a
     bare opaque id; unknown ids render ``(unmapped)``; ids in a drop-category are omitted
     or flagged as ``(artifact)``.
  2. SIGNATURE DE-SATURATION — globally-dominant ("always-on") signatures are excluded
     from the ranked top-K and instead surfaced as compact ``+`` flags, so a constant
     anchor cannot dominate the evidence.
  3. BINNED-AXIS PANEL — only NON-baseline axis bins are shown, capped to the top-N by
     ``|z|``, in a stable canonical order; an optional COMPANION axis is force-rendered
     whenever a designated axis fires (the "never show X without its companion" rule).
  4. CATEGORY + FLAG TOKENS — a dominant categorical token and a thresholded boolean flag.

Every name and threshold is supplied via ``options`` (from ``profile.evidence_options``).
Expected ``evidence_csv`` columns (all configurable):
  ``cluster``, ``dominant_category``, ``dominant_program``, ``signatures_top`` (';'-joined),
  ``companion_bin``, ``flag_metric``, and per-axis ``bin_<axis>`` / ``z_<axis>`` columns.
Optional ``aux_path`` program-map CSV columns: ``program``, ``category``, ``top_genes`` (';'-joined).
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any

import pandas as pd

from ..core.evidence import EvidenceProvider

log = logging.getLogger(__name__)

_DEFAULTS: dict[str, Any] = {
    "axis_order": [],            # canonical axis order; bins read from bin_<axis>
    "max_axes": 6,               # cap on rendered non-baseline axes (ranked by |z|)
    "baseline_bins": ["Mid", "NA", "nan"],
    "companion_axis": None,      # axis rendered WITH an inline companion bin (e.g. "axis=High companion=Low")
    "companion_trigger_axis": None,  # axis whose firing force-renders companion_axis even if capped/baseline
    "companion_bin_col": "companion_bin",
    "bin_prefix": "bin_",
    "z_prefix": "z_",
    "category_col": "dominant_category",
    "program_col": "dominant_program",
    "program_top_genes_n": 5,
    "drop_categories": [],       # program categories treated as artifacts
    "flag_artifact_not_drop": True,  # flag artifact programs vs omit them
    "signatures_col": "signatures_top",
    "always_on_signatures": [],  # excluded from ranked slots, surfaced as flags
    "signature_top_k": 2,
    "flag_metric_col": "flag_metric",
    "flag_threshold": None,      # None ⇒ no flag token
    "flag_labels": ["flagged", "unflagged"],
    "markers_col": "markers_top",
}


class PanelEvidenceProvider(EvidenceProvider):
    """Reference provider; see module docstring. Pure-pandas, offline-testable."""

    def __init__(self, options: dict | None = None, aux_path: str | None = None) -> None:
        merged = dict(_DEFAULTS)
        merged.update(options or {})
        super().__init__(options=merged, aux_path=aux_path)
        self._program_map = self._load_program_map(aux_path)

    # -- program decode (mechanic 1) -------------------------------------
    def _load_program_map(self, path: str | None) -> dict[str, dict]:
        if not path:
            return {}
        p = Path(path)
        if not p.exists():
            log.warning("program-map table missing: %s — programs render (unmapped)", p)
            return {}
        df = pd.read_csv(p)
        out: dict[str, dict] = {}
        n = int(self.options["program_top_genes_n"])
        for _, row in df.iterrows():
            pid = str(row["program"])
            genes = [g.strip() for g in str(row.get("top_genes", "") or "").split(";") if g.strip()][:n]
            out[pid] = {"category": str(row.get("category", "Unknown")), "genes": genes}
        return out

    def _render_program(self, pid: str) -> str | None:
        if not pid or pid in ("NA", "nan", ""):
            return None
        info = self._program_map.get(pid)
        if info is None:
            return f"program={pid}(unmapped)"
        if info["category"] in set(self.options["drop_categories"]):
            if not self.options["flag_artifact_not_drop"]:
                return None
            return f"program={info['category']}:{pid}(artifact)"
        genes = ",".join(info["genes"])
        if genes:
            return f"program={pid}[{info['category']}:{genes}]"
        return f"program={pid}[{info['category']}]"

    # -- signature de-saturation (mechanic 2) ----------------------------
    def _render_signatures(self, raw: Any) -> str:
        if raw is None or (isinstance(raw, float) and pd.isna(raw)):
            raw = ""
        names = [n.strip() for n in str(raw).split(";") if n.strip()]
        always_on = set(self.options["always_on_signatures"])
        k = int(self.options["signature_top_k"])
        ranked = [n for n in names if n not in always_on][:k]
        flags = [f"{n}+" for n in names if n in always_on]
        body = ",".join(ranked) if ranked else "none-discriminative"
        if flags:
            body = f"{body} ({','.join(flags)})"
        return f"signatures={body}"

    # -- binned-axis panel (mechanic 3) ----------------------------------
    def _axis_parts(self, row: pd.Series) -> list[str]:
        order = list(self.options["axis_order"])
        if not order:
            return []
        cap = int(self.options["max_axes"])
        baseline = set(self.options["baseline_bins"])
        bp, zp = self.options["bin_prefix"], self.options["z_prefix"]

        bins = {a: str(row.get(f"{bp}{a}", "NA")) for a in order}
        zs: dict[str, float] = {}
        for a in order:
            try:
                zs[a] = float(row.get(f"{zp}{a}"))
            except (ValueError, TypeError):
                zs[a] = 0.0

        nonbaseline = [a for a in order if bins.get(a, "NA") not in baseline]
        rank = {a: i for i, a in enumerate(order)}
        nonbaseline_sorted = sorted(nonbaseline, key=lambda a: (-abs(zs.get(a, 0.0)), rank.get(a, 999)))
        kept = nonbaseline_sorted[:cap]

        # Companion pairing: when the TRIGGER axis is rendered, force the COMPANION axis in
        # too (so a trigger is never shown without its companion), even if capped/baseline.
        companion = self.options["companion_axis"]
        trigger = self.options["companion_trigger_axis"]
        if companion and companion in order and trigger in kept and companion not in kept:
            kept.append(companion)

        kept_set = set(kept)
        parts: list[str] = []
        comp_col = self.options["companion_bin_col"]
        for a in order:  # stable canonical order
            if a not in kept_set:
                continue
            b = bins.get(a, "Mid")
            if a == companion:
                parts.append(f"{a}={b} companion={row.get(comp_col, 'NA')}")
            else:
                parts.append(f"{a}={b}")
        return parts

    # -- assembly --------------------------------------------------------
    def assemble(self, evidence_csv: str) -> dict[str, str]:
        path = Path(evidence_csv)
        if not path.exists():
            raise FileNotFoundError(f"evidence table missing: {path}")
        df = pd.read_csv(path)
        out: dict[str, str] = {}

        cat_col = self.options["category_col"]
        prog_col = self.options["program_col"]
        sig_col = self.options["signatures_col"]
        flag_col = self.options["flag_metric_col"]
        flag_thr = self.options["flag_threshold"]
        flag_lbls = self.options["flag_labels"]

        for _, row in df.iterrows():
            cid = str(row["cluster"])
            parts: list[str] = []

            if cat_col in df.columns:
                cat = row.get(cat_col)
                parts.append(f"category={'NA' if pd.isna(cat) else cat}")

            parts.extend(self._axis_parts(row))

            if flag_thr is not None and flag_col in df.columns and pd.notna(row.get(flag_col)):
                try:
                    val = float(row.get(flag_col))
                    parts.append(flag_lbls[0] if val > float(flag_thr) else flag_lbls[1])
                except (ValueError, TypeError):
                    pass

            if prog_col in df.columns:
                tok = self._render_program(str(row.get(prog_col, "NA")))
                if tok is not None:
                    parts.append(tok)

            if sig_col in df.columns:
                parts.append(self._render_signatures(row.get(sig_col, "")))

            out[cid] = " | ".join(parts)
        return out
