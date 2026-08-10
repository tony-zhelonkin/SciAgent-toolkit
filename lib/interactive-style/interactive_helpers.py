"""
interactive_helpers.py — the SciAgent-toolkit INTERACTIVE-EXPLORER contract (live-kernel side).
================================================================================================
ONE place that owns the "interactive breakpoint explorer" pattern: at a pipeline inflection
point, a Python `.qmd` drives jscatter (jupyter-scatter) linked panels to brush/lasso cells,
extracts a selection into pandas, characterizes it, and persists selected barcodes + a labelled
matplotlib snapshot. It is the LIVE-KERNEL sibling of the (static, R) `decision-gate-notebook`
skill: it FEEDS the decision gate (produces evidence + a barcode selection); it does NOT replace
it. Nothing here computes biology — it is a read-only projection for visualization + selection.

The explorers load COMPACT per-inflection-point parquet tables (2D coords + a few obs columns +
a handful of gene/score columns) written by a project-owned `export_explorers.py` to
`03_results/interactive/`. That keeps the live kernel light — it never opens the multi-GB `.h5ad`
checkpoints.

LAZY HEAVY-IMPORT DESIGN (important — read before editing):
  This module MUST import cleanly with ONLY the Python standard library + pyyaml. Bare analysis
  boxes have no pandas / numpy / matplotlib / jscatter. Therefore EVERY `import pandas` /
  `import numpy` / `import matplotlib` / `import jscatter` happens LAZILY, inside the function
  that actually needs it — never at module top level. `find_root`, `load_interactive_config`, and
  the config accessors run on stdlib + yaml alone, so they are testable on a bare box. The
  data / plotting / widget functions degrade to a clear ImportError-with-context if their
  backend is absent.

THE OOM LIFECYCLE CONTRACT (the #1 reliability requirement):
  jscatter panels are ipywidgets — every panel a cell builds stays alive in the kernel until
  `.widget.close()` is called. Re-running a grid cell a dozen times while brushing accumulates a
  dozen dead panels and OOMs the kernel. Therefore drive panels ONLY through `grid()` /
  `close_panels()` / `live_panels()`: `grid()` closes any prior panels FIRST, registers the new
  ones, and returns the composed widget. Never hand-roll a raw `_mk` loop in a notebook cell.

Reuse (one import, a few calls):
  from helpers.interactive_style.interactive_helpers import (
      load_interactive_config, load_explorer, grid, first_selection, save_selection, snapshot)
  config = load_interactive_config("02_analysis/config/analysis_config.yaml")
  df = load_explorer("01_qc", config=config)            # compact parquet from 03_results/interactive/
  PANELS = grid(df, ["FOXP3", "pct_counts_mt", "leiden", "tissue"], cat=["leiden", "tissue"])
  sel = first_selection()                               # brushed indices from the live panels
  save_selection(df, sel, label="mt_hi_pocket", subdir="01_qc_explore", config=config)
  snapshot(df, "FOXP3", name="mt_hi_pocket", subdir="01_qc_explore", indices=sel, config=config)
"""
from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Union

import yaml  # stdlib-adjacent (pyyaml); the ONLY non-stdlib top-level import allowed here.

# Default location of the project config, relative to the compartment root.
_DEFAULT_CONFIG_PATH = "02_analysis/config/analysis_config.yaml"

# Fallback FLOORS for the `interactive:` config block, used only when a key is absent. The
# project config is authoritative. `save_selection_cols` deliberately defaults to the two coord
# columns ONLY — project-specific obs/gene columns are opt-in via config, never hardcoded here.
_INTERACTIVE_DEFAULTS: Dict[str, Any] = {
    "save_selection_cols": ["x", "y"],   # columns persisted per selection (extend via config)
    "categorical_obs": [],               # obs columns to ALWAYS treat categorical in grid()
    "summary_stats": [],                 # optional config-driven composition stats (see save_selection)
    "grid_height": 340,                  # per-panel pixel height in grid()
    "grid_rows": 2,                      # rows in the composed panel layout
    "cmap": "magma",                     # continuous colormap for grid() + snapshot()
}


# =============================================================================================
# 0. ROOT + CONFIG — sentinel walk-up + read the `interactive:` block (stdlib + yaml only)
# =============================================================================================
def find_root(sentinel: str = "02_analysis/config/analysis_config.yaml") -> Path:
    """Walk up from CWD to the compartment root (works from any notebook location).

    A live explorer `.qmd` may be opened from `02_analysis/notebooks/<nb>/`; this resolves the
    compartment root by walking up to the sentinel config so `load_explorer`/`save_selection`
    build correct paths regardless of the launch directory.
    """
    d = Path.cwd().resolve()
    while not (d / sentinel).exists():
        if d.parent == d:
            raise FileNotFoundError(
                f"compartment root ({sentinel}) not found above {Path.cwd()}\n"
                "  -> open the explorer from inside the compartment tree.")
        d = d.parent
    return d


def load_interactive_config(path: Union[str, Path, None] = None) -> Dict[str, Any]:
    """Load the project `analysis_config.yaml` and return the full parsed dict.

    Pass the parsed dict to every other function as `config=`. Reading once per notebook and
    threading it keeps these helpers free of global state and import-time side effects. `path`
    defaults to `02_analysis/config/analysis_config.yaml` under the resolved compartment root.
    """
    p = Path(path) if path is not None else (find_root() / _DEFAULT_CONFIG_PATH)
    with open(p, "r") as fh:
        return yaml.safe_load(fh) or {}


def _interactive(config: Optional[Dict[str, Any]]) -> Dict[str, Any]:
    """Return the merged `interactive:` block: project values over `_INTERACTIVE_DEFAULTS` floors."""
    merged = dict(_INTERACTIVE_DEFAULTS)
    if config:
        merged.update(config.get("interactive", {}) or {})
    return merged


def _int_get(config: Optional[Dict[str, Any]], key: str) -> Any:
    """One interactive-key lookup with default fallback."""
    return _interactive(config).get(key, _INTERACTIVE_DEFAULTS.get(key))


def _results_root(config: Optional[Dict[str, Any]]) -> Path:
    """Resolve the results root (`paths.results`, default `03_results/`), anchored at the root."""
    paths = (config or {}).get("paths", {}) or {}
    root = find_root()
    results = paths.get("results", "03_results/")
    p = Path(results)
    return p if p.is_absolute() else root / p


# =============================================================================================
# 1. LOAD — read a compact explorer parquet built by the project's export_explorers.py
# =============================================================================================
def load_explorer(name: str, config: Optional[Dict[str, Any]] = None) -> "pd.DataFrame":
    """Load a compact explorer table, e.g. load_explorer('01_qc').

    Reads `<results_root>/interactive/<name>_explore.parquet` (built by the project's
    `02_analysis/stages/export_explorers.py`). Index = cell barcode. Raises a clear
    FileNotFoundError telling the user to run the exporter first.
    """
    import pandas as pd  # lazy
    path = _results_root(config) / "interactive" / f"{name}_explore.parquet"
    if not path.exists():
        raise FileNotFoundError(
            f"{path} not found — run `python 02_analysis/stages/export_explorers.py` first "
            "to materialize the compact explorer tables.")
    return pd.read_parquet(path)


# =============================================================================================
# 2. PALETTE — a jscatter-friendly categorical color map (glasbey when many levels)
# =============================================================================================
def color_key(df: "pd.DataFrame", col: str, config: Optional[Dict[str, Any]] = None) -> dict:
    """A jscatter categorical palette dict for `col` (glasbey if many levels, else Okabe-Ito).

    Returns {level: hex} suitable for `Scatter.color(by=col, map=color_key(df, col))`. Imports
    jscatter lazily so this module still imports on a box without it.
    """
    import jscatter  # lazy
    levels = sorted(df[col].astype(str).unique())
    pal = list(jscatter.glasbey_dark) if len(levels) > 8 else list(jscatter.okabe_ito)
    return {lvl: pal[i % len(pal)] for i, lvl in enumerate(levels)}


# =============================================================================================
# 3. LIFECYCLE — the OOM-safe panel registry (grid / close_panels / live_panels)
# =============================================================================================
# Module-level registry of the currently-live jscatter panels. jscatter panels are ipywidgets:
# each stays alive in the kernel until `.widget.close()` is called. `grid()` closes the prior
# panels before building new ones, so re-running a grid cell never accumulates dead widgets.
_LIVE_PANELS: List[Any] = []


def _mk(df: "pd.DataFrame", channel: str, *, categorical: bool, height: int,
        cmap: str, config: Optional[Dict[str, Any]],
        tooltip_properties: Optional[Sequence[str]] = None) -> Any:
    """Build ONE jscatter Scatter panel colored by `channel`. Internal — use grid()."""
    import jscatter  # lazy
    s = jscatter.Scatter(data=df, x="x", y="y", height=height)
    # `labeling={'variable': channel}` names the legend. Passing color via a raw `map` leaves
    # jscatter's legend title (`legend_encoding['color']['variable']`) as None, so the legend
    # renders swatches next to bare values with no indication of WHICH channel — unreadable.
    if categorical:
        s.color(by=channel, map=color_key(df, channel, config=config),
                labeling={"variable": channel})
    else:
        lo, hi = float(df[channel].min()), float(df[channel].max())
        s.color(by=channel, map=cmap, norm=[lo, hi if hi > lo else lo + 1e-9],
                labeling={"variable": channel})
    # Legend at top-right (the point cloud's high-mito/low-gene mass sits bottom-left, so a
    # top-left legend would overlap it) and `size="medium"` so it is not an easy-to-miss overlay.
    s.legend(True, position="top-right", size="medium")
    # Hover readout. Without an explicit `.tooltip(True, ...)` jscatter shows NOTHING on hover —
    # `sync_hover` in compose() only propagates WHICH point is hovered, it does not create tooltip
    # content. Default to showing every grid channel for the hovered cell (falls back to this
    # panel's own channel) so a linked hover reads out all panels at once.
    props = list(tooltip_properties) if tooltip_properties else [channel]
    props = [p for p in props if p in df.columns]
    s.tooltip(True, properties=props or [channel])
    return s


def grid(df: "pd.DataFrame", channels: Sequence[str], *, cat: Optional[Sequence[str]] = None,
         rows: Optional[int] = None, height: Optional[int] = None,
         config: Optional[Dict[str, Any]] = None, **compose_kw) -> Any:
    """Build a linked grid of jscatter panels — the OOM-safe entry point.

    Closes any previously-live panels FIRST (so re-running the cell never leaks widgets),
    builds one panel per channel (categorical via `color_key`, continuous via a magma norm),
    registers them in `_LIVE_PANELS`, and returns the composed widget with view/selection/hover
    synced across panels. Brush or lasso ANY panel — the rest highlight the same cells.

    Args:
      channels:   column names to plot, one panel each (2-4 is the sweet spot).
      cat:        channels to force categorical; others are inferred by dtype. Config
                  `interactive.categorical_obs` columns are ALWAYS treated categorical.
      rows/height: layout overrides (default from `interactive.grid_rows` / `grid_height`).
      compose_kw:  extra kwargs forwarded to jscatter.compose (e.g. row_height).

    ALWAYS drive panels through this function — never a raw `_mk` loop in a cell — so the
    lifecycle contract holds. See `close_panels()` for the manual escape hatch.
    """
    import jscatter  # lazy
    close_panels()  # OOM hygiene: drop any panels a prior run of this cell left alive

    always_cat = set(_int_get(config, "categorical_obs") or [])
    forced_cat = set(cat or []) | always_cat
    rows = int(rows if rows is not None else _int_get(config, "grid_rows"))
    height = int(height if height is not None else _int_get(config, "grid_height"))
    cmap = str(_int_get(config, "cmap"))

    import pandas as pd  # lazy — only for the dtype inference
    panels: List[Any] = []
    for c in channels:
        categorical = c in forced_cat or not pd.api.types.is_numeric_dtype(df[c])
        panels.append(_mk(df, c, categorical=categorical, height=height, cmap=cmap,
                          config=config, tooltip_properties=channels))
    _LIVE_PANELS.extend(panels)

    return jscatter.compose(
        [(s, c) for s, c in zip(panels, channels)],
        sync_view=True, sync_selection=True, sync_hover=True, rows=rows, **compose_kw)


def close_panels() -> int:
    """Close every live jscatter panel's widget and clear the registry. Returns the count closed.

    The manual escape hatch for the OOM lifecycle contract — call it from a standalone cell if a
    grid ever lags, or before switching notebooks. `grid()` calls this for you on every re-run.
    """
    n = 0
    for s in _LIVE_PANELS:
        try:
            s.widget.close()
            n += 1
        except Exception:  # already closed / no widget — nothing to free
            pass
    _LIVE_PANELS.clear()
    return n


def live_panels() -> List[Any]:
    """Return the list of currently-live jscatter panels (for first_selection / manual close)."""
    return list(_LIVE_PANELS)


# =============================================================================================
# 4. SELECTION — pull the brushed indices out of the live panels into pandas
# =============================================================================================
def first_selection(scatters: Optional[Sequence[Any]] = None) -> "np.ndarray":
    """Return the brushed row indices from the FIRST panel that has a non-empty selection.

    With `scatters=None` (the common case) it reads the live panels registered by `grid()`, so
    `df.iloc[first_selection()]` works no matter which panel you lassoed. Returns an empty int
    array if nothing is selected.
    """
    import numpy as np  # lazy
    for sc in (scatters if scatters is not None else live_panels()):
        try:
            s = np.asarray(sc.selection())
        except Exception:
            continue
        if s.size:
            return s.astype(int)
    return np.array([], dtype=int)


# =============================================================================================
# 5. PERSIST — durable barcode CSV + a labelled matplotlib snapshot (EDA-tier, under the notebook)
# =============================================================================================
def _eda_dir(subdir: str) -> Path:
    """`02_analysis/notebooks/<subdir>/eda/` — the durable landing spot for a selection."""
    d = find_root() / "02_analysis" / "notebooks" / subdir / "eda"
    d.mkdir(parents=True, exist_ok=True)
    return d


def save_selection(df: "pd.DataFrame", indices, label: str, subdir: str,
                   config: Optional[Dict[str, Any]] = None,
                   cols: Optional[Sequence[str]] = None) -> Path:
    """Persist a selection's BARCODES + key columns for durable reference.

    Writes `02_analysis/notebooks/<subdir>/eda/selection_<label>.csv` (index = cell barcode) and
    appends a one-row summary to `.../eda/selections_index.csv` for quick cross-selection
    comparison. `indices` is what `first_selection()` / `scatter.selection()` returns.

    `cols` defaults to `interactive.save_selection_cols` (DEFAULT ['x', 'y'] — extend via config,
    never hardcoded here). The summary row carries `n_cells` + the mean of each numeric column in
    `cols`, plus any config-driven composition stats from `interactive.summary_stats` — a list of
    {name, column, value?} entries: with `value`, the stat is the fraction of the selection where
    `column == value`; without, it is the modal (most-common) level of `column` and its fraction.
    """
    import numpy as np  # lazy
    import pandas as pd  # lazy
    idx = np.asarray(indices, dtype=int)
    if idx.size == 0:
        raise ValueError("empty selection — brush or lasso a cluster in the grid first")

    cols = [c for c in (cols or _int_get(config, "save_selection_cols")) if c in df.columns]
    sub = df.iloc[idx]
    out = _eda_dir(subdir) / f"selection_{label}.csv"
    sub[cols].to_csv(out)  # index (barcode) is written

    # Summary row: n_cells + numeric-column means + config-driven composition stats.
    summ: Dict[str, Any] = {"label": label, "n_cells": len(sub)}
    for c in cols:
        if c in sub and pd.api.types.is_numeric_dtype(sub[c]):
            summ[c] = round(float(sub[c].mean()), 4)
    for spec in (_int_get(config, "summary_stats") or []):
        name, col = spec.get("name"), spec.get("column")
        if not name or col not in sub:
            continue
        if "value" in spec:                                   # fraction where column == value
            summ[name] = round(float((sub[col].astype(str) == str(spec["value"])).mean()), 3)
        else:                                                 # modal level + its fraction
            vc = sub[col].value_counts(normalize=True)
            if len(vc):
                summ[name] = vc.idxmax()
                summ[f"{name}_frac"] = round(float(vc.max()), 3)

    index_csv = _eda_dir(subdir) / "selections_index.csv"
    prev = pd.read_csv(index_csv) if index_csv.exists() else pd.DataFrame()
    if "label" in prev:
        prev = prev[prev["label"] != label]
    pd.concat([prev, pd.DataFrame([summ])], ignore_index=True).to_csv(index_csv, index=False)
    print(f"saved {len(sub)} barcodes -> {out.relative_to(find_root())}")
    return out


def snapshot(df: "pd.DataFrame", color: str, name: str, subdir: str,
             config: Optional[Dict[str, Any]] = None, indices=None,
             cmap: Optional[str] = None, title: Optional[str] = None):
    """Render a LABELLED matplotlib snapshot of the embedding to `.../<subdir>/eda/<name>.{png,pdf}`.

    jscatter's own viewport export grabs only the points — no legend, axes, or title. This exists
    to capture those. Colors by `color`; if `indices` is given, non-selected cells are greyed and
    the selection drawn on top. EDA-tier (not a 03_results deliverable), so if figure-style's
    `save_figure` is importable it is used, else it saves directly. LAZY matplotlib import.

    The dense point cloud is rasterized (via figure-style's `rasterize_axes`) before saving, so
    the PDF embeds the dots as a crisp raster at `figures.rasterized_dpi` while the title / axes /
    legend stay vector — keeping the file small and fast to open even at 100k+ cells, the same
    contract `save_figure` applies to 03_results figures.
    """
    import matplotlib  # lazy
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd

    cmap = cmap or str(_int_get(config, "cmap"))
    categorical = not pd.api.types.is_numeric_dtype(df[color])
    fig, ax = plt.subplots(figsize=(8, 8))
    xy = df[["x", "y"]].to_numpy()

    if indices is not None and np.asarray(indices).size:
        mask = np.zeros(len(df), bool)
        mask[np.asarray(indices, int)] = True
        ax.scatter(xy[~mask, 0], xy[~mask, 1], s=2, c="lightgrey", linewidths=0)
        base = df.iloc[np.asarray(indices, int)]
        pts = xy[mask]
    else:
        base = df
        pts = xy

    if categorical:
        pal = color_key(df, color, config=config)
        for lvl, hexcol in pal.items():
            m = (base[color].astype(str).to_numpy() == lvl)
            ax.scatter(pts[m, 0], pts[m, 1], s=3, c=hexcol, label=str(lvl), linewidths=0)
        ax.legend(markerscale=4, fontsize=8, loc="best", frameon=True, title=color)
    else:
        sctr = ax.scatter(pts[:, 0], pts[:, 1], s=3, c=base[color].to_numpy(), cmap=cmap,
                          linewidths=0)
        fig.colorbar(sctr, ax=ax, shrink=0.7, label=color)

    n_note = f"  (n={len(base)} selected)" if indices is not None else ""
    ax.set_title(title or f"{color}{n_note}")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel("dim-1")
    ax.set_ylabel("dim-2")
    fig.tight_layout()

    # Rasterize the dense point cloud BEFORE saving: the embedding is thousands-to-millions of
    # scatter dots, which as pure vector make the PDF enormous and slow to open. Marking the dot
    # collections rasterized embeds them as a crisp raster (at figures.rasterized_dpi) in the PDF
    # while title/axes/legend stay vector — the same contract save_figure applies to 03_results
    # figures. Routed through figure-style's rasterize_axes so the idea lives in ONE place.
    _rasterize_dense(ax)

    d = _eda_dir(subdir)
    saved = _save_via_figure_style(fig, subdir, name, config)
    if not saved:
        raster_dpi = _raster_dpi(config)
        fig.savefig(d / f"{name}.png", dpi=200, bbox_inches="tight")
        fig.savefig(d / f"{name}.pdf", dpi=raster_dpi, bbox_inches="tight")
    print(f"saved snapshot -> {(d / f'{name}.png').relative_to(find_root())}")
    return fig


def _raster_dpi(config: Optional[Dict[str, Any]]) -> float:
    """DPI at which rasterized layers embed into the PDF (`figures.rasterized_dpi`, default 600).

    Mirrors save_figure's PDF behaviour so an explorer snapshot and a 03_results figure embed
    their dense dot layer at the same crisp-but-cheap resolution.
    """
    figures = (config or {}).get("figures", {}) or {}
    try:
        return float(figures.get("rasterized_dpi", 600))
    except (TypeError, ValueError):
        return 600.0


def _rasterize_dense(*axes: Any) -> None:
    """Mark each Axes' dense scatter collections rasterized (prefer figure-style's rasterize_axes).

    Owns the fallback so the behaviour holds even on a box where figure-style is not importable:
    when the shared helper is present the idea lives in ONE place (figure_helpers.rasterize_axes),
    otherwise we set `rasterized=True` on the collections directly. Best-effort — never fatal.
    """
    try:
        from helpers.figure_style import rasterize_axes  # per-project shim
    except Exception:
        try:
            from figure_helpers import rasterize_axes  # symlinked lib, no shim
        except Exception:
            rasterize_axes = None
    if rasterize_axes is not None:
        try:
            rasterize_axes(*axes)
            return
        except Exception:
            pass
    for ax in axes:
        if ax is None:
            continue
        for coll in getattr(ax, "collections", []):
            coll.set_rasterized(True)


def _save_via_figure_style(fig, subdir: str, name: str,
                           config: Optional[Dict[str, Any]]) -> bool:
    """Save through figure-style's save_figure into the notebook's eda/ dir if it is importable.

    Returns True if it wrote the files. figure-style resolves output dirs under 03_results/, so
    we save directly into the notebook eda/ dir here (EDA-tier) and only borrow its styling
    ability when present; a plain savefig is the fallback. Best-effort — never fatal.
    """
    try:
        from helpers.figure_style import set_paper_style  # per-project shim
    except Exception:
        try:
            from figure_helpers import set_paper_style  # symlinked lib, no shim
        except Exception:
            return False
    try:
        set_paper_style(config=config)  # apply the legible tier to the already-built figure
        d = _eda_dir(subdir)
        # PNG at screen dpi; PDF at figures.rasterized_dpi so any rasterized dot layer embeds
        # crisp (text/axes stay vector) — matches save_figure's dual-variant contract.
        fig.savefig(d / f"{name}.png", dpi=200, bbox_inches="tight")
        fig.savefig(d / f"{name}.pdf", dpi=_raster_dpi(config), bbox_inches="tight")
        return True
    except Exception:
        return False
