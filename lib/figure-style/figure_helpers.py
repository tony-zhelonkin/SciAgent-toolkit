"""
figure_helpers.py — the SciAgent-toolkit cross-language FIGURE-STYLE CONTRACT (Python side).
============================================================================================
ONE place that owns the project's figure format so phase viz scripts never reinvent it. This
is the Python half of a two-language contract; the R half (`figure_helpers.R`, same directory)
has FUNCTION PARITY — identical public names, equivalent semantics — and both read the SAME
`analysis_config.yaml:figures` block. Centralizing styling here is the load-bearing capability
behind the owner's #1 recurring pain point (figure legibility) and #2 (results placement).

Config keys read (from `analysis_config.yaml:figures`):
  base_size, title_size, axis_title_size, axis_text_size, strip_size, legend_text_size,
  label_size, line_width, point_size, base_size_column, width, height, width_column,
  height_column, dpi, top_n, volcano_label_top, z_clamp, nes_cap, caption_wrap_column,
  variants, by_contrast_dir, overview_dir
Plus, from elsewhere in the config: `paths.results` (results root), `paths.master`
(master-table root), `paths.stage_tables_subdir` / `paths.stage_figures_subdir`.

LAZY HEAVY-IMPORT DESIGN (important — read before editing):
  This module MUST import cleanly with ONLY the Python standard library + pyyaml. Bare analysis
  boxes have no matplotlib / pandas / numpy. Therefore EVERY `import matplotlib` / `import pandas`
  / `import numpy` happens LAZILY, inside the function that actually needs it — never at module
  top level. The path-building, config-reading, caption-writing, master-table, rounding, and
  direction-cue functions all run on stdlib + yaml alone, so they are testable on a bare box.
  Plotting functions (set_paper_style / save_figure / save_overview / style_series) degrade to a
  clear ImportError-with-context if their backend (matplotlib) is absent.

Reuse (one import, a few calls):
  from figure_helpers import set_paper_style, save_overview, overview_path, write_caption
  config = load_figure_config("02_analysis/config/analysis_config.yaml")
  set_paper_style(config=config)                          # once, top of the viz script
  fig = ...                                               # build ONE matplotlib Figure
  save_overview(fig, "04_gsea", "gsea_hallmark_heatmap", table=rows,
                finding="Hallmark IFN-alpha/gamma dominate the ISD90 response.",
                script="02_analysis/scripts/11_gsea_viz.py", fn="save_overview",
                config_kv="figures.nes_cap = 3.5", input="03_results/objects/gsea.rds",
                how_to_read="Rows = pathways; color = NES (orange up / blue down); * = padj<0.05.",
                config=config)                             # figure + sibling table + caption, atomic
"""
from __future__ import annotations

import csv
import math
import os
import textwrap
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Union

import yaml  # stdlib-adjacent (pyyaml); the ONLY non-stdlib top-level import allowed here.

# ---------------------------------------------------------------------------------------------
# Per-variant defaults mirroring the shipped analysis_config.yaml:figures block. These are the
# fallback FLOORS used only when a key is absent from the project config; the project config is
# authoritative. Keep in sync with the R file's `.FIG_DEFAULTS` and the config template.
# ---------------------------------------------------------------------------------------------
_FIG_DEFAULTS: Dict[str, Any] = {
    "base_size": 16,
    "title_size": 18,
    "axis_title_size": 15,
    "axis_text_size": 13,
    "strip_size": 14,
    "legend_text_size": 13,
    "label_size": 5,
    "line_width": 0.8,
    "point_size": 2.0,
    "base_size_column": 9,
    "width": 10,
    "height": 8,
    "width_column": 3.5,
    "height_column": 3.0,
    "dpi": 300,
    "top_n": 20,
    "volcano_label_top": 10,
    "z_clamp": 2.5,
    "nes_cap": 3.5,
    "caption_wrap_column": 70,
    "variants": ["print", "screen"],
    "by_contrast_dir": "by_contrast",
    "overview_dir": "_overview",
}

# Default location of the project config, relative to the project root.
_DEFAULT_CONFIG_PATH = "02_analysis/config/analysis_config.yaml"


# =============================================================================================
# 0. CONFIG — read the figures block once; everything else takes `config=` (no naked literals)
# =============================================================================================
def load_figure_config(path: Union[str, Path, None] = None) -> Dict[str, Any]:
    """Load the project `analysis_config.yaml` and return the full parsed dict.

    Pass the parsed dict to every other function as `config=`. Reading once per script and
    threading it keeps these helpers free of global state and import-time side effects.
    `path` defaults to `02_analysis/config/analysis_config.yaml` under the cwd.
    """
    p = Path(path or _DEFAULT_CONFIG_PATH)
    with open(p, "r") as fh:
        cfg = yaml.safe_load(fh) or {}
    return cfg


def _figures(config: Optional[Dict[str, Any]]) -> Dict[str, Any]:
    """Return the merged `figures:` block: project values over `_FIG_DEFAULTS` floors."""
    fig = {}
    fig.update(_FIG_DEFAULTS)
    if config:
        fig.update(config.get("figures", {}) or {})
    return fig


def _fig_get(config: Optional[Dict[str, Any]], key: str) -> Any:
    """One figures-key lookup with default fallback (mirrors R's `%||%` on CONFIG$figures)."""
    return _figures(config).get(key, _FIG_DEFAULTS.get(key))


def _results_root(config: Optional[Dict[str, Any]]) -> Path:
    """Resolve the results root (`paths.results`, default `03_results/`)."""
    paths = (config or {}).get("paths", {}) or {}
    return Path(paths.get("results", "03_results/"))


def _master_root(config: Optional[Dict[str, Any]]) -> Path:
    """Resolve the master-table root (`paths.master`, default `<results>/master/`)."""
    paths = (config or {}).get("paths", {}) or {}
    return Path(paths.get("master", str(_results_root(config) / "master")))


def _stage_dir(config: Optional[Dict[str, Any]], stage: str, kind: str) -> Path:
    """`03_results/<stage>/<figures|tables>/` — the one place the {figures,tables} subdir is named.

    `kind` ∈ {"figures", "tables"}. Subdir names come from `paths.stage_{figures,tables}_subdir`.
    """
    if kind not in ("figures", "tables"):
        raise ValueError(f"kind must be 'figures' or 'tables', got {kind!r}")
    paths = (config or {}).get("paths", {}) or {}
    subdir_key = "stage_figures_subdir" if kind == "figures" else "stage_tables_subdir"
    subdir = paths.get(subdir_key, kind)
    return _results_root(config) / stage / subdir


# =============================================================================================
# 1. PATHS — per-contrast and cross-contrast; the ONLY sanctioned way to build these dirs
# =============================================================================================
def contrast_path(stage: str, contrast: str, kind: str = "figures",
                  config: Optional[Dict[str, Any]] = None) -> Path:
    """Build + mkdir `03_results/<stage>/<kind>/by_contrast/<contrast>/`; return the dir.

    `kind` ∈ {"figures", "tables"}. The `by_contrast` subdir name comes from
    `figures.by_contrast_dir`. Contrast dir name MUST be the exact config contrast name.
    """
    d = _stage_dir(config, stage, kind) / str(_fig_get(config, "by_contrast_dir")) / contrast
    d.mkdir(parents=True, exist_ok=True)
    return d


def overview_path(stage: str, kind: str = "figures",
                  config: Optional[Dict[str, Any]] = None) -> Path:
    """Build + mkdir `03_results/<stage>/<kind>/_overview/` (cross-contrast); return the dir.

    `kind` ∈ {"figures", "tables"}. The `_overview` subdir name comes from `figures.overview_dir`.
    """
    d = _stage_dir(config, stage, kind) / str(_fig_get(config, "overview_dir"))
    d.mkdir(parents=True, exist_ok=True)
    return d


def _resolve_fig_dir(stage: str, contrast: Optional[str], overview: bool,
                     config: Optional[Dict[str, Any]]) -> Path:
    """Pick the figures dir: by_contrast/<c>/ if `contrast` given, else _overview/ if `overview`."""
    if contrast is not None:
        return contrast_path(stage, contrast, "figures", config)
    if overview:
        return overview_path(stage, "figures", config)
    # Plain stage figures dir (no sub-layout) — e.g. QC/EDA panels that are neither.
    d = _stage_dir(config, stage, "figures")
    d.mkdir(parents=True, exist_ok=True)
    return d


# =============================================================================================
# 2. THEME — the SINGLE style entry point. set_paper_style() is the Python canonical name;
#    project_theme() is the thin same-named alias so a cross-language parity grep finds BOTH
#    contract names in this file (and the R file likewise defines both).
# =============================================================================================
def set_paper_style(base_size: Optional[float] = None, variant: str = "screen",
                    config: Optional[Dict[str, Any]] = None) -> None:
    """Apply matplotlib rcParams from the `figures:` config (LAZY matplotlib import).

    Bold axis titles, no top/right spines, Illustrator-editable text (`pdf.fonttype=42`). Sizes
    come from the config's per-variant font tier (`base_size` for screen, `base_size_column` for
    print) and are enforced as FLOORS (clamped up + warned if a passed `base_size` is below).
    Call ONCE near the top of a viz script. Raises a clear ImportError if matplotlib is absent.

    `variant` ∈ {"screen", "print"} selects which font-tier floor governs.
    """
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError as exc:  # pragma: no cover - exercised only on a backend-present box
        raise ImportError(
            "set_paper_style() needs matplotlib (the plotting backend). Install matplotlib, or "
            "call only the path/caption/table helpers (which need no backend)."
        ) from exc

    floor = float(_variant_base_floor(variant, config))
    bs = _enforce_floor(base_size if base_size is not None else floor, floor, "base_size", variant)
    f = _figures(config)
    plt.rcParams.update({
        "figure.dpi": 100,                                  # on-screen; save dpi set in save_figure
        "savefig.dpi": float(f["dpi"]),
        "font.size": bs,
        "axes.titlesize": float(f["title_size"]),
        "axes.titleweight": "bold",
        "axes.labelsize": float(f["axis_title_size"]),
        "axes.labelweight": "bold",                          # bold axis titles (contract)
        "xtick.labelsize": float(f["axis_text_size"]),
        "ytick.labelsize": float(f["axis_text_size"]),
        "legend.fontsize": float(f["legend_text_size"]),
        "legend.title_fontsize": float(f["legend_text_size"]),
        "lines.linewidth": float(f["line_width"]),
        "lines.markersize": float(f["point_size"]),
        "axes.spines.top": False,                            # remove top spine (declutter)
        "axes.spines.right": False,                          # remove right spine (declutter)
        "figure.autolayout": False,
        "pdf.fonttype": 42,                                  # editable text in Illustrator
        "ps.fonttype": 42,
    })


def project_theme(base_size: Optional[float] = None, legend: bool = True,
                  variant: str = "screen", config: Optional[Dict[str, Any]] = None) -> None:
    """Cross-language alias of set_paper_style() (the R-side canonical name).

    Present so a parity grep finds `project_theme` in BOTH files. `legend` is accepted for R
    signature parity (matplotlib legends are placed per-axes, so it is a no-op here).
    """
    set_paper_style(base_size=base_size, variant=variant, config=config)


# =============================================================================================
# 2b. Per-variant font floors — the legibility contract (clamp up, warn; never silently shrink)
# =============================================================================================
def _variant_base_floor(variant: str, config: Optional[Dict[str, Any]]) -> float:
    """The base-font FLOOR for a variant: base_size (screen) / base_size_column (print)."""
    if variant == "print":
        return float(_fig_get(config, "base_size_column"))
    return float(_fig_get(config, "base_size"))


def _enforce_floor(value: float, floor: float, what: str, variant: str) -> float:
    """Clamp `value` up to `floor` if below it, printing a warning. Legibility never silently lost."""
    value = float(value)
    if value < floor:
        print(f"  [figure-style] WARNING — {what}={value:g} below {variant} floor {floor:g}; "
              f"clamping up to {floor:g} (legibility contract).")
        return floor
    return value


def _variant_geometry(variant: str, config: Optional[Dict[str, Any]]) -> "tuple[float, float]":
    """(width, height) inches for a variant: width/height (screen) / *_column (print)."""
    if variant == "print":
        return float(_fig_get(config, "width_column")), float(_fig_get(config, "height_column"))
    return float(_fig_get(config, "width")), float(_fig_get(config, "height"))


# =============================================================================================
# 3. EXPORT — one plot object, two variant artifacts, from ONE config-driven code path
# =============================================================================================
def save_figure(plot: Any, stage: str, name: str, variant: str = "both",
                contrast: Optional[str] = None, overview: bool = False,
                config: Optional[Dict[str, Any]] = None) -> Dict[str, Path]:
    """Render ONE matplotlib Figure to dual variants from one call (LAZY matplotlib import).

    Variants (default "both"):
      print  → `<stem>.print.pdf`  — vector PDF, column geometry (width_column x height_column),
               print font tier (base_size_column), `pdf.fonttype=42` so text stays Illustrator-
               editable; any `rasterized=True` dense layers stay raster while text/axes vectorize.
      screen → `<stem>.screen.png` — raster PNG @ dpi, screen geometry (width x height),
               screen font tier (base_size).
    `variant` ∈ {"print", "screen", "both"}. The output dir is resolved via contrast_path()
    (if `contrast` given) / overview_path() (if `overview`) / the plain stage figures dir.
    Stale `<name>*.{png,pdf}` are purged first so the run owns its namespace. The same plot
    object is re-sized + re-styled per variant, so one figure yields both print and screen files.

    Returns {variant: Path} of the files written.
    """
    try:
        import matplotlib.pyplot as plt  # noqa: F401  (validates backend presence)
    except ImportError as exc:  # pragma: no cover - exercised only on a backend-present box
        raise ImportError(
            "save_figure() needs matplotlib to render. On a backend-less box call the path / "
            "caption / table helpers instead, which need no plotting backend."
        ) from exc

    variants = _normalize_variants(variant, config)
    out_dir = _resolve_fig_dir(stage, contrast, overview, config)
    purge_figures(stage, name, contrast=contrast, overview=overview, config=config)

    written: Dict[str, Path] = {}
    f = _figures(config)
    for v in variants:
        ext = "pdf" if v == "print" else "png"
        out = out_dir / f"{name}.{v}.{ext}"
        w, h = _variant_geometry(v, config)
        floor = _variant_base_floor(v, config)
        # Re-apply the per-variant font tier + geometry to the SAME figure object.
        try:
            plot.set_size_inches(w, h)
        except AttributeError:
            pass  # not a Figure with set_size_inches (e.g. a seaborn FacetGrid) — size via savefig
        _apply_variant_fonts(plot, floor, v)
        dpi = float(f["dpi"]) if v == "screen" else float(f.get("dpi", 300))
        plot.savefig(out, dpi=dpi, bbox_inches="tight")
        written[v] = out
        print(f"  [figure-style] save_figure: {out.name} ({w:g}x{h:g}in, base>= {floor:g}pt)")
    return written


def _normalize_variants(variant: str, config: Optional[Dict[str, Any]]) -> List[str]:
    """Expand the `variant` arg to a concrete list; "both" => the config's `figures.variants`."""
    if variant == "both":
        vs = list(_fig_get(config, "variants") or ["print", "screen"])
    else:
        vs = [variant]
    for v in vs:
        if v not in ("print", "screen"):
            raise ValueError(f"variant must be print/screen/both, got {v!r}")
    return vs


def _apply_variant_fonts(fig: Any, floor: float, variant: str) -> None:
    """Bump every text artist on the figure up to the per-variant floor (legibility contract)."""
    try:
        for txt in fig.findobj(match=lambda o: hasattr(o, "get_fontsize")):
            try:
                if txt.get_fontsize() < floor:
                    txt.set_fontsize(floor)
            except (TypeError, ValueError):
                continue
    except Exception:  # pragma: no cover - defensive; never let font-bump break a save
        pass


# =============================================================================================
# 3b. SERIES POST-STYLER — fix axis/legend for cross-panel comparability (port of style_running_sum)
# =============================================================================================
def style_series(plot: Any, ylim: Optional[Sequence[float]] = None,
                 config: Optional[Dict[str, Any]] = None) -> Any:
    """Pin a shared y-range + a fixed inside legend so a SERIES of figures stays comparable.

    Ported from 14839's `style_running_sum`: across a family of figures (e.g. one running-sum
    per database) a name-length-dependent outside legend silently resizes the plotting panel,
    so the curves stop being comparable. This clamps the y-axis to a fixed range (zoom, never
    drops data) and moves the legend INSIDE (zero layout width) so every figure in the series
    has identical panel proportions. `ylim` defaults to a symmetric clamp from config when given;
    otherwise the existing axis limits are kept. LAZY matplotlib usage. Returns `plot`.
    """
    try:
        import matplotlib.pyplot as plt  # noqa: F401
    except ImportError as exc:  # pragma: no cover
        raise ImportError("style_series() needs matplotlib.") from exc

    if ylim is None:
        z = _fig_get(config, "z_clamp")
        if z is not None:
            ylim = (-float(z), float(z))
    axes = getattr(plot, "axes", None) or []
    for ax in axes:
        if ylim is not None:
            ax.set_ylim(float(ylim[0]), float(ylim[1]))
        leg = ax.get_legend()
        if leg is not None:
            # Inside, top-right, zero layout width — identical across the series.
            ax.legend(loc="upper right", frameon=True, framealpha=0.9)
    return plot


# =============================================================================================
# 4. PURGE — delete stale figures before a fresh write so a run OWNS its figure namespace
# =============================================================================================
def purge_figures(stage: str, prefix: str, contrast: Optional[str] = None,
                  overview: bool = False, config: Optional[Dict[str, Any]] = None) -> int:
    """Delete `<prefix>*.{png,pdf}` under the resolved stage figures dir (DC semantics).

    Removes orphaned stems a fresh run no longer produces (e.g. a dropped gene panel). Scoped by
    `prefix` so scripts sharing a figures/ dir don't clobber each other. Needs NO plotting
    backend (stdlib only) — safe to call on a bare box. Returns the count removed.
    """
    # Resolve the dir WITHOUT creating it if it already does not exist (avoid empty dir spam).
    if contrast is not None:
        d = _stage_dir(config, stage, "figures") / str(_fig_get(config, "by_contrast_dir")) / contrast
    elif overview:
        d = _stage_dir(config, stage, "figures") / str(_fig_get(config, "overview_dir"))
    else:
        d = _stage_dir(config, stage, "figures")
    n = 0
    if d.is_dir():
        for ext in ("png", "pdf"):
            for fp in d.glob(f"{prefix}*.{ext}"):
                fp.unlink()
                n += 1
    if n:
        print(f"  [figure-style] purge_figures('{prefix}*'): removed {n} stale file(s)")
    return n


# =============================================================================================
# 5. CAPTION — idempotent create/UPDATE of the sibling stage README.md (no plotting backend)
# =============================================================================================
def write_caption(stage: str, filename: str, finding: str, script: str, fn: str,
                  config_kv: str, input: str, how_to_read: str,
                  config: Optional[Dict[str, Any]] = None) -> Path:
    """Idempotently write/replace ONE artifact's caption section in `03_results/<stage>/README.md`.

    The section is a path-qualified heading (`## figures/_overview/<file>` etc.), a one-sentence
    scientific FINDING, a mandatory **How to read** subsection (glyphs / sign convention / claim
    tier), and a `Script | Function | Config | Input` table. Re-running with the SAME `filename`
    REPLACES that file's section in place (idempotent — never duplicates). Prose is wrapped to
    `figures.caption_wrap_column`. Needs NO plotting backend (stdlib only).

    `filename` should be the path-qualified artifact name relative to the stage dir, e.g.
    "figures/_overview/gsea_hallmark_heatmap.screen.png".
    """
    readme = _results_root(config) / stage / "README.md"
    readme.parent.mkdir(parents=True, exist_ok=True)

    wrap = int(_fig_get(config, "caption_wrap_column"))
    heading = f"## {filename}"
    section = _render_caption_section(heading, finding, script, fn, config_kv, input,
                                      how_to_read, wrap)

    if readme.exists():
        existing = readme.read_text()
        new_text = _replace_section(existing, heading, section)
    else:
        new_text = f"# {stage} — artifact captions\n\n{section}"
    readme.write_text(new_text)
    return readme


def _wrap(text: str, width: int) -> str:
    """Wrap one paragraph to `width`, preserving an empty string for falsy input."""
    if not text:
        return ""
    return "\n".join(textwrap.wrap(str(text), width=width)) or str(text)


def _render_caption_section(heading: str, finding: str, script: str, fn: str, config_kv: str,
                            input: str, how_to_read: str, wrap: int) -> str:
    """Render the canonical caption section block for one artifact."""
    lines = [
        heading,
        "",
        _wrap(finding, wrap),
        "",
        "**How to read:** " + _wrap(how_to_read, wrap),
        "",
        "| Script | Function | Config | Input |",
        "|---|---|---|---|",
        f"| `{script}` | `{fn}` | `{config_kv}` | `{input}` |",
        "",
    ]
    return "\n".join(lines)


def _replace_section(text: str, heading: str, new_section: str) -> str:
    """Replace the `## <heading>` section in `text` with `new_section`; append if absent.

    A section runs from its `## ` heading line to the next `## ` (or `# `) heading or EOF. Match
    is exact on the heading line so distinct path-qualified headings never collide.
    """
    lines = text.splitlines()
    out: List[str] = []
    i = 0
    replaced = False
    n = len(lines)
    while i < n:
        if lines[i].strip() == heading.strip():
            # Skip the existing section body up to the next heading.
            i += 1
            while i < n and not (lines[i].startswith("## ") or lines[i].startswith("# ")):
                i += 1
            # Emit the replacement (already terminated with a blank line).
            out.extend(new_section.rstrip("\n").splitlines())
            out.append("")
            replaced = True
            continue
        out.append(lines[i])
        i += 1
    result = "\n".join(out).rstrip("\n") + "\n"
    if not replaced:
        # Append a fresh section, ensuring one blank-line separator.
        result = result.rstrip("\n") + "\n\n" + new_section.rstrip("\n") + "\n"
    return result


# =============================================================================================
# 6. OVERVIEW — the ATOMIC adjacency mechanism: figure + sibling table + caption in ONE call
# =============================================================================================
def save_overview(plot: Any, stage: str, name: str, table: Any, finding: str, script: str,
                  fn: str, config_kv: str, input: str, how_to_read: str,
                  contrast: Optional[str] = None,
                  config: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
    """Write a figure AND its same-stem source table AND its README caption in one call.

    This is the only sanctioned path for an overview/by-contrast figure: you cannot make the
    figure without its neighbor table + caption (the source-table-adjacency + README-adjacency
    contracts, enforced mechanically). Writes:
      figures/_overview/<name>.<variant>.<ext>   (via save_figure; or by_contrast/<c>/ if contrast)
      tables/_overview/<name>.csv                 (the data behind the figure; round_numeric_cols)
      03_results/<stage>/README.md caption        (via write_caption, path-qualified, idempotent)

    Returns {"figures": {variant: Path}, "table": Path, "readme": Path}.
    """
    overview = contrast is None
    figs = save_figure(plot, stage, name, variant="both", contrast=contrast,
                       overview=overview, config=config)

    # Same-stem source table, adjacent under tables/<sub-layout>/.
    if contrast is not None:
        tdir = contrast_path(stage, contrast, "tables", config)
        rel_sub = f"tables/{_fig_get(config, 'by_contrast_dir')}/{contrast}"
    else:
        tdir = overview_path(stage, "tables", config)
        rel_sub = f"tables/{_fig_get(config, 'overview_dir')}"
    table_path = tdir / f"{name}.csv"
    if table is not None:
        _write_table_csv(round_numeric_cols(table), table_path)

    # Path-qualified caption keyed on the SCREEN artifact (the representative deliverable).
    if contrast is not None:
        fig_rel = f"figures/{_fig_get(config, 'by_contrast_dir')}/{contrast}/{name}.screen.png"
    else:
        fig_rel = f"figures/{_fig_get(config, 'overview_dir')}/{name}.screen.png"
    readme = write_caption(stage, fig_rel, finding=finding, script=script, fn=fn,
                           config_kv=config_kv, input=input, how_to_read=how_to_read,
                           config=config)
    print(f"  [figure-style] save_overview: figure + {rel_sub}/{name}.csv + README caption")
    return {"figures": figs, "table": table_path, "readme": readme}


# =============================================================================================
# 7. MASTER TABLE — idempotent, byte-stable cross-stage accumulator append (no plotting backend)
# =============================================================================================
def append_master_table(df_or_rows: Any, database: str, stage: str, name: str,
                        config: Optional[Dict[str, Any]] = None) -> Path:
    """Idempotently append rows to `03_results/master/<name>.csv`, deduped on the `database` column.

    Re-running for the SAME `database` REPLACES those rows (filter-out then append), so the master
    table is a stable accumulator across stages. `round_numeric_cols(sig=9)` is applied for
    byte-stability (re-runs produce identical files). Accepts a list-of-dicts (stdlib csv) OR a
    pandas DataFrame (lazy pandas) so it is testable on a bare box. Needs NO plotting backend.

    `database` labels the rows this call owns; every incoming row gets/keeps that label in the
    `database` column. `stage` is recorded in a `stage` column for provenance.
    """
    rows = _to_rows(df_or_rows)
    # Stamp provenance columns so the master table is self-describing.
    for r in rows:
        r.setdefault("database", database)
        r["database"] = r.get("database", database) or database
        r.setdefault("stage", stage)
    rows = round_numeric_cols(rows)

    out = _master_root(config) / f"{name}.csv"
    out.parent.mkdir(parents=True, exist_ok=True)

    # Read existing rows, drop any belonging to this `database`, then append the fresh ones.
    existing: List[Dict[str, Any]] = []
    if out.exists():
        with open(out, "r", newline="") as fh:
            existing = [dict(r) for r in csv.DictReader(fh)]
    kept = [r for r in existing if str(r.get("database", "")) != str(database)]
    combined = kept + rows

    # Stable union of columns: existing order first, then any new keys in first-seen order.
    fieldnames: List[str] = []
    for r in combined:
        for k in r.keys():
            if k not in fieldnames:
                fieldnames.append(k)
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames)
        w.writeheader()
        for r in combined:
            w.writerow({k: r.get(k, "") for k in fieldnames})
    return out


# =============================================================================================
# 8. ROUNDING — round every numeric column to `sig` significant digits (byte-stable outputs)
# =============================================================================================
def round_numeric_cols(df_or_rows: Any, sig: int = 9) -> Any:
    """Round all numeric values to `sig` significant digits for byte-stable re-runs.

    Accepts a list-of-dicts (stdlib) OR a pandas DataFrame (lazy pandas) and returns the same
    type. Non-numeric values pass through unchanged. Strings that look numeric are NOT coerced
    (a CSV string column stays a string). Needs NO plotting backend.
    """
    # DataFrame path — only if pandas is actually available AND the input is a DataFrame.
    if _looks_like_dataframe(df_or_rows):
        import pandas as pd  # lazy
        df = df_or_rows.copy()
        for col in df.columns:
            if pd.api.types.is_numeric_dtype(df[col]):
                df[col] = df[col].map(lambda x: _round_sig(x, sig))
        return df

    rows = _to_rows(df_or_rows)
    out: List[Dict[str, Any]] = []
    for r in rows:
        nr: Dict[str, Any] = {}
        for k, v in r.items():
            nr[k] = _round_sig(v, sig) if isinstance(v, float) or isinstance(v, int) and not isinstance(v, bool) else v
        out.append(nr)
    return out


def _round_sig(x: Any, sig: int) -> Any:
    """Round a single numeric `x` to `sig` significant digits; pass through non-finite / non-numeric."""
    if isinstance(x, bool):
        return x
    if not isinstance(x, (int, float)):
        return x
    if x == 0 or not math.isfinite(float(x)):
        return x
    return round(float(x), sig - int(math.floor(math.log10(abs(float(x))))) - 1)


# =============================================================================================
# 9. DIRECTION CUE — map a sign to an unambiguous glyph/label (avoid a bare `*`); port of 14839
# =============================================================================================
def direction_cue(value: Any) -> str:
    """Map a signed value to an unambiguous directional glyph/label (never a bare `*`).

    Positive → "up" cue; negative → "down" cue; zero / non-finite / non-numeric → neutral.
    Glyphs are arrows + words so the cue is unambiguous in both color-blind and grayscale views.
    """
    try:
        v = float(value)
    except (TypeError, ValueError):
        return "· n/a"
    if not math.isfinite(v) or v == 0:
        return "· n.s."
    return "↑ up" if v > 0 else "↓ down"


# =============================================================================================
# Internal converters — keep the public API DataFrame-or-rows agnostic (testable without pandas)
# =============================================================================================
def _looks_like_dataframe(obj: Any) -> bool:
    """True only if pandas is importable AND `obj` is a pandas DataFrame (never imports eagerly)."""
    if obj is None or isinstance(obj, (list, tuple, dict)):
        return False
    try:
        import pandas as pd  # lazy
    except ImportError:
        return False
    return isinstance(obj, pd.DataFrame)


def _to_rows(df_or_rows: Any) -> List[Dict[str, Any]]:
    """Normalize input to a list-of-dicts. Accepts list-of-dicts (stdlib) or a DataFrame (lazy)."""
    if df_or_rows is None:
        return []
    if isinstance(df_or_rows, list):
        return [dict(r) for r in df_or_rows]
    if isinstance(df_or_rows, dict):
        return [dict(df_or_rows)]
    if _looks_like_dataframe(df_or_rows):
        return [dict(r) for r in df_or_rows.to_dict(orient="records")]
    raise TypeError(
        f"expected list-of-dicts or a pandas DataFrame, got {type(df_or_rows).__name__}"
    )


def _write_table_csv(rows: Any, path: Union[str, Path]) -> None:
    """Write a list-of-dicts (or DataFrame) to CSV with a stable column union (stdlib csv)."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    rows = _to_rows(rows)
    if not rows:
        path.write_text("")
        return
    fieldnames: List[str] = []
    for r in rows:
        for k in r.keys():
            if k not in fieldnames:
                fieldnames.append(k)
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames)
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in fieldnames})
