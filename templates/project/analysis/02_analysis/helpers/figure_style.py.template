"""
02_analysis/helpers/figure_style.py — per-project figure-style shim.
=====================================================================
ONE import per viz script. Delegates to the SciAgent-toolkit contract
lib (02_analysis/helpers/figure-style/figure_helpers.py, symlinked by
`sciagent activate`) and loads the project analysis_config.yaml.

Usage in any viz script:
    from helpers.figure_style import set_paper_style, save_overview, FIG_CFG
    set_paper_style(config=FIG_CFG)
    save_overview(fig, "04_gsea", "name", table=rows, ..., config=FIG_CFG)
"""
from __future__ import annotations

import importlib.util
import os
import sys
import warnings
from pathlib import Path

# ---------------------------------------------------------------------------
# 1. Resolve the symlinked contract lib and import it (graceful fallback).
# ---------------------------------------------------------------------------
_HERE = Path(__file__).parent
_LIB_DIR = _HERE / "figure-style"
_LIB_FILE = _LIB_DIR / "figure_helpers.py"

_FIGURE_STYLE_LOADED = False

if _LIB_FILE.exists():
    # Inject the lib dir into sys.path so `import figure_helpers` resolves.
    _lib_dir_str = str(_LIB_DIR)
    if _lib_dir_str not in sys.path:
        sys.path.insert(0, _lib_dir_str)
    from figure_helpers import (  # noqa: F401 — re-exported for callers
        load_figure_config,
        project_theme,
        set_paper_style,
        save_figure,
        save_overview,
        contrast_path,
        overview_path,
        style_series,
        purge_figures,
        write_caption,
        append_master_table,
        round_numeric_cols,
        direction_cue,
    )
    _FIGURE_STYLE_LOADED = True
else:
    warnings.warn(
        f"[figure_style] toolkit lib not found at: {_LIB_FILE}\n"
        "Run `sciagent activate` in this repo to link lib/figure-style/. "
        "Falling back to minimal stubs.",
        stacklevel=2,
    )

    # --- MINIMAL FALLBACK stubs (keep scripts from hard-crashing) ----------
    def load_figure_config(path: str = "02_analysis/config/analysis_config.yaml"):  # type: ignore[misc]
        """Fallback: load config from yaml if pyyaml is available."""
        try:
            import yaml
            with open(path) as fh:
                return yaml.safe_load(fh) or {}
        except Exception as exc:
            warnings.warn(f"[figure_style] Could not load config: {exc}", stacklevel=2)
            return {}

    def set_paper_style(base_size=None, variant=None, config=None):  # type: ignore[misc]
        """Fallback: no-op stub (`variant` accepted-but-ignored) — warns if lib absent."""
        warnings.warn("[figure_style] set_paper_style is a no-op stub (toolkit lib absent).",
                      stacklevel=2)

    def project_theme(base_size=None, legend=True, variant=None, config=None):  # type: ignore[misc]
        """Fallback: no-op stub (`variant` accepted-but-ignored)."""
        warnings.warn("[figure_style] project_theme is a no-op stub (toolkit lib absent).",
                      stacklevel=2)


# ---------------------------------------------------------------------------
# 2. Load the project config once (FIG_CFG is the stable project-wide handle).
# ---------------------------------------------------------------------------
try:
    FIG_CFG = load_figure_config("02_analysis/config/analysis_config.yaml")
except Exception as _exc:
    warnings.warn(f"[figure_style] Could not load analysis_config.yaml: {_exc}", stacklevel=1)
    FIG_CFG = {}
