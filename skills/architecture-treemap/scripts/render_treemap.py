#!/usr/bin/env python3
"""
render_treemap.py — Deterministic HTML renderer for architecture-treemap manifests.

Usage:
    python render_treemap.py <components.json> [--output <treemap.html>]

Default output: <input-dir>/treemap.html

The renderer:
  1. Validates components.json via validate_components.py (schema + referential integrity).
     Refuses to render on any schema error.
  2. Injects the validated JSON, D3 v7 (vendored), and the bundle JS into the HTML template.
  3. Writes a single self-contained file openable via file:// (no CDN, no server needed).

Determinism guarantee: same input → same output, modulo no runtime state.
Force-graph simulation seed is in the URL hash (#<integer>) for reproducibility.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

# Resolve paths relative to this script
_SCRIPT_DIR   = Path(__file__).parent
_SKILL_DIR    = _SCRIPT_DIR.parent
_TEMPLATE     = _SKILL_DIR / "assets" / "treemap.template.html"
_BUNDLE       = _SKILL_DIR / "assets" / "treemap.bundle.js"
_D3_VENDOR    = _SKILL_DIR / "assets" / "d3.v7.min.js"
_VALIDATE_MOD = _SCRIPT_DIR / "validate_components.py"

# ── Import validator ──────────────────────────────────────────────────────────
# We import validate_components directly rather than shelling out so that
# errors are richer and there is no subprocess overhead.

import importlib.util

def _load_validator():
    spec = importlib.util.spec_from_file_location("validate_components", _VALIDATE_MOD)
    mod  = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


# ── Rendering ─────────────────────────────────────────────────────────────────

def _load_text(path: Path, label: str) -> str:
    if not path.exists():
        print(f"ERROR: {label} not found at {path}", file=sys.stderr)
        sys.exit(1)
    return path.read_text(encoding="utf-8")


def _build_title(data: dict) -> str:
    project = data.get("project", "Architecture audit")
    date    = data.get("audit_date", "")
    return f"{project} — {date}" if date else project


def _build_frame_badge(data: dict) -> str:
    frame = data.get("frame")
    if not frame:
        return ""
    display = frame.replace("-", " → ")
    return f'<span class="frame">{display}</span>'


def _build_counts(data: dict) -> str:
    n_logical  = len(data.get("logical_components", []))
    n_physical = len(data.get("physical_components", []))
    n_edges    = len(data.get("edges", []))
    n_prune    = len(data.get("prune_candidates", []))
    parts = [f"{n_logical} logical", f"{n_physical} physical", f"{n_edges} edges"]
    if n_prune:
        parts.append(f"{n_prune} prune candidates")
    return " &nbsp;·&nbsp; ".join(parts)


def _build_stamp(data: dict) -> str:
    sha      = data.get("git_sha")
    snap_id  = data.get("snapshot_id") or data.get("audit_date", "")
    parts = []
    if snap_id:
        parts.append(snap_id)
    if sha:
        short = sha[:8] if len(sha) > 8 else sha
        if short != "00000000":  # suppress zero-placeholder shas
            parts.append(f"git:{short}")
    return " &nbsp;·&nbsp; ".join(parts) if parts else ""


def render(components_path: Path, output_path: Path) -> None:
    # ── Load and validate ──────────────────────────────────────────────────────
    validator = _load_validator()
    schema_errors, ref_issues, soft_warnings = validator.run_validation(
        components_path, strict=False
    )

    if schema_errors:
        print(f"RENDER REFUSED — schema errors in {components_path}:")
        for err in schema_errors:
            print(f"  ERROR: {err}")
        print("\nFix the errors above and re-run.")
        sys.exit(1)

    if ref_issues:
        print("Referential-integrity warnings (rendering anyway — use --strict on validate to error):")
        for w in ref_issues:
            print(f"  WARNING: {w}")

    if soft_warnings:
        for w in soft_warnings:
            print(f"  NOTICE: {w}")

    with components_path.open("r", encoding="utf-8") as fh:
        data = json.load(fh)

    # ── Load assets ────────────────────────────────────────────────────────────
    template     = _load_text(_TEMPLATE, "HTML template")
    bundle_js    = _load_text(_BUNDLE,   "treemap bundle JS")
    d3_js        = _load_text(_D3_VENDOR, "D3 v7 vendor JS")

    # ── Build substitutions ────────────────────────────────────────────────────
    # json.dumps with sort_keys ensures same input → same output (dict key order stable).
    components_json = json.dumps(data, ensure_ascii=False, sort_keys=True)

    substitutions = {
        "{{TITLE}}":          _build_title(data),
        "{{PROJECT_NAME}}":   data.get("project", "Architecture audit"),
        "{{FRAME_BADGE}}":    _build_frame_badge(data),
        "{{COUNTS}}":         _build_counts(data),
        "{{STAMP}}":          _build_stamp(data),
        "{{COMPONENTS_JSON}}": components_json,
        "{{D3_BUNDLE}}":      d3_js,
        "{{TREEMAP_BUNDLE}}": bundle_js,
    }

    html = template
    for placeholder, value in substitutions.items():
        html = html.replace(placeholder, value)

    # Sanity-check: no placeholders remain
    remaining = [k for k in substitutions if k in html]
    if remaining:
        print(f"WARNING: template placeholders not replaced: {remaining}", file=sys.stderr)

    # ── Write output ───────────────────────────────────────────────────────────
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(html, encoding="utf-8")
    print(f"OK: rendered {output_path}")
    print(f"    {len(html):,} bytes | open with: file://{output_path.resolve()}")


# ── CLI ───────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Render a components.json manifest into a self-contained treemap.html.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("components_json", help="Path to components.json")
    parser.add_argument(
        "--output", "-o",
        default=None,
        help="Output HTML path. Default: <input-dir>/treemap.html",
    )
    args = parser.parse_args()

    components_path = Path(args.components_json).resolve()
    if not components_path.exists():
        print(f"ERROR: file not found: {components_path}", file=sys.stderr)
        sys.exit(1)

    output_path = Path(args.output).resolve() if args.output else (components_path.parent / "treemap.html")

    render(components_path, output_path)


if __name__ == "__main__":
    main()
