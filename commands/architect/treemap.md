# /architecture-treemap

Render a `components.json` audit manifest into a self-contained interactive
`treemap.html`. Pure transform: no LLM, no network, no synthesis. Same input →
same output.

**Pair with `/diagram` for per-feature Mermaid; use this command for cross-
feature audit-level visualisation.**

## Usage

```
/architecture-treemap <path-to-components.json> [--output <path-to-html>]
```

- `<path-to-components.json>` — required. Path to a valid `components.json`
  that conforms to the architecture-treemap schema.
- `--output <path>` — optional. Default: `<input-dir>/treemap.html`.

If the path is omitted the command looks for the most recent directory under
`docs/_meta/architecture-audit/` matching `YYYY-MM-DD` and uses its
`components.json`.

## Phase 0: Resolve input path

1. If `$ARGUMENTS[0]` is supplied, use it as `<input>`.
2. If omitted, glob `docs/_meta/architecture-audit/????-??-??/components.json`.
   Pick the lexicographically largest directory (ISO date sorts as chronological).
   If none found:

   ```
   No components.json found. Tried:
     - docs/_meta/architecture-audit/YYYY-MM-DD/components.json

   The components.json must conform to the schema at
     skills/architecture-treemap/components.schema.json

   Run /synthesize-audit to generate one, or author it manually.
   ```

   Stop.

## Phase 1: Validate

Run the validator before rendering. Print results and stop on error.

```bash
python skills/architecture-treemap/scripts/validate_components.py <input>
```

On schema error → print errors and stop. Half-rendered treemaps are worse than none.

Referential-integrity warnings are printed but do not block rendering.

## Phase 2: Render

```bash
python skills/architecture-treemap/scripts/render_treemap.py <input> [--output <output>]
```

The renderer:
- Injects the validated JSON into the HTML template.
- Inlines D3 v7 (vendored asset — no CDN required).
- Inlines the `treemap.bundle.js` renderer.
- Writes one self-contained `.html` openable via `file://`.

## Phase 3: Report

Report to the user:

```
Rendered: <absolute-path-to-html>
Open with: file://<absolute-path>

Views available:
  Logical  — treemap by size_estimate_loc, coloured by classification
  Physical — treemap by size_loc, coloured by primary logical-owner classification
  Graph    — force-directed layout with three classification-gravity clusters

Edge overlay (Logical + Graph):
  Click a cell/node → edges for that component light up.
  Shift-click adds. Background-click clears.
  "Show all edges (faint)" toggle in header.

Shift-click nodes in Graph view to compare edge sets.
The force-graph seed is in the URL hash (#<integer>) for reproducible screenshots.
```

If the manifest has no `logical_components` (empty array), note:

```
NOTE: logical_components is empty. The Logical and Graph views will show a
"No logical model yet" notice. Run /synthesize-audit to populate them.
The Physical view renders from physical_components alone.
```
