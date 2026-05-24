# Components-extract

Produce the **deterministic substrate** of an architecture audit: a `components.json` whose every entry traces back to the AST, the git history, or a metric tool — never to judgment. This is the foundation the rest of the `architecture-treemap` stack rests on.

The substrate biases attention toward hotspots (high churn, high fan-in, high complexity) and lets a reviewer sanity-check the synthesis's judgment against ground truth. It emits NO classifications, smells, background-knowledge edges, prune verdicts, or core/seam/removable boundary — those are not statically detectable and are authored later by `/audit-slice` and `/synthesize-audit`.

Pairs with `/architecture-treemap`, which renders a synthesized manifest.

## Phase 0: Parse arguments

- `$ARGUMENTS[0]` — path to the repository root to scan (required). Reject if missing: "Usage: /components-extract <repo-root>".
- `--out <path>` — output path (default: `<repo-root>/components.json`).
- `--since-days <N>` — churn window in days (default: 90).
- `--lang <python>` — source language for the import graph (only `python` is supported today; other values warn and skip the graph).

## Phase 1: Run the extractor

```bash
python skills/architecture-treemap/scripts/extract_components.py <repo-root> \
    [--out <components.json>] [--since-days 90] [--lang python]
```

The script is a pure function of the tree state at the current git HEAD plus the optional tools installed. It prints a one-line summary to stderr (files scanned, python files, static edges, unparseable files skipped, metric tools used, schema-validation result).

## What it emits

- **`physical_components[]`** — one per source/config/test/doc/asset file under the repo. `.gitignore` is honored (via `git ls-files`); vendored and build dirs (`node_modules`, `dist`, `.venv`, `__pycache__`, `*.egg-info`, ...) are skipped. Each carries `id` (collision-free kebab-case of the path), `path`, `size_loc`, `kind` (module / config / test / doc / asset / other, by extension + location), `logical_owners: []` (judgment — left empty), and a `metrics{}` block.
- **`edges[]`** — the Python intra-repo import graph at the physical level. Each edge is `type: "direct-call"`, `evidence_class: "static"` (renders SOLID), with an `evidence` note of the form `path:line -> import module`. Stdlib and third-party imports are dropped; only imports that resolve to a file in the repo become edges. `src/<pkg>/` layout is handled.
- **`logical_components: []`** — intentionally EMPTY. Logical grouping and classification are judgment; `/synthesize-audit` fills this. The schema requires the key but allows an empty array.
- **Top-level provenance** — `schema_version: "1.0"`, `audit_date` (today), `project` (repo dir name), `git_sha` (current HEAD), `snapshot_id`, and an `extractor{}` block recording `tool_version`, `extracted_at`, and which metric tools actually ran.

## The five locked metrics

Attached per physical component in `metrics{}`. All optional in the schema — a metric is omitted (never fabricated) when its tool is absent or it cannot be computed for that file.

| Metric | Source | Notes |
|---|---|---|
| `loc` | stdlib counter | Python: non-blank, non-comment-only lines. Other kinds: raw line count. |
| `fan_out` | `ast` import graph | Distinct intra-repo modules this file imports. |
| `fan_in` | `ast` import graph | Count of repo files importing this one. |
| `cyclomatic` | `radon` (optional) | Mean cyclomatic complexity across blocks. **Omitted gracefully if `radon` is not importable** — there is no hard dependency. |
| `churn_90d` | `git log --since` | Commits touching the file in the window. Skipped when the tree is not a git repo. |
| `test_ratio` | rollup | A repo-level test-LOC / source-LOC ratio attached to test files' metrics. Per-file test_ratio is not statically attributable; the true per-component value is a synth-level rollup once logical grouping exists. |

### Optional dependency

`radon` is the only optional dependency (`pip install radon`). Without it, every other metric still computes and the run records the degraded tool set in `extractor.metric_tools` (e.g. `["ast", "stdlib-loc", "git-log"]` with no `radon`). `cloc`/`tokei` are NOT depended on; the stdlib LOC counter is the deterministic floor.

## Immutability and snapshot semantics

Each run is an **immutable, git-SHA-stamped snapshot**. `git_sha` records the HEAD the substrate was extracted at; `snapshot_id` and `audit_date` datestamp it. Re-running produces a NEW `components.json` — it never mutates a prior snapshot in place. To keep a history, write each run to a dated path (`--out audit/<date>/components.json`).

## Determinism

The extractor is deterministic modulo the wall-clock `extracted_at` field: same tree state + same tool set + same flags produces the same `physical_components`, `edges`, and metric values. File ordering is sorted; id disambiguation suffixes are assigned in sorted order.

## Self-validation

After writing, the script validates its output against `components.schema.json` if `jsonschema` is importable (best-effort; reported on stderr). The authoritative gate is `validate_components.py`, which `/architecture-treemap` runs before rendering and which refuses to render on any schema or referential-integrity error.

## Next steps

```
/audit-slice <hunch>     — author one judgment slice (human-triggered; refuses auto-piped input)
/synthesize-audit        — integrate slices + substrate into the full components.json (fills logical_components, classifications, smells, prune verdicts, core_boundary)
/architecture-treemap <path>  — render the synthesized manifest to a self-contained treemap.html
```

## Rules

1. **Substrate is facts, not assessment.** If asked to emit a classification, smell, or background-knowledge edge, refuse — those are judgment, authored by `/synthesize-audit`. Emitting a guess would defeat the design's trust scaffold.
2. **Omit, never fabricate.** A metric whose tool is absent is left out of the `metrics` block; it is never invented or zero-filled.
3. **One file out.** The extractor writes exactly one `components.json`.
4. **Never crash on a bad file.** A single unparseable source file is skipped and counted in the stderr summary; the run completes.
