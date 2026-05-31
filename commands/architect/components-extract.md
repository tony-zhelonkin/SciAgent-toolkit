# Components-extract

Produce the **deterministic substrate** of an architecture audit: a `components.json` whose every entry traces back to the AST, the git history, or a metric tool — never to judgment. This is the foundation the rest of the `architecture-treemap` stack rests on.

The substrate biases attention toward hotspots (high churn, high fan-in, high complexity) and lets a reviewer sanity-check the synthesis's judgment against ground truth. It emits NO classifications, smells, background-knowledge edges, prune verdicts, or core/seam/removable boundary — those are not statically detectable and are authored later by `/audit-slice` and `/synthesize-audit`.

Pairs with `/architecture-treemap`, which renders a synthesized manifest.

## Phase 0: Parse arguments

- `$ARGUMENTS[0]` — path to the repository root to scan (required). Reject if missing: "Usage: /components-extract <repo-root>".
- `--out <path>` — output path (default: `<repo-root>/components.json`).
- `--since-days <N>` — churn window in days (default: 90).
- `--lang <python>` — source language for the import graph (only `python` is supported today; other values warn and skip the graph).
- `--include-nonsource` — include non-source scaffolding (default: OFF; see below).
- `--include-docs` — include documentation files (default: OFF).
- `--include-audit-output` — include the skill's own prior snapshots (default: OFF).

## Phase 1: Run the extractor

```bash
python skills/architecture-treemap/scripts/extract_components.py <repo-root> \
    [--out <components.json>] [--since-days 90] [--lang python]
```

The script is a pure function of the tree state at the current git HEAD plus the optional tools installed. It prints a one-line summary to stderr (files scanned, python files, static edges, unparseable files skipped, metric tools used, schema-validation result).

## What it emits

- **`physical_components[]`** — one per *architectural* source/config/test file under the repo. `.gitignore` is honored (via `git ls-files`); vendored and build dirs (`node_modules`, `dist`, `.venv`, `__pycache__`, `*.egg-info`, ...) are skipped. Each carries `id` (collision-free kebab-case of the path), `path`, `size_loc`, `kind` (module / config / test / doc / asset / other, by extension + location), `logical_owners: []` (judgment — left empty), a `metrics{}` block, and (if it sits in an import cycle) `cycle_id` + `cycle_size`.
- **`cycles[]`** — import cycles found by Tarjan strongly-connected-components over the static import graph (deterministic). Each is `{id, members[]}` for an SCC of size > 1. Teaches the Acyclic Dependencies Principle.

### Non-source inventory filter (default ON)

`git ls-files` honors `.gitignore` for UNTRACKED files only. This filter is the complement: it drops TRACKED scaffolding that is not architecture and would otherwise drown the signal. **Excluded by default:**

- license files (`LICENSE*`, `COPYING*`), lockfiles (`package-lock.json`, `poetry.lock`, `*.lock`, ...);
- root-level tooling/dotfiles (`.gitignore`, `.mcp.json`, `package.json`, `pyproject.toml`, ...);
- harness/editor/CI dirs (`.claude/`, `.github/`, `.vscode/`, `.idea/`, ...);
- data files by extension (`.csv`, `.tsv`, `.parquet`, `.rds`, `.npy`, `.pkl`, ...);
- **everything under a top-level `docs/` directory, by path-prefix** (this catches `docs/**/*.csv` that an extension-only filter misses).

**Kept:** source modules and tests (including `scripts/verify_*.py` and one-off scripts — those are real); and config/schema/asset files that live UNDER a source root (`src/` or a discovered package dir), e.g. `src/**/plugin_payload.schema.json` — a wire-contract schema belongs with the code it governs. Root-level tooling JSON is dropped; src-resident schema JSON is kept.

Pass `--include-nonsource` to restore the excluded set (e.g. when auditing a tooling repo where the configs ARE the architecture).
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

Two more deterministic metrics are derived from the import graph:

| Metric | Source | Notes |
|---|---|---|
| `instability` | derived | `I = fan_out / (fan_in + fan_out)`, in [0,1]. Undefined (omitted) for isolated nodes. Teaches the Stable-Dependencies Principle. |
| `refactor_pressure` | derived | `churn * cyclomatic * max(fan_in,1) / max(test_ratio, 0.05)`. **Absent when radon is** (cyclomatic missing). |

### Optional dependency — LOUD degradation

`radon` is optional (`pip install radon`) and runs only at extract time — it does not affect the portable HTML. But its absence is announced LOUDLY and honestly, never silently:

- The extractor prints a `WARNING` to stderr naming the unavailable metrics, and stamps them in `extractor.unavailable_metrics` (`["cyclomatic", "refactor_pressure"]`) with a per-metric `degradation_notes`.
- When radon is absent, a **clearly-labelled degraded proxy** `refactor_pressure_loc_proxy` is emitted (LOC/100 substituting for complexity) so the "look here first" lens still ranks something. The renderer shows it marked as a degraded proxy and, if even the proxy is unavailable, greys the lens with "install radon to enable complexity & refactor-pressure" rather than showing an empty lens.
- `jsonschema` is also optional; without it self-validation is SKIPPED with a warning, and the stdlib `validate_components.py` is the authoritative gate.

`cloc`/`tokei` are NOT depended on; the stdlib LOC counter is the deterministic floor.

### Scaling-regime hint

The stderr summary prints a regime hint keyed to component + edge counts: *small* (one `/audit-slice` + `/synthesize-audit` pass holds the whole graph) vs *large* (scope concerns narrowly, expect multiple passes, ~1M tokens). The deterministic substrate does not degrade with size — only the judgment layer does.

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
