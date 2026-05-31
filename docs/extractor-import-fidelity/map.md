---
date: 2026-05-26
feature: extractor-import-fidelity
kind: map
repo_root: /workspaces/DC_hum_verse/01_modules/SciAgent-toolkit
---

# Map: extractor-import-fidelity

## Summary

The architecture-treemap extractor `extract_components.py` resolves relative
imports (`from .X import Y`) against a *package* it derives by stripping the
filename from the module's dotted name. For a regular module `a/b.py` that
stripping is correct (package = `a`). For a package `a/__init__.py` the module
IS package `a`, so stripping yields the *parent* package — one level too
shallow. Every relative re-export out of a package `__init__.py` therefore
mis-resolves: the resolved base falls back through the candidate ladder to the
top package root, collapsing N sibling re-export edges into ONE mis-targeted
edge whose `evidence` string does not match its source line. The defect skews
the static import graph, the fan_in / fan_out / instability / refactor_pressure
metrics derived from it, the logical-edge rollup, the orphan validation, and the
rendered treemap. A separate, related finding: the validator's `--strict`
orphan gate counts a `test -> node` edge as "connected," so a production-orphan
created by this bug passes `--strict`.

All claims below are verified against the live source and against a live
re-run of the resolver over `pathway-explorer` (the empirical witness repo at
`/workspaces/DC_hum_verse/01_modules/pathway-explorer`).

## Blast radius

### Files (touched by this feature)

| File | Role | Nature of involvement |
|------|------|------------------------|
| `skills/architecture-treemap/scripts/extract_components.py:458-480` | `module_name_for_file` | Strips trailing `__init__` (L468-469); produces the dotted module name whose `[:-1]` is taken as the package — the root-cause input. |
| `skills/architecture-treemap/scripts/extract_components.py:493-525` | `collect_imports` | Computes `package_parts = source_module.split(".")[:-1]` (L510); drives `ImportFrom` handling (L517-524). |
| `skills/architecture-treemap/scripts/extract_components.py:528-536` | `resolve_import_from_base` | The relative-import resolver. L533 walks up `level` packages from `package_parts`; the wrong `package_parts` for `__init__.py` makes the base one level too shallow. |
| `skills/architecture-treemap/scripts/extract_components.py:539-550` | `import_from_candidates` | Builds `base.symbol` then bare `base` candidates (L547-549); when `base` is mis-resolved, the bare-`base` fallback later collapses to the package root. |
| `skills/architecture-treemap/scripts/extract_components.py:553-566` | `resolve_to_repo_module` | Trims trailing dotted segments until a repo module matches (L561-565); this is the mechanism that turns a wrong base into an edge to the *top package root* `__init__.py`. |
| `skills/architecture-treemap/scripts/extract_components.py:569-618` | `build_import_graph` | Where edges, `fan_out`, `fan_in` are accumulated (L598-599) and the `evidence` string is formatted (L611-614) from `relative:line -> import target_module`. |
| `skills/architecture-treemap/scripts/extract_components.py:483-490` | `build_module_index` | Maps dotted module -> repo file; `__init__` packages map to their dir name (via `module_name_for_file`). The lookup table the resolver hits. |
| `skills/architecture-treemap/scripts/metric_registry.py:101-128` | `_instability` | Derives `instability = fan_out/(fan_in+fan_out)` from the (now-skewed) fan counts. |
| `skills/architecture-treemap/scripts/metric_registry.py:131-152` | `_refactor_pressure` | Uses `fan_in` (skewed for package-init targets) as a multiplier. |

### Consumers of the emitted edges (downstream blast radius)

| File | Role | Nature of involvement |
|------|------|------------------------|
| `skills/architecture-treemap/scripts/rollup_logical_edges.py:102-159` | `derive_logical_edges` | Reads `edges` where `evidence_class == "static"` (L121) and both endpoints are physical ids (L125); rolls physical import pairs up to logical edges. A lost facade edge never rolls up. |
| `skills/architecture-treemap/scripts/validate_components.py:374-457` | `validate_referential_integrity` | Orphan gate. L421-428 mark a node "connected" if it appears in ANY edge — including a `test -> node` edge (the SECOND finding). |
| `skills/architecture-treemap/scripts/validate_components.py:536-537` | `main` (`--strict`) | `--strict` exits 1 only on `ref_issues`; a test-only-connected node yields no issue, so the run_metadata defect passes `--strict`. |
| `skills/architecture-treemap/scripts/render_treemap.py:94-153` | `render` | Injects `edges` verbatim into the HTML (L127, L135); renders whatever the extractor emitted. |
| `skills/architecture-treemap/assets/treemap.bundle.js:13-26` | edge overlay | Draws outgoing/incoming edges per `evidence_class`; a missing edge simply never lights. |
| `/workspaces/DC_hum_verse/01_modules/pathway-explorer/docs/_meta/architecture-audit/2026-05-25-postrefactor/components.json` | committed snapshot | Contains the mis-resolved edges (verified). Re-extraction will change its edge set — immutability/regeneration implication (see Open questions). |

### Immutable snapshots already on disk (re-extraction changes their edge counts)

In the witness repo `pathway-explorer`, four committed `components.json`
snapshots exist under `docs/_meta/architecture-audit/`:
- `2026-05-21/`, `2026-05-24/`, `2026-05-24-post-refactor/`, `2026-05-25/`,
  `2026-05-25-postrefactor/components.json` (5 dated dirs; the last verified in detail).

These are excluded from re-extraction *as inputs* by
`extract_components.py:333-350` (`is_audit_output`, the `architecture-audit`
path segment). They are not regenerated automatically; a fresh extraction
produces a NEW snapshot. Whether prior snapshots are re-synthesised is a design
decision (Open questions, item c).

### Tests

- No Python unit tests exist for the extractor. `find` over the toolkit returns
  only vendored test files under `mcp_servers/pal/venv/.../site-packages/` and
  one unrelated `docs/.ref/era/implementation/futs_test.py`.
- `skills/architecture-treemap/tests/` does not exist.
- `skills/architecture-treemap/checks/golden_render.sh` — end-to-end check:
  validates + renders the golden fixture, asserts self-containment and
  byte-stability (L1-end). It does NOT exercise import resolution; it renders a
  pre-baked `components.json`.
- `skills/architecture-treemap/references/example/components.json` — the golden
  render fixture (a hand-authored manifest, not extracted from a source tree).
- The toolkit's `tests/run-all.sh` suite (`tests/test_*.sh`) covers the bash
  `sciagent` CLI only — no extractor coverage.

### Configs / data files

- `skills/architecture-treemap/components.schema.json` — the edge schema the
  validator mirrors. `evidence`/`evidence_class` are optional strings/enums
  (`validate_components.py:34-42, 180-210`); the schema does NOT constrain the
  `evidence` string's relationship to its source line, so a mismatched evidence
  string is schema-valid.
- `skills/architecture-treemap/scripts/metric_registry.py:159-273` — the metric
  descriptor table (single source of truth) declaring `fan_in`, `fan_out`,
  `instability`, `refactor_pressure`.

### Out of scope (confirmed NOT touched)

- `skills/architecture-treemap/scripts/pedagogy_registry.py` — explainer/glossary
  text; no edge or import logic.
- Cyclomatic / churn / LOC computation paths (`count_loc`, `cyclomatic_for_file`,
  `churn_for_file`) — independent of the import graph.
- `find_cycles` (`extract_components.py:625-698`) consumes edges but is Tarjan
  SCC over whatever edges exist; it is downstream-affected only insofar as the
  edge set changes (no resolver logic of its own).
- The non-source scaffolding filter (`is_nonsource_scaffolding`,
  `extract_components.py:233-278`) — file inventory, not import resolution.
- Absolute imports from `__init__.py` (e.g.
  `contrast_designs/__init__.py:6 from pathway_explorer.contrast_designs.f1_f4`)
  resolve CORRECTLY (verified) — they take the `node.level == 0` path
  (`resolve_import_from_base:530-531`) and never touch `package_parts`.

## Root cause (verified, not assumed)

The relative-import resolver is `resolve_import_from_base`
(`extract_components.py:528-536`):

```
package_parts = source_module.split(".")[:-1]      # collect_imports:510
base_parts = package_parts[: len(package_parts) - (node.level - 1)]   # L533
```

`source_module` comes from `module_name_for_file`
(`extract_components.py:458-480`), which strips a trailing `__init__`
(L468-469). The consequence:

| File | `module_name_for_file` | `package_parts` (L510) | Correct package |
|------|------------------------|------------------------|-----------------|
| `src/pathway_explorer/html_generator.py` (regular) | `pathway_explorer.html_generator` | `['pathway_explorer']` | `pathway_explorer` ✓ |
| `src/pathway_explorer/api/__init__.py` (package init) | `pathway_explorer.api` | `['pathway_explorer']` | `pathway_explorer.api` ✗ (one short) |
| `src/pathway_explorer/rendering/js/__init__.py` (deep init) | `pathway_explorer.rendering.js` | `['pathway_explorer','rendering']` | `pathway_explorer.rendering.js` ✗ (one short) |

For a regular module, `[:-1]` correctly drops the module's own leaf to leave its
containing package. For an `__init__.py`, `module_name_for_file` already dropped
`__init__`, so the dotted name IS the package; taking `[:-1]` again drops a real
package segment, landing the relative base one level above where it should.

Verified resolution trace for `api/__init__.py` (live re-run of `collect_imports`
+ `resolve_to_repo_module` against the witness repo):

```
line 39 from .channels  -> base 'pathway_explorer.channels'   -> resolve_to_repo_module trims to 'pathway_explorer' -> src/pathway_explorer/__init__.py
line 40 from .payload   -> base 'pathway_explorer.payload'    -> 'pathway_explorer' -> root __init__
line 41 from .run_metadata -> base 'pathway_explorer.run_metadata' -> 'pathway_explorer' -> root __init__
line 42 from .projection-> base 'pathway_explorer.projection' -> 'pathway_explorer' -> root __init__
line 48 from .sink      -> base 'pathway_explorer.sink'       -> 'pathway_explorer' -> root __init__
```

The wrong base (`pathway_explorer.channels`) is not in the module index;
`resolve_to_repo_module` (`L553-566`) then trims trailing segments
(`pathway_explorer.channels` -> `pathway_explorer`), which IS in the index (the
package root `__init__.py`). So all 5 distinct sibling re-exports collapse to a
single edge to the package root, deduped by `seen_edges`
(`build_import_graph:601-604`) to ONE edge. The emitted `evidence` is
`"src/pathway_explorer/api/__init__.py:39 -> import pathway_explorer"`
(`build_import_graph:611-614`) — the line number is the first surviving
relative import (39), but `import pathway_explorer` does not match line 39's
actual source (`from .channels import ...`). Confirmed in the committed snapshot
`docs/_meta/architecture-audit/2026-05-25-postrefactor/components.json`.

### The visible casualty (verified)

`src/pathway_explorer/api/run_metadata.py` — its only correct production importer
is the api facade (`api/__init__.py:41`). With the facade edge lost, the snapshot
shows its only inbound edge is from a test:
`tests/test_run_metadata_contract.py:78 -> import pathway_explorer.api.run_metadata`.
It renders as a test-only / production-orphan node. A human already authored a
compensating `audit-asserted` edge to patch ONE of the five lost edges:
`api/__init__.py:40 -> from .channels import ... (public re-export; the static
extractor did not lift this through the facade).` — i.e. the defect was noticed
and worked around by hand, not at the resolver.

### The four package `__init__.py` files in the witness repo (all affected)

Verified live resolution (`collect_imports` + `resolve_to_repo_module`):

| Package init | relative imports | resolved (buggy) target | correct target |
|--------------|------------------|-------------------------|----------------|
| `src/pathway_explorer/__init__.py` | `.config`, `.data_loader`, `.embedding`, `.api`, `.html_generator`, `.main`, `.similarity` | mis-resolve to siblings of `pathway_explorer` (root has no parent → resolve to root or external) | the named sibling modules |
| `src/pathway_explorer/api/__init__.py` | 5 (`.channels`,`.payload`,`.run_metadata`,`.projection`,`.sink`) | all → `pathway_explorer` root, deduped to 1 edge | `api/channels.py` … `api/sink.py` |
| `src/pathway_explorer/rendering/__init__.py` | 6 (`.chrome`,`.document`,`.index_page`,`.js`,`.styles`,`.worker`) | all → `pathway_explorer` root, deduped to 1 edge | `rendering/chrome.py` … `rendering/worker.py` |
| `src/pathway_explorer/rendering/js/__init__.py` | 9 (`.bootstrap`,`.chart`,…,`.selection`) | all → `pathway_explorer.rendering` (one short), deduped to 1 edge | `rendering/js/bootstrap.py` … `selection.py` |

(There is a 5th package init, `contrast_designs/__init__.py`, but it uses an
ABSOLUTE import and resolves correctly — see Out of scope.)

## Metrics contamination (which metrics, which nodes)

Because the resolver under-counts facade fan-out and mis-attributes inbound edges:

- **`fan_out`** (`build_import_graph:598`, surfaced `extract_components.py:851-852`):
  every package `__init__.py` shows fan_out collapsed toward 1 (or 0) instead of
  its true re-export count. E.g. `api/__init__.py` should have ~5 outbound, the
  snapshot collapses it.
- **`fan_in`** (`build_import_graph:599`, surfaced `L853-854`): the package ROOT
  `__init__.py` gains spurious inbound edges from every sibling init that
  mis-resolved into it; the snapshot shows `pathway_explorer/__init__.py` with
  `fan_in: 2` from two mis-resolved init edges. The re-export TARGETS (e.g.
  `run_metadata.py`) lose their true production `fan_in`.
- **`instability`** (`metric_registry.py:101-128`): `fan_out/(fan_in+fan_out)` is
  recomputed on the wrong counts; verified `run_metadata.py instability: 0.0`
  (pure sink) is computed from `fan_in: 1` (test-only) instead of its true
  production fan_in.
- **`refactor_pressure`** (`metric_registry.py:131-152`): uses `fan_in` as a
  multiplier (`max(fan_in,1)`, L151) — under-counted fan_in for facade targets
  understates pressure.
- **`test_ratio` / coverage** (`compute_module_coverage:902-930`): a module
  "covered" iff a test imports it; this is independent of the resolver bug, BUT
  the renderer's "production orphan vs test-covered" reading is distorted because
  `run_metadata.py` looks test-only.

## Entry points & data flow

1. CLI entry `main` (`extract_components.py:1158`) → `build_manifest`
   (`L947-1047`).
2. File discovery: `list_repo_files` (`L285-303`) → `python_files`
   (`build_manifest:954`).
3. Import roots: `find_python_roots` (`L443-455`) → `module_index =
   build_module_index` (`L483-490`) — dotted module name → repo file. Type:
   `dict[str, str]`.
4. Graph build: `build_import_graph` (`L569-618`) iterates each python file;
   per file → `module_name_for_file` (`L458`, the `__init__`-stripping step) →
   `collect_imports` (`L493`) → for each `ImportFrom`:
   `resolve_import_from_base` (`L528`, **defect site**) → `import_from_candidates`
   (`L539`) → `resolve_to_repo_module` (`L553`, the collapse-to-root step) →
   edge dict with `evidence` (`L605-615`). Type: candidate `list[str]` →
   `target_module: str` → `target_file: path` → `edge: dict`.
5. Fan counts accrue during step 4 (`L598-599`). Metrics derived in
   `build_physical_components` (`L797-899`): base fan_in/fan_out copied
   (`L851-854`), then `metric_registry.compute_derived` (`L874`) adds
   `instability` + `refactor_pressure`.
6. Output `components.json` written in `main` (`L1175-1178`). Type: manifest dict
   with `edges`, `physical_components[].metrics`.
7. Downstream: `rollup_logical_edges.py` reads `edges` → derived logical edges;
   `validate_components.py` reads `edges` → orphan gate; `render_treemap.py`
   reads everything → `treemap.html`.

```mermaid
flowchart LR
    A[extract_components.py:1158 main] --> B[build_manifest L947]
    B --> C[build_module_index L483]
    B --> D[build_import_graph L569]
    C --> D
    D --> E[module_name_for_file L458\nstrips __init__]
    E --> F[collect_imports L510\npackage_parts = name[:-1]]
    F --> G[resolve_import_from_base L528\nDEFECT: base one level short for __init__]
    G --> H[import_from_candidates L539]
    H --> I[resolve_to_repo_module L553\ntrims to package root]
    I --> J[edge + evidence L611\nfan_in/fan_out L598]
    J --> K[compute_derived L874\ninstability, refactor_pressure]
    J --> L[components.json L1175]
    L --> M[rollup_logical_edges.py:121]
    L --> N[validate_components.py:421 orphan gate\n--strict counts test->node]
    L --> O[render_treemap.py:127 -> treemap.html]
    L --> P[committed snapshots docs/_meta/...]
```

## The second finding (related; do NOT conflate)

`validate_referential_integrity` (`validate_components.py:374-457`) builds
`logical_nodes_with_edges` (L421-428) by marking a node connected if it appears
as `from` OR `to` in ANY edge whose endpoint is a logical id — with NO filter on
the edge's *source kind*. The rollup attaches `test -> node` edges as logical
edges (the test suite has logical owners), so a node whose only edge is from a
test counts as connected. `--strict` (`main:536-537`) exits 1 only on
`ref_issues`; the production-orphan never becomes an issue. This is the blind
spot that let the `run_metadata` defect render without failing validation. It is
edge-CLASSIFICATION logic, not import RESOLUTION — a sibling concern, flagged for
the design stage to scope in or out (Open questions, item b).

## Existing patterns to reuse

- **Candidate-ladder resolution** — `import_from_candidates`
  (`extract_components.py:539-550`) + `resolve_to_repo_module` (`L553-566`):
  the established pattern for "offer most-specific dotted name first, trim to the
  package." A corrected `package_parts` feeds straight into this unchanged
  machinery — the analog the fix can reuse without restructuring.
- **`module_name_for_file`** (`L458-480`) already knows a file is an `__init__`
  (it strips it at L468-469); that same predicate is the natural place to derive
  the correct relative-import base — the package-detection signal already exists
  in the file.
- **Registry-driven derived metrics** (`metric_registry.py:306-320`,
  `compute_derived`): metrics auto-recompute from corrected fan counts; no metric
  code changes once edges are right.

## Key files

- `skills/architecture-treemap/scripts/extract_components.py:528-536` — the
  relative-import resolver (the single defect site; fix lives here / in the
  `package_parts` derivation feeding it).
- `skills/architecture-treemap/scripts/extract_components.py:493-525` — where
  `package_parts` is computed and `ImportFrom` is dispatched.
- `skills/architecture-treemap/scripts/extract_components.py:458-480` — the
  `__init__`-stripping module-name derivation, root of the wrong `package_parts`.
- `skills/architecture-treemap/scripts/validate_components.py:374-457` — the
  orphan gate with the test-edge blind spot (second finding).
- `skills/architecture-treemap/scripts/metric_registry.py:101-152` — the metrics
  that inherit the skew.
- `skills/architecture-treemap/checks/golden_render.sh` — the only existing
  automated check; renders a pre-baked fixture, does not test resolution.

## Monolith constraint (where a surgical fix can sit)

`extract_components.py` is 1232 lines, sectioned by comment banners. The import
subsystem is self-contained at `L439-618` ("Python import graph"): six small
pure functions (`find_python_roots`, `module_name_for_file`,
`build_module_index`, `collect_imports`, `resolve_import_from_base`,
`import_from_candidates`, `resolve_to_repo_module`, `build_import_graph`). The
defect is isolated to how `package_parts` is derived for an `__init__.py` and
consumed by `resolve_import_from_base`. The functions are individually testable
(pure, stdlib `ast`), so a well-named resolution helper + targeted unit tests fit
inside the single file without restructuring the monolith.

## Open questions for design

- (a) **Regression-test strategy.** No extractor unit tests exist; the only
  fixture is a hand-authored render manifest
  (`references/example/components.json`) that is never extracted from source.
  Options the design must choose between: unit tests on the pure resolver
  functions (`resolve_import_from_base`, `collect_imports`,
  `resolve_to_repo_module`) with synthetic AST/inputs; OR a golden edge-set
  extracted from a committed fixture package tree (with package inits). No
  fixture source tree currently exists in the skill.
- (b) **Scope of the `--strict` test-edge blind spot** (`validate_components.py:421-428`).
  Same feature or sibling? It is the gate that masked this defect, but it is
  edge-classification, not resolution. Design to ratify in/out.
- (c) **Backward-compat of committed snapshots.** Re-extraction changes edge
  counts and metrics for every consuming repo's committed `components.json`
  (5 dated snapshots in pathway-explorer alone, under
  `docs/_meta/architecture-audit/`). Are existing snapshots frozen historical
  records (left as-is) or regenerated? How is the change surfaced (changelog,
  diff, re-synthesis run)? A human compensating `audit-asserted` edge
  (`api/__init__.py:40`) already exists and would become redundant — does it get
  removed on re-extraction?
- (d) **Multi-dot and bare relative imports.** `from ..X import Y` and
  `from . import X` (no submodule) also flow through `resolve_import_from_base`
  (`L530-536`). The witness repo has no `..`/bare-dot relative imports in its
  inits (verified by grep), so they are UNTESTED against live evidence. The
  level-walk arithmetic (`L533`) interacts with the same wrong `package_parts`
  for `__init__.py`; design must confirm correctness for `level >= 2` and for
  `node.module is None` (`L536`).
- (e) **Does the bug touch non-`__init__` evidence strings anywhere?** Verified
  NO for the witness repo: regular modules (`html_generator.py`) produce correct
  bases and matching evidence (8/8 edges verified). The mismatch between
  `evidence` line number and the imported target is unique to the collapsed
  init edges. Design should confirm this holds for repos with non-src layouts.
- (f) **Resolved base never validated against the import module index before
  candidate-trimming.** `resolve_to_repo_module` silently trims a wrong base
  down to whatever ancestor exists — there is no signal that the original
  most-specific candidate failed. Design may want resolution to fail loudly /
  record a miss rather than trim into the root (mechanism question, not a
  recommendation).

## Proposed phase decomposition (SUGGESTION — for design to ratify, not a plan)

- **Phase A — Resolver fidelity + unit tests.** Correct the `package_parts`
  derivation so a package `__init__.py` resolves single-dot relative imports
  within its own package (and confirm multi-dot / bare-dot behaviour). Add
  resolver unit tests (open question a). Confined to `extract_components.py:439-618`.
- **Phase B — `--strict` test-edge masking decision.** Decide whether the orphan
  gate (`validate_components.py:421-428`) should distinguish production edges
  from test edges so a production-orphan fails `--strict`. In-scope vs sibling
  is open question (b).
- **Phase C — Re-extraction / re-synthesis of affected snapshots.** Re-run the
  extractor over consuming repos, surface the edge-count/metric deltas, and
  decide the fate of committed snapshots and the hand-authored compensating
  edge (open question c).
