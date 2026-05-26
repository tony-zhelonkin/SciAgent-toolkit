---
date: 2026-05-26
feature: extractor-import-fidelity
phase: A
kind: plan
depends_on: [../design/decisions.md, ../map.md]
repo_root: /workspaces/DC_hum_verse/01_modules/SciAgent-toolkit
---

# Phase A — Resolver fidelity (D-1) + robustness harness (D-2)

## Objective

Make the relative-import resolver compute `package_parts` from the file path
(which knows whether a file is an `__init__.py`) rather than from the post-strip
dotted module name, so single-dot/multi-dot/bare re-exports out of a package
`__init__.py` resolve to their true sibling submodules. Simultaneously stand up
the toolkit's first extractor regression suite — a synthetic fixture tree, a
unit layer, and a golden edge-set layer — whose oracle is the CORRECT post-fix
edge set (the fix is behaviour-CHANGING, so byte-equality of an existing
snapshot cannot be the oracle).

All resolver edits sit within the self-contained import subsystem
`skills/architecture-treemap/scripts/extract_components.py:439-618`. No
restructure of the monolith.

---

## Exact steps

All anchors below are verified live in `../design/review.md` ("D-1 — Resolver
correctness [PASS]", anchors section) and pinned in `decisions.md` D-1.

### Step A1 — Add the `containing_package_parts` helper (D-1 Edit 1)

- **File:** `skills/architecture-treemap/scripts/extract_components.py`
- **Location:** insert immediately AFTER `module_name_for_file` (after L480),
  BEFORE `build_module_index` (L483).
- **Change:** add the pure helper exactly as specified in `decisions.md` D-1
  "Implementation contract", Edit 1 (L124-147). The helper returns:
  - `[]` if `source_module` is falsy;
  - `parts` (the full dotted parts) if `relative_path.endswith("/__init__.py")`
    — an `__init__` IS its package;
  - `parts[:-1]` otherwise — a regular module's package is its parent.
- **N-3 (cosmetic, optional, from `review.md`):** the inline comment about a
  top-level `__init__.py` is documenting a clause that is effectively unreachable
  (a top-level `__init__.py` yields `source_module is None`, so the
  `if not source_module: return []` guard fires first). Keep the helper exactly
  as written in `decisions.md` Edit 1; treat the top-level note as defensive
  documentation. No behavioural change either way — do NOT redesign the helper.

### Step A2 — Thread `relative_path` into `collect_imports` (D-1 Edit 2)

- **File:** `skills/architecture-treemap/scripts/extract_components.py`
- **L493 signature:**
  `def collect_imports(absolute_path, source_module):`
  → `def collect_imports(absolute_path, relative_path, source_module):`
- **L510 body:** replace
  `package_parts = source_module.split(".")[:-1] if source_module else []`
  with
  `package_parts = containing_package_parts(relative_path, source_module)`
- Do NOT touch `resolve_import_from_base` (L528-536), `import_from_candidates`
  (L539-550), or `resolve_to_repo_module` (L553-566). They are correct once fed
  a correct `package_parts` (D-1 Decision; the L533 level-walk is already correct
  for `level >= 1`).

### Step A3 — Update the sole caller (D-1 Edit 3)

- **File:** `skills/architecture-treemap/scripts/extract_components.py`
- **L585 (in `build_import_graph`):**
  `imports = collect_imports(absolute, source_module)`
  → `imports = collect_imports(absolute, relative, source_module)`
  (`relative` is already in scope at L582-585).
- **Pre-edit verification:** grep for `collect_imports(` and confirm exactly one
  caller exists (review.md confirms: def at `:493`, single call at `:585`, no
  others). If a second caller is found, STOP and flag — the design assumes one.

### Step A4 — Create the synthetic fixture tree (D-2, "Fixture tree")

- **Location:** `skills/architecture-treemap/tests/fixtures/import_tree/`
  (the `tests/` dir does not exist today — map.md L76, confirmed in review.md).
- **Tree:** create exactly the src-layout tree in `decisions.md` D-2 (L247-264):
  - `src/pkg/__init__.py` — facade: `from .alpha import A`, `from .sub import S`,
    `from . import beta` (single-dot sibling, single-dot subpackage init, bare).
  - `src/pkg/alpha.py` — regular module.
  - `src/pkg/beta.py` — regular module; `from .alpha import A` (the
    regular-module sibling control analog).
  - `src/pkg/sub/__init__.py` — `from .leaf import L` (single-dot) and
    `from ..alpha import A` (multi-dot ancestor).
  - `src/pkg/sub/leaf.py` — deepest regular module.
  - `tests/test_alpha.py` — imports `pkg.alpha` (the test→production edge reused
    by Phase B's gate test).
- Use src-layout deliberately (matches the witness repo and exercises the
  src-beats-root path in `find_python_roots`). Non-src layout is DEF-3, out of
  scope.

### Step A5 — Hand-author the golden oracle BEFORE any extractor run (D-2 + N-1)

- **File:** `skills/architecture-treemap/tests/fixtures/import_tree/expected_edges.json`
- **Content:** the seven physical-path `{from, to}` pairs in `decisions.md` D-2
  "Expected golden edge set" (L290-298):
  ```
  src/pkg/__init__.py        -> src/pkg/alpha.py
  src/pkg/__init__.py        -> src/pkg/beta.py
  src/pkg/__init__.py        -> src/pkg/sub/__init__.py
  src/pkg/beta.py            -> src/pkg/alpha.py
  src/pkg/sub/__init__.py    -> src/pkg/sub/leaf.py
  src/pkg/sub/__init__.py    -> src/pkg/alpha.py
  tests/test_alpha.py        -> src/pkg/alpha.py
  ```
- **N-1 (provenance — authoritative, from `decisions.md` L300-304 / review.md):**
  This set is hand-specified from Python import semantics against the fixture
  source BEFORE any extractor run. The extractor run only CONFIRMS it; it never
  AUTHORS it. Hand-verify each pair against the fixture source you wrote in A4.

### Step A6 — Build the unit layer (D-2 contract, "Location")

- **File:** `skills/architecture-treemap/tests/test_resolver.py`
- stdlib `unittest` or a single-file assertion script (no third-party dep — the
  extractor is stdlib-`ast` only). `importlib`-load `extract_components.py` by
  path for direct helper assertions.
- Assert the resolved dotted module path per `(relative_path, source_module,
  ImportFrom)` input for `containing_package_parts`,
  `resolve_import_from_base`, `import_from_candidates`, `resolve_to_repo_module`.
  Pin the D-1 invariant directly (init keeps full parts; regular module drops
  leaf; level-walk correct for `level>=1`; bare `module is None` path).

### Step A7 — Build the golden edge-set layer (D-2 contract)

- **File:** `skills/architecture-treemap/tests/test_import_graph_golden.py`
- Prefer subprocess invocation (`python3 extract_components.py <fixture_repo_root>`)
  reading the emitted `components.json`, to test the real entry path
  (`extract_components.py` is a script, not a package).
- Normalize emitted `edges` to a sorted set of `(from_path, to_path)` physical
  pairs plus an `evidence`-line spot-check; diff against `expected_edges.json`.

### Step A8 — Add the runner (D-2 contract, "How they run")

- **File:** `skills/architecture-treemap/tests/run.sh`
- Bash entry point that invokes both Python test files, mirroring
  `checks/golden_render.sh`'s shell-check convention so it is discoverable next
  to the existing check.

---

## Gates / acceptance-criteria checklist

Resolver (D-1):

- [ ] Edit 1 helper inserted after L480, before L483, matching `decisions.md`
      Edit 1 verbatim.
- [ ] Edit 2 signature (L493) + body (L510) swapped to use the helper.
- [ ] Edit 3 caller (L585) passes `relative`.
- [ ] Grep confirms exactly one `collect_imports(` caller (else STOP + flag).
- [ ] `resolve_import_from_base` / `import_from_candidates` /
      `resolve_to_repo_module` UNCHANGED.

Harness oracle gates (D-2, in order — the regression-detection proof):

- [ ] **(i) Oracle hand-specified BEFORE any run (N-1).** `expected_edges.json`
      authored from Python import semantics against the fixture source; each of
      the 7 pairs hand-verified. Provenance is NOT the code under test.
- [ ] **(ii) PRE-fix run FAILS the oracle.** With the fixture present but D-1
      Edits 1–3 NOT yet applied (or temporarily reverted), the golden test MUST
      fail against `expected_edges.json`. This proves the harness detects the
      regression (R5).
- [ ] **(iii) POST-fix run PASSES the oracle.** With D-1 Edits 1–3 applied, the
      golden test reproduces the exact 7-pair set and PASSES.
- [ ] **(iv) Regular-module control byte-identical.** Re-extract the witness repo
      (`/workspaces/DC_hum_verse/01_modules/pathway-explorer`) and diff its
      NON-`__init__` (regular-module) edges before vs after the fix. The
      `html_generator.py` 8/8 regular-module edges (map.md L342-346 item e) MUST
      be byte-identical — the helper returns `parts[:-1]` for regular modules,
      identical to the old `[:-1]`.

Unit layer:

- [ ] `test_resolver.py` asserts each pure helper's output, including the
      multi-dot and bare forms (retires R3 — the forms with no live witness).

Suite plumbing:

- [ ] `tests/run.sh` runs both test files and is discoverable beside
      `checks/golden_render.sh`.

---

## Dependencies

- **Upstream:** none. Phase A is self-contained (D-1 Decision: "Self-contained,
  no downstream dependency").
- **Downstream:** Phase B's gate test depends on the `import_tree` fixture from
  Step A4 (create the fixture early to unblock parallel B development). Phase C
  depends on Phase A being verified.

---

## Freeze note

This file is the executable spec for Phase A. It is applied ONLY after the
`sciagent-extension` code-change freeze is cleared AND the branching decision is
made (see `README.md` "IMPLEMENT-TIME ENTRY CONDITIONS"). No edit to
`extract_components.py` or any new `tests/` file occurs before both gates clear.
