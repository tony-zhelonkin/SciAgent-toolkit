---
date: 2026-05-26
feature: extractor-import-fidelity
status: DESIGN
kind: design
depends_on: map.md
repo_root: /workspaces/DC_hum_verse/01_modules/SciAgent-toolkit
---

# Design: extractor-import-fidelity

Architecture-first decision records for the relative-import resolution defect in
the architecture-treemap extractor. Each record is MADR-style (Context /
Decision / Implementation contract / Consequences / Risks) and is written to be
ratifiable and executable by an implementer holding `map.md`, without
re-deriving the investigation.

**Scope is fixed by two prior user decisions** (do not re-litigate):

1. The `--strict` test-edge masking fix is **in scope** for this feature — fix
   both the resolver (D-1) and the gate (D-3).
2. Backward-compat strategy is **regenerate-as-new-dated-snapshot**: after
   implement, affected repos get a NEW dated snapshot; existing snapshots are
   immutable history (D-4).

**Code-change freeze applies now.** This document is design only. The
implementation contracts below cite exact `file:line` anchors but MUST NOT be
applied until the gated implement phase.

The defect, blast radius, metrics contamination, and verified resolution traces
are in `map.md` and are not repeated here except where a decision pins to them.

---

## D-1 — Resolver correctness (the core fix)

### Context

The relative-import resolver derives the importing module's containing package
by stripping the leaf from its dotted module name:

```
package_parts = source_module.split(".")[:-1]      # collect_imports, L510
base_parts = package_parts[: len(package_parts) - (node.level - 1)]   # resolve_import_from_base, L533
```

`source_module` comes from `module_name_for_file` (L458-480), which already
strips a trailing `__init__` segment at L468-469. The verified consequence
(`map.md` "Root cause", table at L127-131):

- **Regular module** `a/b.py` → name `a.b` → `[:-1]` = `['a']` = its package. ✓
- **Package init** `a/__init__.py` → name `a` (already stripped) → `[:-1]` =
  `[]` (one segment short; the real package is `a`). ✗

So `[:-1]` is correct ONLY for regular modules. For an `__init__.py` the dotted
name **already is the containing package**, and `[:-1]` drops a real segment,
landing the relative base one level too shallow. Every single-dot re-export out
of a package `__init__.py` then mis-resolves; `resolve_to_repo_module`
(L553-566) trims the dead base down to the package root, collapsing N sibling
re-exports to one mis-targeted edge whose `evidence` line no longer matches its
import (verified trace, `map.md` L141-159).

The invariant the resolver needs: **`package_parts` must be the dotted name of
the importing module's CONTAINING PACKAGE** — for a regular module that is its
parent package; for an `__init__.py` that is the module itself.

The signal needed to distinguish the two cases (is-this-an-`__init__`) already
exists in the file at the exact site that strips it (`module_name_for_file`
L468-469). The fix threads that signal to the package derivation instead of
re-deriving the package by a leaf-strip that is wrong for inits.

### Decision

Compute `package_parts` from the **file path** (which knows whether the file is
an `__init__.py`), not from the post-strip dotted module name. Introduce one
well-named pure helper, `containing_package_parts(...)`, that states the
invariant in its docstring and is independently unit-testable. The candidate
ladder (`import_from_candidates`, `resolve_to_repo_module`) and the level-walk
arithmetic in `resolve_import_from_base` stay **unchanged** — they are correct
once fed the correct `package_parts`. This is the surgical seam the monolith
gains without a rewrite (the user accepts the monolith; see `map.md` "Monolith
constraint", L302-312).

The fix is **behaviour-CHANGING**: it ADDS correct facade edges that the buggy
resolver dropped or mis-targeted. It MUST NOT alter any currently-correct edge.
The control is the `html_generator.py` regular-module case (8/8 edges verified
correct, `map.md` L342-346, item e): those edges must be byte-identical
before/after on a re-extraction of the witness repo.

#### Relative-import forms the fix must satisfy (Python import semantics)

Let `P` = the importing module's containing package (the new invariant), as a
list of dotted parts.

| Form | `node.level` | `node.module` | Correct base |
|------|--------------|---------------|--------------|
| `from .X import Y` (sibling in current package) | 1 | `"X"` | `P + ["X"]` |
| `from . import X` (bare; submodule of current package) | 1 | `None` | `P` then candidate `P.X` |
| `from ..X import Y` (ancestor package) | 2 | `"X"` | `P[:-1] + ["X"]` |
| `from ...X import Y` | 3 | `"X"` | `P[:-2] + ["X"]` |
| `from .. import X` (bare, ancestor) | 2 | `None` | `P[:-1]` then candidate `P[:-1].X` |
| `from pkg.mod import Y` (absolute) | 0 | `"pkg.mod"` | `node.module` (untouched) |

The current level-walk `base_parts = package_parts[: len(package_parts) - (node.level - 1)]`
(L533) is **already correct for `level >= 1`** GIVEN a correct `package_parts`:
for `level == 1` it keeps all of `P`; for `level == 2` it drops one ancestor;
etc. So the fix is exclusively about supplying the correct `package_parts`. No
change to L530-536 logic is required, only to the value flowing in.

Bare `from . import X` (`node.module is None`, L536) returns `".".join(base_parts)`
= the package dotted name; `import_from_candidates` then offers `P.X` (submodule)
and bare `P` (symbol in `P/__init__.py`) — correct for a bare relative import of
a submodule. This resolves open question (d) for both multi-dot and bare forms:
**both become correct automatically once `package_parts` is right; no
arithmetic change is needed.** D-2 adds fixtures that exercise them so the claim
is enforced, not asserted.

### Implementation contract (exact edits — apply only in implement phase)

All edits in
`skills/architecture-treemap/scripts/extract_components.py`, within the
self-contained import subsystem (L439-618).

**Edit 1 — add the helper.** Insert immediately after `module_name_for_file`
(after L480), before `build_module_index` (L483):

```python
def containing_package_parts(relative_path, source_module):
    """Dotted parts of the package that CONTAINS this module.

    Invariant: relative imports in `relative_path` resolve against the
    returned package. For a regular module `a/b.py` (dotted name `a.b`) the
    containing package is `a`. For a package init `a/__init__.py` (dotted name
    `a`, already stripped of `__init__` by module_name_for_file) the module IS
    the package `a`, so the containing package is `a` itself — NOT its parent.

    `source_module` is module_name_for_file(relative_path, roots) or None.
    Returns a list of dotted parts (possibly empty for a top-level module).
    """
    if not source_module:
        return []
    parts = source_module.split(".")
    if relative_path.endswith("/__init__.py"):
        return parts            # an __init__ IS its package
    # (a top-level "__init__.py" yields source_module=None, returned empty above)
    return parts[:-1]           # a regular module's package is its parent
```

**Edit 2 — use the helper in `collect_imports`.** The current signature
(L493) is `collect_imports(absolute_path, source_module)`. Extend it to also
receive the relative path so it can detect `__init__.py`:

- L493 signature → `def collect_imports(absolute_path, relative_path, source_module):`
- L510 → replace
  `package_parts = source_module.split(".")[:-1] if source_module else []`
  with
  `package_parts = containing_package_parts(relative_path, source_module)`

**Edit 3 — update the sole caller.** In `build_import_graph` (L585):
`imports = collect_imports(absolute, source_module)` →
`imports = collect_imports(absolute, relative, source_module)`
(`relative` is already in scope at L582-585).

No other call sites exist (verify with a grep for `collect_imports(` before
applying). `resolve_import_from_base` (L528-536), `import_from_candidates`
(L539-550), `resolve_to_repo_module` (L553-566) are **unchanged**.

**Self-import guard already present.** A package `__init__.py` that re-exports a
symbol defined in its own `__init__` (bare-`base` candidate resolving to the
init file itself) is dropped by the existing `target_file == relative` guard
(L595-596). The corrected resolver does not introduce spurious self-edges.

### Consequences

- Single-dot re-exports out of every package `__init__.py` now resolve to the
  correct sibling submodules. In the witness repo the four affected inits
  (`map.md` L173-182) gain their true edge sets: `api/__init__.py` → 5 edges,
  `rendering/__init__.py` → 6, `rendering/js/__init__.py` → 9, root
  `__init__.py` → its 7 named siblings.
- `run_metadata.py` regains its production inbound edge from `api/__init__.py:41`;
  it stops rendering as a test-only/production-orphan node.
- The spurious inbound edges that collapsed into the package root `__init__.py`
  disappear (`fan_in` on the root corrects downward; `map.md` L196-199).
- Regular-module edges (the `html_generator.py` 8/8 control) are untouched: the
  helper returns `parts[:-1]` for them, identical to the old `[:-1]`.
- Metrics (`fan_in`/`fan_out`/`instability`/`refactor_pressure`) auto-recompute
  from corrected fan counts via the registry (`map.md` L282-284); no metric code
  changes. The expected metric shift is documented in D-4.

### Risks

- **R1 (path-shape assumption).** The `__init__` test keys on the path suffix
  `"/__init__.py"`. `module_name_for_file` already proves the path ends `.py`
  and uses `/` separators (L464-467), so the suffix test is consistent with the
  existing path model. Mitigation: D-2's deep-subpackage fixture
  (`rendering/js/__init__.py`-shape) exercises a nested init explicitly.
- **R2 (silent over-correction).** A wrong base no longer trims to the root, but
  `resolve_to_repo_module` can still trim a *legitimately* external relative
  import (a re-export of a third-party symbol) down to an existing ancestor.
  This is pre-existing behaviour, not introduced here; open question (f) tracks
  whether resolution should fail loudly instead of trim — see Deferred items.
- **R3 (multi-dot/bare untested against live evidence).** The witness repo has
  no `..`/bare-dot relative imports in its inits (`map.md` L335-341). The fix is
  correct by construction, but the only enforcement is D-2's synthetic fixture.
  This risk is retired by D-2, not by D-1 alone.

---

## D-2 — Robustness harness (the future-proof requirement)

### Context

There is **no extractor test suite** and **no extractable fixture tree**
(`map.md` "Tests", L72-84; open question (a), L316-323). The only fixture is a
hand-authored render manifest (`references/example/components.json`) never
extracted from source, and `checks/golden_render.sh` renders a pre-baked
manifest without exercising import resolution. The fix in D-1 is
behaviour-changing, so the regression oracle cannot be byte-equality of an
existing snapshot; it must be **"the fixture tree produces the correct edge
set."** This feature establishes the missing harness.

### Decision

Establish **two test layers** plus **one synthetic fixture package tree**:

1. **Unit layer** — tests on the pure resolution helpers
   (`containing_package_parts`, `resolve_import_from_base`,
   `import_from_candidates`, `resolve_to_repo_module`) with synthetic inputs.
   These assert the resolved dotted module path per `(relative_path,
   source_module, ImportFrom)` input. They are fast, hermetic, and pin the
   invariant in D-1 directly.
2. **Golden edge-set layer** — run the real extractor over the synthetic
   fixture tree and assert the emitted `edges` (normalized to a sorted set of
   `(from_path, to_path)` pairs, plus an `evidence`-line spot-check) equals a
   committed expected set. This is the regression oracle for the
   behaviour-changing fix: it encodes the CORRECT post-fix edges, not the
   buggy ones.

The fixture is the minimal tree that exercises **every relative-import form** in
the D-1 table, so the multi-dot/bare risk (R3) is enforced by a passing test.

#### Fixture tree (synthetic, committed)

Create under `skills/architecture-treemap/tests/fixtures/import_tree/` a tiny
package using a `src/` layout (mirrors the witness repo so `find_python_roots`
exercises the src-beats-root path):

```
tests/fixtures/import_tree/
  src/
    pkg/
      __init__.py        # package init: sibling re-exports (facade), e.g.
                         #   from .alpha import A          (single-dot sibling)
                         #   from .sub import S            (single-dot subpackage init)
                         #   from . import beta            (bare: submodule of pkg)
      alpha.py           # regular module
      beta.py            # regular module; imports a sibling: from .alpha import A
                         #   (regular-module sibling import — the html_generator control analog)
      sub/
        __init__.py      # nested subpackage init: from .leaf import L (single-dot)
                         #   and a multi-dot ancestor import: from ..alpha import A
        leaf.py          # deepest regular module
  tests/
    test_alpha.py        # imports pkg.alpha (a test->production edge, for D-3)
```

This tree covers, by construction:

- regular-module sibling import (`beta.py: from .alpha`) — the regular-module
  control (`package_parts == ['pkg']`, unchanged behaviour).
- package-`__init__` sibling re-exports (`pkg/__init__.py`) — the api-facade
  analog (must produce `pkg/__init__.py -> alpha.py`, `-> sub/__init__.py`,
  `-> beta.py`).
- a nested subpackage (`sub/__init__.py`) with a single-dot import to `leaf.py`.
- a **multi-dot** import (`sub/__init__.py: from ..alpha`) resolving up to
  `pkg/alpha.py`.
- a **bare** import (`pkg/__init__.py: from . import beta`) resolving to
  `pkg/beta.py`.
- a **test→production** edge (`tests/test_alpha.py -> pkg/alpha.py`) reused by
  D-3's gate test.

#### Expected golden edge set (the oracle)

Committed alongside the fixture as
`tests/fixtures/import_tree/expected_edges.json` — a sorted list of
`{"from": <rel-path>, "to": <rel-path>}` pairs at the **physical-path** level
(id-independent, so it survives id-scheme changes). Authored to the
CORRECT post-fix resolution:

```
src/pkg/__init__.py        -> src/pkg/alpha.py
src/pkg/__init__.py        -> src/pkg/beta.py
src/pkg/__init__.py        -> src/pkg/sub/__init__.py
src/pkg/beta.py            -> src/pkg/alpha.py
src/pkg/sub/__init__.py    -> src/pkg/sub/leaf.py
src/pkg/sub/__init__.py    -> src/pkg/alpha.py
tests/test_alpha.py        -> src/pkg/alpha.py
```

(The pairs above are the authoritative oracle, hand-specified from Python
import semantics against the fixture source **before any extractor run** — they
are NOT authored by running the code under test. The implement-time procedure
only CONFIRMS them: a post-fix run must reproduce this exact set, and a pre-fix
run of the same harness MUST fail it — that is the regression-detection proof.)

### Implementation contract

- **Location.** New directory `skills/architecture-treemap/tests/` (does not
  exist today, `map.md` L76). Layout:
  - `tests/test_resolver.py` — unit layer (stdlib `unittest` or a single-file
    assertion script; no third-party dep — the extractor is stdlib-`ast` only).
  - `tests/test_import_graph_golden.py` — golden edge-set layer; imports
    `extract_components` functions or shells out to `extract_components.py`
    against the fixture tree, then diffs emitted edges to
    `expected_edges.json`.
  - `tests/fixtures/import_tree/...` — the synthetic tree above.
  - `tests/fixtures/import_tree/expected_edges.json` — the oracle.
- **How they run.** Add a runner `skills/architecture-treemap/tests/run.sh`
  that invokes both Python test files (mirroring `checks/golden_render.sh`'s
  shell-check convention so it is discoverable next to the existing check). The
  unit + golden tests are Python; the runner is the bash entry point. The
  fixture's `src/` and `tests/` paths feed straight into the extractor's
  existing `find_python_roots` / `build_module_index` machinery — no extractor
  plumbing changes to make the fixture extractable.
- **Importing the extractor functions.** `extract_components.py` is a script,
  not a package. The golden test should prefer subprocess invocation
  (`python3 extract_components.py <fixture_repo_root>`) reading the emitted
  `components.json`, to test the real entry path; the unit test may
  `importlib`-load the module file by path for direct helper assertions.

### Consequences

- The toolkit gains its first extractor regression suite, scoped to import
  resolution — the seam D-1 created becomes continuously enforced.
- The multi-dot/bare-form correctness claim (D-1, open question (d)) is enforced
  by a test, retiring R3.
- The fixture doubles as D-3's gate fixture (the `test_alpha.py` edge), avoiding
  a second fixture tree (decoupled but not duplicated).

### Risks

- **R4 (fixture drift from real layouts).** A synthetic tree may not capture an
  exotic real layout (namespace packages, non-src roots). Mitigation: the tree
  deliberately uses src-layout to match the witness repo and the only verified
  failure mode; non-src layouts are tracked as open question (e) — see Deferred
  items. The golden test is additive: real-repo layouts can be added as further
  fixtures without redesign.
- **R5 (oracle authored wrong).** The expected edge set is hand-authored. If a
  pair is wrong the test enforces a wrong invariant. Mitigation: implement-time
  procedure requires hand-verifying each pair against fixture source AND
  confirming a *pre-fix* extractor run fails the same oracle (proving the test
  detects the regression).

---

## D-3 — `--strict` gate hardening

### Context

`validate_referential_integrity` (`validate_components.py:374-457`) builds
`logical_nodes_with_edges` (L421-428) by marking a node connected if it appears
as `from` OR `to` in ANY logical-id edge, with **no filter on the edge's source
kind**. The rollup attaches `test -> node` edges as logical edges (the test
suite is the `test-suite` logical node). So a node whose ONLY logical edge is
`test-suite -> node` counts as connected and passes `--strict` (L536-537 exits 1
only on `ref_issues`).

Verified against the witness snapshot: `run-metadata` (classification `seam`,
not `removable`, not `leaf`) has exactly one logical edge,
`test-suite -> run-metadata`, and therefore passes the orphan gate today despite
having no production importer. This is the blind spot that let the D-1 defect
render without failing validation (`map.md` "second finding", L258-269; open
question (b), resolved IN SCOPE per user decision 1).

This is edge-CLASSIFICATION, a sibling concern to D-1's edge-RESOLUTION; it is
in scope by decision, fixed independently in `validate_components.py`.

#### How a "test edge" is identified from the manifest

Verified structurally: the test suite is a single logical node `test-suite`
whose `physical_files` are **all** under `tests/` (witness: 57/57 paths start
`tests/`). The manifest-driven signal, in order of preference:

1. **Primary:** the `from` endpoint is a logical node ALL of whose
   `physical_files[].path` start with a test directory segment (`tests/` or a
   path component named `test`/`tests`). This is the robust definition: a "test
   edge" is one whose source logical node is a test-only owner.
2. **Fallback (physical edges):** for a physical-id `from` endpoint, the
   physical component's `path` matches a test pattern (`tests/` prefix or
   `test_*.py` / `*_test.py` basename). Reused from the same predicate the
   coverage computation already implies (`compute_module_coverage`,
   `map.md` L208-210).

A single shared helper `is_test_source_node(node_id, data)` encodes this so the
definition lives in one place.

### Decision

Redefine the orphan check so a logical node whose **only** logical edges are
test→node edges is treated as a **production-orphan**, NOT as connected.
Introduce a distinct severity: **`production-orphan` is a hard error under
`--strict`** (joins `ref_issues`), with a dedicated message that names the
masking test edge. Rationale: this defect class silently shipped a broken
treemap; a hard error under `--strict` is the gate's purpose. Outside `--strict`
it surfaces as a soft warning (consistent with how non-strict runs never fail).

Preserve the existing `removable` / `leaf:true` exemptions exactly: a node that
is `removable` or `leaf` and test-only-connected stays a soft warning (those
nodes are tolerated standalone by design, L441-449). Do NOT regress a node that
has ANY production logical edge — production connectivity is computed by
EXCLUDING test-source edges, so any single production edge keeps a node
connected.

### Implementation contract (exact edits — apply only in implement phase)

All edits in `validate_components.py`.

**Edit 1 — add the test-source predicate.** Add near the top of
`validate_referential_integrity` (or as a module-level helper above L374):

```python
def is_test_source_node(node_id, logical_meta):
    """True if `node_id` is a logical node whose physical files are ALL tests.

    A test edge is one whose `from` endpoint is such a node; edges from it must
    NOT count toward production connectivity in the orphan gate.
    """
    lc = logical_meta.get(node_id)
    if not lc:
        return False
    paths = [pf.get("path", "") for pf in lc.get("physical_files", [])]
    if not paths:
        return False
    return all(p.startswith("tests/") or "/tests/" in p or
               p.split("/")[-1].startswith("test_") for p in paths)
```

(`logical_meta` is already built at L431-435; move its construction above the
connectivity loop, or pass it in.)

**Edit 2 — split connectivity into production vs any.** Replace the single
`logical_nodes_with_edges` accumulation (L421-428) with two sets:

```python
nodes_with_any_edge: set[str] = set()
nodes_with_production_edge: set[str] = set()
for edge in data.get("edges", []):
    frm = edge.get("from", "")
    to  = edge.get("to",   "")
    from_is_test = is_test_source_node(frm, logical_meta)
    if frm in logical_ids:
        nodes_with_any_edge.add(frm)
    if to in logical_ids:
        nodes_with_any_edge.add(to)
        # `to` is connected to PRODUCTION only if the source is not a test node.
        if not from_is_test:
            nodes_with_production_edge.add(to)
    if frm in logical_ids and not from_is_test:
        nodes_with_production_edge.add(frm)
```

Note: an edge FROM a production node TO a test node still makes the production
node production-connected (it has a real outbound edge); the test node's own
connectivity is governed by its `from`-edges, which are test-sourced. The
asymmetry is intentional and matches Python import direction.

**Edit 3 — orphan classification.** In the per-node loop (L437-455), replace
the `lc_id in logical_nodes_with_edges` test with the two-tier check:

```python
for lc_id in sorted(logical_ids):
    lc = logical_meta.get(lc_id, {})
    is_removable = lc.get("classification") == "removable"
    is_leaf = bool(lc.get("leaf"))
    if lc_id in nodes_with_production_edge:
        continue  # genuinely connected to production
    if is_removable or is_leaf:
        soft_warnings.append(... existing tolerated-orphan message ...)
        continue
    if lc_id in nodes_with_any_edge:
        # connected ONLY by test edges -> production-orphan
        issues.append(
            f"logical_components '{lc_id}': production-orphan — its only logical "
            "edges originate from the test suite; no production component imports "
            "it. Author a production edge or reclassify (removable/leaf)."
        )
    else:
        issues.append(... existing non-leaf-orphan message ...)
```

`--strict` (L536-537) needs **no change** — it already exits 1 on any
`ref_issues`, and `production-orphan` is appended there.

**Edit 4 — gate test.** Add to D-2's harness a
`tests/test_strict_gate.py` (or a section in the golden test) that:

- Builds a manifest where node `X` has only a `test-suite -> X` edge and is
  `classification: seam` (not removable/leaf): asserts `--strict` returns a
  `production-orphan` issue and exit 1.
- Builds the same manifest but with a production edge `Y -> X` added: asserts
  no issue (no regression of legitimately-connected nodes).
- Builds a manifest where `X` is `leaf: true` and test-only: asserts a soft
  warning, not an issue (exemption preserved).

The fixture for this reuses D-2's `import_tree` (the `test_alpha.py -> pkg.alpha`
edge), so the gate test runs against a manifest extracted from the same tree.

### Consequences

- A production-orphan created by the D-1 defect (or any future regression) now
  fails `--strict` instead of rendering silently. After D-1 lands, the witness
  repo's `run-metadata` will be production-connected (via the restored facade
  edge) and will pass — so D-3 does not spuriously fail the corrected repo.
- The `removable`/`leaf` exemptions are unchanged; the only behaviour change is
  for non-exempt nodes connected solely by test edges.

### Risks

- **R6 (test-node identification too narrow/broad).** If a repo splits its tests
  across multiple logical nodes, or names a production node with a `test_`
  prefix, `is_test_source_node` could misclassify. Mitigation: the predicate
  requires ALL physical files of the node to be test files (a production node
  with one `test_`-named file is NOT a test node), and the witness repo's single
  `test-suite` node validates the common case. Multi-test-node repos are a
  Deferred item with a recommended default (extend the predicate, do not block).
- **R7 (interaction with D-1).** D-3 must land in a state where the corrected
  extractor (D-1) does not leave the witness repo's `run-metadata` as a
  production-orphan. Phase ordering (Phase A before re-running validation in
  Phase C) guarantees this; see phase decomposition.

---

## D-4 — Backward-compat + snapshot regeneration

### Context

Re-extraction changes edge counts and metrics for every consuming repo's
committed `components.json`. In the witness repo alone there are 5 dated
snapshots under `docs/_meta/architecture-audit/` (`2026-05-21`, `2026-05-24`,
`2026-05-24-post-refactor`, `2026-05-25`, `2026-05-25-postrefactor`;
verified on disk). These are excluded from re-extraction *as inputs* by
`is_audit_output` (`extract_components.py:333-350`); a fresh run produces a NEW
snapshot, it does not rewrite them (`map.md` L58-69; open question (c)).

A human authored a compensating `audit-asserted` edge to patch ONE of the five
lost facade edges:
`api/__init__.py:40 -> from .channels import ...` (verified in the snapshot as
`public-api-surface -> channel-owners`, `evidence_class: audit-asserted`). Once
D-1 restores the structural edge, this hand-authored edge is **redundant**.

### Decision

Per user decision 2: **regenerate-as-new-dated-snapshot.**

1. After D-1+D-3 land (and are verified by D-2), regenerate the affected repos
   by running the corrected extractor → a NEW dated snapshot directory (e.g.
   `docs/_meta/architecture-audit/2026-MM-DD-import-fidelity/`). The 5 existing
   snapshots are **immutable history** — not rewritten, not deleted.
2. In the regenerated snapshot, the hand-authored compensating `audit-asserted`
   edge at `api/__init__.py:40` is **dropped** — the structural facade edge now
   exists, so the compensation is redundant and keeping it would double-count.
   (The OTHER `audit-asserted` edges — `link-outs -> metadata-assembler`,
   `link-outs -> render-core` — are NOT structural import edges; they stay.
   Only the one that compensated for the resolver bug is dropped. Verify each
   `audit-asserted` edge's evidence against a restored structural edge before
   dropping; drop only those now covered by a real static edge.)
3. Document the intended metrics shift so reviewers expect the diff rather than
   flag it as a regression (below).

#### Intended metrics shift (expected, not a regression)

For package `__init__.py` targets and their re-export targets:

- **`fan_out`** of each package init **increases** to its true re-export count
  (`api/__init__.py` 1→5, `rendering/__init__.py` 1→6,
  `rendering/js/__init__.py` 1→9, root init →7). (`map.md` L191-194.)
- **`fan_in`** of the package ROOT `__init__.py` **decreases** — the spurious
  inbound edges from sibling inits (witness: `fan_in: 2`) disappear.
  (`map.md` L196-199.)
- **`fan_in`** of re-export TARGETS (e.g. `run_metadata.py`) **increases** —
  they regain their true production importers.
- **`instability`** (`fan_out/(fan_in+fan_out)`) recomputes for all the above;
  `run_metadata.py` moves off its bogus test-only `fan_in: 1`.
- **`refactor_pressure`** (uses `max(fan_in,1)` as a multiplier) **increases**
  for facade targets whose fan_in was under-counted.
- **edge count** rises (collapsed 1-edges expand to N); **`evidence` strings**
  on the previously-collapsed init edges now match their source lines.

These are all *corrections*, surfaced via the new snapshot diff. No metric code
changes (registry recomputes from corrected fan counts, `map.md` L282-284).

### Implementation contract

- This is a **process/regeneration** step, not a source edit to the toolkit.
  Executed in Phase C (post-implement, post-verify).
- For each affected consuming repo (pathway-explorer is the verified witness;
  others discovered by searching for `docs/_meta/architecture-audit/`):
  1. Run the corrected `extract_components.py` over the repo → new dated
     snapshot dir.
  2. Run `validate_components.py --strict` on the new snapshot (D-3) — it must
     pass (the restored facade edge resolves the `run-metadata`
     production-orphan).
  3. Run `rollup_logical_edges.py` and `render_treemap.py` to refresh derived
     edges and HTML in the new snapshot.
  4. **N-2 orphan guard — confirm before dropping.** `channel-owners`' ONLY
     edge in the current snapshot is the `audit-asserted` compensating edge
     about to be dropped. After step 3's rollup, verify `channel-owners` now
     carries a *derived static* edge (the restored facade import). Only if that
     structural edge is present may the compensating edge be removed — otherwise
     dropping it orphans `channel-owners` and `--strict` (D-3) will fail. Then
     remove the now-redundant `audit-asserted` compensating edge (item 2 above)
     from the regenerated manifest only.
  5. Commit the new snapshot dir; leave prior snapshots untouched.
- Surface the change: a short changelog/diff note in the regenerated snapshot
  dir (or the repo's audit README) listing the edge-count and metrics deltas, so
  reviewers reading the diff understand it is the fidelity correction.

### Consequences

- History is preserved (prior snapshots immutable); the corrected graph lives in
  a clearly-dated new snapshot.
- The hand-authored compensation is retired at its proper layer (the resolver),
  removing a manual workaround.
- Reviewers have a documented expectation for the metrics movement.

### Risks

- **R8 (other consuming repos unknown).** Only pathway-explorer is verified.
  Mitigation: Phase C begins with a search for all
  `docs/_meta/architecture-audit/` dirs across consuming repos; regenerate each.
- **R9 (dropping the wrong audit-asserted edge).** Mitigation: drop only edges
  whose evidence now maps to a confirmed restored static edge; verify by diffing
  the new static edge set against the audit-asserted set before removal.

---

## Resolution of the map's open questions

- **(a) Regression-test strategy** → **D-2.** Both layers: unit tests on the
  pure helpers AND a golden edge-set extracted from a committed synthetic
  fixture tree (the oracle is the correct post-fix edge set, not byte-equality).
- **(b) `--strict` test-edge blind spot scope** → **in scope** per user decision
  1; designed in **D-3** (hard `production-orphan` error under `--strict`).
- **(c) Backward-compat of committed snapshots** → **D-4.** Regenerate as new
  dated snapshots; old ones immutable; redundant compensating edge dropped on
  regeneration.
- **(d) Multi-dot and bare relative imports** → **D-1.** Both become correct
  automatically once `package_parts` is the correct containing package; the
  level-walk arithmetic (L533) is unchanged and already correct for `level>=1`;
  bare `from . import X` (L536) is correct via the candidate ladder. Enforced
  by D-2's fixture (`from ..alpha`, `from . import beta`).
- **(e) Non-`__init__` evidence-string effects** → confirmed NO for the witness
  repo (regular modules produce correct bases, 8/8 control). The fix preserves
  this (helper returns `parts[:-1]` for regular modules, identical to old).
  Confirming the claim holds for **non-src layouts** is the one residual: the
  helper keys on the `/__init__.py` path suffix, which is layout-independent, so
  the fix is layout-agnostic by construction — but no non-src fixture is
  verified. Recommended default: D-2's fixture uses src-layout (matching the
  only verified case); add a non-src fixture if/when a non-src consuming repo
  appears (Deferred).
- **(f) Resolved base never validated before candidate-trimming** → **partly
  deferred.** D-1 makes the base correct, so the silent-trim-to-root pathology
  no longer fires for the known defect. Whether resolution should *fail loudly*
  on a base that is not in the module index (rather than trim) is a separate
  hardening (R2). Not resolved in this design; see Deferred items.

### Deferred items (with recommended defaults)

- **DEF-1 (open question f): loud-fail on unresolved base.** Recommended
  default: keep the trim-to-ancestor behaviour (it correctly handles
  symbol-vs-submodule and third-party re-exports), but add an optional
  diagnostic count of "relative imports that trimmed past their package" so a
  future run can surface suspicious trims without changing edge output. Do not
  block this feature on it.
- **DEF-2 (R6): multi-test-node repos.** Recommended default: keep the
  all-physical-files-are-tests predicate; extend to a set of test nodes only
  when a real repo needs it.
- **DEF-3 (open question e / R4): non-src layout fixture.** Recommended default:
  src-layout fixture now (matches verified case); add non-src fixture when a
  non-src consuming repo is encountered. The fix is layout-agnostic by
  construction so this is a coverage gap, not a correctness gap.

---

## Proposed phase decomposition

**SUGGESTION for the /plan stage and architect to ratify — not a committed plan.**

- **Phase A — Resolver fidelity + D-2 harness.** Apply D-1 (the
  `containing_package_parts` helper + the two call-site edits in
  `extract_components.py`, all within L439-618). Build the D-2 fixture tree, the
  unit layer, and the golden edge-set layer; confirm a pre-fix run fails the
  oracle and the post-fix run passes. Confined to
  `extract_components.py` + new `tests/`. Self-contained, no downstream
  dependency.

- **Phase B — `--strict` gate hardening.** Apply D-3 (the `is_test_source_node`
  predicate + the two-tier connectivity split + `production-orphan` severity in
  `validate_components.py`) and add the gate test (reusing Phase A's fixture).
  Depends on Phase A only for the shared fixture, not for correctness.

- **Phase C — Snapshot regeneration.** Apply D-4: discover all consuming repos
  with `docs/_meta/architecture-audit/` snapshots, regenerate each as a new
  dated snapshot with the corrected extractor, drop the redundant
  `audit-asserted` compensating edge (only after confirming the restored
  structural edge keeps `channel-owners` connected — N-2), re-validate under
  `--strict`, refresh derived edges/HTML, and write the metrics-shift changelog. Depends on Phases A
  and B landing and being verified. Process step, no toolkit source edits.

Ordering rationale: A must precede C (C re-runs the corrected extractor); B must
precede C (C re-validates under the hardened `--strict`, which must pass for the
corrected repo). A and B are independent except for the shared D-2 fixture and
may be developed in parallel if the fixture is created first.
