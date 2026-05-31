---
date: 2026-05-26
feature: extractor-import-fidelity
phase: B
kind: plan
depends_on: [../design/decisions.md, ../map.md, phase-A.md]
repo_root: /workspaces/DC_hum_verse/01_modules/SciAgent-toolkit
---

# Phase B — `--strict` gate hardening (D-3)

## Objective

Harden the `--strict` orphan gate so a logical node whose ONLY logical edges are
test→node edges is treated as a **production-orphan** (a hard error under
`--strict`, exit 1) rather than counting as "connected." Preserve the existing
`removable` / `leaf:true` exemptions exactly, and never regress a node that has
ANY production logical edge. This closes the blind spot that let the D-1 defect
render `run-metadata` without failing validation.

All edits in `skills/architecture-treemap/scripts/validate_components.py`. This
is edge-CLASSIFICATION, a sibling concern to D-1's edge-RESOLUTION, fixed
independently.

---

## Exact steps

Anchors verified live in `../design/review.md` ("D-3 — `--strict` gate hardening
[PASS]", anchors section) and pinned in `decisions.md` D-3.

### Step B1 — Add the test-source predicate (D-3 Edit 1)

- **File:** `validate_components.py`
- **Location:** module-level helper above `validate_referential_integrity`
  (L374), or near the top of that function.
- **Change:** add `is_test_source_node(node_id, logical_meta)` exactly as in
  `decisions.md` D-3 Edit 1 (L421-435). It returns True iff the node has at least
  one physical file and ALL its `physical_files[].path` are test files
  (`startswith("tests/")` or `"/tests/" in p` or basename `startswith("test_")`).
  The `all()` quantifier is load-bearing (R6 mitigation): a production node with
  a single `test_`-named file is NOT a test node.

### Step B2 — Move `logical_meta` construction up (D-3 Edit 1 note)

- **File:** `validate_components.py`
- `logical_meta` is currently built at L431-435, AFTER the connectivity loop.
  Move its construction ABOVE the connectivity loop (or pass it in) so both the
  predicate (B1) and the split (B3) can use it. Verified necessary in review.md.

### Step B3 — Split connectivity into production vs any (D-3 Edit 2)

- **File:** `validate_components.py`
- Replace the single `logical_nodes_with_edges` accumulation (L421-428) with the
  two-set loop from `decisions.md` D-3 Edit 2 (L444-460):
  - `nodes_with_any_edge` — `from`/`to` in any logical edge (old behaviour).
  - `nodes_with_production_edge` — `to` is production-connected only if the
    `from` endpoint is NOT a test-source node; a `from` endpoint is
    production-connected if it is itself not a test-source node.
- The asymmetry is intentional (`decisions.md` L462-465): an edge FROM a
  production node TO a test node still makes the production node
  production-connected; the test node's own connectivity is governed by its
  `from`-edges, which are test-sourced.

### Step B4 — Two-tier orphan classification (D-3 Edit 3)

- **File:** `validate_components.py`
- In the per-node loop (L437-455), replace the
  `lc_id in logical_nodes_with_edges` test (L438) with the two-tier check from
  `decisions.md` D-3 Edit 3 (L471-489):
  - `lc_id in nodes_with_production_edge` → `continue` (genuinely connected).
  - `is_removable or is_leaf` → soft warning (existing tolerated-orphan message;
    exemption preserved exactly).
  - else `lc_id in nodes_with_any_edge` → append a **`production-orphan`** issue
    with the dedicated message naming the masking test edge.
  - else → existing non-leaf-orphan message.
- `--strict` (L536-537) needs NO change — it already exits 1 on any `ref_issues`,
  and `production-orphan` is appended to `issues` (verified: `issues` is the
  first returned tuple element at L457 and becomes `ref_issues`).

### Step B5 — Gate test (D-3 Edit 4)

- **File:** `skills/architecture-treemap/tests/test_strict_gate.py` (or a section
  in the golden test), added to Phase A's harness.
- **Fixture:** reuse Phase A's `import_tree` (the `tests/test_alpha.py ->
  pkg.alpha` edge), so the gate test runs against a manifest extracted from the
  same tree. This is the ONLY dependency on Phase A.

---

## Gates / acceptance-criteria checklist

Edits:

- [ ] `is_test_source_node` added (D-3 Edit 1); `all()`-quantified over the
      node's own `physical_files`.
- [ ] `logical_meta` moved above the connectivity loop (D-3 Edit 1 note).
- [ ] Connectivity split into `nodes_with_any_edge` + `nodes_with_production_edge`
      (D-3 Edit 2), with the intentional from/to asymmetry.
- [ ] Per-node loop uses the two-tier check; `production-orphan` appended to
      `issues`; `removable`/`leaf` soft-warning path unchanged (D-3 Edit 3).
- [ ] `--strict` (L536-537) UNCHANGED.

Three gate-test cases (D-3 Edit 4, L494-506 — each must be concretely asserted):

- [ ] **Case 1: test-only seam node.** Manifest where node `X` has only a
      `test-suite -> X` edge and `classification: seam` (not removable/leaf):
      `--strict` returns a `production-orphan` issue AND exits 1.
- [ ] **Case 2: production edge added.** Same manifest with a production edge
      `Y -> X` added: NO issue (no regression of legitimately-connected nodes).
- [ ] **Case 3: leaf:true test-only.** Manifest where `X` is `leaf: true` and
      test-only: a SOFT warning, NOT an issue (exemption preserved).

---

## Dependencies

- **Upstream:** Phase A ONLY for the shared `import_tree` fixture (Step A4), NOT
  for D-1 resolver correctness (`decisions.md` D-3 Decision / "Proposed phase
  decomposition": "Depends on Phase A only for the shared fixture, not for
  correctness"). Once the fixture exists, B can be developed in parallel with the
  rest of Phase A.
- **Downstream:** Phase C re-validates regenerated snapshots under this hardened
  gate; B must land and verify before C.

---

## Freeze note

This file is the executable spec for Phase B. It is applied ONLY after the
`sciagent-extension` code-change freeze is cleared AND the branching decision is
made (see `README.md` "IMPLEMENT-TIME ENTRY CONDITIONS"). No edit to
`validate_components.py` or any test file occurs before both gates clear.
