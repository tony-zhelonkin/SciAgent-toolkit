---
date: 2026-05-26
feature: extractor-import-fidelity
kind: plan
status: PLAN
depends_on: [../design/decisions.md, ../design/review.md, ../map.md]
repo_root: /workspaces/DC_hum_verse/01_modules/SciAgent-toolkit
---

# Plan: extractor-import-fidelity

This plan SEQUENCES and DETAILS the architect-PASSED design in
`../design/decisions.md` (gate verdict in `../design/review.md`:
PASS-WITH-NITS; nits N-1/N-2/N-3 already folded into `decisions.md`). It
introduces **no new decisions**. Every step cites the authoritative
`decisions.md` section and the live `file:line` anchors verified in the review.

An implementer holding `decisions.md` can pick up one phase file and execute it
end-to-end without re-deriving the investigation.

---

## The three phases

| Phase | File | Scope | Touches code? |
|-------|------|-------|---------------|
| **A** | `phase-A.md` | Resolver fidelity (D-1) + robustness harness (D-2) | Yes — `extract_components.py` + new `tests/` |
| **B** | `phase-B.md` | `--strict` gate hardening (D-3) | Yes — `validate_components.py` + gate test |
| **C** | `phase-C.md` | Snapshot regeneration (D-4) | No toolkit source — process/regeneration step in consuming repos |

---

## Ordering rationale (ratified from `decisions.md` "## Proposed phase decomposition")

- **A must precede C.** Phase C re-runs the corrected extractor over consuming
  repos; the correction is D-1, which lands in Phase A.
- **B must precede C.** Phase C re-validates each regenerated snapshot under the
  hardened `--strict`; that hardening is D-3, which lands in Phase B. The
  corrected repo must PASS the hardened gate (D-1 restores the facade edge that
  un-orphans `run-metadata`), so B's logic must exist before C validates.
- **A and B are parallel-able once the D-2 fixture exists.** D-3's gate test
  (D-3 Edit 4) reuses Phase A's `import_tree` fixture (the
  `tests/test_alpha.py -> pkg/alpha.py` edge). The dependency is the **fixture
  only**, not D-1's resolver correctness. Once Phase A has created
  `tests/fixtures/import_tree/`, Phase B can be developed in parallel with the
  rest of Phase A.

```
   create D-2 fixture tree (early in Phase A)
            |
   +--------+--------+
   |                 |
 Phase A          Phase B          (parallel-able once fixture exists)
 (D-1 + D-2)      (D-3)
   |                 |
   +--------+--------+
            |
        verify A and B
            |
         Phase C  (D-4 regeneration; depends on A AND B verified)
```

---

## IMPLEMENT-TIME ENTRY CONDITIONS (gate Phase 4 / the /implement stage)

No code is touched until BOTH of the following gates clear:

1. **Code-change freeze cleared.** The `sciagent-extension` code-change freeze is
   in force now (`decisions.md` "Code-change freeze applies now"). This plan,
   and the design it ratifies, are documentation only. Implementation of any
   phase is the **gated next phase** and begins only after the freeze is lifted.

2. **Branching decision made.** The toolkit is currently on
   `feat/sciagent-extension`. The import-fidelity fix is an unrelated concern and
   likely wants its **own branch** (e.g. `fix/extractor-import-fidelity`) cut
   from the appropriate base, not commingled with the extension work on
   `feat/sciagent-extension`. The implementer must make and record this branching
   decision before applying any edit. (Parallel-agent / submodule-branch
   awareness applies: verify the submodule branch and topology before staging.)

Stated plainly: **no edits to `extract_components.py`, `validate_components.py`,
any new `tests/` file, or any consuming-repo snapshot occur until (1) the freeze
is cleared and (2) the branch is chosen.** Each phase file repeats this as a
short "freeze note."

---

## Out of scope for this feature (deferred — do NOT implement)

Carried verbatim from `decisions.md` "### Deferred items (with recommended
defaults)". The implementer must NOT pick these up while executing Phases A–C:

- **DEF-1 (open question f): loud-fail on unresolved base.** Keep the
  trim-to-ancestor behaviour. An optional diagnostic count of "relative imports
  that trimmed past their package" is the recommended default for a FUTURE run —
  do not block this feature on it, and do not change edge output.
- **DEF-2 (R6): multi-test-node repos.** Keep the all-physical-files-are-tests
  predicate (`is_test_source_node`). Extend to a set of test nodes only when a
  real repo needs it. Not in scope now.
- **DEF-3 (open question e / R4): non-src layout fixture.** Use the src-layout
  fixture now (matches the only verified case). The fix is layout-agnostic by
  construction; a non-src fixture is a coverage gap, not a correctness gap. Add
  one only when a non-src consuming repo appears. Not in scope now.

---

## Traceability

Each phase file's gates trace to a specific `decisions.md` decision and anchor:

- Phase A gates → D-1 (Edits 1–3, L118-166) and D-2 (fixture L242-304, oracle
  L283-304, contract L306-329), plus N-1 (oracle provenance).
- Phase B gates → D-3 (Edits 1–4, L414-507), including the three gate-test cases
  in D-3 Edit 4 (L494-506).
- Phase C gates → D-4 (contract L592-616), including the N-2 channel-owners
  guard (L605-612) and the intended metrics shift (L570-590).
