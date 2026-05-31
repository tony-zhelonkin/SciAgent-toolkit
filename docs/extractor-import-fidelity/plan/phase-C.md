---
date: 2026-05-26
feature: extractor-import-fidelity
phase: C
kind: plan
depends_on: [../design/decisions.md, ../map.md, phase-A.md, phase-B.md]
repo_root: /workspaces/DC_hum_verse/01_modules/SciAgent-toolkit
---

# Phase C — Snapshot regeneration (D-4)

## Objective

Regenerate the affected consuming repos with the corrected extractor as NEW
dated snapshots (existing snapshots are immutable history), drop the now-redundant
hand-authored `audit-asserted` compensating edge ONLY after confirming the
restored structural edge keeps the node connected, re-validate each new snapshot
under the hardened `--strict`, refresh derived edges/HTML, and write a
metrics-shift changelog so reviewers expect the diff rather than flag it as a
regression.

This is a **process / regeneration step in the consuming repos — NO toolkit
source edits.** Per user decision 2: regenerate-as-new-dated-snapshot.

---

## Exact steps

From `decisions.md` D-4 "Implementation contract" (L592-616).

### Step C1 — Discover all affected consuming repos (D-4 contract item, R8)

- Search for `docs/_meta/architecture-audit/` directories across consuming
  repos. The verified witness is `pathway-explorer`
  (`/workspaces/DC_hum_verse/01_modules/pathway-explorer`); others are discovered
  here, not assumed. Regenerate EACH discovered repo.

### Step C2 — Re-run the corrected extractor → new dated snapshot (D-4 item 1)

- For each affected repo, run the corrected `extract_components.py` (Phase A) over
  the repo → a NEW dated snapshot directory, e.g.
  `docs/_meta/architecture-audit/2026-MM-DD-import-fidelity/`.
- The 5 existing snapshots in the witness
  (`2026-05-21`, `2026-05-24`, `2026-05-24-post-refactor`, `2026-05-25`,
  `2026-05-25-postrefactor`) are immutable — NOT rewritten, NOT deleted. They are
  excluded as inputs by `is_audit_output` (`extract_components.py:333-350`).

### Step C3 — Refresh derived edges and HTML (D-4 item 3)

- Run `rollup_logical_edges.py` then `render_treemap.py` against the new snapshot
  to refresh derived logical edges and the rendered HTML.

### Step C4 — N-2 orphan guard: confirm BEFORE dropping (D-4 item 4 / N-2)

This is the load-bearing sequencing gate (review.md N-2). The node
`channel-owners` has EXACTLY ONE edge in the current snapshot — the
`audit-asserted` compensating edge `public-api-surface -> channel-owners` that
D-4 is about to drop; it receives no test edge either.

- After C3's rollup, **VERIFY** that `channel-owners` now carries a *derived
  static* edge — the restored facade import `api/__init__.py:39 .channels`
  rolling up to a static `public-api-surface -> channel-owners` logical edge.
- **Only if that structural edge is present** may the redundant `audit-asserted`
  compensating edge be removed from the regenerated manifest. Dropping it WITHOUT
  the structural edge present orphans `channel-owners`, and the hardened
  `--strict` (Phase B) will (correctly) fail.
- Drop ONLY the one compensating edge that D-1 now covers
  (`public-api-surface -> channel-owners`, `evidence_class: audit-asserted`).
  The OTHER two audit-asserted edges — `link-outs -> metadata-assembler`
  (background-knowledge) and `link-outs -> render-core` (shared-state, browser
  state) — are NOT structural import edges and STAY (D-4 item 2; review.md D-4
  table). Verify each audit-asserted edge's evidence against a restored
  structural edge before dropping; drop only those now covered by a real static
  edge (R9).

### Step C5 — Re-validate under hardened `--strict` (D-4 item 2)

- Run `validate_components.py --strict` on the new snapshot. It MUST PASS: D-1
  restored the facade edge that un-orphans `run-metadata`, so the hardened gate
  (Phase B) does not spuriously fail the corrected repo.

### Step C6 — Write the metrics-shift changelog (D-4 item, "Surface the change")

- Write a short changelog/diff note in the regenerated snapshot dir (or the
  repo's audit README) listing the expected deltas from `decisions.md` D-4
  "Intended metrics shift" (L570-590), so reviewers reading the diff understand
  it is the fidelity correction, not a regression:
  - `fan_out` of each package init increases to its true re-export count
    (`api/__init__.py` 1→5, `rendering/__init__.py` 1→6,
    `rendering/js/__init__.py` 1→9, root init →7).
  - `fan_in` of the package ROOT `__init__.py` decreases (spurious inbound edges
    from sibling inits disappear).
  - `fan_in` of re-export TARGETS (e.g. `run_metadata.py`) increases (true
    production importers regained).
  - `instability` recomputes for all the above; `run_metadata.py` moves off its
    bogus test-only `fan_in: 1`.
  - `refactor_pressure` increases for facade targets whose `fan_in` was
    under-counted.
  - edge count rises (collapsed 1-edges expand to N); `evidence` strings on the
    previously-collapsed init edges now match their source lines.
  - No metric code changes — the registry recomputes from corrected fan counts.

### Step C7 — Commit the new snapshot (D-4 item 5)

- Commit the new dated snapshot dir; leave prior snapshots untouched.
  (Branch / submodule awareness per `README.md` entry conditions.)

---

## Gates / acceptance-criteria checklist

- [ ] All consuming repos with `docs/_meta/architecture-audit/` discovered
      (C1); each regenerated.
- [ ] New dated snapshot produced by the CORRECTED extractor (C2); prior
      snapshots untouched / immutable.
- [ ] Rollup + render refreshed (C3).
- [ ] **N-2 guard:** `channel-owners` confirmed to carry the restored DERIVED
      STATIC edge BEFORE the redundant `audit-asserted` compensating edge is
      dropped (C4). Only the one D-1-covered edge dropped; the two non-structural
      audit-asserted edges retained.
- [ ] New snapshot PASSES hardened `--strict` (C5); the witness repo's
      `run-metadata` is now PRODUCTION-connected.
- [ ] Metrics-shift changelog written with the D-4 deltas (C6).
- [ ] New snapshot committed; prior snapshots remain byte-untouched (C7).

---

## Dependencies

- **Upstream:** Phases A AND B must be landed and VERIFIED. A supplies the
  corrected extractor (C2) and the restored facade edge (C4); B supplies the
  hardened `--strict` that C5 validates against (`decisions.md` "Proposed phase
  decomposition": "Depends on Phases A and B landing and being verified").
- **Downstream:** none — this is the terminal phase.

---

## Freeze note

This file is the executable spec for Phase C. It is applied ONLY after the
`sciagent-extension` code-change freeze is cleared AND the branching decision is
made (see `README.md` "IMPLEMENT-TIME ENTRY CONDITIONS"), and only after Phases A
and B are verified. Phase C touches no toolkit source — it regenerates snapshots
in consuming repos — but it still runs only post-freeze-clearance because it
depends on the corrected toolkit code.
