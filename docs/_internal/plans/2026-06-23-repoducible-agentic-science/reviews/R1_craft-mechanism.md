# R1 — Opus review checkpoint: CRAFT block mechanism (P01–P03)

**Date:** 2026-06-23 · **Reviewer:** Opus (orchestrator) · **Verdict:** PASS

## Scope reviewed
- P01 `feat(block): parameterize managed blocks by id (ROLES, CRAFT)` (a31c963)
- P02 `feat(craft): render CRAFT block from craft.yaml in activate/inject/eject` (6140408)
- P03 `feat(craft): author canonical five-convention CRAFT text` (90f8474)
- (baseline) `fix(new): scope research/ namespace to analysis projects` (7013a05)

## Acceptance gate — "code runs + produces artifacts"
End-to-end smoke with the REAL toolkit (`bin/sciagent activate base` in a throwaway project):
- ROLES + CRAFT blocks both render into AGENTS.md; hand-written user content preserved.
- `block_hash_check` returns 0 for ROLES and 0 for CRAFT (independent drift scoping).
- Re-activate is byte-identical (both blocks idempotent).
- `sciagent validate --quiet` rc=0; `sciagent status` reports the managed block hash OK.
- `sciagent deactivate` restores AGENTS.md byte-for-byte.
- `bash tests/run-all.sh` → 75 passed, 0 failed.

## Plan adherence / cleanliness
- block.sh parameterized by block-id with `ROLES` default — existing ROLES tests pass byte-identically; ROLES BEGIN marker is byte-identical (`<!-- BEGIN SCIAGENT:ROLES v1 hash=<40hex> -->`).
- SSOT invariant held: marker framing + `sha1sum` live only in block.sh; `test_block_marker_boundary` updated to guard the new helper names + the `<!-- BEGIN SCIAGENT:` format literal.
- craft.sh is a clean leaf depending only on block.sh; craft.yaml is the single source of truth (floors + body block scalar, `body:` last key).
- `craft_render_and_write` no-ops when craft.yaml is absent → zero blast radius on the fake-toolkit test fixtures; the existing activate/deactivate tests are untouched.
- Canonical CRAFT block renders to 11 lines incl. markers (budget ≤25), all five conventions, all `{{floor}}` tokens substituted.

## Notes / follow-ups (non-blocking)
- `sciagent status` surfaces only the ROLES block hash. CRAFT hand-edit drift is not yet reported there. Fold CRAFT drift/version visibility into the P13 `freshness` check (and optionally a status line).
- craft.yaml `floors:` are *display* values for the always-on rule text; the authoritative numeric figure floors a script reads land in `analysis_config.yaml:figures` (P04). Keep the two in agreement (P04 seeds the same defaults).
