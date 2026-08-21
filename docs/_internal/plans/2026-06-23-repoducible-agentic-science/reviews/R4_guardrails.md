# R4 — Opus review checkpoint: guardrails (P13–P15)

**Date:** 2026-06-23 · **Reviewer:** Opus (orchestrator) · **Verdict:** PASS (with one review-fix applied)

## Scope reviewed
- P13 `feat(validate): add opt-in figure-style/results-layout/captions/provenance/freshness checks` (2a39073)
- P15 `feat(update): add sciagent update verb to re-pin submodule and re-activate` (217c008)
- P14 `feat(hooks): ship PreToolUse no-ephemeral + Stop caption-sweep hooks` (af20086)

## Review-fix applied during R4
P13's freshness check originally grepped a `BEGIN SCIAGENT:CRAFT v[0-9]+` marker literal (a marker-format leak that the boundary test only caught literally because of the missing `<!-- ` prefix) AND compared a version integer that is always `v1` (inert, since block.sh hardcodes the marker version). Replaced with a **hash-staleness** comparison: the repo's stored CRAFT hash (`block_stored_hash AGENTS.md CRAFT`) vs the hash the current `craft.yaml` renders (`_craft_render_body` + block.sh canonicalisation). This (a) removes the marker literal from validate.sh (SSOT honored — marker knowledge stays in block.sh), (b) is actually functional (detects ANY craft.yaml change, not just integer bumps). Added `block.sh craft.sh` to the `validate` verb's module closure; updated the freshness test to the hash mechanism.

## Acceptance gate — "fail closed on violations, soft-warn by default, hooks fire"
On a real scaffolded project:
- `activate base` rc=0 with the full guardrail layer scaffolded: CRAFT block, figure-style symlink, `.claude/hooks/`, `.claude/settings.json` all present. The opt-in invariant holds — `activate`'s internal `validate --quiet` is unchanged; figure/caption/layout findings never block mounting.
- `validate --check all` on a fresh project → rc 0.
- results-layout violation (artifact at `03_results/` root) → default rc 0 with WARN; `--strict` → rc 1 (fail closed).
- no_ephemeral PreToolUse hook: ephemeral `03_results/` write → warn (exit 0) / `SCIAGENT_STRICT=1` block (exit 2); `_scratch/` sanctioned (silent); committed `02_analysis/scripts/*` never blocked even under strict.
- `update --no-pin` → re-activates, CRAFT block still valid.
- `bash tests/run-all.sh` → 82 passed, 0 failed.

## Plan adherence / cleanliness
- Five checks are opt-in via `--check` (repeatable / `all`), soft-warn default, hard-fail under `--strict`; `_scratch/`+`$TMPDIR` always exempt; absent dirs are clean no-ops (software/empty projects pass).
- Hooks ship as `*.template` under `templates/project/_common/.claude/` so `new project` scaffolds them; invoked via `bash "$CLAUDE_PROJECT_DIR/..."` (no exec-bit dependency); in `.claude/settings.json` (separate from toolkit-owned `settings.local.json:outputStyle`); project-scoped; degrade gracefully (validate is the CI mirror).
- `update` re-pins the submodule (best-effort, `--no-pin` to skip) then re-execs `activate` so the fresh lib loads; guides clearly when no active stack.

## Notes / follow-ups (non-blocking)
- The figure-style validate check's grep heuristics (inline theme/ggsave/hex) are conservative and operate on `*_viz.*` scripts; real-world tuning may be needed once run against STING/DC at backfill (P24).
- The CRAFT marker version remains a static `v1`; the freshness signal now rides on the content hash (better) + submodule-commit ancestry, so the static version is no longer load-bearing.
