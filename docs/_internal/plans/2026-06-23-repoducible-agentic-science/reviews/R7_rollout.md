# Rollout — tag v3.2.0, propagation, STING smoke test

**Date:** 2026-06-23 · **By:** Opus (orchestrator)

## Done
- **Tagged `v3.2.0`** (annotated, at dev `1656549`) and **pushed `dev` + `v3.2.0` to origin** (github.com/tony-zhelonkin/SciAgent-toolkit). CHANGELOG updated.
- **Confirmed global runtime already on v3.2.0:** `~/.local/bin/sciagent` → `…/scbio-docker/toolkits/SciAgent-toolkit/bin/sciagent` (the edited toolkit, dev HEAD = v3.2.0). So every project's `sciagent` command already runs v3.2.0 — the runtime switch is automatic.
- **STING re-activated to v3.2.0 + smoke-tested** (`sciagent activate base pathway-signature`):
  - AGENTS.md gained the `SCIAGENT:CRAFT` block (figure rule present); ROLES block re-rendered.
  - `02_analysis/helpers/figure-style` symlink created (relative, resolves to the toolkit lib); figure_helpers imports 6/6 contract functions; `direction_cue(2.1)` → "↑ up".
  - New base capabilities active: `figure-style` skill, `reasoning-trace` skill, `figure-audit` agent.
  - `validate --check all` (rc 0, soft-warn) CAUGHT REAL violations in STING's actual viz scripts: `02_de_limma_trend.R` (saves figure without `project_theme()`); `03_decoupler_tf_viz.R` (inline `theme()` line 77, raw hex line 207, no `project_theme()`); `03e_heat_main_regulators_viz.R` (inline `theme()` line 61). These are exactly the figure-legibility anti-patterns the owner used to hand-steer.

## Deliberately NOT done (protecting the user's uncommitted work)
- **Submodule pin bumps + commits in live project repos were NOT forced.** State at rollout:
  - STING submodule: DIRTY (`M skills/bulk-rnaseq-activity-inference/SKILL.md`, untracked reference).
  - 14839 submodule: DIRTY (`M skills/shinymultiome-uio-host/SKILL.md`).
  - DC submodule: clean; but DC tree has 13 uncommitted files, STING/14839 trees 27/various.
  Checking out v3.2.0 over a dirty submodule, or committing a gitlink bump into a repo full of the user's in-progress work, risks their changes. The global symlink already provides v3.2.0 at runtime, so nothing is lost by deferring the pin bump.
- STING's re-activation changes (`AGENTS.md`, new helper symlink) are left in STING's working tree (uncommitted) for the user to review/commit with their own work — not auto-committed (live repo + dirty submodule).

## Per-project completion recipe (for the user, when ready)
For each project P in {STING-cGAS-GSE329522, DC-nexus/DC_mouse_cancer, 14839-DM-cGAS}:
1. Commit/stash any local edits in `P/01_modules/SciAgent-toolkit` first (STING, 14839 have some).
2. `git -C P/01_modules/SciAgent-toolkit fetch origin --tags && git -C P/01_modules/SciAgent-toolkit checkout v3.2.0`
3. `cd P && sciagent update` (re-pins + re-activates: re-renders ROLES+CRAFT, re-links figure-style) — or `sciagent activate <stack>`.
4. `git -C P add 01_modules/SciAgent-toolkit AGENTS.md 02_analysis/helpers/figure-style && git -C P commit -m "chore: adopt SciAgent-toolkit v3.2.0 (CRAFT + figure-style)"`
5. Backfill to fully use the helpers: add the `figures:` block to `P/02_analysis/config/analysis_config.yaml` (copy from the toolkit template), then refactor the viz scripts `validate --check figure-style --strict` flags to call `project_theme()`/`save_overview()` instead of inline `theme()`/`ggsave()`/hex.
