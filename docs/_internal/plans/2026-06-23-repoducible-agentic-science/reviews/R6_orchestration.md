# R6 — Opus review checkpoint: orchestration + full integration (P21–P23)

**Date:** 2026-06-23 · **Reviewer:** Opus (orchestrator) · **Verdict:** PASS

## Scope reviewed
- P21 `feat(commands): add /add-figure-variant namespaced figure-family command` (1dd1abe)
- P23 `feat(scaffold): reference CRAFT/figure-style/plans in template; supersede dead viz guideline` (5e88610)
- P22 `feat(commands): add /interpret-storm multi-wave interpretation command` (401d66a)

## Acceptance gate — full new-project smoke (scaffold → activate → plan → figure → caption)
On a fresh `sciagent new project --type analysis` + `activate base science-architect`:
- **A figures:** CRAFT block carries the dual-context figure rule ("legible … journal column … back of a room"); `figures.base_size: 16` floor in config; `figure-style` skill + `figure-audit` agent symlinked; `02_analysis/helpers/figure-style` lib symlink resolves; `figure_style.{R,py}` shim rendered.
- **B placement:** stage-based `stages:` + `by_contrast_dir`/`overview_dir` keys in config; template documents `03_results/<stage>/{figures,tables}`.
- **C readme:** CRAFT readme-adjacency rule; `captions` agent symlinked; Stop caption-sweep hook scaffolded.
- **D planning:** AGENTS.md namespace table has the `docs/_internal/plans/` row; all four orchestration commands resolve (`/pipeline-plan`, `/explore-and-plan`, `/add-figure-variant`, `/interpret-storm`); plan templates available.
- **E reproducibility:** CRAFT no-ephemeral rule; `reasoning-trace` skill symlinked; PreToolUse no-ephemeral hook scaffolded; `reasoning/` dir present.
- **Figure pipeline (dependency-free path) actually produces artifacts:** `03_results/02_eda/tables/_overview/demo.csv` (same-stem table neighbor) + `README.md` caption with How-to-read (×1) + `Script|Function|Config|Input` provenance table (×1) + 9-sig byte-stable rounding (`2.12345679`).
- `validate --check all` → rc 0 (soft-warn, non-blocking). `bash tests/run-all.sh` → 83 passed, 0 failed.
- `activate base science-architect` rc 0; both ROLES + CRAFT blocks render; human template prose + managed CRAFT block coexist.

## Plan adherence / cleanliness
- `/add-figure-variant`: compute-first→verify-checkpoint→viz(dual variants via save_overview)→Opus review with a verbatim namespace-isolation grep gate + figure-audit + mandatory captions pass.
- `/interpret-storm`: 4 waves (web → interpret+synthesis(verified/inferred) → design → hand off to /add-figure-variant), every wave persisted before the next ("a chat-only wave is a failure").
- P23 retired the dead viz SSOT: `docs/guidelines/visualization.md` is now a SUPERSEDED redirect to `skills/figure-style` + `lib/figure-style` + `analysis_config.yaml:figures`; the `base_size=12`/`theme_publication` dead-end is gone. AGENTS.md template references the figure-style shim + the toolkit-managed CRAFT block + the plans/ namespace.

## Notes / follow-ups (non-blocking)
- ACTUAL plot rendering (`.print.pdf`/`.screen.png`) remains env-gated (no R/matplotlib here); covered structurally + by code review (see R2). Defer to a sci-stack box.
- `git init` before the project's `.gitignore` is effective can trip `activate`'s pre-flight docs-layout hard-fail (pre-existing, correct behavior — refuses to mount into a repo that would leak `docs/_internal/`); `sciagent new` seeds the gitignore, so the normal flow is unaffected.
- Remaining: P24 backfill (STING/DC/14839 re-activate + gap fixes), tag bump, repoint feeding projects, STING self-test.
