# R2 — Opus review checkpoint: figure stack (P04–P07)

**Date:** 2026-06-23 · **Reviewer:** Opus (orchestrator) · **Verdict:** PASS

## Scope reviewed
- P04 `feat(config): add figures contract with dual-context floors to analysis template` (175d148)
- P05 `feat(figure-style): add cross-language figure_helpers (R+py) contract` (82688bf)
- P06 `feat(activate): link figure-style helper lib into analysis repos` (669a130)
- P07 `feat(figure-style): add concept skill and register in base role` (7489a09)
- (parallel) P09 `feat(tags): add figure, provenance, planning tags` (2d42903)

## Acceptance gate — "code runs + produces artifacts"
On a REAL `sciagent new project demo --type analysis` + `activate base`:
- CRAFT block figure rule present; `figures: base_size 16` floor in config; `02_analysis/helpers/figure-style` symlink resolves to the real helper; shim `figure_style.{R,py}` rendered; `figure-style` skill symlinked into `.claude/skills`.
- Ran the dependency-free figure pipeline (what `save_overview` does minus the plot) via the shim import path. Produced:
  - `03_results/02_eda/tables/_overview/gsea_overview.csv` (stage sub-layout placement)
  - `03_results/02_eda/README.md` (caption: exactly ONE path-qualified heading = idempotent; has **How to read** + `Script|Function|Config|Input` table; second call REPLACED the finding in place → "UPDATED")
  - `03_results/master/gsea_master.csv` (2 data rows after two runs = dedup on `database`)
  - NES rounded to 9 sig figs (`2.3456789`) = byte-stable re-runs.
- `bash tests/run-all.sh` → 77 passed, 0 failed. R↔py parity: all 12 contract names in both files.

## Plan adherence / cleanliness
- `visualization:` dead-end (base_size=12) removed; `figures:` is the single contract block; CRAFT display floors agree with config defaults (16 / 9).
- Helper lib is symlinked (not copied) — RELATIVE symlink for portability (the absolute skill symlinks in DC are already broken by a `/workspaces/`→`/scratch/current/` move; the relative helper symlink avoids that class of bug). Recorded in manifest → torn down on deactivate.
- Python module loads with stdlib+pyyaml only (lazy matplotlib/pandas) → testable + importable on a bare box. `save_overview` is atomic (figure+table+caption inseparable). Tidied a duplicated comment banner in symlinks.sh.
- figure-style skill: scope=concept (119 body lines), tags `[figure, viz]`, explicit named anti-patterns (ggsave-inside-compute, inline theme, truncated labels, ambiguous glyphs, uncapped top-N, figure-without-table).

## Notes / follow-ups (non-blocking)
- ACTUAL plot rendering (`save_figure` → `<stem>.print.pdf` + `<stem>.screen.png` with per-variant font floors) could NOT be executed here: this box has no Rscript / matplotlib / pandas. The variant logic was code-reviewed and is exercised structurally; full render must be validated on a box with the sci stack (defer to rollout STING smoke IF that env has R/matplotlib; otherwise document as env-gated).
- P13 `validate --check figure-style` will add the static guardrail (inline theme/ggsave/hex grep + base_size floor) — the enforcement layer for what P07 documents.
