# R3 — Opus review checkpoint: captions / curator / handoff / figure-audit (P08–P12)

**Date:** 2026-06-23 · **Reviewer:** Opus (orchestrator) · **Verdict:** PASS

## Scope reviewed
- P08 `refactor(scrna-conventions): adopt stage-based results layout` (9af0dca)
- P10 `feat(captions): path-qualified headings, how-to-read, post-figure cleanup` (4be1a62)
- P11 `feat(curator,handoff): tie captions to committed scripts; log scripts/artifacts/decisions` (f81d696)
- P12 `feat(figure-audit): add static-figure legibility reviewer to base role` (0ad0ac1)

## Acceptance gate — consistency across describers (the real risk for prompt-agents is contract DRIFT)
The canonical caption shape (path-qualified `## ` heading + `**How to read:**` + `| Script | Function | Config | Input |` table) is consistent across every place that describes it:
- authored in `lib/figure_helpers.py::_render_caption_section` + `agents/captions.md` (authority),
- referenced by `skills/figure-style/SKILL.md`,
- checked by `doc-curator` C3 (path-qualified + how-to-read + table) and the new C4 (caption `Script:` resolves to a committed `02_analysis/scripts/` path — the C→E provenance tie),
- presence-verified by `figure-audit` checklist item (j).
- `handoff` correctly does NOT author captions; it scans for uncaptioned artifacts and now lists `## Scripts run` (committed paths), `## Artifacts produced`, `## Open decisions → reasoning trace`.
All four agents symlink into a freshly scaffolded project on `activate base`. `bash tests/run-all.sh` → 77/0; `validate --quiet` → 0.

## Plan adherence / cleanliness
- scrna-pipeline-conventions flat layout (`{checkpoints,tables,plots,...}`) DELETED; converged on stage-based `03_results/<stage>/{figures,tables}/` with objects/master/interactive/_scratch at root; defers to figure-style + CRAFT rather than restating a rival layout. config.py sketch updated to `stage_dir()`/`figures(stage)`/`tables(stage)`.
- captions: mandatory-after-figure-change framing + deletion handling + save_overview-backstop note; idempotent replace-in-place documented.
- figure-audit: borrows graphic's Tufte rigor but scoped to RENDERED static figures in base; 10-item D-theme checklist; "what you do NOT do" (review, don't re-render); per-figure verdict table.

## Notes / follow-ups (non-blocking)
- These are LLM-driven agents; their *deterministic* backstops land in P13 (`validate --check captions` mirrors doc-curator C3 harness-agnostically; `--check provenance` mirrors C4) and P14 (Stop hook caption sweep). R3 verified the prompts are mutually consistent; P13/P14 verify the machine enforcement.
- figure-audit's image inspection (opening `*.screen.png`) is only meaningful once real figures render (needs the sci stack) — its file-presence/caption checks work regardless.
