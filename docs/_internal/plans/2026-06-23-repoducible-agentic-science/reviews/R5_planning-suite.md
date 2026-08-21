# R5 — Opus review checkpoint: planning suite (P16–P20)

**Date:** 2026-06-23 · **Reviewer:** Opus (orchestrator) · **Verdict:** PASS

## Scope reviewed
- P16 `feat(templates): add gold-standard plan INDEX + phase brief templates` (a0b37fa)
- P17 `feat(reasoning-trace): add provenance skill; persist bio-interpreter/insight-explorer findings` (78c7cfc)
- P18 `feat(roles): add science-architect orchestration overlay` (4f070cf)
- P19 `feat(commands): add /pipeline-plan orchestration command` (d059e12)
- P20 `feat(commands): add /explore-and-plan research-to-plan command` (012e48c)

## Acceptance gate — "commands produce gold-standard plans; review step runs scripts"
- `/pipeline-plan` encodes the REAL gate verbatim: "actually execute the phase's `02_analysis/scripts/NN_*` script (or confirm its committed/logged run), then assert every artifact the phase's §4 Outputs declares under `03_results/` exists AND is non-empty … this runnable+artifacts assertion is the thing the owner otherwise re-types." Model tiering explicit (planner=Opus, implementer=Sonnet, reviewer=Opus). References the P16 templates; composes with (does not replace) `/plan`/`/implement`/`/verify`.
- `/explore-and-plan` encodes persist-every-wave ("a chat-only finding is a failure"), verified-vs-inferred synthesis tagging, the 14839 "planner never re-explores; phase §3 cites the research note" rule, and hands `_SYNTHESIS.md` to `/pipeline-plan --scope-doc`. Supports `--idempotency-peer` for disjoint namespaces.
- Plan templates render with `_subst` (only the intended `{{DATE}}` token remains); phase table `seq|slug|title|tier|concern|depends_on` + Global notes (idempotency, compute→viz, claim ladder, reasoning traces) present; NN_slug brief has the full fixed section schema incl. acceptance-checks + gotchas.
- `reasoning-trace` skill (scope concept, tag provenance) in base; bio-interpreter + insight-explorer now MUST persist to `research/`/`reasoning/` before returning.
- `science-architect` overlay activates atop base (rc 0); `/pipeline-plan` + `/explore-and-plan` resolve from `commands/science/`.
- `bash tests/run-all.sh` → 83 passed, 0 failed.

## Plan adherence / cleanliness
- Role↔command chicken-egg resolved cleanly: P18 shipped the role activatable with no orchestration commands; P19/P20 each self-wired into `commands:`. (Reordered vs the INDEX's P18-before-P19 listing for correctness — a role can't list a command file that doesn't exist yet.)
- Templates live only under `templates/plan/` (one home); the command docs reference them rather than restating the format.

## Notes / follow-ups (non-blocking)
- `science-architect` lists `commit` (shadowing base's `commit`) — harmless soft shadow; could be dropped since base already provides it. Left as-is.
- Wave 6 adds `/add-figure-variant` (P21) + `/interpret-storm` (P22), which self-wire into this same role, and P23 wires CRAFT/helpers into `sciagent new` + retires the dead viz SSOT.
