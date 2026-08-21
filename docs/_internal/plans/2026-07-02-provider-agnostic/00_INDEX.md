# Provider-agnostic SciAgent — plan index

**Date:** 2026-07-02
**Status:** Draft plan, awaiting decisions (see `50_ADRs.md`)
**Owner:** Anton
**Supersedes/extends:** the "shape B" recommendation in `docs/proposals/ai-research/11-architectural-backbone-and-trace-necessity.md` (recommended, never built); the harness-installers deleted in v3.0.0.

## Why this plan exists

Two triggers, one theme.

1. **Concrete breakage.** (a) Claude Code on the *host* workstation refuses to background a
   conversation / open the agents panel — *"session persistence is disabled."* (b) A freshly
   installed Claude Code *inside a container* comes up with only the default dark theme — none of
   the power-user defaults (vim, effort, thinking, statusline, hooks).

2. **The real goal.** SciAgent-toolkit's thesis is *"the scaffold IS the interface"* — the repo an
   agent wakes up in is the prompt (`docs/_internal/plans/2026-06-04_mindpalace-craft-as-architecture.md:18`).
   Today that interface is **materialized for exactly one harness: Claude Code**. The user now runs
   (or wants to run) `codex`, `pi`, `agy`, and `opencode` too. The injected reproducible-science
   context + coding-style guardrails must reach *whichever* coding agent is rolled out in an analysis
   repo — idempotently, provider-agnostically, without recoupling to any one vendor.

The good news, established by the research in `research/`: **the toolkit's *source* layer is already
provider-neutral** (`roles/*.yaml`, `skills/*/SKILL.md`, `system-prompts/*.md`, `AGENTS.md`,
`craft.yaml`). Only the *materialization* layer is hardcoded to `.claude/`. And the industry has
converged on a portable substrate — **`AGENTS.md` + `.agents/skills/`** — that 4 of the 5 target CLIs
read natively. We are not inventing a standard; we are wiring the toolkit onto one that already exists.

## The tiers (implement gradually, each independently shippable)

| Tier | Name | Effort | Ships | Doc |
|------|------|--------|-------|-----|
| **0** | **Minimal fixes** | ~1 hr | Session-persistence unblock; user-level settings seeding | `10_FIXES_tier0.md` |
| **1** | **Minimal implementation** | ~2–3 days | `AGENTS.md`-first context substrate; `sciagent provision` verb (user-level, cross-provider); devcontainer wiring | `20_provisioning_tier1.md` |
| **2** | **Medium — adapter layer** | ~1–2 wks | Refactor materialization into a provider dispatch table; codex/agy/opencode adapters; `sciagent activate --harness` | `30_adapters_tier2.md` |
| **3** | **Hard — pi role routing** | ~2–4 wks | pi extension that fans a SciAgent *role* out as parallel/chained isolated sub-agents; auto/explicit role selection; guardrail hooks; ADR-011 trace substrate | `40_pi_extension_tier3.md` |

Each tier is a valid stopping point. Tier 0 fixes the pain today. Tiers 1–2 realize "shape B"
(harness-agnostic core + thin per-harness adapters). Tier 3 is the ambitious pi-native fan-out.

## Reading order

1. `01_ARCHITECTURE.md` — the philosophy (what context is, how it propagates, why) and the
   provider-agnostic model (AGENTS.md substrate, source-vs-materialization split, dispatch table).
2. `10_FIXES_tier0.md` → `40_pi_extension_tier3.md` — the tiers, in order.
3. `50_ADRs.md` — **open decisions that gate the tiers.** Read before implementing Tier 1+.
4. `60_blast_radius.md` — every file that changes, per tier; refactor map; test surface.
5. `research/` — source-grounded evidence (CLI capability matrix, bug root causes, prior-art map).

## The one-paragraph summary

Keep the neutral source layer. Make `AGENTS.md` (not `.claude/`) the canonical injection target and
ship the one-line Claude shim the umbrella already uses. Add a **provider dispatch table** so
`activate`/`provision` fan the same role + CRAFT + skills into each harness's own convention
(`.claude/`, `.agents/`, `.pi/`, `~/.codex/`, `.opencode/`). Seed **user-level** defaults in
containers non-destructively. And for pi specifically — the one harness where we can go deeper — ship
a real extension that routes SciAgent roles into parallel/chained sub-agents. Fix the two bugs first;
they are a one-line env change and a ~30-line seeding function.
