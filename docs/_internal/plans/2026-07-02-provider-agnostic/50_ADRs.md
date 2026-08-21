# Open ADRs — decisions that gate this plan

Uses the toolkit's decision taxonomy (`kickoff.md §2`): **[Arch]** hard to reverse, think slowly;
**[Eng]** reversible at refactor cost, pick the cheapest option that meets the bar; **[Taste]** pick
fast, document why.

Two groups: **P-series** = new decisions this plan raises. **Inherited** = pre-existing deferred ADRs
that this plan reactivates or must respect (from `docs/proposals/deffered.md`, `kickoff.md §9`).

---

## New decisions (P-series)

### ADR-P1 — Replace `CLAUDE_CODE_SKIP_PROMPT_HISTORY` privacy flag [Taste]
**Context.** The env var the user set for privacy also disables session persistence + agents (Bug 1).
**Decision needed.** Drop it entirely, or replace with bounded retention (`cleanupPeriodDays: N` in
`~/.claude/settings.json`)? **Recommendation:** drop the env var; if reduced retention is still wanted,
use `cleanupPeriodDays`. Add a `validate`/doctor warning so it can't silently recur.
**Blocks:** Tier 0 Fix 1. **Reversible:** trivially.

### ADR-P2 — User-level settings template: hooks-free variant [Eng]
**Context.** The gold `settings.json.template` references `$CLAUDE_PROJECT_DIR`-scoped hooks + a
project statusline path. Those are wrong at *user* level.
**Decision needed.** Maintain a separate hooks-free user template with a `~/.claude/statusline.sh`
path, or parameterize one template? **Recommendation:** a small dedicated `user/settings.json.template`
(no hooks; user-path statusline). Keep hooks a project-`activate` concern.
**Blocks:** Tier 0 Fix 2 / Tier 1. **Reversible:** yes.

### ADR-P3 — `CLAUDE_CONFIG_DIR` → persistent volume in containers [Eng]
**Context.** Container home is ephemeral; login+settings+sessions are lost on rebuild.
**Decision needed.** Point `CLAUDE_CONFIG_DIR` at a mounted named volume (survives rebuilds, real dir
not symlink), or accept re-seeding each rebuild via `provision`? **Recommendation:** optional; adopt
only if rebuild-churn is painful. Orthogonal to the toolkit (compose/devcontainer change).
**Blocks:** nothing (optional). **Reversible:** yes.

### ADR-P4 — `AGENTS.md` as canonical target; which shims to auto-materialize [Arch]
**Context.** Flipping the default injection target from `.claude/` to `AGENTS.md` is the plan's
load-bearing move. `CLAUDE.md`=`@AGENTS.md` is needed; `GEMINI.md` is usually redundant (agy reads
AGENTS.md natively).
**Decision needed.** Always materialize `CLAUDE.md` shim (yes). Materialize `GEMINI.md` only on request
vs always? **Recommendation:** `CLAUDE.md` always; `GEMINI.md` only if legacy Gemini CLI is a target.
**Blocks:** Tier 1 T1.1 and everything after. **Reversible:** at refactor cost — decide deliberately.

### ADR-P5 — Guardrail portability: "portable floor + per-harness ceiling" [Arch]
**Context.** Layer (c) "unskippable" enforcement differs per harness (Claude settings hooks vs opencode
plugins vs pi `tool_call` vs agy PreToolUse vs codex execpolicy). codex has the weakest hook surface.
**Decision needed.** Accept that the *guaranteed* guardrail is the harness-agnostic `validate --check`
floor (+ CI `check-provenance`), with live hooks as a per-harness best-effort ceiling? **Recommendation:**
yes — this is the only honest cross-provider guarantee. Document per-harness ceiling coverage.
**Blocks:** Tier 2 T2.5, Tier 3 T3.5. **Reversible:** the ceiling is; the floor is the invariant.

### ADR-P6 — `provision` vs `activate` scope boundary [Eng]
**Context.** New `provision` verb (user/global) vs existing `activate` (project/role).
**Decision needed.** Confirm the split: `provision` never activates a role; `activate` never seeds
user-level defaults. **Recommendation:** keep them orthogonal; `provision` at container-create,
`activate` per working session. **Blocks:** Tier 1 T1.2. **Reversible:** yes.

### ADR-P7 — Adapter contract surface + capability honesty [Arch]
**Context.** The `harness_*_ensure_*` contract and the `capabilities` honesty rule (Tier 2).
**Decision needed.** Lock the contract's function set before writing four adapters, so they don't drift.
**Recommendation:** freeze the eight-function contract in `30_adapters_tier2.md:T2.1` after the Claude
extraction proves it. **Blocks:** Tier 2. **Reversible:** at refactor cost.

### ADR-P8 — pi package: vendor vs depend; discovery tree [Eng]
**Context.** Ship SciAgent pi support as a vendored package pinned to pi 0.79.4; point sub-agent
discovery at SciAgent's `.agents/agents/` vs pi's `~/.pi/agent/agents/`.
**Decision needed.** Vendor (yes, per prior-art). Standardize on which agents dir? **Recommendation:**
vendor + pin; generate role agents into pi's native `.pi/agents/`/`~/.pi/agent/agents/` to match the
example's discovery, and treat `.agents/agents/` as the neutral *source* that the extension reads and
renders from. **Blocks:** Tier 3. **Reversible:** yes.

---

## Inherited deferred ADRs (respect / reactivate)

### ADR-011 (proposed) — sciagent trace/result substrate [Arch] — **candidate to promote**
The "missing ADR" (`ai-research/10:220`): ADRs 004/005/006 all depend on a trace substrate that was
never itself an ADR. **Tier 3 pi fan-out is a second trace consumer** → the documented trigger to
un-defer (`ai-research/11:96`). **Recommendation:** promote and decide ADR-011 *before* Tier 3 lands
(sciagent-native JSONL on disk + optional OpenInference exporter; normalizers ~150 LOC/harness from
Claude + pi native JSONL). If Tier 3 is deferred, ADR-011 stays deferred.

### ADR-005 (deferred) — trace recorder + interview distiller
Stays deferred unless ADR-011 is promoted with Tier 3. The lightweight substitute (markdown reasoning
traces + `validate --check provenance`) remains the shipping reproducibility story for Tiers 0–2.

### ADR-004 / ADR-006 (deferred pair) — ERA optimizer + benchmark harness
Not touched by this plan. Note only: if role fan-out (Tier 3) ever wants *measured* role selection, it
reactivates ADR-006. Out of scope here.

### ADR-009 (deferred) — `uv` pinning
Respect as-is (Anton works in Docker with pre-installed envs). A provider-agnostic installer does **not**
reintroduce environment pinning; `provision` seeds *config*, not Python envs.

### Stack-depth-2 cap [Arch, latent]
Flagged as "rubber-banding" (`ai-research/10:226`). Provider-agnosticism does not change it: roles still
stack base+overlay, depth 2, last-wins — per harness. If pi auto-selection ever wants deeper composition,
revisit then, not now.

### Harness-installers-were-removed (v3.0.0) — do not silently resurrect [Arch]
v3.0.0 deleted `install_claude/codex/gemini.sh` deliberately ("roles, agents, skills only"). This plan
does **not** bring back binary installers — the curl-installers stay the user's/devcontainer's job.
`provision` seeds *config + context*, it does not install CLIs. Keep that line bright (record as the
rationale in ADR-P6).

---

## Decision sequencing

Resolve **P1, P2** before Tier 0 (trivial). Resolve **P4** (Arch) before Tier 1 — it's the hinge.
Resolve **P5, P7** before Tier 2. Resolve **P8 + ADR-011** before Tier 3.
