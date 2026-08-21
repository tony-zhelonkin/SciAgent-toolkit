# Code blast radius & refactor map

What changes, per tier, with file:line anchors. Read before touching the mutating core. "New" = new
file; "Edit" = surgical change; "Refactor" = structural move (behavior-preserving).

---

## Tier 0 — Minimal fixes (tiny, non-structural)

| Action | Target | Note |
|--------|--------|------|
| Edit | `~/.bashrc:163` (host, not repo) | Comment `CLAUDE_CODE_SKIP_PROMPT_HISTORY=1`. Not version-controlled. |
| Edit | `lib/sciagent/claude_settings.sh` | Add `claude_settings_ensure_user_defaults()` (~25 LOC), mirrors `ensure_project_defaults` (lines 198–233). |
| New | `templates/user/.claude/settings.json.template` | Hooks-free, user-path statusline (ADR-P2). |
| Edit | `.devcontainer/scripts/postcreate.sh` (umbrella) | Call the user-seed (stopgap) or `provision` (tier 1). |
| Edit | 4× `*/.devcontainer/devcontainer.json` | Add a `postCreateCommand` (they have none today). |
| New | `lib/sciagent/*` doctor check | Warn if `CLAUDE_CODE_SKIP_PROMPT_HISTORY` set. Folds into `validate.sh`. |

**Blast radius:** minimal. No existing test should change behavior. Add
`tests/test_user_settings_seed.sh`.

---

## Tier 1 — Substrate + `provision` (additive, low structural risk)

| Action | Target | Note |
|--------|--------|------|
| New | `lib/sciagent/provision.sh` | The `provision` verb. Reuses `claude_settings.sh` merge helpers + a new `harness_detect()`. |
| Edit | `bin/sciagent` (verbs at :21–53, `VERB_MODULES` :87–98) | Register `provision`. |
| Refactor | `lib/sciagent/status.sh:300–302,631,647` | Extract the existing `.pi/` harness probe into reusable `harness_detect()` covering all 5. |
| New | `templates/global/AGENTS.md.template` | The single-sourced global context file (T1.3). Hand-authored. |
| Edit | `lib/sciagent/activate.sh` (near :264–276) | Guarantee `CLAUDE.md`=`@AGENTS.md` shim always/idempotently (ADR-P4). |
| Edit | devcontainer postcreate (5×) | Run `provision --harness all` at create. |

**Blast radius:** `provision` is a new, isolated verb; the only edit to hot paths is the shim guarantee
in `activate.sh` (idempotent, additive). Tests: `test_provision_*` (detect, seed, idempotency,
non-clobber).

---

## Tier 2 — Adapter layer (the load-bearing refactor)

This is where the real risk lives — it restructures the mutating core.

| Action | Target | Note |
|--------|--------|------|
| Refactor | `lib/sciagent/activate.sh:97–98, 226–251` + `claude_settings.sh:145–233` | **Extract Claude backend → `lib/sciagent/harness/claude.sh`** behind the 8-function contract. Behavior-preserving. |
| New | `lib/sciagent/harness/common.sh` | Shared neutral work: CRAFT/ROLES block render (from `craft.sh`/`block.sh`/`stack.sh`), `.agents/skills` symlink, manifest. |
| New | `lib/sciagent/harness/{codex,agy,opencode,pi}.sh` | ~100–200 LOC each; capability-honest. |
| Edit | `lib/sciagent/activate.sh` | Add `--harness` flag; loop detected/requested adapters. |
| Edit | `lib/sciagent/symlinks.sh:34–200` (manifest schema) | Per-harness manifest sections for clean teardown. |
| Edit | `lib/sciagent/{deactivate,eject}.sh` | Iterate adapters on teardown. |
| Edit | `lib/sciagent/validate.sh` | Per-harness materialization checks + capability reporting. |

**Regression guard (mandatory):** the Claude extraction must keep these green unchanged —
`test_claude_settings_lifecycle.sh`, `test_activate_{solo,stack,idempotent,max_stack,renders_craft}.sh`,
`test_inject_*`, `test_manifest_ownership.sh`, `test_eject_*`. Do the extraction as step 1, prove green,
*then* add adapters. New: `test_activate_harness_{opencode,codex,pi,agy}.sh`, `test_harness_capabilities.sh`.

**Do not touch:** the source layer (`roles/*.yaml`, `skills/`, `system-prompts/`, `craft.yaml`,
`AGENTS.md` templates). It's already neutral; that's the whole point.

---

## Tier 3 — pi extension (new subsystem, self-contained)

| Action | Target | Note |
|--------|--------|------|
| New | `sciagent-pi/` package (vendored) | `package.json` `"pi".extensions`, `src/{index,roles,fanout,guardrails}.ts`, `schema.json`. Modeled on `@aliou/pi-guardrails`; fan-out vendored from pi's `examples/extensions/subagent/`. |
| Edit | `lib/sciagent/harness/pi.sh` | Register the package in `.pi/settings.json:extensions[]`; generate role→`.pi/agents/*.md`. |
| New | role→pi-agent generator | Reads `roles/*.yaml` + `system-prompts/` + CRAFT; emits pi agent markdown. |
| Conditional | ADR-011 trace substrate | Normalizers for Claude + pi native JSONL (~150 LOC each) if promoted. |

**Blast radius:** almost entirely additive and isolated (a TS package + one adapter + a generator).
The only core coupling is the pi adapter, which Tier 2 already introduced. Pin pi 0.79.4; add a version
check. Tests: extension load, single/parallel/chain fan-out, auto vs explicit selection, guardrail block.

---

## Summary — refactor pressure by tier

| Tier | Structural risk | Hot-path edits | New isolated code | Test debt |
|------|----------------|----------------|-------------------|-----------|
| 0 | none | 1 fn + templates + devcontainer | doctor check | 1 test |
| 1 | low | 1 idempotent shim guarantee | `provision.sh`, `harness_detect` | ~4 tests |
| 2 | **high (core extraction)** | `activate`/`claude_settings`/manifest/teardown | 5 adapter files | ~6 tests + regression suite must stay green |
| 3 | medium (new subsystem) | pi adapter only | `sciagent-pi/` package | ~5 tests |

**The one dangerous step is the Tier-2 Claude extraction.** Everything else is additive. Sequence:
extract-and-prove-green *before* any new adapter. If you only ever do Tiers 0–1, you fix the bugs and
get the portable `AGENTS.md` context win for near-zero structural risk.
