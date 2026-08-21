# Implementation progress — resume manifest

**Batch:** Tier 0 (fixes) + Tier 1 (substrate + `provision`).
**Branch:** `feat/provider-agnostic-tier0-1` (off `dev`) in SciAgent-toolkit.
**Started:** 2026-07-03. **Owner:** Anton.

## Scope decisions locked (user away → conservative defaults)
- **Implementer:** Claude subagents (codex NOT installed on this box; install needs interactive login).
  Swap to codex delegation later once installed/authed. Reviewers = opus.
- **Cross-repo:** toolkit-only this batch. Devcontainer wiring (C1) is a tracked follow-up, NOT done here.
- **Host `~/.bashrc` Fix 1:** validate/doctor warning only; the actual dotfile edit is left to the user.
- **Commits:** feature branch, commit per reviewed stage group. No push without approval.
- **Plan docs** live under gitignored `docs/_internal/` — not committed (by repo convention); this
  manifest persists on disk for resume.

## Stage status

Legend: ☐ todo · ◐ in progress · ☑ done · ⊘ deferred/skipped · ⚠ needs attention

### Group A — Tier 0 fixes  ✅ DONE (commit 30726ac)
- ☑ A1  `claude_settings_ensure_user_defaults()` + `_ensure_user_statusline()` in `claude_settings.sh`
- ☑ A2  hooks-free user template `templates/user/.claude/settings.json.template`
- ☑ A3  `_validate_env_hygiene` (`session-persistence` warn) in `validate.sh`
- ☑ R1  opus review → APPROVE-WITH-NITS (no blocking); non-clobber/exit-code/JSON/idempotency verified
- ☑ commit A → `30726ac`
- NOTE: user-level functions are DEAD until wired via `provision`/postcreate (Group B / follow-up C1)
- NOTE: `tests/test_claude_settings_lifecycle.sh` fails on clean `dev` too (tag-vocab gap in `delegate-cli` vs `tags.yaml`) — pre-existing, unrelated

### Group B — Tier 1 core  ✅ DONE (commit 781e963)
- ☑ B1  `lib/sciagent/harness.sh` — `harness_{bin,config_dir,global_context_path,is_present,detect}` (5 harnesses). status.sh left byte-identical (harness.sh is go-forward API; inline probes noted as intentional dup)
- ☑ B2  `lib/sciagent/provision.sh` `cmd_provision` + registered in `bin/sciagent` (usage/VERB_MODULES/VERB_ENTRY). Does NOT call activate/validate. `--user` is the only Tier-1 scope (fwd-compat)
- ☑ B3  `templates/global/AGENTS.md.template` (terse 6-bullet) + fan-out via block.sh `SCIAGENT:CONTEXT` managed block (idempotent, preserves user bytes)
- ☑ B4  `ensure_claude_md_shim` in `activate.sh` (create/preserve `@AGENTS.md`; idempotent)
- ☑ R2  opus review of Group B → APPROVE-WITH-NITS (no blocking; dry-run/source-order/shim/idempotency verified)
- ☑ commit B → `781e963`
- nits (deferred, cosmetic): `provision.sh` `do_user` unused (fwd-compat); `--harness` followed by a flag consumes it
- NOTE: full suite 87 pass / 2 fail; both (`test_claude_settings_lifecycle.sh`, `test_tags_vocabulary.sh`) pre-exist on `dev`

### Group C — wiring + tests  ✅ DONE (commit 1c1a23d)
- ⊘ C1  devcontainer `postCreateCommand` wiring (cross-repo) — **DEFERRED to follow-up (scope decision)**
- ☑ C2  5 hermetic tests: `test_user_settings_seed.sh`, `test_provision_{detect,context,dryrun}.sh`, `test_claude_md_shim.sh` (auto-discovered by `tests/run-all.sh`)
- ☑ R3  opus review → APPROVE-WITH-NITS (anti-vacuity + hermeticity verified; not flaky despite agy/claude on PATH)
- ☑ commit C → `1c1a23d`
- test-coverage nits (non-blocking, optional follow-up): no all-5-present detect case; agy `~/.local/bin` fallback + non-empty `CODEX_HOME` branches untested; doctor-warning (A3) has no dedicated test

## BATCH STATUS: Tier 0 + Tier 1 core COMPLETE ✅
Branch `feat/provider-agnostic-tier0-1` (off `dev`), 3 commits, NOT pushed. Suite 92 pass / 2 fail
(2 pre-existing). Awaiting user review + push decision.
- `30726ac` Group A — user-level settings seeding + session-persistence warning
- `781e963` Group B — harness detection + `provision` verb + global context fan-out + CLAUDE.md shim
- `1c1a23d` Group C — 5 hermetic tests

## Follow-ups parked for next session
- **C1 devcontainer wiring** (deferred by scope): Meta-Aging umbrella `postcreate.sh` + add
  `postCreateCommand` to the 4 dataset `devcontainer.json` → call `si provision --harness all`.
- **Host `~/.bashrc:163` edit** (Fix 1): user to comment `export CLAUDE_CODE_SKIP_PROMPT_HISTORY=1`.
- **PRE-EXISTING BUG (easy fix, out of this batch's scope):** `tests/test_claude_settings_lifecycle.sh`
  + `test_tags_vocabulary.sh` fail on `dev` because the `delegate-cli` skill (commit `e8508b3`) uses tags
  `delegation/multi-model/headless-cli/codex/gemini/orchestration` not present in `tags.yaml`. Fix =
  add those tags to `tags.yaml`. Recommend as a quick separate commit.
- **Test nits** (optional): all-5-present detect case; agy `~/.local/bin` + non-empty `CODEX_HOME`
  branches; a dedicated test for the A3 `CLAUDE_CODE_SKIP_PROMPT_HISTORY` validate warning.
- **codex delegation**: install + auth codex to run future stages through `codex exec` (this batch used
  Claude implementers because codex was absent).
- **Tier 2** (adapter layer, `30_*`) and **Tier 3** (pi extension, `40_*`) — not started.
- **Open ADRs P4/P5** (Arch) confirm before Tier 2; **ADR-011 + P8** before Tier 3.

## Change log (append as stages land)
- 2026-07-03: branch created; manifest initialized.
- 2026-07-03: Group A implemented + opus-reviewed + committed `30726ac`.
- 2026-07-03: Group B implemented + opus-reviewed + committed `781e963`.
- 2026-07-03: Group C (tests) implemented + opus-reviewed + committed `1c1a23d`. Tier 0/1 batch complete.
