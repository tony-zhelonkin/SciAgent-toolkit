# Tracer-Bullet Audit — Consolidated Synthesis

**Date:** 2026-06-02 · **Synthesizer:** Opus 4.8 · **Inputs:** 6 parallel-universe tracer reports (01–06)
**Toolkit:** `/data1/users/antonz/pipeline/scbio-docker/toolkits/SciAgent-toolkit`
**Method:** read all 6 reports → dedup/cluster → verify top claims against source (file:line) → classify → root-cause.

---

## Executive summary (the state of the toolkit)

The **core lifecycle mechanics are sound and trustworthy**: activate → status → inject/eject → deactivate produce clean symlinks, a clean AGENTS.md managed block, and a clean teardown with no orphans for normal usage. 68/69 bash tests pass. The real damage clusters in **one architectural seam — the `requires:` transitive closure** — which is resolved at `activate` but is *invisible* to every machine-readable status surface and *unguarded* in inject/eject. That single seam produces a confirmed **BLOCKER** (ejecting an inject-over-inherited skill silently destroys a live dependency symlink and breaks `tf-footprint-differential-analysis` mid-session), three "missing inherited skills" status bugs, and an unguarded inject. Beyond that seam there is **1 CI-red blocker** (`scope: atomic` lint — pre-existing), **1 confirmed scaffolding bug** (`software-tool` template dir mismatch — user's in-progress rename), a `status --effective` **exit-code-1 bug**, and a layer of **stale docs / orphaned references** (README `roster`, `scrna-atlas`, `planning`; a live real-project manifest pinned to a deleted `planning` role). Net: **~9 real distinct issues**, 1 blocker of each kind (CI-red + runtime-destructive), the rest medium/low. Headline: *the machine works; the requires-closure is the load-bearing crack, and there is no `doctor`/migration path for state that drifts out from under the catalog.*

---

## Convergence table

| # | Finding | Agents | Status | Severity | Category |
|---|---------|--------|--------|----------|----------|
| F1 | `eject` of inject-over-inherited skill destroys the `requires:` symlink, breaks parent skill mid-session | 03 | **CONFIRMED** | BLOCKER | failure / dangling logic |
| F2 | Inherited (`requires:`) skills missing from `--effective`, `--json`, `--source` (text says 39, JSON/effective say 36) | 01, 02, 03 (+04 documents the 3) | **CONFIRMED** | BUG | blank spot / failure |
| F3 | `inject` does not pull transitive `requires:` closure; orchestrator skills land without leaf deps, no warning | 03 | **CONFIRMED** | BUG | edge case / blank spot |
| F4 | `inject` silently accepts an already-inherited skill (guard only reads role YAML, not requires/manifest) | 03 | **CONFIRMED** | BUG | edge case / dangling logic |
| F5 | `software-tool` template dir is `software/`, code renders `$type` = `software-tool/` → silently renders 0 type files | 01, 05 | **CONFIRMED** | BUG | failure (silent) / in-progress rename |
| F6 | `scope: atomic` legacy vocab in `mllmcelltype-consensus-annotation` → 1 test FAIL, CI red | 05 | **CONFIRMED** | BLOCKER (CI) | failure / pre-existing |
| F7 | `status --effective` returns exit 1 when no injected commands (empty-array `:-` loop) | 02 | **CONFIRMED** | BUG | edge case / failure |
| F8 | 1-arg `activate <overlay>` on a 2-deep stack silently nukes the base | 02 | REPORTED (behavior-confirmed by agent) | CONFUSING | edge case / footgun |
| F9 | README Quick Start references nonexistent `sciagent roster` verb | 01 | **CONFIRMED** | ORPHAN | mess from prev session |
| F10 | README line 80 references deleted roles `scrna-atlas` + `planning` | 05 | **CONFIRMED** | ORPHAN | mess from prev session |
| F11 | Real project `DC_hum_verse/.sciagent/manifest.json` stack = `["architect","planning"]`; `planning` no longer exists; no `doctor`/migration | 06 | **CONFIRMED** | ORPHAN | hanging thread / orphaned state |
| F12 | No CLI path from science→skill: `list skills` is names-only; `list tags` doesn't exist though `inject --tag` does | 06 | REPORTED | CONFUSING | blank spot |
| F13 | `pathway-signature` is a strict subset of `base` → `activate base pathway-signature` adds 0 skills | 06 (02/01 corroborate the shadow noise) | REPORTED (counts corroborated) | CONFUSING | dangling logic / taxonomy |
| F14 | AGENTS.md "Active role" template section stays "No role activated yet" after activate | 01 | REPORTED | CONFUSING | edge case |
| F15 | `new project` (no dir) sprays files into CWD + leaks `find: Permission denied` to stdout from `/tmp` | 01 | **CONFIRMED** (code path) | BUG | edge case / footgun |
| F16 | `deactivate _injected` does clean-slate re-activate (drops injected) instead of targeted teardown | 03 | **CONFIRMED** | CONFUSING | edge case |
| F17 | 25 "orphan" skills unreachable except via `inject` | 04, 06 | REPORTED (by design) | ORPHAN (benign) | blank spot |
| F18 | README "context.md" vs actual `docs/_internal/scientific-context.md`; `<scbio-docker>` placeholder unresolved | 01 | **CONFIRMED** | NIT | mess from prev session |
| F19 | 2 skills missing `last-reviewed` (architecture-first-dev, skill-creator) | 05 | REPORTED | NIT | hygiene |
| F20 | Eject error appends `/_injected` to base-only skills; `inject --tag` all-mounted exits 0 silently; misc msg nits | 03 | REPORTED | NIT | edge case |
| F21 | `new project --help` errors (`unknown flag '--help'`); `--title` barely wired into context.md | 06 | REPORTED | NIT | edge case |

---

## Detailed findings (clustered + root-caused)

### CLUSTER A — The `requires:`-closure seam (F1–F4, F7) — the headline

This is **one architectural seam** manifesting as five symptoms. `tf-footprint-differential-analysis` is in `base.yaml` and declares `requires: [tobias-footprint-bindetect, hint-atac-differential-footprint, signac-footprint-visualization]`. `activate.sh` resolves that closure, creates symlinks for all three, and records them in `manifest.symlinks`. But the closure is **only** materialized as on-disk symlinks + manifest entries; it is never represented as a first-class array the other surfaces can read.

**F1 [BLOCKER] — eject destroys an inherited symlink. CONFIRMED.**
- `lib/sciagent/inject.sh:430-456` (`_inject_entry_is_stack_mounted`) only calls `role_load` (role YAML). An inherited skill is *not* in any role YAML, so the guard returns "not mounted" and inject proceeds → creates a redundant `injected` manifest row pointing at the *same* `.claude/skills/<name>` path.
- `lib/sciagent/eject.sh:271-276` (`_eject_remove_symlinks`): `if [[ -L "$claude_path" ]]; then rm "$claude_path"; fi` — unconditional delete, **no check whether activate.sh also owns this symlink**.
- `lib/sciagent/eject.sh:320-329` (`_eject_drop_injected_entry`): rebuilds `manifest.symlinks` dropping that path too.
- Net: the dep symlink is gone from disk *and* manifest; `tf-footprint-differential-analysis` is silently broken for the rest of the session. Agent 03 reproduced this cleanly. This is the most destructive bug found.

**F2 [BUG] — inherited skills absent from `--effective`/`--json`/`--source`. CONFIRMED (3 agents).**
- `status.sh:178-205` (`_status_render_text`) computes `INHERITED_SKILLS` *locally* from `MANIFEST_SKILLS` minus declared/injected, and counts it into the "Skills (N effective)" total.
- `status.sh:411-419` (`_status_render_effective`) iterates only `SKILL_ORDER` + `INJECTED_SKILLS`.
- `status.sh:500-515` (`_status_render_json`) builds `"skills"` from the same two arrays only.
- `status.sh:421-464` (`_status_render_source`) searches only `SKILL_ORDER` + `INJECTED_SKILLS`; for an inherited skill it falls through to `echo "not in active stack: $q" >&2; return 1` (line 462-463) — an **actively wrong** answer, since the skill *is* mounted.
- Root cause: `INHERITED_SKILLS` is a function-local in the text renderer, never hoisted to a global the other three renderers can consume. Result: text=39, JSON/effective=36; any CI/Pi consumer of the machine-readable surface sees an incomplete, self-contradictory manifest.

**F3 [BUG] — inject doesn't resolve the closure. CONFIRMED.**
- `inject.sh` injects exactly one entry and never calls `skill_resolve_transitive` (which `activate.sh` does use). Injecting an orchestrator skill (e.g. `muon-multimodal-analysis`, 8 deps) mounts the orchestrator but leaves 5 leaf skills absent, no warning. The AI then sees a skill referencing tools that aren't mounted.

**F4 [BUG] — inject's stack-mounted guard is closure-blind. CONFIRMED.** Same root cause as F1's first bullet; it is the *entry point* that lets the F1 destruction happen. Fixing the guard to also consult `manifest.symlinks`/closure (warn or refuse) neutralizes F1's trigger.

**F7 [BUG] — `status --effective` exits 1. CONFIRMED.**
- `status.sh:418`: the final loop `for n in "${INJECTED_COMMANDS[@]:-}"; do [[ -n "$n" ]] && echo "$n"; done`. When the array is empty, `${arr[@]:-}` expands to one empty string, the `[[ -n "" ]]` is false, `&&` short-circuits → loop body exits 1 → function (and process) exit 1. Verified independently: `bash -c 'declare -a E=(); for n in "${E[@]:-}"; do [[ -n "$n" ]] && echo "$n"; done; echo $?'` → `Exit: 1`. Any `if sciagent status --effective | grep ...` pipeline sees a spurious failure. This is the *most-likely-to-be-scripted* status mode.

> **The single fix that collapses most of A:** hoist an `INHERITED_SKILLS` (or general "closure-mounted, non-declared") global in `_status_load_state`, feed it to all four renderers, and have inject/eject treat closure-mounted paths as shared (refuse-or-warn on inject; preserve-on-eject). F7 is independent (one-line `return 0` or guard the empty loop).

### CLUSTER B — `software-tool` scaffolding (F5) — CONFIRMED, user's in-progress rename

- `new.sh:170`: `_render_tree "$proj_tpl/$type" "$dir" "$force"` with `$type="software-tool"` → looks for `templates/project/software-tool/`.
- Verified on disk: directory is `templates/project/software/` (contains `AGENTS.md.template`, `README.md.template`, `tool_config.yaml.template`); `templates/project/software-tool/` is **ABSENT**.
- `new.sh:82`: `_render_tree` guards `[[ -d "$tpl_root" ]] || return 0` → silent success, 0 type-specific files rendered. A `software-tool` project gets only `_common/` files (CLAUDE.md, .gitignore, docs/_internal/README.md) — **no AGENTS.md, no tool_config.yaml**.
- `test_new_project_types.sh` passes because it checks directories exist, not that templates rendered — a **test blind spot**.
- This is the documented `software-tool → software` rename-in-progress (CLAUDE.md notes it). Fix = rename dir to `software-tool/` OR add a `$type→dir` map at new.sh:170. Note: report 04 incorrectly called this "Resolved — file name matches role name"; that's about the *role* YAML, not the *template* dir. The template dir is genuinely mismatched.

### CLUSTER C — Stale references & orphaned state (F9, F10, F11) — CONFIRMED

- **F9** README:66 `sciagent roster` — no such verb (CLI has `status`/`list`). First-run reader hits an error.
- **F10** README:80 references `scrna-atlas` + `planning`, both deleted in v1.7.0 (CHANGELOG). Verified present in README.
- **F11 (most interesting orphan):** the *live* real project `/scratch/.../DC_hum_verse/.sciagent/manifest.json` has `"stack": ["architect","planning"]`. `planning` is gone from the catalog. A returning scientist running `status`/`activate` there is pinned to a role that no longer exists, and there is **no `sciagent doctor` / migration path** to detect or repair a manifest that drifted out from under a catalog change. This is the clearest "hanging thread / mess from a previous session" in the whole audit, and no single tracer path was designed to find it — it surfaced because agent 06 grounded against a real project.

### CLUSTER D — Taxonomy / discoverability (F12, F13, F17) — REPORTED, design-level

- **F13** `pathway-signature` ⊂ `base` (all 8 of its skills are among base's 36) → `activate base pathway-signature` yields 39 skills, identical to `activate base`, with every overlay line tagged `(shadows base)`. The taxonomy implies a composition that doesn't compose. Either make `pathway-signature` a lean standalone (use *instead of* base) and say so, or stop base from swallowing it.
- **F12** No science→skill discovery from the CLI: `list skills` (status.sh:741-749) prints names only; `tags.yaml` holds the rich descriptions but `sciagent list tags` doesn't exist even though `inject --tag` is first-class. The tag system is half-surfaced.
- **F17** 25 skills in no role (reachable only via `inject`). Agent 04 and 06 agree this is **by design** (à-la-carte injection), not a bug — but it compounds F12: 25 dormant skills with no CLI-visible way to discover what they do.

### CLUSTER E — UX footguns & nits (F8, F14, F15, F16, F18–F21)

- **F8 [CONFUSING]** 1-arg `activate <overlay>` on a 2-deep stack silently replaces the whole stack (drops base). README says "replaces current stack" — technically true, easy to miss. No `--overlay` shorthand. Real trap for the "just add one more" mental model.
- **F15 [BUG]** `new project` with no dir defaults to `.`; run from `/tmp` it scaffolds into `/tmp` *and* `_seed_gitkeep` (`new.sh:74`) runs `find "$root" -type d -empty` over all of `/tmp`, leaking `find: ... Permission denied` to **stdout** (no `2>/dev/null`). README:54 shows exactly `sciagent new project` with no dir. Two fixes: `2>/dev/null` on the find; README should show an explicit `<dir>`.
- **F16 [CONFUSING]** `deactivate _injected` (deactivate.sh:45-48) calls `cmd_activate base` — a clean-slate re-activate that *drops* injected entries with a warning, rather than cleanly ejecting them. Looks like an activation, not a deactivation.
- **F14 [CONFUSING]** AGENTS.md template "Active role" section stays "No role activated yet" after activate (managed block is appended below, never splices that section). An AI reading top-down hits the contradiction first.
- **F18/F20/F21 [NIT]** README "context.md" vs real `scientific-context.md`; unresolved `<scbio-docker>` placeholder in Next-steps; eject error spuriously appends `/_injected`; `inject --tag` all-mounted exits 0 silent (no stdout); `new project --help` errors; `--title` weakly wired. All low-severity polish.

### Pre-existing CI blocker (F6) — CONFIRMED, not closure-related

- `tests/run-all.sh`: **68 pass, 1 fail** (verified by running it). `test_skill_scope_lint.sh` fails because `skills/mllmcelltype-consensus-annotation/SKILL.md` carries `scope: atomic` (legacy vocab banned by ADR-003) with `last-reviewed: 2026-05-29` (≥ 2026-05-24 cutoff, so no legacy exemption). One-line fix (`atomic → concept|implementation`). This is an independent, pre-existing lint regression, **not** part of the requires-closure story.

---

## Systemic themes (the deep patterns)

1. **The requires-closure is a second-class citizen.** It is resolved at exactly one site (`activate.sh`) and materialized only as side-effects (symlinks + manifest rows). It is **not** a first-class array that status renderers, inject guards, or eject guards can consult. Every symptom in Cluster A (and the F1 blocker) is this one omission. The text renderer recomputes it locally and then *throws it away*. **This is the single highest-leverage architectural fix.** It is also a textbook case of the **orphaned-logic** anti-pattern the user worries about: logic that exists (closure resolution) but whose result is dangling — visible to humans, invisible to machines.

2. **Silent-success-on-absence hides bugs.** `_render_tree`'s `[[ -d ]] || return 0` (F5) and the `find ... -empty` with no stderr redirect (F15) both *succeed loudly-wrong or fail-silently*. The pattern "guard returns 0 when the precondition is missing" is convenient but it converts misconfiguration into invisible no-ops. The test suite inherits this (checks dirs exist, not that content rendered).

3. **No reconciliation between live state and the evolving catalog.** Roles get deleted (v1.7.0 dropped 6); README and *live project manifests* still point at them (F10, F11). There is no `validate`-against-a-manifest, no `doctor`, no migration. State drifts and nothing notices. This is the **mess-from-previous-sessions** category made concrete.

4. **The map from "science" to "harness" is thinner than the mechanics.** Every agent rated the *machinery* 8–10/10 and the *navigation* 3–6/10. `list skills` names-only, no `list tags`, `pathway-signature` ⊂ base, no dedicated scRNA-integration+annotation lane. The CLI can *do* everything but can't *tell you what to do*.

5. **Bash legibility — mostly healthy, one watch item.** `status.sh` is the largest module and concentrates four parallel renderers (text/effective/json/source) that each re-walk the same state slightly differently — this divergence is exactly what produced F2 (text counts inherited, the other three don't). It's not yet an n^n interdependency mess, but it is the one file where "many namespaced blocks doing the almost-same-thing" is starting to drift. Consolidating the four renderers onto one shared state struct would both fix F2 and arrest the drift. Elsewhere modules are small, single-purpose, use `return` not `exit` (test-enforced), and inject/eject deliberately mirror each other with an in-code "must stay in sync" note (eject.sh:210) — good discipline.

**Drift verdict:** Not currently drifting into n^n interdependency. The one real anti-pattern present is **orphaned logic** (the discarded closure result), and the one watch-item is the **four-renderer divergence in status.sh**. Address those two and the architecture stays clean.

---

## Prioritized remediation (max trust-per-effort first)

| Rank | Fix | Effort | Why first | Kind |
|------|-----|--------|-----------|------|
| 1 | **F1 BLOCKER:** make eject preserve a symlink that activate.sh also owns (or refuse to eject closure-mounted skills) | M | Prevents silent mid-session destruction of a live dependency — the only data-loss bug | genuine bug |
| 2 | **F6 CI BLOCKER:** `scope: atomic → concept` in mllmcelltype SKILL.md | XS | One line unblocks green CI | pre-existing |
| 3 | **F2 + F4 (one fix):** hoist `INHERITED_SKILLS` to a global in `_status_load_state`; feed `--effective`/`--json`/`--source`; make inject guard consult it (warn/refuse) | M | Collapses 4 symptoms into one change; closes the F1 trigger; makes machine-readable status truthful | genuine bug |
| 4 | **F7:** guard the empty-array loop / `return 0` in `_status_render_effective` | XS | Unbreaks the most-scripted status mode | genuine bug |
| 5 | **F5:** rename `templates/project/software/ → software-tool/` (or add `$type→dir` map) + a test that asserts type-template *content* rendered | S | Finishes the in-progress rename; software-tool projects become usable | user's in-progress work |
| 6 | **F3:** inject resolves transitive closure (or warns which deps are missing) | M | Orchestrator injects stop silently leaving leaf skills absent | genuine bug |
| 7 | **F15:** `2>/dev/null` on `_seed_gitkeep`'s find; README show explicit `<dir>` | XS | Stops stdout pollution + /tmp spray | genuine bug |
| 8 | **F9/F10/F18:** README — drop `roster`, replace `scrna-atlas`/`planning`, fix `context.md`→`scientific-context.md`, resolve `<scbio-docker>` | S | Stale docs that error on first run | mess from prev session |
| 9 | **F11:** add `sciagent doctor` (or have `status`/`validate` flag a manifest stack referencing a deleted role) + a migration note | M | No path today for drifted live manifests | hanging thread |
| 10 | **F12/F13:** add `sciagent list tags`; enrich `list skills` with one-line desc; decide pathway-signature standalone-vs-overlay and document | M | Moves science→skill from "read the source" to "ask the tool" | design |
| 11 | **F8/F14/F16/F19/F20/F21:** UX polish (overlay-swap guidance, splice "Active role", deactivate-_injected semantics, last-reviewed, message nits) | S each | Low-severity, do opportunistically | mixed |

**Do-not-touch (intentional, per reports 04/06):** the 25 dormant inject-only skills (à-la-carte design), the architect/architecture-treemap allowlisted name collisions, `architecture-treemap` PROBATIONARY sunset, `_archive`/`.deprecated` segregation.

---

## Blank spots & hanging threads (emerged across reports, owned by no single path)

- **Orphaned closure result (the master blank spot).** The inherited-skills array is computed and discarded; no machine surface sees it. Manifests on disk (43 symlinks) disagree with `status --json` (40). Any external consumer is silently misled. (Surfaced from 3 angles: 01 json-count, 02 effective+source, 03 inject/eject.)
- **Live manifest ↔ catalog drift with no reconciler.** Real `DC_hum_verse` manifest pins deleted role `planning`; no `doctor`/migration/validate-manifest exists. Only agent 06 found it, and only because it grounded on a real project — **no tracer path was designed to catch orphaned *live* state.** This is the highest-value untested area.
- **Dead role refs in shipped docs.** README points at `roster` (never existed) and `scrna-atlas`/`planning` (deleted). Docs are not validated against the live catalog.
- **No science→skill discovery surface.** `inject --tag` exists but `list tags` doesn't; `list skills` is names-only; `tags.yaml` descriptions are CLI-invisible. The 25 dormant skills are undiscoverable in practice.
- **Manifest↔filesystem drift is undetectable.** After the F1 eject bug, the manifest and disk diverge; only a full deactivate/reactivate repairs it. Neither `status` nor `validate` checks for or repairs drift (agent 03 #11). A `status` symlink-integrity check exists (status.sh:94-102) but only flags *broken* symlinks, not *missing-vs-manifest* divergence after an over-eager eject already pruned the manifest row.
- **Test blind spots:** `test_new_project_types.sh` asserts directories, not rendered template *content* (missed F5); no test exercises eject-of-inherited (missed F1); no test runs `status --effective` on a no-injected-commands stack (missed F7).

### Where agents disagreed / claims to treat with caution
- Report **04** rated the whole toolkit "healthy, zero fatal inconsistencies" and called software-tool naming "Resolved." That assessment is **too generous**: it ran no inject/eject and no `--effective`/`--source`, so it structurally could not see the Cluster A bugs or F5's content-level breakage. Its *enumeration* (counts, orphan list, requires graph) is solid and corroborated; its *verdict* is path-limited.
- Report **05** labeled F6 a "BLOCKER-CLASS BUG." Accurate for CI-red, but it is a **lint/vocabulary** regression, not a logic bug — distinct in nature from the F1 runtime blocker. Both are real; don't conflate them.
- **F8** (1-arg activate nukes base) is behavior-confirmed by agent 02 but I did not re-run it against source; the README does say "replaces current stack," so it is arguably *documented* behavior — classify CONFUSING (footgun), not BUG.
- **F13** (pathway-signature ⊂ base) — the subset relation is confirmed from the role listings (all 8 ⊂ 36); whether that's "wrong" is a design call, not a defect.
