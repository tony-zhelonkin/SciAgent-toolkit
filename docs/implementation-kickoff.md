# implementation-kickoff — sciagent extension (ARCHIVED — superseded by the three-verb toolkit)

> **STATUS: HISTORICAL RECORD. DO NOT EXECUTE. DO NOT PASTE §7 INTO A SESSION.**
>
> This document planned an *expansion* of the role layer: a `metadata.requires:`
> dependency graph, a `tags.yaml` vocabulary + tag-based `inject --tag`/`eject`
> verb pair, and a `scope: concept|implementation` taxonomy. None of that
> shipped as designed. Instead, the "demolish the role layer" refactor (see
> `docs/architecture.md`, `CHANGELOG.md`) went the opposite direction: it
> **removed** `inject`/`eject` entirely, **removed** the `metadata:`
> frontmatter block and every `requires:`/tag mechanism, and made roles
> pure provenance labels — the full skill/agent/command catalog is mounted
> unconditionally now, so there is nothing left to inject, eject, or gate by
> tag. That checking is `scio lint --check toolkit` today — see the Lint section
> of `docs/architecture.md` — and not the graph-and-tag validator PR 2 below
> designs. The `sciagent` command name and the `validate` verb are both retired.
>
> Every PR plan, file list, code skeleton, and the session-resume prompt in §7
> below describes that abandoned direction. It is kept **only as an audit
> trail** of what was once decided and why the team changed course — read it
> as history, never as a task list. If you are an agent and were pointed at
> this file to "pick up where it left off," stop and re-read
> `docs/architecture.md` instead; nothing here should be implemented.
>
> The imperative mood throughout ("do not soften `skill_deps.sh`", "execute
> PR 1", etc.) is preserved verbatim from the original working document for
> historical accuracy — it does not apply anymore.

---

Written as a working document for executing the sciagent extension after grilling, mirroring the dual purpose of `docs/kickoff.md`:

1. **A guide for the human (Anton)** — what was to be built, in what order, what each PR had to verify, when to stop and ask.
2. **A session-resume prompt for a fresh Claude Code session** — §7 was meant to be pasted into a clean session to pick up the work. (Superseded — do not paste it; see the banner at the top of this file.)

Same content served both because the implementer's reading order and the reviewer's verification order were the same list.

The grilling phase that produced this plan was done at the time. `kickoff.md` §3 was empty; §9 held the 16 resolutions; `docs/proposals/deffered.md` held four wholesale deferrals. The grilling document stays as audit trail; this document was meant to drive PRs, but the plan it describes was abandoned in favor of the opposite direction (see the banner at the top of this file).

---

## 1. Purpose & dual use

This document does not introduce new architectural decisions. Every "should" below traces back to `kickoff.md` §9 or `architecture.md`. When the executing assistant or reviewer finds a real ambiguity, it surfaces as a sub-bullet labelled either **"Decision needed during implementation"** (architectural weight — surface to Anton, do not invent) or **"Implementer-level micro-decisions"** (take the recommendation and proceed; surface at PR review only if uncertain). Those sub-bullets are the only places this document defers downstream; everything else is a transcription of decisions already made.

The reading order in §6 is identical for the human reviewer and for a future assistant: load the same nine files in the same order, then start at PR 1.

---

## 2. What's in scope (in-flight)

The authoritative decision list is `docs/kickoff.md` §9 (entries dated 2026-05-24). Sixteen resolutions land across three ADRs. Compressed one-liners below for verification purposes only — the implementer must read §9 directly before executing each PR.

| # | Decision (kickoff.md §9 anchor) | One-line summary |
|---|---|---|
| 1 | ADR-001 metadata volume | Keep all five fields: `scope`, `requires`, `complementary-skills`, `contraindications`, `tags`. |
| 2 | ADR-001 migration approach | Big-bang codemod in one PR; defaults filled by heuristic; reviewed in same PR. |
| 3 | ADR-001 tag vocabulary location | Single `tags.yaml` at toolkit root; fields `name`, `description`, `since`. |
| 4 | ADR-001 tag vocabulary mutability | PR-only; no `sciagent tags add` verb. |
| 5 | ADR-001 tag vocabulary seed | 10 tags: `trajectory`, `integration`, `annotation`, `de`, `pathway`, `qc`, `viz`, `report`, `architecture`, `tooling`. |
| 6 | ADR-001 inject-by-tag mechanism | Extends existing `inject` verb with `--tag <name>` form. |
| 7 | ADR-001 inject-by-tag × stack-cap | Tag-injection is a separate concept from stack; depth-2 cap stays. |
| 8 | ADR-001 eject verb | Symmetric pair to `inject`: `eject <skill>` and `eject --tag <name>`. |
| 9 | ADR-002 `requires` contract | Warn-and-continue (exit 0); never block activation. |
| 10 | ADR-002 warning surface | STDERR end-of-activate summary block; no per-line noise, no JSON sidecar. |
| 11 | ADR-002 semver on `requires` | Plain names; no version syntax. |
| 12 | ADR-003 scope vocabulary | `concept` / `implementation` (replaces existing `atomic` / `orchestrator` / `foundation`). |
| 13 | ADR-007 validator | Dissolved into a thin `sciagent validate` verb (~30 LOC of checks); `activate` calls it internally. |
| 14 | ADR-008 `lab-loop` | Dropped (depended on deferred ADR-004); see `deffered.md`. |
| 15 | ADR-009 `uv` | Deferred (Docker workflow); see `deffered.md`. |
| 16 | ADR-010 cloud-lab shape | No change. |

Items 1, 6, 8 expand the surface of an existing verb or schema; items 12, 13 are net-new vocabulary or verbs. The codemod (item 2) touches every skill in the repo at once.

---

## 3. What's NOT in scope (do not reopen)

Per `docs/proposals/deffered.md`. These are off the table for this implementation batch:

| Deferral | What it means in code |
|---|---|
| **ADR-005** (trace recording + workflow-skill-creator) | No JSONL recorder. No interview-driven distiller. No conversation-trace infrastructure of any kind. A future docs/code → skill(s) generator is placeholdered as a separate ADR — *not* part of this implementation batch. |
| **ADR-004** (ERA optimizer) + **ADR-006** (benchmark harness) | No `scorable` overlay. No `mutator`/`scorer`/`selector` sub-agents. No FUTS algorithm port. No tier-1/tier-2 benchmark suite. No fitness signal. The sub-agent bisection in `kickoff.md` §4 is illustrative only — moot under this deferral. |
| **ADR-008** (`lab-loop` overlay) | Dropped. Reactivates only when ADR-004 reactivates. Not authored standalone in this batch. |
| **ADR-009** (`uv` pinning) | Not in scope. Skill Python deps stay author-managed. No toolkit-level enforcement. |

If during implementation a tempting "while we're in here…" change touches any of these surfaces, that is scope creep. Push it back to a future ADR and proceed with the planned PR.

---

## 4. PR breakdown

Five PRs, ordered by dependency. PRs 1–3 are the substantive code work; PR 4 closes test coverage; PR 5 is a placeholder for the docs/code → skill generator (separate ADR, not necessarily this batch).

```
PR 1 ───┬──► PR 2 ───► PR 3 ───► PR 4
        │
        └──► (PR 5 standalone, future ADR)
```

| PR | Title | Depends on | Files touched (representative) | Acceptance |
|----|---|---|---|---|
| 1 | ADR-001/003 metadata codemod + `tags.yaml` seed | — | `tags.yaml` (new), `skills/*/SKILL.md` (mass), `templates/skill/SKILL.md`, `tests/test_skill_scope_lint.sh` | All skills carry the five fields; `tags.yaml` present with 10 seeds; lint test updated to new scope vocabulary; `bin/sciagent activate base` green. |
| 2 | `sciagent validate` verb + internal call from `activate` + STDERR warning surface | PR 1 | `bin/sciagent`, `lib/sciagent/validate.sh` (new), `lib/sciagent/activate.sh`, `tests/test_validate_*.sh` (new) | `sciagent validate` exits non-zero on cycle/missing/unknown-tag; `activate` continues on missing `requires` with summary block to STDERR; standalone verb works for debugging. |
| 3 | `sciagent eject` verb + `inject --tag` / `eject --tag` | PR 1, PR 2 | `bin/sciagent`, `lib/sciagent/inject.sh`, `lib/sciagent/eject.sh` (new), `lib/sciagent/symlinks.sh` (manifest entry), `tests/test_inject_tag_*.sh`, `tests/test_eject_*.sh` (new) | Round-trip: `inject X` then `eject X` is a no-op final state; `inject --tag pathway` mounts all skills with that tag; `eject --tag` removes them; `eject` of a stack-mounted skill errors with the right pointer to `deactivate`. |
| 4 | Consolidated test coverage pass | PR 1–3 | `tests/test_*.sh` | All three new verbs/forms have positive, negative, and idempotency tests. `bash tests/run-all.sh` green. |
| 5 | (future ADR) docs/code → skill(s) generator | — | TBD | Out of scope here. Open as a fresh ADR and grill before designing. |

PR 3 depends on PR 1 because `--tag` requires the `tags.yaml` vocabulary and skill `tags:` fields. PR 3 depends on PR 2 because tag-vocab compliance is one of the `validate` checks that `activate` will run before any `--tag` lookup. PR 4 may be folded into PRs 1–3 if it shortens the per-PR diffs; left standalone here for clarity.

---

## 5. Per-PR detailed plan

### PR 1 — Metadata codemod + `tags.yaml` seed

**Scope.** Mechanical, repo-wide. Adds the five `metadata.sciagent.*` (or `metadata.*`) fields to every skill at once and seeds the tag vocabulary file. No CLI changes.

**Files.**

- **New:** `tags.yaml` at toolkit root.
- **Modified (mass):** every `skills/<name>/SKILL.md` that does not already carry all five fields. Per the audit, `scope`, `requires`, `complementary-skills`, `contraindications`, and `tags` already exist on most curated skills; the codemod normalises the ones that are missing fields (e.g. `architecture-first-dev`, `skill-creator` carry only `scope` + `requires`).
- **Modified:** `templates/skill/SKILL.md` — bring the template in line with the canonical shape (the placeholders are already there; align field names with §9's decided vocabulary).
- **Modified:** `tests/test_skill_scope_lint.sh` — change scope vocabulary AND cap values per `kickoff.md` §9, 2026-05-24 cap decision. Specifically: replace the three-way `atomic`/`orchestrator`/`foundation` switch (lines ~79–88) with `concept` ≤ 500 / `implementation` ≤ 350; change default-on-missing from `atomic` to `implementation` (line ~75); bump CUTOFF to `2026-05-24` (line ~19); add a regression-block rule that FAILs on any active skill carrying `scope: atomic|orchestrator|foundation`.

**Codemod behaviour.** One Python (or bash + `yq`) script committed under `scripts/` or run-once and discarded — author's choice. Per skill, it:

1. Parses the frontmatter block.
2. Inserts any missing field under `metadata:` with an empty/default value.
3. Maps existing `scope` values: `foundation` and `orchestrator` → `concept`; `atomic` → `implementation` (per `kickoff.md` §9, 2026-05-24 PR 1 scope-mapping decision). `architecture-treemap` is missing the `scope:` field; assign `implementation` during the codemod. Hand-flip individual skills during code review where the default mapping looks wrong on inspection.
4. Leaves the body untouched.
5. Re-serialises preserving comments where the YAML library supports it; otherwise the diff is reviewed by hand for any prose loss.

**`tags.yaml` content (per kickoff.md §9, 2026-05-24 — "ADR-001 tag vocabulary seed").** Ship the file verbatim from §9:

```yaml
# tags.yaml — sciagent skill tag vocabulary.
# Single source of truth. Adding a tag = edit + PR.
# See docs/proposals/sciagent-extension-design-spec.md ADR-001.

tags:
  - name: trajectory
    description: Time-trajectory inference (pseudotime, RNA velocity, dynamics).
    since: 2026-05-24

  - name: integration
    description: Batch correction and multi-sample harmonization (scVI, scANVI, harmony, MNN).
    since: 2026-05-24

  - name: annotation
    description: Cell-type assignment (marker-based, reference-projection, manual rescue).
    since: 2026-05-24

  - name: de
    description: Differential expression (limma, edgeR, DESeq2, scanpy.rank_genes_groups).
    since: 2026-05-24

  - name: pathway
    description: Pathway analysis — enrichment (fgsea, hypergeometric), activity scoring (GSVA, AUCell, decoupleR), signatures.
    since: 2026-05-24

  - name: qc
    description: Quality control (per-cell, per-sample, batch-level).
    since: 2026-05-24

  - name: viz
    description: Plotting, dimensionality-reduction display, figure composition.
    since: 2026-05-24

  - name: report
    description: Knit/render, publication-figure pipelines, captions, README/Methods generation, skill documentation.
    since: 2026-05-24

  - name: architecture
    description: Software architecture — design patterns, refactoring guidance, system structure.
    since: 2026-05-24

  - name: tooling
    description: General-purpose tool/utility building (CLIs, MCP servers, file-handling helpers, dev workflows).
    since: 2026-05-24
```

**Tests.** Update `tests/test_skill_scope_lint.sh` per the cap decision in `kickoff.md` §9, 2026-05-24:
- `concept` ≤ **500** body lines (replaces the legacy `foundation ≤ 800` / `orchestrator ≤ 250`).
- `implementation` ≤ **350** body lines (replaces `atomic ≤ 300`).
- Default scope on missing frontmatter field: `implementation` (was `atomic`).
- CUTOFF bump: `2026-05-21` → `2026-05-24`.
- New rule: presence of `scope: atomic|orchestrator|foundation` in any active skill is a FAIL (regression block).
- Existing `body_loc_no_fences` counting machinery is unchanged.

The 500/350 numbers catch 100% of today's 62 skills with margin (concept p90=235 → 2.1× headroom; implementation p90=224 → 1.6× headroom). No hand-flips or content refactors are required for any existing skill; `iterative-peak-merging` (317 body lines, 0 companions) was the only would-be casualty of a tighter 300 cap.

Folder-companion file pattern (skill folder may carry reference files referenced by basename from the SKILL.md body) is blessed for **both** scopes — it codifies the existing convention used by 12/62 skills.

Add a new test `tests/test_tags_vocabulary.sh` that verifies every `metadata.tags:` entry across all skills exists in `tags.yaml`.

**Acceptance criteria.**

- `bash tests/run-all.sh` passes.
- `bin/sciagent activate base` succeeds on the host project (`DC_hum_verse`) with no warnings about missing metadata.
- `tags.yaml` parses (valid YAML, ten entries, one entry per line of the seed list above).
- Every existing `metadata.tags:` value references a name in `tags.yaml`.

**Risks.**

- The codemod over ~60 skills may produce YAML that mangles existing comments or quoting. Mitigation: run on a branch; eyeball diff per skill; have the test suite green at the end.
- Scope rename collides with the active `test_skill_scope_lint.sh` semantics. Mitigation: update the lint test in the same PR; do not split.

**All PR 1 decisions resolved 2026-05-24 (see `kickoff.md` §9):**

- ✅ **Scope mapping**: `foundation`/`orchestrator` → `concept`; `atomic` → `implementation`. Hand-flips during code review for skills where the mapping looks wrong on inspection.
- ✅ **Caps**: `concept` ≤ 500 body lines; `implementation` ≤ 350 body lines. Catches 100% of today's 62 skills with margin; no forced refactors.
- ✅ **Folder-relief pattern**: blessed for both scopes (codifies the existing convention used by 12/62 skills).
- ✅ **Default-on-missing**: `implementation` (was `atomic`).
- ✅ **CUTOFF bump**: 2026-05-21 → 2026-05-24.
- ✅ **Regression-block lint rule**: presence of `scope: atomic|orchestrator|foundation` in any active skill = FAIL.
- ✅ **`architecture-treemap`**: assigned `scope: implementation` during the codemod (currently missing the field).

No PR 1 architectural questions remain open. Implementation can proceed.

### PR 2 — `sciagent validate` + activation-time warning surface

**Scope.** A new dispatcher verb that runs cheap graph checks. `activate` calls it internally as a pre-flight; failures of certain checks (cycle, missing required skill — when wired through the existing `skill_deps.sh` resolver) already block activation today, while a *new* class (missing `requires` declared on indirect deps that are themselves absent — see §9 item 9: "Warn-and-continue, visibly") triggers the new STDERR end-of-activate summary block.

**Files.**

- **Modified:** `bin/sciagent` — add a `validate` case branch sourcing `lib/sciagent/validate.sh`.
- **New:** `lib/sciagent/validate.sh` — the verb implementation.
- **Modified:** `lib/sciagent/activate.sh` — call `validate` at the start; defer the STDERR summary print to the end (after the existing activation summary).
- **New tests:** `tests/test_validate_cycle.sh`, `tests/test_validate_missing_requires.sh`, `tests/test_validate_unknown_tag.sh`, `tests/test_activate_warn_on_missing.sh`.

**`lib/sciagent/validate.sh` shape.** Roughly:

```bash
# lib/sciagent/validate.sh — sciagent validate [--quiet]
# Cheap graph checks over skills/ and roles/.
#
# Checks (per kickoff.md §9, ADR-007 resolution 2026-05-24):
#   1. requires resolution    — every named dep exists as a skill
#   2. cycle detection         — DFS over requires graph (skill_deps.sh helper)
#   3. tag-vocab compliance    — every metadata.tags: entry is in tags.yaml
#   4. (optional) skills-ref shell-out, if installed; silently skip otherwise
#
# Exit code:
#   0  — all checks pass
#   1  — any hard check fails (cycle / missing skill / unknown tag)
#
# Called by `sciagent activate` before mounting. Also runnable standalone
# for debugging: `sciagent validate` from the toolkit root.

cmd_validate() {
    local quiet=0
    [[ "${1:-}" == "--quiet" ]] && quiet=1

    local fail=0
    # Check 1+2: walk every skill, attempt transitive resolution.
    # Reuses skill_deps.sh:skill_resolve_transitive which already
    # detects cycles and missing targets.
    # Check 3: parse tags.yaml; verify every metadata.tags entry
    # across skills/ exists in the vocabulary.
    # Check 4 (optional): if `command -v skills-ref` >/dev/null,
    # invoke per skill and surface its exit code.

    # Print a one-line OK summary on success (silenced by --quiet);
    # print a per-failure stanza on failure.
    return $fail
}
```

**Activation-time warning surface (per kickoff.md §9, 2026-05-24 ADR-002 `requires` warning surface).** The block goes to STDERR after the "Activated stack:" lines `activate.sh` already prints. Proposed format:

```
sciagent: activation completed with warnings:

  requires (missing skills):
    - skill-A requires `tool-x`  — not present in skills/
    - skill-B requires `helper-y` — not present in skills/

  Action: install the missing skill(s), drop the `requires:` entries,
  or proceed if the gap is non-breaking for your current work.
```

`activate` already returns 1 on a hard `skill_resolve_transitive` failure (per `activate.sh` Phase B). Per the PR 2 hardness boundary decision in `kickoff.md` §9 (2026-05-24), this is **preserved**: `requires` failures continue to abort before any symlink lands. The new soft-warn applies **only** to `complementary-skills` references that don't resolve — those are buffered during the walk and dumped to STDERR at end-of-activate; exit code stays 0.

**Tests.**

- `test_validate_cycle.sh` — builds a fake toolkit with a `requires:` cycle; asserts `sciagent validate` exits 1 and message names the cycle.
- `test_validate_missing_requires.sh` — `requires: [does-not-exist]`; asserts `sciagent validate` exits 1.
- `test_validate_unknown_tag.sh` — skill with `tags: [bogus]` not in `tags.yaml`; asserts validate exits 1 and names the offending tag + skill.
- `test_activate_warn_on_missing.sh` — soft-warn case (whatever the implementer decides counts as soft); asserts `activate` exit 0 with summary on STDERR.

**Acceptance criteria.**

- `sciagent validate` works standalone from a project root.
- `sciagent activate base` calls validate internally; a misconfigured skill graph aborts before any symlink lands (preserving the existing "pre-mutation contract" from `activate.sh` Phase B).
- Soft warnings appear at end-of-activate on STDERR, formatted, with the exit code remaining 0.

**Risks.**

- The current `activate.sh` already runs `skill_resolve_transitive`. Distinguishing "this is hard-fail in `validate`" from "this is hard-fail in `activate`" is a subtle line. Mitigation: be explicit that `validate` and `activate` share the same resolver code path; the only thing layered on top is the *new* checks (tag-vocab compliance, optional `skills-ref` shell-out) and the STDERR summary buffering.
- The "optional shell-out to `skills-ref`" can produce noisy output for skills that fail upstream schema. Mitigation: silently skip when not installed; surface but don't fail on each finding (warn, not error).

**Decisions resolved 2026-05-24 (see `kickoff.md` §9):**

- ✅ **Hardness boundary**: keep resolver hard-fail on `requires` (preserves the pre-mutation safety invariant); soft-warn only for `complementary-skills` references. STDERR end-of-activate summary block carries only the complementary-skill warnings; `requires` failures still abort.

**Decision still needed during implementation:**

- **STDERR summary exact wording.** The block above is a proposed format. §9 says "formatted, human-readable summary printed at the end of `activate` to STDERR. One place to look" — no exact wording. Implementer's call at PR review time.

### PR 3 — `inject --tag` / `eject` verb pair

**Scope.** Two changes that hang together:

1. Extend `inject` to accept `--tag <name>` in addition to the existing positional skill name.
2. Add a symmetric `eject` verb accepting positional skill name OR `--tag <name>`.

Stack-mounted items (via a role) remain owned by `activate`/`deactivate`. `eject` operates exclusively on items previously placed by `inject` — surfaced via the manifest's `injected[]` array.

**Files.**

- **Modified:** `bin/sciagent` — add `eject` case branch.
- **Modified:** `lib/sciagent/inject.sh` — handle `--tag <name>` form. Resolution: read `tags.yaml`; collect every skill in `skills/*/SKILL.md` whose `metadata.tags:` contains `<name>`; inject each.
- **New:** `lib/sciagent/eject.sh` — mirrors `inject.sh`'s shape.
- **Modified:** `lib/sciagent/symlinks.sh` — the manifest's `injected[]` schema gains a "source" marker (e.g. `{"overlay": "...", "skill": "...", "via": "tag:pathway"}`) so that `eject --tag pathway` knows which skills it owns vs. which were named individually. See "Implementer-level micro-decisions" below for the schema shape recommendation.
- **New tests:** `tests/test_inject_tag_creates_overlay.sh`, `tests/test_eject_named.sh`, `tests/test_eject_tag.sh`, `tests/test_eject_rejects_stack_mounted.sh`, `tests/test_inject_eject_roundtrip.sh`.

**`lib/sciagent/eject.sh` skeleton.**

```bash
# lib/sciagent/eject.sh — sciagent eject <skill> | --tag <name>
# Symmetric counterpart to `inject` (per kickoff.md §9, 2026-05-24
# ADR-001 inject↔eject verb pair).
#
# Removes skills previously placed via `sciagent inject`. Skills
# mounted as part of a role-stack layer are NOT eligible — use
# `sciagent deactivate` for that. The two verbs cover orthogonal
# concepts (stack composition vs. ad-hoc injection).

# shellcheck shell=bash

cmd_eject() {
    if ! manifest_exists; then
        echo "no role active; nothing to eject" >&2
        return 1
    fi

    # Parse: either `eject <skill>` or `eject --tag <name>`.
    local mode="" target=""
    case "${1:-}" in
        --tag)
            mode="tag"; target="${2:-}"
            [[ -z "$target" ]] && { echo "usage: sciagent eject --tag <name>" >&2; return 1; }
            ;;
        "")
            echo "usage: sciagent eject <skill> | --tag <name>" >&2; return 1
            ;;
        *)
            mode="skill"; target="$1"
            ;;
    esac

    # Resolve the list of injected entries to remove.
    # - mode=skill: one entry whose .skill == $target
    # - mode=tag:   all entries whose .via == "tag:$target"
    # Guard: refuse to eject a skill that is stack-mounted (defensive;
    # injected[] should not contain stack-owned skills, but tests
    # exist to enforce it).

    # For each: remove the dual symlinks in .claude/ and .agents/;
    # drop the entry from manifest.injected[]; re-render the AGENTS.md
    # block. If injected[] becomes empty AND the synthetic _injected
    # overlay was the only overlay, collapse stack back to [base].

    # Mirror inject.sh's idempotency: ejecting an absent target prints
    # a one-line "not injected" notice and returns 0.
}
```

**`inject --tag` behaviour.** Treat each tag-matched skill as if it had been named explicitly. The resulting manifest carries `{"overlay": "_injected", "skill": "<name>", "via": "tag:pathway"}` for each. `inject --tag pathway` invoked while `pathway` items are already mounted is a no-op (per existing idempotency contract). A tag-injection that mounts zero skills (no skills carry the tag) is a one-line warning, exit 0.

**Tests.**

- `test_inject_tag_creates_overlay.sh` — `inject --tag pathway` on solo base creates `_injected` overlay and mounts every pathway-tagged skill; manifest entries carry `via: tag:pathway`.
- `test_eject_named.sh` — `inject scvi-basic` then `eject scvi-basic` returns the file tree to pre-inject state.
- `test_eject_tag.sh` — `inject --tag pathway` then `eject --tag pathway` removes only pathway-mounted skills, even if other injected skills remain.
- `test_eject_rejects_stack_mounted.sh` — activate base, attempt `eject <skill-from-base>`; assert error with pointer to `sciagent deactivate`.
- `test_inject_eject_roundtrip.sh` — full cycle: `inject X; inject --tag pathway; eject X; eject --tag pathway` ends at pre-inject state.

**Acceptance criteria.**

- `inject --tag <name>` mounts all matching skills and records `via:` in the manifest.
- `eject <skill>` and `eject --tag <name>` mirror their `inject` counterparts.
- Stack-mounted skills are ejection-rejected with the right error message.
- All five new tests pass; existing `test_inject_creates_overlay.sh` and `test_inject_extends_overlay.sh` still pass.

**Risks.**

- Manifest schema bump (`injected[]` entries gain a `via:` field). Mitigation: schema version 1 → 2; emit a one-time migration on read if old entries are present (or document that the migration is forward-only and pre-existing injections require re-injection — likely fine for a tool only Anton uses).
- A tag-injection that overlaps with skills already present in the active stack risks confusing symlink state. Mitigation: the existing `resolve_canonical` + `ln -sfn` flow already handles overwrite; the new code only needs to track that the manifest stays the source of truth for what `eject --tag` will remove.

**Implementer-level micro-decisions** (not blocking; take the recommendation and proceed; surface at PR review only if uncertain):

- **Manifest schema bump.** Does the existing `injected[]` entry shape get a new `via:` field, or does the implementer route via a parallel `injected_by_tag[]` array? Both work. The single-array-with-`via` choice keeps the schema flat and the eject logic linear; the two-array choice avoids touching existing serialisation. **Recommendation**: single-array with `via:`.
- **Empty tag-injection behaviour.** If `inject --tag <name>` matches zero skills, is that an error (typo in tag name?) or a warning? The tag is in the vocabulary, just unused. **Recommendation**: warn + exit 0.

### PR 4 — Consolidated test coverage

**Scope.** Closure pass on the test suite. If PRs 1–3 each ship with their own tests (recommended above), this PR is small — it adds the cross-cutting tests that don't fit cleanly in any single PR:

- `test_validate_skip_skills_ref_when_absent.sh` — confirms graceful no-op when `skills-ref` CLI is not installed.
- `test_activate_calls_validate_internally.sh` — asserts that breaking a skill graph (planted cycle) causes `activate` to fail without writing symlinks (the existing pre-mutation contract).
- `test_tags_yaml_referenced_by_codemod.sh` — sanity check that every skill's `metadata.tags:` is a subset of `tags.yaml` names; complements `test_tags_vocabulary.sh` from PR 1.

**Acceptance.** `bash tests/run-all.sh` green. If PRs 1–3 covered everything, this PR may be empty and skipped.

### PR 5 — (future ADR) docs/code → skill(s) generator

**Scope.** Not part of this implementation batch. Recorded here so it stays visible in the plan. Per `deffered.md` (ADR-005 entry): a single command/agent that ingests documentation and/or a codebase and emits a skill or family of skills. Treat as a fresh ADR with its own grilling pass before designing.

**Action item.** Open a stub ADR placeholder file at `docs/proposals/adr-docs-to-skill-generator.md` listing the contract sketch from `deffered.md`. Do not design yet.

---

## 6. Reading order for a fresh session (historical — the plan below was never executed as designed)

The plan called for loading these files into context, in this order, before executing PR 1:

1. `docs/kickoff.md` §9 — authoritative resolutions. Every PR's "why" lives here.
2. `docs/proposals/deffered.md` — what is explicitly out of scope.
3. `docs/architecture.md` — current system shape; the cap-at-2 stack and symlink topology that PRs 1–4 must not break.
4. `CLAUDE.md` — project conventions; the "After Modifying Agents, Skills, or Roles" workflow that every PR's smoke test follows.
5. `bin/sciagent` — verb dispatcher; locate where `validate` and `eject` cases will land.
6. `lib/sciagent/activate.sh` — Phase B (transitive `requires:`) is the existing seam; PR 2's validate-call slots in here.
7. `lib/sciagent/inject.sh` — model for PR 3's `eject.sh`.
8. `lib/sciagent/skill_deps.sh` — resolver behaviour shared with `validate`.
9. `lib/sciagent/symlinks.sh` — manifest schema (where the new `via:` field lands).
10. `tests/test_inject_creates_overlay.sh` + `tests/test_skill_scope_lint.sh` + `tests/_lib.sh` — test patterns to mirror.
11. `templates/skill/SKILL.md` and 2–3 representative skill frontmatters (`scvi-basic`, `architecture-first-dev`, `skill-creator`) — current shape vs. target shape for PR 1's codemod.

The grilling artefacts (`docs/proposals/ai-research/*.md`, `docs/proposals/sciagent-extension-design-spec.md`) are optional reading; consult only when a "why was this decided" question arises during a PR.

---

## 7. Session-resume prompt (ARCHIVED — do not paste)

> **This prompt is dead.** It was written to resume the PR 1–3 work above.
> That work was superseded by the demolition refactor (see the banner at the
> top of this file). Pasting it into a session today would misdirect an
> agent into rebuilding `inject`/`eject`/`tags.yaml`/`metadata.requires:` —
> exactly the layer that was deliberately removed. Kept verbatim below for
> audit-trail purposes only.

Paste the body below into a fresh Claude Code session. Assumes working directory `/workspaces/DC_hum_verse/01_modules/SciAgent-toolkit/`.

> ---
>
> Resuming the sciagent extension implementation. Grilling phase is complete: 16 decisions in `docs/kickoff.md` §9 (dated 2026-05-24), four wholesale deferrals in `docs/proposals/deffered.md`. This implementation plan lives at `docs/implementation-kickoff.md` — read it now along with the files §6 lists, in that order.
>
> **State of play.** No code changes yet from the grilling phase. Five PRs planned: (1) metadata codemod + `tags.yaml` seed; (2) `sciagent validate` verb + activation-time STDERR warning surface; (3) `inject --tag` and `eject` verbs; (4) test coverage closure; (5) future ADR placeholder for docs/code → skill generator (not in this batch). PRs 1–3 are the substantive work and must ship in order.
>
> **Hard constraints.**
> - Do NOT reopen any deferred ADR: no JSONL trace recorder, no ERA optimizer, no benchmark harness, no `lab-loop`, no `uv` enforcement, no interview-driven distiller. See `docs/proposals/deffered.md`.
> - Do NOT add metadata fields beyond the five listed in `kickoff.md` §9 (`scope`, `requires`, `complementary-skills`, `contraindications`, `tags`).
> - Do NOT build a "validator subsystem." `sciagent validate` is ~30 LOC of cheap checks per `kickoff.md` §9 ADR-007 dissolution.
> - `tags.yaml` is PR-only — do NOT add a `sciagent tags add` CLI verb.
> - Stack-depth cap stays at 2 (base + overlay). `--tag` injection is a separate concept from stack, not a third stack frame.
> - Voice in code comments, READMEs, and any docs written in this session is strictly impersonal. Chat stays conversational.
>
> **All architectural decisions are resolved** in `kickoff.md` §9 (decisions dated 2026-05-24). Do not invent answers on architectural questions — surface and stop. Key resolutions for PR-1-through-PR-3:
> - PR 1 scope mapping: `foundation`/`orchestrator` → `concept`; `atomic` → `implementation`. Hand-flips during code review for individual outliers.
> - PR 1 caps: `concept` ≤ 500, `implementation` ≤ 350. Folder-relief pattern (companion files in skill folder referenced by basename from SKILL.md body) blessed for both scopes. Default scope on missing field = `implementation`. CUTOFF bumped to 2026-05-24. New regression-block lint rejects `scope: atomic|orchestrator|foundation` on active skills.
> - PR 2 hardness boundary: `requires` keeps hard-fail (preserves the pre-mutation safety invariant); `complementary-skills` warn-and-continue with STDERR end-of-activate summary.
>
> **Implementer-level micro-decisions** (sane-default and proceed; surface at PR review if uncertain):
> - PR 2: exact wording of the STDERR end-of-activate warning summary block. §9 specifies "formatted, human-readable summary printed at the end of `activate` to STDERR. One place to look" — exact prose is the implementer's call.
> - PR 3: manifest schema for tagged injections — add `via:` field to existing `injected[]` entries vs. parallel `injected_by_tag[]` array. Recommended: single-array with `via:` (flatter schema, linear `eject` logic). Both work.
> - PR 3: UX when `inject --tag <name>` matches zero skills — warn (exit 0) vs. error (exit 1). Recommended: warn, since an empty tag is a valid future state.
>
> **Today's task.** Start with PR 1. Eyeball the codemod diff per skill (especially `architecture-treemap`, which needs the new `scope: implementation` field added explicitly). Run `bash tests/run-all.sh` before opening the PR.
>
> Confirm context is loaded and propose the codemod approach (Python with `ruamel.yaml` to preserve comments, or bash + `yq` — whichever is cheaper to get right). Then execute PR 1.
>
> ---

---

## 8. Known risks & mitigations

| Risk | Where it bites | Mitigation |
|---|---|---|
| Codemod mangles YAML comments or quoting across ~60 skills | PR 1 | Use a comment-preserving YAML lib (`ruamel.yaml`); review per-skill diffs; require `bash tests/run-all.sh` green before merge. |
| Scope rename collides with active `test_skill_scope_lint.sh` semantics | PR 1 | Update the lint test in the same PR; do not split. Cap values + side items are decided (`kickoff.md` §9, 2026-05-24 cap decision). |
| `validate` and `activate` share the resolver but want different hardness behaviour | PR 2 | Be explicit in `validate.sh` header comment that the boundary is: `requires` → hard-fail (preserved from `skill_deps.sh`); `complementary-skills` → soft-warn (`kickoff.md` §9, PR 2 hardness boundary decision). |
| Manifest schema v1 → v2 for `via:` field breaks pre-existing `.sciagent/manifest.json` files in active projects | PR 3 | Document forward-only migration; on read, treat missing `via:` as a named-skill injection. Re-inject if Anton needs the `via:` data for any prior injection. |
| `inject --tag` mounted skills overlap with role-stack skills | PR 3 | Existing `ln -sfn` is overwrite-safe; manifest carries the truth for `eject --tag`. Add a test for this collision. |
| `skills-ref` shell-out floods STDERR for skills whose schema upstream is in flux | PR 2 | Silent skip when CLI absent; surface as warning, not error, when present. |
| `validate` walks malformed YAMLs and crashes mid-walk, leaving the user without a usable error | PR 2 | Each per-skill check is wrapped — a single skill's parse failure produces a named-skill error in the summary, not an abort. |

---

## 9. What NOT to do

Explicit don'ts. Lifted from §3, `deffered.md`, and the hard constraints in §5/§7. If any of these surface as "wouldn't it be nice if…" during implementation, push back and proceed with the planned PRs:

- Do **not** add JSONL trace recording or any conversation-trace infrastructure (ADR-005 deferred).
- Do **not** add an ERA-style optimizer overlay, `mutator`/`scorer`/`selector` sub-agents, or FUTS port (ADR-004 deferred).
- Do **not** add a benchmark harness, tier-1/tier-2 task suite, or any fitness signal (ADR-006 deferred).
- Do **not** add `lab-loop` or any other "swap-target stub" overlay (ADR-008 dropped).
- Do **not** add `uv`-based dep pinning or any toolkit-level Python dep enforcement (ADR-009 deferred).
- Do **not** add a `sciagent tags add <name>` CLI verb. Tag vocabulary changes are PR-only (`kickoff.md` §9, ADR-001 mutability).
- Do **not** add metadata fields beyond `scope`, `requires`, `complementary-skills`, `contraindications`, `tags`. The five fields are the complete custom-metadata surface (`kickoff.md` §9, ADR-001 metadata volume).
- Do **not** auto-add tags on first use, infer tags from skill description, or otherwise short-circuit the PR-only flow.
- Do **not** build a "validator subsystem." `sciagent validate` is a thin callable verb, ~30 LOC of cheap checks (`kickoff.md` §9, ADR-007 dissolution).
- Do **not** break the depth-2 stack cap. `inject --tag` is a separate concept from stack — never a third stack frame (`kickoff.md` §9, ADR-001 inject-by-tag × stack-cap).
- Do **not** add semver / version-range syntax to `requires:`. Plain skill names only (`kickoff.md` §9, ADR-002 semver decision).
- Do **not** soften the existing `skill_deps.sh` hard-fail on missing `requires`. Hard-fail is preserved — `activate.sh` Phase B's pre-mutation safety invariant stands (`kickoff.md` §9, PR 2 hardness boundary decision). Warn-and-continue applies only to `complementary-skills` references that don't resolve.
- Do **not** add per-line warnings or JSON sidecars for the `complementary-skills` warning surface. One STDERR summary block at end-of-activate (`kickoff.md` §9, ADR-002 warning surface + PR 2 hardness boundary).
- Do **not** rebrand the deferred ADR-005 distiller as the new docs/code → skill generator. The replacement is a fresh ADR with its own grilling (`deffered.md` ADR-005 entry).
- Do **not** edit the original spec at `docs/proposals/sciagent-extension-design-spec.md`. Where it disagrees with `kickoff.md` §9, §9 wins.

---

*End of implementation kickoff. All architectural decisions were resolved in `kickoff.md` §9; deferrals in `deffered.md`. The plan called for walking into PR 1, reading §5, and executing — but the plan itself was superseded (see the banner at the top of this file); do not execute it.*
