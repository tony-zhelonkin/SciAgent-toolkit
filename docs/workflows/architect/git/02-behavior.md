---
parent: ./README.md
view: behavioral
---

# Behavior — git automation

Behavioral view: when commits fire, what each commit looks like, how
branches and rollback behave, what fails where.

---

## The universal rule

> **Commit after the human gate. Never before. Never bypass.**

Every commit is a bookmark for "the user just approved this artifact".
Rejection is handled at the working-tree layer (`git restore`) — it never
becomes history. Pre-commit hooks are respected; the harness never uses
`--no-verify`, `--amend`, `--no-gpg-sign`, or `-f`.

---

## Per-command commit flow

Each block below describes: (a) the gate that must pass, (b) the subject
template, (c) the paths staged, (d) what happens on rejection.

### `/map {feat}`

- **Gate:** Phase 3 user confirmation (use / regenerate / append).
  If the user picked "use existing", no write happened; no commit.
- **Subject:** `docs({feat}): map codebase ({loc} locs across {files} files)`
  - `{loc}`, `{files}` come from the map summary the command already prints.
- **Body** (optional): one-line open-questions count, e.g.
  `Open questions: 6 (see map.md § Open Questions).`
- **Paths:** `docs/{feat}/map.md`.
- **Rejection:** not applicable — the gate runs *before* writing.

### `/review {feat} --as <set>`

- **Gate:** Phase 4 completion. There is no per-reviewer gate, so once
  all N reviewers return, the set lands atomically.
- **Subject:** `docs({feat}): review [{reviewer_list}]`
  - `{reviewer_list}` = alphabetically-sorted list, e.g. `bioinf, ml, stat`.
- **Body** (optional): one-line headline count per reviewer if present
  in each `review/*.md` frontmatter or top section; otherwise blank.
- **Paths:** every `docs/{feat}/review/*.md` whose mtime is newer than the
  invocation start. (The helper script uses the file list the command
  passes in; mtime is not checked by the script itself.)
- **Rejection:** user can re-run `/review` with a redirect; the prior
  commit stands and the new pass produces a follow-up commit.

### `/synthesize {feat}`

- **Gate:** Phase 3 surface of Headline + P0 list. This is implicit —
  the command writes synthesis.md then reports. Treat the reporting step
  as the gate: if the user chose not to run `/synthesize` (only 1
  reviewer present), nothing happened.
- **Subject:** `docs({feat}): synthesize reviews ({n_convergent} P0, {n_divergent} divergent)`
- **Paths:** `docs/{feat}/synthesis.md`.
- **Rejection:** the user re-runs `/synthesize` with more reviewers; the
  prior synthesis commit stands.

### `/design {feat}` — two gates, two commits

Design is unique: it has two human gates and produces two distinct
artifact sets. One commit per gate.

**Commit A (after Gate 1 — scope approval in Phase 3a):**

- **Subject:** `docs({feat}): design README — scope & acceptance criteria`
- **Paths:** `docs/{feat}/design/README.md` (only this file exists yet).
- **Rationale:** the scope is now pinned. If Gate 2 later falls through,
  the README still represents an agreed scope worth keeping.

**Commit B (after Gate 2 — full bundle approval in Phase 4, post-status-flip):**

- **Subject:** `docs({feat}): design bundle ({n_adrs} ADRs, {n_oqs} OQs) — status APPROVED`
- **Body** (optional): one line per deferred departure, e.g.
  `Deferred: D-013 (opacity-for-degenerate), D-014 (size-channel strategy).`
- **Paths:** `docs/{feat}/design/*.md` (README gets re-staged with the
  status flip; no harm in committing again).

**Commit C (after architect verdict in Phase 5, if verdict = READY):**

- Not a separate commit — `design/review.md` is written by the architect
  agent, then the command reports. We roll architect output into a third
  commit only when the verdict is READY; otherwise surface NEEDS ITERATION
  / NEEDS DISCUSSION and wait (no commit).
- **Subject:** `docs({feat}): architect verdict READY`
- **Paths:** `docs/{feat}/design/review.md`.

**Rejection at either gate:** no commit fires. If the user hand-reverts
with "I changed my mind", `git restore docs/{feat}/design/*.md` reverts
uncommitted edits; the Gate-1 commit (if present) stays as the last
approved state.

### `/architect {feat}` (standalone gate)

- **Gate:** Phase 4 verdict parse.
- **Subject:** `docs({feat}): architect verdict {READY|NEEDS_ITERATION|NEEDS_DISCUSSION}`
- **Paths:** `docs/{feat}/design/review.md`. If the user consented to the
  "flip DRAFT → APPROVED" prompt, also the design/*.md files whose
  frontmatter was flipped.
- **Rejection:** not applicable — architect only reports a verdict.

### `/plan {feat}`

- **Gate:** Phase 4 approval.
- **Subject:** `docs({feat}): plan ({n_phases} phases)`
- **Body** (optional): one line per phase title, e.g.
  `phase-01: config & metadata; phase-02: JS dispatch; …`
- **Paths:** `docs/{feat}/plan/README.md`, `docs/{feat}/plan/phase-*.md`.
- **Rejection:** user edits the plan (hand edit or re-run); the next
  `/plan` invocation writes and commits again.

### `/implement {feat} {N}`

- **Gate:** Phase 4 verification checklist passes. If anything failed,
  **no commit** — the command already halts in Phase 4 and asks the user
  how to proceed. Re-running with the user's fix will re-try verification
  and then commit.
- **Subject:** `feat({feat}): phase {N} — {phase_title}`
  - `{phase_title}` = first-line title from `docs/{feat}/plan/phase-N.md`
    after stripping `# Phase N:`.
  - Use `feat(` — this is real code, not docs. For phases that are pure
    documentation updates (rare but possible, e.g. moving an ADR), use
    `docs(`. The command decides by examining whether any non-`docs/`
    paths were modified.
- **Body:** verification checklist copied as markdown checkboxes,
  pre-filled with the verdicts (`- [x] foo`). This captures "what
  shipped" in-line with the commit.
- **Paths:** every file listed in the phase doc's Files-to-Create and
  Files-to-Modify sections, plus the phase doc itself (because
  `/implement` Phase 3 step 5 marks each file as done in the phase doc).
- **Rejection:** the command already stops at Phase 4 on failure. No
  commit fires. When the user patches and reruns, verification retries
  and the commit happens then.

### `/verify {feat} [phases]`

- **Gate:** none — `/verify` is a read-only drift check that writes
  `verify.md`.
- **Subject:** `docs({feat}): verify — verdict {CLEAN|INCOMPLETE|DRIFT|NEEDS_REVIEW}`
- **Paths:** `docs/{feat}/verify.md`.
- **Note:** committing verify reports is useful for CI or handoff. If a
  user finds the noise excessive, they disable this one command via
  `ARCHITECT_GIT_COMMIT_ON`.

### `/status [slugs]`

- **No commit.** Chat-only output, no file writes.

### `/meta-map [slugs]`

- **Gate:** Phase 4 user surfacing.
- **Subject:** `docs(_meta): map portfolio ({n_features} features, {n_convergent} convergent concerns)`
- **Paths:** `docs/_meta/map.md`.

### `/meta-design [slugs]`

- **Gate:** Phase 4 Gate 1 approval (MADR bundle).
- **Subject:** `docs(_meta): design ({n_madrs} MADRs) — status APPROVED`
- **Body** (optional): per-MADR one-liner, e.g. `MADR-001: shared
  padj-floor convention; MADR-002: opacity priority chain; …`
- **Paths:** `docs/_meta/design.md`, `docs/_meta/deferred.md` if Phase 6.5
  appended rows.

### `/meta-apply [slugs]`

Most intricate case. The compound gate approves some features and rejects
others; we commit the approved ones and `git restore` the rejected ones.

**After Phase 5 option 1 ("approve all OK features"):**

- **One commit per approved feature** (serial, in the order they appeared
  in the compound report).
- **Subject (per feature):** `docs({feat}): apply meta-design ({n_madrs} MADRs)`
- **Body (per feature):** MADR citations, e.g.
  `Applied MADR-001, MADR-003, MADR-006. Deferred rows: D-013..D-014.`
- **Paths (per feature):** the feature's `docs/{feat}/design/*.md` files
  the feature-reviser touched + the feature's slice of
  `docs/_meta/deferred.md` (the script reads the feature's assigned
  ID block and stages only those rows).
- **Second commit for the deferred catalog:** after the per-feature
  commits, one `docs(_meta): deferred — {N} new rows from meta-apply`
  commit collapses any remaining `_meta/deferred.md` diff. In practice
  the per-feature commits already staged each feature's rows; this second
  commit exists only if rows from *all* features were appended in one
  block and per-feature slicing is not clean. The helper script detects
  this and emits one or two commits as appropriate.

**After Phase 5 option 2 ("approve some, reject others"):**

- Commits fire only for approved features.
- Rejected features: `git restore docs/{feat}/design/*.md` and strip the
  feature's deferred rows. Same behavior the command doc already
  describes; now reinforced by the commit-after-gate invariant.

**After Phase 6 architect parallel dispatch:**

- Each architect writes `docs/{feat}/design/review.md`.
- **One combined commit per architect batch:**
  `docs(_meta): architect gate post-meta-apply — {R} READY, {I} NEEDS_ITERATION, {D} NEEDS_DISCUSSION`
- **Paths:** every `docs/{feat}/design/review.md` that was written in
  this batch, plus any frontmatter roll-back edits from Phase 7 (APPROVED
  → DRAFT for NEEDS_ITERATION/NEEDS_DISCUSSION features).

### `/meta-plan [slugs]`

- **Gate:** Phase 4 approval (sequencing).
- **Subject:** `docs(_meta): plan ({n_portfolio_phases} portfolio phases)`
- **Paths:** `docs/_meta/plan.md`.

---

## Commit message conventions

Every subject follows:

```
<type>(<scope>): <imperative summary>
```

- **Type:** `docs` for artifact-only commits, `feat` for `/implement`
  commits touching code. If a phase mixes docs and code, use `feat` — the
  bigger blast radius wins.
- **Scope:** the feature slug, e.g. `color-encoding`, or `_meta` for
  portfolio-level artifacts.
- **Summary:** short, imperative, lowercase except for proper nouns
  (`UMAP`, `MADR`, ADR numbers like `ADR-009`).

Body rules:

- Wrap at 72 columns for readability.
- Bullets with `-` (not `*` or `+`).
- No horizontal rules, no markdown headings, no emojis.
- No trailers (`Co-Authored-By:`, `Signed-off-by:`, `Reviewed-by:`,
  `Generated-with:` …). The commit author is whatever `git config` says;
  nothing else gets appended.

### What must never appear in any commit message

The following strings are explicitly forbidden and must be grep-able on
a post-implementation audit:

```
Co-Authored-By: Claude
Co-Authored-By: Gemini
Co-Authored-By: Codex
Generated with Claude
Generated with Gemini
Generated with Codex
🤖
Claude Code
Anthropic
```

These are not exhaustive; the policy is **no LLM / no provider
references, ever**. See ADR-001 in `03-decisions.md`.

### Worked examples

```
docs(color-encoding): map codebase (69 locs across 5 files)
```

```
docs(color-encoding): review [bioinf, graphic, ml, stat]
```

```
docs(color-encoding): synthesize reviews (7 P0, 3 divergent)
```

```
docs(color-encoding): design README — scope & acceptance criteria
```

```
docs(color-encoding): design bundle (8 ADRs, 2 OQs) — status APPROVED

Deferred: D-013 (opacity-for-degenerate), D-014 (size-channel strategy).
```

```
docs(color-encoding): architect verdict READY
```

```
docs(color-encoding): plan (6 phases)

phase-01: config & metadata foundation
phase-02: JS dispatch wiring
phase-03: CSS rewrite + legend
phase-04: dot-count rework
phase-05: dirColor runtime resolution
phase-06: literal-forbidden assertions
```

```
feat(color-encoding): phase 1 — config & metadata foundation

- [x] COLORS_DIVERGING_LUT has length 101
- [x] ENTITY_PROFILES['Pathway'] has 11 keys
- [x] ENTITY_PROFILES['TF']['native_score_axis_max'] == 15.0
- [x] generate_dashboard emits entity_profiles / palette_lut / score_saturation
- [x] COLORS_DIVERGING (9-stop) still present for CSS fallback
- [x] ENTITY_TYPES/ENTITY_SHAPES still work for legacy callers
```

```
docs(_meta): design (4 MADRs) — status APPROVED

MADR-001: shared padj-floor convention across features
MADR-002: opacity priority chain (entity-level then score-level)
MADR-003: portfolio-coupled gene-set parse for umap + biological-workflow
MADR-004: ENTITY_PROFILES schema shared between color-encoding and biological-workflow
```

```
docs(biological-workflow): apply meta-design (3 MADRs)

Applied MADR-001, MADR-002, MADR-004.
ADR-007 revised (opacity priority chain per MADR-002).
ADR-011 added (ENTITY_PROFILES consumption per MADR-004).
Deferred rows: D-093..D-099.
```

---

## Branch behavior

### Default (`ARCHITECT_GIT_BRANCHES=false`)

- Everything commits to the current branch.
- The user decides when to branch, push, or merge.
- Works for solo research workflows where `main` is the working surface.

### Opt-in (`ARCHITECT_GIT_BRANCHES=true`)

- On the **first** commit for a feature slug, the helper script creates
  `${ARCHITECT_GIT_BRANCH_PREFIX}{slug}` (default: `architect/{slug}`) if
  it does not exist, and `git switch` to it.
- Subsequent commits for the same slug commit to that branch.
- Branch creation uses `git switch -c` from the current HEAD. If
  uncommitted changes exist, the switch fails; the script surfaces the
  error and exits 5. The user stashes/commits manually and re-runs.
- The harness never deletes, rebases, or merges a branch. Users merge
  when they choose, by whatever tool they prefer (git, PR, etc.).

### Meta-layer branches

- If `ARCHITECT_GIT_META_BRANCH=current`, meta commands commit wherever
  the user is. This is correct when meta work rides the same branch as
  the features it touches.
- If `ARCHITECT_GIT_META_BRANCH=portfolio`, meta commands commit on
  `${ARCHITECT_GIT_BRANCH_PREFIX}_meta` (default: `architect/_meta`). This
  gives portfolio-level artifacts a clean lane when many feature
  branches are in flight. Crossing the streams (`/meta-apply` editing a
  feature's `design/*.md`) in this mode requires the user to be on the
  feature branch first; the script refuses to commit a feature's
  `docs/{feat}/design/*.md` onto `architect/_meta` and exits 5 with a
  clear error.

### Worktrees

**Not supported.** Parallel agent dispatch writes to distinct paths;
worktrees add checkout overhead and force the user to juggle working
directories. Users who want them can `git worktree add` manually and run
the harness inside that worktree — the helper script will honour whatever
working tree it is invoked in.

---

## Rollback behavior

### Before a commit (working tree only)

- `git restore` on the artifact paths reverts the uncommitted edits. This
  is what `/meta-apply` Phase 5 option 2 already does; the
  commit-after-gate rule preserves that contract.

### After a commit (history exists)

The harness does not rewrite history. User options:

- `git revert <sha>` — creates a new commit that inverts the one to
  undo. History stays linear, no force-push needed. This is the
  recommended path.
- `git reset --hard <sha-before>` on an unpushed branch — fast but
  destructive; user's call.
- `/architect-undo` (not proposed here) could be added later as a
  convenience wrapper around `git revert HEAD`. Out of scope for v1.

### Pre-commit hook rejection

- Exit code 5 from the helper script.
- Calling command surfaces the hook's output verbatim to the user.
- No retry, no `--no-verify`. The user fixes the underlying issue and
  re-runs the command.

---

## Failure modes and safeguards

| Failure | Behavior |
|---------|----------|
| Not in a git repo | Script exits 3 with a one-line message. The command logs a warning and continues (other phases already ran; commit is the last step). |
| `.claude/architect.git` absent | Script exits 0 with "disabled" notice. No noise. |
| `ARCHITECT_GIT_AUTOCOMMIT=false` | Same as absent — exit 0, skip. |
| Command not in `ARCHITECT_GIT_COMMIT_ON` | Exit 0, skip. |
| Nothing staged (no paths written) | Exit 0, skip. Avoids empty commits. |
| Pre-commit hook fails | Exit 5 with verbatim hook output. User fixes and re-runs. |
| GPG signing fails (`user.signingkey` set but no key) | Exit 5. User disables signing or fixes config. |
| Branch create conflicts with uncommitted changes | Exit 5 with guidance to stash/commit first. |
| Path outside git root | Exit 4 (misconfiguration; command bug). Surface to user. |
| Detached HEAD | Exit 5 with "refuse to commit on detached HEAD; `git switch -c` first". |
| Submodule touched | The helper script commits in the git worktree it runs in. Submodule commits stay in the submodule; super-repo commits (if any) are the user's job. The harness does not recurse. |

---

## Log-shape after a full feature walkthrough

After `/map → /review (all) → /synthesize → /design → /plan → /implement 1..6`
for a feature `color-encoding`, the expected `git log --oneline` is:

```
feat(color-encoding): phase 6 — literal-forbidden assertions
feat(color-encoding): phase 5 — dirColor runtime resolution
feat(color-encoding): phase 4 — dot-count rework
feat(color-encoding): phase 3 — CSS rewrite + legend
feat(color-encoding): phase 2 — JS dispatch wiring
feat(color-encoding): phase 1 — config & metadata foundation
docs(color-encoding): plan (6 phases)
docs(color-encoding): architect verdict READY
docs(color-encoding): design bundle (8 ADRs, 2 OQs) — status APPROVED
docs(color-encoding): design README — scope & acceptance criteria
docs(color-encoding): synthesize reviews (7 P0, 3 divergent)
docs(color-encoding): review [bioinf, graphic, ml, stat]
docs(color-encoding): map codebase (69 locs across 5 files)
```

13 commits for a full run. Granular, bisectable, rollback-friendly.
Users who find this too noisy drop cheap commands from
`ARCHITECT_GIT_COMMIT_ON`:

```bash
ARCHITECT_GIT_COMMIT_ON=design,plan,implement,meta-design,meta-apply
```

yielding just 10 commits (design, plan, six phases, meta-design,
meta-apply) for the same workflow.
