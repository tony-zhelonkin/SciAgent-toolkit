---
parent: ./README.md
view: behavioral
version: 2
revised: 2026-04-29
---

# Behavior — git automation (v2)

Behavioral view: when commits fire, what each commit looks like, how
branches and rollback behave, what fails where, what gets logged.

---

## The universal rule

> **Commit after the human gate. Never before. Never bypass.**

Every commit is a bookmark for "the user just approved this artifact".
Rejection is handled at the working-tree layer (`git restore`, with the
v2 dirty-path safety check) — it never becomes history. Pre-commit
hooks are respected; the harness never uses `--no-verify`, `--amend`,
`--no-gpg-sign`, or `-f`.

---

## Per-command commit flow

Each block below describes: (a) the gate that must pass, (b) the subject
the renderer composes, (c) the paths staged, (d) what happens on
rejection. Rendered subjects are now produced by `scripts/_render.py`,
not by command-file heredocs.

### `/map {feat}`

- **Gate:** Phase 3 user confirmation (use / regenerate / append). If the
  user picked "use existing", no write happened; the helper logs
  `skipped | empty` and exits 0.
- **Rendered subject:** `docs({feat}): map codebase ({loc} locs across {files} files)`
- **Rendered body** (optional): `Open questions: {n_oq} (see map.md § Open Questions).`
- **Paths:** `docs/{feat}/map.md`.
- **Rejection:** not applicable — gate runs before write.

### `/review {feat} --as <set>`

- **Gate:** Phase 4 completion. All N reviewers return; the set lands
  atomically.
- **Rendered subject:** `docs({feat}): review [{sorted_reviewer_list}]`
- **Paths:** every `docs/{feat}/review/*.md` written this invocation.
- **Rejection:** user re-runs `/review`; prior commit stands and the new
  pass produces a follow-up commit.

### `/synthesize {feat}`

- **Gate:** Phase 3 surface of Headline + P0 list.
- **Rendered subject:** `docs({feat}): synthesize reviews ({n_convergent} P0, {n_divergent} divergent)`
- **Paths:** `docs/{feat}/synthesis.md`.
- **Rejection:** user re-runs with more reviewers; prior synthesis commit stands.

### `/design {feat}` — **single commit (v2 simplification)**

v1 specified two or three commits (scope, bundle, verdict). v2 collapses
to one commit fired at Gate 2 (full bundle approval, post-status-flip).
The architect verdict, if READY in the same invocation, is included in
this commit; if NEEDS_ITERATION/NEEDS_DISCUSSION, the command halts
before the commit step and no commit fires.

- **Gate:** Phase 4 Gate 2 approval AND (in the same invocation) Phase 5
  architect verdict = READY. If Gate 2 passes but architect returns
  NEEDS_ITERATION/NEEDS_DISCUSSION, halt without committing.
- **Rendered subject:** `docs({feat}): design ({n_adrs} ADRs, {n_oqs} OQs) — APPROVED`
- **Rendered body** (optional): `Deferred: {deferred_ids_csv}.`
- **Paths:** `docs/{feat}/design/*.md` (the full bundle, README + 01–04
  + review.md). One SHA, one approval.
- **Rejection at Gate 1 or Gate 2:** no commit; user edits/re-runs;
  `git restore docs/{feat}/design/*.md` if user wants the working tree
  reset (subject to ADR-016 dirty-path check).

See ADR-014 in `03-decisions.md` for why v2 collapses /design's commits.

### `/architect {feat}` (standalone gate)

- **Gate:** Phase 4 verdict parse.
- **Rendered subject:** `docs({feat}): architect verdict {READY|NEEDS_ITERATION|NEEDS_DISCUSSION}`
- **Paths:** `docs/{feat}/design/review.md`. If the user consented to
  the "flip DRAFT → APPROVED" prompt, also the design/*.md files whose
  frontmatter was flipped.
- **Rejection:** not applicable — architect only reports a verdict.

### `/plan {feat}`

- **Gate:** Phase 4 approval.
- **Rendered subject:** `docs({feat}): plan ({n_phases} phases)`
- **Rendered body** (optional): one line per phase title from each
  `phase-NN.md` frontmatter `focus:`.
- **Paths:** `docs/{feat}/plan/README.md`, `docs/{feat}/plan/phase-*.md`.

### `/implement {feat} {N}`

- **Gate:** Phase 4 verification checklist passes. Anything failed →
  command halts at Phase 4; commit phase unreachable.
- **Rendered subject:** `feat({feat}): phase {N} — {phase_title}`
  - `{phase_title}` ← `phase-N.md` frontmatter `focus:` field.
  - `feat(` if any non-`docs/` path was modified; otherwise `docs(`.
- **Rendered body:** verification checklist as `- [x]` bullets, sourced
  from the phase doc's body.
- **Paths:** every file in `phase-N.md` frontmatter `files_touched:` +
  the phase doc itself (because `/implement` Phase 3 step 5 stamps the
  frontmatter).
- **Rejection:** command stops at Phase 4 on failure; no commit fires.

### `/verify {feat} [phases]`

- **Gate:** none — `/verify` is a read-only drift check that writes
  `verify.md`.
- **Rendered subject:** `docs({feat}): verify — verdict {CLEAN|INCOMPLETE|DRIFT|NEEDS_REVIEW}`
- **Paths:** `docs/{feat}/verify.md`.
- **Note:** users who run `/verify` frequently as a sanity check drop
  it from `ARCHITECT_GIT_COMMIT_ON`.

### `/status [slugs]`

- **No commit.** Chat-only output, no file writes. The L4 trigger
  carries one line: `This command does not commit.`

### `/meta-map [slugs]`

- **Gate:** Phase 4 user surfacing.
- **Rendered subject:** `docs(_meta): map portfolio ({n_features} features, {n_convergent} convergent concerns)`
- **Paths:** `docs/_meta/map.md`.

### `/meta-design [slugs]`

- **Gate:** Phase 4 Gate 1 approval (MADR bundle).
- **Rendered subject:** `docs(_meta): design ({n_madrs} MADRs) — APPROVED`
- **Rendered body** (optional): per-MADR one-liner.
- **Paths:** `docs/_meta/design.md`, `docs/_meta/deferred.md` if Phase
  6.5 appended rows.

### `/meta-apply [slugs]` — **single atomic commit (v2 simplification)**

v1 specified one commit per approved feature plus an architect-batch
commit afterward. v2 collapses to **one atomic commit per `/meta-apply`
invocation**, scope `_meta`, listing all approved features in the body.

- **Gate:** Phase 5 compound gate. After option 1 (approve all) or
  option 2 (approve some, reject others), the helper fires once.
- **Rendered subject:** `docs(_meta): apply meta-design ({n_madrs} MADRs) → {n_features} features`
- **Rendered body:** one line per approved feature listing applied
  MADRs and deferred-row IDs:
  ```
  Applied to:
  - biological-workflow: MADR-001, MADR-002, MADR-004 (deferred D-093..D-099)
  - color-encoding:      MADR-001, MADR-004 (deferred D-100..D-102)
  Rejected: theme-bundles (Gate 5 option 2)
  ```
- **Paths:** every approved feature's `docs/{feat}/design/*.md` files
  the feature-reviser touched + the assigned slice of
  `docs/_meta/deferred.md`. **All staged together** before one
  `git commit`. Hook failure leaves the working tree intact (atomic).
- **Phase 6 architect batch:** the architect re-gate writes
  `docs/{feat}/design/review.md` files. These are NOT committed by
  `/meta-apply`. They land in the working tree; the user runs
  `/architect {feat}` per feature (or accepts them in place) on their
  cadence. (Justified in ADR-014.)
- **Rejected features:** `git restore docs/{feat}/design/*.md` and
  strip the feature's deferred rows. Subject to ADR-016 dirty-path
  check; if hand-edits would be clobbered, the helper exits 4 and the
  user stashes first.

### `/meta-plan [slugs]`

- **Gate:** Phase 4 approval.
- **Rendered subject:** `docs(_meta): plan ({n_portfolio_phases} portfolio phases)`
- **Paths:** `docs/_meta/plan.md`.

---

## Commit message conventions

Every subject (rendered by `_render.py`) follows:

```
<type>(<scope>): <imperative summary>
```

- **Type:** `docs` for artifact-only commits, `feat` for `/implement`
  commits touching code.
- **Scope:** the feature slug (e.g., `color-encoding`) or `_meta` for
  portfolio-level artifacts.
- **Summary:** short, imperative, lowercase except for proper nouns
  (`UMAP`, `MADR`, ADR numbers).

Body rules:

- Wrap at 72 columns.
- Bullets with `-`.
- No horizontal rules, no markdown headings, no emojis.
- No trailers (`Co-Authored-By:`, `Signed-off-by:` unless user's git
  config produces it, etc.).

### What must never appear in any commit message

The renderer's `guard()` filter rejects any text containing any of:

```
Co-Authored-By: Claude
Co-Authored-By: Gemini
Co-Authored-By: Codex
Co-Authored-By: GPT
Co-Authored-By: OpenAI
Co-Authored-By: Anthropic
Co-Authored-By: Google
Generated with Claude
Generated with Gemini
Generated with Codex
🤖
Claude Code
Anthropic
<noreply@anthropic.com>
<noreply@google.com>
<noreply@openai.com>
```

These patterns are encoded in `scripts/_render.py` as `FORBIDDEN_PATTERNS`
and unit-tested in `scripts/_render_test.py`. See ADR-001.

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
docs(color-encoding): design (8 ADRs, 2 OQs) — APPROVED

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
docs(_meta): design (4 MADRs) — APPROVED

MADR-001: shared padj-floor convention across features
MADR-002: opacity priority chain (entity-level then score-level)
MADR-003: portfolio-coupled gene-set parse for umap + biological-workflow
MADR-004: ENTITY_PROFILES schema shared between color-encoding and biological-workflow
```

```
docs(_meta): apply meta-design (3 MADRs) → 2 features

Applied to:
- biological-workflow: MADR-001, MADR-002, MADR-004 (deferred D-093..D-099)
- color-encoding:      MADR-001, MADR-004 (deferred D-100..D-102)
```

---

## Branch behavior

### Default (`ARCHITECT_GIT_BRANCHES=false`)

- Everything commits to the current branch.
- The user decides when to branch, push, or merge.

### Opt-in (`ARCHITECT_GIT_BRANCHES=true`)

- On the **first** commit for a feature slug, the helper creates
  `${ARCHITECT_GIT_BRANCH_PREFIX}{slug}` (default: `architect/{slug}`)
  if it does not exist, and `git switch` to it.
- **Branch creation/switch is announced visibly:**
  - Chat: `architect-commit: switched to branch architect/{slug} (was main)`
  - Log: `... | switched | architect/{slug} (was main)`
- Subsequent commits for the same slug stay on that branch.
- If uncommitted changes prevent the switch, exit 5 with guidance to
  stash first. The user resolves and re-runs.
- The harness never deletes, rebases, or merges a branch.

### Meta-layer branches

Meta commands always commit on the current branch. (v1's
`ARCHITECT_GIT_META_BRANCH=portfolio` mode and its cross-stream-refusal
logic were removed in v2 — see ADR-015.)

### Worktrees

**Not supported.** Same as v1. Users who want them run the harness
inside a manually-created worktree.

---

## Rollback behavior

### Before a commit (working tree only) — with v2 dirty-path safety

- `git restore` on the artifact paths reverts uncommitted edits.
- **v2 addition:** before running `git restore`, the helper runs
  `git status --porcelain` against the target paths. If any path has
  uncommitted hand-edits the user might lose, the helper exits 4 with:
  ```
  architect-commit: would clobber uncommitted edits in:
    docs/foo/design/01-architecture.md (modified)
    docs/foo/design/03-decisions.md (modified)
  Stash or commit first; then re-run.
  ```
  The user stashes/commits and re-runs the rollback.
- See ADR-016.

### After a commit (history exists)

- `git revert <sha>` is the recommended path.
- `git reset --hard <sha-before>` on unpushed branches — destructive,
  user's call.
- The harness does not rewrite history.

### Pre-commit hook rejection

- Exit code 5; log line: `... | hook-failed | <hook-name-if-known>`.
- Calling command surfaces the hook's verbatim output to the user.
- No retry, no `--no-verify`. User fixes and re-runs the command.

---

## Transparency log (L0) format

Every helper invocation appends one line to `.claude/architect.git.log`:

```
<ISO-8601 timestamp> | <command> | <feature> | <verdict> | <sha-or-dash> | <subject-or-reason>
```

Verdict vocabulary:

| Verdict | When | sha field |
|---------|------|-----------|
| `committed` | `git commit` succeeded | short SHA |
| `skipped` | Disabled, not in COMMIT_ON, nothing staged, `--dry-run` | `-` |
| `refused` | Bad input, not in worktree, path traversal, dirty paths | `-` |
| `hook-failed` | `git commit` exit non-zero | `-` |
| `switched` | Branch switch (precedes the commit's own log line) | `-` |
| `dry-run` | `--dry-run` flag set | `-` |

The user audits with:

```bash
tail -20 .claude/architect.git.log
grep ' | hook-failed ' .claude/architect.git.log    # historical hook failures
grep " | umap | committed " .claude/architect.git.log | wc -l   # commits to umap
```

The log is never read by the harness — only written. It is the user's
audit surface, not state input.

---

## Failure modes and safeguards

| Failure | Behavior | Log verdict |
|---------|----------|-------------|
| Not in a git repo | Exit 3, one-line message | `refused | not-in-worktree` |
| `.claude/architect.git` absent | Exit 0 with "disabled" notice | `skipped | disabled` |
| `ARCHITECT_GIT_AUTOCOMMIT=false` | Same as absent | `skipped | disabled` |
| Command not in `ARCHITECT_GIT_COMMIT_ON` | Exit 0, skip | `skipped | not-in-allowlist` |
| Nothing staged (no paths written, or paths empty after `git add`) | Exit 0, skip | `skipped | empty` |
| Renderer raises `ValueError` (forbidden string) | Exit 5 with renderer's stderr | `refused | guard` |
| Pre-commit hook fails | Exit 5 with verbatim hook output | `hook-failed | <hook>` |
| GPG signing fails | Exit 5 | `hook-failed | gpg` |
| Branch create conflicts with uncommitted changes | Exit 5 with guidance | `refused | dirty-tree` |
| Path outside git root | Exit 4 (command bug) | `refused | path-traversal` |
| Detached HEAD | Exit 5 with "refuse to commit on detached HEAD" | `refused | detached-head` |
| Dirty paths would be clobbered by rollback | Exit 4 with file list | `refused | dirty-rollback` |
| Submodule touched | Helper commits in the worktree it runs in; super-repo is the user's job | `committed` (in the submodule's log) |

---

## Log-shape after a full feature walkthrough

After `/map → /review (all) → /synthesize → /design → /plan → /implement 1..6`
for a feature `color-encoding` (with default `ARCHITECT_GIT_COMMIT_ON`),
the expected `git log --oneline` is:

```
feat(color-encoding): phase 6 — literal-forbidden assertions
feat(color-encoding): phase 5 — dirColor runtime resolution
feat(color-encoding): phase 4 — dot-count rework
feat(color-encoding): phase 3 — CSS rewrite + legend
feat(color-encoding): phase 2 — JS dispatch wiring
feat(color-encoding): phase 1 — config & metadata foundation
docs(color-encoding): plan (6 phases)
docs(color-encoding): design (8 ADRs, 2 OQs) — APPROVED
docs(color-encoding): synthesize reviews (7 P0, 3 divergent)
docs(color-encoding): review [bioinf, graphic, ml, stat]
docs(color-encoding): map codebase (69 locs across 5 files)
```

**11 commits** for a full run (v1's 13 minus the two `/design` extras).
Granular, bisectable, rollback-friendly. Users who find this too noisy
drop cheap commands from `ARCHITECT_GIT_COMMIT_ON`:

```bash
ARCHITECT_GIT_COMMIT_ON=design,plan,implement,meta-design,meta-apply
```

yielding 9 commits (design, plan, six phases, and any meta work) for the
same workflow.
