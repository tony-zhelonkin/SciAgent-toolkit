---
parent: ./README.md
view: how-to
---

# Integration guide — wiring git automation into the harness

The concrete steps to take this proposal from spec to shipped. Read
`01-architecture.md` and `03-decisions.md` first; this file assumes
those decisions are accepted.

---

## Files to create

| Path | Purpose |
|------|---------|
| `scripts/architect-commit.sh` | The helper script. All git logic. |
| `templates/architect.git.template` | Config template. Copied into `.claude/architect.git` by `activate-role.sh --with-git`. |
| `docs/workflows/architect/git/examples/sample-log.md` | Worked `git log --oneline` example from a synthetic feature. |

## Files to modify

| Path | What changes |
|------|--------------|
| `commands/map.md` | + Phase N: Commit (see template below) |
| `commands/review.md` | + Phase N: Commit |
| `commands/synthesize.md` | + Phase N: Commit |
| `commands/design.md` | + Phase N: Commit A (Gate 1), Commit B (Gate 2), Commit C (architect READY) |
| `commands/architect.md` | + Phase N: Commit |
| `commands/plan.md` | + Phase N: Commit |
| `commands/implement.md` | + Phase 5 → Phase 6: Commit (shift the hand-off phase to Phase 6) |
| `commands/verify.md` | + Phase N: Commit |
| `commands/meta-map.md` | + Phase N: Commit |
| `commands/meta-design.md` | + Phase N: Commit |
| `commands/meta-apply.md` | + Phase 5 amendments (commit after option-1/option-2 resolution) + Phase 6 architect-batch commit |
| `commands/meta-plan.md` | + Phase N: Commit |
| `commands/status.md` | + "This command does not commit." note |
| `scripts/activate-role.sh` | + `--with-git` flag that copies the template |
| `docs/workflows/architect/README.md` | + link to `git/README.md` in document index |

---

## The helper script — `scripts/architect-commit.sh`

Specification (not yet code). Implement in bash, pass shellcheck,
executable bit set.

```bash
#!/usr/bin/env bash
# architect-commit.sh — provider-agnostic commit helper for the architect harness.
# Usage: see 01-architecture.md § 1.
set -euo pipefail

# --- Parse flags ---
# --command, --feature, --phase, --subject, --body, --body-file, --paths,
# --dry-run, --verbose. Unknown flag => exit 2.

# --- Load config (if present) ---
CONFIG=".claude/architect.git"
if [[ ! -f "$CONFIG" ]]; then
  echo "architect-commit: disabled (no $CONFIG); skipping"
  exit 0
fi
# shellcheck disable=SC1090
source "$CONFIG"
: "${ARCHITECT_GIT_AUTOCOMMIT:=false}"
: "${ARCHITECT_GIT_COMMIT_ON:=}"
: "${ARCHITECT_GIT_BRANCHES:=false}"
: "${ARCHITECT_GIT_BRANCH_PREFIX:=architect/}"
: "${ARCHITECT_GIT_META_BRANCH:=current}"

# Hardcoded safety rails — config values are ignored.
ARCHITECT_GIT_ATTRIBUTION=false
ARCHITECT_GIT_AUTOPUSH=false

# --- Early exits ---
[[ "$ARCHITECT_GIT_AUTOCOMMIT" == "true" ]] || { echo "architect-commit: AUTOCOMMIT=false; skipping"; exit 0; }
[[ -z "$ARCHITECT_GIT_COMMIT_ON" || ",${ARCHITECT_GIT_COMMIT_ON}," == *",${COMMAND},"* ]] \
  || { echo "architect-commit: command '$COMMAND' not in COMMIT_ON; skipping"; exit 0; }

git rev-parse --is-inside-work-tree >/dev/null 2>&1 \
  || { echo "architect-commit: not inside a git worktree" >&2; exit 3; }

# --- Branch management ---
if [[ "$ARCHITECT_GIT_BRANCHES" == "true" ]]; then
  # Meta commands may pin to a portfolio branch.
  if [[ "$FEATURE" == "_meta" && "$ARCHITECT_GIT_META_BRANCH" == "portfolio" ]]; then
    TARGET_BRANCH="${ARCHITECT_GIT_BRANCH_PREFIX}_meta"
  elif [[ "$FEATURE" != "_meta" ]]; then
    TARGET_BRANCH="${ARCHITECT_GIT_BRANCH_PREFIX}${FEATURE}"
  fi
  if [[ -n "${TARGET_BRANCH:-}" ]]; then
    CURRENT_BRANCH=$(git rev-parse --abbrev-ref HEAD)
    if [[ "$CURRENT_BRANCH" != "$TARGET_BRANCH" ]]; then
      if git show-ref --verify --quiet "refs/heads/$TARGET_BRANCH"; then
        git switch "$TARGET_BRANCH" \
          || { echo "architect-commit: failed to switch to $TARGET_BRANCH" >&2; exit 5; }
      else
        git switch -c "$TARGET_BRANCH" \
          || { echo "architect-commit: failed to create $TARGET_BRANCH" >&2; exit 5; }
      fi
    fi
  fi
fi

# --- Stage paths (explicit list only) ---
# Validate each path is inside the worktree; abort on traversal.
REPO_ROOT=$(git rev-parse --show-toplevel)
for p in "${PATH_LIST[@]}"; do
  ABS=$(readlink -f -- "$p" 2>/dev/null || realpath -- "$p")
  [[ "$ABS" == "$REPO_ROOT"/* ]] \
    || { echo "architect-commit: path '$p' outside repo root" >&2; exit 4; }
done
git add -- "${PATH_LIST[@]}"

# --- Skip empty commits ---
if git diff --cached --quiet; then
  echo "architect-commit: nothing staged; skipping"
  exit 0
fi

# --- Compose commit message ---
# Subject + optional body. NO TRAILERS. NO LLM MENTIONS.
MSG_FILE=$(mktemp)
trap 'rm -f "$MSG_FILE"' EXIT
{
  echo "$SUBJECT"
  if [[ -n "${BODY:-}" ]]; then
    echo ""
    echo "$BODY"
  fi
} >"$MSG_FILE"

# --- Commit (hooks respected; never --no-verify, --amend, or -f) ---
if [[ "$DRY_RUN" == "true" ]]; then
  echo "architect-commit: [dry-run] would commit with:"
  cat "$MSG_FILE"
  exit 0
fi

git commit -F "$MSG_FILE" \
  || { echo "architect-commit: git commit failed (see output above)" >&2; exit 5; }

echo "architect-commit: committed $(git rev-parse --short HEAD) on $(git rev-parse --abbrev-ref HEAD)"
```

### Self-test mode

```bash
./scripts/architect-commit.sh --selftest
```

runs:

1. A `--dry-run` invocation with fake args, confirms no git state change.
2. A config-absent path, confirms graceful skip.
3. A forbidden-string audit: `grep -RnE "(Co-Authored-By|Generated with|🤖|Anthropic|Claude Code)"` against the script and template, expecting **zero hits**.

Failure of any sub-test exits non-zero and prints which sub-test broke.

---

## The config template — `templates/architect.git.template`

Exact contents (copy verbatim into `.claude/architect.git` by
`activate-role.sh --with-git`):

```bash
# Architect git automation — local config
# Delete this file to disable auto-commit entirely.

ARCHITECT_GIT_AUTOCOMMIT=true
ARCHITECT_GIT_COMMIT_ON=map,review,synthesize,design,architect,plan,implement,verify,meta-map,meta-design,meta-apply,meta-plan
ARCHITECT_GIT_BRANCHES=false
ARCHITECT_GIT_BRANCH_PREFIX=architect/
ARCHITECT_GIT_META_BRANCH=current

# HARDCODED false — the script ignores any other value.
ARCHITECT_GIT_ATTRIBUTION=false
ARCHITECT_GIT_AUTOPUSH=false
```

Note the comments. Users will edit this file; comments are the nearest
docs when they're three weeks deep in a project.

---

## Command integration — the template

Paste this block near the end of each `commands/*.md` that writes an
artifact. Replace the three placeholders per command.

```markdown
## Phase N: Commit (if enabled)

After the gate approves the artifact(s) above, invoke the commit helper:

```
scripts/architect-commit.sh \
  --command=<command-name> \
  --feature={slug} \
  --subject="<type>({slug}): <templated-subject>" \
  --body-file=<(cat <<BODY
<optional multi-line body>
BODY
) \
  --paths=<comma-separated path list>
```

The helper no-ops if `.claude/architect.git` is absent or if this
command is not in `ARCHITECT_GIT_COMMIT_ON`. It commits on the current
branch unless `ARCHITECT_GIT_BRANCHES=true`.

**Do NOT append any authorship trailers, `Co-Authored-By` lines,
`Signed-off-by: <provider>` markers, or "Generated with …" lines. The
commit author is whatever git config produces; the harness is
provider-agnostic.**

If the helper exits non-zero (hook failure, branch conflict, etc.),
surface the error verbatim to the user and stop. Do NOT retry; do NOT
pass `--no-verify`.
```

---

## Per-command integration — exact subjects & paths

### `commands/map.md` — new Phase 4

```
--command=map
--feature={slug}
--subject="docs({slug}): map codebase ({loc} locs across {files} files)"
--body-file=<(printf 'Open questions: %d (see map.md § Open Questions).\n' {n_oq})
--paths=docs/{slug}/map.md
```

Skip the commit when Phase 1 resolved to "use existing" (nothing was
written).

### `commands/review.md` — new Phase 5

```
--command=review
--feature={slug}
--subject="docs({slug}): review [{sorted_reviewer_list}]"
--paths=docs/{slug}/review/<reviewer-1>.md,docs/{slug}/review/<reviewer-2>.md,...
```

Path list = files actually written this invocation (not pre-existing).

### `commands/synthesize.md` — new Phase 4

```
--command=synthesize
--feature={slug}
--subject="docs({slug}): synthesize reviews ({n_convergent} P0, {n_divergent} divergent)"
--paths=docs/{slug}/synthesis.md
```

Skip when Phase 1 early-exited (0 or 1 reviewer files).

### `commands/design.md` — three commit points

**After Phase 3a Gate 1:**

```
--command=design
--feature={slug}
--subject="docs({slug}): design README — scope & acceptance criteria"
--paths=docs/{slug}/design/README.md
```

**After Phase 4 Gate 2 (post status-flip to APPROVED):**

```
--command=design
--feature={slug}
--subject="docs({slug}): design bundle ({n_adrs} ADRs, {n_oqs} OQs) — status APPROVED"
--body-file=<(printf 'Deferred: %s.\n' "$deferred_ids_csv")
--paths=docs/{slug}/design/README.md,docs/{slug}/design/01-architecture.md,docs/{slug}/design/02-behavior.md,docs/{slug}/design/03-decisions.md[,...]
```

**After Phase 5 architect verdict (only if READY):**

```
--command=design
--feature={slug}
--subject="docs({slug}): architect verdict READY"
--paths=docs/{slug}/design/review.md
```

If verdict is NEEDS_ITERATION or NEEDS_DISCUSSION: no commit. The
command already halts and asks the user.

### `commands/architect.md` — new Phase 5

```
--command=architect
--feature={slug}
--subject="docs({slug}): architect verdict {VERDICT}"
--paths=docs/{slug}/design/review.md[,docs/{slug}/design/*.md if status flipped]
```

### `commands/plan.md` — new Phase 5

```
--command=plan
--feature={slug}
--subject="docs({slug}): plan ({n_phases} phases)"
--body-file=<(printf 'phase-%02d: %s\n' ...)
--paths=docs/{slug}/plan/README.md,docs/{slug}/plan/phase-01.md,...,docs/{slug}/plan/phase-{N}.md
```

### `commands/implement.md` — renumbered Phase 5 → Phase 6 (new Phase 5: Commit)

Insert BEFORE the current Phase 5 (which becomes Phase 6: Hand off to
next phase):

```
## Phase 5: Commit (if enabled)

Run ONLY if Phase 4 verification passed cleanly (no ❌ or ⚠️ marks).
If verification failed, the command halts at Phase 4 — this phase is
unreachable until the user fixes and re-runs.

Commit template:

--command=implement
--feature={slug}
--phase={N}
--subject="feat({slug}): phase {N} — {phase_title}"        # or docs({slug}): if no non-docs paths
--body-file=<( ... verification checklist as - [x] bullets ... )
--paths=<Files-to-Create from phase-N.md>,<Files-to-Modify from phase-N.md>,docs/{slug}/plan/phase-{N}.md
```

Only use `feat(` when `--paths` contains at least one path not starting
with `docs/`. Otherwise `docs(`.

Renumber the existing Phase 5 (hand-off) to Phase 6 and leave its body
unchanged. The WAIT-for-user instruction at the end still applies —
commit fires *before* the WAIT, so the committed state is what the user
sees when they decide to proceed.

### `commands/verify.md` — new Phase 10 (current Phase 9 stays as "Write verify.md")

```
--command=verify
--feature={slug}
--subject="docs({slug}): verify — verdict {CLEAN|INCOMPLETE|DRIFT|NEEDS_REVIEW}"
--paths=docs/{slug}/verify.md
```

### `commands/status.md` — add a one-line note

Replace any commit-adjacent content with:

```
## Commit behaviour

`/status` is chat-only and writes no files; this command never commits.
```

### `commands/meta-map.md` — new Phase 5

```
--command=meta-map
--feature=_meta
--subject="docs(_meta): map portfolio ({n_features} features, {n_convergent} convergent concerns)"
--paths=docs/_meta/map.md
```

### `commands/meta-design.md` — new Phase 7

```
--command=meta-design
--feature=_meta
--subject="docs(_meta): design ({n_madrs} MADRs) — status APPROVED"
--body-file=<( per-MADR one-liner listing )
--paths=docs/_meta/design.md[,docs/_meta/deferred.md if appended]
```

### `commands/meta-apply.md` — amend Phase 5 and add Phase 6.5 (architect-batch commit)

Amend Phase 5 option 1 ("approve all OK features"):

> After flipping `status: DRAFT → APPROVED` on each approved feature,
> invoke the commit helper ONCE PER approved feature (serially), before
> dispatching the architect batch in Phase 6:
>
> ```
> --command=meta-apply
> --feature={feat}
> --subject="docs({feat}): apply meta-design ({n_madrs} MADRs)"
> --body-file=<( applied-MADR list + deferred-row IDs )
> --paths=docs/{feat}/design/*.md,<slice of docs/_meta/deferred.md for this feature's ID block>
> ```

Amend Phase 5 option 2 ("approve some, reject others"):

> For REJECTED features: run `git restore docs/{feat}/design/*.md` and
> strip their deferred rows. **Do not commit rejected features.** For
> APPROVED features: proceed as option 1.

Add Phase 6.5 after architect batch returns (between Phase 6 and
Phase 7):

```
--command=meta-apply
--feature=_meta
--subject="docs(_meta): architect gate post-meta-apply — {R} READY, {I} NEEDS_ITERATION, {D} NEEDS_DISCUSSION"
--paths=docs/{feat-1}/design/review.md,docs/{feat-2}/design/review.md,...[,any design/*.md files where APPROVED was rolled back to DRAFT in Phase 7]
```

Phase 7's "roll back APPROVED → DRAFT for NEEDS_ITERATION features"
edits are included in this Phase 6.5 commit.

### `commands/meta-plan.md` — new Phase 7

```
--command=meta-plan
--feature=_meta
--subject="docs(_meta): plan ({n_portfolio_phases} portfolio phases)"
--paths=docs/_meta/plan.md
```

---

## `activate-role.sh` change

Add one flag and one block:

```bash
# New flag: --with-git
# ...in the arg parsing loop:
    --with-git) WITH_GIT=true ;;
# ...

# Near the end, after agent/skill symlinking:
if [[ "${WITH_GIT:-false}" == "true" ]]; then
  TEMPLATE="$TOOLKIT_ROOT/templates/architect.git.template"
  TARGET="$PROJECT_DIR/.claude/architect.git"
  if [[ -f "$TARGET" ]]; then
    echo "architect-commit: $TARGET already exists; leaving it alone"
  else
    cp "$TEMPLATE" "$TARGET"
    echo "architect-commit: installed $TARGET (edit to customise)"
  fi
fi
```

---

## Rollout order

Do not do all of this in one PR. The safe order:

1. **PR 1 — helper script + template (no command changes).**
   Adds `scripts/architect-commit.sh` and
   `templates/architect.git.template` with `--selftest` passing.
   Run `./scripts/architect-commit.sh --selftest` in CI (if CI exists).
   Run the forbidden-string audit:
   ```
   grep -RnE "(Co-Authored-By|Generated with|🤖|Anthropic|Claude Code|Gemini|Codex)" \
     scripts/architect-commit.sh templates/architect.git.template
   ```
   Expect zero hits.
2. **PR 2 — activate-role.sh `--with-git` flag.** Copies the template
   into `.claude/architect.git` on opt-in. No command changes yet; the
   file exists but no command invokes the script.
3. **PR 3 — wire `/map` and `/verify` only.** Small, low-blast-radius
   commands. Run a full walkthrough on a test feature; confirm
   `git log` shape matches `02-behavior.md § Log-shape`.
4. **PR 4 — wire `/review`, `/synthesize`, `/plan`.** Artifact-only
   commands, no code touches.
5. **PR 5 — wire `/design` and `/architect`.** Two-commit design is the
   subtle case; test both gates fire correctly.
6. **PR 6 — wire `/implement`.** Code-touching command. Test with a
   feature where pre-commit hooks fail intentionally to confirm the
   halt-and-surface flow.
7. **PR 7 — wire meta commands (`/meta-map`, `/meta-design`,
   `/meta-apply`, `/meta-plan`).** `/meta-apply`'s per-feature commit
   loop is the largest change.
8. **PR 8 — docs update.** Add `docs/workflows/architect/git/` to the
   main README index, update `05-meta-approach.md` with cross-refs to
   the meta-layer commit behaviour.

Between PRs, the harness is usable: users without `.claude/architect.git`
see no behavioural change at all; users who opted in see progressively
more commands auto-committing.

---

## Acceptance checklist

A v1 rollout is done when:

- [ ] `scripts/architect-commit.sh --selftest` passes.
- [ ] Forbidden-string grep across the script and all touched
      `commands/*.md` returns zero hits.
- [ ] Walkthrough on a throwaway feature produces the log shape in
      `02-behavior.md § Log-shape`.
- [ ] With `.claude/architect.git` absent, the harness behaves exactly
      as it did before (no new commits, no script invocations that
      fail).
- [ ] With `ARCHITECT_GIT_BRANCHES=true`, a second feature creates a
      distinct `architect/{slug}` branch without clobbering the first.
- [ ] Intentionally failing a pre-commit hook halts the harness with
      the hook output visible to the user; no `--no-verify` is ever
      invoked.
- [ ] `/meta-apply` with a mix of approved/rejected features produces
      exactly one commit per approved feature and no commit (plus
      `git restore`) per rejected feature.

When all boxes are checked, update this doc's status in the parent
README from "proposal" to "implemented".

---

## After-implementation maintenance

- When a new command is added to the harness (e.g. a hypothetical
  `/deprecate` or `/release`), add it to `ARCHITECT_GIT_COMMIT_ON`'s
  default list and add a "Phase N: Commit" block to the command file.
- When a new reviewer is added (e.g. `security`), no commit changes are
  needed — `/review` already globs the review directory.
- If pre-commit hooks become a common source of friction, consider
  adding a `--skip-hooks=<specific-hook-name>` passthrough. Still NEVER
  a blanket `--no-verify`.
- Keep the forbidden-string grep in CI. It is the single strongest
  guard against attribution drift as Claude Code / Gemini / Codex
  update their defaults.
