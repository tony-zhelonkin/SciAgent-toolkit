---
parent: ./README.md
view: how-to
version: 2
revised: 2026-04-29
---

# Integration guide — wiring git automation into the harness (v2)

The concrete steps to take this proposal from spec to shipped. The adjacent
`01-architecture.md` and `03-decisions.md` preserve the v2 design rationale;
this guide's file map and project-opt-in section record the current
router/reference integration points after the three-verb demolition.

---

## Files to create

| Path | Purpose |
|------|---------|
| `scripts/architect-commit.sh` | L2 git operator. ~120 lines bash. |
| `scripts/_render.py` | L3 subject/body renderer + `guard()` filter. ~150 lines Python. |
| `scripts/_render_test.py` | pytest: 13 workflow operations × happy-path + forbidden-string rejection. ~80 lines. |
| `templates/architect.git.template` | L1 config template. ~15 lines. |
| `docs/workflows/architect/git/examples/sample-log.md` | Worked `git log --oneline` from one feature run. |
| `docs/workflows/architect/git/examples/sample-architect.git.log` | Worked transparency log from the same run. |

## Files to modify

| Path | What changes |
|------|--------------|
| `skills/architecture-first-dev/references/map.md` | + 3-line "Phase N: Commit" trigger |
| `skills/architecture-first-dev/references/review.md` | + 3-line trigger |
| `skills/architecture-first-dev/references/synthesize.md` | + 3-line trigger |
| `skills/architecture-first-dev/references/design.md` | + 3-line trigger (single commit at Gate 2 — see ADR-014) |
| `skills/architecture-first-dev/references/architect.md` | + 3-line trigger |
| `skills/architecture-first-dev/references/plan.md` | + 3-line trigger |
| `commands/implement.md` | + 3-line trigger as Phase 5 (current Phase 5 hand-off becomes Phase 6) |
| `skills/architecture-first-dev/references/verify.md` | + 3-line trigger |
| `skills/architecture-first-dev/references/portfolio.md` | + triggers for `meta-map`, `meta-design`, `meta-apply`, and `meta-plan` |
| `skills/architecture-first-dev/references/status.md` | + one-line note: "This route does not commit." |
| `docs/workflows/architect/git/README.md` | + explicit project opt-in instructions |
| `docs/workflows/architect/README.md` | + link to `git/README.md` in document index |

Total docs/script delta from v1: **~325 lines fewer in stage specifications**
(no bash heredocs); **~230 lines added in scripts** (renderer + tests +
helper trim); net ≈ 95 lines smaller AND testable AND single-source-of-truth.

---

## L2 — `scripts/architect-commit.sh`

Specification (not yet code). Implement in bash, pass `shellcheck`,
executable bit set.

```bash
#!/usr/bin/env bash
# architect-commit.sh — provider-agnostic commit helper for the architect harness.
# Usage: see git/01-architecture.md § L2.
set -euo pipefail

# --- Constants ---
LOG_FILE=".claude/architect.git.log"
CONFIG=".claude/architect.git"
RENDERER="scripts/_render.py"

# --- Parse flags ---
# --command, --feature, --paths, --phase, --dry-run, --verbose, --selftest
# Unknown flag => exit 2.

# --- log_line: append one structured line to L0 transparency log ---
log_line() {
  local verdict="$1" sha="${2:-}" detail="${3:-}"
  printf '%s | %-12s | %-7s | %-12s | %-7s | %s\n' \
    "$(date -u +%Y-%m-%dT%H:%M:%S+00:00)" \
    "$COMMAND" "$FEATURE" "$verdict" "${sha:--}" "$detail" \
    >>"$LOG_FILE" 2>/dev/null || true
}

# --- Selftest mode ---
if [[ "${1:-}" == "--selftest" ]]; then
  # 1. Dry-run with fake args, confirm exit 0 and no git state change.
  # 2. Config-absent path, confirm graceful skip.
  # 3. Renderer reachability check (python3 + scripts/_render.py present).
  # Failure of any sub-test exits non-zero with which sub-test broke.
  exit 0
fi

# --- Load config ---
if [[ ! -f "$CONFIG" ]]; then
  log_line "skipped" "-" "disabled (no $CONFIG)"
  echo "architect-commit: disabled (no $CONFIG); skipping"
  exit 0
fi
# shellcheck disable=SC1090
source "$CONFIG"
: "${ARCHITECT_GIT_AUTOCOMMIT:=false}"
: "${ARCHITECT_GIT_COMMIT_ON:=}"
: "${ARCHITECT_GIT_BRANCHES:=false}"
: "${ARCHITECT_GIT_BRANCH_PREFIX:=architect/}"

# Hardcoded safety rails — config values are ignored.
ARCHITECT_GIT_ATTRIBUTION=false
ARCHITECT_GIT_AUTOPUSH=false

# --- Early exits ---
if [[ "$ARCHITECT_GIT_AUTOCOMMIT" != "true" ]]; then
  log_line "skipped" "-" "AUTOCOMMIT=false"
  echo "architect-commit: AUTOCOMMIT=false; skipping"; exit 0
fi
if [[ -n "$ARCHITECT_GIT_COMMIT_ON" && ",${ARCHITECT_GIT_COMMIT_ON}," != *",${COMMAND},"* ]]; then
  log_line "skipped" "-" "command not in COMMIT_ON"
  exit 0
fi
git rev-parse --is-inside-work-tree >/dev/null 2>&1 \
  || { log_line "refused" "-" "not-in-worktree"; echo "architect-commit: not in a git worktree" >&2; exit 3; }

# --- Branch management (visible switch) ---
if [[ "$ARCHITECT_GIT_BRANCHES" == "true" && "$FEATURE" != "_meta" ]]; then
  TARGET="${ARCHITECT_GIT_BRANCH_PREFIX}${FEATURE}"
  CURRENT=$(git rev-parse --abbrev-ref HEAD)
  if [[ "$CURRENT" != "$TARGET" ]]; then
    if git show-ref --verify --quiet "refs/heads/$TARGET"; then
      git switch "$TARGET" || { log_line "refused" "-" "branch-switch-failed"; exit 5; }
    else
      git switch -c "$TARGET" || { log_line "refused" "-" "branch-create-failed"; exit 5; }
    fi
    log_line "switched" "-" "$TARGET (was $CURRENT)"
    echo "architect-commit: switched to branch $TARGET (was $CURRENT)"
  fi
fi

# --- Stage paths (explicit list only; reject traversal) ---
REPO_ROOT=$(git rev-parse --show-toplevel)
for p in "${PATH_LIST[@]}"; do
  ABS=$(readlink -f -- "$p" 2>/dev/null || realpath -- "$p")
  [[ "$ABS" == "$REPO_ROOT"/* ]] \
    || { log_line "refused" "-" "path-outside-repo: $p"; echo "architect-commit: path '$p' outside repo root" >&2; exit 4; }
done
git add -- "${PATH_LIST[@]}"

# --- Skip empty commits ---
if git diff --cached --quiet; then
  log_line "skipped" "-" "empty staging"
  echo "architect-commit: nothing staged; skipping"
  exit 0
fi

# --- Render subject + body via L3 ---
RENDER_OUT=$(python3 "$RENDERER" \
  --command="$COMMAND" --feature="$FEATURE" \
  ${PHASE:+--phase=$PHASE} \
  --paths="$(IFS=,; echo "${PATH_LIST[*]}")" 2>&1) || {
  log_line "refused" "-" "renderer-error"
  echo "architect-commit: renderer error:" >&2
  echo "$RENDER_OUT" >&2
  exit 5
}

# RENDER_OUT format: subject\n\nbody (body may be empty).
SUBJECT=$(printf '%s\n' "$RENDER_OUT" | sed -n '1p')
BODY=$(printf '%s\n' "$RENDER_OUT" | sed -n '3,$p')

# --- Compose commit message file (NO TRAILERS) ---
MSG_FILE=$(mktemp)
trap 'rm -f "$MSG_FILE"' EXIT
{
  echo "$SUBJECT"
  if [[ -n "$BODY" ]]; then
    echo ""
    echo "$BODY"
  fi
} >"$MSG_FILE"

# --- Dry-run short circuit ---
if [[ "$DRY_RUN" == "true" ]]; then
  log_line "dry-run" "-" "$SUBJECT"
  echo "architect-commit: [dry-run] would commit:"
  cat "$MSG_FILE"
  exit 0
fi

# --- Commit (hooks respected; never --no-verify, --amend, or -f) ---
if ! git commit -F "$MSG_FILE"; then
  log_line "hook-failed" "-" "git-commit-rc=$?"
  echo "architect-commit: git commit failed (see output above)" >&2
  exit 5
fi

SHA=$(git rev-parse --short HEAD)
BRANCH=$(git rev-parse --abbrev-ref HEAD)
log_line "committed" "$SHA" "$SUBJECT"
echo "architect-commit: ✓ $SUBJECT ($SHA) on $BRANCH"
```

---

## L3 — `scripts/_render.py`

Specification (not yet code).

```python
#!/usr/bin/env python3
"""Subject/body renderer for architect-commit.sh.

Reads canonical frontmatter from docs/{feature}/... and produces a
deterministic (subject, body) tuple. Applies guard() to reject any
LLM/vendor attribution string before returning.
"""
import argparse
import re
import sys
from pathlib import Path

FORBIDDEN_PATTERNS = [
    re.compile(p, re.IGNORECASE) for p in [
        r"co-authored-by:\s*(claude|gemini|codex|gpt|openai|anthropic|google)",
        r"generated with (claude|gemini|codex|openai)",
        r"claude code", r"anthropic", r"openai", r"<noreply@(anthropic|google|openai)\.com>",
        r"🤖", r"🧠", r"✨",
    ]
]


class GuardError(ValueError):
    pass


def guard(text: str) -> str:
    for pat in FORBIDDEN_PATTERNS:
        if pat.search(text):
            raise GuardError(f"forbidden pattern matched: {pat.pattern!r}")
    return text


def _read_frontmatter(path: Path) -> dict:
    """Parse simple YAML frontmatter from the top of a markdown file."""
    text = path.read_text()
    m = re.match(r"^---\s*\n(.*?\n)---\s*\n", text, re.DOTALL)
    if not m:
        return {}
    fm: dict = {}
    for line in m.group(1).splitlines():
        if ":" in line and not line.startswith("  "):
            k, _, v = line.partition(":")
            fm[k.strip()] = v.strip()
    return fm


def _render_map(feature, **_):
    fm = _read_frontmatter(Path(f"docs/{feature}/map.md"))
    loc = fm.get("loc", "?")
    files = fm.get("files", "?")
    n_oq = fm.get("open_questions", "0")
    subject = f"docs({feature}): map codebase ({loc} locs across {files} files)"
    body = f"Open questions: {n_oq} (see map.md § Open Questions)." if n_oq != "0" else ""
    return subject, body


def _render_review(feature, paths, **_):
    reviewers = sorted({Path(p).stem for p in paths if "/review/" in p})
    subject = f"docs({feature}): review [{', '.join(reviewers)}]"
    return subject, ""


def _render_implement(feature, phase, **_):
    fm = _read_frontmatter(Path(f"docs/{feature}/plan/phase-{phase:02d}.md"))
    title = fm.get("focus", f"phase {phase}")
    subject = f"feat({feature}): phase {phase} — {title}"
    body = ""  # body composed from `## Verification` checkboxes if present
    return subject, body


# ... 10 more _render_<command> functions
COMMAND_HANDLERS = {
    "map": _render_map,
    "review": _render_review,
    "synthesize": lambda **kw: ...,
    "design": lambda **kw: ...,
    "architect": lambda **kw: ...,
    "plan": lambda **kw: ...,
    "implement": _render_implement,
    "verify": lambda **kw: ...,
    "meta-map": lambda **kw: ...,
    "meta-design": lambda **kw: ...,
    "meta-apply": lambda **kw: ...,
    "meta-plan": lambda **kw: ...,
}


def render_subject_body(command, feature, paths, phase=None):
    fn = COMMAND_HANDLERS[command]
    subject, body = fn(feature=feature, paths=paths, phase=phase)
    return guard(subject), guard(body)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--command", required=True)
    ap.add_argument("--feature", required=True)
    ap.add_argument("--paths", required=True, help="comma-separated")
    ap.add_argument("--phase", type=int, default=None)
    args = ap.parse_args()
    paths = [p for p in args.paths.split(",") if p]
    try:
        subject, body = render_subject_body(args.command, args.feature, paths, args.phase)
    except GuardError as e:
        print(f"_render: guard rejected: {e}", file=sys.stderr)
        sys.exit(2)
    except Exception as e:
        print(f"_render: error: {e}", file=sys.stderr)
        sys.exit(1)
    print(subject)
    print()
    if body:
        print(body)


if __name__ == "__main__":
    main()
```

---

## Pytest — `scripts/_render_test.py`

Coverage:

- 13 workflow operations × happy-path: feed a fixture `docs/` tree, assert
  `(subject, body)` shape matches the catalog in `02-behavior.md §
  Worked examples`.
- 7 forbidden-string injection cases (one per pattern class): assert
  `guard()` raises `GuardError`.
- Frontmatter-missing graceful failure: render returns a generic-but-
  valid subject when `phase-N.md` lacks fields.
- Idempotency: same inputs → same outputs.

Run with:

```bash
pytest scripts/_render_test.py -v
```

---

## L1 — `templates/architect.git.template`

Exact contents. A project opts in by copying this template to
`.claude/architect.git` explicitly:

```bash
# Architect git automation — local config
# Delete this file to disable auto-commit entirely.

ARCHITECT_GIT_AUTOCOMMIT=true

# Which workflow operations may auto-commit. Comma-separated. Default: all.
ARCHITECT_GIT_COMMIT_ON=map,review,synthesize,design,architect,plan,implement,verify,meta-map,meta-design,meta-apply,meta-plan

# Create a feature branch (prefix + slug) on first commit for that
# feature. Leave false for single-branch workflows.
ARCHITECT_GIT_BRANCHES=false
ARCHITECT_GIT_BRANCH_PREFIX=architect/

# HARDCODED false — the script ignores any other value.
ARCHITECT_GIT_ATTRIBUTION=false
ARCHITECT_GIT_AUTOPUSH=false
```

Note: `ARCHITECT_GIT_META_BRANCH` was removed in v2 (see ADR-015).

---

## L4 — Stage-level integration template

Add this block to each router reference or external command that writes an
artifact. The renderer's `guard()` enforces the attribution constraint.

```markdown
## Phase N: Commit (if enabled)

Invoke `scripts/architect-commit.sh --command=<this-command-name> \
--feature={slug} --paths=<comma-separated paths just written>`. The
helper composes the message via `scripts/_render.py` and no-ops if
`.claude/architect.git` is absent.

If the helper exits non-zero (hook failure, branch conflict, etc.),
surface the error verbatim to the user and stop. Do NOT retry; do NOT
pass `--no-verify`.
```

`/status` carries instead:

```markdown
## Commit behaviour

`/status` is chat-only and writes no files; this command never commits.
```

---

## Per-operation path lists

The `<paths>` substitution per operation (the renderer derives subject
and body from frontmatter, so paths are the only command-specific
context here):

| Command | `--paths` value |
|---------|-----------------|
| `/map` | `docs/{slug}/map.md` |
| `/review` | comma-separated list of review files written this invocation |
| `/synthesize` | `docs/{slug}/synthesis.md` |
| `/design` | `docs/{slug}/design/README.md,docs/{slug}/design/01-architecture.md,docs/{slug}/design/02-behavior.md,docs/{slug}/design/03-decisions.md,docs/{slug}/design/review.md` (full bundle, single commit) |
| `/architect` | `docs/{slug}/design/review.md` + any `design/*.md` whose status flipped |
| `/plan` | `docs/{slug}/plan/README.md` + every `docs/{slug}/plan/phase-*.md` |
| `/implement` | every path in `phase-N.md` frontmatter `files_touched:` + the phase doc itself; `--phase=N` flag also passed |
| `/verify` | `docs/{slug}/verify.md` |
| `/meta-map` | `docs/_meta/map.md`; `--feature=_meta` |
| `/meta-design` | `docs/_meta/design.md`[,`docs/_meta/deferred.md`]; `--feature=_meta` |
| `/meta-apply` | every approved feature's `docs/{feat}/design/*.md` + relevant `docs/_meta/deferred.md` slice (single atomic commit); `--feature=_meta` |
| `/meta-plan` | `docs/_meta/plan.md`; `--feature=_meta` |

---

## Project opt-in

Git automation remains project-owned configuration. After `scio link`, a
project enables this proposal explicitly:

```bash
mkdir -p .claude
cp /path/to/SciAgent-toolkit/templates/architect.git.template \
  .claude/architect.git
```

An existing `.claude/architect.git` stays project-owned. `scio link`
manages catalog links and guardrail hooks, so this opt-in file remains outside
its write set.

---

## Rollout — 5 PRs (compressed from v1's 8)

The reduced v2 scope merits a tighter sequence. Between PRs the harness
remains usable; users without `.claude/architect.git` see no
behavioural change.

1. **PR 1 — L0 + L1 + L2 + L3 + tests.**
   Adds `scripts/architect-commit.sh`, `scripts/_render.py`,
   `scripts/_render_test.py`, and `templates/architect.git.template`.
   Acceptance: `pytest scripts/_render_test.py -v` passes;
   `./scripts/architect-commit.sh --selftest` passes.
2. **PR 2 — document explicit project opt-in.** A project copies the template
   into `.claude/architect.git`; catalog binding stays unchanged.
3. **PR 3 — wire cheap routes.** `/map`, `/verify`, `/status`. Run a
   walkthrough on a test feature; confirm `git log` shape matches
   `02-behavior.md § Log-shape` and `.claude/architect.git.log` reflects
   each invocation.
4. **PR 4 — wire artifact and design routes.** `/review`,
   `/synthesize`, `/plan`, `/design`, `/architect`. Test single-commit
   `/design` (Gate 2 only).
5. **PR 5 — wire `/implement` and portfolio routes.** Code-touching
   command with intentional pre-commit hook failure to verify
   halt-and-surface. Then `/meta-map`, `/meta-design`, `/meta-apply`
   (single atomic commit), `/meta-plan`.

---

## Acceptance checklist

A v2 rollout is done when:

- [ ] `pytest scripts/_render_test.py -v` passes (13 workflow operations × happy
      path + 7 forbidden-string rejections + frontmatter-missing).
- [ ] `./scripts/architect-commit.sh --selftest` passes.
- [ ] Walkthrough on a throwaway feature produces the 11-commit log
      shape in `02-behavior.md § Log-shape`.
- [ ] `.claude/architect.git.log` after the walkthrough has one line
      per command invocation, no missing entries.
- [ ] With `.claude/architect.git` absent, the harness behaves exactly
      as it did before (no commits, helper exits 0 quietly).
- [ ] With `ARCHITECT_GIT_BRANCHES=true`, branch creation emits the
      visible chat line AND a `switched` log line.
- [ ] Intentionally failing a pre-commit hook halts the harness with
      the hook output; `hook-failed` line in the log.
- [ ] `/meta-apply` with mixed approved/rejected features produces
      exactly one atomic commit covering all approved features (not N).
- [ ] `git restore` against a path with uncommitted hand-edits is
      refused with the file list (ADR-016).
- [ ] Forbidden-string injection attempts in any rendered field raise
      `GuardError` and surface `refused | renderer-error` in the log.

When all boxes are checked, update `git/README.md` status from
"proposal v2" to "implemented".

---

## After-implementation maintenance

- **New route added to the workflow** (e.g. hypothetical `/deprecate`):
  add an entry to `COMMAND_HANDLERS` in `_render.py`, add a happy-path
  test in `_render_test.py`, add the command to
  `ARCHITECT_GIT_COMMIT_ON`'s default list, paste the L4 trigger block
  into the new reference or external command file.
- **New reviewer added to `/review`** (e.g. `security`): no changes
  needed. The renderer derives the reviewer list from the `--paths`
  list automatically.
- **Pre-commit hooks become friction:** never blanket `--no-verify`.
  If a specific hook is the issue, consider a future
  `--skip-hooks=<specific-hook-name>` passthrough; surface the
  motivation as a deferred decision under ADR-019.
- **Forbidden-string drift** (a provider invents a new attribution
  format): add the new pattern to `FORBIDDEN_PATTERNS` in
  `_render.py`; add a test case in `_render_test.py`. Single-file
  change.
- **Log file gets large:** `truncate -s 0 .claude/architect.git.log`
  manually. Don't add rotation logic to the harness.
