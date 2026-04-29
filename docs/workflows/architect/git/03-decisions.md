---
parent: ./README.md
view: decision
version: 2
revised: 2026-04-29
---

# Decisions — git automation (v2)

Decision view: the architectural choices behind this proposal and the
alternatives explicitly rejected. One ADR per choice. v1 ADRs preserved
verbatim where unchanged; ADR-006 and ADR-012 superseded with a trail;
new ADRs (ADR-014, ADR-015, ADR-016, ADR-017, ADR-018) introduce the
v2-specific decisions.

---

## ADR-001 — No LLM / provider attribution in any commit, ever

**Status:** ACCEPTED (hard constraint from the user). Unchanged in v2.

**Context.** The user works across Claude Code, Gemini CLI, and Codex
CLI. Each CLI's built-in commit prompts append vendor-specific
attribution trailers and "Generated with …" lines. Adopting any of
these would (1) leak which LLM produced which commit, (2) fragment
history with multiple `Co-Authored-By:` trailers, (3) bind a user-facing
convention to a single vendor, and (4) pollute `git shortlog` with
noreply addresses unrelated to the human doing the work.

**Decision.** No trailers, no LLM mentions, ever. Two enforcement
points (v2):

1. `guard()` filter in `scripts/_render.py` rejects forbidden patterns
   before the message reaches `git commit`.
2. `scripts/architect-commit.sh` composes messages with no trailers and
   ignores any `ARCHITECT_GIT_ATTRIBUTION` value other than `false`.

Forbidden patterns (encoded in `_render.py`'s `FORBIDDEN_PATTERNS`):

- `Co-Authored-By: Claude`, `Gemini`, `Codex`, `GPT`, `OpenAI`,
  `Anthropic`, `Google`, or any variant
- `Generated with Claude`, `Generated with Gemini`, `Generated with Codex`
- `🤖 Generated with [Claude Code]` and markdown-link variants
- `Anthropic`, `Google`, `OpenAI` (bare in trailers or body)
- `<noreply@anthropic.com>`, `<noreply@google.com>`, `<noreply@openai.com>`
- Emoji-only lines (`🤖`, `🧠`, `✨`)

**v1 vs v2.** v1 placed a "Do NOT append authorship" reminder in every
command file's commit section and relied on a CI grep for backstop. v2
removes the per-command reminders (token noise), keeps the script's
no-trailer composition, and adds the renderer's `guard()` as the
primary enforcement (covered by pytest). See ADR-017.

**Consequences.**

- Every `commands/*.md` file's commit trigger is ≤3 lines; no
  attribution boilerplate.
- `pytest scripts/_render_test.py` exercises the forbidden patterns on
  every commit-class subject and body. Drift caught locally, not in CI.
- If a pre-commit hook in the user's repo tries to add a trailer, that's
  the user's repo policy; the harness does not fight it.

**Alternatives rejected.**

- *Per-provider attribution.* Rejected: user works across providers.
- *Opt-in attribution.* Rejected: foot-gun that drifts between provider
  defaults over time.
- *Single neutral trailer like `Co-Authored-By: Architect Harness`.*
  Rejected: fake authorship is still fake.

---

## ADR-002 — Commit after the human gate, never before

**Status:** ACCEPTED. Unchanged in v2.

**Context.** Each command has human gates. Commits can land before or
after the gate. Before-gate commits pollute history when work is
rejected; after-gate commits leave the working tree as the rejection
surface (`git restore`).

**Decision.** Commit after the gate. Generalises `/meta-apply` Phase 5's
existing rollback contract to every command.

**Consequences.** Rejection never pollutes history. Partial work within
a command does not auto-commit; the user handles it in the working tree
and reruns. v2 adds dirty-path safety to the rollback step (ADR-016).

**Alternatives rejected.** Commit-before-gate-revert-on-reject (history
pollution, harder bisect). Shadow-branch-fast-forward-on-approve
(over-engineered for a solo workflow).

---

## ADR-003 — Helper script, not Claude Code hooks

**Status:** ACCEPTED. Unchanged in v2.

**Context.** Hooks fire on tool-use events, not on command completion. A
slash command is many tool calls; a hook cannot easily tell "command
done, gate passed" from "intermediate tool call finished". Hooks are
also Claude-Code-specific.

**Decision.** Explicit helper script invoked from each command file.
Cross-provider portable. Centralised logic. Trivial to debug
(`--dry-run --verbose`).

**Consequences.** Each command file gains ~3 lines (v2 reduction from
v1's ~30; see ADR-012's supersession). The script is portable to
Gemini CLI and Codex CLI without changes.

---

## ADR-004 — Default branch = current; feature branches opt-in

**Status:** ACCEPTED. Unchanged in v2.

**Context.** Single-branch (solo, linear) vs feature-branch
(multi-feature parallel) workflows are both valid; binding to either is
wrong for the other.

**Decision.** Default `ARCHITECT_GIT_BRANCHES=false` — current branch.
Multi-feature users flip to `true`; the harness creates
`architect/{slug}` on the first commit for each feature.

**v2 addition.** Branch switches are now visibly announced to the user
(chat line + L0 log entry). See ADR-018 (hidden-state mitigations).

**Consequences.** Single-branch users see no new branches unless they
opt in. Multi-feature users get automatic isolation with one flag.
Branch-creation surface is no longer silent.

---

## ADR-005 — No worktrees

**Status:** ACCEPTED. Unchanged in v2.

**Context.** Worktrees give each feature an isolated checkout.
Superficially appealing but adds checkout juggling and stale-ref
edge cases.

**Decision.** Skip. Parallel agent dispatch already writes to distinct
paths. Users who want worktrees create them manually.

---

## ADR-006 — One commit per command invocation (v1)

**Status:** SUPERSEDED by ADR-014 (2026-04-29).

**Original v1 decision.** One commit per command invocation, with
`/design` producing two (one per gate) and `/meta-apply` producing one
per approved feature plus a final architect-batch commit.

**Why superseded.** Architectural review identified the multi-commit
exceptions as over-engineering for the harness's actual use cases:

- `/design`'s 2-or-3 commits buy nothing in a solo workflow; bisect
  granularity is already at the command level.
- `/meta-apply`'s per-feature loop creates a hidden-state failure mode:
  a hook failure mid-loop leaves working tree half-committed.

See ADR-014 for the v2 single-commit decision.

---

## ADR-007 — Never `--no-verify`, never `--amend`, never `-f`

**Status:** ACCEPTED. Unchanged in v2.

**Context.** Pre-commit hooks, gpg signing, and existing history
protection are repo-level policy.

**Decision.** Never bypass. Hook failure → exit 5, surface verbatim.

---

## ADR-008 — Stage only paths the command wrote

**Status:** ACCEPTED. Unchanged in v2.

**Context.** `git add .` and `git add -A` are dangerous in a working
tree with unrelated edits.

**Decision.** Explicit `--paths` list, built by the calling command.
Never `git add .` or `-A`.

---

## ADR-009 — No autopush, no autopull, no autotag

**Status:** ACCEPTED. Unchanged in v2.

**Decision.** Never. Push has remote side-effects (CI, webhooks); pull
can reset state; tags are human-curated. The user invites those
consequences when ready.

---

## ADR-010 — Config lives in `.claude/architect.git` as key=value

**Status:** ACCEPTED. Unchanged in v2.

**Decision.** Shell `KEY=value` form. Sourceable by the bash helper.

**v2 modification.** The key set is reduced — `ARCHITECT_GIT_META_BRANCH`
removed (see ADR-015). Otherwise unchanged.

---

## ADR-011 — `/status` never commits; `/verify` does (by default)

**Status:** ACCEPTED. Unchanged in v2.

**Decision.** `/status` is chat-only. `/verify` commits by default;
users who run it as a frequent sanity check drop it from
`ARCHITECT_GIT_COMMIT_ON`.

---

## ADR-012 — Commit subject templates live in commands, not config (v1)

**Status:** SUPERSEDED by ADR-017 (2026-04-29).

**Original v1 decision.** Templates live in command files. The helper
script receives the fully-rendered subject as `--subject=…`.

**Why superseded.** Architectural review identified ~325 lines of
bash-in-markdown (template + heredoc body composition) spread across
13 command files as a major source of token bloat AND a
single-source-of-truth violation for the `guard()` policy.

See ADR-017 for the v2 renderer-layer decision.

---

## ADR-013 — `activate-role.sh` offers optional git opt-in

**Status:** ACCEPTED. Unchanged in v2.

**Decision.** `activate-role.sh architect …` gains `--with-git` flag.
When present, copies `templates/architect.git.template →
.claude/architect.git` (never overwriting). When absent, nothing
changes.

---

## ADR-014 — Single commit per `/design` and per `/meta-apply` (v2)

**Status:** ACCEPTED (supersedes ADR-006's multi-commit exceptions).

**Context.** v1 prescribed 2–3 commits for `/design` (scope, bundle,
verdict) and N+1 commits for `/meta-apply` (one per approved feature
plus a final architect-batch). Architectural review identified two
problems:

1. **Token overhead per command file.** Multi-commit ceremonies require
   command files to hold multiple subject templates and conditional
   prose for each branch.
2. **Hidden-state risk in `/meta-apply`.** A pre-commit hook failure on
   the third per-feature commit leaves history with features 1–2
   committed and feature 3 in the working tree. The user inherits a
   half-committed mess that the helper cannot atomically reverse.

**Decision.** One commit per command invocation, full stop:

- `/design` fires one commit at Gate 2 (full bundle approval). If the
  architect verdict in the same invocation is READY, it's included; if
  NEEDS_ITERATION/NEEDS_DISCUSSION, the command halts before the commit
  step and no commit fires.
- `/meta-apply` fires one atomic commit covering all approved features.
  Phase 6 architect re-gate writes `review.md` files into the working
  tree but does NOT auto-commit them; the user runs `/architect {feat}`
  per feature on their cadence.

**Consequences.**

- Fewer commits per workflow (11 instead of 13 for a full feature run;
  meta-apply drops from N+1 to 1).
- Bisect granularity preserved at the command level — every
  human-approved unit of work maps to one SHA.
- Hook failure leaves working tree intact; no half-committed state.
- v2 log shape (see `02-behavior.md § Log-shape`) is the new contract.

**Alternatives rejected.**

- *Keep v1's multi-commit scheme.* Hidden-state failure mode unresolved.
- *Per-MADR commits inside `/meta-apply`.* Even worse granularity
  fragmentation; user couldn't bisect the apply as a unit.
- *Pre-`/meta-apply` stash + atomic apply.* Stash recovery complexity
  exceeds the value; ADR-014 is simpler.

---

## ADR-015 — Meta-layer always commits on current branch (v2)

**Status:** ACCEPTED.

**Context.** v1 introduced `ARCHITECT_GIT_META_BRANCH=portfolio` mode,
which commits meta-layer artifacts on a shared `architect/_meta` branch
when `ARCHITECT_GIT_BRANCHES=true`. The mode required a cross-stream
refusal: the helper would refuse to commit a feature's
`docs/{feat}/design/*.md` onto `architect/_meta`, exiting 5 with an
error directing the user to switch to the feature branch first.

**Decision.** Remove the mode. Meta-layer commands (`/meta-map`,
`/meta-design`, `/meta-apply`, `/meta-plan`) always commit on the
current branch.

**Rationale.**

- The cross-stream refusal logic adds bash complexity and a non-obvious
  error condition for a workflow nobody runs by default.
- Solo and multi-feature workflows alike commit meta work on whichever
  branch the user is currently on; this matches the "meta work rides
  the same branch as the features it touches" reality.
- If a user genuinely wants a portfolio branch, they `git switch -c
  architect/_meta` manually before invoking a meta command. The harness
  does not need to express this as a config knob.

**Consequences.**

- One fewer config key (`ARCHITECT_GIT_META_BRANCH`).
- One fewer error path in `architect-commit.sh`.
- `02-behavior.md § Branch behavior` no longer needs a "Meta-layer
  branches" subsection.

**Alternatives rejected.**

- *Keep the mode.* Edge case for ~zero users; complexity tax for everyone.
- *Make portfolio mode the default.* Breaks the "current branch"
  default for users who don't want `architect/_meta` to spring into
  existence.

---

## ADR-016 — Refuse `git restore` if dirty paths would be clobbered (v2)

**Status:** ACCEPTED.

**Context.** v1's `/meta-apply` rollback contract on rejection runs
`git restore docs/{slug}/design/*.md`. If the user has hand-edits to
those files mid-session (e.g., they were typing notes in a design doc
while reviewing), `git restore` silently discards those edits. This is
a hidden-state risk.

**Decision.** Before any `git restore` triggered by the harness, run
`git status --porcelain` against the target paths. If any path has
uncommitted hand-edits the user might lose, exit 4 with a verbatim
list:

```
architect-commit: would clobber uncommitted edits in:
  docs/foo/design/01-architecture.md (modified)
  docs/foo/design/03-decisions.md (modified)
Stash or commit first; then re-run.
```

The user stashes/commits their edits and re-runs the rollback step.

**Consequences.**

- One additional `git status --porcelain` check before any `git
  restore`. Negligible cost.
- Users who genuinely want to discard their edits stash with
  `git stash -u`, run the rollback, then drop the stash.

**Alternatives rejected.**

- *Silent restore (v1 behavior).* Hidden-state risk; data loss.
- *Auto-stash before restore, auto-pop on success.* Stash recovery
  complexity (what if the pop conflicts?). Defer until concrete
  motivation exists.

---

## ADR-017 — Subject/body composition lives in a Python renderer (v2)

**Status:** ACCEPTED (supersedes ADR-012).

**Context.** v1 placed commit-message templates inside each command
file's "Phase N: Commit" section as bash heredocs:

```bash
--body-file=<(printf 'phase-%02d: %s\n' ...)
```

This required:
- ~30 lines of boilerplate per command file × 13 commands = ~400 lines.
- Bash syntax embedded in markdown command specs.
- A separate "Do NOT append authorship" reminder in each file.
- A CI grep job to detect drift.

**Decision.** Pull all subject/body composition into
`scripts/_render.py`. Command files invoke
`architect-commit.sh --command=<X> --feature=<Y> --paths=<list>`; the
helper calls the renderer to obtain `(subject, body)`. Renderer reads
canonical frontmatter (post-Tier-1 phase-doc schema) for every field.

**Renderer module shape:**

```python
def render_subject_body(command, feature, paths, phase=None) -> (str, str):
    fn = COMMAND_HANDLERS[command]
    subject, body = fn(feature=feature, paths=paths, phase=phase)
    return guard(subject), guard(body)

# 13 _render_<command> functions, each ~10 lines.
# guard() filter rejects forbidden patterns.
```

**Tested with `scripts/_render_test.py` (pytest):**
- 13 commands × happy-path against fixture `docs/` tree.
- Forbidden-string rejection: 7 attempted-injection cases.
- Frontmatter-missing graceful failure.

**Consequences.**

- Command files shrink to ~3 lines for the commit trigger.
- All forbidden-string enforcement in one `guard()` function.
- CI grep replaced by local pytest run; faster feedback.
- Renderer testable without git.
- New command added to harness = new entry in `COMMAND_HANDLERS` +
  ~10 lines of renderer code + ~3 lines of trigger in command file.

**Net change:** ~400 lines of bash-in-markdown removed. ~150 lines of
Python (renderer) + ~80 lines of pytest (tests) added. Net **~170
lines** of code, replacing a sprawling, untestable, drift-prone
template surface with a single testable module.

**Alternatives rejected.**

- *Keep templates in command files (v1).* Token-bloated; not
  unit-testable; encourages drift.
- *Templates in config (key=value).* Mini-language to maintain;
  fragmentation between config and renderer logic.
- *Bash-only renderer in `architect-commit.sh`.* Hard to test; mixes
  policy with execution.

---

## ADR-018 — Append-only transparency log (v2)

**Status:** ACCEPTED.

**Context.** v1 surfaced helper outcomes only via chat output. Chat
scrolls away; the user has no durable record of "what did the harness
just do to git?" Branch switches under
`ARCHITECT_GIT_BRANCHES=true` were silent. This is the largest
hidden-state risk in v1.

**Decision.** Every helper invocation appends one line to
`.claude/architect.git.log`:

```
<ISO-8601 timestamp> | <command> | <feature> | <verdict> | <sha-or-dash> | <subject-or-reason>
```

Verdict vocabulary: `committed | skipped | refused | hook-failed |
switched | dry-run`. Append-only. No rotation, no parsing, no harness
read-back. The user audits with `tail`/`grep`.

**Consequences.**

- The user can answer "what did the harness commit yesterday?" with
  one shell command, even after weeks of session churn.
- Branch switches are no longer silent — they get a `switched | <new>
  (was <old>)` line in addition to the chat announcement.
- Hook failures, refusals, and skips are all traceable.
- File lives under `.claude/`, which the toolkit's `.gitignore` already
  excludes — log is local, not pushed.

**Alternatives rejected.**

- *Chat output only (v1).* Hidden-state risk; user cannot audit history.
- *Structured JSON log.* Heavier; harder to `tail`/`grep`. The pipe-
  delimited format is the right tradeoff.
- *Auto-rotate at N MB.* Premature; the file grows ~80 bytes per
  invocation; ~10 invocations/day = ~300 KB/year. Truncate manually
  if it ever matters.
- *Read the log back as harness state.* Rejected: the log is for the
  user, not for the harness. Reading it back creates the hidden state
  it's meant to expose.

---

## ADR-019 — Out of scope for v2

These remain plausible future features deliberately deferred:

- **`/architect-undo`** — convenience wrapper around `git revert HEAD`.
  `git revert` is already one command.
- **Pre-phase stash** — auto-stash before `/implement`, pop on rollback.
  Adds stash-recovery complexity. ADR-016's refuse-on-dirty path is
  simpler and equally safe.
- **Per-commit co-signing** — already handled by user's git config.
- **Pre-commit message hooks** (`prepare-commit-msg`) — out of scope.
- **Automatic `git tag`-ing of `/implement` final-phase commits.**
  Tagging is human-curated.
- **Automatic PR creation via `gh pr create`.** Vendor-locked.
- **Log rotation / dashboards on `.claude/architect.git.log`.** User's
  job; the log is plain text by design.
- **Provider-detection auto-disable.** ADR-001 covers all providers
  uniformly; no detection needed.

All of the above can be added post-v2 without changing the core
architecture if the motivation becomes concrete.
