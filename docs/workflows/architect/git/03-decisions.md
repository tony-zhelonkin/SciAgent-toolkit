---
parent: ./README.md
view: decision
---

# Decisions — git automation

Decision view: the architectural choices behind this proposal and the
alternatives explicitly rejected. One ADR per choice.

---

## ADR-001 — No LLM / provider attribution in any commit, ever

**Status:** ACCEPTED (hard constraint from the user).

**Context.** The user works across Claude Code, Gemini CLI, and Codex
CLI. Claude Code's built-in commit prompts append:

```
Co-Authored-By: Claude Sonnet 4.6 <noreply@anthropic.com>
🤖 Generated with [Claude Code](https://claude.com/claude-code)
```

Gemini CLI and Codex CLI have analogous conventions. If we adopted
provider-specific attribution we would:

1. Leak which LLM produced which commit, which is irrelevant to the work.
2. Fragment history — a commit partly by Claude and partly by Gemini
   would accrete multiple `Co-Authored-By:` trailers.
3. Bind a user-facing convention to a single vendor. If the user
   switches providers, every existing commit suddenly misattributes.
4. Pollute `git shortlog` / `git log --pretty=format:'%an %ae'` with
   noreply addresses that have no relationship to the human doing the
   work.

**Decision.** The helper script hard-codes "no trailers, no LLM
mentions". The config key `ARCHITECT_GIT_ATTRIBUTION` exists only for
discoverability — the script ignores any value other than `false`.
Explicitly forbidden strings (tested post-implementation with a
grep-based audit):

- `Co-Authored-By: Claude`, `Gemini`, `Codex`, or any variant
- `Generated with Claude`, `Generated with Gemini`, `Generated with Codex`
- `🤖 Generated with [Claude Code]` and the markdown-link variants
- `Anthropic`, `Google`, `OpenAI` in trailers or body
- Emoji-only lines (`🤖`, `🧠`, `✨`) in subject or body

Each command file that drafts commit messages has this explicit note:

> Do NOT append any authorship trailers, Co-Authored-By lines, or
> "Generated with …" markers. The commit author is whatever git config
> produces; the harness is provider-agnostic.

**Consequences.**

- Every `commands/*.md` file gets a one-line reminder in its "Phase N:
  Commit" section. The helper script is the backstop.
- A post-implementation `grep` audit is part of the rollout checklist
  (`04-integration-guide.md § Rollout`).
- If a pre-commit hook in the user's repo tries to add a trailer, that
  is the user's repo policy; the harness does not fight it.

**Alternatives rejected.**

- *Per-provider attribution.* Rejected: user works across providers.
- *Opt-in attribution.* Rejected: the config key would become a foot-gun
  that drifts between provider defaults over time.
- *A single neutral trailer like `Co-Authored-By: Architect Harness
  <noreply@local>`.* Rejected: fake authorship is still fake. If the
  user wants a trailer they add it via git config `trailer.*`; the
  harness does not.

---

## ADR-002 — Commit after the human gate, never before

**Status:** ACCEPTED.

**Context.** Each command has human gates. Commits can live on either
side of a gate:

- *Before the gate:* artifacts land in history regardless of outcome.
  Rejection requires `git revert` or history rewrite. Noise accrues.
- *After the gate:* the working tree holds unapproved edits; `git
  restore` reverts on rejection. Only approved work becomes history.

**Decision.** Commit after the gate. This matches the existing
`/meta-apply` Phase 5 option 2 behaviour (`git restore
docs/{slug}/design/*.md` on rejection), which only works if there is no
intervening commit. Generalising: every command's commit step is the
*final* phase, after the last user approval.

**Consequences.**

- Rejection never pollutes history.
- Partial work within a command (e.g. architect verdict NEEDS_ITERATION)
  does not auto-commit; the user handles it in the working tree and
  reruns.
- `/meta-apply` keeps its existing rollback contract unchanged.

**Alternatives rejected.**

- *Commit before the gate, revert on reject.* History pollution from
  rejected work; bisect becomes harder.
- *Commit to a shadow branch, fast-forward on approve.* Over-engineered
  for a solo research workflow. The shadow-branch pattern belongs to CI
  systems, not editing helpers.

---

## ADR-003 — Helper script, not Claude Code hooks

**Status:** ACCEPTED.

**Context.** Claude Code supports hooks (PostToolUse, Stop,
UserPromptSubmit, etc.) that run shell commands in response to events.
An alternative architecture is: a hook fires after every assistant
turn, detects which slash command just ran, and commits.

**Decision.** Use an explicit helper script invoked from each command
file. No hooks.

**Rationale.**

- Hooks fire on *tool use events*, not on *command completion*. A slash
  command is many tool calls (agent dispatches, file writes, user
  prompts). A hook cannot easily tell "the command is done, gate
  passed" vs "one intermediate tool call finished".
- Hooks would have to scrape the recent conversation to recover the
  command name, feature slug, and gate outcome. Fragile and expensive.
- Command files already orchestrate phase-by-phase; adding a terminal
  "Phase N: Commit" phase is a natural extension that keeps all the
  context local.
- Hooks are provider-specific (Claude Code's feature). We want the
  harness portable across Claude / Gemini / Codex. A helper script is
  the smallest cross-provider unit.
- Debugging hooks is painful (they run out-of-band); debugging a script
  invocation is trivial (`--dry-run --verbose`).

**Consequences.**

- Each command file gains ~15 lines (one "Phase N: Commit" section).
  Centralised logic stays in `scripts/architect-commit.sh`.
- The script is portable to Gemini CLI and Codex CLI with no changes.
  Each CLI's corresponding command files gain the same phase.

**Alternatives rejected.**

- *PostToolUse hook that matches a tool-call pattern.* Fragile; see
  above.
- *Stop hook that commits on every assistant turn.* Wildly wrong — many
  turns are mid-phase and nothing should commit.
- *UserPromptSubmit hook.* Irrelevant; this fires at the wrong moment.

---

## ADR-004 — Default branch = current; feature branches opt-in

**Status:** ACCEPTED.

**Context.** Two common workflows:

1. *Single-branch* (often `main`): solo researcher, one feature at a
   time, linear history.
2. *Feature-branch*: multi-person or multi-feature-parallel work.

Binding the harness to either is wrong for the other.

**Decision.** Default `ARCHITECT_GIT_BRANCHES=false` — commit to current
branch. Users with parallel workflows flip it to `true`; the harness
then creates `architect/{slug}` on first commit for each feature.

**Rationale.**

- The existing DC_hum_verse repo (and pathway-explorer inside it) works
  on `main`. Forcing branches would create friction out of the gate.
- Feature-branch users know they want branches; one config flag is
  enough.
- Branch creation on *first commit* (not first command) means
  "`/map` that resulted in no write" doesn't churn branches for nothing.

**Consequences.**

- Single-branch users see no new branches unless they opt in.
- Multi-feature users get automatic isolation once they set one flag.
- Meta-layer branches are a separate knob (`ARCHITECT_GIT_META_BRANCH`)
  because portfolio work often rides on a specific feature branch.

**Alternatives rejected.**

- *Always branch.* Friction; surprises single-branch users.
- *Branch on `/map` regardless of whether anything was written.* Churn
  for no-op runs.
- *Branch per phase within a feature.* Absurd. Phases are commit units,
  not branch units.

---

## ADR-005 — No worktrees

**Status:** ACCEPTED.

**Context.** Worktrees give each feature an isolated checkout.
Superficially appealing for parallel agent dispatch.

**Decision.** Skip. Worktrees are a user-level concern, not a harness
concern.

**Rationale.**

- Parallel agent dispatch writes to distinct paths (`docs/{feat-A}/…`
  vs `docs/{feat-B}/…`). No file collision.
- `/meta-apply` parallel dispatch writes to per-feature subdirs in
  parallel — again, no collision.
- Worktrees force the user to `cd` into the right checkout before every
  command. The harness gets invoked in one place; multiple worktrees
  would mean multiple sessions.
- If a user really wants worktrees, nothing in the harness stops them:
  `git worktree add ../feat-x` + run the harness there. The helper
  script operates on whatever working tree it finds.

**Consequences.**

- One fewer moving piece.
- One fewer set of edge cases (stale refs, pruning, locking).

**Alternatives rejected.**

- *Worktree per feature.* See above.
- *Worktree per meta-apply run.* Creates a temp worktree, runs the
  parallel dispatch, merges. Over-engineered.

---

## ADR-006 — One commit per command invocation (two for `/design`)

**Status:** ACCEPTED.

**Context.** Granularity choices:

- *One commit per file written.* Very granular, potentially noisy.
- *One commit per phase within a command.* Matches internal phases but
  does not match the gate cadence.
- *One commit per command invocation.* Matches the "one user gate per
  command" rhythm. `/design` is the exception with two gates.
- *One commit per milestone* (e.g. map + review + synthesize batched).
  Loses bisect granularity.

**Decision.** One commit per command invocation, with `/design`
producing two (one per gate) and `/meta-apply` producing one per
approved feature plus a final architect-batch commit.

**Rationale.**

- Matches the harness's existing "one gate per command" invariant.
- Gives bisectable history: a regression introduced by a specific
  design or implement phase is isolated to one commit.
- `/design`'s two-gate exception is justified by the scope-vs-bundle
  separation. If scope is approved but bundle is later rejected, the
  scope commit still represents agreed state.

**Consequences.**

- 10–15 commits per full-feature walkthrough (see `02-behavior.md §
  Log-shape`).
- Users who find it noisy drop cheap commands from
  `ARCHITECT_GIT_COMMIT_ON`.

**Alternatives rejected.**

- *Squash on feature completion.* Lossy; precludes phase-level rollback.
- *Per-phase commits inside `/design`.* Phases inside `/design` are not
  individually gated; would violate the commit-after-gate rule.

---

## ADR-007 — Never `--no-verify`, never `--amend`, never `-f`

**Status:** ACCEPTED.

**Context.** Pre-commit hooks, gpg signing, and existing history
protection are repo-level policy. The harness could bypass them for
"smoother automation".

**Decision.** Absolutely not. The helper script never passes
`--no-verify`, `--no-gpg-sign`, `--amend`, or `-f`. On hook failure it
surfaces the error verbatim and exits 5.

**Rationale.**

- Hook bypass corrupts trust. Any linter, secret scanner, or format
  check the user put in place is there for a reason.
- `--amend` rewrites history, which can clobber signed commits or
  in-flight reviews.
- `-f` on push is out of scope anyway (we do not push), but even
  `git commit -f` is never the right answer.

**Consequences.**

- Hook failures halt the harness. User fixes the issue, reruns.
- The "create a NEW commit rather than amending" rule from git best
  practices is honored trivially — we never amend.

**Alternatives rejected.**

- *Opt-in bypass.* Rejected: would become a shortcut when hooks are
  annoying rather than wrong.

---

## ADR-008 — Stage only paths the command wrote

**Status:** ACCEPTED.

**Context.** `git add .` and `git add -A` are dangerous in a working
tree with unrelated edits (experimental code, scratch files, sensitive
data).

**Decision.** The helper script takes an explicit `--paths` list. The
calling command builds this list from the artifacts it just wrote. Never
`git add .` or `-A`.

**Rationale.**

- Users often have unrelated edits in progress.
- Sensitive paths (`.env`, credentials) are typically `.gitignore`'d,
  but a harness should not rely on gitignore alone.
- Explicit staging makes commits clean: a commit titled `docs(foo):
  plan` contains only `docs/foo/plan/**`.

**Consequences.**

- `02-behavior.md § Paths touched` becomes the authoritative table.
- `/implement` computes its path list from the phase doc's
  Files-to-Create + Files-to-Modify sections — no surprise staging.

**Alternatives rejected.**

- *`git add -u` (only modified tracked files).* Still too broad; stages
  unrelated modifications.
- *`git add` with a pathspec pattern like `docs/{feat}/`.* Good for
  artifact commands but fails `/implement`, where code can be anywhere.

---

## ADR-009 — No autopush, no autopull, no autotag

**Status:** ACCEPTED.

**Context.** Autopush on every commit is a common "productivity" bolt-on.

**Decision.** Never. `ARCHITECT_GIT_AUTOPUSH` exists for discoverability
only; the script ignores any value other than `false`.

**Rationale.**

- Push has remote side-effects: triggers CI, notifies reviewers, kicks
  webhooks. The user decides when to invite those consequences.
- Pull / fetch can reset working state unexpectedly.
- Tags are a human-curated concept; `/implement` finishing a phase is
  not a release.

**Consequences.**

- After a full feature walkthrough, the user runs `git push` manually
  when they're ready.
- CI integration, if any, is triggered by the user's push — not by the
  harness.

**Alternatives rejected.**

- *Autopush at phase N.* Opinionated and wrong for many workflows.
- *Autopush only on branch creation.* Still wrong; user may not want
  the branch public yet.

---

## ADR-010 — Config lives in `.claude/architect.git` as key=value

**Status:** ACCEPTED.

**Context.** Config format options: YAML (`.claude/architect.git.yaml`),
TOML, JSON, a shell-sourceable `.env`-style file.

**Decision.** `.claude/architect.git` in shell `KEY=value` form.

**Rationale.**

- The helper script is bash; `source .claude/architect.git` is one line.
- The file is tiny (<10 keys); YAML/TOML would be overkill.
- Provider-agnostic: Gemini CLI and Codex CLI wrapper scripts can read
  the same file.
- Easy to comment out / flip values without a schema.

**Consequences.**

- One config file, one format, one loader.
- No schema validation — the script treats unknown keys as ignored and
  missing required keys as "off by default".

**Alternatives rejected.**

- *Extend `settings.local.json`.* Would couple the harness's git policy
  to Claude Code's settings format. Against the provider-agnostic
  principle.
- *YAML.* Adds a dependency on a YAML parser in the script. Over-kill.

---

## ADR-011 — `/status` never commits; `/verify` does (by default)

**Status:** ACCEPTED.

**Context.** `/status` is chat-only with no file writes — trivially no
commit. `/verify` writes `docs/{feat}/verify.md`, which is a useful
artifact but also noisy if run frequently.

**Decision.**

- `/status`: hard-coded no commit (nothing to commit anyway).
- `/verify`: commits by default (goes in `ARCHITECT_GIT_COMMIT_ON`).
  Users who run verify as a frequent sanity check drop it from the
  allow-list.

**Rationale.**

- `/verify`'s output is a milestone artifact useful for handoff and CI
  pickup. Versioning it is useful.
- Users who find it noisy have a single-line config edit to silence it.

**Consequences.**

- `/verify` commits add ~1 entry per run. Grep-based dashboards
  (`git log --grep='verify'`) become a lightweight audit trail.

**Alternatives rejected.**

- *Default `/verify` off.* Fine, but users who want the audit trail
  would have to know to enable it.
- *Hardcode `/verify` on, not configurable.* Unnecessary rigidity.

---

## ADR-012 — Commit subject templates live in commands, not config

**Status:** ACCEPTED.

**Context.** Commit subject templates could live:

- Inside each command file (next to the phase that computes the values).
- Inside the config (`.claude/architect.git`).
- Inside the helper script (defaults per command name).

**Decision.** Templates live in the command files. The helper script
receives the fully-rendered subject as `--subject=…`.

**Rationale.**

- The command file is the only place that knows the relevant values
  (`{n_phases}`, `{n_adrs}`, `{reviewer_list}`, …) because it just
  computed them for its own report.
- Templating in config would require the helper script to parse values
  from arbitrary inputs — a new mini-language to maintain.
- The command file already has user-visible report text; reusing its
  phrasing in the commit message ensures chat-output and history
  agree.

**Consequences.**

- Each command file carries its subject templates in the "Phase N:
  Commit" section.
- Changing the template requires editing the command file — the right
  place, since the template reuses values unique to that command.

**Alternatives rejected.**

- *Config-driven templates with variable substitution.* New
  mini-language. Out of scope.
- *Helper-script defaults.* Loses per-command context.

---

## ADR-013 — `activate-role.sh` offers optional git opt-in

**Status:** ACCEPTED.

**Context.** The config file `.claude/architect.git` must come into
existence somehow. Options:

- Created manually by the user.
- Copied by `activate-role.sh` unconditionally.
- Offered by `activate-role.sh` behind a flag.

**Decision.** `activate-role.sh architect …` gains an optional
`--with-git` flag. When present, it copies
`templates/architect.git.template → .claude/architect.git` (never
overwriting). When absent, nothing changes — git automation stays off.

**Rationale.**

- Opt-in by design (see Motivation in `README.md`).
- Existing activations don't regress.
- Users who want it toggle one flag.

**Consequences.**

- One-line change to `activate-role.sh`.
- Users can also hand-copy the template file; the flag is a convenience,
  not a gate.

**Alternatives rejected.**

- *Create on every activation.* Violates opt-in.
- *Never create; users copy manually.* Friction; most users won't
  discover it.

---

## ADR-014 — Out of scope for v1

These are plausible features deliberately deferred:

- **`/architect-undo`** — convenience wrapper around `git revert HEAD`.
  Useful but not load-bearing; `git revert` is already one command.
- **Pre-phase stash** — auto-stash before `/implement`, pop on
  rollback. Attractive for safety but adds stash-recovery complexity.
  Defer until a concrete user incident argues for it.
- **Per-commit co-signing** — some repos require gpg signing or
  ssh-signing. Already handled by the user's git config; the harness
  just respects it.
- **Pre-commit message hooks** (`prepare-commit-msg`) — out of scope;
  users who want them wire them in git config.
- **Automatic `git tag`-ing of `/implement` final-phase commits.**
  Tagging is a human-curated concept; don't guess.
- **Automatic PR creation via `gh pr create`.** Vendor-locked to
  GitHub; not portable to GitLab/Bitbucket. User's choice.
- **`--dry-run` mode across the whole command.** The helper script has
  `--dry-run`; extending that to the command-level is possible but
  invites scope-creep.

All of the above can be added post-v1 without changing the core
architecture if the motivation becomes concrete.
