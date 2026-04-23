---
parent: ../README.md
---

# Git automation — plan for the architect harness

**Status:** proposal (pre-implementation). This directory holds the design
documents for tying git commits into the architect workflow without
over-engineering and without binding the harness to any particular LLM
provider.

---

## Problem

The architect harness produces many small, high-value artifacts —
`map.md`, `review/*.md`, `synthesis.md`, `design/*.md`, `plan/phase-NN.md`,
code under `/implement`. Today none of them are auto-committed. The only
git touch anywhere in the harness is one aspirational `git restore` line in
`/meta-apply` Phase 5 (rollback on rejection), and `Shell(git:*)` is
permitted in `.claude/settings.local.json` — but nothing invokes it.

Consequences:

- Users commit manually after "milestones they're happy with", which is
  fine for slow cadence but noisy for active multi-phase work.
- Rollback is working-tree-level only. `git restore` assumes uncommitted
  edits; once the user commits, the harness has no opinions.
- Sub-agents writing in parallel (e.g. `/review --as all`, `/meta-apply`)
  land atomically from the harness's view but as one big diff from git's
  view — no natural boundaries.
- Across providers (Claude, Gemini, Codex), commit-message conventions
  differ. Claude Code's built-in prompts inject `Co-Authored-By: Claude`
  and `🤖 Generated with Claude Code`. The user works across providers
  and does **not** want provider attribution baked into history.

## Non-goals

- Replacing the human gate. Commits are bookmarks for approved state, not
  a substitute for approval.
- Branching per phase. Features and (optionally) portfolio work earn
  branches; individual phases do not.
- Worktrees. Parallel dispatch already uses distinct paths
  (`docs/{feature-A}/…` vs `docs/{feature-B}/…`) — a worktree buys
  nothing.
- Pushing, pulling, merging, or opening PRs. That remains the user's
  choice and tool.
- "Smart" commit messages that summarise the diff. Templated messages from
  command context, full stop.

---

## TL;DR of the plan

A single provider-agnostic helper script + a small "Phase N: Commit"
addendum in each command file + a one-file config.

```
scripts/architect-commit.sh           ← one helper, invoked with command context
.claude/architect.git                  ← opt-in config (key=value; absent = off)
commands/*.md                          ← each gains a terminal "Phase N: Commit" step
```

Commit semantics:

| When | What | Why |
|------|------|-----|
| After a human gate approves an artifact | Commit the artifact paths just written | Commits represent approved state |
| After `/implement` phase verification passes | Commit the phase's code + the updated phase doc | Phase = commit unit |
| After `/meta-apply` compound gate approves some features | One commit per approved feature | Clean rollback if one is later reverted |
| When the user chose "use existing" at a pre-gate | No commit | Nothing changed |
| When the user rejected at a gate | No commit; `git restore` the paths | Rejection means working tree never becomes history |

Branch policy: default is **current branch**. Opt-in mode `architect/{slug}`
creates a feature branch on the first `/map` and suggests merge after the
last phase.

Attribution policy: **none, ever**. No `Co-Authored-By:`, no trailer
referencing Claude, Gemini, Codex, or any LLM. No emojis in commit
subjects. No `Generated with …` line in the body. The commit author is
whatever `git config user.{name,email}` returns — i.e. the human.

---

## Document index

| File | Purpose |
|------|---------|
| [01-architecture.md](./01-architecture.md) | Components: the helper script, the config file, command integration points, what each piece does and does not do. |
| [02-behavior.md](./02-behavior.md) | Per-command commit flow, commit-message catalog, branch & rollback semantics, failure modes. |
| [03-decisions.md](./03-decisions.md) | ADRs: provider-agnostic attribution, commit-after-gate, no-branch default, script-not-hook, one-commit-per-command, no worktrees, no autopush, etc. |
| [04-integration-guide.md](./04-integration-guide.md) | Concrete diffs to apply to each `commands/*.md` file. The helper-script spec with interface, exit codes, and flag contract. Rollout order. |

---

## Decision tree for adopting this

```
Do you want the harness to auto-commit?
├── No  → stop. Keep manual commits. (Do nothing.)
└── Yes
    ├── Do you work across LLM providers (Claude, Gemini, Codex)?
    │   └── The attribution-off policy is mandatory regardless.
    │       Every provider must produce the same history.
    │
    ├── Do you run multiple features in parallel?
    │   ├── Yes → enable `ARCHITECT_GIT_BRANCHES=true`
    │   └── No  → leave branches off; commit to current branch
    │
    └── Do you want granular per-artifact commits or batched milestones?
        ├── Granular (default) → one commit per command invocation
        └── Batched            → set `ARCHITECT_GIT_COMMIT_ON=implement,meta-apply`
                                 to only auto-commit expensive phases; use
                                 manual commits for cheap artifacts (map,
                                 review, synthesize).
```

---

## What "done" looks like after implementation

- `scripts/architect-commit.sh` exists, is executable, and passes a
  smoke test (`./scripts/architect-commit.sh --selftest`).
- `.claude/architect.git.template` exists with documented keys; activation
  copies it to `.claude/architect.git` only if the user opts in.
- Every command file under `commands/` has a terminal "Phase N: Commit"
  section (or an explicit "This command does not commit" note for the
  handful that don't produce artifacts, e.g. `/status`).
- A worked example exists under `docs/workflows/architect/git/examples/`
  showing the resulting `git log --oneline` for a small feature taken
  through `/map → /review → /synthesize → /design → /plan → /implement 1`.
- The integration guide in `04-integration-guide.md` has been applied and
  the diff is reviewable; no command file gained more than ~30 lines.

---

## Status

Proposal — written 2026-04-18. No code or command changes have been made
yet. Review, revise, then execute the plan in `04-integration-guide.md`.
