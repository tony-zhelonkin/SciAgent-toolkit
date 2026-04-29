---
parent: ../README.md
version: 2
revised: 2026-04-29
---

# Git automation — plan for the architect harness

**Status:** proposal v2 (pre-implementation). v1 architectural review identified
over-engineering and hidden-state risks; this revision trims the scope and
introduces a renderer layer + transparency log. Diff vs v1 summarised at the
bottom of this file.

This directory holds the design documents for tying git commits into the
architect workflow without over-engineering and without binding the harness
to any particular LLM provider.

---

## Problem (unchanged from v1)

The architect harness produces many small, high-value artifacts —
`map.md`, `review/*.md`, `synthesis.md`, `design/*.md`, `plan/phase-NN.md`,
code under `/implement`. Today none of them are auto-committed. Users
commit manually after "milestones they're happy with", which is fine for
slow cadence but noisy for active multi-phase work. Rollback is
working-tree-level only. Across providers (Claude, Gemini, Codex),
commit-message conventions differ and would inject vendor-specific
authorship trailers if left to default.

## Non-goals (unchanged from v1)

- Replacing the human gate. Commits are bookmarks for approved state, not
  a substitute for approval.
- Branching per phase. Features earn branches (opt-in); individual phases
  do not.
- Worktrees. Parallel dispatch already uses distinct paths
  (`docs/{feature-A}/…` vs `docs/{feature-B}/…`).
- Pushing, pulling, merging, or opening PRs. That remains the user's
  choice and tool.
- "Smart" commit messages that summarise the diff.

## Additional non-goals (new in v2)

- Per-command attribution-policy reminders in markdown. The `guard()`
  filter in the renderer is the single enforcement point; reminders in
  13 command files are token noise.
- Multi-commit ceremonies. `/design` ships one commit (Gate 2), not two
  or three. `/meta-apply` ships one atomic commit per invocation, not N
  per feature. Bisect granularity stays at the command level.

---

## TL;DR of the v2 plan

A 4-layer + 1 cross-cut design. All git logic in one bash script; all
message composition in one Python module with unit tests; one config file;
one transparency log; ~3-line trigger blocks in command files.

```
scripts/architect-commit.sh           ← L2: git operator (bash, no policy)
scripts/_render.py                    ← L3: subject/body renderer + guard() filter
scripts/_render_test.py               ← pytest: 13-command coverage + forbidden-string
templates/architect.git.template      ← L1: opt-in config; absent = off
commands/*.md                         ← L4: terminal "Phase N: Commit" trigger (~3 lines)
.claude/architect.git                 ← L1: runtime config (copied from template)
.claude/architect.git.log             ← L0: append-only transparency log
```

Commit semantics:

| When | What | Why |
|------|------|-----|
| After a human gate approves an artifact | Commit the artifact paths just written | Commits represent approved state |
| After `/implement` phase verification passes | Commit the phase's code + the updated phase doc | Phase = commit unit |
| After `/meta-apply` compound gate approves features | One atomic commit covering all approved features | Hook failure leaves working tree intact |
| When the user chose "use existing" at a pre-gate | No commit | Nothing changed |
| When the user rejected at a gate | No commit; `git restore` the paths after dirty-path safety check | Rejection means working tree never becomes history |

Branch policy: default is **current branch**. Opt-in mode `architect/{slug}`
creates a feature branch on the first `/map` and announces the switch
visibly. Meta-layer always commits on current branch (no `META_BRANCH=portfolio`
mode in v2).

Attribution policy: **none, ever**. No `Co-Authored-By:`, no trailer
referencing Claude, Gemini, Codex, or any LLM. No emojis in commit
subjects. No `Generated with …` line in the body. The commit author is
whatever `git config user.{name,email}` returns. Enforced at two
independent points: (a) `guard()` filter in the renderer, (b) hard-coded
no-trailer composition in the helper script.

Transparency: every helper invocation appends one line to
`.claude/architect.git.log` (timestamp | command | feature | verdict |
sha-or-reason). The user can `tail -20 .claude/architect.git.log` to
audit everything the harness ever did to git.

---

## Document index

| File | Purpose |
|------|---------|
| [01-architecture.md](./01-architecture.md) | The 4 layers + cross-cut. What each piece does and does not do. |
| [02-behavior.md](./02-behavior.md) | Per-command commit flow, commit-message catalog, branch & rollback semantics, transparency log format, failure modes. |
| [03-decisions.md](./03-decisions.md) | ADRs. v1 ADR-006 and ADR-012 superseded; new ADRs for the renderer layer, transparency log, atomic /meta-apply, dirty-path refusal. |
| [04-integration-guide.md](./04-integration-guide.md) | Concrete diffs (≤3 lines per command file). Renderer module spec. 5-PR rollout. |

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
    ├── Do you run multiple features in parallel on shared branches?
    │   ├── Yes → enable ARCHITECT_GIT_BRANCHES=true
    │   └── No  → leave branches off; commit to current branch
    │
    └── Do you want granular per-artifact commits or batched milestones?
        ├── Granular (default) → one commit per command invocation
        └── Batched            → set ARCHITECT_GIT_COMMIT_ON=implement,meta-apply
                                 to only auto-commit expensive phases
```

---

## What "done" looks like after implementation

- `scripts/architect-commit.sh` exists, executable, passes `--selftest`.
- `scripts/_render.py` exists; `pytest scripts/_render_test.py` passes
  (covers 13 commands × happy-path + forbidden-string-rejection).
- `templates/architect.git.template` exists with documented keys;
  activation copies it to `.claude/architect.git` only on `--with-git`.
- Every command file under `commands/` has a ≤3-line "Phase N: Commit"
  trigger (or an explicit "This command does not commit" note for `/status`).
- Worked-example log under `docs/workflows/architect/git/examples/`
  showing the resulting `git log --oneline` for one full feature run.
- `.claude/architect.git.log` accumulates one line per invocation; user
  can audit harness git activity without reading source.
- No command file gained more than ~3 lines.

---

## Diff vs v1 (for reviewers)

| Concern | v1 disposition | v2 disposition | Reason |
|---|---|---|---|
| Subject/body composition | In each command file (heredoc with `<(printf ...)`) | In `scripts/_render.py` reading frontmatter | ~325 lines of bash-in-markdown removed; testable; single chokepoint for `guard()` |
| `/design` commits | A (scope) + B (bundle) + C (verdict) = 3 | One at Gate 2 (subsumes A and C) | Bisect granularity already at command level; 3 SHAs for one approval is noise |
| `/meta-apply` commits | N per-feature commits + 1 architect-batch | One atomic commit per invocation | Hook failure mid-loop in v1 leaves half-committed state; v2 = all-or-nothing |
| `ARCHITECT_GIT_META_BRANCH=portfolio` | Supported with cross-stream refusal logic | Removed; meta always on current branch | Edge case; complexity for ~zero users |
| Forbidden-string CI grep | CI job over scripts + command files | Pytest unit tests on renderer output | No CI assumed; tests run locally; faster feedback |
| Per-command attribution reminder | One line per command file | Removed — `guard()` is the enforcement | 13 lines saved; single source of truth |
| Skip-detection ("Skip when nothing was written") | Repeated in command file prose | Single check in helper (empty staging → exit 0) | Don't ask command authors to re-derive |
| Branch switch on first commit | Silent | Visible announcement: `architect-commit: switched to architect/{slug} (was main)` + log line | Hidden-state mitigation |
| `git restore` on rejection | Direct | Refused if target paths have uncommitted hand-edits; user must stash first | Hidden-state mitigation |
| Transparency log (L0) | — | New: `.claude/architect.git.log` append-only | Hidden-state mitigation |
| Rollout PR count | 8 | 5 | Reduced scope merits compression |

The architectural bones from v1 (single helper script, opt-in config,
attribution hard-off, after-gate commits, explicit `--paths`, no
`--no-verify`, no autopush) are preserved verbatim.

---

## Status

Proposal v2 — written 2026-04-18, revised 2026-04-29. No code or command
changes have been made yet. Review, revise, then execute the plan in
`04-integration-guide.md`.
