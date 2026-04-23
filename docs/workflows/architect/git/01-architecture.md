---
parent: ./README.md
view: structural
---

# Architecture — git automation

Structural view: the components, where they live, what each touches, and
what each explicitly does not touch.

---

## Component map

```mermaid
graph TD
    subgraph "SciAgent-toolkit/ (source)"
        SCRIPT[scripts/architect-commit.sh<br/>single helper, all commit logic]
        TEMPLATE[templates/architect.git.template<br/>opt-in config defaults]
        subgraph "commands/"
            CMD_MAP[map.md]
            CMD_REV[review.md]
            CMD_SYN[synthesize.md]
            CMD_DES[design.md]
            CMD_ARCH[architect.md]
            CMD_PLN[plan.md]
            CMD_IMP[implement.md]
            CMD_VRF[verify.md]
            CMD_STAT[status.md]
            CMD_MM[meta-map.md]
            CMD_MD[meta-design.md]
            CMD_MA[meta-apply.md]
            CMD_MP[meta-plan.md]
        end
    end

    subgraph "Project/.claude/ (runtime)"
        CFG[architect.git<br/>key=value config<br/>absent = off]
        SET[settings.local.json<br/>Shell(git:*) already allowed]
    end

    subgraph "Project/ (runtime, git-tracked)"
        DOCS[docs/*/map.md, review/, design/, plan/]
        CODE[src code written by /implement]
    end

    CMD_MAP -->|invokes at Phase N| SCRIPT
    CMD_REV --> SCRIPT
    CMD_SYN --> SCRIPT
    CMD_DES --> SCRIPT
    CMD_ARCH --> SCRIPT
    CMD_PLN --> SCRIPT
    CMD_IMP --> SCRIPT
    CMD_VRF --> SCRIPT
    CMD_MM --> SCRIPT
    CMD_MD --> SCRIPT
    CMD_MA --> SCRIPT
    CMD_MP --> SCRIPT

    SCRIPT -.reads.-> CFG
    SCRIPT -.stages+commits.-> DOCS
    SCRIPT -.stages+commits.-> CODE

    CMD_STAT -.no commit.-> CFG

    TEMPLATE -.copied into.-> CFG
```

Everything new is in one script, one config template, and one "Phase N:
Commit" block per command file. No hooks, no agents, no wrappers.

---

## Components

### 1. `scripts/architect-commit.sh` — the only executable

**Responsibility:** stage a set of paths, compose a templated commit
message from command context, and run `git commit`. One script per
operation. No loops, no interactive prompts, no git config writes.

**Inputs (CLI contract):**

```
architect-commit.sh \
  --command=<map|review|synthesize|design|architect|plan|implement|verify|meta-map|meta-design|meta-apply|meta-plan> \
  --feature=<slug>                 # or __meta for meta-layer commands
  --phase=<N>                      # /implement only; omit otherwise
  --subject=<free-text subject>    # command-specific, see 02-behavior.md
  --body=<free-text body>          # optional, multi-line OK (stdin alternative: --body-file)
  --paths=<p1,p2,…>                # explicit list; if absent, derive from --command + --feature
  [--dry-run]                      # print what would happen, don't touch git
  [--verbose]
```

**Behaviour:**

1. Read `.claude/architect.git` if present; otherwise treat as "off" and
   print a single-line skip notice (`architect-commit: disabled (no
   .claude/architect.git); skipping`) and exit 0.
2. Refuse to run outside a git worktree (`git rev-parse --is-inside-work-tree`).
   Exit non-zero with a clear message.
3. If `ARCHITECT_GIT_BRANCHES=true` and the current branch does not match
   `${ARCHITECT_GIT_BRANCH_PREFIX}{feature}`, create and switch to it *only
   on the first commit for that feature* (see §4 in `02-behavior.md`).
   Otherwise leave branch alone.
4. Stage the requested paths. **Never `git add .` or `git add -A`.** Only
   the explicit list, which is always inside `docs/{feature}/` or (for
   `/implement`) inside the phase's declared Files-to-Create /
   Files-to-Modify list. For meta-layer commands, scope is
   `docs/_meta/**` plus (for `/meta-apply`) the touched per-feature
   design files.
5. If nothing is staged after step 4, exit 0 with `architect-commit:
   nothing to commit; skipping`. This handles "user chose 'use existing'"
   and partial no-op scenarios.
6. Compose the commit message from `--subject` and `--body`. Append no
   trailers. No `Co-Authored-By`, no `Signed-off-by` unless the user's
   git config itself produces one (we do not add to user config).
7. `git commit` using a here-doc or `-F -`. Respect pre-commit hooks.
   **Never** pass `--no-verify`, `--no-gpg-sign`, or `--amend`.
8. On commit failure (hook rejection, gpg, etc.), do NOT retry; print the
   git output and exit non-zero so the calling command can surface it to
   the user.

**Exit codes:**

| Code | Meaning |
|------|---------|
| 0    | Commit made, or intentionally skipped (disabled / nothing to stage / --dry-run) |
| 2    | Bad input (missing required flag, unknown command) |
| 3    | Not in a git worktree |
| 4    | Staging failed (path outside workspace, etc.) |
| 5    | `git commit` itself failed (hook, signing, etc.) — surface output |

**What the script does not do:**

- Does not push. Remote state is the user's responsibility.
- Does not pull, fetch, merge, or rebase.
- Does not create tags, stashes, or worktrees.
- Does not write to `~/.gitconfig`, `.git/config`, or `.gitignore`.
- Does not parse diffs to "summarise changes" for the message. The calling
  command supplies the subject and optional body; that is the only source
  of message text.
- Does not prompt interactively. Stdin is only used for `--body-file -`.

---

### 2. `.claude/architect.git` — single config file

**Format:** `KEY=value` lines, `#` comments, bash-sourceable. Same form
as a `.env` file so `scripts/architect-commit.sh` can `source` it.

**Template** (`templates/architect.git.template`, copied into
`.claude/architect.git` only if the user opts in during
`activate-role.sh`):

```bash
# Architect git automation — local config
# Delete this file to disable auto-commit entirely.

# Master switch. If false or absent, the helper script no-ops.
ARCHITECT_GIT_AUTOCOMMIT=true

# Which commands may auto-commit. Comma-separated. Default: all.
# Drop expensive-to-review commands here to get batched manual commits
# for them (e.g. "implement,meta-apply" only).
ARCHITECT_GIT_COMMIT_ON=map,review,synthesize,design,architect,plan,implement,verify,meta-map,meta-design,meta-apply,meta-plan

# Create a feature branch (prefix + slug) on first commit for that
# feature. Leave false for single-branch workflows.
ARCHITECT_GIT_BRANCHES=false
ARCHITECT_GIT_BRANCH_PREFIX=architect/

# Meta-layer branch behaviour when ARCHITECT_GIT_BRANCHES=true.
# "portfolio" = single shared branch "architect/_meta"; "current" = commit on current branch.
ARCHITECT_GIT_META_BRANCH=current

# Safety rails — do not turn off without reading 03-decisions.md.
ARCHITECT_GIT_ATTRIBUTION=false   # HARDCODED false — the script ignores any other value.
ARCHITECT_GIT_AUTOPUSH=false      # Same — ignored if flipped.
```

**What the config is not:**

- Not a per-command customisation of commit messages. Message templates
  live inside each command file and are passed to the script as
  `--subject`/`--body`. Keeping templates in commands (where context
  lives) avoids an out-of-band second source of truth.
- Not a rules engine. The config only answers "commit or not" and "branch
  or not". Anything more complex is a smell.

---

### 3. Command-level integration — "Phase N: Commit"

Each `commands/*.md` that writes an artifact gains a terminal phase:

```markdown
## Phase N: Commit (if enabled)

After the human gate approves this command's output, invoke:

`scripts/architect-commit.sh --command=<name> --feature={slug} --subject=... --paths=...`

The script no-ops if `.claude/architect.git` is absent or the command is
not listed in `ARCHITECT_GIT_COMMIT_ON`. See
`docs/workflows/architect/git/02-behavior.md` for the per-command commit
subject and path list.
```

The integration is intentionally thin:

- The command file contains the subject template and the path expression.
- Every substitution (`{slug}`, `{phase_name}`, `{adr_count}`) is resolved
  from data the command already computed for its own report step.
- No command invokes `git` directly. Only via the helper script.

See `04-integration-guide.md` for the exact phrasing to paste into each
command file.

---

### 4. What already exists and is reused

| Existing | Purpose | What we use it for |
|----------|---------|---------------------|
| `.claude/settings.local.json` with `Shell(git:*)` | Permission to run any git subcommand | Permits the helper script without new grants. |
| `/meta-apply` Phase 5 option 2: `git restore docs/{slug}/design/*.md` | Rollback on rejection | Unchanged — still applies because we commit *after* the gate. |
| Phase-file "Verification" checklists in `plan/phase-NN.md` | Acceptance criteria per phase | The `/implement` Phase 4 verdict is the gate for the Phase 5 commit step. |
| `activate-role.sh` symlink installer | Role/agent/skill activation | Gains an optional "copy `.claude/architect.git.template → .claude/architect.git`" step behind a `--with-git` flag. |

No existing agent, MCP server, or skill is modified.

---

## Paths touched, by command

The helper script stages only the paths listed below. All are relative to
the git root.

| Command | Paths staged |
|---------|-------------|
| `/map {feat}` | `docs/{feat}/map.md` |
| `/review {feat} --as <set>` | `docs/{feat}/review/*.md` (all files mtime-newer than the command start) |
| `/synthesize {feat}` | `docs/{feat}/synthesis.md` |
| `/design {feat}` (Gate 1) | `docs/{feat}/design/README.md` |
| `/design {feat}` (Gate 2) | `docs/{feat}/design/*.md` (the full bundle, status=APPROVED post-flip) |
| `/architect {feat}` | `docs/{feat}/design/review.md` + any `design/*.md` whose status was flipped |
| `/plan {feat}` | `docs/{feat}/plan/README.md`, `docs/{feat}/plan/phase-*.md` |
| `/implement {feat} {N}` | All Files-to-Create + Files-to-Modify from `docs/{feat}/plan/phase-N.md` + the phase doc itself (for the "Completed" marker update) |
| `/verify {feat}` | `docs/{feat}/verify.md` |
| `/meta-map` | `docs/_meta/map.md` |
| `/meta-design` | `docs/_meta/design.md`, `docs/_meta/deferred.md` (if new rows appended) |
| `/meta-apply` | **Per approved feature:** `docs/{feat}/design/*.md` + relevant slice of `docs/_meta/deferred.md`. One commit per feature — see `02-behavior.md § /meta-apply`. |
| `/meta-plan` | `docs/_meta/plan.md` |
| `/status` | **No commit.** Chat-only output. |

---

## Invariants

1. **The helper script is the only place that runs `git`.** Command files
   never shell out to git directly; agents never run `git add`/`commit`.
2. **Attribution is off, always.** `ARCHITECT_GIT_ATTRIBUTION` is listed
   in the config for discoverability only — the script ignores it and
   hard-codes "no trailers, no LLM mentions". See ADR-001 in
   `03-decisions.md`.
3. **Commit happens after the gate.** The human approval precedes the
   commit in every command. Rejection means the working tree is reverted
   (never a post-factum `git revert`).
4. **Paths staged ⊆ paths the command just wrote.** No "while we're here"
   staging of unrelated files.
5. **No --no-verify, ever.** Pre-commit hook failures surface to the
   user; the harness does not bypass them.
6. **No empty commits.** If nothing was written, nothing is committed.
