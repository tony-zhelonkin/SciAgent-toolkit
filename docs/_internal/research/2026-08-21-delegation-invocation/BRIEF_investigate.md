# Investigation brief — why claude→codex delegation works in one repo and fights in another

Read-only. Modify nothing. Your output is the report in §5.

## 1. What happened

The owner drives implementation by having Claude Code author prompts for codex
workers, then running them himself with a `!` prefix. In **JR-MC** (created
2026-08-21, pinned at the current toolkit tip) this failed repeatedly. In
**Meta-Aging** it works. Both use the same toolkit, the same container image,
the same operator.

Observed JR-MC failures, in order:

1. A long inline invocation **wrapped in the terminal**, splitting `-m
   gpt-5.6-sol` so codex reported *"a value is required for `--model <MODEL>`"*
   and bash then tried to execute `gpt-5.6-sol` as a command.
2. `--search` was passed and **does not exist**. Verified: absent from
   `codex exec --help` in both 0.147.0 (host) and 0.149.0 (the JR-MC container).
   Web search is config-gated (`tools.web_search`, `--enable
   web_search_request`), not a flag.
3. Recovery required several rounds of `strings` on the codex binary to discover
   the real capability surface.

The operator's own read: *"The line wrapped and broke `-m gpt-5.6-sol`."*

## 2. What the field built in response — three different conventions

| Project | Prompts live | Output lives | Invocation |
|---|---|---|---|
| **JR-MC** | `docs/_internal/codex/{00_house_rules.md,<unit>.md}` | `logs/codex/<unit>.{final.md,run.log}` | `.devcontainer/scripts/codex_unit.sh <unit> [--web]` — a tracked-looking wrapper pinning model and effort |
| **Meta-Aging** | `integration/03_results/_scratch/codex_handoff/<UNIT>_prompt.md` | `<UNIT>_run_last.md` beside it | long inline `codex exec … - < prompt.md` |
| **14782-DM** | `/tmp/claude-<pid>/…/scratchpad/codex/` | same | `probe.sh` + `launch.sh` in that scratchpad |

And the toolkit's own `skills/delegate-cli/SKILL.md` prescribes a **fourth**:
`/tmp/<unit>_prompt.md`, `/tmp/<unit>_last.md`, `/tmp/<unit>_run.log`.

Two facts already established, do not re-derive:

- JR-MC's `codex_unit.sh` and its whole `docs/_internal/codex/` tree are
  **untracked** (`git ls-files` returns nothing for either). `docs/_internal/`
  is gitignored by the toolkit's own `SCIO:GITIGNORE` block, so the house rules
  and every unit prompt are invisible to git by construction. The fix the
  project invented cannot survive the container.
- Meta-Aging has **no wrapper script at all**.

## 3. Your assignment

Investigate **only the project paths given in your prompt**. Establish, with
evidence:

1. **Why the invocation failed or held.** Command length and shape; whether
   prompts are passed by stdin, by `-` , or inline; absolute vs relative paths;
   whether a wrapper absorbs the length. Measure actual character counts of the
   invocations you find in logs, scripts, and tracked docs — a wrapped line is a
   length problem, so quantify it.
2. **Version and capability drift.** Which codex version this project's
   container runs, which flags its docs/scripts/prompts assume, and where those
   disagree. Check `~/.codex/config.toml`, any `.codex/` in the project, and any
   project doc naming a codex flag.
3. **What the project invented, and whether it survives.** For each artifact:
   is it tracked, gitignored, or in `/tmp`? Would it exist after a container
   recreate? (Only `/workspaces` — the bind mount — survives; home and `/tmp`
   do not.)
4. **Whether `delegate-cli` was reachable.** Does `.claude/skills/` exist? Is
   it a symlink, and does it resolve *from the container's perspective* as well
   as the host's? Count reachable skills against the count in the project's own
   pinned `01_modules/*/skills/`. Note: the six category links are written
   **absolute**, so a container-bound project's links dangle on the host and
   vice versa — say which side yours resolve on.
5. **What the operator had to discover by execution** that a reachable,
   accurate skill would have told them.

## 4. Constraints

- **Read-only.** No writes, no `git add`, no `--go`/`--apply`/`--fix`.
- Bound the effort: `ls`, `find`, `wc`, `grep`, `git ls-files`,
  `git check-ignore`, `readlink`, `codex exec --help`. Sample; do not read whole
  trees.
- Report absence explicitly. "No wrapper exists" is a finding.
- **Do not design the fix.** A synthesis step does that. If you have a design
  opinion, put it in one short closing section limited to what your evidence
  supports.

## 5. Deliverable

1. A table: artifact → path → tracked/ignored/tmp → survives container
   recreate?
2. The invocation failure mechanism for your project, with measured command
   lengths.
3. Version/capability drift: what is installed vs what the project's own text
   assumes.
4. Skill reachability numbers, and which side the category links resolve on.
5. **What a correct, reachable `delegate-cli` would have prevented** — be
   specific about which failure.
6. Anything that surprised you; anything you could not determine and why.

State plainly what you verified versus inferred.
