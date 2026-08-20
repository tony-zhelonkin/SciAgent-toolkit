# Inspection brief — what `docs/_internal/` actually looks like in the field

You are **surveying, not designing and not fixing**. Read only. Do not modify,
create, or delete any file in any repository. Your entire output is the report
described in "Deliverable".

## Why this survey exists

`scio` (the AI harness toolkit vendored into each analysis project) ships an
always-on instruction in every project's `AGENTS.md` saying:

> Durable project memory lives in tracked files alone — `AGENTS.md`, committed
> stages, decisions in `docs/_internal/reasoning/`, handoffs in
> `docs/_internal/sessions/`. Harness auto-memory is off and untracked harness
> state counts as no record.

It simultaneously writes `docs/_internal/` into the `SCIO:GITIGNORE` managed
block of every project's `.gitignore`, and its `docs-layout` lint check warns
when `docs/_internal/` is **not** gitignored. So the mandated memory location is
guaranteed untracked, and by the instruction's own definition counts as no
record. Spot-checked before this survey: 14761-DM has 1,161 files and 13 MB
under `docs/_internal/`, zero tracked. DC-nexus 4,169 files, 2 tracked.

The owner intends to stabilize this into one skeleton that every project shares
and that any harness (Claude Code, Codex, others) can rely on. **Before**
designing that, we need to know what actually accumulated — the real layout, the
real names, the real drift.

## What to inspect

For **each project path given in your assignment**, report on:

1. **`docs/_internal/` layout.** Immediate subdirectory names and the file count
   in each. Note depth: does content sit directly in the subdir, or nested
   further (e.g. `plans/<date-slug>/NN_*.md`)?
2. **Naming drift.** Across projects, the same concept appears under different
   spellings — `handoff/` vs `handoffs/`, `session/` vs `sessions/`,
   `reasoning/` vs `decisions/` vs `reports/`. Enumerate every variant you see
   and which project uses it. This is the single most important output: it
   quantifies the drift the owner wants to stop.
3. **Tracked vs untracked.** `git -C <proj> ls-files docs/_internal | wc -l`
   against the on-disk file count. Also check what `git check-ignore -v` says,
   so we know whether the ignore came from the managed `SCIO:GITIGNORE` block or
   was hand-written.
4. **What agents actually wrote.** Sample 3–6 files per project (do not read
   them all). For each sample: is it a durable decision record, a session
   handoff, a plan, a one-off scratch probe, or a stale artifact superseded long
   ago? Give a rough proportion per project. We need to know the signal-to-noise
   ratio, because a memory system that ships everything is a landfill.
5. **Scratch behaviour.** Look for `_scratch/` at the project root and for
   evidence of agents writing to `/tmp/` instead (references to `/tmp/` paths
   inside `docs/_internal/` files, `02_analysis/`, or notebooks). The owner has
   observed agents writing probe scripts to
   `/tmp/claude-<pid>/<hashed-workspace>/<uuid>/scratchpad/` — outside the
   project entirely. Report any trace of that pattern.
6. **Harness footprints.** Which of `.claude/`, `.codex/`, `.agents/`,
   `.gemini/`, `AGENTS.md`, `CLAUDE.md` exist. Whether any harness-specific
   directory contains files an agent authored (as opposed to scio-managed
   mounts, which are symlinks). The goal is one skeleton that is harness-neutral,
   so harness-specific accumulation is a finding.
7. **Anything that surprised you.** Content in an unexpected place, evidence of a
   convention nobody documented, or a project that solved something well.

## Constraints

- **Read-only.** Never write. Do not run `git add`, `git commit`, or any command
  with `--go`, `--apply`, `--fix`.
- **Bound your effort.** These trees are large. Use `find`, `ls`, `wc`,
  `git ls-files`, `du`. Sample files; never read a whole subtree.
- **Report absence explicitly.** "No `docs/_internal/` at all" is a real and
  useful finding. So is "directory exists but is empty".
- Do not propose a design. A separate consultation does that. If you have a
  design opinion, put it in one short "Observations for the designer" section at
  the end and keep it to what your evidence supports.

## Deliverable

Markdown. Lead with a table: one row per project, columns = `docs/_internal`
file count, tracked count, subdirectory names present, harness dirs present,
`_scratch/` present. Then per-project detail only where that project differs
from the pattern. Then a consolidated **naming-variant table** (concept → every
spelling observed → projects using it). Then signal-to-noise estimates. Then
"Observations for the designer".

State plainly what you could not determine and why. An honest gap is worth more
than a confident guess.
