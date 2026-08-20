# Consultation — where delegation artifacts live, and who normalizes it

Consulting on **placement and ownership**. Produce a written recommendation.
Modify nothing.

Read first: `FINDINGS_field.md`, `CONSULT_enforcement.report.md`, and
`../../plans/2026-08-20-memory-and-seams/00_INDEX.md` (the decisions already
taken). Then read `skills/delegate-cli/SKILL.md` in full — 416 lines — because
three findings below are inside it.

## 1. The workflow in question

The owner drives implementation by having Claude Code author **runnable script
prompts** for codex workers, which the owner then inspects and executes with a
`!` prefix. A real instance, from 14782-DM:

- `scratchpad/codex/probe.sh` — capability check: does codex run, which model
  answered, is tool use live (it verified by reading the real first line of
  `AGENTS.md` rather than answering from memory), exit status.
- `scratchpad/codex/launch.sh` — blocks until three workers exit, prints each
  worker's byte counts.
- Three worker prompt files — a shared spec (no git writes, no network, no
  dependency changes, `00_data/` read-only, never open the 18.8 GB panel) plus
  **disjoint path ownership** per worker, plus per-worker traps (one worker's
  do-not-touch list names `install_local_packages.sh`, which looks exactly like
  the file to edit and must stay byte-identical), plus acceptance commands the
  orchestrator ran itself first.

All of it under `/tmp/claude-<pid>/<hashed-workspace>/<uuid>/scratchpad/codex/`.
The reasoning that produced it — a venv-from-venv `.pth` fix, `pathlib`'s `**`
descending into dot-directories and polluting a test glob, a released tag
existing only on GitHub, and a live decision to park a frozen-label-space
question rather than fold it into this branch — exists only in a chat
transcript. The owner's question: does this belong in a durable place, and is
normalizing it scio's job?

## 2. Three things already true — verify each, they reframe the question

**(a) The skill already contains the pattern the field believes it invented.**
`skills/delegate-cli/SKILL.md:152-160`, "Three rules that hold for every
implementer, learned the expensive way": put exact numeric targets in the prompt
and say *report a deviation rather than editing the artifact to match*; forbid
weakening or deleting a test to make a change fit, and ask for the disposition of
every test removed, by name; review is not optional. The owner described the
first of these as one of "the two rules that earn their keep." Line 211 also
already documents the `!`-prefix division of labour: *"Anton runs in the `!`
prefix; I prepare the prompt file and monitor output."*

**(b) scio itself prescribes `/tmp`.** Lines 104, 139-140: `/tmp/<unit>_prompt.md`,
`/tmp/<unit>_last.md`, `/tmp/<unit>_run.log`. The disappearance the owner is
complaining about is not agent drift — it is the toolkit's own documented
convention. Assess this directly: is `/tmp` right for some of these artifacts and
wrong for others, or wrong for all?

**(c) The Meta-Aging umbrella has ZERO skills mounted.** `.claude/skills` does
not exist there, while 86 skills sit in its own pinned toolkit checkout. So its
`AGENTS.md` grew a hand-written "Delegating to codex / agy" section that
restates the skill's container-sandbox findings — not because the field invented
a convention, but because it could reach nothing and **reconstructed** one.

That section also appears to hold knowledge the skill lacks: that `agy --model`
takes an exact label or slug from `agy models` while `--effort` is rejected for
every model, and that `agy` takes its prompt as the *value* of `-p` so flags
must follow it because a bare `--print` swallows the next flag. Verify whether
the skill has these. If field knowledge is genuinely newer than the catalog,
say what should carry it upward.

## 3. The hypothesis to attack

The delegation run produces three artifact classes, and they already have homes
under decisions taken — **no new location needs naming**:

| Artifact | Nature | Proposed home |
|---|---|---|
| `probe.sh` | near-identical every run; a capability check | a **toolkit asset** under `skills/delegate-cli/assets/`, mounted via the catalog |
| `launch.sh` | boilerplate driver + a per-run worker/path table | asset for the driver; the table is per-run |
| worker prompt files | the spec, path partition, traps, acceptance commands — **this is the intent record for the change** | `docs/_internal/<stage-stem>/reasoning/` (the memory skeleton) |
| run logs, byte counts, `-o` final messages | evidence, superseded by the next run | disposable; at most a line in `session.md` |

If that split is right, the whole intervention is: ship `probe.sh`/`launch.sh`
as skill assets, change three `/tmp` lines in the skill to name the in-tree
route for the prompt, and leave the logs where they are. That is a skill edit
plus two files — not a mechanism.

Attack it. Specifically: is a worker prompt really *memory*, or is it a
disposable input that merely looks valuable? Does putting prompts in-tree create
a landfill of the kind already measured (one project's `docs/_internal/` is
568 MB)? Does a prompt referencing absolute container paths and a specific
model version have any value a week later, or is the durable residue only the
*decisions* it encoded?

## 4. Constraints — these are settled, do not reopen

- **No new hook.** Not Claude, not Codex, not git. Enforcement ranking already
  decided: populated filesystem > `scio lint` > git hook > harness hooks.
- **Never scaffold empty.** A directory arrives with its first real file.
  `.gitkeep` was ignored in every project that shipped it.
- **Do not name a location nothing creates.** The `_scratch/` claim is being
  deleted precisely because no project has one. If you propose a path, say what
  makes it exist.
- **No hand-maintained metadata.** Rejected by the owner with reasons.
- **The seam:** follows the repo → scio; follows the human → dev-env. dev-env
  may name a path for a tool, never for an agent.
- **Simplicity is the binding constraint.** The owner: *"I really don't want
  this make more complex."* Removing a convention counts as progress.

## 5. Questions

1. Verify §2(a), (b), (c). Correct anything wrong.
2. Take a position on §3's split, per artifact class.
3. **Whose job is normalization** — scio (catalog skill + lint), dev-env
   (user-global habit), or neither? Note that "delegation replaces the typing,
   never the verification; re-run the tests yourself, confirm which files the
   agent touched, confirm no git write reached any repo it could see — a
   submodule is its own repo" is a *stance*, not a flag reference. Where does a
   stance live under the seam?
4. **Is any of it enforceable?** Say plainly what cannot be, rather than
   proposing a check that would need semantic judgment.
5. **The duplication.** Umbrella `AGENTS.md` and `delegate-cli` now both carry
   container-sandbox findings. One claim, one home — which home, and what
   happens to the other copy? Consider that the umbrella's copy exists because
   the catalog was unreachable, and that the fleet sweep will make it reachable.
6. **First move**, doable in one bounded edit, that improves this without
   depending on the rest.

Ground every claim in what you read. Distinguish verified from inferred.
