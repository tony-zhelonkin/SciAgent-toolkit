# Consultation — redesign `delegate-cli` into durable, normalised behaviour

Consulting on **the design of one skill and the assets it ships**. Produce a
written recommendation. Modify nothing.

Read: `skills/delegate-cli/SKILL.md` in full (416 lines), the three
`investigate-*.report.md` in this directory, and
`../../plans/2026-08-20-memory-and-seams/07_delegation-seam.md` (which already
commits to shipping `probe.sh`/`launch.sh` as skill assets and to a retention
contract).

## 1. Five verified failures of the current skill

**(a) It prescribes a flag that does not exist.** `SKILL.md:58` and `:246` tell
the reader to pass `--search` for web search. Verified absent from
`codex exec --help` in **both** 0.147.0 (host) and 0.149.0 (JR-MC container).
Web search is config-gated: `-c tools.web_search=true --enable
web_search_request`. An operator burned several rounds running `strings` on the
codex binary to discover this. The skill did not merely fail to help — **it
caused the failure.**

**(b) It documents flags against a moving target.** Its devcontainer section is
dated "Verified 2026-07-28". Two codex versions are live in the fleet right now
(0.147.0 host, 0.149.0 container) and the container's is installed by an
*unpinned* `setup_ai_env.sh`. A prose skill asserting CLI flags is a claim whose
only enforcement is the installed binary's version — which scio does not
control.

**(c) Its backgrounding advice does not apply to how the owner runs it.**
`SKILL.md:168` says *"Run in background if it's real implementation
(`run_in_background: true`)"* — a **Claude Code tool parameter**. The owner runs
delegations himself with a `!` prefix, where no such parameter exists. Observed:
`! .devcontainer/scripts/codex_unit.sh u1_build` blocked the session for 1m42s
with no way to detach. Meta-Aging's inline form gets backgrounded by accident
(Claude's 120s timeout moves it) or by a manual `&`.

**(d) It says nothing about concurrency, and the gap corrupted a run.** Two
`u1_build` runs overlapped for ~2.5 minutes writing the same paths. The
orchestrator's own conclusion: *"I'll treat its self-report as unreliable and
verify the files against disk."* The duplicate was killed (exit 143). No lock,
no already-running detection, no per-unit output uniqueness. Note the owner's
framing: parallel workers are **desirable** for disjoint units — *"not
inherently bad and good for some workflows but again would need some
regulation."* So the answer is not "forbid concurrency" but "make it explicit
and safe."

**(e) Long inline invocations wrap and break.** A wrapped line split `-m
gpt-5.6-sol`, so codex reported *"a value is required for `--model <MODEL>`"*
and bash then tried to execute `gpt-5.6-sol`. Measured: JR-MC's wrapper reduces
the typed command from a 373-441 char expanded pipeline to **44-55 typed
characters**, with a 9.7-10.5 KB prompt through stdin.

## 2. The asymmetry to explain

Meta-Aging works. JR-MC fought. Same toolkit, same image, same operator.

The investigation found Meta-Aging holds because the operator **consistently**
uses external prompt files through stdin (4.9-7.7 KB each), and because its
tasks forbade network access so `--search` never came up. But: *"The exact
executed command was not preserved in a script, tracked document, or run log."*

So Meta-Aging's success rests on **operator habit, not mechanism**. JR-MC
reached for an inline command first and hit every sharp edge. That is the thing
to fix: a habit that works is one keystroke from not working.

## 3. Three field conventions, none of them scio's

| Project | Prompts | Outputs | Invocation |
|---|---|---|---|
| JR-MC | `docs/_internal/codex/{00_house_rules.md,<unit>.md}` | `logs/codex/<unit>.{final.md,run.log}` | `codex_unit.sh <unit> [--web]` |
| Meta-Aging | `03_results/_scratch/codex_handoff/<UNIT>_prompt.md` | `<UNIT>_run_last.md` beside it | inline `codex exec … - < prompt.md` |
| 14782-DM | `/tmp/…/scratchpad/codex/` | same | `probe.sh` + `launch.sh` there |

`delegate-cli` prescribes a **fourth**: `/tmp/<unit>_prompt.md`.

**And JR-MC's fix cannot survive.** `codex_unit.sh` is untracked;
`docs/_internal/codex/` is untracked *and* gitignored by the toolkit's own
`SCIO:GITIGNORE`. A fresh clone reconstructs none of it. The best answer in the
fleet is the least durable.

## 4. Questions

1. **Capability drift (a, b).** How should a skill relate to a CLI that
   version-drifts underneath it? Options: probe at run time and report; pin the
   CLI version in the container; ship a capability table keyed by version;
   assert nothing and read `--help`. Recommend one, and say what the skill
   should stop claiming. Consider that "verified on <date>" prose has already
   failed once.
2. **The shipped asset's contract.** Phase 07 commits to `probe.sh` and
   `launch.sh` as `skills/delegate-cli/assets/`. Specify them concretely.
   `launch.sh` must at minimum resolve: backgrounding (c), a per-unit lock with
   explicit opt-in for parallel disjoint units (d), short invocation (e), unique
   output paths, exit-status and byte-count reporting, failure when an expected
   output is absent. Say what it must **not** do.
3. **Where does project-local delegation state live?** House rules and unit
   prompts are project-specific and currently land in a gitignored tree.
   Reconcile with decisions already taken: memory is stage-keyed under
   `docs/_internal/<stage-stem>/`; delegation *prompts are disposable* while
   their decisions are durable; nothing may be scaffolded empty; and no location
   may be named that nothing creates. Does a house-rules file belong in the
   *project* at all, or is it the skill's own content?
4. **Normalising four conventions into one.** Which survives, and what makes the
   other three converge on it? The fleet sweep is the only delivery event
   available. Be concrete about what an existing project has to do.
5. **What is enforceable.** Which of these can `scio lint` check, which can the
   asset enforce mechanically, and which are irreducibly the operator's
   judgement? Do not propose a check needing semantic judgement.
6. **The verification stance.** *"Delegation replaces the typing, never the
   verification. Re-run the tests yourself, confirm which files the agent
   touched, confirm no git write reached any repo it could see — a submodule is
   its own repo."* Where does that live: skill prose, asset behaviour, or the
   user-global habit layer? Note the seam already decided: follows the repo →
   scio; follows the human → dev-env.

## 5. Constraints — settled, do not reopen

- **No new hook.** Not Claude, not Codex, not git.
- **Never scaffold empty**; never name a location nothing creates.
- **No hand-maintained metadata.**
- **Bash only** in anything scio ships — no `jq`, `yq`, `python3`.
- Skill descriptions cap at 350 chars; `lint --check toolkit` enforces it.
- **Source shape must equal mount shape** — a directory symlink cannot filter,
  so anything added under `skills/delegate-cli/` is mounted into every project.
  `assets/` is therefore reachable with no new delivery mechanism, and also
  unavoidable.
- Simplicity binds: *"I really don't want this make more complex."* Deleting
  counts as progress. The current skill is 416 lines; say what should be cut.

## 6. Deliverable

1. The capability-drift answer (Q1), and the specific lines/claims to delete.
2. `probe.sh` and `launch.sh` contracts (Q2) — behaviour, arguments, exit codes,
   what they refuse. Pseudocode or a shell sketch is welcome; it must be
   POSIX-ish bash with no external deps.
3. Where project-local delegation state lives (Q3), consistent with §5.
4. The convergence plan (Q4).
5. The enforcement map (Q5) — asset / lint / operator-judgement, and what is
   *not* enforceable.
6. Ownership of the stance (Q6).
7. What to delete from the 416 lines.
8. Where you disagree with this brief.

Ground everything in what you read. Distinguish verified from inferred.
