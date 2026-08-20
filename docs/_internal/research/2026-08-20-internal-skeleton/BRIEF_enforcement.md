# Consultation — harness-agnostic enforcement of project structure

You are consulting on **mechanism and ownership**, not on the directory layout
(that is settled enough — see §2). Produce a written recommendation. Modify
nothing.

Read first: `FINDINGS_field.md` (what 24 projects actually contain) and
`CONSULT_skeleton.report.md` (the layout recommendation this builds on).

## 1. The question

The owner, in their words:

> "I'm thinking though how to enforce `_internal` and its structure in an
> agent-agnostic manner, meaning no matter codex or claude code all would keep
> their thoughts organised and actions intent-full. What kind of hooks, whose
> responsibility would they be, would encourage, enforce, and guardrail the
> agents into a clean tidy structure? I already mostly see them follow the
> `01_` `02_analysis` `03_results` and its decomposed folder structures with
> adjacent README files and tables and figures derived from the tables living
> beside each other, but though **mostly through example not through AGENTS.md
> encouragement or hooks** or anything. Should this be responsibility of scio or
> `/data1/users/antonz/pipeline/dev-env`?"

Three questions, in priority order:

1. **What is the enforcement mechanism**, given that it must work identically for
   Claude Code, Codex, and a harness that does not exist yet?
2. **Whose responsibility** — `scio` (project-scoped, versioned per project) or
   `dev-env` (user-global, follows the human across machines)? Or split, and
   where exactly is the seam?
3. **The example-versus-instruction puzzle** in §3, which is the most
   interesting thing in the owner's message and may invalidate the framing of
   question 1.

## 2. The layout, assumed settled

Memory under `docs/_internal/` mirrors the analysis structure, keyed by stage
stem, so intent is filed where the work happened:

```
02_analysis/stages/30_grn.R  →  03_results/30_grn/{figures,tables}/  →  docs/_internal/30_grn/
```

with `_project/` for memory that spans or precedes stages. Two durable forms
only: `reasoning/` (topic files) and `session.md` (updated in place; git history
is the archive). The owner also raised mirroring `02_analysis/helpers/` — assess
whether that helps or is a category error.

Do not redesign this. Assume it and answer how it gets enforced.

## 3. The puzzle at the centre — resolve this first

The owner observes that agents **already** comply with `01_`/`02_analysis`/
`03_results`, adjacent READMEs, and figure-beside-its-source-table — and that
this compliance came from **example, not instruction**. No hook enforces it. The
always-on text mentions it, but the owner's read is that the tree itself is what
teaches.

Now put that beside two facts from the field survey:

- **Empty scaffold is active harm.** Every STING child project ships
  `handoffs/` containing only `.gitkeep`. Agents wrote their handoffs to
  `sessions/` instead, in every case. An empty directory did not teach; it
  misdirected. Ten spellings of "handoff" exist across the fleet.
- **A populated tree taught correctly.** `03_results/<stage>/{figures,tables}/`
  with adjacent READMEs is followed nearly everywhere, and it is the part that
  always has real content in it.

So the hypothesis to test: **structure propagates by example only when the
example carries content.** An empty directory is not an example, it is an
unsupported claim, and agents correctly ignore it. If that is right, the
mechanism for `docs/_internal/` is not a hook and not a scaffold — it is
*ensuring the first real entry exists and is good*, then letting the pattern
replicate.

Take a position on this. If it is right, say what follows for mechanism and for
migration into the 10 projects that have no `docs/_internal/` at all. If it is
wrong, say why and what actually drives the observed compliance.

## 4. Enforcement surfaces — establish the facts, then choose

**Claude Code.** Rich hooks via `.claude/settings.json`. scio currently ships
exactly two, materialized by `link` under a hash-and-cede ownership discipline:
`no_ephemeral.sh` (PreToolUse, matcher `Bash`) and `caption_sweep.sh` (Stop).
So the PreToolUse-veto and end-of-turn-sweep pattern is already proven here.

**Codex.** Less clear, and you should determine it rather than assume.
`codex-cli 0.147.0` exposes `--dangerously-bypass-hook-trust` ("run enabled
hooks without requiring persisted hook trust"), so a hook system with a trust
model exists. `~/.codex/` contains `skills/`, `rules/`, `policy/`, `plugins/`,
`memories/`, `sessions/`. Inspected: `rules/default.rules` and
`policy/default.codexpolicy` are command-approval allowlists
(`prefix_rule(pattern=[...], decision="allow")`), not structural guardrails.
`memories/` is empty. **Find the actual hook configuration format and event
surface** — read the binary's help output, any bundled docs, config schema
errors, `~/.codex/plugins/`, whatever resolves it. Report what you could not
determine.

Critical field fact: **`.codex/` exists in zero of 24 projects.** Codex agents
leave no project-level footprint today, while Claude leaves settings, locks,
hooks and status lines. Whatever is built must not require per-harness project
directories to multiply.

**Git hooks.** The only enforcement point that fires regardless of which agent
acted, because every harness eventually runs `git`. Consider `pre-commit` /
`commit-msg` seriously, including their weaknesses: they are not installed by
clone, they are trivially bypassed with `--no-verify`, they fire only at commit
time, and they cannot see work that is never committed — which is exactly the
`/tmp` scratch case.

**The filesystem.** A path that exists and contains a good example. No
execution, no configuration, no harness awareness. Weakest enforcement,
strongest reach. §3 argues this may be the real answer.

**`scio lint`.** Already the house answer for "a convention we can measure."
Harness-blind by construction — it reads the tree. But it only runs when
somebody runs it.

Rank these by (reach × strength ÷ complexity) for this specific job, and say
which combination you recommend. Name what each surface *cannot* enforce.

## 5. Ownership — scio or dev-env

Facts:

- **`scio`** is vendored per project as a git submodule, versioned with the
  project, and owns "project tree, analysis config, docs namespaces, AI harness
  (agents/skills/commands)". Three verbs: `link` (availability, convergent),
  `craft` (always-on standing text, idempotent, ≤25 rendered lines with 8
  spare), `lint` (enforcement, read-only). Bash only, no `jq`/`yq`/`python3`.
- **`dev-env`** is `nvim` + `tmux` + dotfiles + a notebook layer (inline figures,
  Jupyter kernels). Rolled out per machine or mounted read-only into a container
  at `/opt/dev-env`. It owns **the human's tools**, user-global, following the
  person rather than the project. It is clean at `0b5a4b0`.
- An 8-file bundle already sits staged at
  `docs/_internal/handoff/dev-env/` — user settings, global `AGENTS.md`
  template, statusline, provisioning scripts — intended as the *user-global
  harness* layer, adopted nowhere yet.
- The house ownership test, from `scbio-docker/AGENTS.md`: *does the content
  change when the **image** changes → scbio-docker; when the **project / AI
  harness** changes → SciAgent-toolkit.*
- **ADR-D5 is decided but unimplemented**: extract project binding behind a
  *common* layer plus *harness adapters*. That is precisely the
  harness-agnostic question, already decided in principle. Say whether this work
  should be the thing that finally forces D5, or whether D5 is a distraction
  here.

Apply the ownership test. If the answer is "both", give the seam precisely: what
is user-global (travels with the human, any project) versus project-scoped
(travels with the repo, any human). Consider explicitly: an agent's *habits*
arguably belong to the user; a project's *structure* clearly belongs to the
project. Where does "write your reasoning to `docs/_internal/<stage>/`" fall?

## 6. Constraints

- **`I really don't want this make more complex.`** Removing a convention counts
  as progress. Every added moving part needs justification.
- The owner has **rejected** intent metadata that humans or agents must maintain
  by hand — specifically config-file knobs (`decisions.<stage>`,
  `status: APPROVED`). Do not re-propose a scheme with a field to fill.
- 10 of 24 projects have no `docs/_internal/`. Absence is the majority state.
- Always-on text is scarce: 8 spare lines, and the survey shows deletions are
  available (the root `_scratch/` claim names a directory that exists nowhere).
- The owner is not interested in remote/hub syncing right now. Durability of
  *structure and intent* is the subject; where bytes get pushed is not.

## 7. Deliverable

1. **Your position on §3**, with what follows from it.
2. **The mechanism** — concretely, which surfaces, which events, what each one
   does. If you recommend a hook, give its trigger, its verdict semantics
   (advisory versus blocking), and what it inspects. If you recommend no hook,
   say so plainly and defend it.
3. **Harness parity** — how Codex gets the same behaviour as Claude Code, given
   `.codex/` exists nowhere today. Include what you determined about Codex's
   hook surface and what you could not.
4. **The ownership verdict** — scio, dev-env, or a named seam. One sentence
   stating the test you applied.
5. **What this makes unnecessary.** Existing hooks, always-on lines, or
   conventions that this design lets you delete.
6. **The first move.** One concrete change, doable this week, that improves
   things and does not depend on the rest being right.

Ground everything in what you can verify by reading. Distinguish what you
checked from what you inferred.
