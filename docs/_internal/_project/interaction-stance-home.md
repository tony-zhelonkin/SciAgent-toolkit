# Where the interaction stance lives — it is not scio's to hold

**Date:** 2026-08-24 · **Closes:** task #33

## Scope

The stance retired with `system-prompts/architect-mentor.md`: the user drives the
thought process, the assistant extends their reach and asks for the judgement only
the user can supply. It has been active nowhere since — `mentor-mode` is
invoked-only, and `docs/_internal/handoff/user-level-background.md` is staged and
unadopted.

## Decision

**The stance is user-level context and belongs in `dev-env`, which is what the
staged file already declares as its destination.** scio should stop looking like a
candidate home for it.

## Evidence

The question was framed as a choice between a `craft.yaml` bullet, the dev-env
user-level file, the conventions router skill, or a combination. Scope settles it
without appeal to taste:

- **`craft.yaml` is wrong by scope.** That block is per-analysis-project and
  concerns code, results, figures and memory — properties of the work. The stance
  concerns how an assistant collaborates with one person, and it holds while
  reading a paper or debugging a container, not only inside an analysis repo.
  Putting it there would propagate a personal preference into ~25 project
  `AGENTS.md` files where it describes nothing about the project.
- **The router skill is wrong by reachability.** `analysis-code-conventions`
  routes judgement calls about code. A stance nobody thinks to ask for is a stance
  that never applies, which is exactly `mentor-mode`'s current condition.
- **dev-env is right by reach.** Harness-neutral user-level context loads every
  session, everywhere, for every repository — which matches what the stance
  governs.

## The one part that IS scio's, and is already done

The staged file's **Language caution** — abstract "X, not Y" constructions obscure
the intended operation, so prefer concrete actors, actions, failure cases and
remedies — is a property of prose this repository emits. It is already active in
two places that bind: `craft.yaml`'s always-on Prose convention, and AGENTS.md
rule 6 for agents editing the toolkit itself. No further action.

## What remains, and whose it is

Adopting `user-level-background.md` into `/data1/users/antonz/pipeline/dev-env`
is an **owner action** in another repository. Until then the file sits in scio's
memory under `handoff/dev-env/`, which is part of the rehoming task #52 lists:
those are assets staged for elsewhere, not memory.

Nothing in scio changes for this decision. That is the finding.
