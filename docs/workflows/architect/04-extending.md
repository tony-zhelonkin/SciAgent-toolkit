---
parent: ./README.md
view: how-to
---

# Extending the architecture-first workflow

The workflow is a router plus complete stage references. Extend that structure
at the layer that owns the new behavior.

---

## Choose the extension type

| Need | Add or edit |
|---|---|
| Another lens within `review` | `agents/<name>.md` and the roster in `references/review.md` |
| Another workflow stage | `references/<stage>.md` and a router row in `SKILL.md` |
| Another portfolio operation | The route table and `references/portfolio.md`, or a new reference if the contract is independently large |
| Another retrospective audit operation | `references/architecture-treemap-audit.md`, or a separate reference when it has an independent contract |
| Independently invocable execution workflow | `commands/<name>.md` plus a router boundary row when architecture-first work routes to it |

Every project linked to the toolkit sees the full `skills/`, `agents/`, and
`commands/` trees. Catalog availability therefore follows the source tree;
selective role manifests are absent from this model.

## Add a reviewer lens

Example: add a `security` reviewer.

1. Copy a reviewer agent and edit its frontmatter and body:

   ```bash
   cp agents/stat.md agents/security.md
   ```

2. Give the agent a unique `name`, a usage-focused `description`, its model and
   tools, a precise lens, the `map.md` input contract, and the single output
   path `docs/{feature}/review/security.md`.

3. Update
   [`references/review.md`](../../../skills/architecture-first-dev/references/review.md):

   - add `security` to the canonical roster;
   - define when `--as security` is useful;
   - include it in `--as all` only if every full-panel review should pay for
     that lens;
   - preserve the round-1 and iterate-round contracts.

4. Update the reviewer tables in this documentation and add focused tests when
   the new roster creates behavior that can be checked mechanically.

5. Run the toolkit catalog check and full suite:

   ```bash
   ./bin/scio lint --check toolkit
   bash tests/run-all.sh
   ```

An already-correct whole-tree link exposes the new agent immediately. A newly
bound project receives it through `scio link`.

## Add a workflow stage

Example: add an `audit` stage that writes `docs/{feature}/audit.md`.

### 1. Define the contract

Write down:

- the user intent that selects the stage;
- positional arguments, flags, defaults, rejections, and usage text;
- required upstream artifacts;
- tools and sub-agents;
- the exact output path and format;
- human gates, stop conditions, and next-route guidance.

### 2. Create the complete stage reference

Add:

```text
skills/architecture-first-dev/references/audit.md
```

Use a neighboring reference as the structural template. Keep the full phase
logic in this file so the router can load one authoritative contract after it
selects the stage.

### 3. Add the router row

Edit `skills/architecture-first-dev/SKILL.md` in both places that define
reachability:

- add the intent branch to the routing decision tree;
- add a route-table row linking `references/audit.md` and naming its argument
  entry.

If the stage leads to or consumes an external command, add that boundary to the
router contract as well.

### 4. Add an agent only when the stage needs one

Create `agents/auditor.md` when a distinct persona, model, or tool boundary is
part of the stage design. A main-agent stage can stay entirely inside the
reference.

### 5. Keep the canonical docs aligned

Update:

- `00-quickstart.md` for cadence and after-stage tips;
- `01-architecture.md` for component and artifact flow;
- `02-behavior.md` for exact stage behavior;
- `03-decisions.md` when the new shape carries a durable trade-off.

### 6. Verify

Run `./bin/scio lint --check toolkit` and `bash tests/run-all.sh`. Exercise
the route against a scratch project when its file or dispatch behavior is new.

## Extend the portfolio routes

`references/portfolio.md` owns `meta-map`, `meta-design`, `meta-apply`, and
`meta-plan` as one connected contract. Add a portfolio operation there when it
shares their artifacts and invariants. Split it into a new reference when its
inputs, phases, or safety boundary stand on their own; then add a router row for
the new reference.

Portfolio operations preserve these boundaries:

- `meta-architect` reads project artifacts and writes under `docs/_meta/`;
- `feature-reviser` applies accepted MADRs to one existing feature design;
- judgment-heavy redesign returns to the per-feature `design` route;
- compound edits retain their human gate.

## Add an external command

External commands are appropriate for execution workflows that remain useful
outside the architecture router. `implement` and `decompose` are the shipped
examples.

1. Add `commands/<name>.md` with its complete argument and execution contract.
2. Add a router boundary row when architecture-first work should hand off to
   it.
3. Link to the command with a repository-relative path from `SKILL.md`.
4. Update `docs/commands.md` and the relevant workflow docs.
5. Add tests for argument parsing, mutations, refusals, and idempotency as the
   command requires.

The complete `commands/` directory is already mounted by `scio link`.

## Domain specialization

For project-specific emphasis, pass the domain context through the feature
description and let `map.md` carry it downstream. A durable catalog-wide lens
belongs in a new reviewer agent and the review roster. Editing an existing
reviewer's semantics changes every linked consumer and should be treated as a
catalog change.

## Provider portability

The Markdown artifacts and stage contracts are portable. Harnesses differ in
agent frontmatter, explicit skill invocation, and parallel-dispatch mechanics.
Keep provider-specific mechanics at those boundaries while preserving route
arguments, artifact paths, gates, and stop conditions.

`scio link` exposes the same catalog through `.claude/` and `.agents/`.
Provider support beyond those discovery paths belongs in that provider's own
adapter or extension.

## Existing project artifacts

The workflow writes its canonical files alongside earlier documentation. A
new `map` or `review` route can therefore begin without moving a historical
feature folder. When consolidating older files, use `git mv` and classify each
file by its actual role: reviewer output belongs under `review/`, consensus at
`synthesis.md`, and accepted architecture under `design/`.

## Debug checklist

| Symptom | Check |
|---|---|
| Router does not select the new stage | The intent branch and route-table row both exist in `SKILL.md` |
| Router selects the stage but cannot load it | The linked `references/<stage>.md` path and filename match exactly |
| A reviewer name is rejected | The agent file and `references/review.md` roster agree |
| A consumer cannot see new catalog content | Its category path is a whole-tree link to the intended pinned toolkit |
| A route writes an unexpected file | The reference's output contract and the agent body name the same path |
| Portfolio work crosses feature boundaries incorrectly | Re-check `portfolio.md` and the selected agent's write scope |

The route reference is the executable specification; these workflow docs are
its navigation and rationale layer.
