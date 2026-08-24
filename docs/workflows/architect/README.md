# Architecture-first workflow documentation

An architecture-first software workflow with composable expert review, staged
artifacts, and human judgment gates.

---

## TL;DR

The mounted entry point is
[`skills/architecture-first-dev/SKILL.md`](../../../skills/architecture-first-dev/SKILL.md).
It selects one route and loads that route's complete specification from
`skills/architecture-first-dev/references/`.

```text
map → review → synthesize → design → plan → implement → verify
```

`review` and `synthesize` are conditional. `architect`, `diagram`, and `status`
are utility routes. The portfolio routes coordinate two or more related
features, and the architecture-treemap route supports retrospective audits.

Each route produces or inspects the durable Markdown artifacts described in
[00-quickstart.md](./00-quickstart.md). Implementation and large-campaign
decomposition remain external commands at
[`commands/implement.md`](../../../commands/implement.md) and
[`commands/decompose.md`](../../../commands/decompose.md).

## Document index

| File | View | What it answers |
|---|---|---|
| [00-quickstart.md](./00-quickstart.md) | Start here | Cadence, route selection, stage transitions, and worked examples |
| [01-architecture.md](./01-architecture.md) | Structural | Router, references, agents, commands, mounts, and artifact flow |
| [02-behavior.md](./02-behavior.md) | Behavioral | Exact stage behavior and data flow |
| [03-decisions.md](./03-decisions.md) | Decisions | Why the workflow has this shape and where its gates live |
| [04-extending.md](./04-extending.md) | How-to | Adding a reviewer, stage reference, or router route |
| [05-meta-approach.md](./05-meta-approach.md) | How-to | Cross-feature map, design, apply, and planning |

## Availability

Bind the toolkit catalog into a project from that project's pinned checkout:

```bash
./01_modules/scio/bin/scio link
```

`link` mounts the complete `skills/`, `agents/`, and `commands/` trees into
both `.claude/` and `.agents/`. The architecture-first router, its references,
the reviewer agents, and the external implementation commands then travel as
one pinned catalog.

Invoke the `architecture-first-dev` skill explicitly or describe the
architecture task in terms covered by its frontmatter. The router preserves
the selected route's arguments and loads its full reference before acting.

## Common routes

```text
map <feature>                         → docs/{feature}/map.md
review <feature> --as <spec>          → docs/{feature}/review/<reviewer>.md
synthesize <feature>                  → docs/{feature}/synthesis.md
design <feature>                      → docs/{feature}/design/*.md
plan <feature>                        → docs/{feature}/plan/phase-NN.md
implement <feature> <phase>           → code through commands/implement.md
verify <feature> [phases]             → docs/{feature}/verify.md
status [feature,...]                  → chat-only portfolio snapshot
```

The optional portfolio routes are `meta-map`, `meta-design`, `meta-apply`, and
`meta-plan`; one shared reference specifies all four.

## Review invocation shapes

```text
review <feature> --as all
review <feature> --as all --but stat,divergent
review <feature> --as bioinf,ml
review <feature>
```

The last form is interactive. An iterate round adds `--iterate` and preserves
the preceding reviews under `.history/`.

## Scope

This workflow is intended for non-trivial features, refactors, redesigns,
cross-cutting portfolio decisions, and retrospective architecture audits.
One-line fixes, throwaway scripts, and routine data analysis usually proceed
directly through the relevant coding or scientific skill.

## Shipped pieces

| Artifact | Location | Purpose |
|---|---|---|
| Router | `skills/architecture-first-dev/SKILL.md` | Selects a route and loads its complete contract |
| Stage specifications | `skills/architecture-first-dev/references/*.md` | Own phases, arguments, gates, outputs, and stop conditions |
| Reviewer and pipeline agents | `agents/*.md` | Provide cartography, expert lenses, synthesis, gates, and portfolio operations |
| External commands | `commands/implement.md`, `commands/decompose.md` | Preserve implementation and campaign-decomposition boundaries |
| Workflow rationale | `docs/workflows/architect/` | Canonical sequence, behavior, decisions, and extension guide |
| Runtime artifacts | `docs/{feature}/`, `docs/_meta/` | Durable feature and portfolio state in consumer projects |

See [01-architecture.md](./01-architecture.md) for the component map.
