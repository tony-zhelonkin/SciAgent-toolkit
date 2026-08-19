# Helper-Family APIs

Use this guide to choose between a focused flat helper and a cohesive helper
family, then design the boundary a stage calls.

## Why the named boundary matters

A stage reads in story order because implementation chunks sit behind a named
entry point. The reviewer can mentally outsource “here we do this, then that,”
continue through the analytical narrative, and open the machinery only when a
specific implementation question makes that necessary.

## Flat helper or family

A single focused helper stays a flat file. A subfolder is earned when several
functions or files form a cohesive family by sharing one or more of:

- vocabulary that has a stable domain meaning;
- invariants enforced across multiple operations;
- state or resources with a common lifecycle;
- tests that describe one behavioral unit; and
- one small entry point that can represent the family to its callers.

File count alone supplies no evidence of cohesion. Two unrelated helpers do not
become a family by sharing a directory, and a long focused implementation does
not automatically need a folder.

## Design the provisional API

Start from the call site the stage should read, then work inward:

1. Name the analytical operation in project vocabulary.
2. Pass consequential inputs and configuration explicitly.
3. Return the smallest result that expresses the stage handoff.
4. State guarantees and refusals at the entry point.
5. Keep orchestration, adapters, and low-level transformations internal.

“Provisional” means the boundary may evolve with the project. It should still be
deliberate: callers depend on its behavior even before it belongs to a package.

## Public and internal surfaces

The entry point exposes domain inputs, consequential choices, a named result
contract, and failures the stage can act on. Internal files own adapters,
tool-specific calls, shared state, intermediate representations, batching, and
the invariant checks used across operations.

Avoid a façade that simply mirrors every internal function. A small API earns its
value by reducing how much machinery the stage must understand.

## Configuration boundary

Keep project values in the project's configuration surface, validate them near
loading, and pass relevant values into the helper API. Config loaders may derive
paths and provide typed accessors. Hidden mutable configuration obscures tests
and restart behavior.

Defaults belong in code only when they are stable implementation defaults.
Scientific thresholds, dataset identities, and choices a reviewer may revisit
belong in project configuration or an explicit decision record.

## Example shapes

One focused operation can remain direct:

```text
helpers/validate_contrasts.R
```

A cohesive capability can earn a family:

```text
helpers/contrast_model/
├── fit.R
├── schema.R
├── diagnostics.R
└── contrast_model.R    # small entry point sourced by stages
```

The façade filename is a project choice until observed practice supports a
stronger convention. Do not invent a universal `api.R` or `__init__.py` rule.

## Tests and callers

Test the public behavior through the entry point, with focused internal tests
for risky invariants. Search callers before changing the boundary. A family is
successful when stage code imports or sources the small surface and tests can
exercise the capability without reconstructing the stage.

## Promotion route

When the family becomes load-bearing, important, or reused across projects,
read [`promotion-readiness.md`](promotion-readiness.md). Cohesion earns a local
subfolder; broader responsibility and a stable user contract may earn a repository.
