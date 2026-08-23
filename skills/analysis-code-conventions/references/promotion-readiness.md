# Promotion Readiness

Use this guide when a cohesive helper family may have outgrown one analysis
project.

## Telos

A helper family that becomes load-bearing, important, or reused across projects
earns its own repository with a cleanly designed user API. 

Promotion gives
- capability an explicit audience, 
- ownership model, 
- test surface, 
- and release life.

## Signals that justify the move

Reinforcing signals:
- two or more real projects need the capability without project-specific edits;
- callers rely on a stable conceptual API rather than copying internals;
- failures would block important analyses or corrupt consequential outputs;
- independent tests can express the behavior outside any one project fixture;
- dependencies and supported environments can be stated clearly;
- someone owns review, releases, compatibility, and issue response; and
- the capability needs a cadence different from the originating analysis.

Importance alone can justify promotion when the family is load-bearing and
independent stewardship materially improves reliability. Reuse alone is weaker
when every caller still requires a different interface.

## Signals to keep it in the project

Keep the family local while:
- its vocabulary and schema remain specific to one dataset or claim;
- the API changes whenever the analysis question changes;
- tests require most of the originating project to construct a fixture;
- there is one caller and little expected reuse;
- dependencies are inseparable from the project's environment; or
- nobody can own releases and compatibility.
Local code is a valid destination. 

## Readiness review

Before moving code, write down:
1. the intended users and their smallest useful entry points and workflow behaviour;
2. guarantees, refusals, schemas, and error behavior;
3. the project assumptions that must become arguments or adapters;
4. the minimum independent test matrix;
5. ownership and release expectations; and
6. the migration path for existing callers.

Design the user API from those needs. 
The local helper façade is evidence, not an API that must be exported unchanged.

## Decision

Suggest to promote when the capability has a coherent audience and independent lifecycle.
Keep it local when extraction would mainly relocate project assumptions or create
maintenance ceremony without a durable user contract. 
Record the decision where the project's architecture decisions live and revisit it when the signals change.
