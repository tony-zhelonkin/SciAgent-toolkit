# Contributing

## Before you start

Read `AGENTS.md` — it covers the critical rules (idempotency, bash-only, manifest ownership, test gate).

## Running tests

```bash
bash tests/run-all.sh
```

All 13 tests must pass. Tests cover activate/deactivate/inject round-trips, block hash drift detection, and manifest ownership.

## Commit messages

Imperative mood, subject ≤72 chars. No AI attribution trailers.

```
# good
add divergent-thinking agent to base role
fix inject when overlay already holds the same skill name

# bad
Added new agent
Update
```

## Pull requests

- One logical change per PR
- Include the `tests/run-all.sh` output in the PR description
- Update `agents/README.md` or `skills/README.md` if you add content to those directories
