# Contributing

## Before you start

Read `AGENTS.md` — it covers the critical rules (idempotency, bash-only, manifest ownership, test gate).

## Running tests

```bash
bash tests/run-all.sh
```

All tests must pass (`bash tests/run-all.sh`). Tests cover activate/deactivate/inject round-trips, block hash drift detection, and manifest ownership.

## Name collisions

CI fails when a new skill, agent, or command basename collides with an existing entry in another namespace.

```
error: name 'clash-name' collides across skill and command. either rename one, or add to tests/collision-allowlist.txt with rationale.
```

Fix: rename the new entry, or — if the overlap is intentional (e.g., a command that invokes a same-named skill) — add a line to `tests/collision-allowlist.txt`. The file's header documents the line format and kind-ordering rule.

If an already-allowlisted name later acquires a new namespace, CI also fails:

```
error: name 'clash-name' collides across skill,agent,command but allowlist records 'skill,command'. update tests/collision-allowlist.txt to reflect the new kinds, or rename one.
```

Fix: update the allowlist row to include the new kind.

Manifest schema (for the curious): `docs/architecture.md` §9.

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
