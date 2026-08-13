# Contributing

## Before you start

Read `AGENTS.md` — it covers the critical rules (idempotency, bash-only, link ownership, test gate).

## Running tests

```bash
bash tests/run-all.sh
```

All tests must pass (`bash tests/run-all.sh`). Tests cover link convergence and ownership, block hash drift detection, linting, and catalog validation. Run `bin/scio lint --check toolkit --strict --quiet` for the catalog gate alone.

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

## Skill lifecycle

Skills move from active use to `_attic` by judgment. The directory location is
the lifecycle state. See `docs/skill-lifecycle.md`.

## Error handling (lib/scio)

Library functions in `lib/scio/*.sh` only ever `return <code>` — never `exit` (only `bin/scio`, at the top dispatch level, may exit; awk/subshell `exit` is fine since it tears down the awk/subshell, not the caller's shell). Every side-effecting call (`block_write`, `ln -sfn`, `mkdir -p`, …) is checked: `cmd || { echo "scio <verb>: <message>" >&2; return 1; }`. User-facing errors use the prefix `scio <verb>: <message>` on stderr. Reserve `|| true` for genuinely best-effort, non-state operations and annotate each with a `# best-effort: <reason>` comment. `set -e`/`set -o pipefail` are deliberately off — the codebase relies on explicit return-code dispatch.

## Commit messages

Imperative mood, subject ≤72 chars. No AI attribution trailers.

```
# good
add divergent-thinking agent
preserve private directories during link

# bad
Added new agent
Update
```

## Pull requests

- One logical change per PR
- Include the `tests/run-all.sh` output in the PR description
- Update `docs/agents.md` or `docs/skills.md` if you add content to those directories
