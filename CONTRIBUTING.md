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

## Skill lifecycle

Skills move through a natural `experimental → stable → deprecated → _attic`
lifecycle, driven by judgment (not hooks). The optional `metadata.status:` field
(default `stable`) is surfaced softly in `sciagent list skills`; retired skills
move to `_attic/` as reference-only. See `docs/skill-lifecycle.md`.

## Error handling (lib/sciagent)

Library functions in `lib/sciagent/*.sh` only ever `return <code>` — never `exit` (only `bin/sciagent`, at the top dispatch level, may exit; awk/subshell `exit` is fine since it tears down the awk/subshell, not the caller's shell). Every side-effecting call (`block_write`, `manifest_*`, `ln -sfn`, `mkdir -p`, …) is checked: `cmd || { echo "sciagent <verb>: <message>" >&2; return 1; }`. User-facing errors use the prefix `sciagent <verb>: <message>` on stderr. Reserve `|| true` for genuinely best-effort, non-state operations and annotate each with a `# best-effort: <reason>` comment. `set -e`/`set -o pipefail` are deliberately off — the codebase relies on explicit return-code dispatch.

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
- Update `docs/agents.md` or `docs/skills.md` if you add content to those directories
