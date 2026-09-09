# Scio toolkit architecture

## Purpose

Scio supplies project-local computational-biology context and guardrails
to AI coding harnesses. Git is the version authority. A consumer normally pins
the toolkit as a submodule and binds that exact checkout into its harness
discovery paths.

## Source layout

```text
scio/
├── bin/scio
├── lib/scio/          # common.sh first, then the modules one verb needs
├── skills/<name>/SKILL.md
├── agents/<name>.md
├── commands/<name>.md
├── craft.yaml
├── templates/
├── tests/
├── docs/              # this tree; docs/_internal/ is its own repository
└── _attic/            # retired skills, each still a SKILL.md
```

The three mount directories contain mountable content only. Their source shape
is identical to the consumer shape: skill directories and flat Markdown files
for agents and commands.

## CLI

The CLI has three verbs:

- `link` binds the catalog, materializes hooks, merges settings, and refreshes
  the managed gitignore block.
- `craft` renders the managed CRAFT block.
- `lint` runs project checks and the explicit toolkit catalog check.

The dispatcher resolves its own real path, derives `$SCIO_TOOLKIT`, sources
`common.sh` plus the dependency closure for one verb, and calls that verb's
`cmd_*` entrypoint. Module order is decided there rather than by modules
sourcing each other.

## Link topology

`scio link [--project-dir D]` creates six directory symlinks:

```text
D/.claude/skills   -> $SCIO_TOOLKIT/skills
D/.claude/agents   -> $SCIO_TOOLKIT/agents
D/.claude/commands -> $SCIO_TOOLKIT/commands
D/.agents/skills   -> $SCIO_TOOLKIT/skills
D/.agents/agents   -> $SCIO_TOOLKIT/agents
D/.agents/commands -> $SCIO_TOOLKIT/commands
```

The operation is convergent. It sweeps legacy toolkit-owned child links,
preserves user-owned links, and silently keeps correct bindings. A locality
check protects submodule pins.

A category whose directory holds anything the toolkit does not own — a skill
the project wrote, or one a third-party installer put there — keeps that
directory and receives one link per catalog entry beside the project's own.
The category symlink is the better binding because it needs no refresh when
the catalog changes, but it makes the mount point resolve into the toolkit
checkout, so an installer writing to `.claude/skills/` writes into the vendored
submodule, where the result is untracked and the next checkout there deletes
it. A project entry named like a catalog entry keeps loading, and `link`
reports the shadow on every run. Remove the last project entry and the
category converges back to the single symlink.

Those six paths are the door. The same files are also reachable through the
vendored tree, and that path is the fallback: a skill reached through a mount
discloses progressively, while one named by vendor path arrives whole.

The vendored directory is named `scio`, and `SciAgent-toolkit` for copies
predating the ADR-D9 rename. Discovery, ownership, freshness and the two
executable templates read one list, `lib/scio/common.sh::_SCIO_TOOLKIT_DIRS`,
preferring the new name — so a project may sit at either during the migration.

Analysis projects also receive shared helper-library links and import shims.
`link` refreshes the `SCIO:GITIGNORE` block because those bindings and the
harness/state paths are project-local artifacts.

## Guardrail hooks

`link` materializes `no_ephemeral.sh` and `caption_sweep.sh`, then merges their
registrations into `.claude/settings.json`. The settings merge preserves every
unrelated project setting.

Hook bodies follow hash-and-cede ownership. Recognized toolkit bytes refresh;
user-owned bytes remain in place and receive a one-time warning.
`templates/PROVENANCE.sha1` carries historical hashes for release archives and
shallow clones.

## Managed blocks

`block.sh` owns marker-framed updates in shared text files. `craft` writes the
`SCIO:CRAFT` block in `AGENTS.md`; a hand-edited body requires `--force`.
Text outside the markers and file permissions are preserved.

Historical `SCIAGENT:CRAFT` blocks are drift-checked and rewritten in place
with the `SCIO` prefix. `SCIAGENT:ROLES` and `SCIO:ROLES` blocks remain
readable and removable through `block_read` and `block_remove`. The writer
accepts CRAFT alone.

## Lint

`lint` runs thirteen selectable checks: `figure-style`, `results-layout`,
`captions`, `provenance`, `freshness`, `hooks`, `docs-layout`, `stage-thinness`,
`comment-intent`, `stage-layout`, `harness-links`, `internal-memory` and
`toolkit`. Findings warn by default and become failures under `--strict`.

`internal-memory` and `toolkit` are **opt-in**: neither is a member of `all`,
because each audits a subject a project may legitimately not have yet. Every
check treats an absent subject as silence rather than a finding.

`scio lint --check toolkit` is the catalog gate the release builder uses. It
checks skill frontmatter shape and description length, reports basename
collisions across skills, agents and commands, holds the CRAFT body to
`craft.yaml`'s `max_lines`, requires every `_attic/` entry to be a retired skill,
refuses a `.gitkeep` or an empty directory under `skills/` and `templates/`, and
requires delegate-cli's three assets to be present, executable and parseable —
that skill's instructions are executed rather than read.

## State and removal

Catalog links have no state file. `.scio/hook_state/` and
`.scio/hook_settings.state` are ownership records for materialized hook
artifacts. On the first rebranded `link`, `.sciagent/` is moved to `.scio/`
before hook ownership is evaluated.

Project removal is manual: unlink the six category paths, remove complete
`SCIO` or legacy `SCIAGENT` managed blocks from `AGENTS.md`, and edit
materialized hook settings when the project retires enforcement.

## Project memory

`docs/_internal/` is its own Git repository in every project, always (ADR-D10).
That is what makes rewriting a record safe and what lets the owner publish the
reasoning on its own schedule; publication embeds it as a submodule.

Scope is the only structure: a directory per stage stem, plus `_project/` for
what spans or precedes stages, holding flat topic notes and `plans/`. `session.md`
is the live record and always carries that name; what it superseded lives in
`session-history/<UTC timestamp>.md`, so a plain sort is chronological order.

See [propagation.md](propagation.md) for the two-hop delivery model.
