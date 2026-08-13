# SciAgent toolkit architecture

## Purpose

SciAgent supplies project-local computational-biology context and guardrails
to AI coding harnesses. Git is the version authority. A consumer normally pins
the toolkit as a submodule and binds that exact checkout into its harness
discovery paths.

## Source layout

```text
SciAgent-toolkit/
├── bin/sciagent
├── lib/sciagent/
├── skills/<name>/SKILL.md
├── agents/<name>.md
├── commands/<name>.md
├── craft.yaml
├── templates/
└── tests/
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

The dispatcher resolves its own real path, derives `$SCIAGENT_TOOLKIT`, sources
the dependency closure for one verb, and calls that verb's `cmd_*` entrypoint.

## Link topology

`sciagent link [--project-dir D]` creates six directory symlinks:

```text
D/.claude/skills   -> $SCIAGENT_TOOLKIT/skills
D/.claude/agents   -> $SCIAGENT_TOOLKIT/agents
D/.claude/commands -> $SCIAGENT_TOOLKIT/commands
D/.agents/skills   -> $SCIAGENT_TOOLKIT/skills
D/.agents/agents   -> $SCIAGENT_TOOLKIT/agents
D/.agents/commands -> $SCIAGENT_TOOLKIT/commands
```

The operation is convergent. It sweeps legacy toolkit-owned child links,
preserves populated directories and user-owned links, and silently keeps
correct bindings. A locality check protects submodule pins.

Analysis projects also receive shared helper-library links and import shims.
`link` refreshes the `SCIAGENT:GITIGNORE` block because those bindings and the
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
`SCIAGENT:CRAFT` block in `AGENTS.md`; a hand-edited body requires `--force`.
Text outside the markers and file permissions are preserved.

Historical `SCIAGENT:ROLES` blocks remain readable and removable through
`block_read` and `block_remove`. The writer accepts CRAFT alone.

## Lint

`lint` contains project guardrail checks for figure style, results layout,
captions, provenance, freshness, stage structure, comment intent,
documentation layout, and registered hook existence. Findings warn by default
and become failures under `--strict`.

`sciagent lint --check toolkit` checks skill frontmatter shape and reports
basename collisions across skills, agents, and commands. The release builder
uses this check as its catalog gate.

## State and removal

Catalog links have no state file. `.sciagent/hook_state/` and
`.sciagent/hook_settings.state` are ownership records for materialized hook
artifacts.

Project removal is manual: unlink the six category paths, remove complete
SCIAGENT managed blocks from `AGENTS.md`, and edit materialized hook settings
when the project retires enforcement.

See [propagation.md](propagation.md) for the two-hop delivery model.
