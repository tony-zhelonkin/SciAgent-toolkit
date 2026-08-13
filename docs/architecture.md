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
├── roles/<name>.yaml
├── craft.yaml
├── templates/
└── tests/
```

The three mount directories contain mountable content only. Their source shape
is identical to the consumer shape: skill directories and flat Markdown files
for agents and commands.

Roles are retained data during the demolition. They do not select catalog
membership. The transitional `validate` verb reads role basenames only when it
checks cross-namespace name collisions.

## CLI

The target surface has three verbs:

- `link` binds the catalog and ensures project hooks.
- `craft` renders the managed CRAFT block.
- `lint` runs project checks.

Two transitional verbs remain in this step: `validate` for the definition-of-
done gate and `gitignore` for its existing managed block.

The dispatcher resolves its own real path, derives `$SCIAGENT_TOOLKIT`, sources
the small dependency closure for one verb, and calls that verb's `cmd_*`
entrypoint.

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

The operation is convergent:

1. Sweep legacy toolkit-owned child links from known mount directories.
2. Keep a correct category link silently.
3. Replace an outside-pointing category link and report the change.
4. Collapse an emptied legacy directory into the category link.
5. Preserve and refuse a real directory containing user entries.

The sweep recognizes links resolving inside the current toolkit and dangling
legacy paths under a `SciAgent-toolkit` source. It covers retired output-style
and analysis-helper links. Regular files and outside-pointing links survive.

A locality check protects submodule pins. When a project carries its own
`SciAgent-toolkit`, an external installation cannot bind that project.

## Guardrail hooks

`link` materializes `no_ephemeral.sh` and `caption_sweep.sh`, then merges their
registrations into `.claude/settings.json`. The settings merge preserves every
unrelated project setting.

Hook bodies follow a hash-and-cede rule. Current or historically shipped bytes
are toolkit-owned and refreshable. Unrecognized bytes are user-owned, retained,
and ceded after one warning. `templates/PROVENANCE.sha1` carries historical
hashes for release archives and shallow clones.

## Managed blocks

`block.sh` owns marker-framed updates in shared text files. A block begins with
an HTML marker containing its id, marker version, and SHA1 of the body. Readers
verify the stored hash against the body itself.

`craft` writes the `SCIAGENT:CRAFT` block in `AGENTS.md`. A consistent older
block can be refreshed. A hand-edited body fails its stored hash and requires
`--force`. Text outside the markers and file permissions are preserved.

## Lint and validation

`lint` contains project guardrail checks, including figure style, results
layout, captions, provenance, freshness, stage structure, comment intent,
documentation layout, and registered hook existence. Findings warn by default
and become failures under `--strict`.

`validate` checks toolkit skill frontmatter and cross-namespace collisions. It
remains available while the demolition definition of done calls it directly.

## State and removal

Catalog links have no state file. Their target is the complete binding.
`.sciagent/hook_state/` and `.sciagent/hook_settings.state` are ownership
records for materialized hook artifacts.

Project removal is manual: unlink the six category paths, then remove complete
SCIAGENT managed blocks from `AGENTS.md`. Materialized hooks and their settings
registrations are removed separately when the project retires enforcement.

See [propagation.md](propagation.md) for the two-hop delivery model.
