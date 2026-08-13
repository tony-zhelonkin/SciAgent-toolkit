# sciagent

SciAgent is a per-project catalog of computational-biology skills, sub-agents,
commands, CRAFT conventions, and executable guardrails for AI coding harnesses.

The toolkit keeps its mount sources in the same shape the harnesses consume:
84 skill directories, 21 flat agent files, and 7 flat command files. A project
binds each category with one directory symlink.

## Install

Most analysis projects vendor this repository at
`01_modules/SciAgent-toolkit/`. Use a relative alias so each project invokes its
own pinned copy:

```bash
alias si='./01_modules/SciAgent-toolkit/bin/sciagent'
```

`sciagent link` refuses an external toolkit when the project contains its own
`SciAgent-toolkit` checkout. This keeps the binding aligned with the project's
submodule pin.

For one checkout, a PATH symlink is convenient:

```bash
ln -sf /absolute/path/to/SciAgent-toolkit/bin/sciagent ~/.local/bin/sciagent
sciagent --help
```

### From a release tarball

The offline installer accepts a local release artifact and checksum. The
release artifact is named `scio`; the installed command is `sciagent`.

```bash
./install.sh --archive  scio-0.1.0-<short-sha>.tar.gz \
             --checksum scio-0.1.0-<short-sha>.tar.gz.sha256 \
             --prefix   "$HOME/.local"
```

It installs a content-addressed version under
`~/.local/share/scio/versions/<full-git-sha>/`, links
`~/.local/bin/sciagent`, and writes a receipt. The installer has no network
path. Maintainers build an artifact with
`scripts/build-release.sh <ref>` from a clean tree.

## Quick start

```bash
# Bind the catalog and materialize the project guardrail hooks.
sciagent link

# Render the shared computational-biology conventions into AGENTS.md.
sciagent craft

# Run project checks.
sciagent lint
```

Every verb accepting a project path defaults to the current directory:

```bash
sciagent link --project-dir /path/to/project
sciagent craft --project-dir /path/to/project
sciagent lint --project-dir /path/to/project
```

## Verbs

| Verb | Description |
|---|---|
| `link [--project-dir D]` | Bind the six catalog trees and ensure guardrail hooks |
| `craft [--project-dir D] [--force] [--quiet]` | Render or refresh `SCIAGENT:CRAFT` in `AGENTS.md` |
| `lint [--project-dir D] [--check <name>...] [--strict] [--quiet]` | Run project guardrail checks |
| `validate [--quiet]` | Validate toolkit skill frontmatter and namespace collisions |
| `gitignore [<path>]` | Refresh the `SCIAGENT:GITIGNORE` block |

`validate` and `gitignore` are transitional surfaces during the CLI
demolition. Run `sciagent --help` for the terse reference.

## What `link` writes

```text
project/
├── .claude/
│   ├── skills   -> <toolkit>/skills
│   ├── agents   -> <toolkit>/agents
│   ├── commands -> <toolkit>/commands
│   ├── settings.json
│   └── hooks/
│       ├── no_ephemeral.sh
│       └── caption_sweep.sh
└── .agents/
    ├── skills   -> <toolkit>/skills
    ├── agents   -> <toolkit>/agents
    └── commands -> <toolkit>/commands
```

`link` is convergent. A missing category link is created, a correct link is a
silent no-op, and a link pointing elsewhere is replaced with a message. It
also sweeps legacy toolkit-owned child mounts, retired output-style links, and
dangling absolute mounts from older container paths.

A real populated category directory is preserved. `link` names every entry in
the refusal and asks the user to relocate it before retrying. Private skills,
agents, and commands should live outside these six category paths.

Hook bodies use a hash-and-cede ownership discipline. An unchanged body from
any shipped toolkit version can be refreshed. A user-edited body is preserved,
reported once, and ceded from future management. Hook registrations are merged
into `.claude/settings.json` while unrelated project settings remain intact.

Claude Code treats project settings as replacements for user settings except
for permission rules. SciAgent writes only the two project guardrail hook
registrations; editor, model, memory, attribution, thinking, and statusline
preferences remain user-owned.

## Manual removal

The `deactivate` verb has retired, so removing a project binding is manual.
The catalog and managed context come out in two operations:

```bash
rm .claude/{skills,agents,commands} .agents/{skills,agents,commands}
# Edit AGENTS.md and remove the complete BEGIN/END SCIAGENT managed blocks.
```

One link per category makes the filesystem portion a single explicit `rm`
line. The materialized hooks and their settings registrations remain project
guardrails; remove those files and registrations separately when retiring the
enforcement layer too.

## Validation and linting

`sciagent validate` checks every skill's `name` and `description` frontmatter.
The name must match the directory and the description must fit the configured
length cap. Cross-namespace collisions among skills, agents, commands, and
roles are warnings.

`sciagent lint` runs project checks for figure style, results layout, captions,
provenance, freshness, stage structure, comment intent, documentation layout,
and registered hook existence. Findings warn by default and become failures
under `--strict`.

## Harness support

Claude Code reads `.claude/{skills,agents,commands}`. Pi reads
`.agents/skills`; extensions can consume `.agents/agents` and
`.agents/commands`. Both trees point at the same toolkit sources.

## Documentation

- [Propagation model](docs/propagation.md)
- [Architecture](docs/architecture.md)
- [Skills catalog](docs/skills.md)
- [Agent catalog](docs/agents.md)
- [Command catalog](docs/commands.md)
