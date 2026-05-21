# sciagent

Per-project context manager for AI coding harnesses. Activates a role — a bundle of skills, sub-agents, and slash commands — for a session.

## Why

Different work modes need different context. Reviewing code and running bioinformatics analysis call for different agents, different skills, different commands. Roles let you swap that context in one command. Everything is per-project — nothing touches your home directory.

## Quick start

```bash
# Bootstrap a new project directory with AGENTS.md, CLAUDE.md, context.md
sciagent new project

# Activate a role — symlinks agents, skills, commands into .claude/ and .agents/
sciagent activate base

# Add one skill on top of the current stack
sciagent inject simplify

# Show the active stack, effective tables, and block/symlink health
sciagent status

# Tear down: remove symlinks and the managed block from AGENTS.md
sciagent deactivate
```

## The RPG model

A project has at most two active roles: a `base` (the foundation) and an optional `overlay` (the specialization). Last-wins on name collisions; `sciagent status` shows what got shadowed.

```bash
sciagent activate base reviewer
```

`base` provides bioinformatics context; `reviewer` overlays code-review agents and commands. Stack depth is capped at 2 to stay inspectable — Claude Code's own three-tier resolution already makes "where did this come from?" painful enough.

## Verbs

| Verb | Description |
|------|-------------|
| `activate <base> [overlay]` | Activate role(s); replaces current stack |
| `deactivate [<role>]` | Tear down the stack or remove one role |
| `inject <skill>` | Add one skill to the current overlay |
| `status [--json\|--effective\|--source <name>]` | Report active stack and effective tables |
| `list [roles\|skills\|agents\|commands]` | List available content in the toolkit |
| `new project\|role\|skill\|agent [args]` | Scaffold from templates |

Run `sciagent --help` for the terse reference. `si` is available as an alias if you symlink `bin/sciagent` as `si` in your PATH.

## What it writes

`sciagent activate` creates symlinks and a managed block in `AGENTS.md`:

```
project/
├── AGENTS.md                         # your file; sciagent appends a managed block
├── CLAUDE.md                         # 1-line @AGENTS.md shim (from template)
├── .claude/
│   ├── skills/<name>  →  toolkit/skills/<name>
│   ├── agents/<name>.md  →  toolkit/agents/<name>.md
│   ├── commands/<name>.md  →  toolkit/commands/<name>.md
│   └── output-styles/<name>.md  →  toolkit/system-prompts/<name>.md
├── .agents/
│   ├── skills/<name>  →  toolkit/skills/<name>
│   ├── agents/<name>.md  →  toolkit/agents/<name>.md
│   └── commands/<name>.md  →  toolkit/commands/<name>.md
└── .sciagent/manifest.json           # machine-readable state for safe teardown
```

The managed block is delimited by HTML comments (`<!-- BEGIN SCIAGENT:ROLES v1 hash=... -->`), invisible in rendered markdown. On each run, sciagent recomputes the hash and warns if you've edited inside the block.

`deactivate` removes the block and removes only symlinks it owns (tracked via `manifest.json`).

## Harness support

Reads the native discovery directories of Claude Code (`.claude/skills/`, `.claude/agents/`, `.claude/commands/`) and Pi (`.agents/skills/`). Sub-agents and slash commands are Claude-specific today. Pi extension for `.agents/agents/` and `.agents/commands/` is the user's job — symlinks are already there.

## Architecture

See [docs/architecture.md](docs/architecture.md) for the full design spec.
