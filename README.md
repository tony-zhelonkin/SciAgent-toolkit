# sciagent

Per-project context manager for AI coding harnesses. 
Activates a role — a bundle of skills, sub-agents, and slash commands — for a session.

## Premise

I worked as computational biologist across different projects for two collaborating labs, 
with occasional collaborations on the side. Somehow I needed to manage and perform. 
Here came ChatGPT -> then Claude Code -> vibe coding -> agentic coding -> and we\`re
still falling down to see how deep the rabbit hole goes. 

Since 2024 I\`ve been exploring how I could make this at least somewhat more reproducible in my own hands.
This project is an attempt at sharpening the blunt tools.

## Why

Different projects make you put on different hats. 
Just as we have to switch contexts, so do LLM.
Different work modes need different context. 
Different data modalities, different goals, might require different context.
- different agents, 
- different skills, 
- different commands. 
Roles let you swap that context in one command. 
Everything is by design per-project folder.

Even with increased context window size I don\`t personally believe that 
universal agents are universally good. So I decided my agents need to have specialization.
But requirements would constantly change, so I needed something that allow me to control the
context management process. So this tool has been an exploration of the process

## Install

The CLI is one bash script at `bin/sciagent`. To call `sciagent` (or its short alias `si`) from anywhere, symlink it into a directory on your `PATH`.

`PATH` is the colon-separated list shell searches when a command is typed — `echo $PATH` prints it. 
`~/.local/bin` is the conventional spot for user-installed binaries and is already on `PATH` in most modern shells; if it isn't, add `export PATH="$HOME/.local/bin:$PATH"` to your shell rc (`~/.bashrc`, `~/.zshrc`). 

```bash
# Replace /absolute/path/to/SciAgent-toolkit with wherever you cloned (or submoduled) the toolkit
ln -sf /absolute/path/to/SciAgent-toolkit/bin/sciagent ~/.local/bin/sciagent
ln -sf /absolute/path/to/SciAgent-toolkit/bin/sciagent ~/.local/bin/si

sciagent --help     # verify it resolves
```

`ln -sf` is idempotent — re-run to repoint at a different checkout. Uninstall with `rm ~/.local/bin/sciagent ~/.local/bin/si`; the toolkit itself is untouched. 

## Quick start

```bash
# Bootstrap a new project directory with AGENTS.md, CLAUDE.md, docs/_internal/scientific-context.md
sciagent new project

# Activate a role — symlinks agents, skills, commands into .claude/ and .agents/
sciagent activate base

# Add one skill on top of the current stack
sciagent inject simplify

# Show the active stack, effective tables, and block/symlink health (add --json for a machine-readable manifest)
sciagent status

# List the roles available to activate, with per-role skill/agent/command counts
sciagent list roles

# Tear down: remove symlinks and the managed block from AGENTS.md
sciagent deactivate
```

## The RPG model

A project has at most two active roles: 
- a `base` (the foundation) 
- and an optional `overlay` (the specialization). 
Layering runs bottom-to-top — the overlay's entries shadow matching ones from the base. `sciagent status` shows what got shadowed.

Any role can occupy either slot; there's no enforced base/overlay typing. 
Toggle whichever combination fits the session — `base` + `pathway-signature` for downstream interpretation, `base` + `scatac-regulatory` to layer a chromatin stack on top of the scRNA foundation, `architect` solo for design sessions. `sciagent list roles` enumerates what's available.

```bash
sciagent activate base pathway-signature
```

Stack depth is for now capped at 2 to stay inspectable.

## Verbs

| Verb | Description |
|------|-------------|
| `activate <base> [overlay]` | Activate role(s); replaces current stack |
| `deactivate [<role>]` | Tear down the stack or remove one role |
| `inject <name>` | Add one skill / agent / command (auto-detect; `--skill` / `--agent` / `--command` for explicit) |
| `eject <name>` | Remove one injected entry (symmetric to `inject`) |
| `validate [--quiet]` | Check toolkit integrity (requires-graph, tags, refs, name collisions) |
| `status [--json\|--effective\|--source <name>]` | Report active stack and effective tables |
| `list [roles\|skills\|agents\|commands]` | List available content in the toolkit |
| `new project\|role\|skill\|agent [args]` | Scaffold from templates |

Run `sciagent --help` for the terse reference. `si` is available as an alias if you symlink `bin/sciagent` as `si` in your PATH.

## inject · eject · validate

`inject <name>` auto-detects whether `<name>` is a skill, agent, or command:

```
$ sciagent inject extra-skill
injected: extra-skill (into _injected)

$ sciagent inject extra-agent
injected: extra-agent (into _injected)

$ sciagent inject extra-command
injected: extra-command (into _injected)
```

Ambiguous names hard-fail; the explicit flags resolve them:

```
$ sciagent inject dual-name
error: ambiguous — 'dual-name' exists as both skill and command. use --skill <name>, --agent <name>, or --command <name>

$ sciagent inject --command dual-name
injected: dual-name (into _injected)
note: companion skill 'dual-name' available — `inject --skill dual-name` to add
```

The companion-skill note is informational — the skill is not auto-mounted. Unknown names hard-fail:

```
$ sciagent inject definitely-does-not-exist
error: 'definitely-does-not-exist' not found as skill, agent, or command
```

`eject <name>` is symmetric. Ambiguous when the same name was injected as 2+ kinds:

```
$ sciagent eject extra-agent
ejected: extra-agent (agent)

$ sciagent eject dual-name
error: ambiguous — 'dual-name' is injected as both command and skill. use --skill <name>, --agent <name>, or --command <name>
```

`validate` checks toolkit integrity. Cross-namespace name collisions are soft-warns (mounting both is supported); other failures are hard. `--quiet` suppresses all output (exit code only):

```
$ sciagent validate
validate: warning — name 'architect' appears as agent, command, and role (mounting both is supported; ensure the overlap is intentional)
validate: warning — name 'architecture-treemap' appears as both skill and command (mounting both is supported; ensure the overlap is intentional)
sciagent validate: all checks passed

$ sciagent validate --quiet
```

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
