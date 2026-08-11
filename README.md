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
Roles used to gate which of these got mounted; they don't anymore — the whole
catalog is always mounted (see "The RPG model" below), and activating a role now
labels provenance. The output-style is a separate choice at activation time
(`--output-style <name>`, else `craft.yaml`'s `output_style:`), not a property of
the role.
Everything is by design per-project folder.

Even with increased context window size I don\`t personally believe that 
universal agents are universally good. So I decided my agents need to have specialization.
But requirements would constantly change, so I needed something that allow me to control the
context management process. So this tool has been an exploration of the process

## Install

The CLI is one bash script at `bin/sciagent`. How you put `sciagent` (or its short alias
`si`) on your `PATH` depends on whether you work against one checkout or many.

### Working across multiple projects (recommended default)

Most projects vendor their own copy at `01_modules/SciAgent-toolkit/` (submodule pin) so
activation is reproducible. If you work across several such projects — or inside an
**umbrella** that nests several sibling projects, each with its own toolkit copy, under
one workspace (`docs/architecture.md` §13) — do **not** symlink `si` to one fixed
checkout. `bin/sciagent` derives `$SCIAGENT_TOOLKIT` from its own file location when the
env var isn't set, so a PATH symlink bakes a single absolute checkout as the target for
every invocation, no matter which project's directory you're actually in. Run it from a
different project and mutating verbs get refused by `_guard_toolkit_locality`
("refusing to activate against an external toolkit") — or, with
`--allow-external-toolkit`, silently activate that project against the wrong copy.

Define `si` as a relative-path shell alias instead:

```bash
alias si='./01_modules/SciAgent-toolkit/bin/sciagent'
```

Relative to `$PWD`, this always resolves to whichever project's own in-repo toolkit
you're currently `cd`'d into — correct standalone, and correct inside an umbrella
container that vendors several sibling copies at once. No `SCIAGENT_TOOLKIT` export, no
per-project symlink, and it can't mutate the wrong project's `.claude`/`.agents`.

### Single checkout, single project

If you only ever work against one toolkit checkout, a `PATH` symlink is simpler:

`PATH` is the colon-separated list shell searches when a command is typed — `echo $PATH` prints it. 
`~/.local/bin` is the conventional spot for user-installed binaries and is already on `PATH` in most modern shells; if it isn't, add `export PATH="$HOME/.local/bin:$PATH"` to your shell rc (`~/.bashrc`, `~/.zshrc`). 

```bash
# Replace /absolute/path/to/SciAgent-toolkit with wherever you cloned (or submoduled) the toolkit
ln -sf /absolute/path/to/SciAgent-toolkit/bin/sciagent ~/.local/bin/sciagent
ln -sf /absolute/path/to/SciAgent-toolkit/bin/sciagent ~/.local/bin/si

sciagent --help     # verify it resolves
```

`ln -sf` is idempotent — re-run to repoint at a different checkout. Uninstall with `rm ~/.local/bin/sciagent ~/.local/bin/si`; the toolkit itself is untouched. Once a second vendored project enters the picture, switch to the alias form above.

## Quick start

```bash
# Bootstrap a new project directory with AGENTS.md, CLAUDE.md, docs/_internal/scientific-context.md
sciagent new project

# Activate a role — symlinks the whole skill/agent/command catalog into .claude/ and .agents/;
# the role only decides provenance labels (output-style is a separate --output-style choice)
sciagent activate base

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
Layering runs bottom-to-top for *provenance*: a role's `skills:`/`agents:`/`commands:`
lists decide which role a mounted entry is attributed to (and which entry wins the
attribution on a name collision), not whether it gets mounted — the full catalog is
always mounted regardless of the active stack. `sciagent status` shows what got shadowed.

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
| `validate [--quiet]` | Check skill frontmatter shape + cross-namespace name collisions |
| `lint [--project-dir D] [--check <name>...]` | Opt-in PROJECT guardrail checks against an analysis repo |
| `status [--json\|--effective\|--source <name>]` | Report active stack and effective tables |
| `list [roles\|skills\|agents\|commands]` | List available content in the toolkit |
| `new project\|role\|skill\|agent [args]` | Scaffold from templates |
| `craft [--project-dir D]` | Render/refresh the SCIAGENT:CRAFT block in AGENTS.md |
| `gitignore [<path>]` | Add/update the SCIAGENT:GITIGNORE block in .gitignore |
| `update [--to <ref>]` | Re-pin the toolkit submodule + re-activate the current stack |
| `provision [--harness <csv\|all>]` | Seed user-level / global context + settings per harness |

Run `sciagent --help` for the terse reference.

## validate

There is no `inject`/`eject` verb anymore — the whole catalog is always mounted, so
there is nothing left to add or remove on top of a role's stack (see "The RPG model").

`validate` checks the frontmatter shape of every skill in the toolkit (matching
`name:`, a `description:` under the length cap). Cross-namespace name collisions
are soft-warns (mounting both is supported); frontmatter failures are hard.
`--quiet` suppresses all output (exit code only):

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
├── CLAUDE.md                         # created, or `@AGENTS.md` prepended to yours
├── .claude/
│   ├── skills/<name>  →  toolkit/skills/<name>
│   ├── agents/<name>.md  →  toolkit/agents/<name>.md
│   ├── commands/<name>.md  →  toolkit/commands/<name>.md
│   ├── output-styles/<name>.md  →  toolkit/system-prompts/<name>.md
│   ├── settings.json                 # created, or missing keys backfilled
│   ├── statusline.sh                 # created, +x
│   └── hooks/{no_ephemeral,caption_sweep}.sh   # created, +x, registered in settings.json
├── .agents/
│   ├── skills/<name>  →  toolkit/skills/<name>
│   ├── agents/<name>.md  →  toolkit/agents/<name>.md
│   └── commands/<name>.md  →  toolkit/commands/<name>.md
├── 02_analysis/helpers/{figure-style,interactive-style}  →  toolkit/lib/<name>
└── .sciagent/manifest.json           # records the active stack + the links created
```

The managed block is delimited by HTML comments (`<!-- BEGIN SCIAGENT:ROLES v1 hash=... -->`), invisible in rendered markdown.

`deactivate` removes both managed blocks and every symlink whose target resolves
inside the toolkit — ownership is derived from the link target, so a lost or
stale `manifest.json` cannot strand a mount. Anything that is not such a symlink
is left alone, which is why a file or symlink of your own inside `.claude/skills/`
survives untouched.

**`deactivate` is not a full inverse of `activate`.** The four `.claude/` entries
above marked *created* — `settings.json`, `statusline.sh`, both hooks — and the
`@AGENTS.md` line in `CLAUDE.md` are **not** removed, and the hooks stay
registered and live. There is currently no verb that unwinds them; remove them by
hand if you want a clean tree. Only the output-style key in
`settings.local.json` is reversed, and only when an output-style was actually
applied.

## Harness support

Reads the native discovery directories of Claude Code (`.claude/skills/`, `.claude/agents/`, `.claude/commands/`) and Pi (`.agents/skills/`). Sub-agents and slash commands are Claude-specific today. Pi extension for `.agents/agents/` and `.agents/commands/` is the user's job — symlinks are already there.

## Architecture

See [docs/architecture.md](docs/architecture.md) for the full design spec.
