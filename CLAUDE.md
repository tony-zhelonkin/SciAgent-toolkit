# CLAUDE.md

This file provides guidance to Claude Code when working with the SciAgent-toolkit codebase.

---

## Repository Overview

**SciAgent-toolkit** is a per-project context manager for AI coding harnesses. It activates a role — a bundle of skills, sub-agents, and slash commands — into a project's `.claude/` and `.agents/` directories.

**Integration:** Used as a submodule at `01_modules/SciAgent-toolkit/` in analysis projects.

---

## Key Directories

| Path | Purpose |
|------|---------|
| `bin/sciagent` | CLI dispatcher |
| `lib/sciagent/` | Internal bash modules (`activate.sh`, `deactivate.sh`, `inject.sh`, `status.sh`, `new.sh`, `block.sh`, `roles.sh`, `symlinks.sh`, `stack.sh`) |
| `roles/` | Role definitions (YAML) |
| `agents/` | Canonical sub-agent definitions (`.md` files) |
| `skills/` | Canonical skill definitions (`<name>/SKILL.md` format) |
| `commands/` | Canonical slash command definitions (`.md` files) |
| `system-prompts/` | System-prompt / output-style definitions (`.md` files; Claude Code consumer mounts into `.claude/output-styles/`) |
| `templates/` | Project scaffolding templates (`AGENTS.md.template`, `CLAUDE.md.template`, `context.md.template`) |
| `tests/` | Bash test suite (`run-all.sh`) |

---

## Common Commands

```bash
# Activate base role (symlinks agents/skills/commands into .claude/ and .agents/)
bin/sciagent activate base

# Run the test suite
bash tests/run-all.sh

# Scaffold a new role
bin/sciagent new role my-role

# Check active stack
bin/sciagent status
```

---

## Role System

Role YAML schema:

```yaml
name: base
description: Default bioinformatics analysis role

skills:
  - skill-name
agents:
  - agent-name
commands:
  - command-name
# output_style: architect-mentor   # optional, Claude-specific
```

Stack is capped at 2 (base + optional overlay). Last-wins on name collisions.

### Adding a New Agent

1. Create `agents/<name>.md` with YAML frontmatter (`name`, `description`, `model`, `color`)
2. Add to relevant `roles/*.yaml`
3. Test: `bin/sciagent activate base && echo OK`
4. Update `agents/README.md`

### Adding a New Role

1. `bin/sciagent new role <name>` — scaffolds `roles/<name>.yaml`
2. Edit the YAML to add skills/agents/commands
3. Test: `bin/sciagent activate <name>`

### Adding a New Skill

1. `bin/sciagent new skill <name>` — copies `skills/_TEMPLATE/` to `skills/<name>/`
2. Edit `skills/<name>/SKILL.md`
3. Add to relevant `roles/*.yaml`

---

## After Modifying Agents, Skills, or Roles

Test activation from the toolkit directory:

```bash
bin/sciagent activate base && echo "OK" || echo "FAILED"
ls -la .claude/agents/ .claude/skills/
```

If this toolkit is a submodule, re-activate from the parent:

```bash
./01_modules/SciAgent-toolkit/bin/sciagent activate base
```

---

## Architecture Reference

See [docs/architecture.md](docs/architecture.md) for the full design spec.
