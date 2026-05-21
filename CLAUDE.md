# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

---

## Repository Overview

**SciAgent-toolkit** is a role-based agent and skills framework for AI-assisted bioinformatics projects.

**Integration:** Used as a submodule in [scbio-docker](https://github.com/tony-zhelonkin/scbio-docker) at `toolkits/SciAgent-toolkit/` and added to new projects at `01_modules/SciAgent-toolkit/`.

---

## Key Directories and Files

### Scripts
- `sciagent new project` — **Primary entry point** for project setup (templates + role activation)
- `sciagent activate` — **Role activator** for agent/skill symlinks

### Role System
- `roles/` — Role definitions (YAML files)
- `roles/base.yaml` — Default bioinformatics analysis role
- `agents/` — Canonical location for custom Claude agents (flat `.md` files)
- `skills/` — Canonical location for custom Claude skills (directory format)
  - `skills/<name>/SKILL.md` — canonical entry point for each skill
  - `skills/<name>/{references,scripts,checks,assets}/` — optional supporting dirs
  - `skills/_TEMPLATE/` — starter template (copy to author new skills)
  - `skills/skill-creator/` — reference Rich-tier implementation
- `.claude/agents/` — Symlinked agents (populated by `activate-role.sh`)
- `.claude/skills/` — Symlinked skills (populated by `activate-role.sh`)
  - Directory-format skills are symlinked as directories (preserves supporting dirs)
  - Legacy flat-format skills are still supported; see `skills/README.md` for the format guide

### Template System
- `templates/vendor/` — AI context templates installed by `setup-ai.sh`
  - `CLAUDE.md.template` — Claude Code project instructions
  - `AGENTS.md.template` — Universal AI rules for all agents
  - `context.md.template` — Scientific project context
  - `analysis_config.yaml.template` — Analysis parameters for `02_analysis/config/`

### Docker Testing Infrastructure
- `docker/test/` — CI/CD test Dockerfiles and test suite

---

## Common Development Commands

### Setup

```bash
# Full project setup (templates + role activation)
sciagent new project

# View options
sciagent new --help
```

### Role Activation

```bash
# Activate base role (default for bioinformatics)
sciagent activate base

# List available roles
ls roles/*.yaml
```

### Docker Testing

```bash
cd docker/test
./test-all.sh
```

---

## Role System

The role system provides a declarative way to configure agents and skills per project.

### Role Definition (`roles/base.yaml`)

```yaml
name: base
description: Default bioinformatics analysis role with full agent suite

agents:
  - bioinf-librarian
  - bio-research-visualizer
  - rnaseq-insight-explorer
  - rnaseq-methods-writer
  - figure-caption-generator
  - repo-doc-curator
  - refactor-stage-reviewer
  - handoff

skills: []
```

### How Role Activation Works (`activate-role.sh`)

1. **Read role YAML**: Parses `roles/<role>.yaml`
2. **Create directories**: `.claude/agents/` and `.claude/skills/`
3. **Clear existing symlinks**: Removes old role configuration
4. **Symlink agents**: Links each agent from `agents/<name>.md` to `.claude/agents/`
5. **Symlink skills** (directory-first, flat fallback):
   - If `skills/<name>/SKILL.md` exists → symlinks the whole directory to `.claude/skills/<name>`
   - Otherwise if `skills/<name>.md` exists → symlinks the file to `.claude/skills/<name>.md`

### Creating Custom Roles

```yaml
# roles/my-custom-role.yaml
name: my-custom-role
description: Custom role for specific workflow

agents:
  - bioinf-librarian
  - my-custom-agent

skills:
  - my-custom-skill
```

---

## Custom Agents

### Agent File Structure

Agents are Markdown files with YAML frontmatter in the `agents/` directory:

```markdown
---
name: "agent-identifier"
description: "When to use this agent (with examples)"
model: "sonnet" | "opus" | "haiku"
color: "yellow" | "blue" | "green"
---

# Agent Identity
...

# Methodology
...
```

### Creating New Agents

1. Create `agents/new-agent-name.md` following the structure above
2. Add to relevant role(s) in `roles/*.yaml`
3. Test activation: `sciagent activate base`
4. Document in `agents/README.md`

---

## After Modifying Agents or Roles

**IMPORTANT**: After making ANY changes to files in `agents/`, `skills/`, or `roles/`, test role activation:

```bash
# One-liner to test from toolkit directory
sciagent activate base && echo "Role activation: OK" || echo "Role activation: FAILED"
```

Verify symlinks:
```bash
ls -la .claude/agents/
ls -la .claude/skills/
```

If a parent project uses this toolkit as a submodule, re-activate from the parent directory:
```bash
./01_modules/SciAgent-toolkit/bin/sciagent activate base
```

---

## Guidelines Reference

The `docs/guidelines/` directory contains modular methodology documentation for bioinformatics projects.
