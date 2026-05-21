# AGENTS.md - SciAgent-toolkit Codebase Instructions

**Version:** 2.1.0
**Last Updated:** 2026-05-20
**Purpose:** Universal AI agent instructions for working with the SciAgent-toolkit codebase

---

> **This file is for AI agents working on the SciAgent-toolkit codebase itself.**
> For project-specific analysis methodology, see `templates/vendor/AGENTS.md.template`.

---

## Repository Overview

**SciAgent-toolkit** is a role-based agent and skills framework for AI-assisted bioinformatics projects.

### Core Purpose

1. **Provide reusable agents** for bioinformatics workflows
2. **Provide a skills library** for single-cell and multi-omics analysis
3. **Manage roles** for different analysis contexts
4. **Template project setup** for consistent AI-assisted analysis

---

## Architecture

```
SciAgent-toolkit/
├── agents/          # Canonical agent definitions (.md files)
├── skills/          # Canonical skill definitions (directory format)
├── roles/           # Role definitions (YAML)
├── commands/        # Slash command definitions
├── templates/vendor/ # Project context templates
├── docs/guidelines/  # Modular methodology guidelines
├── scripts/
│   ├── setup-ai.sh          # Primary entry point
│   └── activate-role.sh     # Role activator
└── docker/test/             # CI/CD test infrastructure
```

### Key Directories

| Directory | Purpose |
|-----------|---------|
| `agents/` | Canonical agent definitions |
| `skills/` | Canonical skill definitions |
| `roles/` | Role definitions (YAML) |
| `templates/vendor/` | Project context templates |
| `docs/guidelines/` | Modular methodology guidelines |
| `docker/test/` | CI/CD test infrastructure |

---

## Critical Rules for Toolkit Development

### 1. Separation of Concerns

**Root-level files** describe the **toolkit codebase**:
- `AGENTS.md` (this file) — Universal instructions for toolkit development
- `CLAUDE.md` — Claude Code context for the toolkit
- `README.md` — Project overview

**Template files** (`templates/vendor/`) are for **user projects**:
- `AGENTS.md.template` — Comprehensive project methodology
- `CLAUDE.md.template` — Claude Code project context
- `context.md.template` — Scientific question scaffold

### 2. Template Placeholders

Use double-brace placeholders in templates. The `setup-ai.sh` script substitutes them:

| Placeholder | Substituted With |
|-------------|------------------|
| `{{PROJECT_ID}}` | Basename of project directory |
| `{{PROJECT_TITLE}}` | Same as PROJECT_ID |
| `{{DATE}}` | Current date (YYYY-MM-DD) |
| `{{SPECIES}}` | Default: "Mus musculus" |
| `{{EXPERIMENTAL_DESIGN}}` | Default: "TBD" |

### 3. Idempotent Scripts

All scripts MUST be idempotent:
- Check for existing state before acting
- Safe to run multiple times
- No destructive operations

### 4. Agent File Structure

Agents are Markdown files with YAML frontmatter:

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

---

## Development Workflow

### Adding a New Agent

1. Create `agents/new-agent-name.md` following the structure above
2. Add to relevant role(s) in `roles/*.yaml`
3. Test by activating the role: `sciagent activate <role>`
4. Document in `agents/README.md`

### Adding a New Role

1. Create `roles/new-role.yaml`:
   ```yaml
   name: new-role
   description: Role description

   agents:
     - agent-name-1
     - agent-name-2

   skills: []
   ```
2. Test activation: `sciagent activate new-role`

---

## Testing

```bash
cd docker/test
./test-all.sh
```

---

## File Reference

### Scripts

| Script | Purpose |
|--------|---------|
| `setup-ai.sh` | Primary entry point for project setup |
| `activate-role.sh` | Role activation (symlinks agents/skills) |

### Documentation

| File | Purpose |
|------|---------|
| `docs/guidelines/*.md` | Modular methodology guidelines |
| `docs/workflows/architect/` | Architect role workflow docs |

---

## Guidelines Reference

The `docs/guidelines/` directory contains modular methodology documentation:

| Module | Content |
|--------|---------|
| `core_architecture.md` | Phased workflow, directory structure |
| `data_processing.md` | filterByExpr, normalization, DE |
| `gsea_analysis.md` | GSEA patterns, msigdbr usage |
| `checkpoint_caching.md` | load_or_compute pattern |
| `master_tables.md` | CSV schema standardization |
| `visualization.md` | Colors, themes, plots |
| `code_style.md` | R/Python conventions |

These guidelines are referenced from project templates but maintained here as the single source of truth.

---

## Version History

- **2.1.0** (2026-05-20): Removed MCP infrastructure and harness installer docs; toolkit now covers roles/agents/skills only
- **2.0.0** (2025-12-16): Restructured as toolkit codebase documentation; project methodology moved to templates
- **1.0.0** (2025-12-10): Initial version
