# AGENTS.md - SciAgent-toolkit Codebase Instructions

**Version:** 2.0.0
**Last Updated:** 2025-12-16
**Purpose:** Universal AI agent instructions for working with the SciAgent-toolkit codebase

---

> **This file is for AI agents working on the SciAgent-toolkit codebase itself.**
> For project-specific analysis methodology, see `templates/vendor/AGENTS.md.template`.

---

## Repository Overview

**SciAgent-toolkit** is a modular MCP (Model Context Protocol) infrastructure orchestrator that integrates AI assistants (Claude Code, Gemini CLI, Codex CLI) with specialized scientific research tools for bioinformatics projects.

### Core Purpose

1. **Install and configure AI CLI tools** (Claude Code, Gemini CLI, Codex CLI)
2. **Orchestrate MCP servers** (PAL, Sequential Thinking, Context7, ToolUniverse, Serena)
3. **Provide reusable agents** for bioinformatics workflows
4. **Manage roles and profiles** for different analysis contexts
5. **Template project setup** for consistent AI-assisted analysis

---

## Architecture

### Four-Tier System

```
┌─────────────────────────────────────────────────────────┐
│  User Interfaces: Claude Code, Gemini CLI, Codex CLI    │
├─────────────────────────────────────────────────────────┤
│  Role & Profile System: activate-role.sh, switch-mcp.sh │
├─────────────────────────────────────────────────────────┤
│  MCP Servers: PAL, Sequential Thinking, Context7,       │
│               ToolUniverse, Serena                      │
├─────────────────────────────────────────────────────────┤
│  External APIs: ChEMBL, UniProt, PubMed, DrugBank, etc. │
└─────────────────────────────────────────────────────────┘
```

### Key Directories

| Directory | Purpose |
|-----------|---------|
| `scripts/` | Installation and configuration scripts |
| `scripts/mcp_servers/` | Individual MCP server installers |
| `templates/vendor/` | Project context templates (AGENTS.md, CLAUDE.md, etc.) |
| `templates/mcp-profiles/` | MCP configuration profiles |
| `agents/` | Canonical agent definitions |
| `skills/` | Canonical skill definitions |
| `roles/` | Role definitions (YAML) |
| `docs/guidelines/` | Modular methodology guidelines |
| `docker/test/` | CI/CD test infrastructure |

---

## Critical Rules for Toolkit Development

### 1. Separation of Concerns

**Root-level files** describe the **toolkit codebase**:
- `AGENTS.md` (this file) - Universal instructions for toolkit development
- `CLAUDE.md` - Claude Code context for the toolkit
- `README.md` - Project overview

**Template files** (`templates/vendor/`) are for **user projects**:
- `AGENTS.md.template` - Comprehensive project methodology
- `CLAUDE.md.template` - Claude Code project context
- `GEMINI.md.template` - Gemini project context
- `context.md.template` - Scientific question scaffold

### 2. Absolute Paths in MCP Configurations

**Critical**: `.mcp.json` files MUST use absolute paths because Claude Code doesn't guarantee working directory.

```json
{
  "mcpServers": {
    "tooluniverse": {
      "command": "uv",
      "args": ["--directory", "/absolute/path/to/tooluniverse-env", "run", "tooluniverse-mcp"]
    }
  }
}
```

### 3. Template Placeholders

Use double-brace placeholders in templates. The `setup-ai.sh` script substitutes them:

| Placeholder | Substituted With |
|-------------|------------------|
| `{{PROJECT_ID}}` | Basename of project directory |
| `{{PROJECT_TITLE}}` | Same as PROJECT_ID |
| `{{DATE}}` | Current date (YYYY-MM-DD) |
| `{{SPECIES}}` | Default: "Mus musculus" |
| `{{EXPERIMENTAL_DESIGN}}` | Default: "TBD" |

### 4. Idempotent Scripts

All installation scripts MUST be idempotent:
- Check for existing installations before installing
- Safe to run multiple times
- No destructive operations
- Support `--force` flag for reinstall

### 5. Agent File Structure

Agents are Markdown files with YAML frontmatter:

```markdown
---
name: "agent-identifier"
description: "When to use this agent (with examples)"
model: "sonnet" | "opus" | "haiku"
color: "yellow" | "blue" | "green"
---

# Agent Identity
Your agent's role and capabilities...

# Methodology
How the agent works...

# Examples
Example use cases...
```

---

## Development Workflow

### Adding a New Agent

1. Create `agents/new-agent-name.md` following the structure above
2. Add to relevant role(s) in `roles/*.yaml`
3. Test by activating the role: `./scripts/activate-role.sh <role>`
4. Document in `agents/README.md`

### Adding a New MCP Server

1. Create installer: `scripts/mcp_servers/setup_newserver.sh`
2. Add to orchestrator: `scripts/setup_mcp_infrastructure.sh`
3. Update config generator: `scripts/configure_mcp_servers.sh`
4. Create MCP profile(s): `templates/mcp-profiles/*.json`
5. Update documentation

### Adding a New Role

1. Create `roles/new-role.yaml`:
   ```yaml
   name: new-role
   description: Role description
   mcp_profile: coding  # or research-lite, etc.

   agents:
     - agent-name-1
     - agent-name-2

   skills: []
   ```
2. Test activation: `./scripts/activate-role.sh new-role`

### Adding a New MCP Profile

1. Create `templates/mcp-profiles/new-profile.mcp.json`
2. Add Gemini mapping if needed: `templates/gemini-profiles/new-profile.json`
3. Update `switch-mcp-profile.sh` if special handling needed

---

## Testing

### Architecture Tests

```bash
cd docker/test
./test-all.sh
```

### Individual Component Tests

```bash
# Test ToolUniverse
./scripts/test_tooluniverse.sh

# Test full installation
./scripts/test_installation.sh

# Validate MCP config
python3 -m json.tool .mcp.json
```

### Docker Test Images

| Dockerfile | Purpose |
|------------|---------|
| `Dockerfile.architecture-test` | Validates roles, templates, profiles |
| `Dockerfile.tooluniverse-test` | Base image with ToolUniverse |
| `Dockerfile.claude-test` | Claude Code integration |
| `Dockerfile.codex-test` | Codex CLI integration |
| `Dockerfile.gemini-test` | Gemini CLI integration |

---

## File Reference

### Scripts

| Script | Purpose |
|--------|---------|
| `setup-ai.sh` | **Primary entry point** for project setup |
| `setup_mcp_infrastructure.sh` | MCP server installation orchestrator |
| `activate-role.sh` | Role activation (symlinks agents/skills) |
| `switch-mcp-profile.sh` | Profile switcher for context management |
| `configure_mcp_servers.sh` | MCP configuration generator |
| `install_claude.sh` | Claude Code installer |
| `install_gemini.sh` | Gemini CLI installer |
| `install_codex.sh` | Codex CLI installer |

### Configuration Files

| File | Scope | Purpose |
|------|-------|---------|
| `.mcp.json` | Project-local | Claude Code MCP configuration |
| `.gemini/settings.json` | Project-local | Gemini CLI configuration |
| `~/.codex/config.toml` | User-global | Codex CLI configuration |
| `.env` | Project-local | API keys (git-ignored) |

### Documentation

| File | Purpose |
|------|---------|
| `docs/ARCHITECTURE.md` | System design |
| `docs/CONFIGURATION.md` | Advanced configuration |
| `docs/INSTALLATION.md` | Setup instructions |
| `docs/TROUBLESHOOTING.md` | Common issues |
| `docs/guidelines/*.md` | Modular methodology guidelines |

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

- **2.0.0** (2025-12-16): Restructured as toolkit codebase documentation; project methodology moved to templates
- **1.0.0** (2025-12-10): Initial version with mixed toolkit/project content
