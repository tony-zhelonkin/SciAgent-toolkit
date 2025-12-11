# SciAgent Toolkit - Complete Installation Guide

Complete setup and configuration guide for creating customized AI scientist agents using Claude Code, Gemini CLI, Codex CLI, and scientific MCP servers.

## Prerequisites

Your system should have:
- Linux or macOS (WSL supported on Windows)
- Python 3.10 or later
- Internet connection
- sudo access (for some system package installations on Linux)

## Installation Overview

The toolkit will install:
1. **Claude Code** (optional) - IDE/CLI interface with advanced reasoning
2. **Gemini CLI** (optional) - Google's advanced multimodal CLI agent
3. **Codex CLI** (optional) - Terminal-based scientific research interface
4. **Serena MCP Server** - Semantic code search and editing
5. **Sequential Thinking MCP Server** - Structured reasoning for complex decisions
6. **ToolUniverse MCP Server** - 600+ scientific research tools
7. **PAL MCP Server** - Multi-model collaboration and planning
8. **PubMed Plugin** - Biomedical literature access (manual setup)

## Full Installation (Recommended)

```bash
# Make the script executable (if needed)
chmod +x scripts/setup-ai.sh

# Run the installer
./scripts/setup-ai.sh
```

The installer will:
1. Check for required dependencies (Python, Node.js)
2. Setup environment variables (copying `.env.template` if needed)
3. Install selected AI CLIs (Claude, Gemini, Codex)
4. Install and configure all MCP servers (ToolUniverse, Serena, PAL, etc.)
5. Apply the default "coding" profile to all agents

## Selective Installation Options

### Minimal Setup (Fastest)
Skips Serena (saves ~15m compilation time) and Codex/Gemini.
```bash
./scripts/setup-ai.sh --minimal
```

### Skip Specific Components
```bash
./scripts/setup-ai.sh --skip-codex --skip-serena
```

## Component Details

### 1. Claude Code
**Purpose**: IDE/CLI interface with advanced reasoning.
**Installation**: Native installer via https://claude.ai/install.sh
**Location**: `~/.local/bin/claude`
**Verification**: `claude --version`

### 2. Gemini CLI
**Purpose**: Google's advanced multimodal CLI agent with large context window.
**Installation**: via npm package `@google/gemini-cli`
**Location**: `~/.npm-global/bin/gemini`
**Verification**: `gemini --version`

### 3. Codex CLI
**Purpose**: Terminal-based scientific research interface.
**Installation**: npm or Homebrew
**Location**: Global npm package or Homebrew installation
**Verification**: `codex --version`

### 4. MCP Servers (Core & Scientific)
- **Serena**: Semantic code search (requires `uvx`).
- **Sequential Thinking**: Structured reasoning (requires `npx`).
- **ToolUniverse**: 600+ scientific tools (Drug Discovery, Genomics, Literature).
- **PAL**: Multi-model collaboration (orchestrate Gemini from Claude).

## Post-Installation Configuration

### 1. API Keys (Crucial)
The setup script creates a `.env` file from a template if one doesn't exist. **You must edit this file** to add your API keys:

```bash
nano .env
```
Required keys for full functionality:
- `GEMINI_API_KEY`: For Gemini CLI and PAL
- `OPENAI_API_KEY`: For Codex and PAL fallback
- `CONTEXT7_API_KEY`: For Context7 (optional but recommended)

### 2. Profile Selection
We use a unified profile system to manage context usage across all agents.

```bash
# Default Coding Profile (Lightweight)
./scripts/switch-mcp-profile.sh coding

# Research Profile (ToolUniverse enabled)
./scripts/switch-mcp-profile.sh research-lite

# Hybrid Profile (Claude=Coding, Gemini=Research)
./scripts/switch-mcp-profile.sh hybrid-research
```

See [MCP-CONTEXT-MANAGEMENT.md](MCP-CONTEXT-MANAGEMENT.md) for details.

## Automated MCP Configuration

The installation uses `switch-mcp-profile.sh` to configure all agents simultaneously.

### What Gets Configured
- **Claude Code**: `.mcp.json` (Project-specific)
- **Gemini CLI**: `.gemini/settings.json` (Project-specific)
- **Codex CLI**: `~/.codex/config.toml` (User-level)

**Important**: These files contain absolute paths and API keys injected from your `.env`. Do not commit them to version control.

## Verification

### Check Installations
```bash
# Verify CLIs
claude --version
gemini --version
codex --version

# Verify Configuration
ls -l .mcp.json .gemini/settings.json
```

### Test Connectivity
```bash
# Start Gemini and list tools
gemini "List available MCP tools"

# Start Claude and check MCP status
claude
/mcp
```

## File Locations

### Configuration Files
```
SciAgent-toolkit/
├── .env                         # API Keys (Git-ignored)
├── .mcp.json                    # Claude Code config
├── .gemini/
│   └── settings.json            # Gemini CLI config
└── tooluniverse-env/            # Python venv for tools
```

### User-Level Configurations
```
~/.codex/config.toml             # Codex CLI config
~/.npm-global/bin/               # CLI binaries (gemini, codex)
```

## Usage Examples

### Claude Code (General Coding)
```bash
claude
"Refactor this function to be more efficient."
"Use PAL to ask Gemini for a second opinion."
```

### Gemini CLI (Deep Research)
```bash
# Switch to research profile first
./scripts/switch-mcp-profile.sh research-lite

gemini "Search PubMed for recent papers on CRISPR off-target effects."
```

### Hybrid Workflow (Orchestration)
```bash
# Switch to hybrid profile
./scripts/switch-mcp-profile.sh hybrid-research

claude
"Write a Python script to analyze these sequences. Use PAL to run Gemini and search for the latest protocols first."
```

## Uninstallation

### Remove CLIs
```bash
rm ~/.local/bin/claude
npm uninstall -g @google/gemini-cli @openai/codex
```

### Remove Toolkit
```bash
rm -rf tooluniverse-env
rm .mcp.json .env
rm -rf .gemini
```

## Troubleshooting

For common issues, see [TROUBLESHOOTING.md](TROUBLESHOOTING.md).

