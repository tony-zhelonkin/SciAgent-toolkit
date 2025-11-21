# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

---

## Repository Overview

**SciAgent-toolkit** is a modular MCP (Model Context Protocol) infrastructure orchestrator that integrates AI assistants (Claude Code, Codex CLI) with specialized scientific research tools. This is an **installation and configuration framework**, not a standalone application.

### Core Architecture

Three-tier system:
1. **User Interfaces**: Claude Code CLI, Codex CLI (optional)
2. **MCP Server Layer**: Sequential Thinking, ToolUniverse (600+ scientific tools), Serena (code intelligence), PubMed (optional plugin)
3. **External Data Sources**: ChEMBL, UniProt, DrugBank, FDA, PubMed, ClinicalTrials.gov, Europe PMC, etc.

---

## Key Directories and Files

### Installation Scripts
- `scripts/setup_mcp_infrastructure.sh` - **Main orchestrator**, runs all installation steps in correct order
- `scripts/install_claude.sh` - Installs Claude Code to `~/.local/bin/claude`
- `scripts/install_codex.sh` - Installs Codex CLI (optional alternative interface)
- `scripts/mcp_servers/setup_serena.sh` - Installs Serena MCP (code intelligence, uvx-based)
- `scripts/mcp_servers/setup_sequential_thinking.sh` - Installs Sequential Thinking MCP (npx-based)
- `scripts/mcp_servers/setup_tooluniverse.sh` - Installs ToolUniverse MCP (600+ scientific tools, Python-based)
- `scripts/mcp_servers/setup_pubmed.sh` - Provides manual PubMed plugin installation instructions
- `scripts/configure_mcp_servers.sh` - Auto-generates `.mcp.json` and `~/.codex/config.toml` configurations
- `scripts/test_installation.sh` - Validates full installation

### Configuration Files
- `.mcp.json` - **Project-local** MCP server configuration for Claude Code (JSON format)
- `~/.codex/config.toml` - **User-global** MCP server configuration for Codex CLI (TOML format)
- `.claude/settings.local.json` - Claude Code settings (pre-approved commands, enabled MCP servers)

### Virtual Environments
- `tooluniverse-env/` - Python virtual environment for ToolUniverse (managed by uv)
- `scripts/.mcp_setup/` - MCP configuration backups and staging area

### Agents
- `agents/` - Canonical location for custom Claude agents
- `.claude/agents/` - Symlinked to `agents/` for auto-discovery by Claude Code
- `agents/bioinf-librarian.md` - Expert agent for finding bioinformatics tools/docs
- `agents/rnaseq-methods-writer.md` - Auto-generates publication Methods sections from RNA-seq code

### Docker Testing Infrastructure
- `docker/test/Dockerfile.tooluniverse-test` - Base image with ToolUniverse
- `docker/test/Dockerfile.claude-test` - Extends base with Claude Code
- `docker/test/Dockerfile.codex-test` - Extends base with Codex CLI
- `docker/test/test-all.sh` - Sequential build and validation of all images

---

## Common Development Commands

### Installation and Setup

```bash
# Full installation (Claude Code + all MCP servers)
./scripts/setup_mcp_infrastructure.sh

# Claude Code only (skip Codex CLI)
./scripts/setup_mcp_infrastructure.sh --skip-codex

# MCP servers only (assumes Claude/Codex already installed)
./scripts/setup_mcp_infrastructure.sh --mcp-only

# Skip specific components
./scripts/setup_mcp_infrastructure.sh --skip-tooluniverse --skip-pubmed

# View all options
./scripts/setup_mcp_infrastructure.sh --help
```

### Configuration Management

```bash
# Regenerate MCP configurations
./scripts/configure_mcp_servers.sh --project-dir "$(pwd)"

# Validate existing configuration
python3 -m json.tool .mcp.json

# List configured MCP servers
jq '.mcpServers | keys' .mcp.json

# Add MCP server manually (Claude Code)
claude mcp add server-name --scope local -- command arg1 arg2
```

### Testing and Verification

```bash
# Verify Claude Code installation
claude --version
claude doctor

# Verify Codex CLI installation (if installed)
codex --version

# Test ToolUniverse installation
./scripts/test_tooluniverse.sh

# Test individual MCP server
npx -y @modelcontextprotocol/server-sequential-thinking --help
uvx --from git+https://github.com/oraios/serena serena --help
uv --directory ./tooluniverse-env run tooluniverse-mcp --version

# Comprehensive installation test
./scripts/test_installation.sh
```

### Docker Testing

```bash
# Build and test all images sequentially
cd docker/test
./test-all.sh

# Build individual images
docker build -f Dockerfile.tooluniverse-test -t scibio-tooluniverse:test .
docker build -f Dockerfile.claude-test -t scibio-claude:test .
docker build -f Dockerfile.codex-test -t scibio-codex:test .
```

---

## Architecture Deep Dive

### Installation Flow

The main orchestrator (`setup_mcp_infrastructure.sh`) follows this sequence:

1. **Parse command-line arguments** (--skip-*, --*-only flags)
2. **Install Claude Code** (if not skipped): Downloads from https://claude.ai/install.sh, installs to `~/.local/bin/claude`, updates PATH
3. **Install Codex CLI** (if not skipped): Attempts npm global install, falls back to Homebrew on macOS
4. **Install Base MCP Servers**:
   - **Serena**: Installs uv/uvx, uses `git+https://github.com/oraios/serena` (HTTPS only, no SSH required)
   - **Sequential Thinking**: Installs Node.js/npm if needed, uses npx with `@modelcontextprotocol/server-sequential-thinking`
5. **Install Scientific MCP Servers**:
   - **PubMed**: Provides manual plugin installation instructions (Claude Code marketplace only)
   - **ToolUniverse**: Installs uv, creates `./tooluniverse-env/`, auto-detects command name (`tooluniverse-mcp` or `tooluniverse-smcp-stdio`)
6. **Generate Configurations**: Runs `configure_mcp_servers.sh` to create `.mcp.json` and `~/.codex/config.toml`
7. **Verification**: Checks installed binaries and validates configuration files

### Configuration Generation Strategy

`configure_mcp_servers.sh` implements **runtime detection**:

1. **Detect Sequential Thinking**: Checks for `npx` availability
2. **Detect ToolUniverse**: Checks for `tooluniverse-env/` directory in both `scripts/` and project root
3. **Detect Serena**: Checks for `uvx` with 5-second timeout (skips if not pre-cached to avoid 5-15 min build time)
4. **Generate `.mcp.json`**: Project-local JSON configuration for Claude Code
5. **Generate `~/.codex/config.toml`**: User-global TOML configuration for Codex CLI
6. **Validate**: Uses `python3 -m json.tool` to check JSON syntax

**Critical**: MCP configurations MUST use **absolute paths** because Claude Code doesn't guarantee working directory.

**Environment-Specific Configuration**: The `.mcp.json` file is environment-specific and contains absolute paths. It is:
- **Generated automatically** during installation
- **Not committed to git** (in `.gitignore`)
- **Must be regenerated** when moving between environments (Dev containers, different machines)

To regenerate for a new environment:
```bash
./scripts/configure_mcp_servers.sh --force
```

### MCP Server Details

#### Sequential Thinking MCP
- **Purpose**: Structured reasoning for complex decisions
- **Command**: `npx -y @modelcontextprotocol/server-sequential-thinking`
- **Requirements**: Node.js 18+, npm/npx
- **Features**: Step-by-step analysis, decision trees, multi-step reasoning chains

#### ToolUniverse MCP
- **Purpose**: 600+ scientific research tools
- **Command**: `uv --directory ./scripts/tooluniverse-env run tooluniverse-mcp` (or `tooluniverse-smcp-stdio`)
- **Requirements**: Python 3.10+, uv package manager
- **Environment**: Local venv in `./scripts/tooluniverse-env/` or `./tooluniverse-env/`
- **Tool Categories**: Drug discovery (ChEMBL, DrugBank, FDA), Genomics (UniProt, protein interactions), Literature (PubMed, Europe PMC, Semantic Scholar), Clinical trials (ClinicalTrials.gov)
- **Important**: Command name varies by version; setup script auto-detects correct name and location

#### Serena MCP
- **Purpose**: Semantic code search and editing
- **Command**: `uvx --from git+https://github.com/oraios/serena serena start-mcp-server`
- **Requirements**: uvx (automatically installed)
- **Build Time**: 5-15 minutes on first run (one-time Rust compilation), seconds for cached builds
- **Features**: Symbol-level analysis, cross-reference tracking, intelligent refactoring

#### PubMed Plugin (Claude Code only)
- **Purpose**: Biomedical literature access (36M+ articles)
- **Installation**: Manual via Claude Code plugin marketplace
- **Commands**:
  ```bash
  claude
  /plugin marketplace add anthropics/life-sciences
  /plugin install pubmed@life-sciences
  # Restart Claude Code
  ```
- **No API Key Required**: Uses public NCBI E-utilities

### Tool Filtering for Context Management

ToolUniverse provides 600+ tools which can overflow Claude's context window. Use filtering:

#### Include Specific Tools
```json
{
  "mcpServers": {
    "tooluniverse": {
      "args": [
        "--directory", "/absolute/path/to/tooluniverse-env",
        "run", "tooluniverse-mcp",
        "--include-tools", "EuropePMC_search_articles,ChEMBL_search_similar_molecules"
      ]
    }
  }
}
```

#### Exclude Tool Categories
```json
{
  "args": [
    "--directory", "/absolute/path/to/tooluniverse-env",
    "run", "tooluniverse-mcp",
    "--exclude-tool-types", "PackageTool"
  ]
}
```

#### Multiple Specialized Instances
```json
{
  "mcpServers": {
    "tooluniverse-literature": {
      "command": "uv",
      "args": ["--directory", "/path/to/tooluniverse-env", "run", "tooluniverse-mcp", "--include-tools", "EuropePMC_search_articles,PubMed_search"]
    },
    "tooluniverse-drugs": {
      "command": "uv",
      "args": ["--directory", "/path/to/tooluniverse-env", "run", "tooluniverse-mcp", "--include-tools", "ChEMBL_search_similar_molecules,FDA_get_drug_info"]
    }
  }
}
```

### Optional: Azure OpenAI Summarization

ToolUniverse can auto-summarize long outputs using Azure OpenAI:

```bash
# Add to ~/.bashrc or ~/.zshrc
export AZURE_OPENAI_API_KEY="your-api-key"
export AZURE_OPENAI_ENDPOINT="https://your-resource.openai.azure.com"

# Re-run ToolUniverse setup to enable summarization
./scripts/mcp_servers/setup_tooluniverse.sh
```

Configuration example:
```json
{
  "mcpServers": {
    "tooluniverse": {
      "env": {
        "AZURE_OPENAI_API_KEY": "${AZURE_OPENAI_API_KEY}",
        "AZURE_OPENAI_ENDPOINT": "${AZURE_OPENAI_ENDPOINT}"
      },
      "args": ["--directory", "/path/to/tooluniverse-env", "run", "tooluniverse-mcp", "--hook-type", "SummarizationHook"]
    }
  }
}
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
Your agent's role and capabilities...

# Methodology
How the agent works...

# Examples
Example use cases...
```

### Creating New Agents

1. Create `agents/new-agent-name.md` following the structure above
2. Claude Code auto-discovers agents from `agents/` directory
3. Test agent by referencing it in prompts
4. Document in `agents/README.md`

### Existing Agents

#### Bioinformatics Research Librarian
- **File**: `agents/bioinf-librarian.md`
- **Purpose**: Find bioinformatics tools, documentation, resources
- **Methodology**: Prioritizes GitHub repos → official docs → peer-reviewed literature
- **Output**: Structured notes in `web_notes.md`

#### RNA-seq Methods Writer
- **File**: `agents/rnaseq-methods-writer.md`
- **Purpose**: Auto-generate publication Methods sections from RNA-seq analysis code
- **Input**: Analysis scripts and notebooks
- **Output**: Formal scientific writing with statistical models

---

## Workflow Patterns

### Literature-Driven Discovery
```
PubMed Search → Sequential Thinking (analyze) → ToolUniverse (validate) → Synthesis
```

**Example**: Search papers on topic → analyze methodologies → check clinical trials → recommend next steps

### Drug Discovery Pipeline
```
ChEMBL → FDA Database → ClinicalTrials.gov → PubMed → Sequential Analysis
```

**Example**: Find similar molecules → check safety profiles → find ongoing trials → review recent literature

### Genomics Research
```
UniProt → Protein Interaction DBs → PubMed → Sequential Analysis
```

**Example**: Get protein info → find interactions → search recent literature → analyze functional implications

---

## Troubleshooting

### MCP Servers Not Loading

**Symptom**: `/mcp` in Claude Code shows "No MCP servers configured"

**Diagnosis**:
```bash
# Check if .mcp.json exists
ls -la .mcp.json

# Validate JSON syntax
python3 -m json.tool .mcp.json

# Check Claude Code diagnostics
claude doctor
```

**Solutions**:
1. Regenerate configuration: `./scripts/configure_mcp_servers.sh`
2. Manually add servers: `claude mcp add server-name --scope local -- command args`
3. Restart Claude Code

### ToolUniverse Command Not Found

**Symptom**: `tooluniverse-mcp: command not found` or similar

**Diagnosis**:
```bash
# Check which command is installed
uv --directory ./tooluniverse-env run tooluniverse-mcp --version
uv --directory ./tooluniverse-env run tooluniverse-smcp-stdio --version
```

**Solution**: Update `.mcp.json` with correct command name detected above

### Serena Installation Takes Too Long

**Symptom**: First Serena run takes 5-15 minutes

**Explanation**: Serena is a Rust project compiled on first run (one-time build)

**Solutions**:
1. **Wait it out**: First build is one-time, subsequent runs are instant
2. **Skip Serena**: Use `--skip-serena` flag (not currently exposed but can modify scripts)
3. **Pre-cache**: Run `uvx --from git+https://github.com/oraios/serena serena --help` before configuration

### No .mcp.json After Installation

**Symptom**: MCP servers installed but not configured

**Solution**:
```bash
# Run configuration script manually
./scripts/configure_mcp_servers.sh --project-dir "$(pwd)"

# Verify it worked
cat .mcp.json
```

### PubMed Not Available

**Symptom**: PubMed doesn't appear in `/mcp` list

**Explanation**: PubMed is a **Claude Code plugin**, not a standalone MCP server

**Solution**:
```bash
claude
/plugin marketplace add anthropics/life-sciences
/plugin install pubmed@life-sciences
# Press Ctrl+D to exit, then restart Claude Code
```

### Git SSH Errors During Serena Install

**Symptom**: `Permission denied (publickey)` or `ssh_askpass: exec failed`

**Explanation**: Fixed in v1.1.0+ by using HTTPS instead of SSH

**Solution**:
```bash
# Update to latest version
git pull origin main

# Force Git to use HTTPS globally (if needed)
git config --global url."https://github.com/".insteadOf git@github.com:

# Re-run Serena setup
./scripts/mcp_servers/setup_serena.sh
```

---

## Important Implementation Details

### Idempotent Installation
- All scripts check for existing installations
- Safe to run multiple times
- Incremental updates supported
- No destructive operations

### Absolute Paths Required
- **Critical**: `.mcp.json` MUST use absolute paths
- Reason: Claude Code doesn't guarantee working directory
- Generated configs automatically use absolute paths

### Version-Dependent Command Names
- ToolUniverse command name varies: `tooluniverse-mcp` or `tooluniverse-smcp-stdio`
- Setup script auto-detects and uses correct name
- Configuration generation handles both variants

### Virtual Environment Management
- ToolUniverse creates isolated venv in `./tooluniverse-env/`
- Managed by uv package manager (faster than pip/conda)
- **Critical**: Create venv BEFORE installing packages (fixed in v1.1.0+)

### PackageTool Exclusion
- ToolUniverse's PackageTool excluded by default: `--exclude-tool-types PackageTool`
- Reason: Not useful for research workflows, reduces context usage

### Platform Differences
- **macOS**: Claude Code installed to `~/.local/bin/claude`, Codex CLI via Homebrew or npm
- **Linux**: Claude Code installed to `~/.local/bin/claude`, Codex CLI via npm (may require sudo)
- **Windows**: Use WSL (Windows Subsystem for Linux)

---

## File Locations Reference

### System-Level Installations
```
~/.local/bin/claude              # Claude Code binary
~/.cargo/bin/uv                  # UV package manager
/usr/local/bin/codex             # Codex CLI (Homebrew on macOS)
~/.npm-global/bin/codex          # Codex CLI (npm)
```

### Project-Level Files
```
SciAgent-toolkit/
├── .mcp.json                    # Claude Code MCP configuration (project-local)
├── tooluniverse-env/            # ToolUniverse Python environment
│   └── .venv/                   # Managed by uv
├── test_tooluniverse.sh         # Quick test script (auto-generated)
├── agents/                      # Canonical agent definitions
├── .claude/
│   ├── settings.local.json      # Claude Code settings
│   └── agents/                  # Symlinked to agents/
└── scripts/
    ├── .mcp_setup/              # MCP configuration staging
    └── mcp_servers/             # Individual MCP setup scripts
```

### User-Level Configurations
```
~/.bashrc                        # Updated with PATH additions
~/.codex/config.toml             # Codex CLI MCP configuration (user-global)
```

---

## Adding New MCP Servers

### Step 1: Create Setup Script

Create `scripts/mcp_servers/setup_newserver.sh`:

```bash
#!/usr/bin/env bash
set -euo pipefail

echo "Installing NewServer MCP..."

# Check prerequisites
if ! command -v required_command &>/dev/null; then
    echo "Installing required_command..."
    # Installation logic
fi

# Install MCP server
echo "Setting up NewServer..."
# Installation commands here

# Test installation
if command -v newserver &>/dev/null; then
    echo "NewServer MCP installed successfully"
else
    echo "NewServer MCP installation failed"
    exit 1
fi
```

### Step 2: Add to Main Orchestrator

Edit `scripts/setup_mcp_infrastructure.sh`:

```bash
# Add flag for new server
INSTALL_NEWSERVER=true

# Add to installation section
if [ "$INSTALL_NEWSERVER" = true ]; then
    separator "Installing NewServer MCP"
    bash "${SCRIPT_DIR}/mcp_servers/setup_newserver.sh" || {
        log_error "NewServer MCP setup failed"
        exit 1
    }
fi
```

### Step 3: Update Configuration Generator

Edit `scripts/configure_mcp_servers.sh`:

```bash
# Add detection logic
if command -v newserver &>/dev/null; then
    NEWSERVER_AVAILABLE=true
fi

# Add to JSON generation
if [ "$NEWSERVER_AVAILABLE" = true ]; then
    cat >> .mcp.json <<EOF
    "newserver": {
      "command": "newserver",
      "args": ["start"]
    }
EOF
fi
```

### Step 4: Document

- Add to `docs/INSTALLATION.md`
- Update `README.md` features list
- Add usage examples to `docs/QUICKSTART.md`

---

## Contributing Changes

### Before Making Changes

1. Read `CONTRIBUTING.md` for guidelines
2. Test changes on clean environment (use Docker)
3. Validate all shell scripts with `shellcheck`
4. Update documentation for new features

### Testing Changes

```bash
# Test individual components
./scripts/mcp_servers/setup_serena.sh

# Test full installation
./scripts/setup_mcp_infrastructure.sh

# Test configuration generation
./scripts/configure_mcp_servers.sh

# Validate JSON
python3 -m json.tool .mcp.json

# Run comprehensive tests
./scripts/test_installation.sh

# Test in Docker
cd docker/test
./test-all.sh
```

### Commit Message Format

```bash
# Good examples
git commit -m "Add support for custom ToolUniverse filters"
git commit -m "Fix PATH configuration in install_claude.sh"
git commit -m "Update INSTALLATION.md with macOS M1 instructions"

# Use imperative mood
# Keep first line under 72 characters
# Add detailed description if needed
```

---

## Quick Reference

### Check Installation Status
```bash
claude --version                 # Claude Code
codex --version                  # Codex CLI
cat .mcp.json                    # MCP configuration
claude doctor                    # Diagnostics
```

### Rebuild Configuration
```bash
./scripts/configure_mcp_servers.sh --project-dir "$(pwd)"
```

### Test MCP Servers
```bash
./scripts/test_tooluniverse.sh
./scripts/test_installation.sh
```

### Launch Interfaces
```bash
claude                           # Start Claude Code
codex                            # Start Codex CLI
```

### Common Claude Code Commands
```bash
/mcp                             # List MCP servers
/plugin marketplace list         # List available plugins
/plugin list                     # List installed plugins
```

---

## Resources

- [Claude Code Documentation](https://docs.anthropic.com/en/docs/claude-code)
- [MCP Protocol](https://modelcontextprotocol.io/)
- [ToolUniverse GitHub](https://github.com/greedyai/tool-universe)
- [Serena GitHub](https://github.com/oraios/serena)
- [Sequential Thinking MCP](https://github.com/modelcontextprotocol/servers/tree/main/src/sequentialthinking)
- [PubMed E-utilities](https://www.ncbi.nlm.nih.gov/books/NBK25501/)
