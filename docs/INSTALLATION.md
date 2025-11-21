# SciAgent Toolkit - Complete Installation Guide

Complete setup and configuration guide for creating customized AI scientist agents using Claude Code, Codex CLI, and scientific MCP servers.

## Prerequisites

Your system should have:
- Linux or macOS (WSL supported on Windows)
- Python 3.10 or later
- Internet connection
- sudo access (for some system package installations on Linux)

## Installation Overview

The toolkit will install:
1. **Claude Code** (optional) - IDE/CLI interface with advanced reasoning
2. **Codex CLI** (optional) - Terminal-based scientific research interface
3. **Serena MCP Server** - Semantic code search and editing
4. **Sequential Thinking MCP Server** - Structured reasoning for complex decisions
5. **ToolUniverse MCP Server** - 600+ scientific research tools
6. **PubMed Plugin** - Biomedical literature access (manual setup)

## Full Installation (Recommended)

```bash
# Make the script executable (if needed)
chmod +x scripts/setup_mcp_infrastructure.sh

# Run the installer
./scripts/setup_mcp_infrastructure.sh
```

The installer will:
1. Install Claude Code (if not present)
2. Install Codex CLI (if not present)
3. Install required dependencies (uv, Node.js)
4. Set up Serena MCP server
5. Set up Sequential Thinking MCP server
6. Set up ToolUniverse MCP server
7. Create configuration files
8. Run verification checks

## Selective Installation Options

### Claude Code Only
```bash
./scripts/setup_mcp_infrastructure.sh --skip-codex
```

### Codex CLI Only
```bash
./scripts/setup_mcp_infrastructure.sh --skip-claude
```

### MCP Servers Only
```bash
# Assumes Claude/Codex already installed
./scripts/setup_mcp_infrastructure.sh --mcp-only
```

### Individual Components

Install components separately:

```bash
# Install only Claude Code
./scripts/install_claude.sh

# Install only Codex CLI
./scripts/install_codex.sh

# Install only Serena MCP
./scripts/mcp_servers/setup_serena.sh

# Install only Sequential Thinking MCP
./scripts/mcp_servers/setup_sequential_thinking.sh

# Install only ToolUniverse MCP
./scripts/mcp_servers/setup_tooluniverse.sh
```

## Component Details

### 1. Claude Code

**Purpose**: IDE/CLI interface with advanced reasoning

**Installation**: Native installer via https://claude.ai/install.sh

**Location**: `~/.local/bin/claude`

**Verification**:
```bash
claude --version
claude doctor
```

### 2. Codex CLI

**Purpose**: Terminal-based scientific research interface

**Installation**: npm or Homebrew

**Location**: Global npm package or Homebrew installation

**Authentication**: Sign in with ChatGPT or OpenAI API key

**Verification**:
```bash
codex --version
```

### 3. Serena MCP Server

**Purpose**: Semantic code search and editing

**Requirements**: uvx (automatically installed)

**Installation**: Via uvx from GitHub

**Features**:
- Symbol-level code analysis
- Cross-reference tracking
- Intelligent refactoring
- Context-aware editing

### 4. Sequential Thinking MCP Server

**Purpose**: Structured reasoning for complex decisions

**Requirements**: npx (automatically installed via Node.js)

**Installation**: Via npx from npm

**Features**:
- Step-by-step analysis
- Decision trees
- Complex reasoning chains
- Multi-step technical decisions

### 5. ToolUniverse MCP Server

**Purpose**: 600+ scientific research tools

**Requirements**: Python 3.10+, uv

**Installation**: Local environment in `./tooluniverse-env/`

**Features**:
- **Drug Discovery**: ChEMBL, DrugBank, FDA databases
- **Genomics**: UniProt, protein interactions, pathways
- **Literature**: PubMed, Europe PMC, Semantic Scholar, OpenAlex
- **Clinical Trials**: ClinicalTrials.gov
- **FDA Data**: Approvals, safety data, adverse events
- **Molecular Biology**: Protein structures, gene expression

### 6. PubMed Plugin

**Purpose**: Biomedical literature access (36+ million articles)

**Requirements**: Claude Code

**Installation**: Via Claude Code plugin marketplace (manual)

**Manual Setup Required**:
```bash
# Start Claude Code
claude

# In the Claude Code interface, run:
/plugin marketplace add anthropics/life-sciences
/plugin install pubmed@life-sciences

# Restart Claude Code
# Press Ctrl+D to exit, then run 'claude' again

# Verify installation
/mcp
```

**Features**:
- Search 36M+ PubMed citations
- Access 8M+ PMC full-text articles
- Article metadata and citations
- Related article discovery
- ID format conversion (PMID, PMC ID, DOI)

## Post-Installation Configuration

### Optional: Azure OpenAI for Summarization

If you want ToolUniverse to automatically summarize long outputs:

```bash
# Add to your ~/.bashrc or ~/.zshrc
export AZURE_OPENAI_API_KEY="your-api-key-here"
export AZURE_OPENAI_ENDPOINT="https://your-resource.openai.azure.com"

# Reload your shell configuration
source ~/.bashrc  # or source ~/.zshrc
```

Then re-run the ToolUniverse setup to enable summarization:

```bash
./scripts/mcp_servers/setup_tooluniverse.sh
```

## Automated MCP Configuration

The installation script automatically creates a `.mcp.json` configuration file for Claude Code with all installed MCP servers.

### What Gets Configured Automatically

The following servers are automatically detected and configured:

1. **Sequential Thinking** - If Node.js/npx is available
2. **ToolUniverse** - If `tooluniverse-env/` directory exists
3. **Serena** - If uvx is available

### Configuration File Location

```bash
# Project-specific configuration (Claude Code)
SciAgent-toolkit/.mcp.json

# User-level configuration (Codex CLI)
~/.codex/config.toml
```

**Important**: The `.mcp.json` file is **environment-specific** and contains absolute paths. It is automatically generated and should not be committed to version control (already in `.gitignore`).

### Dev Container and Multi-Environment Setup

When using Dev containers or multiple environments, regenerate the configuration after installation:

```bash
# After setting up in a new environment (Dev container, different machine, etc.)
./scripts/configure_mcp_servers.sh --force
```

This ensures paths are correct for your current environment. The ToolUniverse installation can be in either:
- `scripts/tooluniverse-env/` (default location)
- `tooluniverse-env/` (project root)

The configuration script automatically detects the correct location.

### Manual Configuration (If Needed)

If automatic configuration fails, you can manually configure MCP servers:

```bash
# Run the configuration script separately
./scripts/configure_mcp_servers.sh

# Or manually add servers via Claude Code
claude mcp add sequential-thinking --scope local -- \
  npx -y @modelcontextprotocol/server-sequential-thinking

claude mcp add tooluniverse --scope local -- \
  uv --directory "$(pwd)/tooluniverse-env" run tooluniverse-mcp

claude mcp add serena --scope local -- \
  uvx --from git+https://github.com/oraios/serena serena start-mcp-server
```

### Verifying Configuration

```bash
# Check configuration file exists
ls -la .mcp.json

# Validate JSON syntax
python3 -m json.tool .mcp.json

# List configured servers
jq '.mcpServers | keys' .mcp.json

# Start Claude and verify
claude
/mcp  # Should show all configured servers
```

## Verification

### Check Claude Code Installation

```bash
# Verify Claude Code is installed
claude --version

# Run diagnostics
claude doctor

# Start Claude Code and check MCP servers
claude
# In Claude: /mcp
```

### Check Codex CLI Installation

```bash
# Verify Codex CLI is installed
codex --version

# Start Codex and check MCP servers
codex
# In Codex: /mcp
```

### Check MCP Configuration Files

```bash
# For Claude Code (project-specific)
cat .mcp.json

# For Codex CLI (user-level)
cat ~/.codex/config.toml
```

### Test ToolUniverse

```bash
# Auto-generated test script
./test_tooluniverse.sh
```

## File Locations

### System-Level Installations

```
~/.local/bin/claude              # Claude Code binary
~/.cargo/bin/uv                  # UV package manager
/usr/local/bin/codex             # Codex CLI (if using Homebrew)
# OR
<npm-global>/codex               # Codex CLI (if using npm)
/usr/bin/node                    # Node.js (Linux)
/usr/local/bin/node              # Node.js (macOS via Homebrew)
```

### Project-Level Installations

```
SciAgent-toolkit/
├── .mcp.json                    # Claude Code MCP configuration
├── tooluniverse-env/            # ToolUniverse Python environment
│   └── .venv/                   # Managed by uv
└── test_tooluniverse.sh         # Quick test script
```

### User-Level Configurations

```
~/.bashrc                        # Updated with PATH additions
~/.codex/config.toml             # Codex CLI MCP configuration
```

## Usage Examples

### Claude Code Examples

```bash
# Start Claude Code
claude

# Example scientific queries:
"What scientific tools are available from ToolUniverse?"
"Find recent studies about immunotherapy for melanoma using PubMed"
"Search ChEMBL for molecules similar to aspirin"
"Get protein information for BRCA1 from UniProt"
"Find clinical trials for diabetes treatment"
"Search Europe PMC for CRISPR gene editing papers from the last year"
```

### Codex CLI Examples

```bash
# Start Codex CLI
codex

# Sign in (first time only)
# Select "Sign in with ChatGPT"

# Example scientific queries:
"List all available ToolUniverse tools"
"Find FDA approval information for metformin"
"Search for protein-protein interactions for TP53"
"Get clinical trial data for Alzheimer's disease treatments"
```

## Common Research Workflows

### 1. Drug Discovery Workflow

```
Step 1: "Find similar molecules to [compound name] in ChEMBL"
Step 2: "Search for clinical trials involving [similar compounds]"
Step 3: "Get FDA safety information for [compounds]"
Step 4: "Search PubMed for recent studies on [compounds]"
Step 5: "Summarize findings and identify promising candidates"
```

### 2. Genomics Research Workflow

```
Step 1: "Get protein information for [gene] from UniProt"
Step 2: "Find protein-protein interactions for [gene]"
Step 3: "Identify pathways involving [gene]"
Step 4: "Search literature for [gene] mutations using PubMed"
Step 5: "Analyze functional implications"
```

### 3. Literature Review Workflow

```
Step 1: "Search PubMed for recent papers on [topic]"
Step 2: "Get full-text articles from PMC for top results"
Step 3: "Extract key findings and methodologies"
Step 4: "Find related articles in Semantic Scholar"
Step 5: "Identify research gaps and opportunities"
```

## Uninstallation

### Remove Claude Code

```bash
rm ~/.local/bin/claude
# Remove from PATH in ~/.bashrc if added
```

### Remove Codex CLI

```bash
# If installed via npm
npm uninstall -g @openai/codex

# If installed via Homebrew
brew uninstall codex
```

### Remove MCP Servers

```bash
# Remove ToolUniverse environment
rm -rf ./tooluniverse-env

# Remove configuration files
rm .mcp.json
rm ~/.codex/config.toml

# Remove test scripts
rm test_tooluniverse.sh
```

## Next Steps

1. **Configure PubMed**: Follow manual setup instructions above
2. **Learn the Tools**: Explore available tools via `/mcp` command
3. **Read Advanced Configuration**: Check [CONFIGURATION.md](CONFIGURATION.md)
4. **Try Examples**: Run the example workflows above
5. **Customize**: Create tool-specific ToolUniverse instances

## Resources

- [Claude Code Documentation](https://docs.anthropic.com/en/docs/claude-code)
- [MCP Protocol](https://modelcontextprotocol.io/)
- [ToolUniverse GitHub](https://github.com/greedyai/tool-universe)
- [PubMed E-utilities](https://www.ncbi.nlm.nih.gov/books/NBK25501/)
- [Serena GitHub](https://github.com/oraios/serena)

## Common Installation Issues

### Issue: ToolUniverse Installation Fails

**Error:** `No virtual environment found; run 'uv venv' to create an environment`

**Cause:** Older versions of the script (before v1.1.0) had a bug where the virtual environment wasn't created before installation.

**Solution:** Update to the latest version:
```bash
git pull origin main
./scripts/mcp_servers/setup_tooluniverse.sh
```

### Issue: Serena Installation Fails with SSH Error

**Error:** `Permission denied (publickey)` or `ssh_askpass: exec failed`

**Cause:** Git attempting to use SSH instead of HTTPS for GitHub access.

**Solution:** The latest version (v1.1.0+) automatically uses HTTPS. Update your installation:
```bash
git pull origin main
./scripts/mcp_servers/setup_serena.sh
```

**Manual workaround (if needed):**
```bash
# Force Git to use HTTPS
git config --global url."https://github.com/".insteadOf git@github.com:

# Test Serena installation
uvx --from git+https://github.com/oraios/serena serena --help
```

### Issue: No .mcp.json Created

**Error:** MCP servers installed but not configured, `/mcp` shows "No MCP servers configured"

**Cause:** Configuration not automated in versions before v1.1.0.

**Solution:** Run the configuration script:
```bash
./scripts/configure_mcp_servers.sh
```

### Issue: ToolUniverse Command Not Found

**Error:** `tooluniverse-mcp: command not found` or similar

**Cause:** The command name varies between ToolUniverse versions.

**Solution:** The latest script auto-detects the correct command. Check which is installed:
```bash
uv --directory ./tooluniverse-env run tooluniverse-mcp --version
# OR
uv --directory ./tooluniverse-env run tooluniverse-smcp-stdio --version
```

Update `.mcp.json` with the correct command name if needed.

### Issue: Claude Code Shows "MCP Server Failed to Start"

**Solution:** Test each server independently:
```bash
# Test Sequential Thinking
npx -y @modelcontextprotocol/server-sequential-thinking --help

# Test ToolUniverse
./test_tooluniverse.sh

# Test Serena
uvx --from git+https://github.com/oraios/serena serena --help

# Check Claude logs
claude doctor
```

## Additional Troubleshooting

For more common issues and solutions, see [TROUBLESHOOTING.md](TROUBLESHOOTING.md).

For frequently asked questions, see [FAQ.md](FAQ.md).

---

Happy researching with your AI scientist agents!
