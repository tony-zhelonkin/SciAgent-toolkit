# Building AI Scientists with MCP Infrastructure

Complete setup and configuration guide for creating customized AI scientist agents using Claude Code, Codex CLI, and scientific MCP servers.

## Quick Start

### Full Installation (Recommended)

```bash
# Install everything: Claude Code, Codex CLI, and all MCP servers
./setup_mcp_infrastructure.sh
```

### Selective Installation

```bash
# Install only Claude Code with MCP servers (skip Codex)
./setup_mcp_infrastructure.sh --skip-codex

# Install only Codex CLI with MCP servers (skip Claude)
./setup_mcp_infrastructure.sh --skip-claude

# Install only MCP servers (assumes Claude/Codex already installed)
./setup_mcp_infrastructure.sh --mcp-only
```

## What Gets Installed

### 1. Claude Code (Optional)
- **Purpose**: IDE/CLI interface with advanced reasoning
- **Installation**: Native installer via https://claude.ai/install.sh
- **Location**: `~/.local/bin/claude`

### 2. Codex CLI (Optional)
- **Purpose**: Terminal-based scientific research interface
- **Installation**: npm or Homebrew
- **Location**: Global npm package or Homebrew installation

### 3. MCP Servers

#### Serena (Automatic)
- **Purpose**: Semantic code search and editing
- **Requirements**: uvx (automatically installed)
- **Installation**: Via uvx from GitHub

#### Sequential Thinking (Automatic)
- **Purpose**: Structured reasoning for complex decisions
- **Requirements**: npx (automatically installed via Node.js)
- **Installation**: Via npx from npm

#### ToolUniverse (Automatic)
- **Purpose**: 600+ scientific research tools
- **Requirements**: Python 3.10+, uv
- **Installation**: Local environment in `./tooluniverse-env/`
- **Features**:
  - Drug discovery and development
  - Genomics and molecular biology
  - Literature research (PubMed, Semantic Scholar, Europe PMC)
  - Clinical trials (ClinicalTrials.gov)
  - FDA approvals and safety data
  - Protein analysis (UniProt, ChEMBL)

#### PubMed (Manual Setup Required)
- **Purpose**: Biomedical literature access (36+ million articles)
- **Requirements**: Claude Code
- **Installation**: Via Claude Code plugin marketplace (manual)
- **Setup Instructions**: See below

## Installation Instructions

### Prerequisites

Your system should have:
- Linux or macOS (WSL supported on Windows)
- Python 3.10 or later
- Internet connection
- sudo access (for some system package installations on Linux)

### Step 1: Run the Main Installer

```bash
# Make the script executable (if needed)
chmod +x setup_mcp_infrastructure.sh

# Run the installer
./setup_mcp_infrastructure.sh
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

### Step 2: Configure PubMed (Claude Code only)

PubMed must be installed manually via the Claude Code plugin system:

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

### Step 3: Optional - Set Up Azure OpenAI for Summarization

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
./mcp_servers/setup_tooluniverse.sh
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

## Advanced Configuration

### Tool Filtering for Better Performance

ToolUniverse provides 600+ tools, which can overwhelm the context window. You can filter tools to only include what you need:

#### For Claude Code (edit `.mcp.json`):

```json
{
  "mcpServers": {
    "tooluniverse-research": {
      "type": "stdio",
      "command": "uv",
      "args": [
        "--directory",
        "/path/to/GVDRP1_prj/tooluniverse-env",
        "run",
        "tooluniverse-smcp-stdio",
        "--include-tools",
        "EuropePMC_search_articles,ChEMBL_search_similar_molecules,search_clinical_trials"
      ]
    }
  }
}
```

#### For Codex CLI (edit `~/.codex/config.toml`):

```toml
[mcp_servers.tooluniverse-research]
command = "uv"
args = [
  "--directory",
  "/path/to/GVDRP1_prj/tooluniverse-env",
  "run",
  "tooluniverse-smcp-stdio",
  "--include-tools",
  "EuropePMC_search_articles,ChEMBL_search_similar_molecules,search_clinical_trials"
]
```

### Multiple ToolUniverse Instances

You can run multiple ToolUniverse instances with different tool sets:

- `tooluniverse-research`: Literature research tools only
- `tooluniverse-analysis`: Drug discovery and genomics tools
- `tooluniverse-clinical`: Clinical trials and FDA tools

See `config/README.md` for detailed configuration examples.

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

## Troubleshooting

### Issue: Claude Code or Codex not found after installation

**Solution:**
```bash
# Reload your shell configuration
source ~/.bashrc  # or source ~/.zshrc

# Or restart your terminal
```

### Issue: MCP servers not loading

**Solution:**
```bash
# For Claude Code
claude doctor

# Check MCP configuration
python3 -m json.tool .mcp.json

# Test individual MCP servers
uvx --from git+https://github.com/oraios/serena serena --help
npx -y @modelcontextprotocol/server-sequential-thinking --help
uv --directory ./tooluniverse-env run tooluniverse-smcp-stdio --help
```

### Issue: ToolUniverse tools not executing

**Solution:**
- Verify Python version: `python3 --version` (should be 3.10+)
- Check uv installation: `uv --version`
- Test ToolUniverse: `./test_tooluniverse.sh`
- Check API keys if using summarization features

### Issue: PubMed not available in Claude Code

**Solution:**
PubMed is installed separately via the plugin system:
```bash
claude
/plugin marketplace add anthropics/life-sciences
/plugin install pubmed@life-sciences
# Restart Claude Code
```

### Issue: Permission denied errors

**Solution:**
```bash
# Make scripts executable
chmod +x setup_mcp_infrastructure.sh
chmod +x install_claude.sh
chmod +x install_codex.sh
chmod +x mcp_servers/*.sh
chmod +x config/*.sh
```

## File Structure

```
GVDRP1_prj/
├── setup_mcp_infrastructure.sh      # Main orchestrator (START HERE)
├── install_claude.sh                # Claude Code installer
├── install_codex.sh                 # Codex CLI installer
├── test_tooluniverse.sh             # ToolUniverse test script (auto-generated)
├── .mcp.json                        # Claude Code MCP config (auto-generated)
├── tooluniverse-env/                # ToolUniverse installation (auto-generated)
│   └── .venv/                       # Python virtual environment
├── mcp_servers/
│   ├── setup_serena.sh              # Serena MCP installer
│   ├── setup_sequential_thinking.sh # Sequential thinking MCP installer
│   ├── setup_tooluniverse.sh        # ToolUniverse MCP installer
│   └── setup_pubmed.sh              # PubMed setup instructions
├── config/
│   ├── README.md                    # Detailed configuration guide
│   └── merge_mcp_configs.sh         # Config merger utility
└── MCP_SETUP_README.md              # This file
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

1. **Learn the Tools**: Explore available tools via `/mcp` command
2. **Read Documentation**: Check `config/README.md` for advanced configuration
3. **Try Examples**: Run the example workflows above
4. **Customize**: Create tool-specific ToolUniverse instances
5. **Integrate**: Build custom research pipelines

## Resources

- [Claude Code Documentation](https://docs.anthropic.com/en/docs/claude-code)
- [MCP Protocol](https://modelcontextprotocol.io/)
- [ToolUniverse GitHub](https://github.com/greedyai/tool-universe)
- [PubMed E-utilities](https://www.ncbi.nlm.nih.gov/books/NBK25501/)
- [Serena GitHub](https://github.com/oraios/serena)

## Support

For issues or questions:
1. Check the troubleshooting section above
2. Review `config/README.md` for detailed configuration
3. Run diagnostics: `claude doctor` or `DEBUG=true codex`
4. Check individual setup script logs

---

**Happy researching with your AI scientist agents!**
