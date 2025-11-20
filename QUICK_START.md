# Quick Start Guide - MCP Infrastructure

## TL;DR - Get Started in 5 Minutes

### Step 1: Run the Installer (2 minutes)

```bash
cd /data1/users/antonz/projects/GVDRP1_prj

# Full installation (recommended)
./setup_mcp_infrastructure.sh

# OR selective installation
./setup_mcp_infrastructure.sh --skip-codex  # Claude only
```

### Step 2: Configure PubMed (1 minute, Claude Code only)

```bash
claude
/plugin marketplace add anthropics/life-sciences
/plugin install pubmed@life-sciences
# Press Ctrl+D to exit, then restart
```

### Step 3: Verify (1 minute)

```bash
claude
/mcp

# You should see:
# - serena
# - sequential-thinking
# - tooluniverse
# - pubmed (after plugin install)
```

### Step 4: Try It Out (1 minute)

```bash
# In Claude Code, try:
"What scientific tools are available from ToolUniverse?"
"Find recent studies about cancer immunotherapy on PubMed"
"Search ChEMBL for molecules similar to aspirin"
```

Done! You now have a complete AI scientist agent.

---

## Installation Options Reference

| Command | Installs | Use Case |
|---------|----------|----------|
| `./setup_mcp_infrastructure.sh` | Everything | Full setup |
| `./setup_mcp_infrastructure.sh --skip-codex` | Claude + MCP | Claude only |
| `./setup_mcp_infrastructure.sh --skip-claude` | Codex + MCP | Codex only |
| `./setup_mcp_infrastructure.sh --mcp-only` | MCP servers only | Add to existing |
| `./mcp_servers/setup_tooluniverse.sh` | ToolUniverse only | Scientific tools |

---

## Quick Command Reference

### Installation
```bash
./setup_mcp_infrastructure.sh              # Full install
./setup_mcp_infrastructure.sh --help       # Show options
```

### Verification
```bash
claude --version                           # Check Claude
codex --version                            # Check Codex
claude doctor                              # Run diagnostics
python3 -m json.tool .mcp.json            # Validate config
```

### Testing
```bash
./test_tooluniverse.sh                     # Test ToolUniverse
```

### Usage
```bash
claude                                     # Start Claude Code
codex                                      # Start Codex CLI
```

---

## Common First Queries

### Check Available Tools
```
"What MCP servers are available?"
"List all ToolUniverse scientific tools"
"/mcp"  # Show MCP server status
```

### Test Scientific Capabilities

**Literature Search:**
```
"Find recent papers about CRISPR gene editing on PubMed"
"Search Europe PMC for immunotherapy studies from 2024"
```

**Drug Discovery:**
```
"Search ChEMBL for molecules similar to ibuprofen"
"Find FDA approval information for metformin"
```

**Genomics:**
```
"Get protein information for TP53 from UniProt"
"Find protein-protein interactions for BRCA1"
```

**Clinical Research:**
```
"Search ClinicalTrials.gov for diabetes treatments"
"Find recent Phase III trials for cancer drugs"
```

---

## File Structure at a Glance

```
GVDRP1_prj/
├── QUICK_START.md                ← You are here
├── MCP_SETUP_README.md           ← Full user guide
├── ARCHITECTURE.md               ← System architecture
├── INSTALLATION_SUMMARY.md       ← What was installed
│
├── setup_mcp_infrastructure.sh   ← MAIN INSTALLER
├── install_claude.sh
├── install_codex.sh
│
├── mcp_servers/
│   ├── setup_serena.sh
│   ├── setup_sequential_thinking.sh
│   ├── setup_tooluniverse.sh
│   └── setup_pubmed.sh
│
├── config/
│   ├── README.md                 ← Configuration guide
│   └── merge_mcp_configs.sh
│
├── .mcp.json                     ← Auto-generated
├── tooluniverse-env/             ← Auto-generated
└── test_tooluniverse.sh          ← Auto-generated
```

---

## What Gets Installed

### Required (Always)
- ✅ Claude Code CLI
- ✅ Serena MCP (code intelligence)
- ✅ Sequential Thinking MCP (reasoning)
- ✅ uv/uvx (for Serena)
- ✅ Node.js/npm (for Sequential)

### Scientific Research (Recommended)
- ✅ ToolUniverse MCP (600+ tools)
- ✅ PubMed Plugin (36M+ articles)

### Optional
- ⚪ Codex CLI (alternative interface)
- ⚪ Azure OpenAI (for summarization)

---

## Troubleshooting Cheat Sheet

### Problem: Command not found after install
```bash
source ~/.bashrc
# OR restart terminal
```

### Problem: MCP servers not loading
```bash
claude doctor
python3 -m json.tool .mcp.json
```

### Problem: ToolUniverse not working
```bash
./test_tooluniverse.sh
uv --version
python3 --version  # Should be 3.10+
```

### Problem: PubMed not available
```bash
# PubMed requires manual plugin install
claude
/plugin marketplace add anthropics/life-sciences
/plugin install pubmed@life-sciences
```

---

## Next Steps After Quick Start

1. **Read the docs**: `MCP_SETUP_README.md`
2. **Explore configuration**: `config/README.md`
3. **Understand architecture**: `ARCHITECTURE.md`
4. **Try example workflows**: See MCP_SETUP_README.md
5. **Customize**: Filter tools, add instances

---

## Scientific Research Workflow Examples

### Example 1: Drug Discovery (5 minutes)
```
1. "Find molecules similar to aspirin in ChEMBL"
2. "Check FDA approval status for similar NSAIDs"
3. "Search clinical trials for pain management"
4. "Find recent publications about COX-2 inhibitors on PubMed"
5. "Summarize findings and safety profiles"
```

### Example 2: Genomics Research (5 minutes)
```
1. "Get protein information for BRCA1 from UniProt"
2. "Find protein-protein interactions for BRCA1"
3. "Search PubMed for BRCA1 mutations"
4. "Identify pathways involving BRCA1"
5. "Analyze functional implications"
```

### Example 3: Literature Review (5 minutes)
```
1. "Search PubMed for immunotherapy meta-analyses"
2. "Get full-text articles from PMC for top 3 results"
3. "Extract key findings and methodologies"
4. "Find related articles in Semantic Scholar"
5. "Identify research gaps"
```

---

## Configuration Quick Tips

### Reduce Context Usage
Edit `.mcp.json` to filter ToolUniverse tools:
```json
{
  "mcpServers": {
    "tooluniverse": {
      "args": [
        "--directory", "/path/to/tooluniverse-env",
        "run", "tooluniverse-smcp-stdio",
        "--include-tools", "EuropePMC_search_articles,ChEMBL_search_similar_molecules"
      ]
    }
  }
}
```

### Enable Summarization
```bash
# Add to ~/.bashrc
export AZURE_OPENAI_API_KEY="your-key"
export AZURE_OPENAI_ENDPOINT="https://your-resource.openai.azure.com"

# Re-run ToolUniverse setup
./mcp_servers/setup_tooluniverse.sh
```

### Multiple Tool Sets
Create specialized instances:
- `tooluniverse-lit`: Literature research only
- `tooluniverse-drug`: Drug discovery only
- `tooluniverse-genomics`: Genomics only

See `config/README.md` for examples.

---

## Support & Documentation

| Question | Resource |
|----------|----------|
| How do I install? | This file (QUICK_START.md) |
| What can I do? | MCP_SETUP_README.md |
| How do I configure? | config/README.md |
| How does it work? | ARCHITECTURE.md |
| What was installed? | INSTALLATION_SUMMARY.md |

---

## Remember

1. **PubMed requires manual plugin install** (Claude Code only)
2. **ToolUniverse has 600+ tools** - filter them for better performance
3. **Both Claude and Codex are supported** - use either or both
4. **All scripts are idempotent** - safe to run multiple times
5. **Check /mcp in Claude/Codex** to verify servers are loaded

---

**Happy researching! 🔬🧬💊**

For detailed documentation, see: `MCP_SETUP_README.md`
