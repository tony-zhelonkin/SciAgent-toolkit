# SciAgent Toolkit

A modular infrastructure for building AI-powered scientific research agents using the Model Context Protocol (MCP). This toolkit enables integration of powerful AI assistants (Claude Code, Gemini CLI, Codex CLI) with specialized scientific tools for drug discovery, genomics, literature research, and more.

## Features

- **600+ Scientific Tools** via ToolUniverse MCP server (ChEMBL, UniProt, DrugBank, FDA databases)
- **36M+ Biomedical Articles** via PubMed integration
- **Code Intelligence** with Serena MCP server
- **Structured Reasoning** via Sequential Thinking MCP server
- **Multi-Agent Support** for Claude Code, Gemini CLI, and Codex CLI
- **Smart Context Management** via profile switching (Coding vs Research modes)
- **Modular Installation** with selective component installation
- **Production-Ready** with comprehensive error handling and validation

## Quick Start

```bash
# Clone the repository
git clone https://github.com/[your-username]/SciAgent-toolkit
cd SciAgent-toolkit

# 1. Setup Environment
cp templates/.env.template .env
# Edit .env and add your API keys (GEMINI_API_KEY, OPENAI_API_KEY, etc.)

# 2. Run the automated installer
./scripts/setup_mcp_infrastructure.sh

# 3. Switch to a profile (e.g., for coding)
./scripts/switch-mcp-profile.sh coding

# Verify installation
claude
/mcp  # Should show: sequential-thinking, pal, context7
```

## Configuration & Secrets

SciAgent-toolkit uses a `.env` file for secure API key management. A `.env.template` is provided in the `templates/` directory.

1.  **Single Source of Truth**: The `.env` file (git-ignored) holds all your API keys.
2.  **Auto-Injection**: The `switch-mcp-profile.sh` script automatically reads keys from `.env` and injects them into the generated tool configurations (e.g., `.gemini/settings.json`).
3.  **Security**: Generated configuration files containing keys are automatically git-ignored to prevent accidental exposure.

**What gets installed automatically:**
- **AI CLIs:** Claude Code, Gemini CLI, Codex CLI (if selected)
- **Core MCPs:** Sequential Thinking, PAL, Context7
- **Research MCPs:** ToolUniverse (600+ tools), Serena (Code Search)
- **Smart Config:** Unified profile switcher for all CLIs

For detailed installation instructions, see [Installation Guide](docs/INSTALLATION.md).

## Context & Profile Management

SciAgent-toolkit includes a powerful **Profile Switcher** to manage context window usage across Claude, Gemini, and Codex.

Instead of loading all tools at once (which would consume 500k+ tokens), you select a profile optimized for your current task:

- `coding`: Lightweight, fast. (PAL + Context7)
- `research-lite`: Targeted scientific tools. (ToolUniverse core)
- `hybrid-research`: **Claude** for coding, **Gemini** for research.

👉 **Read the full guide:** [MCP Context Management](docs/MCP-CONTEXT-MANAGEMENT.md)

## Docker/Container Deployment

**For local development:** Use the installation scripts above.

**For production containerized deployments:** Use [scbio-docker](https://github.com/tony-zhelonkin/scbio-docker), which:
- Includes SciAgent-toolkit as a Git submodule
- Handles build-time MCP server installation
- Generates runtime MCP configurations
- Provides full scientific computing environment (R 4.5, Python 3.10, bioinformatics packages)

The Docker test images in `docker/test/` are for CI/CD validation only.

## What Gets Installed

- **Claude Code** or **Codex CLI** (AI interfaces)
- **PAL MCP Server** - Collaboration, planning, and code analysis
- **Serena MCP Server** - Code intelligence and semantic search
- **Sequential Thinking MCP Server** - Structured problem solving
- **ToolUniverse MCP Server** - 600+ scientific research tools
- **PubMed Plugin** - Biomedical literature access (manual setup)

## Use Cases

### Drug Discovery
Search molecular databases, check FDA approvals, find clinical trials, and review literature.

### Genomics Research
Query protein databases, analyze interactions, identify pathways, and search recent studies.

### Literature Review
Search PubMed/PMC, access full-text articles, extract findings, and map research gaps.

## Documentation

- [Installation Guide](docs/INSTALLATION.md) - Complete setup instructions
- [Quick Start](docs/QUICKSTART.md) - Get running in 5 minutes
- [Architecture](docs/ARCHITECTURE.md) - System design and data flow
- [Configuration](docs/CONFIGURATION.md) - Advanced customization
- [Troubleshooting](docs/TROUBLESHOOTING.md) - Common issues and solutions
- [FAQ](docs/FAQ.md) - Frequently asked questions

## Custom Agents

This toolkit includes pre-configured Claude agents for bioinformatics workflows:
- **Bioinformatics Research Librarian** - Find tools, docs, and resources
- **RNA-seq Methods Writer** - Auto-generate publication methods sections

See [agents/README.md](agents/README.md) for details.

## Requirements

**Minimum:**
- Python 3.10+
- Node.js 18+
- 2GB free disk space
- Internet connection
- macOS or Linux (WSL supported on Windows)

**Recommended:**
- Python 3.11+
- Node.js 20+
- 5GB free disk space (for ToolUniverse with all dependencies)
- Stable internet connection

## Known Limitations

- **Serena**: Optional component, requires 5-15 min build from source (skipped by default to save time)
- **PubMed**: Requires manual plugin installation via Claude Code marketplace
- **ToolUniverse**: Large installation (~300MB with dependencies)
- **SSH Keys**: Not required - all installations use HTTPS

### Installing Optional Serena

If you want code intelligence features, install Serena manually:

```bash
# This will take 5-15 minutes on first run (one-time build)
uvx --from git+https://github.com/oraios/serena serena start-mcp-server --help

# Then add to Claude Code
claude mcp add serena --scope local -- \
  uvx --from git+https://github.com/oraios/serena serena start-mcp-server
```

## License

This project is licensed under the MIT License - see [LICENSE](LICENSE) file for details.

**Important:** This toolkit is a configuration and installation framework. It does not own or claim any rights to the AI models, APIs, or external services it integrates with:
- Claude Code and related services are provided by Anthropic
- OpenAI Codex and related services are provided by OpenAI
- All MCP servers are separate open-source projects
- External APIs (PubMed, ChEMBL, UniProt, etc.) are provided by their respective organizations
- All external services require their own subscriptions, API keys, or authentication

## Contributing

Contributions are welcome! Please see [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines.

## Acknowledgments

This project integrates and builds upon:
- [Serena MCP Server](https://github.com/oraios/serena)
- [ToolUniverse](https://github.com/mims-harvard/ToolUniverse)
- [Sequential Thinking MCP Server](https://github.com/modelcontextprotocol/servers/tree/main/src/sequentialthinking)
- [Model Context Protocol](https://modelcontextprotocol.io)
