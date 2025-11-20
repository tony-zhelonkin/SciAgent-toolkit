# SciAgent Toolkit

A modular infrastructure for building AI-powered scientific research agents using the Model Context Protocol (MCP). This toolkit enables integration of powerful AI assistants (Claude Code, Codex CLI) with specialized scientific tools for drug discovery, genomics, literature research, and more.

## Features

- **600+ Scientific Tools** via ToolUniverse MCP server (ChEMBL, UniProt, DrugBank, FDA databases)
- **36M+ Biomedical Articles** via PubMed integration
- **Code Intelligence** with Serena MCP server
- **Structured Reasoning** via Sequential Thinking MCP server
- **Dual Interface Support** for Claude Code and Codex CLI
- **Modular Installation** with selective component installation
- **Production-Ready** with comprehensive error handling and validation

## Quick Start

```bash
# Clone the repository
git clone https://github.com/[your-username]/SciAgent-toolkit
cd SciAgent-toolkit

# Run the installer
./scripts/setup_mcp_infrastructure.sh

# Verify installation
claude
/mcp  # Should show all installed MCP servers
```

For detailed installation instructions, see [Installation Guide](docs/INSTALLATION.md).

## What Gets Installed

- **Claude Code** or **Codex CLI** (AI interfaces)
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

- Python 3.10+
- Node.js 18+ (for Sequential Thinking server)
- macOS or Linux (WSL supported on Windows)
- Internet connection

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
- [ToolUniverse](https://github.com/greedyai/tool-universe)
- [Sequential Thinking MCP Server](https://github.com/modelcontextprotocol/server-sequential-thinking)
- [Model Context Protocol](https://modelcontextprotocol.io/)
