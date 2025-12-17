# SciAgent Toolkit

A modular infrastructure for building AI-powered scientific research agents using the Model Context Protocol (MCP). This toolkit enables integration of powerful AI assistants (Claude Code, Gemini CLI, Codex CLI) with specialized scientific tools for drug discovery, genomics, literature research, and more.

## Features

- **Multi-Agent Support** for [Claude Code](https://github.com/anthropics/claude-code), [Gemini CLI](https://github.com/google-gemini/gemini-cli), and [Codex CLI](https://github.com/openai/codex)
- **Structured Reasoning** via [Sequential Thinking MCP server](https://github.com/modelcontextprotocol/servers/tree/main/src/sequentialthinking) for Claude code
- **Multi-Agent Consensus and Brain-storming** via [PAL MCP](https://github.com/BeehiveInnovations/pal-mcp-server)
- **Role-Based Configuration** via declarative YAML role definitions
- **Smart Context Management** via profile switching (Coding vs Research modes)
- **Template System** for consistent AI context files across projects
- **Modular Installation** with selective component installation
- **Scientific Tools** via [ToolUniverse MCP](https://github.com/mims-harvard/ToolUniverse) server (ChEMBL, UniProt, DrugBank, FDA databases)
- **Biomedical Articles** via PubMed integration through Claude code [life-sciences marketplace](https://github.com/anthropics/life-sciences)
- **Code Intelligence** with [Serena MCP](https://github.com/oraios/serena) server (I\`m thinking of deprecating it)
- **Production-Ready** with comprehensive Docker testing (11 tests)





## Quick Start

```bash
# Clone the repository
git clone https://github.com/tony-zhelonkin/SciAgent-toolkit
cd SciAgent-toolkit

# 1. Setup Environment
cp templates/.env.template .env
# Edit .env and add your API keys (GEMINI_API_KEY, OPENAI_API_KEY, etc.)

# 2. Run the automated installer (creates templates, activates role, configures MCP)
./scripts/setup-ai.sh

# 3. Switch to a profile (e.g., for coding)
./scripts/switch-mcp-profile.sh coding

# Verify installation
claude
/mcp  # Should show: sequential-thinking, pal, context7
```

### Minimal Setup (Faster)

```bash
# Skip Serena, Codex, and Gemini for faster setup (2-3 minutes)
./scripts/setup-ai.sh --minimal
```

## Configuration & Secrets

SciAgent-toolkit uses a `.env` file for secure API key management. A `.env.template` is provided in the `templates/` directory.

1.  **Single Source of Truth**: The `.env` file (git-ignored) holds all your API keys.
2.  **Auto-Injection**: The `switch-mcp-profile.sh` script automatically reads keys from `.env` and injects them into the generated tool configurations (e.g., `.gemini/settings.json`).
3.  **Security**: Generated configuration files containing keys are automatically git-ignored to prevent accidental exposure.

**What gets installed automatically:**
- **AI CLIs:** Claude Code, Gemini CLI, Codex CLI (if selected)
- **Core MCPs:** Sequential Thinking, PAL, [Context7](https://context7.com/dashboard)
- **Research MCPs:** ToolUniverse (600+ tools), Serena (Code Search)
- **Smart Config:** Unified profile switcher for all CLIs

For detailed installation instructions, see [Installation Guide](docs/INSTALLATION.md).

## Context & Profile Management

SciAgent-toolkit includes a **Profile Switcher** to manage context window usage across Claude, Gemini, and Codex.

Instead of loading all tools at once (which would consume 500k+ tokens), you select a profile optimized for your current task:

- `coding`: Lightweight, fast. (PAL + Context7)
- `research-lite`: Targeted scientific tools. (ToolUniverse core)
- `hybrid-research`: **Claude** for coding, **Gemini** for research.

👉 **Read the full guide:** [MCP Context Management](docs/MCP-CONTEXT-MANAGEMENT.md)

## Docker/Container Deployment

**For local development:** Use the installation scripts above.

**For production containerized deployments:** Use [scbio-docker](https://github.com/tony-zhelonkin/scbio-docker), which:
- Set up to be the base container for AI tooling
- Handles build-time MCP server installation
- Generates runtime MCP configurations
- Provides full scientific computing environment (R 4.5, Python 3.10, bioinformatics packages)

The Docker test images in `docker/test/` are for CI/CD validation only.

## What Gets Installed

### AI CLIs
- **Claude Code** - Primary coding interface
- **Gemini CLI** - High-context coding and research research agent (optional)
- **Codex CLI** - Alternative terminal interface (optional)

### MCP Servers
- **PAL MCP Server** - Multi-model collaboration and planning
- **Context7** - Library documentation lookup
- **Sequential Thinking** - Structured problem solving
- **ToolUniverse** - 600+ scientific research tools
- **Serena** - Code intelligence and semantic search (optional, 5-15 min build)

### Project Files (via `setup-ai.sh`)
- **CLAUDE.md, GEMINI.md, AGENTS.md** - AI context files
- **context.md** - Scientific project context
- **.claude/agents/, .claude/skills/** - Role-activated agents/skills
- **.mcp.json** - MCP server configuration

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

## Role System & Custom Agents

I plan to implement a **role-based system** to configure agents and skills per project or switch roles at runtime during project analysis:

```bash
# Activate a role (symlinks agents/skills to .claude/)
./scripts/activate-role.sh base --project-dir /path/to/project
```

### Pre-configured Agents (Base Role)

**Research & Documentation**
- **Bioinformatics Research Librarian** - Find tools, docs, and resources via web research
- **Bio-Research Visualizer** - Deep biological mechanism research + visualization recommendations

**Data Exploration & Analysis**
- **RNA-seq Insight Explorer** - Explore RNAseq results with scientific skepticism

**Publication & Documentation**
- **RNA-seq Methods Writer** - Auto-generate publication methods sections from code
- **Figure Caption Generator** - Publication-quality captions for figures and tables (fire-and-forget)
- **Repo Doc Curator** - Audit and consolidate repository documentation

**Code Review & Quality**
- **Refactor Stage Reviewer** - Peer review of refactored analysis stages

**Session Management**
- **Handoff** - Timestamped session handoff documentation

### Creating Custom Roles
Create `roles/my-role.yaml`:
```yaml
name: my-role
description: Custom workflow role
mcp_profile: research-lite
agents:
  - bioinf-librarian
  - my-custom-agent
skills: []
```

See [agents/README.md](agents/README.md) for details.

## Requirements

**Minimum:**
- Python 3.10+
- Node.js 18+
- 2GB free disk space
- Internet connection
- macOS or Linux (WSL probably supported on Windows, but not tested)

**Recommended:**
- Python 3.11+
- Node.js 20+
- 5GB free disk space (for ToolUniverse with all dependencies)
- Stable internet connection

## Known Limitations

- **Serena**: Optional component, skipped by default
- **PubMed**: Requires manual plugin installation via Claude Code marketplace
- **ToolUniverse**: Large installation (~300MB with dependencies)

## License

This project is licensed under the MIT License - see [LICENSE](LICENSE) file for details.

**Important:** This toolkit is a configuration and installation framework. It does not own or claim any rights to the AI models, APIs, or external services it integrates with:
- Claude Code and related services are provided by Anthropic
- OpenAI Codex and related services are provided by OpenAI
- All MCP servers are separate open-source projects
- External APIs (PubMed, ChEMBL, UniProt, etc.) are provided by their respective organizations
- All external services require their own subscriptions, API keys, or authentication

## Contributing

Contributions are welcome! See [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines.

## Acknowledgments

This project integrates and builds upon:
- [Serena MCP Server](https://github.com/oraios/serena)
- [ToolUniverse](https://github.com/mims-harvard/ToolUniverse)
- [Sequential Thinking MCP Server](https://github.com/modelcontextprotocol/servers/tree/main/src/sequentialthinking)
- [Model Context Protocol](https://modelcontextprotocol.io)
