# Frequently Asked Questions

## General Questions

### What is the SciAgent Toolkit?

The SciAgent Toolkit is a modular infrastructure that connects AI assistants (Claude Code, Codex CLI) with scientific research tools through the Model Context Protocol (MCP). It provides access to 600+ scientific tools, 36M+ biomedical articles, and specialized reasoning capabilities.

### Do I need both Claude Code and Codex CLI?

No. You can install either one or both:
- **Claude Code**: Best for IDE-integrated workflows and PubMed access
- **Codex CLI**: Best for terminal-based workflows
- **Both**: Maximum flexibility

Use `--skip-claude` or `--skip-codex` flags to install selectively.

### Is this free to use?

The toolkit itself is free and open source (MIT License). However:
- **Claude Code** requires an Anthropic account
- **Codex CLI** requires ChatGPT or OpenAI API key
- **External APIs** may have usage limits or require subscriptions
- **Azure OpenAI** (optional summarization) requires Azure account

### What operating systems are supported?

- macOS (Intel and Apple Silicon)
- Linux (Ubuntu, Debian, Fedora, etc.)
- Windows via WSL (Windows Subsystem for Linux)

## Installation Questions

### How long does installation take?

- **Quick setup**: 5-10 minutes
- **Full setup**: 10-20 minutes (depending on internet speed)

Most time is spent downloading dependencies.

### Can I install components separately?

Yes! Each component can be installed independently:

```bash
./scripts/install_claude.sh                      # Claude Code only
./scripts/install_codex.sh                       # Codex CLI only
./scripts/mcp_servers/setup_serena.sh            # Serena only
./scripts/mcp_servers/setup_tooluniverse.sh      # ToolUniverse only
```

### What if I already have Claude Code installed?

The installer detects existing installations and skips them. It's safe to run multiple times.

### Can I uninstall everything?

Yes. See the [Uninstallation section in INSTALLATION.md](INSTALLATION.md#uninstallation).

## MCP Server Questions

### What is an MCP server?

An MCP (Model Context Protocol) server is a standardized way to extend AI assistants with custom tools and capabilities. Think of them as plugins or extensions.

### Which MCP servers are included?

1. **Serena**: Code intelligence and semantic search
2. **Sequential Thinking**: Structured reasoning
3. **ToolUniverse**: 600+ scientific research tools
4. **PubMed** (manual): Biomedical literature access

### Can I add my own MCP servers?

Yes! See [CONFIGURATION.md - Custom MCP Server Integration](CONFIGURATION.md#custom-mcp-server-integration).

### Why doesn't PubMed install automatically?

PubMed is distributed through Claude Code's plugin marketplace, which requires manual installation via the `/plugin` commands. This is a Claude Code-specific feature.

## ToolUniverse Questions

### What tools are included in ToolUniverse?

600+ scientific tools across:
- Drug discovery (ChEMBL, DrugBank, FDA)
- Genomics (UniProt, protein interactions)
- Literature (Europe PMC, Semantic Scholar, OpenAlex)
- Clinical trials (ClinicalTrials.gov)
- Molecular biology (protein structures, pathways)

### Why so many tools?

ToolUniverse aims to be comprehensive. However, you can filter to only the tools you need using `--include-tools` or `--exclude-tool-types`. See [CONFIGURATION.md](CONFIGURATION.md).

### Do I need Azure OpenAI?

No. Azure OpenAI is optional and only used for summarizing long API responses. ToolUniverse works fine without it.

### Can I use multiple ToolUniverse instances?

Yes! Create specialized instances for different research domains. See [CONFIGURATION.md - Multiple ToolUniverse Instances](CONFIGURATION.md#multiple-tooluniverse-instances).

## Usage Questions

### How do I know which tools are available?

```bash
# In Claude Code or Codex CLI
/mcp

# Or ask
"What ToolUniverse tools are available?"
"List all scientific tools"
```

### Can I use this for commercial research?

Yes, but:
- Check the terms of service for each external API (PubMed, ChEMBL, etc.)
- Comply with Claude Code or Codex CLI license agreements
- Some databases have restrictions on commercial use

### How do I cite this in publications?

If you use this toolkit in research:
1. Cite the specific tools you used (ChEMBL, PubMed, UniProt, etc.)
2. Mention the AI assistant used (Claude Code or Codex CLI)
3. Optionally acknowledge this toolkit in your methods section

Example:
```
Literature searches were performed using PubMed (NCBI) accessed through
the SciAgent Toolkit (https://github.com/[your-repo]) via Claude Code
(Anthropic). Molecular similarity searches were conducted using ChEMBL...
```

### Can I share my configurations with collaborators?

Yes! The `.mcp.json` file can be committed to version control and shared. Just ensure paths are adjusted for each user's system or use relative paths where possible.

## Performance Questions

### Why is startup slow?

Initial startup loads all MCP servers. To speed up:
1. Use tool filtering to reduce the number of tools
2. Comment out unused MCP servers
3. Increase timeout settings (see [CONFIGURATION.md](CONFIGURATION.md))

### Why is memory usage high?

600+ tools create a large context. Solutions:
1. Filter tools with `--include-tools`
2. Use multiple specialized instances
3. Disable unused MCP servers

### Can I run this on a low-spec machine?

Minimum requirements:
- 4GB RAM (8GB recommended)
- 2GB free disk space
- Python 3.10+
- Node.js 18+

For low-spec systems, use aggressive tool filtering.

## Troubleshooting Questions

For detailed troubleshooting, see **[TROUBLESHOOTING.md](TROUBLESHOOTING.md)**.

Quick diagnostic:
```bash
claude mcp list              # Check server status
claude doctor                # Run diagnostics
python3 -m json.tool .mcp.json  # Validate config
```

## Configuration Questions

### What's the difference between .mcp.json and config.toml?

- **`.mcp.json`**: Claude Code configuration (project-specific)
- **`~/.codex/config.toml`**: Codex CLI configuration (user-level)

Choose based on which tool(s) you're using.

### Can I use environment-specific configs?

Yes! Maintain separate configurations:
- Development: All tools, verbose logging
- Production: Filtered tools, essential only

See [CONFIGURATION.md - Development vs Production](CONFIGURATION.md#development-vs-production-configs).

### How do I filter tools?

```json
{
  "args": [
    "--include-tools",
    "tool1,tool2,tool3"
  ]
}
```

Or exclude by type:
```json
{
  "args": [
    "--exclude-tool-types",
    "type1,type2"
  ]
}
```

See [CONFIGURATION.md](CONFIGURATION.md) for details.

## Research Workflow Questions

### What's the best workflow for drug discovery?

See [examples/drug-discovery-workflow.md](../examples/drug-discovery-workflow.md) for a detailed example.

General pattern:
1. ChEMBL: Find similar molecules
2. FDA: Check safety profiles
3. ClinicalTrials.gov: Find ongoing trials
4. PubMed: Review literature
5. Sequential Thinking: Synthesize findings

### How do I perform literature reviews efficiently?

1. Start broad: PubMed search
2. Get full text: PMC access
3. Find related: Semantic Scholar
4. Cross-reference: Europe PMC
5. Synthesize: Sequential Thinking

See [examples/literature-review-workflow.md](../examples/literature-review-workflow.md).

### Can I combine multiple data sources?

Absolutely! That's the power of this toolkit. Example:

```
1. "Find protein information for TP53 in UniProt"
2. "Search PubMed for TP53 mutations in cancer"
3. "Find clinical trials for TP53-targeted therapies"
4. "Analyze and synthesize all findings"
```

## Contributing Questions

### Can I contribute new features?

Yes! See [CONTRIBUTING.md](../CONTRIBUTING.md) for guidelines.

Contributions welcome:
- New MCP servers
- Documentation improvements
- Bug fixes
- Example workflows
- Configuration templates

### How do I report bugs?

1. Check [TROUBLESHOOTING.md](TROUBLESHOOTING.md)
2. Search existing GitHub issues
3. Create a new issue with:
   - Operating system
   - Python/Node.js versions
   - Full error messages
   - Steps to reproduce

### Can I add new MCP servers to the toolkit?

Yes! Follow these steps:
1. Create setup script in `scripts/mcp_servers/`
2. Test thoroughly
3. Document in [INSTALLATION.md](INSTALLATION.md)
4. Submit pull request

See [CONTRIBUTING.md](../CONTRIBUTING.md).

## Advanced Questions

### Can I run MCP servers remotely?

MCP servers run locally by default. Remote execution requires custom MCP server implementation. See [MCP Protocol documentation](https://modelcontextprotocol.io/).

### Can I create custom agents?

Yes! See `agents/` directory for examples:
- Bioinformatics Research Librarian
- RNA-seq Methods Writer

Create new agents by following the existing format.

### How do I optimize for specific research domains?

1. Create domain-specific ToolUniverse instances
2. Filter to relevant tools only
3. Configure custom agents for your domain
4. Document workflows in `examples/`

### Can I integrate with Jupyter notebooks?

While the toolkit is CLI-focused, you can:
1. Call CLI tools from notebooks via subprocess
2. Use APIs directly in notebooks
3. Use MCP servers programmatically (advanced)

## Still Have Questions?

- Check [INSTALLATION.md](INSTALLATION.md) for setup details
- See [TROUBLESHOOTING.md](TROUBLESHOOTING.md) for common issues
- Review [CONFIGURATION.md](CONFIGURATION.md) for advanced options
- Read [ARCHITECTURE.md](ARCHITECTURE.md) to understand how it works
- Open a GitHub issue for specific problems
