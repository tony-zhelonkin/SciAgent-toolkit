# Advanced Configuration Guide

This guide covers advanced configuration options for the SciAgent Toolkit.

## Role System Configuration

The role system provides declarative configuration for agents and skills per project.

### Role Definition Format

```yaml
# roles/my-role.yaml
name: my-role
description: Description shown during activation
mcp_profile: coding  # Recommended MCP profile

agents:
  - bioinf-librarian
  - rnaseq-methods-writer

skills:
  - my-custom-skill
```

### Activating Roles

```bash
# Activate for current directory
./scripts/activate-role.sh base

# Activate for specific project
./scripts/activate-role.sh base --project-dir /path/to/project

# List available roles
ls roles/*.yaml
```

### What Role Activation Does

1. Creates `.claude/agents/` and `.claude/skills/` directories
2. Clears existing symlinks (enables role switching)
3. Symlinks agents from `agents/<name>.md` to `.claude/agents/`
4. Symlinks skills from `skills/<name>.md` to `.claude/skills/`
5. Displays recommended MCP profile

---

## Profile System Configuration

The profile system manages MCP server configuration across all AI CLIs.

### Available Profiles

| Profile | MCP Servers | Context | Use Case |
|---------|-------------|---------|----------|
| `minimal` | context7, sequential-thinking | ~3k | Fastest startup |
| `coding` | + pal | ~25k | General coding |
| `codebase` | + serena | ~75k | Code analysis |
| `research-lite` | + tooluniverse (core) | ~30k | Targeted research |
| `research-full` | + tooluniverse (14 tools) | ~50k | Scientific research |
| `full` | all servers | ~100k | Maximum capability |

### Switching Profiles

```bash
# Basic usage
./scripts/switch-mcp-profile.sh coding

# With explicit project directory
./scripts/switch-mcp-profile.sh research-lite --project-dir /path/to/project
```

### API Key Substitution

The profile switcher injects API keys from environment variables:

```bash
# Set in .env file or environment
export GEMINI_API_KEY="your-gemini-key"
export OPENAI_API_KEY="your-openai-key"
export CONTEXT7_API_KEY="your-context7-key"

# Keys are substituted in:
# - .gemini/settings.json
# - PAL MCP configuration
# - Context7 configuration
```

### Profile Validation

Before switching, the script validates:
- Required MCP servers are installed
- Required environment variables are set
- ToolUniverse environment exists (for research profiles)

---

## Tool Filtering for Better Performance

ToolUniverse provides 600+ tools, which can overwhelm the context window. You can filter tools to only include what you need.

### For Claude Code

Edit `.mcp.json` in your project directory.

#### Option 1: Direct Binary Execution (Recommended)

Point directly to the executable in the virtual environment. This method is simpler and avoids argument confusion.

```json
{
  "mcpServers": {
    "tooluniverse-research": {
      "type": "stdio",
      "command": "/path/to/SciAgent-toolkit/tooluniverse-env/bin/tooluniverse-smcp-stdio",
      "args": [
        "--include-tools",
        "EuropePMC_search_articles,ChEMBL_search_similar_molecules,search_clinical_trials"
      ]
    }
  }
}
```

#### Option 2: Using `uv`

Run the server via `uv`. Ensure you use the correct `uv` arguments (`--directory`, `run`).

```json
{
  "mcpServers": {
    "tooluniverse-research": {
      "type": "stdio",
      "command": "uv",
      "args": [
        "--directory",
        "/path/to/SciAgent-toolkit/tooluniverse-env",
        "run",
        "tooluniverse-smcp-stdio",
        "--include-tools",
        "EuropePMC_search_articles,ChEMBL_search_similar_molecules,search_clinical_trials"
      ]
    }
  }
}
```

### For Codex CLI

Edit `~/.codex/config.toml`:

```toml
[mcp_servers.tooluniverse-research]
command = "uv"
args = [
  "--directory",
  "/path/to/SciAgent-toolkit/tooluniverse-env",
  "run",
  "tooluniverse-smcp-stdio",
  "--include-tools",
  "EuropePMC_search_articles,ChEMBL_search_similar_molecules,search_clinical_trials"
]
```

### Filtering Options

#### Include Specific Tools

```bash
--include-tools "tool1,tool2,tool3"
```

#### Exclude Tool Types

```bash
--exclude-tool-types "type1,type2"
```

## Multiple ToolUniverse Instances

You can run multiple ToolUniverse instances with different tool sets for better organization and performance.

### Example: Specialized Instances

#### Instance 1: Literature Research

```json
{
  "mcpServers": {
    "tooluniverse-literature": {
      "type": "stdio",
      "command": "uv",
      "args": [
        "--directory",
        "/path/to/tooluniverse-env",
        "run",
        "tooluniverse-smcp-stdio",
        "--include-tools",
        "EuropePMC_search_articles,PubMed_search,Semantic_Scholar_search,OpenAlex_search"
      ]
    }
  }
}
```

#### Instance 2: Drug Discovery

```json
{
  "mcpServers": {
    "tooluniverse-drugs": {
      "type": "stdio",
      "command": "uv",
      "args": [
        "--directory",
        "/path/to/tooluniverse-env",
        "run",
        "tooluniverse-smcp-stdio",
        "--include-tools",
        "ChEMBL_search_similar_molecules,DrugBank_search,FDA_search_drugs,search_clinical_trials"
      ]
    }
  }
}
```

#### Instance 3: Genomics & Proteomics

```json
{
  "mcpServers": {
    "tooluniverse-genomics": {
      "type": "stdio",
      "command": "uv",
      "args": [
        "--directory",
        "/path/to/tooluniverse-env",
        "run",
        "tooluniverse-smcp-stdio",
        "--include-tools",
        "UniProt_search,protein_interaction_search,pathway_analysis"
      ]
    }
  }
}
```

## Azure OpenAI Configuration

Enable automatic summarization of long outputs from ToolUniverse.

### Setup

```bash
# Add to ~/.bashrc or ~/.zshrc
export AZURE_OPENAI_API_KEY="your-api-key"
export AZURE_OPENAI_ENDPOINT="https://your-resource.openai.azure.com"
export AZURE_OPENAI_DEPLOYMENT="your-deployment-name"  # Optional

# Reload configuration
source ~/.bashrc  # or source ~/.zshrc
```

### Re-run ToolUniverse Setup

```bash
./scripts/mcp_servers/setup_tooluniverse.sh
```

## Performance Tuning

### Startup Timeout

For slower systems, increase the MCP server startup timeout:

#### Codex CLI (`~/.codex/config.toml`)

```toml
[mcp_servers.tooluniverse]
command = "uv"
args = ["--directory", "/path/to/tooluniverse-env", "run", "tooluniverse-smcp-stdio"]
startup_timeout_sec = 120  # Increase from default 60
```

### Context Window Management

Tips for managing context window usage:

1. **Use Tool Filtering**: Reduce the number of available tools
2. **Multiple Specialized Instances**: Create focused tool sets
3. **Disable Unused Servers**: Comment out MCP servers you don't need

```json
{
  "mcpServers": {
    "serena": { ... },
    "sequential-thinking": { ... }
    // "tooluniverse": { ... }  // Commented out when not needed
  }
}
```

## Development vs Production Configs

### Development Configuration

Full tool access for exploration:

```json
{
  "mcpServers": {
    "tooluniverse-dev": {
      "type": "stdio",
      "command": "uv",
      "args": [
        "--directory",
        "/path/to/tooluniverse-env",
        "run",
        "tooluniverse-smcp-stdio"
      ]
    }
  }
}
```

### Production Configuration

Filtered, focused tool sets:

```json
{
  "mcpServers": {
    "tooluniverse-prod": {
      "type": "stdio",
      "command": "uv",
      "args": [
        "--directory",
        "/path/to/tooluniverse-env",
        "run",
        "tooluniverse-smcp-stdio",
        "--include-tools",
        "specific,required,tools,only"
      ]
    }
  }
}
```

## Custom MCP Server Integration

### Adding a New MCP Server

1. **Create a setup script**:

```bash
# scripts/mcp_servers/setup_newserver.sh
#!/usr/bin/env bash

echo "Installing NewServer MCP..."
# Installation logic here

echo "NewServer MCP installed successfully"
```

2. **Add to orchestrator**:

Edit `scripts/setup_mcp_infrastructure.sh` to include your new server.

3. **Update configuration**:

Add server to `.mcp.json` or `~/.codex/config.toml`.

4. **Document**:

Update [INSTALLATION.md](INSTALLATION.md) with server details.

## Environment-Specific Configurations

### Project-Specific (`.mcp.json`)

Best for:
- Project-specific tool requirements
- Collaborative projects with shared configs
- Different tool sets per project

Location: Project root directory

### User-Level (`~/.codex/config.toml`)

Best for:
- Personal preferences
- Consistent tools across all projects
- Codex CLI users

Location: `~/.codex/config.toml`

### Project-Level Codex Configuration

**Finding:** The Codex CLI currently *only* reads configuration from the global file `~/.codex/config.toml`. It does not natively support a project-local `.codex/config.toml` or `.codex.toml`.

**Recommended Strategy:**
To maintain project-specific settings (like specific tools or strict context limits):

1.  **Store Config in Project:** Create `.codex/config.toml` in your project root.
2.  **Sync on Context Switch:** Use a script (like `switch-mcp-profile.sh`) to copy/symlink this file to `~/.codex/config.toml` when you start working.

*Example Sync Script Snippet:*
```bash
# In your setup script
if [ -f "${PROJECT_ROOT}/.codex/config.toml" ]; then
    echo "Applying project-specific Codex config..."
    cp "${PROJECT_ROOT}/.codex/config.toml" "$HOME/.codex/config.toml"
fi
```

## Configuration Validation

### Validate Claude Code Configuration

```bash
# Check JSON syntax
python3 -m json.tool .mcp.json

# Run Claude diagnostics
claude doctor
```

### Validate Codex CLI Configuration

```bash
# Codex will report errors on startup
codex

# Check config file exists
cat ~/.codex/config.toml
```

## Troubleshooting Configuration Issues

### Issue: MCP Server Not Loading

**Check configuration syntax**:
```bash
python3 -m json.tool .mcp.json
```

**Verify paths are absolute**:
```json
{
  "args": [
    "--directory",
    "/absolute/path/to/tooluniverse-env",  // ✓ Absolute
    // NOT: "./tooluniverse-env"           // ✗ Relative
  ]
}
```

### Issue: Tool Filtering Not Working

**Verify tool names**:
```bash
# List all available tools
./test_tooluniverse.sh
```

**Check spelling and capitalization**:
Tool names are case-sensitive.

### Issue: Azure OpenAI Not Working

**Check environment variables**:
```bash
echo $AZURE_OPENAI_API_KEY
echo $AZURE_OPENAI_ENDPOINT
```

**Re-run setup**:
```bash
./scripts/mcp_servers/setup_tooluniverse.sh
```

## Best Practices

1. **Start Broad, Then Narrow**: Begin with all tools, then filter based on usage
2. **Use Multiple Instances**: Organize by research domain
3. **Document Custom Configs**: Keep notes on why you included/excluded tools
4. **Version Control Configs**: Commit `.mcp.json` for project sharing
5. **Test After Changes**: Always verify with `/mcp` command
6. **Keep Backups**: Save working configs before experimenting

## Example Complete Configuration

### Claude Code (`.mcp.json`)

```json
{
  "mcpServers": {
    "serena": {
      "type": "stdio",
      "command": "uvx",
      "args": [
        "--from",
        "git+https://github.com/oraios/serena",
        "serena",
        "start-mcp-server",
        "--base-dir",
        "."
      ]
    },
    "sequential-thinking": {
      "type": "stdio",
      "command": "npx",
      "args": [
        "-y",
        "@modelcontextprotocol/server-sequential-thinking"
      ]
    },
    "tooluniverse-literature": {
      "type": "stdio",
      "command": "uv",
      "args": [
        "--directory",
        "/Users/tony/Bioinf GDrive/Projects/SciAgent-toolkit/tooluniverse-env",
        "run",
        "tooluniverse-smcp-stdio",
        "--include-tools",
        "EuropePMC_search_articles,Semantic_Scholar_search"
      ]
    },
    "tooluniverse-drugs": {
      "type": "stdio",
      "command": "uv",
      "args": [
        "--directory",
        "/Users/tony/Bioinf GDrive/Projects/SciAgent-toolkit/tooluniverse-env",
        "run",
        "tooluniverse-smcp-stdio",
        "--include-tools",
        "ChEMBL_search_similar_molecules,FDA_search_drugs"
      ]
    }
  }
}
```

### Codex CLI (`~/.codex/config.toml`)

```toml
[mcp_servers.serena]
command = "uvx"
args = [
  "--from",
  "git+https://github.com/oraios/serena",
  "serena",
  "start-mcp-server",
  "--base-dir",
  "."
]

[mcp_servers.sequential-thinking]
command = "npx"
args = ["-y", "@modelcontextprotocol/server-sequential-thinking"]

[mcp_servers.tooluniverse]
command = "uv"
args = [
  "--directory",
  "/Users/tony/Bioinf GDrive/Projects/SciAgent-toolkit/tooluniverse-env",
  "run",
  "tooluniverse-smcp-stdio"
]
startup_timeout_sec = 60
```

---

For more information, see:
- [INSTALLATION.md](INSTALLATION.md) - Installation guide
- [TROUBLESHOOTING.md](TROUBLESHOOTING.md) - Common issues
- [ARCHITECTURE.md](ARCHITECTURE.md) - System architecture
