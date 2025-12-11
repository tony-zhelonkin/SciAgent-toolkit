# MCP Server Setup Guide

**Last Updated:** 2025-12-10
**Status:** ✅ ALL ISSUES RESOLVED

---

## Quick Start

### Inside a Dev Container

```bash
# 1. Run AI setup (installs Claude CLI, configures MCP servers)
./toolkits/SciAgent-toolkit/scripts/setup-ai.sh

# 2. Verify MCP servers
claude mcp list

# Expected output (without API keys):
# ✓ sequential-thinking: Connected
# ✓ context7: Connected
# ✓ tooluniverse: Connected
# ✓ serena: Connected
# ✗ pal: Failed (needs API keys)
```

### Configuring PAL at Runtime

PAL requires at least one AI provider API key. To enable it:

```bash
# 1. Edit your .env file
nano .devcontainer/.env

# 2. Add at least ONE of these keys:
GEMINI_API_KEY=your-key-here     # Get from: https://aistudio.google.com/apikey
OPENAI_API_KEY=your-key-here     # Get from: https://platform.openai.com/api-keys
# XAI_API_KEY=your-key-here      # Get from: https://console.x.ai/

# 3. Re-run configuration
./toolkits/SciAgent-toolkit/scripts/configure_mcp_servers.sh --force

# 4. Verify PAL is now connected
claude mcp list
```

---

## How MCP Setup Works

### Script Hierarchy

```
setup-ai.sh                          # Entry point (runs everything)
  └── setup_mcp_infrastructure.sh    # Installs Claude CLI, MCP servers
       ├── install_claude.sh         # Installs ~/.local/bin/claude
       ├── setup_tooluniverse.sh     # Creates tooluniverse-env/
       ├── setup_serena.sh           # Installs via uvx
       ├── setup_sequential_thinking.sh
       └── configure_mcp_servers.sh  # Generates .mcp.json (can run standalone)
```

### When to Use Each Script

| Script | Use Case |
|--------|----------|
| `setup-ai.sh` | First-time setup in new container |
| `configure_mcp_servers.sh --force` | After changing API keys or fixing issues |
| `debug_tooluniverse.sh` | Diagnosing ToolUniverse connection problems |

### Files Created

```
project/
├── .mcp.json                    # MCP server configuration (gitignored)
├── .claude/
│   ├── settings.json            # Enabled servers list
│   └── settings.local.json      # Local overrides
├── tooluniverse-env/            # ToolUniverse Python venv (gitignored)
└── .devcontainer/
    └── .env                     # API keys (gitignored)
```

---

## MCP Servers Reference

| Server | Purpose | Requirements |
|--------|---------|--------------|
| **sequential-thinking** | Structured reasoning for complex decisions | npx (Node.js) |
| **context7** | Up-to-date library documentation | npx; API key optional |
| **tooluniverse** | 600+ scientific tools (ChEMBL, UniProt, PubMed, etc.) | uv, Python 3.10+ |
| **serena** | Code intelligence and semantic search | uvx |
| **pal** | Multi-model AI collaboration (Gemini, GPT, Grok) | uvx + API key required |

### Context7 (No API Key Required)

Context7 works without an API key for basic usage. Add a key only for:
- Higher rate limits
- Private repository access

```bash
# Get optional key from: https://context7.com/dashboard
CONTEXT7_API_KEY=your-key-here
```

### PAL Multi-Model Collaboration

PAL enables Claude to consult other AI models for second opinions. Requires at least one key:

```bash
# Choose one or more:
GEMINI_API_KEY=xxx    # Google's Gemini models
OPENAI_API_KEY=xxx    # GPT-4, etc.
XAI_API_KEY=xxx       # Grok
OPENROUTER_API_KEY=xxx # Access to many models
```

---

## Troubleshooting

### ToolUniverse shows "connecting..." or "Failed"

```bash
# Run debug script
./toolkits/SciAgent-toolkit/scripts/debug_tooluniverse.sh --verbose

# Common fixes:
# 1. Reinstall
rm -rf tooluniverse-env
./toolkits/SciAgent-toolkit/scripts/mcp_servers/setup_tooluniverse.sh
./toolkits/SciAgent-toolkit/scripts/configure_mcp_servers.sh --force

# 2. Check the command being used (should be tooluniverse-smcp-stdio)
cat .mcp.json | grep tooluniverse
```

### PAL shows "Failed to connect"

This is **expected** if no API keys are set. See "Configuring PAL at Runtime" above.

### npm/nvm Prefix Warning

```
nvm is not compatible with the npm config "prefix" option
```

This is a **cosmetic warning** - MCP servers still work. To fix:
```bash
nvm use --delete-prefix $(node --version)
```

### Configuration Not Taking Effect

```bash
# 1. Force regenerate
./toolkits/SciAgent-toolkit/scripts/configure_mcp_servers.sh --force

# 2. Restart Claude Code
exit  # or Ctrl+D
claude

# 3. Verify
claude mcp list
```

---

## Technical Details

### ToolUniverse Transport Modes

ToolUniverse 1.0.14+ has two commands:
- `tooluniverse-mcp` → HTTP mode (default, NOT compatible with Claude Code)
- `tooluniverse-smcp-stdio` → stdio mode (required for Claude Code)

The configuration script automatically prefers `tooluniverse-smcp-stdio`.

### API Key Detection

The script sources `.env` files in this order:
1. `${PROJECT_DIR}/.env`
2. `${PROJECT_DIR}/.devcontainer/.env`

Keys must be exported or the script uses `set -a` to auto-export.

### Version Compatibility Notes

| Component | Version | Notes |
|-----------|---------|-------|
| ToolUniverse | 1.0.14+ | `--exclude-tool-types` removed; use `--compact-mode` |
| Claude CLI | 2.0.64+ | Required for `claude mcp add` |
| Node.js | 18+ | Required for npx-based servers |
| Python | 3.10+ | Required for ToolUniverse |

---

## Resolution History

**2025-12-10:** All issues resolved

| Issue | Root Cause | Fix |
|-------|------------|-----|
| PAL missing from `/mcp` | Script skipped if no API keys | Always configure, show helpful warnings |
| Context7 not configured | API key required but optional | Removed requirement |
| ToolUniverse "connecting..." | Wrong command (`tooluniverse-mcp` uses HTTP) | Prefer `tooluniverse-smcp-stdio` |
| ToolUniverse CLI error | `--exclude-tool-types` removed in 1.0.14+ | Removed flag from all scripts |
| Empty env vars in config | PAL env included empty strings | Only include non-empty values |
