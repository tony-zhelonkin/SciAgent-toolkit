# Troubleshooting

Common issues and solutions for SciAgent-toolkit MCP server configuration.

---

## Quick Diagnostic Commands

```bash
# Check MCP server status
claude mcp list

# Validate configuration file
python3 -m json.tool .mcp.json

# Run diagnostics
claude doctor

# Debug ToolUniverse specifically
./scripts/debug_tooluniverse.sh --verbose
```

---

## MCP Server Connection Issues

### PAL shows "Failed to connect"

**Symptoms:**
- `/mcp` shows PAL with `✕ failed` status
- Error about missing environment variables

**Cause:** PAL requires at least one AI provider API key.

**Fix:**
```bash
# 1. Edit your .env file
nano .devcontainer/.env

# 2. Add at least ONE of these keys:
GEMINI_API_KEY=your-key-here     # Get from: https://aistudio.google.com/apikey
OPENAI_API_KEY=your-key-here     # Get from: https://platform.openai.com/api-keys

# 3. Re-run configuration
./toolkits/SciAgent-toolkit/scripts/configure_mcp_servers.sh --force

# 4. Restart Claude Code
```

---

### ToolUniverse shows "connecting..." or "Failed"

**Symptoms:**
- `/mcp` shows ToolUniverse stuck on "connecting..."
- Logs show `ModuleNotFoundError: No module named 'tooluniverse'`

**Cause:** Virtual environment is corrupted or `uv` is detecting the wrong Python.

**Fix:**

1. **Reinstall the environment:**
   ```bash
   rm -rf tooluniverse-env
   uv venv tooluniverse-env
   uv pip install --python tooluniverse-env/bin/python tooluniverse
   ```

2. **Verify installation:**
   ```bash
   ./tooluniverse-env/bin/tooluniverse-smcp-stdio --help
   ```

3. **Regenerate configuration:**
   ```bash
   ./toolkits/SciAgent-toolkit/scripts/configure_mcp_servers.sh --force
   ```

**Alternative Fix (Direct Binary):**

If `uv run` keeps failing, edit `.mcp.json` to use the direct executable:

```json
"tooluniverse": {
  "command": "/absolute/path/to/tooluniverse-env/bin/tooluniverse-smcp-stdio",
  "args": ["--include-tools", "..."]
}
```

---

### Context7 not connecting

**Symptoms:** Context7 shows connection errors.

**Note:** Context7 works without an API key for basic usage. Only add a key for higher rate limits or private repository access.

```bash
# Optional key from: https://context7.com/dashboard
CONTEXT7_API_KEY=your-key-here
```

---

### Configuration not taking effect

**Symptoms:** Changes to `.mcp.json` aren't reflected in Claude Code.

**Fix:**
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

## ToolUniverse Specific Issues

### Issue: ToolUniverse installation fails with "Permission denied"

**Symptoms:**
- `setup_tooluniverse.sh` fails
- Error: `error: failed to remove file ... Permission denied`

**Fix:** Use explicit Python targeting:
```bash
uv pip install --python tooluniverse-env/bin/python tooluniverse
```

---

### Issue: ToolUniverse fails with invalid arguments

**Symptoms:**
- Errors about `--directory` or `--exclude-tool-types` being unrecognized

**Cause:** Mixing direct binary execution with `uv` arguments, or using deprecated CLI flags.

**Fix:** Use the `research-lite` profile or edit `.mcp.json`:
```bash
./scripts/switch-mcp-profile.sh research-lite
```

Or manually edit `.mcp.json` to remove conflicting arguments.

---

### Issue: Gemini - ToolUniverse tools not appearing

**Symptoms:**
- `find_tools` returns empty list
- Tools listed in `.gemini/settings.json` but not available

**Cause:** Tool names passed as single comma-separated string instead of separate arguments.

**Incorrect:**
```json
"args": ["--include-tools", "Tool1,Tool2"]
```

**Correct:**
```json
"args": ["--include-tools", "Tool1", "Tool2"]
```

---

### Issue: Context window overflow with ToolUniverse

**Symptoms:**
- Claude becomes slow or unresponsive
- Context usage shows 100%+

**Cause:** Loading all 600+ ToolUniverse tools exceeds context window.

**Fix:** Use a filtered profile:
```bash
./scripts/switch-mcp-profile.sh research-lite  # 6 essential tools
```

Or use `hybrid` profile (Claude=coding, Gemini=research).

---

## Codex CLI Issues

### Issue: Installation fails with `npm error code EACCES`

**Symptoms:**
- `npm install -g @openai/codex` fails
- Error: `Error: EACCES: permission denied, mkdir '/opt/nvm/versions/node/...'`

**Cause:** npm trying to write to system directory.

**Fix:**
```bash
# 1. Create user-writable directory
mkdir -p "$HOME/.npm-global"

# 2. Configure npm
npm config set prefix "$HOME/.npm-global"

# 3. Add to PATH
export PATH="$HOME/.npm-global/bin:$PATH"
echo 'export PATH="$HOME/.npm-global/bin:$PATH"' >> "$HOME/.bashrc"

# 4. Retry installation
npm install -g @openai/codex
```

---

### Issue: MCP Handshake Failure (Stdio Noise)

**Symptoms:**
- Codex shows `⚠ MCP client for 'tooluniverse' failed to start`
- Error: `handshaking with MCP server failed` or `JSON parse error`
- Claude/Gemini connect fine to the same server

**Cause:** Codex CLI enforces strict MCP over Stdio compliance. Some servers print logs to stdout, corrupting JSON-RPC messages.

**Workarounds:**
- Use the server's quiet mode if available
- Use direct binary instead of `uv run`
- Exclude the server from Codex until patched

---

## Installation Issues

### Issue: Command not found after install

```bash
source ~/.bashrc
# OR restart terminal
```

---

### Issue: No .mcp.json file created

```bash
# Configuration script wasn't run or failed
./scripts/configure_mcp_servers.sh

# Verify
cat .mcp.json
```

---

### Issue: Serena installation failed (SSH error)

**Cause:** Git SSH configuration issue.

**Fix:**
```bash
git pull origin main
./scripts/mcp_servers/setup_serena.sh

# Test
uvx --from git+https://github.com/oraios/serena serena --help
```

---

### Issue: npm/nvm Prefix Warning

```
nvm is not compatible with the npm config "prefix" option
```

This is a **cosmetic warning** - MCP servers still work.

**Optional fix:**
```bash
nvm use --delete-prefix $(node --version)
```

---

## Profile Switching Issues

### Issue: PAL fails after profile switch

**Symptoms:** PAL worked, then fails after `switch-mcp-profile.sh`.

**Cause:** Some profiles use a wrapper script, others use direct `uvx`. The wrapper may not exist.

**Fix:**
```bash
# Use a profile with direct uvx invocation
./scripts/switch-mcp-profile.sh hybrid

# Or install the wrapper
./scripts/mcp_servers/setup_pal.sh
```

---

### Issue: API keys not applied after profile switch

**Cause:** Profile templates use placeholders that require shell environment variables to be exported.

**Fix:**
```bash
# Export keys in your shell
export GEMINI_API_KEY="your-key"
export OPENAI_API_KEY="your-key"

# Re-run profile switch
./scripts/switch-mcp-profile.sh coding
```

---

### Issue: Python IndentationError during profile switch

**Status:** Fixed in versions after 2025-12-16.

**Fix:**
```bash
cd toolkits/SciAgent-toolkit
git pull origin main
```

---

## Technical Reference

### ToolUniverse Transport Modes

ToolUniverse 1.0.14+ has two commands:
- `tooluniverse-mcp` → HTTP mode (NOT compatible with Claude Code)
- `tooluniverse-smcp-stdio` → stdio mode (required for Claude Code)

The configuration scripts automatically prefer `tooluniverse-smcp-stdio`.

### API Key Detection Order

Scripts source `.env` files in this order:
1. `${PROJECT_DIR}/.env`
2. `${PROJECT_DIR}/.devcontainer/.env`

### Version Compatibility

| Component | Version | Notes |
|-----------|---------|-------|
| ToolUniverse | 1.0.14+ | Use `tooluniverse-smcp-stdio` |
| Claude CLI | 2.0.64+ | Required for `claude mcp add` |
| Node.js | 18+ | Required for npx-based servers |
| Python | 3.10+ | Required for ToolUniverse |

---

## Getting More Help

- [INSTALLATION.md](INSTALLATION.md) - Full installation guide
- [CONFIGURATION.md](CONFIGURATION.md) - Advanced configuration
- [ISSUES.md](ISSUES.md) - Known architectural issues
- [FAQ.md](FAQ.md) - Frequently asked questions
