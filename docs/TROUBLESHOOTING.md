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

## PAL clink to Gemini CLI Issues

### Issue: PAL clink fails with "YOLO mode is disabled by settings"

**Symptoms:**
- `mcp__pal__clink(cli_name: "gemini", ...)` returns error
- Error message: `Cannot start in YOLO mode when it is disabled by settings`
- Return code 52

**Cause:** PAL's embedded Gemini CLI config (`conf/cli_clients/gemini.json`) adds `--yolo` flag, but user's Gemini settings have `disableYoloMode: true`.

**Fix:**
```bash
# Enable YOLO mode in Gemini settings
cat ~/.gemini/settings.json
# Change "disableYoloMode": true → "disableYoloMode": false
```

Or manually edit:
```json
{
  "security": {
    "disableYoloMode": false
  }
}
```

---

### Issue: PAL clink fails with "Could not find child token in parent raw content"

**Symptoms:**
- Gemini CLI crashes during startup
- Error: `[ERROR] [ImportProcessor] Could not find child token in parent raw content`
- Often references content from `GEMINI.md` or `AGENTS.md`

**Cause:** Gemini CLI's markdown parser has trouble with certain constructs, especially:
- Blockquotes (`>`) containing bold (`**`) text
- Complex nested markdown formatting

**Example problematic content:**
```markdown
> **Full methodology:** See `AGENTS.md` in this directory...
```

**Fix:** Simplify the markdown in your project's `GEMINI.md`:
```markdown
# Change FROM:
> **Full methodology:** See `AGENTS.md` in this directory for guidelines.

# Change TO:
Full methodology is documented in `AGENTS.md`. This file contains Gemini-specific context.
```

---

### Issue: PAL clink fails with "GEMINI_API_KEY environment variable" error

**Symptoms:**
- Gemini CLI starts but immediately exits
- Error: `When using Gemini API, you must specify the GEMINI_API_KEY environment variable`

**Cause:** The `.env` file containing API keys is not loaded into the shell environment. Claude Code runs in a separate process from your terminal, so manual `source .env` doesn't propagate.

**Root Cause:** Unlike PAL's internal API calls (which can read `.env` files), `clink` spawns an external CLI process that inherits environment variables from the shell. If `GEMINI_API_KEY` isn't exported, the Gemini CLI can't authenticate.

**Workaround Options:**

1. **Add to shell profile** (recommended):
   ```bash
   echo 'source /path/to/project/.devcontainer/.env' >> ~/.bashrc
   source ~/.bashrc
   # Restart Claude Code
   ```

2. **Export before starting Claude Code:**
   ```bash
   export GEMINI_API_KEY="your-key-here"
   claude
   ```

3. **Use PAL's direct model calls instead of clink:**
   ```
   mcp__pal__chat(prompt: "...", model: "gemini-2.5-pro")
   ```
   PAL's direct API calls can source `.env` files internally.

**Note:** This is an architectural limitation - see ISSUES.md for ISSUE-018.

---

## Project Analysis Troubleshooting

### Memory Issues

**Issue: R session crashes during large operations**

**Solution:**
```r
# Check memory usage
gc()
pryr::mem_used()

# Force garbage collection between operations
rm(large_object); gc()

# Use disk-backed operations for very large data
options(future.globals.maxSize = 8 * 1024^3)  # 8GB limit
```

### Peak Calling Issues

**Issue: MACS3 can't find fragment files**

**Diagnosis:**
```r
# Check fragments path
Fragments(seurat_obj)

# Verify file exists
file.exists(Fragments(seurat_obj)[[1]]@path)
```

**Solution:** Use original fragment paths, not symlinks or copied files.

### Checkpoint Loading Issues

**Issue: "cannot open connection" or "object not found"**

**Diagnosis:**
```bash
# Check file exists and size
ls -lh /path/to/checkpoint.rds

# Check file integrity
R -e 'x <- readRDS("/path/to/checkpoint.rds"); class(x)'
```

**Solution:** Re-run the stage that creates the checkpoint, or restore from backup.

### Quick Health Check Script

Run this to diagnose common issues:

```bash
#!/bin/bash
echo "=== DC_Dictionary Health Check ==="
echo ""
echo "1. Environment Variables"
for var in GEMINI_API_KEY OPENAI_API_KEY CONTEXT7_API_KEY; do
    [ -n "${!var}" ] && echo "   $var: SET" || echo "   $var: NOT SET"
done
echo ""
echo "2. Gemini Settings"
[ -f ~/.gemini/settings.json ] && cat ~/.gemini/settings.json | grep -E "(disableYoloMode|auth)" || echo "   No settings file"
echo ""
echo "3. R Environment"
which R && R --version | head -1 || echo "   R not found"
echo ""
echo "4. Working Directory"
pwd
ls -la *.md 2>/dev/null | head -5
```

---

## Getting More Help

- [INSTALLATION.md](INSTALLATION.md) - Full installation guide
- [CONFIGURATION.md](CONFIGURATION.md) - Advanced configuration
- [ISSUES.md](ISSUES.md) - Known architectural issues
- [FAQ.md](FAQ.md) - Frequently asked questions