# Troubleshooting

Common issues and solutions for SciAgent-toolkit MCP server configuration.

**Related docs**: [ISSUES.md](ISSUES.md) (detailed issue tracking) | [ARCHITECTURE_REVIEW.md](ARCHITECTURE_REVIEW.md) (strategic refactoring plan)

---

## Environment Variable Loading (Common Root Cause)

Many issues below stem from API keys not being available in the shell. The `.devcontainer/.env` file is NOT auto-sourced into interactive shells, and `~/.bashrc` resets on every container rebuild.

### Permanent Fix: `source_env.sh` (Recommended)

Place an idempotent injection script at `.devcontainer/scripts/source_env.sh` and call it from `postCreateCommand` in `devcontainer.json`. This survives container rebuilds.

**`.devcontainer/scripts/source_env.sh`:**
```bash
#!/usr/bin/env bash
# source_env.sh - Inject .env auto-sourcing into shell profile
# Called by postCreateCommand to survive container rebuilds.
set -euo pipefail

MARKER="_source_project_env"
TARGET="${HOME}/.bashrc"

# Idempotent: skip if already present
if grep -q "$MARKER" "$TARGET" 2>/dev/null; then
    exit 0
fi

cat >> "$TARGET" << 'EOF'

# Auto-source project .env files (API keys for Gemini, PAL, etc.)
_source_project_env() {
    local project_dir=""
    if [ -d "/workspaces" ]; then
        for d in /workspaces/*/; do
            [ -f "${d}.devcontainer/.env" ] && project_dir="${d%/}" && break
        done
    fi
    [ -z "$project_dir" ] && project_dir="${PWD}"
    [ -f "${project_dir}/.devcontainer/.env" ] && { set -a; source "${project_dir}/.devcontainer/.env"; set +a; }
    [ -f "${project_dir}/.env" ] && { set -a; source "${project_dir}/.env"; set +a; }
}
_source_project_env
EOF
```

**`devcontainer.json` hook:**
```json
"postCreateCommand": "bash .devcontainer/scripts/source_env.sh && echo 'Container up and running'"
```

**Why this works:** `postCreateCommand` runs once per container creation, injecting the function into the fresh `~/.bashrc`. Every subsequent shell (new terminal, tmux pane, `bash -l`) inherits the exported keys.

**First-time setup for an existing container** (if you haven't rebuilt yet):
```bash
bash .devcontainer/scripts/source_env.sh
```

### Quick one-liner (current session only):
```bash
set -a && source .devcontainer/.env && set +a
```

### Why manual `.bashrc` edits alone don't persist

Dev Containers recreate `~/.bashrc` from the base image on every rebuild. Adding the function manually works until the next rebuild. The `source_env.sh` + `postCreateCommand` approach re-injects it automatically.

**Verify:** `echo "Key set: ${GEMINI_API_KEY:+yes}"`

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

### PAL fails with "ModuleNotFoundError: No module named 'pandas'" (ISSUE-019)

**Symptoms:**
- `/mcp` shows PAL with `✕ failed` status
- Error traceback mentions your project's `config.py` file
- Error: `from config import ...` followed by `ModuleNotFoundError`

**Cause:** Your project has a `config.py` file in `PYTHONPATH`, which shadows PAL's internal `config.py` module. This is a Python namespace collision.

**Fix:**
```bash
# 1. Add PYTHONPATH clearing to PAL's env block in .mcp.json
# Edit .mcp.json and ensure the PAL section looks like:
{
  "pal": {
    "env": {
      "PYTHONPATH": "",
      "GEMINI_API_KEY": "${GEMINI_API_KEY}",
      "OPENAI_API_KEY": "${OPENAI_API_KEY}"
    }
  }
}

# 2. Or re-run profile switcher (templates now include this fix)
source .devcontainer/.env
./01_Scripts/SciAgent-toolkit/scripts/switch-mcp-profile.sh coding

# 3. Restart Claude Code
```

**Note:** This fix was added to all profile templates in 2026-01-18. If using older templates, manually add `"PYTHONPATH": ""` to PAL's env block.

---

### API keys show as "${GEMINI_API_KEY}" literal in .mcp.json (ISSUE-020)

**Symptoms:** `.mcp.json` contains literal `${GEMINI_API_KEY}` instead of actual key value.

**Fix:** See [Environment Variable Loading](#environment-variable-loading-common-root-cause) above, then re-run profile switcher:
```bash
./01_Scripts/SciAgent-toolkit/scripts/switch-mcp-profile.sh coding
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

## Gemini CLI Issues

### Issue: Gemini CLI says "must specify GEMINI_API_KEY" despite key in .env

**Symptoms:** `gemini` shows "must specify GEMINI_API_KEY" even though key is in `.devcontainer/.env`.

**Fix:** See [Environment Variable Loading](#environment-variable-loading-common-root-cause) above. The `.env` file is NOT auto-exported to shell sessions.

---

### Issue: PAL shows "Disconnected" in Gemini CLI (ISSUE-019 variant)

**Symptoms:**
- Running `gemini` and typing `/mcp` shows:
  ```
  🟢 sequential-thinking - Ready
  🟢 context7 - Ready
  🔴 pal - Disconnected
  ```
- PAL works fine in Claude Code but fails in Gemini

**Cause:** Same PYTHONPATH collision issue as Claude Code. Your project has a `config.py` file that shadows PAL's internal module. The `.gemini/settings.json` file needs `PYTHONPATH: ""` in the PAL env block.

**Fix:**
```bash
# Edit .gemini/settings.json and add PYTHONPATH to PAL env block:
{
  "pal": {
    "env": {
      "PYTHONPATH": "",        # <-- Add this line
      "GEMINI_API_KEY": "...",
      "OPENAI_API_KEY": "..."
    }
  }
}
```

Or re-run the profile switcher (templates now include this fix):
```bash
set -a && source .devcontainer/.env && set +a
./01_Scripts/SciAgent-toolkit/scripts/switch-mcp-profile.sh coding
```

Then restart Gemini CLI.

---

### Issue: API keys not loaded in new terminal sessions (Devcontainer)

**Symptoms:** New terminal tabs lose API keys even though first terminal worked.

**Fix:** Use the `source_env.sh` script approach described in [Environment Variable Loading](#environment-variable-loading-common-root-cause) above. This injects the auto-sourcing function into `~/.bashrc` via `postCreateCommand`, so it survives container rebuilds and applies to every new shell.

**Alternative (requires container rebuild):** Add `env_file: [.env]` to your `docker-compose.yml` service definition, then rebuild via **Ctrl+Shift+P** → "Dev Containers: Rebuild Container".

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

### Issue: PAL setup fails with "Broken symlink" error

**Symptoms:**
- `setup-ai.sh` fails during PAL MCP installation
- Error message:
  ```
  error: Failed to inspect Python interpreter from active virtual environment at `venv/bin/python3`
    Caused by: Broken symlink at `venv/bin/python3`, was the underlying Python interpreter removed?
  hint: Consider recreating the environment (e.g., with `uv venv`)
  [ERROR] Failed to install pal-mcp-server
  ```

**Cause:** A stale `venv` directory exists from a previous setup on a different machine or container. The Python symlinks inside point to a non-existent path (e.g., another user's home directory or a removed Python installation).

Example broken symlink:
```
venv/bin/python -> /home/other-user/.local/share/uv/python/cpython-3.11.14-linux-x86_64-gnu/bin/python3.11
```

The setup script checks `if [ ! -d "venv" ]` and skips venv creation if the directory exists, even if it's broken.

**Fix:**
```bash
# Remove the stale venv directory
rm -rf /path/to/SciAgent-toolkit/mcp_servers/pal/venv

# Re-run setup
./scripts/setup-ai.sh --force
```

**Prevention:** The `venv/` pattern is in `.gitignore`, so these directories should not be committed. If you're sharing the toolkit across environments (e.g., via a mounted volume or synced folder), ensure venv directories are excluded.

**Verification:**
```bash
# Check if venv symlinks are valid
ls -la mcp_servers/pal/venv/bin/python*
# If the target path doesn't exist on this machine, delete and recreate
```

---

### Issue: PAL setup reports failure but actually works

**Symptoms:**
- `setup-ai.sh` shows `[ERROR] Failed to install pal-mcp-server` and `[ERROR] PAL MCP setup failed`
- However, PAL MCP is actually functional when tested

**Cause:** The error occurs when `uv` tries to inspect an existing venv directory during the setup script. If the venv was created successfully moments before (in a previous step), `uv` may report a misleading error about the Python interpreter while the installation actually completed.

**Diagnosis:**
```bash
# Test if PAL actually works
/path/to/SciAgent-toolkit/mcp_servers/pal/venv/bin/pal-mcp-server --help

# Check if Python in venv works
/path/to/SciAgent-toolkit/mcp_servers/pal/venv/bin/python --version
```

**Resolution:** If the above tests pass, PAL is installed correctly despite the error message. You can safely ignore the error and verify PAL works in Claude Code with `/mcp`.

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

1. **Use `source_env.sh`** (recommended, persists across rebuilds):
   ```bash
   bash .devcontainer/scripts/source_env.sh
   # Restart Claude Code
   ```
   See [Environment Variable Loading](#environment-variable-loading-common-root-cause) for full setup.

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

## MCP Addon Tools Denied in Subagents

### Symptoms

- MCP addon tools (like `nb_search`, `nb_read_cell`) are denied when called from Claude Code subagents (Task tool)
- The main conversation can use the tools (after interactive approval), but subagents auto-deny them
- Only tools that the user has manually approved at least once work; the rest are silently blocked

### Cause

`manage-addon.sh` (before the fix) added MCP server entries to `.mcp.json` but never added `mcp__<server>__<tool>` entries to `permissions.allow` in `.claude/settings.local.json`. Claude Code subagents run autonomously and can't prompt interactively, so unapproved MCP tools are auto-denied.

### What was changed (4 files)

| File | Change |
|------|--------|
| `templates/mcp-addons/notebook-tools.addon.json` | Added `tool_permissions` array with all 11 tool names |
| `templates/mcp-addons/jupyter.addon.json` | Added empty `tool_permissions` array (ready for future population) |
| `scripts/manage-addon.sh` | `update_settings_local()` now reads `tool_permissions` from all enabled addon templates, strips stale `mcp__*` entries, and merges current addon permissions into `permissions.allow` |
| `scripts/switch-mcp-profile.sh` | Settings generation block now reads `active-addons.json` + addon templates to include tool permissions alongside default permissions |

### How the permission flow works

**Addon template** (e.g. `notebook-tools.addon.json`):
```json
{
  "_meta": {
    "tool_permissions": [
      "mcp__notebook-tools__nb_overview",
      "mcp__notebook-tools__nb_read_cell",
      "..."
    ]
  }
}
```

**Enable flow:** `cmd_enable()` → `update_settings_local()` reads all enabled addons' `tool_permissions`, strips old `mcp__*` entries from `permissions.allow`, appends current ones (deduplicated).

**Disable flow:** `cmd_disable()` sets `enabled: false` first → `update_settings_local()` rebuilds from only enabled addons, so disabled addon's permissions are removed.

**Profile switch flow:** `switch-mcp-profile.sh` reads `active-addons.json` + addon templates to include addon tool permissions in the fresh `settings.local.json` it writes, so permissions survive profile changes.

### Verification

After enabling an addon:
```bash
# Check permissions were added
jq '.permissions.allow' .claude/settings.local.json | grep mcp__notebook-tools

# Test: disable should remove them
./scripts/manage-addon.sh disable notebook-tools --project-dir .
jq '.permissions.allow' .claude/settings.local.json | grep mcp__notebook-tools  # should return empty

# Test: enable should add them back
./scripts/manage-addon.sh enable notebook-tools --project-dir .
jq '.permissions.allow' .claude/settings.local.json | grep mcp__notebook-tools  # should show 11 entries

# Test: profile switch preserves addon permissions
./scripts/switch-mcp-profile.sh coding
jq '.permissions.allow' .claude/settings.local.json | grep mcp__notebook-tools  # should still show 11 entries
```

### For addon authors

When creating a new addon template, always include a `tool_permissions` array in `_meta` listing every `mcp__<server>__<tool>` that the addon provides. Without this, subagents will be unable to use the addon's tools:

```json
{
  "_meta": {
    "name": "my-addon",
    "tool_permissions": [
      "mcp__my-addon__tool_one",
      "mcp__my-addon__tool_two"
    ]
  }
}
```

---

## Getting More Help

- [INSTALLATION.md](INSTALLATION.md) - Full installation guide
- [CONFIGURATION.md](CONFIGURATION.md) - Advanced configuration
- [ISSUES.md](ISSUES.md) - Known architectural issues
- [FAQ.md](FAQ.md) - Frequently asked questions