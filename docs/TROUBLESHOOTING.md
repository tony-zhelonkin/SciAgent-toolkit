# Troubleshooting Guide

Common issues and solutions for the SciAgent Toolkit.

## Quick Diagnostics

Run these commands for automated diagnostics:

```bash
# Check installation status
claude doctor

# Verify MCP configuration
python3 -m json.tool .mcp.json

# Test individual servers
npx -y @modelcontextprotocol/server-sequential-thinking --help
./test_tooluniverse.sh
uvx --from git+https://github.com/oraios/serena serena --help
```

---

## Critical Bug Fixes (v1.1.0+)

**If you installed before v1.1.0, please update:**

```bash
cd SciAgent-toolkit
git pull origin main

# Re-run installation to get fixes
./scripts/setup_mcp_infrastructure.sh --mcp-only
```

**What was fixed:**
- ✅ ToolUniverse virtual environment creation
- ✅ Serena SSH/HTTPS authentication
- ✅ Automated MCP configuration
- ✅ Command name auto-detection

---

## Installation Issues

### Issue: Claude Code or Codex not found after installation

**Symptoms**:
```bash
bash: claude: command not found
bash: codex: command not found
```

**Solution 1: Reload shell configuration**
```bash
source ~/.bashrc  # or source ~/.zshrc

# Or restart your terminal
```

**Solution 2: Check PATH**
```bash
# For Claude Code
echo $PATH | grep ".local/bin"

# Manually add if needed
export PATH="$HOME/.local/bin:$PATH"

# Add to ~/.bashrc permanently
echo 'export PATH="$HOME/.local/bin:$PATH"' >> ~/.bashrc
```

### Issue: Permission denied errors

**Symptoms**:
```bash
./setup_mcp_infrastructure.sh: Permission denied
```

**Solution**:
```bash
# Make scripts executable
chmod +x scripts/setup_mcp_infrastructure.sh
chmod +x scripts/install_*.sh
chmod +x scripts/mcp_servers/*.sh
```

### Issue: Python version too old

**Symptoms**:
```bash
Error: Python 3.10+ required
```

**Solution (macOS)**:
```bash
# Install Python 3.10+ via Homebrew
brew install python@3.11

# Verify
python3 --version
```

**Solution (Linux)**:
```bash
# Ubuntu/Debian
sudo apt update
sudo apt install python3.11

# Fedora
sudo dnf install python3.11
```

### Issue: Node.js not found

**Symptoms**:
```bash
Error: Node.js 18+ required
```

**Solution (macOS)**:
```bash
brew install node
```

**Solution (Linux - Ubuntu/Debian)**:
```bash
curl -fsSL https://deb.nodesource.com/setup_20.x | sudo -E bash -
sudo apt-get install -y nodejs
```

### Issue: ToolUniverse Virtual Environment Error

**Symptoms**:
```bash
error: No virtual environment found
error: run `uv venv` to create an environment
[ERROR] Failed to install ToolUniverse
```

**Cause**: Bug in v1.0.x - missing `uv venv` command before installation

**Solution (v1.1.0+)**:
```bash
# Update to latest version
git pull origin main
./scripts/mcp_servers/setup_tooluniverse.sh
```

**Manual Fix**:
```bash
cd SciAgent-toolkit
uv venv tooluniverse-env
uv --directory tooluniverse-env pip install tooluniverse

# Verify installation
uv --directory tooluniverse-env run python -c "import tooluniverse; print(tooluniverse.__version__)"
```

### Issue: Serena SSH Authentication Failure

**Symptoms**:
```bash
Permission denied (publickey)
ssh_askpass: exec(/path/to/ssh-askpass): No such file or directory
git@github.com: Permission denied (publickey)
failed to fetch commit
```

**Cause**: Git defaulting to SSH instead of HTTPS for GitHub

**Solution (v1.1.0+ - Automatic)**:
```bash
git pull origin main
./scripts/mcp_servers/setup_serena.sh
```

**Manual Fix**:
```bash
# Force Git to use HTTPS
git config --global url."https://github.com/".insteadOf git@github.com:

# Test Serena installation
uvx --from git+https://github.com/oraios/serena serena --help
```

**Alternative - Set up SSH keys** (if you prefer SSH):
```bash
# Generate SSH key
ssh-keygen -t ed25519 -C "your_email@example.com"

# Add to GitHub
cat ~/.ssh/id_ed25519.pub
# Copy and add to https://github.com/settings/keys
```

### Issue: MCP Configuration Not Created

**Symptoms**:
- No `.mcp.json` file in project directory
- `/mcp` shows "No MCP servers configured"
- Manual configuration required

**Cause**: Missing in v1.0.x - config merger not implemented

**Solution (v1.1.0+)**:
```bash
# Run configuration script
./scripts/configure_mcp_servers.sh

# Verify it worked
ls -la .mcp.json
python3 -m json.tool .mcp.json

# Start Claude and check
claude
/mcp  # Should show all servers
```

**Manual Configuration**:
```bash
# Add servers one by one
claude mcp add sequential-thinking --scope local -- \
  npx -y @modelcontextprotocol/server-sequential-thinking

claude mcp add tooluniverse --scope local -- \
  uv --directory "$(pwd)/tooluniverse-env" run tooluniverse-mcp

claude mcp add serena --scope local -- \
  uvx --from git+https://github.com/oraios/serena serena start-mcp-server
```

### Issue: ToolUniverse Command Not Found

**Symptoms**:
```bash
tooluniverse-mcp: command not found
# OR
tooluniverse-smcp-stdio: command not found
```

**Cause**: Command name varies between ToolUniverse versions

**Solution (v1.1.0+ - Auto-detected)**:
The latest script automatically detects the correct command name.

**Manual Check**:
```bash
# Try both possible commands
uv --directory ./tooluniverse-env run tooluniverse-mcp --version
# OR
uv --directory ./tooluniverse-env run tooluniverse-smcp-stdio --version

# Update .mcp.json with the working command
```

## MCP Server Issues

### Issue: MCP servers not loading

**Symptoms**:
- `/mcp` shows no servers
- Server status: "disconnected"

**Diagnostic Steps**:

**Step 1: Check configuration syntax**
```bash
# For Claude Code
python3 -m json.tool .mcp.json

# For Codex CLI
cat ~/.codex/config.toml
```

**Step 2: Run diagnostics**
```bash
# For Claude Code
claude doctor
```

**Step 3: Test individual servers**
```bash
# Test Serena
uvx --from git+https://github.com/oraios/serena serena --help

# Test Sequential Thinking
npx -y @modelcontextprotocol/server-sequential-thinking --help

# Test ToolUniverse
uv --directory ./tooluniverse-env run tooluniverse-smcp-stdio --help
```

**Step 4: Check paths are absolute**
```json
{
  "args": [
    "--directory",
    "/absolute/path/to/tooluniverse-env"  // Must be absolute, not relative
  ]
}
```

**Solution**: Fix configuration and restart Claude/Codex.

### Issue: ToolUniverse tools not executing

**Symptoms**:
- Tools appear in `/mcp` but fail when used
- Timeout errors

**Diagnostic Steps**:

**Step 1: Verify Python version**
```bash
python3 --version  # Should be 3.10+
```

**Step 2: Check uv installation**
```bash
uv --version
```

**Step 3: Test ToolUniverse directly**
```bash
./test_tooluniverse.sh
```

**Step 4: Check for Azure OpenAI errors (if configured)**
```bash
echo $AZURE_OPENAI_API_KEY
echo $AZURE_OPENAI_ENDPOINT
```

**Solution**:
```bash
# Reinstall ToolUniverse
./scripts/mcp_servers/setup_tooluniverse.sh

# If Azure OpenAI issues, unset variables
unset AZURE_OPENAI_API_KEY
unset AZURE_OPENAI_ENDPOINT

# Re-run setup
./scripts/mcp_servers/setup_tooluniverse.sh
```

### Issue: Serena MCP not starting

**Symptoms**:
- Serena appears disconnected
- "Failed to start Serena" errors

**Solution**:
```bash
# Check uvx installation
uvx --version

# Reinstall Serena
./scripts/mcp_servers/setup_serena.sh

# Test manually
uvx --from git+https://github.com/oraios/serena serena start-mcp-server --help
```

### Issue: Sequential Thinking MCP not starting

**Symptoms**:
- Sequential Thinking disconnected
- npx errors

**Solution**:
```bash
# Check Node.js and npx
node --version
npx --version

# Reinstall Node.js if needed (see above)

# Test manually
npx -y @modelcontextprotocol/server-sequential-thinking --help

# Reinstall
./scripts/mcp_servers/setup_sequential_thinking.sh
```

## PubMed Plugin Issues

### Issue: PubMed not available in Claude Code

**Symptoms**:
- `/mcp` doesn't show PubMed
- "PubMed plugin not found"

**Solution**:
```bash
# PubMed requires manual installation
claude

# In Claude Code:
/plugin marketplace add anthropics/life-sciences
/plugin install pubmed@life-sciences

# Restart Claude Code (Ctrl+D, then 'claude' again)

# Verify
/mcp
```

**Note**: PubMed is a Claude Code exclusive feature and not available in Codex CLI.

## Configuration Issues

### Issue: Tool filtering not working

**Symptoms**:
- `--include-tools` is ignored
- All 600+ tools still appear

**Solution**:

**Step 1: Verify tool names**
```bash
# List all available tools
./test_tooluniverse.sh

# Or check in Claude/Codex
"List all available ToolUniverse tools"
```

**Step 2: Check tool name spelling**
Tool names are case-sensitive. Common mistakes:
- `europepmc_search` ✗ → `EuropePMC_search_articles` ✓
- `chembl_search` ✗ → `ChEMBL_search_similar_molecules` ✓

**Step 3: Verify configuration syntax**
```json
{
  "args": [
    "--include-tools",
    "tool1,tool2,tool3"  // Comma-separated, no spaces
  ]
}
```

### Issue: Multiple instances interfering

**Symptoms**:
- Duplicate tools appearing
- Unexpected tool availability

**Solution**:

Ensure each instance has a unique name:
```json
{
  "mcpServers": {
    "tooluniverse-literature": { ... },    // ✓ Unique name
    "tooluniverse-drugs": { ... },         // ✓ Unique name
    // NOT:
    "tooluniverse": { ... },               // ✗ Duplicate name
    "tooluniverse": { ... }                // ✗ Duplicate name
  }
}
```

### Issue: Changes to .mcp.json not taking effect

**Symptoms**:
- Modified configuration not reflected
- Old servers still loading

**Solution**:
```bash
# Validate JSON syntax
python3 -m json.tool .mcp.json

# Completely restart Claude/Codex
# Press Ctrl+D to exit
claude  # Start again

# Verify changes
/mcp
```

## Performance Issues

### Issue: Slow MCP server startup

**Symptoms**:
- Long delays when starting Claude/Codex
- Timeout errors

**Solution 1: Increase timeout (Codex CLI)**
```toml
[mcp_servers.tooluniverse]
command = "uv"
args = [ ... ]
startup_timeout_sec = 120  # Increase from default 60
```

**Solution 2: Reduce tool count**
```json
{
  "args": [
    "--include-tools",
    "fewer,specific,tools"  // Less overhead
  ]
}
```

### Issue: High memory usage

**Symptoms**:
- System slowing down
- Out of memory errors

**Solution**: Use tool filtering to reduce context:
```json
{
  "args": [
    "--exclude-tool-types",
    "unused_category1,unused_category2"
  ]
}
```

## Platform-Specific Issues

### macOS: "xcrun: error" during installation

**Symptoms**:
```bash
xcrun: error: invalid active developer path
```

**Solution**:
```bash
# Install Xcode Command Line Tools
xcode-select --install
```

### Linux: Package installation failures

**Symptoms**:
```bash
E: Unable to locate package python3.11
```

**Solution (Ubuntu/Debian)**:
```bash
sudo apt update
sudo apt upgrade

# Add deadsnakes PPA for newer Python
sudo add-apt-repository ppa:deadsnakes/ppa
sudo apt update
sudo apt install python3.11
```

### WSL: PATH not updating

**Symptoms**:
- Commands not found in WSL
- PATH changes don't persist

**Solution**:
```bash
# Check which shell you're using
echo $SHELL

# Edit correct config file
# For bash:
nano ~/.bashrc

# For zsh:
nano ~/.zshrc

# Add PATH permanently
export PATH="$HOME/.local/bin:$PATH"

# Reload
source ~/.bashrc  # or source ~/.zshrc
```

## Data/API Issues

### Issue: API rate limiting

**Symptoms**:
- "Rate limit exceeded" errors
- Failed API calls

**Solution**:
- Wait before retrying
- For PubMed: Use NCBI API key (see [PubMed E-utilities](https://www.ncbi.nlm.nih.gov/books/NBK25497/))
- For other APIs: Check provider documentation

### Issue: Network connectivity errors

**Symptoms**:
- Timeout errors
- "Failed to connect" messages

**Solution**:
```bash
# Test internet connectivity
ping google.com

# Check firewall settings
# Ensure ports are open for:
# - GitHub (git+https://)
# - npm registry
# - Python package index (PyPI)
```

## Debugging Tips

### Enable Verbose Logging

**Claude Code**:
```bash
# Run with debug output
DEBUG=* claude
```

**Codex CLI**:
```bash
# Run with debug output
DEBUG=true codex
```

### Check System Logs

```bash
# macOS
tail -f ~/Library/Logs/Claude/mcp.log

# Linux
journalctl -f | grep -i mcp
```

### Test Components Individually

```bash
# Test Python
python3 -c "import sys; print(sys.version)"

# Test Node.js
node -e "console.log(process.version)"

# Test uv
uv --version

# Test uvx
uvx --version

# Test npx
npx --version
```

## Getting Help

If you're still experiencing issues:

1. **Check existing issues**: Review [FAQ.md](FAQ.md)
2. **Check logs**: Look for error messages in terminal output
3. **Minimal reproduction**: Try installing on a clean system
4. **Gather information**:
   - Operating system and version
   - Python version: `python3 --version`
   - Node.js version: `node --version`
   - Error messages (full output)
5. **Report issues**: Create a GitHub issue with details above

## Common Error Messages

### "command not found: claude"
→ See [Claude Code or Codex not found after installation](#issue-claude-code-or-codex-not-found-after-installation)

### "Python version 3.10+ required"
→ See [Python version too old](#issue-python-version-too-old)

### "MCP server failed to start"
→ See [MCP servers not loading](#issue-mcp-servers-not-loading)

### "Tool execution timeout"
→ See [Slow MCP server startup](#issue-slow-mcp-server-startup)

### "Invalid JSON in .mcp.json"
→ Run: `python3 -m json.tool .mcp.json`

### "Permission denied"
→ See [Permission denied errors](#issue-permission-denied-errors)

---

For more information:
- [INSTALLATION.md](INSTALLATION.md) - Installation guide
- [CONFIGURATION.md](CONFIGURATION.md) - Configuration options
- [FAQ.md](FAQ.md) - Frequently asked questions
