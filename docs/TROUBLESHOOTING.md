# Troubleshooting

## ToolUniverse Issues

### Issue: "Connection failed for 'tooluniverse'" or "ModuleNotFoundError"

**Symptoms:**
- Gemini/Claude shows `✕ Error during discovery for server 'tooluniverse'`
- Logs show `ModuleNotFoundError: No module named 'tooluniverse'`
- `uv` fails to detect the environment or defaults to the system Python.

**Cause:**
- The `tooluniverse-env` virtual environment may be corrupted or created with incorrect permissions.
- `uv run` might be picking up the wrong Python interpreter.

**Fix:**

1.  **Reinstall the environment manually:**

    ```bash
    # Remove existing environment
    rm -rf 01_scripts/SciAgent-toolkit/scripts/tooluniverse-env
    rm -rf tooluniverse-env

    # Create new venv in project root
    uv venv tooluniverse-env

    # Install tooluniverse using the specific python binary to avoid permission issues
    uv pip install --python tooluniverse-env/bin/python tooluniverse
    ```

2.  **Verify the installation:**

    ```bash
    # Check if the executable works directly
    ./tooluniverse-env/bin/tooluniverse-smcp-stdio --help
    ```

3.  **Update `.mcp.json` to use the direct executable:**

    Edit `.mcp.json` and change the `tooluniverse` section to point directly to the binary instead of using `uv run`.

    *From:*
    ```json
    "tooluniverse": {
      "command": "uv",
      "args": ["--directory", "...", "run", "tooluniverse-smcp-stdio", ...]
    }
    ```

    *To:*
    ```json
    "tooluniverse": {
      "command": "/absolute/path/to/tooluniverse-env/bin/tooluniverse-smcp-stdio",
      "args": ["--include-tools", "..."]
    }
    ```

4.  **Restart the CLI:**
    Exit and restart Gemini/Claude.

---

### Issue: ToolUniverse installation fails with "Permission denied"

**Symptoms:**
- `setup_tooluniverse.sh` fails.
- Error: `error: failed to remove file ... Permission denied`

**Fix:**
- Use the manual installation steps above, ensuring you use `uv pip install --python ...` to explicitly target the venv's interpreter.

---

### Issue: ToolUniverse fails in "full" or "research-full" profiles

- Switching to `full` or `research-full` profile using `switch_mcp` command results in broken ToolUniverse functionality.
- Errors regarding invalid arguments passed to `tooluniverse-smcp-stdio` (e.g., `unrecognized arguments: --directory`).
- This happens when the `.mcp.json` configuration mixes the direct binary path with `uv` arguments.

**Fix:**

- Use the `research-lite` profile which uses a cleaner configuration:
  ```bash
  ./scripts/switch-mcp-profile.sh research-lite
  ```
- Or manually edit `.mcp.json` to remove `--directory ... run ...` arguments from the `tooluniverse` args list, keeping only the tool-specific arguments (e.g., `--include-tools ...`).

---

### Issue: Gemini: ToolUniverse tools not appearing

**Symptoms:**
- After switching to a research profile, `find_tools` returns an empty list `[]` or only default tools.
- `.gemini/settings.json` appears to have the tools listed in the `args` array.

**Cause:**
- The `--include-tools` argument in the JSON configuration was passed as a single comma-separated string (e.g., `"Tool1,Tool2"`).
- The Python `argparse` library (used by ToolUniverse) expects separate arguments for `nargs='+'` or `nargs='*'` when receiving a list from `subprocess`.

**Fix:**
- Ensure that in `.gemini/settings.json` (and the `templates/gemini-profiles/research.json`), each tool name is a separate string in the `args` array.

**Incorrect:**
```json
"args": ["--include-tools", "Tool1,Tool2"]
```

**Correct:**
```json
"args": ["--include-tools", "Tool1", "Tool2"]
```

---

## Codex CLI Issues

### Issue: Installation fails with `npm error code EACCES`

**Symptoms:**
- Running `npm install -g @openai/codex` or `install_codex.sh` fails.
- Error log shows `Error: EACCES: permission denied, mkdir '/opt/nvm/versions/node/...'`

**Cause:**
- The global Node.js installation directory (often `/opt/nvm/...` in dev containers) is read-only for the current user.
- npm attempts to write to this system directory by default.

**Fix:**
- Configure npm to use a user-writable directory (e.g., `~/.npm-global`) for global packages.

```bash
# 1. Create directory
mkdir -p "$HOME/.npm-global"

# 2. Configure npm
npm config set prefix "$HOME/.npm-global"

# 3. Add to PATH (add to ~/.bashrc for persistence)
export PATH="$HOME/.npm-global/bin:$PATH"
echo 'export PATH="$HOME/.npm-global/bin:$PATH"' >> "$HOME/.bashrc"

# 4. Retry installation
npm install -g @openai/codex
```

### Issue: MCP Handshake Failure (Stdio Noise)

**Symptoms:**
- Codex CLI starts but shows `⚠ MCP client for 'tooluniverse' failed to start`.
- Error: `handshaking with MCP server failed: connection closed: initialize response` or `JSON parse error`.
- Claude/Gemini may connect fine to the same server.

**Cause:**
- The Codex CLI enforces strict [MCP over Stdio](https://modelcontextprotocol.io/docs/concepts/transports#stdio) compliance.
- The `tooluniverse` server (and some python scripts) prints initialization logs (e.g., `🚀 Starting ToolUniverse...`) to **stdout**.
- This extra text corrupts the JSON-RPC messages required by Codex.

**Fix / Workaround:**
- **Patch the Server:** Modify the MCP server to print all logs to `stderr` only.
- **Use Quiet Mode:** If the server supports it, use a flag to suppress stdout logs.
- **Exclude from Codex:** In `~/.codex/config.toml`, remove the offending server until it is patched.
- **Use Direct Binary:** Sometimes using the direct binary (e.g., `tooluniverse-env/bin/python`) instead of `uv run` helps if `uv` itself is adding noise.

---