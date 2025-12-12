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