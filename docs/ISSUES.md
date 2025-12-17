# SciAgent-toolkit Architecture Issues

**Created**: 2025-12-16
**Version**: Current (as of v1.x)
**Status**: Active issue tracking for MCP configuration system

This document catalogs architectural issues, design inconsistencies, and technical debt in the SciAgent-toolkit MCP configuration system. Issues are prioritized by severity and impact.

---

## Critical Issues

### ISSUE-001: PAL Server Definition Inconsistency Across Profiles

**Severity**: Critical
**Impact**: PAL MCP fails on initial setup with `--force` flag
**Reproducible**: Yes
**Status**: ✅ Resolved (2025-12-16) - All profiles now use direct uvx approach

#### Description

PAL MCP server is defined using **two incompatible approaches** across different profile templates:

| Profile | PAL Command | Requires Setup |
|---------|-------------|----------------|
| `coding.mcp.json` | `${TOOLKIT_ROOT}/mcp_servers/pal/start-pal.py` | Yes (`setup_pal.sh`) |
| `codebase.mcp.json` | `uvx --from git+...` (direct) | No |
| `full.mcp.json` | `uvx --from git+...` (direct) | No |
| `hybrid-research.mcp.json` | `uvx --from git+...` (direct) | No |

#### Root Cause

1. `coding.mcp.json` expects a wrapper script at `${TOOLKIT_ROOT}/mcp_servers/pal/start-pal.py`
2. This wrapper is **only created** when `setup_pal.sh` is run
3. The wrapper's purpose is to load `.env` files and inject API keys at runtime
4. Direct `uvx` invocation (other profiles) requires environment variables to be **already set** in the shell or passed via `env:` block

#### Failure Scenario

```bash
./setup-ai.sh --force
# setup_mcp_infrastructure.sh runs
# PAL is installed (setup_pal.sh creates wrapper)
# configure_mcp_servers.sh applies "coding" profile → uses wrapper → WORKS

./switch-mcp-profile.sh hybrid
# Applies hybrid-research.mcp.json → uses direct uvx → FAILS
# Because env vars (GEMINI_API_KEY, OPENAI_API_KEY) are not substituted

claude
# PAL shows "Missing environment variables: GEMINI_API_KEY, OPENAI_API_KEY"
```

#### Recommended Fix

**Option A (Consistent Wrapper)**: Standardize on the wrapper approach for all profiles:
```json
{
  "pal": {
    "type": "stdio",
    "command": "${TOOLKIT_ROOT}/mcp_servers/pal/start-pal.py",
    "args": []
  }
}
```
- Wrapper handles environment loading consistently
- Requires `setup_pal.sh` to have been run

**Option B (Direct uvx + Substitution)**: Use direct uvx with proper env var substitution:
```json
{
  "pal": {
    "type": "stdio",
    "command": "uvx",
    "args": ["--from", "git+https://github.com/BeehiveInnovations/pal-mcp-server.git", "pal-mcp-server"],
    "env": {
      "GEMINI_API_KEY": "${GEMINI_API_KEY}",
      "OPENAI_API_KEY": "${OPENAI_API_KEY}"
    }
  }
}
```
- Requires `switch-mcp-profile.sh` to substitute `${GEMINI_API_KEY}` placeholders
- Currently **NOT implemented** in the substitution logic

---

### ISSUE-002: Environment Variable Substitution Incomplete

**Severity**: Critical
**Impact**: API keys not injected into .mcp.json when using direct uvx profiles
**Reproducible**: Yes

#### Description

The Python substitution logic in `switch-mcp-profile.sh` only handles:
- `${TOOLUNIVERSE_ENV}` → path substitution
- `${TOOLKIT_ROOT}` → path substitution

It does **NOT** handle:
- `${GEMINI_API_KEY}` → API key substitution
- `${OPENAI_API_KEY}` → API key substitution
- `${CONTEXT7_API_KEY}` → API key substitution (for Context7)

#### Affected Code

```python
# switch-mcp-profile.sh, lines 115-134
content = content.replace('${TOOLUNIVERSE_ENV}', tool_env)
content = content.replace('${TOOLKIT_ROOT}', toolkit_root)
# Missing:
# content = content.replace('${GEMINI_API_KEY}', gemini_key)
# content = content.replace('${OPENAI_API_KEY}', openai_key)
```

#### Impact

When `hybrid-research.mcp.json` is applied:
```json
{
  "env": {
    "GEMINI_API_KEY": "${GEMINI_API_KEY}",
    "OPENAI_API_KEY": "${OPENAI_API_KEY}"
  }
}
```

The literal strings `${GEMINI_API_KEY}` are written to `.mcp.json` instead of the actual values.

#### Recommended Fix

Update the Python substitution block to include API keys:

```python
# In switch-mcp-profile.sh
content = content.replace('${TOOLUNIVERSE_ENV}', tool_env)
content = content.replace('${TOOLKIT_ROOT}', toolkit_root)
content = content.replace('${GEMINI_API_KEY}', os.environ.get('GEMINI_API_KEY', ''))
content = content.replace('${OPENAI_API_KEY}', os.environ.get('OPENAI_API_KEY', ''))
content = content.replace('${CONTEXT7_API_KEY}', os.environ.get('CONTEXT7_API_KEY', ''))
```

---

## High Priority Issues

### ISSUE-003: Two-Stage Environment Binding Creates Stale Configuration

**Severity**: High
**Impact**: Configuration drift between Claude/Gemini and PAL
**Root Cause**: Architectural

#### Description

Environment variables are bound at different times:

| Tool | Binding Time | Source |
|------|--------------|--------|
| Claude Code (.mcp.json) | Profile switch | Static substitution |
| Gemini (.gemini/settings.json) | Profile switch | Static substitution |
| PAL wrapper (start-pal.py) | MCP startup | Dynamic from .env |

#### Problem

1. User edits `.env` to add/change `GEMINI_API_KEY`
2. PAL wrapper picks up new key immediately (runtime loading)
3. Claude/Gemini configs still have old/missing key (static substitution)
4. User must re-run `switch-mcp-profile.sh` to update Claude/Gemini

#### Recommended Fix

**Option A (Document)**: Clearly document that profile switch must be re-run after `.env` changes

**Option B (Runtime Loading)**: Create wrappers for all servers that load `.env` at runtime

**Option C (Symlink/Include)**: Use Claude Code's ability to read env vars from shell environment instead of hardcoding in `.mcp.json`

---

### ISSUE-004: Profile Templates Missing Installation Validation

**Severity**: High
**Impact**: Profiles reference uninstalled servers
**Reproducible**: Yes

#### Description

`switch-mcp-profile.sh` applies profile templates without checking if referenced servers are installed:

- `codebase.mcp.json` references `serena` → may not be installed (skip-serena flag)
- `full.mcp.json` references `tooluniverse` → may not be installed (skip-tooluniverse flag)
- `coding.mcp.json` references `${TOOLKIT_ROOT}/mcp_servers/pal/start-pal.py` → may not exist

#### Current Behavior

Profile applies successfully, but Claude Code fails at startup:
```
[pal] Status: ✘ failed
Command: /path/to/start-pal.py (file not found)
```

#### Recommended Fix

Add validation before applying profile:

```bash
# In switch-mcp-profile.sh, after reading profile template
validate_profile_requirements() {
    local profile_content="$1"

    # Check for serena
    if echo "$profile_content" | grep -q '"serena"'; then
        command -v uvx &>/dev/null || {
            log_error "Profile requires serena but uvx not found"
            return 1
        }
    fi

    # Check for PAL wrapper
    if echo "$profile_content" | grep -q 'start-pal.py'; then
        [ -f "${TOOLKIT_ROOT}/mcp_servers/pal/start-pal.py" ] || {
            log_error "Profile requires PAL wrapper but it's not installed"
            log_info "Run: ./scripts/mcp_servers/setup_pal.sh"
            return 1
        }
    fi

    # ... similar checks for other servers
}
```

---

### ISSUE-005: ToolUniverse Tool Count Inconsistency

**Severity**: Medium
**Impact**: Token estimates inaccurate, potential context overflow
**Reproducible**: Yes
**Status**: ✅ Resolved - research-full.mcp.json now has --include-tools with 14 tools

#### Description

Different profiles specify different tool subsets, but one profile is misconfigured:

| Profile | Stated Tokens | Tool Specification |
|---------|---------------|-------------------|
| `research-lite.mcp.json` | ~30k | `--include-tools` with 6 tools |
| `research-full.mcp.json` | ~50k | `--hook-type SummarizationHook` **but no --include-tools** |
| `full.mcp.json` | ~100k | `--include-tools` with 14 tools |

#### Problem

`research-full.mcp.json` without `--include-tools` loads **all 600+ tools**, consuming 500k+ tokens, vastly exceeding the stated ~50k estimate and the 200k context limit.

#### Current research-full.mcp.json

```json
{
  "tooluniverse": {
    "type": "stdio",
    "command": "${TOOLUNIVERSE_ENV}/bin/python",
    "args": [
      "-m", "tooluniverse.smcp_server",
      "--transport", "stdio",
      "--hook-type", "SummarizationHook"
    ]
  }
}
```

#### Recommended Fix

Add `--include-tools` to `research-full.mcp.json`:

```json
{
  "args": [
    "-m", "tooluniverse.smcp_server",
    "--transport", "stdio",
    "--hook-type", "SummarizationHook",
    "--include-tools", "EuropePMC_search_articles,openalex_literature_search,pubmed_search,semantic_scholar_search,ChEMBL_search_similar_molecules,ChEMBL_get_molecule,ChEMBL_search_targets,DrugBank_search,UniProt_search_proteins,UniProt_get_protein,search_clinical_trials,get_clinical_trial_details,FDA_adverse_events,FDA_drug_labels"
  ]
}
```

---

## Medium Priority Issues

### ISSUE-006: Multiple ToolUniverse Installation Locations

**Severity**: Medium
**Impact**: Confusion about which venv is used
**Root Cause**: Design

#### Description

ToolUniverse venv can exist in two locations:

1. `${PROJECT_DIR}/tooluniverse-env/` (preferred, project-local)
2. `${SCRIPT_DIR}/tooluniverse-env/` (fallback, toolkit-local)

The fallback was added for cases where setup wasn't run from project directory, but it creates confusion:
- Which version is being used?
- Are they in sync?
- Which one gets updates?

#### Recommended Fix

Remove the fallback location; always require project-local installation:

```bash
# In switch-mcp-profile.sh
TOOLUNIVERSE_ENV=""
if [ -d "${PROJECT_DIR}/tooluniverse-env" ]; then
    TOOLUNIVERSE_ENV="${PROJECT_DIR}/tooluniverse-env"
else
    log_warn "ToolUniverse not found at ${PROJECT_DIR}/tooluniverse-env"
    log_info "Run: ./scripts/mcp_servers/setup_tooluniverse.sh"
fi
# Remove the fallback: elif [ -d "${SCRIPT_DIR}/tooluniverse-env" ]
```

---

### ISSUE-007: Inconsistent Error Handling in Python Blocks

**Severity**: Medium
**Impact**: Silent failures during profile switch
**Example**: ISSUE-001 root cause (fixed 2025-12-16)

#### Description

Python code embedded in bash scripts has inconsistent error handling:

1. Some blocks use `try/except` with `sys.exit(1)` on error
2. Some blocks silently ignore errors
3. Some blocks have syntax errors (fixed example: missing `import stat`, wrong indentation)

#### Examples

**Good pattern** (in Gemini configuration, after fix):
```python
try:
    with open(settings_path, 'w') as f:
        f.write(content)
    os.chmod(settings_path, stat.S_IRUSR | stat.S_IWUSR)
except Exception as e:
    sys.stderr.write(f'Error writing settings: {e}\n')
    sys.exit(1)
```

**Bad pattern** (in PAL wrapper .env loading):
```python
except Exception as e:
    # Silently ignore errors reading .env, just proceed
    pass
```

#### Recommended Fix

Standardize on explicit error handling with clear messages:

```python
except Exception as e:
    sys.stderr.write(f'Warning: Could not load {env_file}: {e}\n')
    # Continue with default values or fail explicitly depending on severity
```

---

### ISSUE-008: Bash Script Code Comments as Documentation

**Severity**: Low
**Impact**: Orphaned comments create confusion
**Example**: Line 143 in switch-mcp-profile.sh

#### Description

Leftover comment from refactoring:
```bash
# ... (previous code)
```

This comment doesn't describe what code does or should do - it's a placeholder that was never cleaned up.

#### Recommended Fix

Remove orphaned comments or replace with meaningful ones:
```bash
# --------------------------
# Gemini Configuration
# --------------------------
```

---

## Architectural Recommendations

### Recommendation 1: Unified Server Definition Format

Create a single, consistent server definition approach:

```json
{
  "pal": {
    "type": "stdio",
    "command": "python3",
    "args": ["-m", "pal_mcp_server"],
    "cwd": "${TOOLKIT_ROOT}/mcp_servers/pal/venv",
    "env_files": ["${PROJECT_DIR}/.env", "${PROJECT_DIR}/.devcontainer/.env"]
  }
}
```

This would require either:
- Claude Code to support `env_files` directive (feature request)
- Or a universal wrapper pattern for all servers

### Recommendation 2: Validate-on-Switch Pattern

Profile switching should:
1. Parse profile to identify required servers
2. Check each server is installed/available
3. Warn about missing servers (non-fatal)
4. Substitute all placeholders including API keys
5. Write configuration
6. Verify configuration is valid JSON

### Recommendation 3: Single Source of Truth for Environment

Choose ONE approach and use consistently:

| Approach | Pros | Cons |
|----------|------|------|
| Static substitution | Simple, predictable | Stale after .env changes |
| Runtime loading (wrapper) | Always current | Extra complexity, startup cost |
| Shell environment | Standard Unix pattern | Requires user to export vars |

Current codebase mixes all three, causing confusion.

### Recommendation 4: Profile Dependency Declaration

Add metadata to profiles declaring their requirements:

```json
{
  "_meta": {
    "name": "hybrid-research",
    "requires": ["uvx", "npx", "node>=18"],
    "optional": ["serena"],
    "env_vars": ["GEMINI_API_KEY", "OPENAI_API_KEY"]
  },
  "mcpServers": { ... }
}
```

The switcher can then validate requirements before applying.

---

## Low Priority Issues

### ISSUE-009: Claude Code Installation Hangs After Permission Failure

**Severity**: Low
**Impact**: First-run setup gets stuck, requires container rebuild
**Reproducible**: Yes (intermittent, after certain failure states)
**Status**: ✅ Mitigated (2025-12-16) - Root cause addressed in scbio-docker v0.5.3

#### Description

When `setup-ai.sh` runs after a previous partial failure (e.g., permission denied during toml install), the Claude Code installation process can hang indefinitely at:

```
[INFO] Installing Claude Code native build latest...
```

#### Reproduction Steps

1. Run `setup-ai.sh --force` in a fresh container
2. If toml installation fails with permission error, script continues
3. Claude Code installation begins but never completes
4. Process hangs at "Installing Claude Code native build latest..."
5. Subsequent manual runs of `install_claude.sh` report directory not writable

#### Observed Behavior

```
[INFO] Installing python package: toml...
ERROR: Could not install packages due to an OSError: [Errno 13] Permission denied
[WARN] Standard install failed (read-only venv?). Retrying with --user...
ERROR: Can not perform a '--user' install. User site-packages are not visible in this virtualenv.
[WARN] Failed to install 'toml'. Codex config updates may fail.
...
[INFO] Installing Claude Code native build latest...
# HANGS HERE
```

#### Workaround

Rebuild container from scratch: `Dev Containers: Rebuild Container`

#### Investigation Notes

- May be related to corrupted state in `~/.local/` after partial install
- May be related to npm/node state from previous nvm installation
- Fresh container build succeeds on first attempt

#### Mitigation

**scbio-docker v0.5.3** pre-installs the dependencies that were causing permission failures:
- Node.js 20 LTS (npm, npx) - eliminates nvm permission issues
- uv/uvx - eliminates permission issues with MCP server installation
- `toml` Python package - eliminates permission errors in base venv

With these dependencies pre-installed, the cascade of permission failures that led to installation hangs no longer occurs. Users on older image versions should rebuild with v0.5.3.

---

### ISSUE-010: toml Package Installation Fails in Read-Only Base Venv

**Severity**: Low
**Impact**: Codex config.toml updates may fail
**Reproducible**: Yes
**Status**: ✅ Resolved (2025-12-16)

#### Description

The setup script attempts to install the `toml` Python package for Codex configuration, but fails because:

1. `/opt/venvs/base` is a read-only system venv (by design in scbio-docker)
2. `--user` fallback fails because user site-packages aren't visible in the virtualenv

#### Error Output

```
[INFO] Installing python package: toml...
ERROR: Could not install packages due to an OSError: [Errno 13] Permission denied: '/opt/venvs/base/lib/python3.10/site-packages/toml'
Check the permissions.

[WARN] Standard install failed (read-only venv?). Retrying with --user...
ERROR: Can not perform a '--user' install. User site-packages are not visible in this virtualenv.
[WARN] Failed to install 'toml'. Codex config updates may fail.
```

#### Root Cause

The scbio-docker base image has a two-tier library architecture:
- System venv `/opt/venvs/base` is read-only (owned by root)
- User packages should go to user library, but virtualenv isolation prevents `--user` installs

#### Potential Solutions (Not Implemented)

1. Pre-install `toml` in the base Docker image during build
2. Create a separate user-writable venv for setup scripts
3. Use a different TOML library that's already installed (e.g., Python 3.11+ has `tomllib`)
4. Rewrite Codex config generation to not require Python toml parsing

#### Resolution

**Fix implemented**: Added `toml` package to `docker/requirements/base.txt` in scbio-docker.

The `toml` package is now pre-installed in the base venv during Docker image build (Stage 1: Builder), then copied to the runtime image (Stage 2). This ensures the package is available when `setup-ai.sh` runs for Codex configuration.

**File changed**: `docker/requirements/base.txt`
```
# Configuration parsing (for SciAgent-toolkit Codex config)
toml
```

**Note**: Requires rebuilding the Docker image to take effect:
```bash
scripts/build.sh  # or docker build -f docker/base/Dockerfile
```

---

### ISSUE-011: Agent Activation Script Parses Markdown Headers as Agent Names

**Severity**: Low
**Impact**: Cosmetic warnings during role activation, no functional impact
**Reproducible**: Yes
**Status**: ✅ Resolved (2025-12-16)

#### Description

When `activate-role.sh` activates the base role, it produces numerous warnings about non-existent agents. These appear to be parsing errors where markdown header content and description text are being interpreted as agent names.

#### Example Output

```
[OK]   Agent: bioinf-librarian
[WARN]   Agent not found: # (expected at .../agents/#.md)
[WARN]   Agent not found: Tool/documentation (expected at .../agents/Tool/documentation.md)
[WARN]   Agent not found: research (expected at .../agents/research.md)
[OK]   Agent: bio-research-visualizer
[WARN]   Agent not found: # (expected at .../agents/#.md)
[WARN]   Agent not found: Biological (expected at .../agents/Biological.md)
[WARN]   Agent not found: mechanism (expected at .../agents/mechanism.md)
```

#### Pattern Analysis

The warnings follow a pattern suggesting the script is parsing the agent file's YAML frontmatter `description:` field:
- `#` - Markdown header marker
- `Tool/documentation` - Part of "Tool/documentation research" description
- `Biological mechanism research` - Split into individual words

#### Root Cause (Suspected)

The `roles/base.yaml` file or the parsing logic in `activate-role.sh` may be:
1. Reading agent file headers in addition to names
2. Splitting description text on spaces/special characters
3. Treating each word as a separate agent reference

#### Functional Impact

**None** - The actual agents are correctly symlinked despite the warnings. Only cosmetic output is affected.

#### Resolution

**Fixed** in `scripts/activate-role.sh`: Updated `get_yaml_list()` function to properly strip inline YAML comments.

The parsing pipeline now:
1. Extracts the YAML section for the key
2. Finds list items (lines with `- item`)
3. Removes leading dash and whitespace
4. **NEW**: Strips inline comments (everything after `#`)
5. **NEW**: Removes trailing whitespace
6. Filters out empty lines

**Before (broken):**
```yaml
- bioinf-librarian            # Tool/documentation research
```
Was parsed as multiple items: `bioinf-librarian`, `#`, `Tool/documentation`, `research`

**After (fixed):**
```yaml
- bioinf-librarian            # Tool/documentation research
```
Is correctly parsed as single item: `bioinf-librarian`

---

### ISSUE-012: PAL MCP Not Connected with Default 'coding' Profile

**Severity**: Low
**Impact**: PAL unavailable until profile switch
**Reproducible**: Yes
**Status**: ✅ Resolved (2025-12-16)

#### Description

After running `setup-ai.sh`, PAL MCP is not automatically connected in Claude Code when using the default `coding` profile. It only becomes available after switching to `hybrid` profile via `switch-mcp-profile.sh`.

#### Observed Behavior

1. Run `setup-ai.sh --force` → completes successfully
2. Start `claude` → `/mcp` shows only context7 and sequential-thinking (no PAL)
3. Run `./switch-mcp-profile.sh hybrid` → PAL now connected
4. Restart `claude` → `/mcp` shows context7, sequential-thinking, AND pal

#### Investigation Notes

- Setup completion shows `coding` profile applied
- The `coding` profile template should include PAL
- May be related to ISSUE-001 (PAL definition inconsistency)
- May require explicit `switch-mcp-profile.sh coding` after initial setup

#### Workaround

Run `./switch-mcp-profile.sh hybrid` (or `coding` explicitly) after setup completes.

#### Resolution

**Root Cause**: The `coding.mcp.json` profile template was using the PAL wrapper approach (`${TOOLKIT_ROOT}/mcp_servers/pal/start-pal.py`) which requires `setup_pal.sh` to create the wrapper script. This script was not being called consistently.

**Fix**: Updated `templates/mcp-profiles/coding.mcp.json` to use the direct `uvx` approach (same as hybrid-research profile):

```json
{
  "pal": {
    "type": "stdio",
    "command": "uvx",
    "args": ["--from", "git+https://github.com/BeehiveInnovations/pal-mcp-server.git", "pal-mcp-server"],
    "env": {
      "GEMINI_API_KEY": "${GEMINI_API_KEY}",
      "OPENAI_API_KEY": "${OPENAI_API_KEY}"
    }
  }
}
```

This change also addresses **ISSUE-001** (PAL Server Definition Inconsistency) by standardizing all profiles on the direct uvx approach. Combined with the scbio-docker v0.5.3 pre-installation of `uvx`, PAL now works out-of-the-box with all profiles.

**Note**: Users still need to set `GEMINI_API_KEY` and/or `OPENAI_API_KEY` in their `.env` file for PAL to authenticate with the model providers.

---

### ISSUE-013: ToolUniverse MCP Fails to Connect in Claude Code

**Severity**: Low
**Impact**: ToolUniverse scientific tools unavailable in Claude Code
**Reproducible**: Yes
**Status**: 🔍 Needs Investigation - works in Gemini, suggests config format difference

#### Description

ToolUniverse MCP server consistently fails to connect in Claude Code regardless of which MCP profile is selected (coding, hybrid, research-lite, full, etc.). Other MCP servers (context7, pal, sequential-thinking, serena) connect successfully.

#### Observed Behavior

```
> /mcp
  ⎿  Failed to reconnect to tooluniverse.

> /mcp
 5 servers

 ❯ 1. context7             ✔ connected
   2. pal                  ✔ connected
   3. sequential-thinking  ✔ connected
   4. serena               ✔ connected
   5. tooluniverse         ✘ failed
```

#### Investigation Notes

- Same ToolUniverse installation works in Gemini CLI (connects successfully)
- Suggests configuration format difference between Claude Code and Gemini
- May be path resolution issue in `.mcp.json`
- May be startup timeout issue (ToolUniverse can be slow to initialize)
- Check logs at: `~/.cache/claude-cli-nodejs/-workspaces-<project>/`

#### Debugging Steps

```bash
# Check ToolUniverse can start manually
uv --directory /path/to/tooluniverse-env run tooluniverse-mcp --help

# Check Claude Code logs
cat ~/.cache/claude-cli-nodejs/-workspaces-*/mcp-*.log

# Run Claude Code with debug output
claude --debug
```

---

### ISSUE-014: ToolUniverse in Gemini Shows Limited Tool Set

**Severity**: Low
**Impact**: Missing research tools (PubMed, etc.) in Gemini despite successful connection
**Reproducible**: Yes
**Status**: 🔍 Needs Investigation - may be ToolUniverse version or profile mismatch

#### Description

When ToolUniverse connects successfully in Gemini CLI, it only shows 5 tools instead of the expected research tool set:

```
🟢 tooluniverse - Ready (5 tools)
  Tools:
  - ChEMBL_search_similar_molecules
  - EuropePMC_search_articles
  - find_tools
  - openalex_literature_search
  - search_clinical_trials
```

Missing expected tools:
- PubMed search
- UniProt protein search
- DrugBank search
- FDA adverse events
- Clinical trial details
- Semantic Scholar search

#### Investigation Notes

- May be `--include-tools` filter in profile template being too restrictive
- May be ToolUniverse version difference in available tools
- Check profile template being used (research-lite vs research-full)
- `find_tools` is available - can be used to discover what tools ToolUniverse actually provides

#### Related

- ISSUE-005: ToolUniverse Tool Count Inconsistency (profile misconfiguration)

---

### ISSUE-015: Codex CLI Does Not Detect MCP Servers

**Severity**: Low
**Impact**: Codex CLI MCP integration non-functional
**Reproducible**: Yes
**Status**: 🔍 Needs Investigation - may be Codex CLI version or config format issue

#### Description

Despite `setup-ai.sh` creating `~/.codex/config.toml` with MCP server configurations, Codex CLI does not detect or connect to any MCP servers.

#### Observed Behavior

- `~/.codex/config.toml` exists and contains valid ToolUniverse configuration
- Codex CLI starts successfully
- No MCP servers appear in Codex session
- Gemini CLI (same setup process) successfully connects to MCPs

#### config.toml Content (Created by Setup)

```toml
# ToolUniverse MCP Server Configuration
[mcp_servers.tooluniverse]
command = "/workspaces/.../tooluniverse-env/bin/tooluniverse-mcp"
args = []
startup_timeout_sec = 60
```

#### Investigation Notes

- May be a Codex CLI version compatibility issue
- May require different config format than documented
- May require explicit MCP enable command within Codex
- Gemini uses different config format (JSON) and works correctly

#### Gemini Comparison (Working)

```json
// .gemini/settings.json - this works
{
  "mcpServers": {
    "tooluniverse": {
      "command": "...",
      "args": [...]
    }
  }
}
```

---

### ISSUE-016: Project CLAUDE.md Points to Wrong AGENTS.md (Template Version Mismatch)

**Severity**: Medium
**Impact**: AI agents receive toolkit codebase instructions instead of analysis methodology
**Reproducible**: Yes (projects created before commit 36fae5c)
**Status**: Partially Resolved (template fixed, but existing projects affected)

#### Description

Projects created using SciAgent-toolkit submodule version `16ab338` or earlier have incorrect references in their CLAUDE.md and AGENTS.md files:

**Incorrect (old template):**
```markdown
> **This is a thin wrapper.** Full methodology is in `01_modules/SciAgent-toolkit/AGENTS.md`.
```

**Correct (new template, commit 36fae5c):**
```markdown
> **Full methodology:** See `AGENTS.md` in this directory for comprehensive analysis guidelines.
```

#### The Problem

The old template directed AI agents to `01_modules/SciAgent-toolkit/AGENTS.md`, which is the **toolkit's codebase documentation** (how to develop SciAgent-toolkit), NOT the analysis methodology.

The correct reference should be to the **project's own `AGENTS.md`** (created from `AGENTS.md.template`), which contains the RNAseq analysis methodology, critical rules, and workflow guidelines.

#### File Purpose Clarification

| File | Location | Purpose |
|------|----------|---------|
| `AGENTS.md` | Project root | **Analysis methodology** for AI agents |
| `CLAUDE.md` | Project root | Claude Code context + references AGENTS.md |
| `01_modules/SciAgent-toolkit/AGENTS.md` | Submodule | **Toolkit development docs** (NOT for analysis) |
| `01_modules/SciAgent-toolkit/CLAUDE.md` | Submodule | Toolkit codebase overview |

#### Root Cause

The template fix (commit `36fae5c`) was made AFTER some projects were created using the old submodule version. The submodule in scbio-docker or init-project.sh may point to an older commit.

#### Affected Projects

Projects created with SciAgent-toolkit before commit `36fae5c` (Dec 16, 2025 17:01).

#### Resolution for Existing Projects

**Manual fix** - Update project's CLAUDE.md line 3:

```markdown
# Change FROM:
> **This is a thin wrapper.** Full methodology is in `01_modules/SciAgent-toolkit/AGENTS.md`.

# Change TO:
> **Full methodology:** See `AGENTS.md` in this directory for comprehensive analysis guidelines.
```

Or regenerate context files:
```bash
# Update submodule first
cd 01_modules/SciAgent-toolkit && git pull origin main

# Regenerate (backup existing files first!)
./01_modules/SciAgent-toolkit/scripts/setup-ai.sh --force
```

#### Prevention

Ensure scbio-docker's SciAgent-toolkit submodule is updated to at least commit `36fae5c`:

```bash
# In scbio-docker repo
cd toolkits/SciAgent-toolkit
git fetch origin
git checkout 36fae5c  # or later
cd ../..
git add toolkits/SciAgent-toolkit
git commit -m "Update SciAgent-toolkit to fix template references"
```

---

## Change Log

| Date | Change | Author |
|------|--------|--------|
| 2025-12-16 | Initial document created | Claude Code |
| 2025-12-16 | Fixed Python IndentationError in switch-mcp-profile.sh (ISSUE-007 example) | Claude Code |
| 2025-12-16 | Added ISSUE-009 through ISSUE-011 (Claude install hang, toml permission, agent parsing) | Claude Code |
| 2025-12-16 | Added ISSUE-012 (PAL not connected with coding profile) | Claude Code |
| 2025-12-16 | Added ISSUE-013 through ISSUE-014 (ToolUniverse connection failures) | Claude Code |
| 2025-12-16 | Added ISSUE-015 (Codex MCP detection) | Claude Code |
| 2025-12-16 | Added ISSUE-016 (CLAUDE.md/AGENTS.md template version mismatch) | Claude Code |
| 2025-12-16 | **Resolved ISSUE-010**: Added toml to docker/requirements/base.txt | Claude Code |
| 2025-12-16 | **Mitigated ISSUE-009**: Added Node.js 20 + uv to Dockerfile (scbio-docker v0.5.3) | Claude Code |
| 2025-12-16 | **Resolved ISSUE-011**: Fixed YAML inline comment parsing in activate-role.sh | Claude Code |
| 2025-12-16 | **Resolved ISSUE-001, ISSUE-012**: Standardized PAL on direct uvx in coding.mcp.json | Claude Code |
| 2025-12-16 | **Resolved ISSUE-005**: Verified research-full.mcp.json has --include-tools | Claude Code |
| 2025-12-16 | **Triaged ISSUE-013, 014, 015**: Marked as needs investigation | Claude Code |

---

## Related Documentation

- [ARCHITECTURE.md](ARCHITECTURE.md) - System architecture overview
- [CONFIGURATION.md](CONFIGURATION.md) - Configuration details
- [TROUBLESHOOTING.md](TROUBLESHOOTING.md) - Common issues and solutions
