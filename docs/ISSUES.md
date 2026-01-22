# SciAgent-toolkit Architecture Issues

**Created**: 2025-12-16
**Version**: Current (as of v2.0)
**Status**: Active issue tracking for MCP configuration system

This document catalogs architectural issues, design inconsistencies, and technical debt in the SciAgent-toolkit MCP configuration system. Issues are prioritized by severity and impact.

**Related docs**: [TROUBLESHOOTING.md](TROUBLESHOOTING.md) (user-facing quick fixes) | [ARCHITECTURE_REVIEW.md](ARCHITECTURE_REVIEW.md) (strategic refactoring plan)

---

## Critical Issues

### ISSUE-001: PAL Server Definition Inconsistency Across Profiles

**Severity**: Critical | **Status**: ✅ Resolved (2025-12-16) | **See also**: ISSUE-012

**Resolution**: All profiles now standardized on direct `uvx` approach. Wrapper script (`start-pal.py`) removed from templates. PAL env block uses `${GEMINI_API_KEY}` / `${OPENAI_API_KEY}` placeholders substituted at profile switch time.

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

**Severity**: Medium | **Status**: ✅ Resolved

**Resolution**: `research-full.mcp.json` now includes `--include-tools` with 14 specific tools, preventing 600+ tool context overflow. All profiles now have explicit tool specifications.

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

**Severity**: Low | **Status**: ✅ Mitigated in scbio-docker v0.5.3

**Resolution**: Pre-installing Node.js 20 LTS, uv/uvx, and `toml` in the Docker image eliminates the cascade of permission failures that caused installation hangs. Users on older images should rebuild.

---

### ISSUE-010: toml Package Installation Fails in Read-Only Base Venv

**Severity**: Low | **Status**: ✅ Resolved (2025-12-16)

**Resolution**: Added `toml` to `docker/requirements/base.txt` in scbio-docker. Pre-installed during image build.

---

### ISSUE-011: Agent Activation Script Parses Markdown Headers as Agent Names

**Severity**: Low | **Status**: ✅ Resolved (2025-12-16)

**Resolution**: Fixed `get_yaml_list()` in `activate-role.sh` to strip inline YAML comments (everything after `#`), preventing description words from being interpreted as agent names.

---

### ISSUE-012: PAL MCP Not Connected with Default 'coding' Profile

**Severity**: Low | **Status**: ✅ Resolved (2025-12-16) | **See also**: ISSUE-001

**Resolution**: All profile templates (including `coding.mcp.json`) now use direct `uvx` for PAL instead of wrapper script. Combined with scbio-docker v0.5.3 pre-installing `uvx`, PAL works out-of-the-box. Users must set API keys in `.env`.

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

### ISSUE-017: nvm Conflicts with Pre-installed Node.js + Custom npm Prefix

**Severity**: High
**Impact**: setup-ai.sh fails when Node.js is pre-installed with custom npm prefix
**Reproducible**: Yes
**Status**: 🔍 Design Decision Required - Option C recommended

#### Description

When running `setup-ai.sh` in an environment where Node.js was pre-installed (e.g., from NodeSource) AND has a custom npm prefix configured (`npm config set prefix ~/.npm-global`), the setup fails because nvm is incompatible with existing npm prefix configurations.

#### Error Output

```
[INFO] Installing nvm...
=> You currently have modules installed globally with `npm`. These will no
=> longer be linked to the active version of Node when you install a new node
=> with `nvm`; and they may (depending on how you construct your `$PATH`)
=> override the binaries of modules installed with `nvm`:

/home/devuser/.npm-global/lib
├── @anthropic-ai/claude-code@2.0.71
├── ts-node@10.9.2
└── typescript@5.9.3

[INFO] Installing/using Node.js LTS via nvm...
Installing latest LTS version.
Downloading and installing node v24.12.0...
...
Your user's .npmrc file (${HOME}/.npmrc)
has a `globalconfig` and/or a `prefix` setting, which are incompatible with nvm.
Run `nvm use --delete-prefix v24.12.0` to unset it.
[ERROR] Setup failed. Check the output above for errors.
```

#### Root Cause

**Two conflicting Node.js management strategies:**

| Component | Strategy | Mechanism |
|-----------|----------|-----------|
| **scbio-docker v0.5.3** | System-level Node.js | NodeSource PPA (`/usr/bin/node`) |
| **SciAgent-toolkit** | User-level Node.js | nvm (`~/.nvm/versions/node/...`) |
| **Custom prefix pattern** | User-writable global packages | `npm config set prefix ~/.npm-global` |

The conflict occurs when:
1. Node.js is pre-installed at system level (scbio-docker Dockerfile lines 229-235)
2. User configures `npm config set prefix ~/.npm-global` to avoid sudo for global packages
3. This creates `~/.npmrc` with a `prefix` setting
4. `setup-ai.sh` calls `ensure_nvm()` which installs nvm
5. nvm installs fine, but `nvm use` fails because nvm doesn't work with custom npm prefixes

#### Affected Code Paths

**scbio-docker Dockerfile** (docker/base/Dockerfile:229-235):
```dockerfile
# Install Node.js (for MCP servers: npx required)
RUN curl -fsSL https://deb.nodesource.com/setup_20.x | bash - && \
    apt-get install -y nodejs && \
    rm -rf /var/lib/apt/lists/* && \
    node --version && npm --version && npx --version
```

**SciAgent-toolkit common.sh** (scripts/common.sh:144-187):
```bash
ensure_nvm() {
    export NVM_DIR="$HOME/.nvm"
    # ... installs nvm if not present ...
    if ! command -v node &>/dev/null || [[ $(which node) == "/usr/bin/node" ]]; then
        log_info "Installing/using Node.js LTS via nvm..."
        nvm install --lts  # <-- This fails when .npmrc has prefix
        nvm use --lts
    fi
}
```

**The trigger condition** (common.sh:172):
```bash
[[ $(which node) == "/usr/bin/node" ]]
```
This condition explicitly detects system Node.js and tries to switch to nvm-managed Node.js, which triggers the conflict.

#### Why This Design Exists

1. **scbio-docker pre-installs Node.js**: Ensures MCP servers work immediately without installation delays
2. **SciAgent-toolkit prefers nvm**: Allows user-controlled Node.js versions without sudo
3. **Custom npm prefix pattern**: Common workaround for containers where `/usr/lib/node_modules` is read-only

The three patterns are individually valid but mutually incompatible.

#### Recommended Solution: Option C - Clean Separation of Concerns

**Don't pre-install Node.js in Docker image. Let SciAgent-toolkit manage it via nvm.**

This is the cleanest architectural approach:

| Responsibility | Owner |
|----------------|-------|
| Base OS + build tools + R + Python venvs | scbio-docker Dockerfile |
| Node.js + npm + AI CLI tools | SciAgent-toolkit (via nvm) |

**Changes required in scbio-docker:**

1. **Remove Node.js installation from Dockerfile** (lines 229-235):
   ```dockerfile
   # REMOVE THIS SECTION:
   # Install Node.js (for MCP servers: npx required)
   # RUN curl -fsSL https://deb.nodesource.com/setup_20.x | bash - && \
   #     apt-get install -y nodejs && \
   #     rm -rf /var/lib/apt/lists/* && \
   #     node --version && npm --version && npx --version
   ```

2. **Keep only curl** (already present for other downloads)

3. **Update CLAUDE.md** to reflect that Node.js is installed at runtime by SciAgent-toolkit

**Benefits:**
- Clean separation: scbio-docker = science stack, SciAgent-toolkit = AI stack
- No conflicting Node.js installations
- User gets nvm for version management
- Consistent behavior across all environments (not just scbio-docker)

**Tradeoffs:**
- First-run setup takes ~1-2 minutes longer (nvm + Node.js download)
- Requires network access during setup

#### Alternative Options (Not Recommended)

**Option A: Skip nvm when Node.js is pre-installed**
- Would require modifying `ensure_nvm()` to detect and accept system Node.js
- Violates SciAgent-toolkit's design of user-controlled Node.js
- Creates implicit dependency on host environment

**Option B: Remove conflicting .npmrc before nvm installation**
- Destructive: removes user's npm configuration
- May break other tools that depend on that configuration

**Option D: Document and provide manual workaround**
- Poor UX: user must troubleshoot and fix manually
- Doesn't solve root cause

#### User Workarounds (Until Fix is Implemented)

Users encountering this error can:

1. **Remove the conflicting .npmrc**:
   ```bash
   rm ~/.npmrc
   ./modules/SciAgent-toolkit/scripts/setup-ai.sh
   ```

2. **Remove only the prefix line**:
   ```bash
   sed -i '/^prefix/d' ~/.npmrc
   ./modules/SciAgent-toolkit/scripts/setup-ai.sh
   ```

3. **Rebuild container without pre-installed Node.js**:
   Edit Dockerfile to remove Node.js installation, rebuild image

#### Related Issues

- ISSUE-009: Claude Code Installation Hangs After Permission Failure
- ISSUE-010: toml Package Installation Fails in Read-Only Base Venv

---

### ISSUE-018: PAL clink to Gemini CLI Fails Due to Multiple Configuration Conflicts

**Severity**: High
**Impact**: PAL's `clink` tool cannot delegate to Gemini CLI out-of-the-box
**Reproducible**: Yes
**Status**: 🔍 Documented - Multiple root causes identified

#### Description

When using `mcp__pal__clink(cli_name: "gemini", ...)` from Claude Code, the call fails with one or more errors:

1. **YOLO mode conflict**: Return code 52, "Cannot start in YOLO mode when it is disabled by settings"
2. **Markdown parsing error**: `[ERROR] [ImportProcessor] Could not find child token in parent raw content`
3. **Missing API key**: "When using Gemini API, you must specify the GEMINI_API_KEY environment variable"

#### Root Causes

**Issue A: YOLO Mode Hardcoded in PAL**

PAL's embedded CLI client config (`conf/cli_clients/gemini.json`) contains:
```json
{
  "additional_args": ["--yolo"]
}
```

This conflicts with Gemini's default user settings (`~/.gemini/settings.json`):
```json
{
  "security": {
    "disableYoloMode": true
  }
}
```

**Problem**: PAL assumes YOLO mode is available, but Gemini CLI's default security settings disable it.

**Issue B: Gemini CLI Markdown Parser Sensitivity**

Gemini CLI's `ImportProcessor` fails to parse certain markdown constructs in project context files (`GEMINI.md`, `AGENTS.md`):

```markdown
# This causes parsing failure:
> **Full methodology:** See `AGENTS.md` in this directory...
```

The blockquote (`>`) with nested bold (`**`) triggers a token parsing error.

**Issue C: Environment Variable Inheritance Gap**

PAL's `clink` tool spawns external CLI processes that inherit environment variables from the shell. Unlike PAL's internal API calls (which can read `.env` files), external CLIs require:

| Requirement | PAL Internal API | PAL clink (External CLI) |
|-------------|------------------|--------------------------|
| API key source | Reads `.env` files | Inherits from shell env |
| User action needed | None (automatic) | Must export vars manually |

If `GEMINI_API_KEY` isn't in the shell environment when Claude Code starts, the spawned Gemini CLI cannot authenticate.

#### Affected Components

| Component | Issue | Owner |
|-----------|-------|-------|
| PAL MCP (`conf/cli_clients/gemini.json`) | Hardcoded `--yolo` | PAL upstream |
| Gemini CLI | `disableYoloMode: true` default | Google |
| Gemini CLI | Markdown ImportProcessor | Google |
| SciAgent-toolkit templates | Problematic markdown patterns | SciAgent-toolkit |
| Shell environment | `.env` not auto-loaded | User/devcontainer config |

#### Current Workarounds

See [TROUBLESHOOTING.md](TROUBLESHOOTING.md#pal-clink-to-gemini-cli-issues) for user-facing solutions:

1. **YOLO mode**: Edit `~/.gemini/settings.json` to set `disableYoloMode: false`
2. **Markdown parsing**: Simplify `GEMINI.md` to avoid blockquotes with nested formatting
3. **API key**: Add `source .env` to `~/.bashrc` or export before starting Claude Code

#### Recommended Fixes

**Short-term (SciAgent-toolkit):**

1. Update `templates/vendor/GEMINI.md.template` line 3 to avoid problematic markdown patterns:
   ```markdown
   # Change FROM:
   > **Full methodology:** See `AGENTS.md` in this directory...

   # Change TO:
   Full methodology is documented in `AGENTS.md`. This file contains Gemini-specific context.
   ```
2. Document the YOLO mode requirement in setup scripts
3. Add `.env` sourcing to devcontainer `postCreateCommand`

**Medium-term (PAL upstream):**

1. Make `--yolo` flag configurable or check if YOLO is enabled before passing
2. Pass environment variables from `.env` files to spawned CLI processes
3. Add CLI compatibility checks before spawning

**Long-term (Gemini CLI upstream):**

1. Improve ImportProcessor robustness for common markdown patterns
2. Consider making YOLO mode opt-out rather than opt-in for programmatic usage

#### Related Issues

- ISSUE-003: Two-Stage Environment Binding Creates Stale Configuration (same `.env` inheritance problem)
- ISSUE-016: Project CLAUDE.md Points to Wrong AGENTS.md (template quality)

---

### ISSUE-019: PAL MCP Fails Due to PYTHONPATH Collision with Project config.py

**Severity**: Critical | **Status**: ✅ Resolved (2026-01-18)

**Problem**: PAL's `from config import ...` resolves to project's `config.py` when PYTHONPATH includes project scripts directory, causing `ModuleNotFoundError`.

**Resolution**: Added `"PYTHONPATH": ""` to PAL's env block in all profile templates (7 files: `templates/mcp-profiles/*.mcp.json` + `templates/gemini-profiles/*.json`). Re-run profile switcher to apply.

---

### ISSUE-020: API Keys Not Auto-Sourced from .env in Shell Environment

**Severity**: High
**Impact**: MCP profile switching fails to substitute API keys, resulting in `${GEMINI_API_KEY}` literal strings in .mcp.json
**Reproducible**: Yes
**Status**: ✅ Documented - User must add .env sourcing to bashrc

#### Description

The `.devcontainer/.env` file containing API keys is not automatically sourced into the shell environment. The `switch-mcp-profile.sh` script substitutes `${GEMINI_API_KEY}` etc. from the shell environment, but if those variables aren't set, the literal placeholder strings end up in `.mcp.json`.

#### Root Cause

1. Docker Compose reads `.devcontainer/.env` during container creation
2. This sets environment variables for the container lifecycle, NOT for interactive shells
3. When user opens a terminal, shell environment is separate from `.env` file contents
4. `switch-mcp-profile.sh` calls `load_env()` which sources `.env`, but only for that script's session
5. The script then uses Python to substitute, but Python inherits the original shell's environment

#### Solution

Add auto-sourcing to `~/.bashrc`:

```bash
# Auto-source project .env files for API keys
_source_project_env() {
    local project_dir="${PWD}"
    if [ -f "${project_dir}/.devcontainer/.env" ]; then
        set -a
        source "${project_dir}/.devcontainer/.env"
        set +a
    fi
    if [ -f "${project_dir}/.env" ]; then
        set -a
        source "${project_dir}/.env"
        set +a
    fi
}
_source_project_env
```

#### Related Issues

- ISSUE-003: Two-Stage Environment Binding Creates Stale Configuration
- ISSUE-018: PAL clink to Gemini CLI Fails Due to Multiple Configuration Conflicts

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
| 2025-12-16 | Added ISSUE-017 (nvm conflicts with pre-installed Node.js + custom npm prefix) | Claude Code |
| 2025-12-29 | Added ISSUE-018 (PAL clink to Gemini CLI configuration conflicts) | Claude Code |
| 2026-01-18 | Added ISSUE-019 (PYTHONPATH collision with project config.py) - **Resolved** | Claude Code |
| 2026-01-18 | Added ISSUE-020 (API keys not auto-sourced from .env) - **Documented** | Claude Code |
| 2026-01-18 | Fixed profile templates: Added PYTHONPATH clearing to PAL env block | Claude Code |

---

## Related Documentation

- [ARCHITECTURE.md](ARCHITECTURE.md) - System architecture overview
- [CONFIGURATION.md](CONFIGURATION.md) - Configuration details
- [TROUBLESHOOTING.md](TROUBLESHOOTING.md) - Common issues and solutions
