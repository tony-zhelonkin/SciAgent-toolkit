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

## Change Log

| Date | Change | Author |
|------|--------|--------|
| 2025-12-16 | Initial document created | Claude Code |
| 2025-12-16 | Fixed Python IndentationError in switch-mcp-profile.sh (ISSUE-007 example) | Claude Code |

---

## Related Documentation

- [ARCHITECTURE.md](ARCHITECTURE.md) - System architecture overview
- [CONFIGURATION.md](CONFIGURATION.md) - Configuration details
- [TROUBLESHOOTING.md](TROUBLESHOOTING.md) - Common issues and solutions
