# Security: API Key Management

> **Last Updated:** 2024-12-19
> **Version:** 1.0.0

This document describes the security architecture for managing API keys in SciAgent-toolkit and explains the design decisions behind it.

---

## Table of Contents

1. [Overview](#overview)
2. [The Problem](#the-problem)
3. [Architecture](#architecture)
4. [Per-CLI Behavior](#per-cli-behavior)
5. [Security Validation](#security-validation)
6. [Best Practices](#best-practices)
7. [Troubleshooting](#troubleshooting)

---

## Overview

SciAgent-toolkit integrates multiple AI CLI tools (Claude Code, Codex CLI, Gemini CLI) that each require API keys for MCP servers like PAL. This document explains:

- Why API keys appear in config files
- How to prevent accidental exposure
- The security model and its trade-offs

### Key Principle

```
+--------------------------------------------------+
|           SINGLE SOURCE OF TRUTH                 |
|                                                  |
|    .devcontainer/.env  OR  ~/.env               |
|    (Always gitignored)                          |
+--------------------------------------------------+
                      |
      +---------------+---------------+
      |               |               |
      v               v               v
+----------+   +-----------+   +------------+
| Claude   |   | Codex CLI |   | Gemini CLI |
| .mcp.json|   | config.toml|  | settings.json
| (git-ign)|   | (env_vars) |   | (git-ign)  |
+----------+   +-----------+   +------------+
```

**All files containing actual API keys MUST be gitignored.**

---

## The Problem

### Why Do API Keys Appear in Config Files?

During investigation (December 2024), we found that API keys were being written to:
- `.mcp.json` (Claude Code)
- `.gemini/settings.json` (Gemini CLI)

**Root Cause:** The `switch-mcp-profile.sh` script intentionally substitutes `${GEMINI_API_KEY}` placeholders with actual values from environment variables.

```python
# From switch-mcp-profile.sh lines 180-183
# API key substitutions (ISSUE-002 fix)
content = content.replace('${GEMINI_API_KEY}', os.environ.get('GEMINI_API_KEY', ''))
content = content.replace('${OPENAI_API_KEY}', os.environ.get('OPENAI_API_KEY', ''))
```

This was implemented as a workaround for [Claude Code Issue #1254](https://github.com/anthropics/claude-code/issues/1254) where environment variable expansion in the `env` block of MCP configs was unreliable.

### Why Can't We Just Use ${VAR} References?

Each CLI tool handles environment variables differently:

| CLI Tool | `${VAR}` Expansion | Secure Alternative |
|----------|-------------------|-------------------|
| Claude Code | Partial (bugs in `env` block) | Substitution required |
| Codex CLI | Not supported | `env_vars` whitelist |
| Gemini CLI | Not supported | Hardcoded (by design) |

---

## Architecture

### File Locations and Gitignore Status

| File | Contains Keys? | Must be Gitignored? | Notes |
|------|---------------|---------------------|-------|
| `.devcontainer/.env` | Yes (source) | Yes | Single source of truth |
| `.env` | Yes (if used) | Yes | Alternative source |
| `.mcp.json` | Yes (substituted) | **YES** | Claude Code config |
| `.gemini/settings.json` | Yes (substituted) | Yes | Gemini CLI config |
| `.claude/settings.local.json` | No (usually) | Yes (precaution) | Claude settings |
| `~/.codex/config.toml` | No (uses env_vars) | N/A (user-global) | Codex config |

### Required .gitignore Entries

Your project's `.gitignore` MUST contain:

```gitignore
# Environment files (contain API keys)
.env
.devcontainer/.env
.devcontainer/.env.local

# AI CLI configurations (contain API keys)
.gemini/
.mcp.json
.claude/settings.local.json
```

---

## Per-CLI Behavior

### Claude Code

**Config File:** `.mcp.json` (project-local)

**How Keys Are Handled:**
1. Profile templates contain `${GEMINI_API_KEY}` placeholders
2. `switch-mcp-profile.sh` substitutes with actual values from environment
3. Resulting `.mcp.json` contains actual keys
4. File MUST be gitignored

**Why Substitution?**
Claude Code's `${VAR}` expansion has known issues in the `env` block of MCP server configs ([Issue #1254](https://github.com/anthropics/claude-code/issues/1254)).

**Security:** The script now warns if `.mcp.json` is not gitignored.

### Codex CLI

**Config File:** `~/.codex/config.toml` (user-global)

**How Keys Are Handled:**
1. Config uses `env_vars` to whitelist shell environment variables
2. **No actual keys stored in config**
3. Keys are forwarded from shell at runtime

**Example (Secure Pattern):**
```toml
[mcp_servers.pal]
type = "stdio"
command = "uvx"
args = ["--from", "git+https://github.com/BeehiveInnovations/pal-mcp-server.git", "pal-mcp-server"]
env_vars = ["GEMINI_API_KEY", "OPENAI_API_KEY", "OPENROUTER_API_KEY"]
```

**Security:** This is the most secure pattern - keys never touch config files.

### Gemini CLI

**Config File:** `.gemini/settings.json` (project-local)

**How Keys Are Handled:**
1. Gemini CLI stores actual API key values in settings.json
2. This is by design - Gemini CLI does not support `${VAR}` expansion
3. File MUST be gitignored

**Security:** The `.gemini/` directory must always be in `.gitignore`.

### PAL MCP Server

**How Keys Are Read:**
PAL MCP Server reads API keys directly from environment variables at runtime:
- `GEMINI_API_KEY`
- `OPENAI_API_KEY`
- `OPENROUTER_API_KEY`
- `XAI_API_KEY`

The keys in `.mcp.json` are passed to PAL via the `env` block configuration.

---

## Security Validation

### Automatic Validation

Security validation runs automatically when using these scripts:

| Script | When Validation Runs |
|--------|---------------------|
| `setup-ai.sh` | At the end of setup |
| `switch-mcp-profile.sh` | After switching profiles |
| `activate-role.sh` | After role activation |
| `configure_mcp_servers.sh` | After configuration |

You don't need to run validation manually - it's built into the workflow.

### validate-secrets.sh

A validation script is also available for manual checks:

```bash
# Run from project root
./01_modules/SciAgent-toolkit/scripts/validate-secrets.sh

# Or specify project directory
./01_modules/SciAgent-toolkit/scripts/validate-secrets.sh /path/to/project
```

**What It Checks:**
1. All files that may contain API keys are gitignored
2. No API keys in staged git changes
3. Templates use placeholders, not actual keys

### Pre-commit Hook

Install as a git hook for automatic checking:

```bash
# Create symlink
ln -s ../../01_modules/SciAgent-toolkit/scripts/validate-secrets.sh .git/hooks/pre-commit

# Or copy the script
cp 01_modules/SciAgent-toolkit/scripts/validate-secrets.sh .git/hooks/pre-commit
chmod +x .git/hooks/pre-commit
```

### Manual Verification

```bash
# Check if .mcp.json is gitignored
git check-ignore -v .mcp.json

# Check for API key patterns in a file
grep -E 'sk-proj-|AIzaSy|sk-or-v1-' .mcp.json

# Check what would be committed
git diff --cached --name-only
```

---

## Best Practices

### DO

1. **Keep keys in `.devcontainer/.env` only** - This is your single source of truth
2. **Always verify gitignore** before committing
3. **Use `validate-secrets.sh`** before pushing
4. **Rotate keys immediately** if accidentally exposed
5. **Use environment-specific keys** for dev/staging/prod

### DON'T

1. **Never commit API keys** to version control
2. **Never share `.mcp.json` or `.gemini/settings.json`** between machines
3. **Never put keys in templates** - use `${VAR}` placeholders only
4. **Never disable the gitignore warnings** in `switch-mcp-profile.sh`

### Setting Up a New Environment

When setting up a new development environment:

```bash
# 1. Create .env file from template
cp 01_modules/SciAgent-toolkit/templates/.env.template .devcontainer/.env

# 2. Edit and add your API keys
nano .devcontainer/.env

# 3. Load environment
source .devcontainer/.env

# 4. Run profile switcher (will substitute keys)
./01_modules/SciAgent-toolkit/scripts/switch-mcp-profile.sh coding

# 5. Verify security
./01_modules/SciAgent-toolkit/scripts/validate-secrets.sh
```

---

## Troubleshooting

### "API keys but is NOT git-ignored!"

**Problem:** `switch-mcp-profile.sh` or `validate-secrets.sh` warns that a config file is not gitignored.

**Solution:** Add the file to `.gitignore`:
```bash
echo ".mcp.json" >> .gitignore
echo ".gemini/" >> .gitignore
```

### Keys Not Working in PAL

**Problem:** PAL MCP server can't authenticate.

**Check:**
1. Environment variables are set: `echo $GEMINI_API_KEY`
2. Keys were loaded before profile switch: `source .devcontainer/.env`
3. `.mcp.json` contains the keys (should, if substitution worked)

### Codex CLI Not Finding Keys

**Problem:** Codex CLI's MCP servers fail to authenticate.

**Check:**
1. Environment variables are set in your shell
2. Config uses `env_vars` not `env`:
   ```toml
   # Correct
   env_vars = ["GEMINI_API_KEY"]

   # Wrong
   env = { GEMINI_API_KEY = "..." }
   ```

### Accidentally Committed API Keys

**Immediate Steps:**
1. **Rotate all exposed keys immediately** - Generate new keys from provider dashboards
2. Remove from git history (if not pushed):
   ```bash
   git reset HEAD~1
   git checkout -- .mcp.json
   ```
3. If already pushed, contact your security team and rotate keys

---

## References

- [Claude Code Settings Documentation](https://code.claude.com/docs/en/settings)
- [Claude Code MCP Bug #1254](https://github.com/anthropics/claude-code/issues/1254)
- [Codex CLI MCP Documentation](https://developers.openai.com/codex/mcp)
- [PAL MCP Server Configuration](https://github.com/BeehiveInnovations/pal-mcp-server/blob/main/docs/configuration.md)
- [12-Factor App: Config](https://12factor.net/config)

---

## Changelog

### 1.0.0 (2024-12-19)
- Initial security documentation
- Added `validate-secrets.sh` script
- Updated Codex config generation to use `env_vars` pattern
- Added gitignore warnings to `switch-mcp-profile.sh`
