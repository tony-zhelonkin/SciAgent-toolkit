# SciAgent-toolkit Security Model

This document describes API key handling in SciAgent-toolkit, identifies security concerns, and provides recommendations.

---

## 1. Current Security Architecture

### 1.1 Key Flow Overview

```
┌─────────────────────────────────────────────────────────────────────────┐
│                         API KEY LIFECYCLE                                │
├─────────────────────────────────────────────────────────────────────────┤
│                                                                          │
│  STORAGE                                                                 │
│  ────────                                                                │
│  .env.template ──(copy)──> .env ──(user edits)──> Keys added            │
│       │                      │                                           │
│       │                      ├── /.env (project root)                    │
│       │                      └── /.devcontainer/.env (container)         │
│       │                                                                  │
│  INJECTION                                                               │
│  ─────────                                                               │
│  switch-mcp-profile.sh                                                   │
│       │                                                                  │
│       ├── Loads .env via load_env()                                      │
│       ├── Substitutes ${KEY} → actual_value                              │
│       └── Writes to .mcp.json  ←── KEYS NOW IN FILE                      │
│                                                                          │
│  CONSUMPTION                                                             │
│  ───────────                                                             │
│  Claude Code reads .mcp.json                                             │
│       │                                                                  │
│       ├── PAL MCP: loads .env via wrapper (runtime)                      │
│       ├── Context7: reads from .mcp.json env section                     │
│       └── ToolUniverse: reads from .mcp.json env section                 │
│                                                                          │
└─────────────────────────────────────────────────────────────────────────┘
```

### 1.2 Key Storage Locations

| Location | Purpose | Git Status |
|----------|---------|------------|
| `/.env` | Project-level secrets | `.gitignore` (should be) |
| `/.devcontainer/.env` | Container secrets | `.gitignore` (should be) |
| `/.mcp.json` | MCP server config with keys | `.gitignore` (critical!) |
| `~/.gemini/settings.json` | Gemini CLI auth | User home (not in repo) |
| `~/.codex/config.toml` | Codex CLI auth | User home (not in repo) |

### 1.3 Template Substitution System

**Templates with placeholders** (`templates/mcp-profiles/*.mcp.json`):
```json
{
  "pal": {
    "env": {
      "GEMINI_API_KEY": "${GEMINI_API_KEY}",
      "OPENAI_API_KEY": "${OPENAI_API_KEY}"
    }
  }
}
```

**After `switch-mcp-profile.sh`** (`.mcp.json` in project root):
```json
{
  "pal": {
    "env": {
      "GEMINI_API_KEY": "AIzaSy...(actual key)",
      "OPENAI_API_KEY": "sk-proj-...(actual key)"
    }
  }
}
```

---

## 2. User Workflow (Current State)

### 2.1 Initial Setup

After running `setup-ai.sh`, the user must:

| Step | Command/Action | Notes |
|------|----------------|-------|
| 1 | `setup-ai.sh` runs | Creates `.env` from template |
| 2 | Edit `.env` manually | Add API keys by hand |
| 3 | `switch-mcp-profile.sh coding` | Substitutes keys into `.mcp.json` |
| 4 | Source `.env` in shell | `source .devcontainer/.env` |
| 5 | Verify with `validate-secrets.sh` | Optional but recommended |

### 2.2 Each New Shell Session

```bash
# Must run manually each time (not automated)
source /workspaces/DC_Dictionary/.devcontainer/.env
```

### 2.3 After Adding/Changing Keys

```bash
# 1. Edit .env file
vim .devcontainer/.env

# 2. Re-run profile switch to update .mcp.json
./01_scripts/SciAgent-toolkit/scripts/switch-mcp-profile.sh coding

# 3. Restart Claude Code to pick up new .mcp.json
```

---

## 3. Security Issues Identified

### 3.1 Critical Issues

| Issue | Severity | File | Lines | Description |
|-------|----------|------|-------|-------------|
| **Keys written to .mcp.json** | Critical | `switch-mcp-profile.sh` | 157-191 | After substitution, actual API keys exist in `.mcp.json` file on disk |

**Risk:** If `.mcp.json` is accidentally committed or shared, all API keys are exposed.

**Current Mitigation:** `.gitignore` includes `.mcp.json`

**Weakness:** Relies on `.gitignore` being correct and not bypassed (`git add -f`)

### 3.2 Medium Issues

| Issue | Severity | File | Lines | Description |
|-------|----------|------|-------|-------------|
| **No auto-sourcing** | Medium | `devcontainer.json` | 46 | `.env` not sourced automatically in shell |
| **Passive validation** | Medium | `validate-secrets.sh` | 73-177 | Warns about issues but doesn't block commits |

### 3.3 Low Issues

| Issue | Severity | File | Lines | Description |
|-------|----------|------|-------|-------------|
| **Multiple .env locations** | Low | `common.sh` | 228-252 | `/.env` vs `/.devcontainer/.env` causes confusion |
| **No encryption** | Low | `.env` files | N/A | Plaintext storage relies on file permissions |
| **Manual hook setup** | Low | `validate-secrets.sh` | 17 | Pre-commit hook must be manually linked |

---

## 4. Protection Mechanisms

### 4.1 .gitignore Patterns

The following should be in `.gitignore`:
```gitignore
# Secrets
.env
.devcontainer/.env
.mcp.json

# CLI configs
.gemini/
.codex/
```

**Verify with:**
```bash
git check-ignore .env .mcp.json .devcontainer/.env
```

### 4.2 Validation Script

`validate-secrets.sh` performs four checks:

1. **Gitignore verification** - Confirms sensitive files are ignored
2. **Staged file scan** - Checks git staged changes for key patterns
3. **Content scan** - Scans files for API key regex patterns
4. **Template validation** - Ensures templates have `${PLACEHOLDER}` not actual keys

**API Key Patterns Detected:**
```
sk-proj-[A-Za-z0-9_-]{20,}     # OpenAI project
sk-[A-Za-z0-9_-]{40,}          # OpenAI legacy
AIzaSy[A-Za-z0-9_-]{30,}       # Google/Gemini
sk-or-v1-[a-f0-9]{64}          # OpenRouter
xai-[A-Za-z0-9_-]{20,}         # X.AI/Grok
ghp_[A-Za-z0-9]{36}            # GitHub PAT
gho_[A-Za-z0-9]{36}            # GitHub OAuth
```

### 4.3 PAL Runtime Loading

PAL MCP server uses a wrapper that loads `.env` at runtime:
- **File:** `scripts/mcp_servers/setup_pal.sh` (lines 114-180)
- **Advantage:** Keys loaded from `.env`, not read from `.mcp.json`
- **Fallback:** Still reads from `.mcp.json` if `.env` unavailable

---

## 5. Recommendations

### 5.1 Immediate User Actions

#### A. Auto-source .env in Shell Profile

Add to `~/.bashrc` or `~/.zshrc`:
```bash
# Auto-load project environment
if [ -f "/workspaces/DC_Dictionary/.devcontainer/.env" ]; then
    set -a
    source /workspaces/DC_Dictionary/.devcontainer/.env
    set +a
fi
```

#### B. Use Single Canonical Location

Choose ONE location and stick to it:
- **Recommended:** `/.devcontainer/.env` (works in containers)
- Delete or symlink the other

#### C. Enable Pre-commit Hook

```bash
cd /workspaces/DC_Dictionary
ln -sf ../../01_scripts/SciAgent-toolkit/scripts/validate-secrets.sh .git/hooks/pre-commit
chmod +x .git/hooks/pre-commit
```

#### D. Set File Permissions

```bash
chmod 600 .devcontainer/.env
chmod 600 .mcp.json
```

### 5.2 Future Toolkit Improvements

| Improvement | Priority | Description |
|-------------|----------|-------------|
| **Runtime-only injection** | High | Pass keys via environment, never write to `.mcp.json` |
| **Auto-source in devcontainer** | High | Add to `postStartCommand` in `devcontainer.json` |
| **Enforce validation** | Medium | Block commits if keys detected (exit 1) |
| **Pre-commit auto-install** | Medium | Install hook automatically in `setup-ai.sh` |
| **Secrets manager** | Low | Integrate with 1Password, HashiCorp Vault, etc. |

---

## 6. Security Checklist

Run this verification after setup:

```bash
#!/bin/bash
echo "=== SciAgent-toolkit Security Check ==="

# 1. Gitignore
echo -e "\n1. Checking .gitignore..."
for f in .env .mcp.json .devcontainer/.env; do
    git check-ignore -q "$f" 2>/dev/null && echo "   [OK] $f ignored" || echo "   [FAIL] $f NOT ignored!"
done

# 2. File permissions
echo -e "\n2. Checking file permissions..."
for f in .devcontainer/.env .mcp.json; do
    if [ -f "$f" ]; then
        perms=$(stat -c %a "$f" 2>/dev/null || stat -f %Lp "$f")
        [ "$perms" = "600" ] && echo "   [OK] $f (600)" || echo "   [WARN] $f ($perms) - should be 600"
    fi
done

# 3. Pre-commit hook
echo -e "\n3. Checking pre-commit hook..."
[ -x .git/hooks/pre-commit ] && echo "   [OK] Pre-commit hook installed" || echo "   [WARN] No pre-commit hook"

# 4. Environment variables
echo -e "\n4. Checking environment..."
[ -n "$GEMINI_API_KEY" ] && echo "   [OK] GEMINI_API_KEY set" || echo "   [WARN] GEMINI_API_KEY not in environment"
[ -n "$OPENAI_API_KEY" ] && echo "   [OK] OPENAI_API_KEY set" || echo "   [INFO] OPENAI_API_KEY not set (optional)"

echo -e "\n=== Check complete ==="
```

---

## 7. Quick Reference

### Key Files

| File | Location | Purpose |
|------|----------|---------|
| `.env.template` | `01_scripts/SciAgent-toolkit/templates/` | Template with empty placeholders |
| `common.sh` | `01_scripts/SciAgent-toolkit/scripts/` | `load_env()` function (lines 228-252) |
| `switch-mcp-profile.sh` | `01_scripts/SciAgent-toolkit/scripts/` | Substitutes keys (lines 157-191) |
| `validate-secrets.sh` | `01_scripts/SciAgent-toolkit/scripts/` | Validation checks |
| `setup_pal.sh` | `01_scripts/SciAgent-toolkit/scripts/mcp_servers/` | PAL wrapper (lines 114-180) |

### Commands

```bash
# Validate security
./01_scripts/SciAgent-toolkit/scripts/validate-secrets.sh

# Switch MCP profile (re-inject keys)
./01_scripts/SciAgent-toolkit/scripts/switch-mcp-profile.sh coding

# Check if files are git-ignored
git check-ignore .env .mcp.json

# Source environment
source .devcontainer/.env
```

---

## Version

- **Document version:** 1.0.0
- **Created:** 2025-12-29
- **Based on:** SciAgent-toolkit analysis

```