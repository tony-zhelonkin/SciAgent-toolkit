# SciAgent-toolkit Architecture Review

**Created**: 2026-01-18
**Status**: Active - Complexity Analysis & Improvement Plan
**Version**: 1.0.0

This document analyzes the current architecture of SciAgent-toolkit, identifies sources of accidental complexity, and proposes a plan to reduce the interoperation surface area while maintaining functionality.

---

## Executive Summary

The SciAgent-toolkit has grown organically, resulting in **accidental complexity** that increases the maintenance burden and bug surface area. This review identifies the core issues and proposes a phased approach to simplify the architecture.

**Key Finding**: The current architecture has **too many scripts with overlapping responsibilities** and **implicit dependencies** that create fragile interoperation chains.

---

## Current Architecture Analysis

### Script Dependency Graph

```
setup-ai.sh (Entry Point)
    │
    ├── common.sh (Utilities)
    │       └── load_env() ←── .devcontainer/.env (NOT auto-sourced to shell!)
    │
    ├── setup_mcp_infrastructure.sh
    │       ├── install_claude.sh
    │       ├── install_codex.sh
    │       ├── install_gemini.sh
    │       └── mcp_servers/
    │               ├── setup_pal.sh
    │               ├── setup_sequential_thinking.sh
    │               ├── setup_tooluniverse.sh
    │               └── setup_serena.sh
    │
    ├── configure_mcp_servers.sh (DEPRECATED - delegates to switch-mcp-profile.sh)
    │       └── switch-mcp-profile.sh (The REAL configuration engine)
    │
    └── activate-role.sh
            └── roles/*.yaml
```

### Identified Complexity Issues

#### ISSUE A: Two-Layer Configuration System (Critical)

**Problem**: `configure_mcp_servers.sh` and `switch-mcp-profile.sh` overlap in responsibility.

| Script | Original Intent | Current Reality |
|--------|-----------------|-----------------|
| `configure_mcp_servers.sh` | Generate MCP configs from detected servers | Now just delegates to `switch-mcp-profile.sh` |
| `switch-mcp-profile.sh` | Switch between pre-defined profiles | Does ALL the work: .mcp.json, .gemini/settings.json, ~/.codex/config.toml |

**Impact**:
- Users confused about which script to run
- Extra indirection adds no value
- Documentation references both, creating confusion

**Recommendation**:
1. Mark `configure_mcp_servers.sh` as deprecated
2. Have `setup-ai.sh` call `switch-mcp-profile.sh` directly
3. Remove `configure_mcp_servers.sh` in v3.0

---

#### ISSUE B: Environment Variable Injection Gap (Critical)

**Problem**: API keys flow through a fragmented path with multiple failure points.

```
.devcontainer/.env (Source of Truth)
        │
        ▼ (NOT automatically exported to shell!)
Shell Environment (MISSING KEYS)
        │
        ▼ (load_env() only works within script session)
switch-mcp-profile.sh
        │
        ▼ (Python subprocess inherits ORIGINAL shell env, not script env!)
.mcp.json (Contains ${LITERAL_PLACEHOLDERS} if keys missing)
```

**Root Cause**: Docker Compose reads `.env` files for container lifecycle, but this does NOT export variables to interactive shells.

**Impact**: ISSUE-019, ISSUE-020, and likely others trace back to this gap.

**Recommendation**:
1. **Immediate**: Add `.env` sourcing to `~/.bashrc` (done in this session)
2. **Short-term**: Have `switch-mcp-profile.sh` validate key substitution succeeded
3. **Long-term**: Consider runtime key injection (never write keys to files)

---

#### ISSUE C: PYTHONPATH Collision (Critical - Now Fixed)

**Problem**: PAL MCP uses relative imports (`from config import ...`), which collide with project `config.py` files when PYTHONPATH includes project directories.

**Root Cause**: uvx runs PAL in the project's environment context, inheriting PYTHONPATH.

**Fix Applied**: Added `"PYTHONPATH": ""` to PAL env block in all profile templates.

**Recommendation**:
1. Document this as a pattern for any MCP server using common module names
2. Consider adding PYTHONPATH clearing to the profile template generation logic

---

#### ISSUE D: Script Naming Inconsistency (Medium)

**Problem**: Scripts use mixed naming conventions.

| Convention | Examples |
|------------|----------|
| `verb-noun.sh` | `setup-ai.sh`, `check-role.sh`, `switch-mcp-profile.sh` |
| `verb_noun.sh` | `configure_mcp_servers.sh`, `install_claude.sh` |
| `test_noun.sh` | `test_installation.sh` |

**Impact**: Cognitive load, tab completion issues, documentation inconsistency.

**Recommendation**: Standardize on `verb-noun.sh` (kebab-case) in v3.0.

---

#### ISSUE E: Multiple Entry Points (Medium)

**Problem**: Users can enter the system through multiple scripts:
- `setup-ai.sh` (full setup)
- `setup_mcp_infrastructure.sh` (MCP only)
- `configure_mcp_servers.sh` (deprecated but still exists)
- `switch-mcp-profile.sh` (profile switching)
- `activate-role.sh` (role activation)

**Impact**:
- No single source of truth for "what do I run?"
- Each entry point has slightly different behavior
- Documentation must cover all entry points

**Recommendation**:
1. Implement unified `sci` CLI (see CLI_IMPROVEMENT_PLAN.md)
2. Keep existing scripts for backward compatibility but deprecate direct usage
3. Document primary workflow as: `sci setup` → `sci role activate` → `sci profile switch`

---

#### ISSUE F: Implicit vs Explicit Dependencies (Medium)

**Problem**: Scripts have implicit dependencies that aren't validated upfront.

| Script | Implicit Dependencies |
|--------|----------------------|
| `switch-mcp-profile.sh` | `uvx`, `npx`, API keys in env, Python 3 |
| `setup_pal.sh` | `uvx`, `git+https` network access |
| `setup_tooluniverse.sh` | `uv`, Python 3.10+ |

**Impact**: Failures occur deep in execution rather than at validation time.

**Recommendation**:
1. Add `validate_dependencies()` function to `common.sh`
2. Call at script entry point, fail fast with actionable message
3. Document all dependencies in script header comments

---

#### ISSUE G: Profile Template Duplication (Low)

**Problem**: Profile templates duplicate common patterns.

Example: All profiles with PAL have the same PAL block:
```json
"pal": {
  "type": "stdio",
  "command": "uvx",
  "args": ["--from", "git+https://...", "pal-mcp-server"],
  "env": { "PYTHONPATH": "", "GEMINI_API_KEY": "...", "OPENAI_API_KEY": "..." }
}
```

**Impact**: When PAL configuration changes, ALL profile templates must be updated (done 4 times in this session).

**Recommendation**:
1. Create `templates/mcp-servers/*.json` with individual server definitions
2. Profile templates reference server definitions
3. Single source of truth for each server configuration

---

## Proposed Architecture (v3.0)

### Simplified Dependency Graph

```
sci (Unified CLI Entry Point)
    │
    ├── lib/common.sh (Utilities + Dependency Validation)
    │
    ├── lib/env-loader.sh (Explicit .env Loading)
    │       └── Validates keys loaded into environment
    │
    ├── lib/mcp-config.sh (MCP Configuration Engine)
    │       ├── templates/mcp-servers/*.json (Server Definitions)
    │       └── templates/mcp-profiles/*.json (Profile Compositions)
    │
    └── lib/role-manager.sh (Role Management)
            └── roles/*.yaml (Role Definitions)
```

### Key Architectural Changes

#### 1. Single Entry Point
- `sci` wrapper calls appropriate library functions
- Existing scripts become thin wrappers for backward compatibility
- Clear workflow: `sci setup` → `sci status` → `sci role` → `sci profile`

#### 2. Explicit Environment Loading
- `lib/env-loader.sh` explicitly loads and validates `.env` files
- Returns error if required keys missing
- No implicit "maybe it's in the environment" behavior

#### 3. Composable Server Definitions
```
templates/
├── mcp-servers/
│   ├── pal.json           # PAL server definition
│   ├── context7.json      # Context7 server definition
│   └── ...
└── mcp-profiles/
    ├── coding.json        # References: ["context7", "pal", "sequential-thinking"]
    └── ...
```

#### 4. Fail-Fast Validation
- All scripts validate dependencies at entry
- Clear error messages with fix instructions
- No deep-stack failures

---

## Implementation Roadmap

### Phase 1: Stabilization (Immediate)
- [x] Fix PYTHONPATH collision in profile templates
- [x] Add `.env` sourcing to bashrc
- [x] Document ISSUE-019, ISSUE-020 in ISSUES.md
- [x] Update TROUBLESHOOTING.md with new cases
- [ ] Add validation in `switch-mcp-profile.sh` that API keys were substituted

### Phase 2: Deprecation (v2.x)
- [ ] Mark `configure_mcp_servers.sh` as deprecated
- [ ] Add deprecation warnings to deprecated scripts
- [ ] Update documentation to use preferred entry points
- [ ] Create migration guide for existing users

### Phase 3: Refactoring (v3.0)
- [ ] Implement `sci` unified CLI (see CLI_IMPROVEMENT_PLAN.md)
- [ ] Refactor server definitions into composable templates
- [ ] Add dependency validation to all scripts
- [ ] Standardize naming convention (kebab-case)

### Phase 4: Simplification (v3.1+)
- [ ] Remove deprecated scripts
- [ ] Consolidate to single profile engine
- [ ] Implement runtime key injection (no keys in files)
- [ ] Add comprehensive test suite for interoperation

---

## Metrics for Success

### Complexity Reduction
| Metric | Current | Target (v3.0) |
|--------|---------|---------------|
| Number of entry-point scripts | 5+ | 1 (`sci`) |
| Lines to update for PAL config change | 4 files | 1 file |
| Implicit dependencies per script | 3-5 | 0 (all validated) |
| Documented troubleshooting cases | 20+ | <10 (prevented by design) |

### User Experience
- [ ] New user can get started with 2 commands
- [ ] Error messages include fix instructions
- [ ] No "works on my machine" issues
- [ ] Tab completion for all commands

---

## Related Documentation

- [CLI_IMPROVEMENT_PLAN.md](CLI_IMPROVEMENT_PLAN.md) - Unified CLI proposal
- [ISSUES.md](ISSUES.md) - Tracked issues
- [TROUBLESHOOTING.md](TROUBLESHOOTING.md) - Common problems and solutions
- [SECURITY.md](SECURITY.md) - API key handling security model

---

## Appendix: Scripts to Consolidate

### Scripts to Keep (Renamed)
| Current | Proposed | Reason |
|---------|----------|--------|
| `switch-mcp-profile.sh` | `lib/mcp-config.sh` | Core configuration engine |
| `activate-role.sh` | `lib/role-manager.sh` | Role management engine |
| `common.sh` | `lib/common.sh` | Shared utilities |

### Scripts to Deprecate
| Script | Replacement | Timeline |
|--------|-------------|----------|
| `configure_mcp_servers.sh` | `sci profile switch` | v2.x deprecation, v3.0 removal |
| `setup-ai.sh` | `sci setup` | v2.x deprecation, v3.0 thin wrapper |
| `setup_mcp_infrastructure.sh` | `sci setup mcp` | v2.x deprecation, v3.0 thin wrapper |

### Scripts to Archive
| Script | Reason | Action |
|--------|--------|--------|
| Individual `install_*.sh` | Rarely run standalone | Move to `lib/installers/` |
| Individual `setup_*.sh` | Called by orchestrator | Move to `lib/mcp-setup/` |

---

## Version History

| Version | Date | Changes |
|---------|------|---------|
| 1.0.0 | 2026-01-18 | Initial architecture review following ISSUE-019/020 resolution |
