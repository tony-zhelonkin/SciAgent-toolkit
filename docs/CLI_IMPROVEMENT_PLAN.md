# CLI Improvement Plan

**Created**: 2025-12-30
**Status**: Proposal / Design Document
**Version**: 1.0

This document outlines a comprehensive plan to improve the SciAgent-toolkit CLI interface, making it more streamlined, intuitive, and user-friendly.

---

## Executive Summary

The current SciAgent-toolkit CLI consists of 10+ standalone bash scripts with inconsistent interfaces, making it difficult for users to discover functionality and understand the system state. This plan proposes a **unified CLI wrapper** (`sci`) that provides a consistent, discoverable interface while preserving backward compatibility with existing scripts.

---

## Current State Analysis

### Existing CLI Scripts

| Script | Purpose | Naming Convention | Arguments Style |
|--------|---------|-------------------|-----------------|
| `setup-ai.sh` | Full project setup | verb-noun | `--flags` |
| `setup_mcp_infrastructure.sh` | MCP installation | verb_noun_noun | `--flags` |
| `activate-role.sh` | Role activation | verb-noun | positional + `--flags` |
| `check-role.sh` | Check active role | verb-noun | `--flags` |
| `switch-mcp-profile.sh` | Profile switching | verb-noun-noun | positional + `--flags` |
| `toggle-mcp.sh` | MCP enable/disable | verb-noun | subcommand + positional |
| `validate-secrets.sh` | Security check | verb-noun | positional |
| `test_installation.sh` | Verify install | verb_noun | none |
| `install_claude.sh` | Install Claude | verb_noun | none |
| `configure_mcp_servers.sh` | Generate configs | verb_noun_noun | `--flags` |

### Identified Pain Points

1. **Discovery Problem**
   - Users must know script names to use them
   - No central help or command listing
   - Tab completion not available

2. **Inconsistent Interface**
   - Mixed naming: `setup-ai.sh` vs `setup_mcp_infrastructure.sh`
   - Mixed arguments: positional vs flags vs subcommands
   - Different output formats and verbosity

3. **Missing Functionality**
   - No unified status command
   - No way to list available profiles/roles
   - No update mechanism
   - No version reporting

4. **Cognitive Load**
   - Users must remember 10+ script names
   - Path awareness required (`./scripts/` vs `./`)
   - No clear workflow guidance

5. **State Visibility**
   - Current role/profile not easily discoverable
   - MCP server status requires starting Claude
   - No single dashboard view

---

## Proposed Solution: Unified `sci` CLI

### Design Principles

1. **Single Entry Point**: One command (`sci`) for all operations
2. **Discoverable**: `sci help` and `sci <cmd> --help` for all commands
3. **Consistent**: Same argument patterns across all subcommands
4. **Progressive Disclosure**: Simple commands for common tasks, flags for advanced use
5. **Backward Compatible**: Existing scripts continue to work

### Command Structure

```
sci <command> [subcommand] [options] [arguments]
```

### Proposed Command Hierarchy

```
sci
├── setup                    # Installation & setup
│   ├── full                 # Full setup (default)
│   ├── minimal              # Minimal setup (skip optional components)
│   └── mcp                  # MCP servers only
│
├── status                   # Show current state (role, profile, servers)
│
├── role                     # Role management
│   ├── list                 # List available roles
│   ├── show [name]          # Show role details
│   ├── activate <name>      # Activate a role
│   └── current              # Show current role (default)
│
├── profile                  # MCP profile management
│   ├── list                 # List available profiles
│   ├── show [name]          # Show profile details
│   ├── switch <name>        # Switch to a profile
│   └── current              # Show current profile (default)
│
├── mcp                      # MCP server management
│   ├── list                 # List configured servers
│   ├── enable <server>      # Enable a server
│   ├── disable <server>     # Disable a server
│   └── test [server]        # Test server connectivity
│
├── check                    # Validation & diagnostics
│   ├── secrets              # Validate secret protection
│   ├── install              # Verify installation
│   └── all                  # Run all checks (default)
│
├── update                   # Update toolkit components
│   ├── self                 # Update SciAgent-toolkit
│   └── mcp                  # Update MCP servers
│
├── doctor                   # Diagnose and fix common issues
│
├── version                  # Show version information
│
└── help [command]           # Show help
```

---

## Detailed Command Specifications

### `sci status` - Dashboard View

**Purpose**: Single command to see everything at a glance

**Output Example**:
```
SciAgent-toolkit v2.1.0
══════════════════════════════════════════════════════════════

Role:     dc-dictionary
Profile:  hybrid-research (~35k tokens)

MCP Servers:
  ✔ context7            connected
  ✔ pal                 connected
  ✔ sequential-thinking connected
  ✘ tooluniverse        failed

Agents: 7 active   Skills: 12 active

Location: /workspaces/DC_Dictionary
══════════════════════════════════════════════════════════════

Run 'sci doctor' to diagnose issues
```

**Flags**:
- `--json`: Output as JSON (for scripting)
- `--quiet`: Only show issues

---

### `sci role` - Role Management

**Subcommands**:

```bash
# List available roles
sci role list
# Output:
#   base           Default bioinformatics analysis role
#   planning       Planning + cross-model consensus
#   sc-atac        Single-cell ATAC-seq role
# * dc-dictionary  DC Dictionary project role (active)

# Show role details
sci role show base
# Output:
#   Name: base
#   Description: Default bioinformatics analysis role
#   MCP Profile: coding
#   Agents (8):
#     • bioinf-librarian
#     • bio-research-visualizer
#     ...
#   Skills (0): none

# Activate a role
sci role activate base
# Output:
#   ✔ Activated role: base
#   ✔ MCP profile: coding
#   → Restart Claude Code to apply changes

# Show current role (default)
sci role
# Output:
#   dc-dictionary
```

---

### `sci profile` - Profile Management

**Subcommands**:

```bash
# List profiles with token estimates
sci profile list
# Output:
#   minimal        ~3k tokens    Sequential Thinking + Context7
#   coding         ~25k tokens   + PAL
#   codebase       ~75k tokens   + Serena
# * hybrid-research ~35k tokens  Claude=coding, Gemini=research (active)
#   research-lite  ~30k tokens   + ToolUniverse (6 tools)
#   research-full  ~50k tokens   + ToolUniverse (14 tools)
#   full           ~100k tokens  All servers

# Show profile details
sci profile show hybrid-research
# Output:
#   Name: hybrid-research
#   Context: ~35k tokens (17.5% of 200k)
#
#   Claude Code servers:
#     • context7
#     • pal
#     • sequential-thinking
#
#   Gemini servers:
#     • tooluniverse (research tools)

# Switch profile
sci profile switch coding
# Output:
#   ✔ Switched to profile: coding
#   ✔ Updated .mcp.json
#   ✔ Updated .gemini/settings.json
#   → Restart Claude Code to apply changes

# Show current profile
sci profile
# Output:
#   hybrid-research
```

---

### `sci mcp` - MCP Server Management

**Subcommands**:

```bash
# List servers and status
sci mcp list
# Output:
#   Server              Status      Profile
#   context7            enabled     all
#   pal                 enabled     coding+
#   sequential-thinking enabled     all
#   serena              disabled    codebase+
#   tooluniverse        enabled     research+

# Enable/disable servers
sci mcp enable serena
sci mcp disable tooluniverse

# Test server
sci mcp test pal
# Output:
#   Testing PAL MCP server...
#   ✔ uvx available
#   ✔ Server starts successfully
#   ✔ API keys configured (GEMINI_API_KEY)
#   ✔ Ready to use

# Test all servers
sci mcp test
```

---

### `sci setup` - Installation

**Subcommands**:

```bash
# Full setup (default)
sci setup
# Equivalent to: ./setup-ai.sh

# Minimal setup (fast)
sci setup minimal
# Equivalent to: ./setup-ai.sh --minimal

# MCP only (no templates/roles)
sci setup mcp
# Equivalent to: ./setup_mcp_infrastructure.sh

# With options
sci setup --skip-codex --skip-gemini
sci setup --force
```

---

### `sci check` - Validation

**Subcommands**:

```bash
# Security validation
sci check secrets
# Equivalent to: ./validate-secrets.sh

# Installation verification
sci check install
# Equivalent to: ./test_installation.sh

# All checks (default)
sci check
# Output:
#   Security Validation
#   ✔ .mcp.json is gitignored
#   ✔ .gemini/settings.json is gitignored
#   ✔ No API keys in staged files
#
#   Installation Verification
#   ✔ Claude Code installed (v2.0.71)
#   ✔ MCP config valid JSON
#   ✔ ToolUniverse environment exists
#   ⚠ Serena not installed (optional)
#
#   4 passed, 0 failed, 1 warning
```

---

### `sci doctor` - Diagnostics & Auto-fix

**Purpose**: Identify and optionally fix common issues

```bash
sci doctor
# Output:
#   Diagnosing SciAgent-toolkit...
#
#   ✔ Claude Code installed
#   ✔ MCP configuration valid
#   ✘ PAL server failing
#     → GEMINI_API_KEY not set in environment
#     → Fix: Add to .env and run 'sci profile switch hybrid-research'
#
#   ⚠ ToolUniverse connection timeout
#     → May need longer startup time
#     → Try: sci mcp test tooluniverse --timeout 60
#
#   Found 1 error, 1 warning
#
#   Run 'sci doctor --fix' to attempt automatic fixes

sci doctor --fix
# Attempts to fix issues automatically where possible
```

---

### `sci version` - Version Information

```bash
sci version
# Output:
#   SciAgent-toolkit  v2.1.0
#   Claude Code       v2.0.71
#   Node.js           v20.12.0
#   uv                v0.5.4
#
#   Active MCP Servers:
#     context7            @upstash/context7-mcp@latest
#     pal                 pal-mcp-server (git)
#     sequential-thinking @modelcontextprotocol/server-sequential-thinking@latest

sci version --short
# Output:
#   2.1.0
```

---

## Implementation Phases

### Phase 1: Foundation (Week 1)

1. **Create `sci` wrapper script**
   - Bash script at `scripts/sci`
   - Command routing to existing scripts
   - Help system infrastructure

2. **Implement core commands**
   - `sci status` (new)
   - `sci version` (new)
   - `sci help` (new)

3. **Add aliases for existing functionality**
   - `sci setup` → `setup-ai.sh`
   - `sci role activate` → `activate-role.sh`
   - `sci role current` → `check-role.sh`
   - `sci profile switch` → `switch-mcp-profile.sh`

### Phase 2: Enhanced Commands (Week 2)

1. **Role management**
   - `sci role list` - Parse `roles/*.yaml`
   - `sci role show` - Display role details

2. **Profile management**
   - `sci profile list` - Parse `templates/mcp-profiles/*.mcp.json`
   - `sci profile show` - Display profile details

3. **MCP management**
   - `sci mcp list` - Parse `.mcp.json` and settings
   - `sci mcp test` - Test server connectivity

### Phase 3: Diagnostics & Polish (Week 3)

1. **Doctor command**
   - Issue detection logic
   - Auto-fix capabilities

2. **Check command consolidation**
   - Unified validation runner

3. **Output formatting**
   - Consistent color scheme
   - JSON output option
   - Quiet mode

### Phase 4: Quality of Life (Week 4)

1. **Tab completion**
   - Bash completion script
   - Zsh completion script

2. **Update mechanism**
   - `sci update self`
   - `sci update mcp`

3. **Documentation**
   - Man page generation
   - Interactive tutorial (`sci tutorial`)

---

## Backward Compatibility

### Strategy

1. **Preserve all existing scripts** - They continue to work unchanged
2. **`sci` is a wrapper** - Calls existing scripts internally
3. **Deprecation warnings** - When scripts are called directly (optional)
4. **Migration path** - Document how to update workflows

### Script Mapping

| Old Command | New Command |
|-------------|-------------|
| `./scripts/setup-ai.sh` | `sci setup` |
| `./scripts/setup-ai.sh --minimal` | `sci setup minimal` |
| `./scripts/activate-role.sh base` | `sci role activate base` |
| `./scripts/check-role.sh` | `sci role` or `sci role current` |
| `./scripts/switch-mcp-profile.sh coding` | `sci profile switch coding` |
| `./scripts/toggle-mcp.sh enable pal` | `sci mcp enable pal` |
| `./scripts/validate-secrets.sh` | `sci check secrets` |
| `./scripts/test_installation.sh` | `sci check install` |

---

## User Experience Improvements

### Discoverability

1. **Contextual help**: `sci help <command>` shows examples
2. **Command suggestions**: "Did you mean `sci role activate`?"
3. **Tab completion**: Complete commands, roles, profiles, servers

### Feedback

1. **Progress indicators**: Spinners for long operations
2. **Color coding**: Green=success, Yellow=warning, Red=error
3. **Actionable messages**: Always tell user what to do next

### Error Handling

1. **Clear error messages**: What went wrong
2. **Suggested fixes**: How to resolve
3. **Links to docs**: Where to learn more

### Examples

```bash
# Good error message
$ sci role activate nonexistent
Error: Role 'nonexistent' not found

Available roles:
  • base
  • planning
  • sc-atac
  • dc-dictionary

Run 'sci role list' to see all roles with descriptions.

# Good success message
$ sci profile switch coding
✔ Switched to profile: coding

  Enabled servers:
    • context7
    • pal
    • sequential-thinking

  Context usage: ~25k tokens (12.5% of 200k)

→ Restart Claude Code to apply: exit && claude
```

---

## Configuration

### Global Settings File

Location: `~/.config/sciagent/config.yaml`

```yaml
# User preferences
defaults:
  project_dir: auto           # 'auto' or specific path
  output_format: pretty       # 'pretty', 'json', or 'quiet'
  color: auto                 # 'auto', 'always', or 'never'

# Command defaults
setup:
  skip_codex: false
  skip_gemini: false

# Aliases
aliases:
  s: status
  r: role
  p: profile
```

### Environment Variables

| Variable | Purpose | Default |
|----------|---------|---------|
| `SCI_PROJECT_DIR` | Override project directory | Current directory |
| `SCI_OUTPUT` | Output format | `pretty` |
| `SCI_COLOR` | Color mode | `auto` |
| `SCI_CONFIG` | Config file path | `~/.config/sciagent/config.yaml` |

---

## Technical Considerations

### Language Choice

**Recommendation**: Bash for initial implementation

**Rationale**:
- Consistent with existing scripts
- No additional dependencies
- Easy integration with existing code
- Portable across Unix systems

**Future consideration**: Rewrite in Python or Go for:
- Better argument parsing
- Easier testing
- Cross-platform support (Windows)

### Script Location

```
scripts/
├── sci                      # Main entry point (NEW)
├── sci-completion.bash      # Bash completion (NEW)
├── sci-completion.zsh       # Zsh completion (NEW)
├── lib/                     # Shared functions (NEW)
│   ├── output.sh            # Output formatting
│   ├── config.sh            # Configuration handling
│   └── utils.sh             # Common utilities
├── setup-ai.sh              # Existing (unchanged)
├── activate-role.sh         # Existing (unchanged)
└── ...                      # Other existing scripts
```

### Installation

```bash
# Add to PATH (in .bashrc)
export PATH="$PATH:/path/to/SciAgent-toolkit/scripts"

# Or create symlink
ln -s /path/to/SciAgent-toolkit/scripts/sci ~/.local/bin/sci

# Enable completion
source /path/to/SciAgent-toolkit/scripts/sci-completion.bash
```

---

## Success Metrics

### Usability

- [ ] User can discover available commands without docs
- [ ] Common workflows require fewer keystrokes
- [ ] Error messages are actionable
- [ ] New users can get started with `sci setup && sci status`

### Adoption

- [ ] Existing users can migrate without breaking workflows
- [ ] Documentation updated to use new commands
- [ ] Examples in QUICKSTART.md use `sci` commands

### Maintainability

- [ ] Adding new commands is straightforward
- [ ] Shared code reduces duplication
- [ ] Tests cover critical paths

---

## Open Questions

1. **Should `sci` be the default name?**
   - Alternatives: `sag` (SciAgent), `sat` (SciAgent-toolkit), `agent`
   - `sci` is short, memorable, and fits "scientific"

2. **How to handle project vs global scope?**
   - Some commands are project-specific (role, profile)
   - Some are global (update self, version)
   - Proposal: Auto-detect, with `--global` flag override

3. **Should we support Windows?**
   - Current: Unix-only (bash scripts)
   - Future: Consider Python rewrite for cross-platform

4. **Interactive mode?**
   - Consider `sci interactive` for guided setup
   - Menu-driven interface for complex workflows

---

## Related Documentation

- [ARCHITECTURE.md](ARCHITECTURE.md) - System architecture
- [QUICKSTART.md](QUICKSTART.md) - Getting started guide
- [ISSUES.md](ISSUES.md) - Known issues
- [TROUBLESHOOTING.md](TROUBLESHOOTING.md) - Common problems

---

## Appendix: Full Command Reference

```
sci - SciAgent-toolkit unified CLI

USAGE:
    sci <command> [subcommand] [options] [arguments]

COMMANDS:
    setup [variant]         Install and configure toolkit
        full                Full setup with all components (default)
        minimal             Skip optional components (Codex, Gemini, Serena)
        mcp                 MCP servers only

    status                  Show current state (role, profile, servers)
        --json              Output as JSON
        --quiet             Only show issues

    role <subcommand>       Manage roles
        list                List available roles
        show [name]         Show role details
        activate <name>     Activate a role
        current             Show current role (default)

    profile <subcommand>    Manage MCP profiles
        list                List available profiles
        show [name]         Show profile details
        switch <name>       Switch to a profile
        current             Show current profile (default)

    mcp <subcommand>        Manage MCP servers
        list                List configured servers
        enable <server>     Enable a server
        disable <server>    Disable a server
        test [server]       Test server connectivity

    check [type]            Run validation checks
        secrets             Validate secret protection
        install             Verify installation
        all                 Run all checks (default)

    doctor                  Diagnose and fix issues
        --fix               Attempt automatic fixes

    update <target>         Update components
        self                Update SciAgent-toolkit
        mcp                 Update MCP servers

    version                 Show version information
        --short             Version number only

    help [command]          Show help

GLOBAL OPTIONS:
    --project-dir <dir>     Override project directory
    --quiet, -q             Minimal output
    --json                  JSON output (where supported)
    --no-color              Disable colored output
    --version, -v           Show version
    --help, -h              Show help

EXAMPLES:
    sci setup                       # Full installation
    sci status                      # Dashboard view
    sci role activate base          # Switch to base role
    sci profile switch coding       # Use coding profile
    sci mcp test                    # Test all MCP servers
    sci doctor --fix                # Auto-fix issues

ENVIRONMENT:
    SCI_PROJECT_DIR     Override project directory
    SCI_OUTPUT          Output format (pretty|json|quiet)
    SCI_COLOR           Color mode (auto|always|never)

For more information: https://github.com/tony-zhelonkin/SciAgent-toolkit
```
