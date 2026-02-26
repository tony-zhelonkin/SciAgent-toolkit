# Changelog

All notable changes to the SciAgent-toolkit project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Fixed
- **MCP addon tool permissions**: `manage-addon.sh` now reads `tool_permissions` from addon templates and merges them into `settings.local.json` `permissions.allow`. Previously, addon tools were denied in subagents (Task tool) because subagents can't prompt interactively for permission.
  - `manage-addon.sh`: `update_settings_local()` strips stale `mcp__*` entries and rebuilds from enabled addons
  - `switch-mcp-profile.sh`: Settings generation includes addon tool permissions (survives profile switches)
  - `notebook-tools.addon.json`: Added `tool_permissions` array with all 11 tool names
  - `jupyter.addon.json`: Added empty `tool_permissions` array for future use

## [2.0.0] - 2025-12-16

### Added
- **Role System**: New declarative role-based configuration for agents and skills
  - `roles/base.yaml` - Default bioinformatics analysis role
  - `scripts/activate-role.sh` - Role activator script (symlinks agents/skills to `.claude/`)
  - `skills/` directory for custom skills
- **Template System**: Centralized AI context templates in `templates/vendor/`
  - `CLAUDE.md.template` - Claude Code project instructions
  - `GEMINI.md.template` - Gemini CLI project instructions
  - `AGENTS.md.template` - Universal AI rules for all agents
  - `context.md.template` - Scientific project context
  - `analysis_config.yaml.template` - Analysis parameters
- **Enhanced Profile System** (`switch-mcp-profile.sh`):
  - API key substitution (`${GEMINI_API_KEY}`, `${OPENAI_API_KEY}`, `${CONTEXT7_API_KEY}`)
  - `validate_profile()` function for dependency checking
  - Profile validation before switching
- **setup-ai.sh Enhancements**:
  - Template installation with placeholder substitution
  - Automatic role activation (`activate-role.sh base`)
  - Creates `02_analysis/config/analysis_config.yaml`
- **Architecture Tests** (`docker/test/Dockerfile.architecture-test`):
  - Validates role system, templates, profile switching
  - 8 sub-tests covering all modularization changes
- **Gemini CLI Test** (`docker/test/Dockerfile.gemini-test`)
- **Full test suite** now includes 11 tests (up from 10)

### Changed
- **scbio-docker Integration**: Directory renamed from `01_scripts/` to `01_modules/`
- **Template References**: All templates now reference `01_modules` (not `01_scripts`)
- **research-full.mcp.json**: Added explicit `--include-tools` with 14 curated tools to prevent context overflow

### Fixed
- **MCP Configuration**: Fixed incorrect arguments in `full.mcp.json` and `research-full.mcp.json` templates that caused `switch_mcp full` to fail.
  - Removed incorrect `uv` arguments (`--directory`, `run`) that were passed to the direct binary executable.
  - This ensures `tooluniverse` works correctly in `full` and `research-full` profiles, consistent with `research-lite`.
- **Dev Container npm EACCES**: Fixed npm global install failures in dev containers where nvm is installed to a read-only location (e.g., `/opt/nvm/` owned by root).
  - Added `ensure_npm_writable_prefix()` function to detect read-only npm prefix and auto-configure user-local prefix (`~/.npm-global`)
  - Added `ensure_npm_global_path()` function to ensure previously installed npm binaries are in PATH
  - Updated `install_gemini.sh` and `install_codex.sh` to call these functions before npm install
  - PATH updates are now persisted to `~/.bashrc` automatically
- **Prefix Conflict**: Fixed an issue where `nvm` usage was conflicting with a hardcoded `npm` prefix setting in `common.sh`, causing warnings and potential environment issues.
  - Removed `configure_npm_prefix` call from `ensure_nvm` in `scripts/common.sh`.
  - Removed `prefix` setting from `~/.npmrc`.
  - Removed `~/.npm-global/bin` from `~/.bashrc` to prevent shadowing of `nvm` managed binaries.
- **Claude Code Settings**: Fixed invalid permission glob patterns in `.claude/settings.local.json` generation.
  - Updated `scripts/switch-mcp-profile.sh` to use the correct `Bash(for f in :*)` wildcard syntax instead of invalid `Bash(for f in :*.ext)` patterns.
  - This resolves "Invalid Settings" warnings in `claude doctor` and ensures proper file globbing permissions.
- **MCP Profile Switching**: Fixed `switch-mcp-profile.sh` to correctly generate the `permissions` block in `settings.local.json`, ensuring MCP servers load correctly in Claude Code.

### Changed
- **Dependencies**: Updated `scripts/common.sh` to respect `nvm` environment management and avoid forcing global `npm` prefixes when `nvm` is active.
