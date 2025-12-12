# Changelog

All notable changes to the SciAgent-toolkit project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Fixed
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
