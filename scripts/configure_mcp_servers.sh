#!/usr/bin/env bash
#
# MCP Server Configuration Script
#
# Automatically configures all installed MCP servers using the unified profile system.
# Delegates completely to 'switch-mcp-profile.sh' to ensure consistency across Claude, Gemini, and Codex.
#
# Usage:
#   ./configure_mcp_servers.sh [OPTIONS]
#
# Options:
#   --project-dir DIR      Project directory for .mcp.json (default: current dir)
#   --force                Force re-application of the default profile
#   --dry-run              (No-op, retained for compatibility)
#   --help                 Show this help message

set -euo pipefail

# Get script directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="${PWD}"
FORCE=false

# Source common utilities
if [ -f "${SCRIPT_DIR}/common.sh" ]; then
    source "${SCRIPT_DIR}/common.sh"
else
    echo "Error: common.sh not found in ${SCRIPT_DIR}"
    exit 1
fi

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --project-dir)
            PROJECT_DIR="$2"
            shift 2
            ;;
        --force)
            FORCE=true
            shift
            ;;
        --dry-run)
            log_info "Dry run: Would apply default profile to ${PROJECT_DIR}"
            exit 0
            ;;
        --help)
            grep '^#' "$0" | grep -v '#!/usr/bin/env' | sed 's/^# //' | sed 's/^#//'
            exit 0
            ;;
        *)
            log_error "Unknown option: $1"
            exit 1
            ;;
    esac
done

log_info "Configuring MCP servers..."
log_info "Project directory: ${PROJECT_DIR}"

# ----------------------------------------------------------------------------
# Delegate to Profile Switcher
# ----------------------------------------------------------------------------
SWITCH_SCRIPT="${SCRIPT_DIR}/switch-mcp-profile.sh"

if [ ! -f "${SWITCH_SCRIPT}" ]; then
    log_error "Critical: Profile switcher not found at ${SWITCH_SCRIPT}"
    exit 1
fi

# Load env variables safely
load_env "${PROJECT_DIR}"

# Determine which profile to apply
# If we are forcing, or if no config exists, we apply the default 'coding' profile.
# This ensures that a fresh install gets a working "coding" setup (Claude+Gemini+Codex).

if [ "$FORCE" = true ] || [ ! -f "${PROJECT_DIR}/.mcp.json" ]; then
    log_info "Applying default 'coding' profile..."
    
    if "${SWITCH_SCRIPT}" coding --project-dir "${PROJECT_DIR}"; then
        log_ok "Configuration complete!"
        log_info "Applied 'coding' profile to Claude, Gemini, and Codex."
    else
        log_error "Failed to apply profile."
        exit 1
    fi
else
    log_info "Configuration already exists. Skipping default profile application."
    log_info "Use --force to overwrite with default 'coding' profile."
fi