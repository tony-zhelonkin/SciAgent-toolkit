#!/usr/bin/env bash
#
# PAL MCP Server Setup Script
#
# Sets up the PAL MCP server for collaboration, planning, and code analysis.
# Works with Claude Code.
#
# Requirements:
#   - uvx (from uv package manager)

set -euo pipefail

# Get script directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Source common utilities
if [ -f "${SCRIPT_DIR}/../common.sh" ]; then
    source "${SCRIPT_DIR}/../common.sh"
else
    echo "Error: common.sh not found in ${SCRIPT_DIR}/.."
    exit 1
fi

log_info "Setting up PAL MCP Server..."

# ============================================================================
# 1. Check prerequisites
# ============================================================================
log_info "Checking prerequisites..."

# Check for git
if ! command -v git &>/dev/null; then
    log_error "git not found. Installing git..."

    OS_TYPE="$(uname -s)"
    case "$OS_TYPE" in
        Linux*)
            if sudo apt-get update && sudo apt-get install -y git &>/dev/null; then
                log_ok "git installed successfully"
            else
                log_error "Failed to install git"
                exit 1
            fi
            ;;
        Darwin*)
            if command -v brew &>/dev/null && brew install git &>/dev/null; then
                log_ok "git installed via Homebrew"
            else
                log_error "git not available. Please install git manually"
                log_info "On macOS, install Xcode Command Line Tools: xcode-select --install"
                exit 1
            fi
            ;;
        *)
            log_error "git is required but not installed"
            log_info "Please install git manually for your platform"
            exit 1
            ;;
    esac
else
    GIT_VERSION=$(git --version 2>/dev/null || echo "unknown")
    log_ok "git already installed (${GIT_VERSION})"
fi

# Check for uv/uvx
if ! command -v uvx &>/dev/null; then
    log_error "uvx not found. Installing uv package manager..."

    # Detect OS
    OS_TYPE="$(uname -s)"
    case "$OS_TYPE" in
        Linux*)     PLATFORM=Linux;;
        Darwin*)    PLATFORM=macOS;;
        *)          PLATFORM=Unknown;;
    esac

    if [ "$PLATFORM" = "macOS" ]; then
        # On macOS, use Homebrew or pip3 --user
        if command -v brew &>/dev/null; then
            if brew install uv &>/dev/null; then
                log_ok "uv installed via Homebrew"
            else
                log_error "Failed to install uv via Homebrew"
                exit 1
            fi
        elif pip3 install --user uv &>/dev/null; then
            log_ok "uv installed via pip3 --user"
            export PATH="$HOME/Library/Python/$(python3 -c 'import sys; print(f"{sys.version_info.major}.{sys.version_info.minor}")')/bin:$PATH"
        else
            log_error "Failed to install uv"
            exit 1
        fi
    else
        # On Linux, use system pip3 with sudo
        if sudo pip3 install uv &>/dev/null; then
            log_ok "uv package installed successfully"
        else
            log_error "Failed to install uv package"
            exit 1
        fi
    fi

    # Verify installation
    if ! command -v uvx &>/dev/null; then
        log_error "uvx command not found after uv installation"
        exit 1
    fi
else
    UVX_VERSION=$(uvx --version 2>/dev/null | head -1 || echo "unknown")
    log_ok "uvx already installed (${UVX_VERSION})"
fi

# ============================================================================
# 2. Test PAL MCP server
# ============================================================================
log_info "Testing PAL MCP server..."

# Configure Git to use HTTPS instead of SSH for GitHub
# This ensures users without SSH keys can still install PAL
export GIT_TERMINAL_PROMPT=0
export GIT_SSH_COMMAND="ssh -o BatchMode=yes -o StrictHostKeyChecking=no"

# Temporarily configure Git to prefer HTTPS over SSH
git config --global --add url."https://github.com/".insteadOf git@github.com: 2>/dev/null || true
git config --global --add url."https://".insteadOf ssh:// 2>/dev/null || true

log_info "First-time download may take a moment..."

if timeout 30 uvx --from git+https://github.com/BeehiveInnovations/pal-mcp-server.git pal-mcp-server --help &>/dev/null; then
    log_ok "PAL MCP server can start"
else
    log_warn "Could not verify PAL startup (may need first-time download or API keys)"
    log_info "This is normal on first run. PAL will download automatically when first used."
fi

# Clean up temporary environment variables
unset GIT_TERMINAL_PROMPT
unset GIT_SSH_COMMAND

# Revert Git config changes (only for this session)
git config --global --unset-all url."https://github.com/".insteadOf 2>/dev/null || true
git config --global --unset-all url."https://".insteadOf 2>/dev/null || true

# ============================================================================
# Summary
# ============================================================================
log_ok "PAL MCP Server setup complete"
log_info ""
log_info "PAL provides:"
log_info "  - Collaboration & Planning (chat, thinkdeep, planner)"
log_info "  - Code Analysis & Quality (debug, codereview)"
log_info ""
log_info "PAL will be configured in the project .mcp.json file"
