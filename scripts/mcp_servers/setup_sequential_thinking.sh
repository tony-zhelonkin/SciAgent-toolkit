#!/usr/bin/env bash
#
# Sequential Thinking MCP Server Setup Script
#
# Sets up the Sequential Thinking MCP server for structured reasoning.
# Works with Claude Code.
#
# Requirements:
#   - npx (from Node.js/npm)

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

log_info "Setting up Sequential Thinking MCP Server..."

# ============================================================================
# 1. Check prerequisites
# ============================================================================
log_info "Checking prerequisites..."

# Check for Node.js/npx
if ! command -v npx &>/dev/null; then
    log_error "npx not found. Installing Node.js..."

    # Detect OS
    OS_TYPE="$(uname -s)"
    case "$OS_TYPE" in
        Linux*)     PLATFORM=Linux;;
        Darwin*)    PLATFORM=macOS;;
        *)          PLATFORM=Unknown;;
    esac

    if [ "$PLATFORM" = "macOS" ]; then
        # On macOS, use Homebrew
        if command -v brew &>/dev/null; then
            if brew install node@20 &>/dev/null || brew install node &>/dev/null; then
                log_ok "Node.js installed via Homebrew"
            else
                log_error "Failed to install Node.js via Homebrew"
                exit 1
            fi
        else
            log_error "Homebrew not found. Please install from https://brew.sh"
            exit 1
        fi
    else
        # On Linux, use NodeSource
        log_info "Installing Node.js 20.x from NodeSource..."
        if curl -fsSL https://deb.nodesource.com/setup_20.x | sudo -E bash - &>/dev/null && \
           sudo apt-get install -y nodejs &>/dev/null; then
            log_ok "Node.js installed successfully"
        else
            log_error "Failed to install Node.js"
            exit 1
        fi
    fi

    # Verify installation
    if ! command -v npx &>/dev/null; then
        log_error "npx command not found after Node.js installation"
        exit 1
    fi
else
    NODE_VERSION=$(node --version 2>/dev/null || echo "unknown")
    NPM_VERSION=$(npm --version 2>/dev/null || echo "unknown")
    log_ok "Node.js already installed (node: ${NODE_VERSION}, npm: ${NPM_VERSION})"
fi

# ============================================================================
# 2. Test Sequential Thinking MCP server
# ============================================================================
log_info "Testing Sequential Thinking MCP server..."

if timeout 10 npx -y @modelcontextprotocol/server-sequential-thinking --help &>/dev/null; then
    log_ok "Sequential Thinking MCP server can start"
else
    log_warn "Could not verify server startup (may need first-time download)"
fi

# ============================================================================
# Summary
# ============================================================================
log_ok "Sequential Thinking MCP Server setup complete"
log_info ""
log_info "Sequential Thinking provides:"
log_info "  - Structured reasoning for complex decisions"
log_info "  - Step-by-step problem solving"
log_info "  - Multi-step technical analysis"
log_info ""
log_info "Sequential Thinking will be configured in the project .mcp.json file"
