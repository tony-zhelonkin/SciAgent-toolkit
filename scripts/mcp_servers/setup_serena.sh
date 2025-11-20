#!/usr/bin/env bash
#
# Serena MCP Server Setup Script
#
# Sets up the Serena MCP server for semantic code search and editing.
# Works with Claude Code.
#
# Requirements:
#   - uvx (from uv package manager)

set -euo pipefail

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

log_info()  { echo -e "${BLUE}[INFO]${NC} $*"; }
log_ok()    { echo -e "${GREEN}[OK]${NC} $*"; }
log_warn()  { echo -e "${YELLOW}[WARN]${NC} $*"; }
log_error() { echo -e "${RED}[ERROR]${NC} $*"; }

log_info "Setting up Serena MCP Server..."

# ============================================================================
# 1. Check prerequisites
# ============================================================================
log_info "Checking prerequisites..."

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
# 2. Test Serena MCP server
# ============================================================================
log_info "Testing Serena MCP server..."

if timeout 10 uvx --from git+https://github.com/oraios/serena serena --help &>/dev/null; then
    log_ok "Serena MCP server can start"
else
    log_warn "Could not verify Serena startup (may need first-time download)"
fi

# ============================================================================
# Summary
# ============================================================================
log_ok "Serena MCP Server setup complete"
log_info ""
log_info "Serena provides:"
log_info "  - Semantic code search and editing"
log_info "  - Symbol-level understanding"
log_info "  - Intelligent refactoring"
log_info ""
log_info "Serena will be configured in the project .mcp.json file"
