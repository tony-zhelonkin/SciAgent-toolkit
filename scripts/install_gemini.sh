#!/usr/bin/env bash
#
# Google Gemini CLI Installation Script
#
# Installs Google Gemini CLI using npm.
# Works on macOS, Linux, and WSL.
#

set -euo pipefail

# Get script directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Source common utilities
if [ -f "${SCRIPT_DIR}/common.sh" ]; then
    source "${SCRIPT_DIR}/common.sh"
else
    echo "Error: common.sh not found in ${SCRIPT_DIR}"
    exit 1
fi

# Detect OS
PLATFORM=$(detect_os)

log_info "Platform: $PLATFORM"

# Check if running as root (should not be)
if [ "$(id -u)" -eq 0 ]; then
    log_error "This script should not be run as root."
    exit 1
fi

log_info "Installing Google Gemini CLI..."

# First, ensure any previously installed npm-global binaries are in PATH
ensure_npm_global_path

# Check if already installed (after ensuring PATH includes npm-global)
if command -v gemini &>/dev/null; then
    GEMINI_VERSION=$(gemini --version 2>/dev/null | head -1 || echo "unknown")
    log_ok "Gemini CLI already installed (version: ${GEMINI_VERSION})"
    log_ok "Gemini CLI installation verified"
    exit 0
fi

# Try npm (required)
ensure_nvm

# Ensure npm prefix is writable (handles read-only container environments)
ensure_npm_writable_prefix

if command -v npm &>/dev/null; then
    log_info "Installing Gemini CLI via npm..."
    if npm install -g @google/gemini-cli; then
        log_ok "Gemini CLI installed via npm"
        # Ensure the new binary is in PATH
        ensure_npm_global_path
    else
        log_error "Failed to install Gemini CLI via npm"
        exit 1
    fi
else
    log_error "Cannot install Gemini CLI: npm required"
    log_info "Please install Node.js/npm first:"
    log_info "  - On Linux: curl -fsSL https://deb.nodesource.com/setup_20.x | sudo -E bash - && sudo apt-get install -y nodejs"
    log_info "  - On macOS: brew install node"
    exit 1
fi

# Verify installation
if command -v gemini &>/dev/null; then
    GEMINI_VERSION=$(gemini --version 2>/dev/null | head -1 || echo "installed")
    log_ok "Gemini CLI version: ${GEMINI_VERSION}"
else
    log_error "Gemini CLI installation failed - command not found"
    exit 1
fi

log_ok "Gemini CLI installation complete"
log_info ""
log_info "Next steps:"
log_info "  1. Run 'gemini' to start Gemini CLI"
log_info "  2. Sign in with your Google account"
log_info "  3. Configure MCP servers (check 'gemini mcp --help')"
