#!/usr/bin/env bash
#
# Claude Code Installation Script
#
# Installs Claude Code using the native installer method.
# Works on macOS, Linux, and WSL.

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

# Check if running as root (should not be)
if [ "$(id -u)" -eq 0 ]; then
    log_error "This script should not be run as root."
    exit 1
fi

log_info "Installing Claude Code..."

# Check if already installed
if command -v claude &>/dev/null; then
    CLAUDE_VERSION=$(claude --version 2>/dev/null | head -1 || echo "unknown")
    log_ok "Claude Code already installed (version: ${CLAUDE_VERSION})"

    # Ensure PATH is set even if already installed
    if ! grep -q 'export PATH="$HOME/.local/bin:$PATH"' "$HOME/.bashrc" 2>/dev/null; then
        echo 'export PATH="$HOME/.local/bin:$PATH"' >> "$HOME/.bashrc"
        log_info "Added ~/.local/bin to PATH in .bashrc"
    fi
    export PATH="$HOME/.local/bin:$PATH"

    log_ok "Claude Code installation verified"
    exit 0
fi

# Install using native installer
log_info "Downloading and installing Claude Code (latest version)..."
log_info "This may take 1-2 minutes (downloading ~50MB binary)..."

# Use timeout to prevent hanging indefinitely
# 300 seconds (5 minutes) should be enough for most connections
INSTALL_TIMEOUT=300

set +e  # Temporarily allow command failure
timeout $INSTALL_TIMEOUT bash -c 'curl -fsSL https://claude.ai/install.sh | bash -s latest'
INSTALL_EXIT_CODE=$?
set -e  # Re-enable exit on error

if [ $INSTALL_EXIT_CODE -eq 0 ]; then
    log_ok "Claude Code installed successfully"

    # Add ~/.local/bin to PATH permanently for bash
    if ! grep -q 'export PATH="$HOME/.local/bin:$PATH"' "$HOME/.bashrc"; then
        echo 'export PATH="$HOME/.local/bin:$PATH"' >> "$HOME/.bashrc"
        log_info "Added ~/.local/bin to PATH in .bashrc"
    fi

    # Ensure ~/.local/bin is in PATH for this session
    export PATH="$HOME/.local/bin:$PATH"

    # Verify installation
    if command -v claude &>/dev/null; then
        CLAUDE_VERSION=$(claude --version 2>/dev/null | head -1 || echo "installed")
        log_ok "Claude Code version: ${CLAUDE_VERSION}"

        # Skip diagnostics - claude doctor can hang in container environments
        # Run manually with: claude doctor
        log_ok "Skipping diagnostics (run 'claude doctor' manually if needed)"
    else
        log_error "Claude Code installation failed - command not found"
        exit 1
    fi
else
    if [ $INSTALL_EXIT_CODE -eq 124 ]; then
        log_error "Claude Code installation timed out after ${INSTALL_TIMEOUT} seconds"
        log_info "This might be a network issue. Try again or install manually:"
        log_info "  curl -fsSL https://claude.ai/install.sh | bash -s latest"
    else
        log_error "Failed to install Claude Code (exit code: $INSTALL_EXIT_CODE)"
    fi
    log_info "You can also try installing via npm: npm install -g @anthropic-ai/claude-code"
    exit 1
fi

log_ok "Claude Code installation complete"
log_info ""
log_info "Next steps:"
log_info "  1. Restart your terminal or run: source ~/.bashrc"
log_info "  2. Run 'claude' to start Claude Code"
log_info "  3. Sign in with your Anthropic account"
