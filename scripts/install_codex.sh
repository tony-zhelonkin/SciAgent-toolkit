#!/usr/bin/env bash
#
# GPT Codex CLI Installation Script
#
# Installs GPT Codex CLI using npm or Homebrew.
# Works on macOS, Linux, and WSL.

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

# Detect OS
OS_TYPE="$(uname -s)"
case "$OS_TYPE" in
    Linux*)     PLATFORM=Linux;;
    Darwin*)    PLATFORM=macOS;;
    *)          PLATFORM=Unknown;;
esac

log_info "Platform: $PLATFORM"

# Check if running as root (should not be)
if [ "$(id -u)" -eq 0 ]; then
    log_error "This script should not be run as root."
    exit 1
fi

log_info "Installing GPT Codex CLI..."

# Check if already installed
if command -v codex &>/dev/null; then
    CODEX_VERSION=$(codex --version 2>/dev/null | head -1 || echo "unknown")
    log_ok "Codex CLI already installed (version: ${CODEX_VERSION})"
    log_ok "Codex CLI installation verified"
    exit 0
fi

# Try npm first (recommended)
if command -v npm &>/dev/null; then
    log_info "Installing Codex CLI via npm..."
    if npm install -g @openai/codex &>/dev/null; then
        log_ok "Codex CLI installed via npm"
    else
        log_error "Failed to install Codex CLI via npm"
        exit 1
    fi

# Try Homebrew on macOS
elif [ "$PLATFORM" = "macOS" ] && command -v brew &>/dev/null; then
    log_info "Installing Codex CLI via Homebrew..."
    if brew install codex &>/dev/null; then
        log_ok "Codex CLI installed via Homebrew"
    else
        log_error "Failed to install Codex CLI via Homebrew"
        exit 1
    fi

# No installation method available
else
    log_error "Cannot install Codex CLI: npm or Homebrew required"
    log_info "Please install Node.js/npm first:"
    log_info "  - On Linux: curl -fsSL https://deb.nodesource.com/setup_20.x | sudo -E bash - && sudo apt-get install -y nodejs"
    log_info "  - On macOS: brew install node"
    exit 1
fi

# Verify installation
if command -v codex &>/dev/null; then
    CODEX_VERSION=$(codex --version 2>/dev/null | head -1 || echo "installed")
    log_ok "Codex CLI version: ${CODEX_VERSION}"
else
    log_error "Codex CLI installation failed - command not found"
    exit 1
fi

log_ok "Codex CLI installation complete"
log_info ""
log_info "Next steps:"
log_info "  1. Run 'codex' to start Codex CLI"
log_info "  2. Sign in with ChatGPT (recommended) or use OPENAI_API_KEY"
log_info "  3. Configure MCP servers in ~/.codex/config.toml"
log_info ""
log_info "Authentication options:"
log_info "  - Sign in with ChatGPT: codex (then select sign-in option)"
log_info "  - Use API key: export OPENAI_API_KEY='your-key-here' && codex"
