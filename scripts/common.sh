#!/usr/bin/env bash
#
# Common utilities for SciAgent-toolkit scripts
#

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
CYAN='\033[0;36m'
MAGENTA='\033[0;35m'
NC='\033[0m' # No Color

# Logging functions
log_info()  { echo -e "${BLUE}[INFO]${NC} $*"; }
log_ok()    { echo -e "${GREEN}[OK]${NC} $*"; }
log_warn()  { echo -e "${YELLOW}[WARN]${NC} $*"; }
log_error() { echo -e "${RED}[ERROR]${NC} $*"; }
log_step()  { echo -e "${CYAN}[STEP]${NC} $*"; }
separator() { echo -e "${MAGENTA}======================================${NC}"; echo -e "${MAGENTA}==== $* ${NC}"; echo -e "${MAGENTA}======================================${NC}"; }

# OS Detection
detect_os() {
    local os_type
    os_type="$(uname -s)"
    case "$os_type" in
        Linux*)     echo "Linux";;
        Darwin*)    echo "macOS";;
        *)          echo "Unknown";;
    esac
}

# Dependency Checks
ensure_command() {
    local cmd="$1"
    local msg="${2:-$cmd is required but not installed.}"
    if ! command -v "$cmd" &>/dev/null; then
        log_error "$msg"
        return 1
    fi
    return 0
}

ensure_node() {
    if ! command -v node &>/dev/null; then
        log_warn "Node.js not found. Attempting to install..."
        # This is a placeholder. Real installation might be complex.
        # For now, we just warn and return error if strict.
        return 1
    fi
    local version
    version=$(node --version)
    log_info "Node.js version: $version"
}

# Ensure NVM is installed and loaded
ensure_nvm() {
    export NVM_DIR="$HOME/.nvm"

    # Try to load nvm if not loaded
    if ! command -v nvm &>/dev/null; then
        if [ -s "$NVM_DIR/nvm.sh" ]; then
            set +u  # Disable nounset for nvm script
            \. "$NVM_DIR/nvm.sh"
            set -u  # Re-enable nounset
        fi
    fi

    # Install if still missing
    if ! command -v nvm &>/dev/null; then
        log_info "Installing nvm..."
        curl -o- https://raw.githubusercontent.com/nvm-sh/nvm/v0.40.3/install.sh | bash
        if [ -s "$NVM_DIR/nvm.sh" ]; then
            set +u
            \. "$NVM_DIR/nvm.sh"
            set -u
        fi
    fi

    # Verify nvm is loaded
    if command -v nvm &>/dev/null; then
        log_ok "nvm is available"

        # Install/use Node LTS if node not found or if we want to ensure user-managed node
        if ! command -v node &>/dev/null || [[ $(which node) == "/usr/bin/node" ]]; then
            log_info "Installing/using Node.js LTS via nvm..."
            set +u
            nvm install --lts
            nvm use --lts
            set -u
        fi
    else
        log_error "Failed to install/load nvm"
        return 1
    fi

    # Configure npm to use user-writable prefix for global packages
    # This fixes EACCES errors when node is in a read-only location (e.g., /opt/nvm/)
    configure_npm_prefix
}

# Configure npm to use a user-writable prefix for global packages
# Prevents EACCES errors when system node is in read-only locations like /opt/nvm/
configure_npm_prefix() {
    local npm_prefix="$HOME/.npm-global"

    # Check if current npm global directory is writable
    local current_prefix
    current_prefix=$(npm config get prefix 2>/dev/null || echo "")

    if [ -n "$current_prefix" ] && [ -w "$current_prefix/lib" ] 2>/dev/null; then
        # Current prefix is writable, no change needed
        return 0
    fi

    # Create user-writable npm global directory
    if [ ! -d "$npm_prefix" ]; then
        log_info "Creating user-writable npm prefix at $npm_prefix"
        mkdir -p "$npm_prefix/bin" "$npm_prefix/lib"
    fi

    # Configure npm to use user prefix
    npm config set prefix "$npm_prefix" 2>/dev/null || true

    # Add to PATH for this session
    export PATH="$npm_prefix/bin:$PATH"

    # Persist PATH in .bashrc if not already present
    if ! grep -q 'npm-global/bin' "$HOME/.bashrc" 2>/dev/null; then
        echo 'export PATH="$HOME/.npm-global/bin:$PATH"' >> "$HOME/.bashrc"
        log_info "Added ~/.npm-global/bin to PATH in .bashrc"
    fi

    log_ok "npm configured to use user-writable prefix: $npm_prefix"
}
