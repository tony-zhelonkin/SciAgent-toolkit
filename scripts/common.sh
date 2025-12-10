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
}
