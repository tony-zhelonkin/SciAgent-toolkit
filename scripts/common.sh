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

# Configure npm to use user-writable prefix for global packages
# This is called automatically by ensure_npm_writable_prefix() only when needed
configure_npm_prefix() {
    # If npm is not available, do nothing
    if ! command -v npm &>/dev/null; then
        return 0
    fi

    local npm_global_dir="$HOME/.npm-global"

    # Create directory if it doesn't exist
    if [[ ! -d "$npm_global_dir" ]]; then
        mkdir -p "$npm_global_dir"
        mkdir -p "$npm_global_dir/bin"
        mkdir -p "$npm_global_dir/lib"
    fi

    # Set prefix
    npm config set prefix "$npm_global_dir"
    log_info "Configured npm global prefix to $npm_global_dir"

    # Add to PATH for current session
    if [[ ":$PATH:" != *":$npm_global_dir/bin:"* ]]; then
        export PATH="$npm_global_dir/bin:$PATH"
    fi

    # Add to PATH permanently in .bashrc if not already there
    if [ -f "$HOME/.bashrc" ] && ! grep -q 'export PATH="$HOME/.npm-global/bin:$PATH"' "$HOME/.bashrc" 2>/dev/null; then
        echo 'export PATH="$HOME/.npm-global/bin:$PATH"' >> "$HOME/.bashrc"
        log_info "Added ~/.npm-global/bin to PATH in .bashrc"
    fi
}

# Check if npm global prefix is writable; if not, configure user-local prefix
# This handles read-only nvm installations in containers without conflicting
# with writable nvm setups
ensure_npm_writable_prefix() {
    # If npm is not available, do nothing
    if ! command -v npm &>/dev/null; then
        return 0
    fi

    # Get current prefix
    local current_prefix
    current_prefix=$(npm config get prefix 2>/dev/null)

    # Check if the node_modules directory under prefix is writable
    local node_modules_dir="$current_prefix/lib/node_modules"

    # If node_modules doesn't exist, check if we can create it
    if [[ ! -d "$node_modules_dir" ]]; then
        # Try to check if parent is writable
        if [[ ! -w "$current_prefix/lib" ]] && [[ ! -w "$current_prefix" ]]; then
            log_info "npm prefix $current_prefix is read-only, switching to user-local prefix"
            configure_npm_prefix
            return 0
        fi
    elif [[ ! -w "$node_modules_dir" ]]; then
        log_info "npm node_modules ($node_modules_dir) is read-only, switching to user-local prefix"
        configure_npm_prefix
        return 0
    fi

    # Prefix is writable, no changes needed
    return 0
}

# Ensure ~/.npm-global/bin is in PATH if it exists and contains binaries
# Call this to pick up previously installed npm global packages
ensure_npm_global_path() {
    local npm_global_bin="$HOME/.npm-global/bin"

    # If directory exists and has content, ensure it's in PATH
    if [[ -d "$npm_global_bin" ]] && [[ -n "$(ls -A "$npm_global_bin" 2>/dev/null)" ]]; then
        if [[ ":$PATH:" != *":$npm_global_bin:"* ]]; then
            export PATH="$npm_global_bin:$PATH"
        fi

        # Also ensure it's in .bashrc for persistence
        if [ -f "$HOME/.bashrc" ] && ! grep -q 'export PATH="$HOME/.npm-global/bin:$PATH"' "$HOME/.bashrc" 2>/dev/null; then
            echo 'export PATH="$HOME/.npm-global/bin:$PATH"' >> "$HOME/.bashrc"
            log_info "Added ~/.npm-global/bin to PATH in .bashrc"
        fi
    fi
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
    # configure_npm_prefix
}

# ... (previous code)

# Ensure Python package is installed
ensure_pip_package() {
    local package="$1"
    if ! python3 -c "import ${package}" &>/dev/null; then
        log_info "Installing python package: ${package}..."
        if command -v pip3 &>/dev/null; then
            # Check if inside a virtualenv
            if [ -n "${VIRTUAL_ENV:-}" ] || python3 -c "import sys; sys.exit(0 if sys.prefix != sys.base_prefix else 1)" 2>/dev/null; then
                # Inside virtualenv: try direct install, fallback to --user on failure (e.g. read-only venv)
                if ! pip3 install "${package}"; then
                    log_warn "Standard install failed (read-only venv?). Retrying with --user..."
                    pip3 install --user "${package}" || return 1
                fi
            else
                # Outside virtualenv: try --user first
                pip3 install --user "${package}" || {
                    pip3 install "${package}" || return 1
                }
            fi
        elif command -v pip &>/dev/null; then
             if [ -n "${VIRTUAL_ENV:-}" ] || python3 -c "import sys; sys.exit(0 if sys.prefix != sys.base_prefix else 1)" 2>/dev/null; then
                if ! pip install "${package}"; then
                    log_warn "Standard install failed (read-only venv?). Retrying with --user..."
                    pip install --user "${package}" || return 1
                fi
            else
                pip install --user "${package}" || pip install "${package}" || return 1
            fi
        else
            log_error "pip not found. Cannot install ${package}."
            return 1
        fi
    fi
    return 0
}

# Securely load environment variables
load_env() {
    local project_dir="${1:-$PWD}"
    
    # Priority 1: .env in project root
    if [ -f "${project_dir}/.env" ]; then
        log_info "Loading secrets from ${project_dir}/.env"
        set -a
        source "${project_dir}/.env"
        set +a
    fi
    
    # Priority 2: .devcontainer/.env (often used in devcontainers)
    if [ -f "${project_dir}/.devcontainer/.env" ]; then
        log_info "Loading secrets from ${project_dir}/.devcontainer/.env"
        set -a
        source "${project_dir}/.devcontainer/.env"
        set +a
    fi
    
    # Warn if critical keys are missing
    if [ -z "${GEMINI_API_KEY:-}" ] && [ -z "${OPENAI_API_KEY:-}" ]; then
        log_warn "No API keys (GEMINI_API_KEY or OPENAI_API_KEY) found in environment."
        log_info "You may need to create a .env file with your API keys."
    fi
}

