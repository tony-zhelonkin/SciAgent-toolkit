#!/usr/bin/env bash
#
# Notebook Tools MCP Server Setup Script
#
# Installs the notebook-tools-mcp package for lightweight notebook
# read/search/navigate capabilities. No Jupyter server required.
#
# Requirements:
#   - python3
#   - pip

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

if [ -f "${SCRIPT_DIR}/../common.sh" ]; then
    source "${SCRIPT_DIR}/../common.sh"
else
    echo "Error: common.sh not found in ${SCRIPT_DIR}/.."
    exit 1
fi

# ============================================================================
# Parse arguments
# ============================================================================
PROJECT_DIR="${PWD}"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --project-dir)
            PROJECT_DIR="$2"
            shift 2
            ;;
        *)
            log_error "Unknown argument: $1"
            exit 1
            ;;
    esac
done

REPO_URL="https://github.com/tony-zhelonkin/notebook-tools-mcp.git"

log_info "Setting up Notebook Tools MCP Server..."

# ============================================================================
# 1. Check prerequisites
# ============================================================================
PYTHON_CMD="python"
if ! command -v python &>/dev/null; then
    if command -v python3 &>/dev/null; then
        PYTHON_CMD="python3"
    else
        log_error "python not found. Please install Python 3."
        exit 1
    fi
fi

# ============================================================================
# 2. Install package
# ============================================================================
if ${PYTHON_CMD} -c "import notebook_tools_mcp" 2>/dev/null; then
    log_ok "notebook-tools-mcp is already installed"
else
    log_info "Installing notebook-tools-mcp from ${REPO_URL}..."

    # Prefer local clone if present (submodule), otherwise install from git
    LOCAL_DIR="${PROJECT_DIR}/01_modules/notebook-tools-mcp"
    if [ -d "${LOCAL_DIR}" ] && [ -f "${LOCAL_DIR}/pyproject.toml" ]; then
        log_info "Found local clone at ${LOCAL_DIR}, installing editable"
        pip install -e "${LOCAL_DIR}"
    else
        pip install "git+${REPO_URL}"
    fi
fi

# ============================================================================
# 3. Verify
# ============================================================================
if ${PYTHON_CMD} -c "import notebook_tools_mcp" 2>/dev/null; then
    log_ok "notebook-tools-mcp installed and verified"
else
    log_error "notebook-tools-mcp import failed after installation"
    exit 1
fi

log_info "Enable with: ./scripts/manage-addon.sh enable notebook-tools --project-dir ${PROJECT_DIR}"
