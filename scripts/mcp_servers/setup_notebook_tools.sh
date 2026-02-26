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
LOCAL_DIR="${PROJECT_DIR}/01_modules/notebook-tools-mcp"
NEEDS_INSTALL=false

if ${PYTHON_CMD} -c "import notebook_tools_mcp" 2>/dev/null; then
    # Check if local submodule version is newer than installed version
    if [ -d "${LOCAL_DIR}" ] && [ -f "${LOCAL_DIR}/pyproject.toml" ]; then
        INSTALLED_VER=$(${PYTHON_CMD} -c "import notebook_tools_mcp; print(notebook_tools_mcp.__version__)" 2>/dev/null || echo "0.0.0")
        SOURCE_VER=$(${PYTHON_CMD} -c "
import configparser, pathlib
p = pathlib.Path('${LOCAL_DIR}/pyproject.toml')
# Simple TOML version extraction (avoids tomllib dep on py3.10)
for line in p.read_text().splitlines():
    if line.strip().startswith('version'):
        print(line.split('=')[1].strip().strip('\"'))
        break
" 2>/dev/null || echo "0.0.0")
        if [ "${INSTALLED_VER}" != "${SOURCE_VER}" ]; then
            log_info "Version mismatch: installed=${INSTALLED_VER}, source=${SOURCE_VER}"
            NEEDS_INSTALL=true
        else
            log_ok "notebook-tools-mcp ${INSTALLED_VER} is already installed and up to date"
        fi
    else
        log_ok "notebook-tools-mcp is already installed"
    fi
else
    NEEDS_INSTALL=true
fi

if [ "${NEEDS_INSTALL}" = true ]; then
    log_info "Installing notebook-tools-mcp..."

    # Prefer local clone if present (submodule), otherwise install from git
    if [ -d "${LOCAL_DIR}" ] && [ -f "${LOCAL_DIR}/pyproject.toml" ]; then
        log_info "Found local clone at ${LOCAL_DIR}, installing editable"
        pip install -e "${LOCAL_DIR}"
    else
        log_info "Installing from ${REPO_URL}"
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
