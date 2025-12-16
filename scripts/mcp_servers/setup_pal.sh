#!/usr/bin/env bash
#
# PAL MCP Server Setup Script
#
# Sets up the PAL MCP server for collaboration, planning, and code analysis.
# Works with Claude Code.
#
# Requirements:
#   - git
#   - python3
#   - pip

set -euo pipefail

# Get script directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TOOLKIT_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
INSTALL_DIR="${TOOLKIT_ROOT}/mcp_servers/pal"

# Source common utilities
if [ -f "${SCRIPT_DIR}/../common.sh" ]; then
    source "${SCRIPT_DIR}/../common.sh"
else
    echo "Error: common.sh not found in ${SCRIPT_DIR}/.."
    exit 1
fi

log_info "Setting up PAL MCP Server..."

# ============================================================================
# 1. Check prerequisites
# ============================================================================
log_info "Checking prerequisites..."

# Check for git
if ! command -v git &>/dev/null; then
    log_error "git not found. Please install git manually."
    exit 1
fi

# Check for python3
if ! command -v python3 &>/dev/null; then
    log_error "python3 not found. Please install python3 manually."
    exit 1
fi

# Check for pip3
if ! command -v pip3 &>/dev/null; then
    log_warn "pip3 not found. Trying python3 -m pip..."
    if ! python3 -m pip --version &>/dev/null; then
        log_error "pip not found. Please install python3-pip."
        exit 1
    fi
fi

# Check for uv (optional but recommended for faster installs)
if ! command -v uv &>/dev/null; then
    log_info "uv not found. Installing uv for faster setup..."
    
    # Install uv user-local
    if python3 -m pip install --user uv &>/dev/null; then
        log_ok "uv installed via pip --user"
        # Add to PATH for this session if needed (heuristic for linux/mac)
        export PATH="$HOME/.local/bin:$PATH"
    else
        log_warn "Failed to install uv. Will fall back to standard pip."
    fi
fi

# ============================================================================
# 2. Install PAL MCP server
# ============================================================================
log_info "Installing PAL MCP server locally to ${INSTALL_DIR}..."

mkdir -p "${INSTALL_DIR}"
cd "${INSTALL_DIR}"

# Create venv if not exists
if [ ! -d "venv" ]; then
    if command -v uv &>/dev/null; then
        uv venv venv
    else
        python3 -m venv venv
    fi
    log_ok "Created virtual environment in ${INSTALL_DIR}/venv"
fi

# Install dependencies and package
log_info "Installing pal-mcp-server package..."

PIP_CMD="${INSTALL_DIR}/venv/bin/python -m pip"
if command -v uv &>/dev/null; then
    PIP_CMD="uv pip"
    export VIRTUAL_ENV="${INSTALL_DIR}/venv"
fi

# Install from git
# We use the HTTPS URL to avoid SSH key requirements
REPO_URL="git+https://github.com/BeehiveInnovations/pal-mcp-server.git"

if $PIP_CMD install --upgrade "${REPO_URL}"; then
    log_ok "Installed pal-mcp-server"
else
    log_error "Failed to install pal-mcp-server"
    exit 1
fi

# ============================================================================
# 3. Create Wrapper Script
# ============================================================================
WRAPPER_PATH="${INSTALL_DIR}/start-pal.py"
log_info "Creating startup wrapper at ${WRAPPER_PATH}..."

cat > "${WRAPPER_PATH}" << 'EOF'
#!/usr/bin/env python3
import os
import sys
import subprocess
from pathlib import Path

def main():
    # 1. Determine Project Root
    # We assume the current working directory is the project root when this script is called by MCP.
    project_root = Path.cwd()
    
    # 2. Load .env variables
    # We look for .env in the project root.
    env_files = [project_root / ".env", project_root / ".devcontainer" / ".env"]
    
    # Current environment
    env = os.environ.copy()
    
    for env_file in env_files:
        if env_file.exists():
            try:
                with open(env_file, 'r') as f:
                    for line in f:
                        line = line.strip()
                        if not line or line.startswith('#'): continue
                        if '=' in line:
                            k, v = line.split('=', 1)
                            # Remove surrounding quotes
                            v = v.strip()
                            if len(v) >= 2 and ((v[0] == '"' and v[-1] == '"') or (v[0] == "'" and v[-1] == "'")):
                                v = v[1:-1]
                            # Only set if not already in env (let shell override)? 
                            # actually for secrets, the .env usually holds the truth if shell is empty.
                            if k not in env:
                                env[k] = v
            except Exception as e:
                # Silently ignore errors reading .env, just proceed
                pass

    # 3. Path to the server executable in the venv
    # This script is located at .../mcp_servers/pal/start-pal.py
    # The venv is at .../mcp_servers/pal/venv
    script_dir = Path(__file__).parent.absolute()
    
    # We want to run the python interpreter from the venv, calling the module
    venv_python = script_dir / "venv" / "bin" / "python"
    
    if not venv_python.exists():
        sys.stderr.write(f"Error: Python interpreter not found at {venv_python}\n")
        sys.exit(1)

    # 4. Construct command
    # python -m pal_mcp_server [args...]
    cmd = [str(venv_python), "-m", "pal_mcp_server"] + sys.argv[1:]
    
    # 5. Execute
    try:
        # Replace current process
        os.execve(str(venv_python), cmd, env)
    except OSError as e:
        sys.stderr.write(f"Error executing PAL server: {e}\n")
        sys.exit(1)

if __name__ == "__main__":
    main()
EOF

chmod +x "${WRAPPER_PATH}"
log_ok "Created wrapper script"

# ============================================================================
# 4. Verify
# ============================================================================
log_info "Verifying installation..."

if timeout 10 "${WRAPPER_PATH}" --help &>/dev/null; then
    log_ok "PAL server verification successful"
else
    # It might fail if it strictly requires API keys even for help, but usually help is free.
    # We'll warn but not fail the script, as the keys might be provided later by the user.
    log_warn "PAL server installed, but --help check didn't exit cleanly."
    log_info "This might be normal if API keys are missing. The wrapper will load them from .env at runtime."
fi

log_ok "PAL MCP Server setup complete"