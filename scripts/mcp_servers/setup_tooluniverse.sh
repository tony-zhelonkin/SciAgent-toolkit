#!/usr/bin/env bash
#
# ToolUniverse MCP Server Setup Script
#
# Sets up ToolUniverse MCP server with 600+ scientific tools.
# Works with both Claude Code and Codex CLI.
#
# Requirements:
#   - uv package manager
#   - Python 3.10+
#
# Optional environment variables:
#   - AZURE_OPENAI_API_KEY: For summarization hooks
#   - AZURE_OPENAI_ENDPOINT: Azure OpenAI endpoint URL

set -euo pipefail

# Get script directory and project directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"

# Source common utilities
if [ -f "${SCRIPT_DIR}/../common.sh" ]; then
    source "${SCRIPT_DIR}/../common.sh"
else
    echo "Error: common.sh not found in ${SCRIPT_DIR}/.."
    exit 1
fi



# ToolUniverse installation directory
TOOLUNIVERSE_ENV="${PROJECT_DIR}/tooluniverse-env"

log_info "Setting up ToolUniverse MCP Server..."
log_info "Installation directory: ${TOOLUNIVERSE_ENV}"

# ============================================================================
# 1. Check prerequisites
# ============================================================================
log_info "Checking prerequisites..."

# Check for uv
if ! command -v uv &>/dev/null; then
    log_error "uv package manager not found"
    log_info "Installing uv..."

    if curl -LsSf https://astral.sh/uv/install.sh | sh; then
        log_ok "uv installed successfully"

        # Add to PATH for current session (UV installs to ~/.local/bin)
        export PATH="$HOME/.local/bin:$PATH"

        # Verify installation
        if ! command -v uv &>/dev/null; then
            log_error "uv installation failed - command not found after install"
            exit 1
        fi
    else
        log_error "Failed to install uv"
        log_info "Please install manually: https://github.com/astral-sh/uv"
        exit 1
    fi
else
    UV_VERSION=$(uv --version 2>/dev/null || echo "unknown")
    log_ok "uv already installed (${UV_VERSION})"
fi

# Check Python version
if command -v python3 &>/dev/null; then
    PYTHON_VERSION=$(python3 --version 2>&1 | awk '{print $2}')
    log_ok "Python ${PYTHON_VERSION} found"
else
    log_error "Python 3 not found. Please install Python 3.10 or later."
    exit 1
fi

# ============================================================================
# 2. Install ToolUniverse
# ============================================================================
log_info "Installing ToolUniverse..."

# Create working directory if it doesn't exist
mkdir -p "${TOOLUNIVERSE_ENV}"

# Create virtual environment first
log_info "Creating virtual environment..."
if uv venv "${TOOLUNIVERSE_ENV}" &>/dev/null; then
    log_ok "Virtual environment created"
else
    log_error "Failed to create virtual environment"
    exit 1
fi

# Install ToolUniverse using uv
log_info "Installing ToolUniverse package (this may take a few minutes)..."
if uv --directory "${TOOLUNIVERSE_ENV}" pip install tooluniverse &>/dev/null; then
    log_ok "ToolUniverse installed successfully"
else
    log_error "Failed to install ToolUniverse"
    exit 1
fi

# Verify installation
if uv --directory "${TOOLUNIVERSE_ENV}" run python -c "import tooluniverse; print('ToolUniverse installed successfully')" &>/dev/null; then
    log_ok "ToolUniverse installation verified"
else
    log_error "ToolUniverse installation verification failed"
    exit 1
fi

# ============================================================================
# 3. Detect correct MCP command name
# ============================================================================
log_info "Detecting ToolUniverse MCP command..."

# Try different command names to find which one works
TOOLUNIVERSE_CMD=""
if uv --directory "${TOOLUNIVERSE_ENV}" run tooluniverse-mcp --help &>/dev/null; then
    TOOLUNIVERSE_CMD="tooluniverse-mcp"
    log_ok "Using command: tooluniverse-mcp"
elif uv --directory "${TOOLUNIVERSE_ENV}" run python -m tooluniverse.smcp_server --help &>/dev/null; then
    TOOLUNIVERSE_CMD="python -m tooluniverse.smcp_server"
    log_ok "Using command: python -m tooluniverse.smcp_server"
else
    log_error "Could not detect ToolUniverse MCP command"
    log_error "Tried: tooluniverse-mcp, python -m tooluniverse.smcp_server"
    exit 1
fi

# ============================================================================
# 4. Test MCP server
# ============================================================================
log_info "Testing ToolUniverse MCP server..."

if timeout 10 uv --directory "${TOOLUNIVERSE_ENV}" run ${TOOLUNIVERSE_CMD} --help &>/dev/null; then
    log_ok "ToolUniverse MCP server can start"
else
    log_warn "Could not verify MCP server startup (may need first-time initialization)"
fi

# ============================================================================
# 5. Configure for Claude Code
# ============================================================================
if command -v claude &>/dev/null; then
    log_info "Configuring ToolUniverse for Claude Code..."

    # Determine configuration options
    CLAUDE_MCP_ARGS="--scope local"

    # Check for Azure OpenAI credentials (optional)
    # Note: --exclude-tool-types removed in ToolUniverse 1.0.14+
    # Use --compact-mode for minimal tool set if needed
    if [ -n "${AZURE_OPENAI_API_KEY:-}" ] && [ -n "${AZURE_OPENAI_ENDPOINT:-}" ]; then
        log_info "Azure OpenAI credentials detected - enabling SummarizationHook"
        CLAUDE_MCP_ARGS="${CLAUDE_MCP_ARGS} --env AZURE_OPENAI_API_KEY=${AZURE_OPENAI_API_KEY}"
        CLAUDE_MCP_ARGS="${CLAUDE_MCP_ARGS} --env AZURE_OPENAI_ENDPOINT=${AZURE_OPENAI_ENDPOINT}"
        TOOLUNIVERSE_CMD_ARGS="--hook-type SummarizationHook"
    else
        log_warn "Azure OpenAI credentials not set - SummarizationHook disabled"
        log_info "Set AZURE_OPENAI_API_KEY and AZURE_OPENAI_ENDPOINT to enable"
        TOOLUNIVERSE_CMD_ARGS=""
    fi

    # Add MCP server using claude mcp add command
    log_info "Adding ToolUniverse to Claude Code MCP servers..."

    # Note: This is the recommended approach from the documentation
    log_warn "Claude Code MCP server must be added manually or via CLI"
    log_info ""
    log_info "To add ToolUniverse to Claude Code, run:"
    log_info "  claude mcp add tooluniverse ${CLAUDE_MCP_ARGS} -- \\"
    log_info "    uv --directory ${TOOLUNIVERSE_ENV} run ${TOOLUNIVERSE_CMD} ${TOOLUNIVERSE_CMD_ARGS}"
    log_info ""
    log_info "To verify: claude mcp list"
    log_info ""
else
    log_warn "Claude Code not installed - skipping Claude-specific configuration"
fi

# ============================================================================
# 6. Configure for Codex CLI
# ============================================================================
if command -v codex &>/dev/null; then
    log_info "Configuring ToolUniverse for Codex CLI..."

    # Create Codex config directory if needed
    mkdir -p "$HOME/.codex"

    CODEX_CONFIG="$HOME/.codex/config.toml"

    # Check if config exists
    if [ -f "$CODEX_CONFIG" ]; then
        log_warn "Codex config already exists: ${CODEX_CONFIG}"
        log_info "Append the following to your config manually:"
    else
        log_info "Creating Codex config: ${CODEX_CONFIG}"
    fi

    # Generate config snippet
    cat << EOF

# ToolUniverse MCP Server Configuration
# Add this to ${CODEX_CONFIG}

[mcp_servers.tooluniverse]
command = "uv"
args = [
  "--directory",
  "${TOOLUNIVERSE_ENV}",
  "run",
  "${TOOLUNIVERSE_CMD}"
]
startup_timeout_sec = 60

# Optional: Research-focused ToolUniverse instance with specific tools
# [mcp_servers.tooluniverse-research]
# command = "uv"
# args = [
#   "--directory",
#   "${TOOLUNIVERSE_ENV}",
#   "run",
#   "${TOOLUNIVERSE_CMD}",
#   "--include-tools",
#   "EuropePMC_search_articles,ChEMBL_search_similar_molecules,search_clinical_trials"
# ]
# startup_timeout_sec = 60

EOF

    # Optionally write config if it doesn't exist
    if [ ! -f "$CODEX_CONFIG" ]; then
        cat > "$CODEX_CONFIG" << EOF
# Codex CLI Configuration
# Generated by setup_tooluniverse.sh

[mcp_servers.tooluniverse]
command = "uv"
args = [
  "--directory",
  "${TOOLUNIVERSE_ENV}",
  "run",
  "${TOOLUNIVERSE_CMD}"
]
startup_timeout_sec = 60
EOF
        log_ok "Created Codex config with ToolUniverse"
    fi
else
    log_warn "Codex CLI not installed - skipping Codex-specific configuration"
fi

# ============================================================================
# 7. Create helper scripts
# ============================================================================
log_info "Creating helper scripts..."

# Create a quick test script
cat > "${PROJECT_DIR}/test_tooluniverse.sh" << EOF
#!/usr/bin/env bash
# Quick test of ToolUniverse MCP server
SCRIPT_DIR="\$(cd "\$(dirname "\${BASH_SOURCE[0]}")" && pwd)"
TOOLUNIVERSE_ENV="\${SCRIPT_DIR}/tooluniverse-env"

echo "Testing ToolUniverse MCP server..."
uv --directory "\${TOOLUNIVERSE_ENV}" run ${TOOLUNIVERSE_CMD} --help
EOF

chmod +x "${PROJECT_DIR}/test_tooluniverse.sh"
log_ok "Created test script: ${PROJECT_DIR}/test_tooluniverse.sh"

# ============================================================================
# Summary
# ============================================================================
log_ok "ToolUniverse MCP Server setup complete"
log_info ""
log_info "Installation details:"
log_info "  Location: ${TOOLUNIVERSE_ENV}"
log_info "  Command: uv --directory ${TOOLUNIVERSE_ENV} run ${TOOLUNIVERSE_CMD}"
log_info ""
log_info "Features available:"
log_info "  - 600+ scientific tools across multiple domains"
log_info "  - Drug discovery and development"
log_info "  - Genomics and molecular biology"
log_info "  - Literature research (PubMed, Semantic Scholar)"
log_info "  - Clinical research (ClinicalTrials.gov, FDA)"
log_info ""
log_info "Next steps:"
log_info "  For Claude Code:"
log_info "    Run the command shown above to add ToolUniverse"
log_info "    Then: claude mcp list"
log_info ""
log_info "  For Codex CLI:"
if [ -n "${CODEX_CONFIG:-}" ]; then
    log_info "    Config file: ${CODEX_CONFIG}"
else
    log_info "    Config file: ~/.codex/config.toml (not created)"
fi
log_info "    Run: codex"
log_info "    In chat: /mcp"
log_info ""
log_info "Optional configuration:"
log_info "  Set AZURE_OPENAI_API_KEY and AZURE_OPENAI_ENDPOINT"
log_info "  to enable SummarizationHook for long outputs"
