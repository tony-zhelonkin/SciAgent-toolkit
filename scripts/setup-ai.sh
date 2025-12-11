#!/usr/bin/env bash
#
# setup-ai.sh - Runtime AI tools installer for scbio-docker
#
# Run this script inside the container to set up AI tools (Claude Code, MCP servers).
# First run takes 5-15 minutes due to Serena compilation. Subsequent runs are fast.
#
# Usage:
#   ./setup-ai.sh [OPTIONS]
#
# Options:
#   --minimal         Fast setup: Skip Serena and Codex (2-3 min instead of 15 min)
#   --skip-serena     Skip Serena MCP server (saves 5-15 min on first run)
#   --skip-codex      Skip Codex CLI installation
#   --skip-gemini     Skip Gemini CLI installation
#   --force           Force reinstall even if already configured
#   --help            Show this help message
#
# Creates:
#   - .mcp.json in current directory (Claude Code MCP configuration)
#   - tooluniverse-env/ in current directory (ToolUniverse Python environment)
#   - ~/.local/bin/claude (Claude Code CLI, if not present)
#
# Example:
#   # Full setup (5-15 minutes first time)
#   ./setup-ai.sh
#
#   # Fast setup (2-3 minutes)
#   ./setup-ai.sh --minimal
#

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
log_error() { echo -e "${RED}[ERROR]${NC} $*" >&2; }

# Default settings
PROJECT_DIR="$(pwd)"
# Determine script directory to find sibling scripts
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCIAGENT_SCRIPTS="${SCRIPT_DIR}"
FORCE=false
MINIMAL=false
SKIP_SERENA=false
SKIP_CODEX=false
SKIP_GEMINI=false

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --minimal)
            MINIMAL=true
            SKIP_SERENA=true
            SKIP_CODEX=true
            SKIP_GEMINI=true
            shift
            ;;
        --skip-serena)
            SKIP_SERENA=true
            shift
            ;;
        --skip-codex)
            SKIP_CODEX=true
            shift
            ;;
        --skip-gemini)
            SKIP_GEMINI=true
            shift
            ;;
        --force)
            FORCE=true
            shift
            ;;
        --help)
            grep '^#' "$0" | grep -v '#!/usr/bin/env' | sed 's/^# //' | sed 's/^#//'
            exit 0
            ;;
        *)
            log_error "Unknown option: $1"
            log_info "Use --help for usage information"
            exit 1
            ;;
    esac
done

echo ""
echo -e "${BLUE}╔══════════════════════════════════════════╗${NC}"
echo -e "${BLUE}║       scbio-docker AI Setup              ║${NC}"
echo -e "${BLUE}╚══════════════════════════════════════════╝${NC}"
echo ""

# Check if SciAgent-toolkit is available
if [ ! -d "$SCIAGENT_SCRIPTS" ]; then
    log_error "SciAgent-toolkit scripts not found at $SCIAGENT_SCRIPTS"
    log_error "This script is part of SciAgent-toolkit and expects sibling scripts in the same directory."
    exit 1
fi

# Check if already configured
if [ -f "${PROJECT_DIR}/.mcp.json" ] && [ "$FORCE" = false ]; then
    log_ok "AI tools already configured (.mcp.json exists)"
    log_info "Run with --force to reinstall"

    # Quick status check
    if command -v claude &> /dev/null; then
        log_ok "Claude Code: $(claude --version 2>/dev/null || echo 'installed')"
    else
        log_warn "Claude Code not found in PATH"
    fi

    if [ -d "${PROJECT_DIR}/tooluniverse-env" ]; then
        log_ok "ToolUniverse: installed"
    fi

    echo ""
    log_info "To start Claude Code: claude"
    log_info "To check MCP servers: /mcp (inside claude)"
    exit 0
fi

# Build setup arguments
SETUP_ARGS=""
if [ "$SKIP_CODEX" = true ]; then
    SETUP_ARGS="$SETUP_ARGS --skip-codex"
fi
if [ "$SKIP_GEMINI" = true ]; then
    SETUP_ARGS="$SETUP_ARGS --skip-gemini"
fi

# Timing information
if [ "$MINIMAL" = true ]; then
    log_info "Running minimal setup (skipping Serena, Codex, Gemini)"
    log_info "Estimated time: 2-3 minutes"
else
    log_info "Running full setup"
    if [ "$SKIP_SERENA" = true ]; then
        log_info "Estimated time: 3-5 minutes (Serena skipped)"
    else
        log_info "Estimated time: 5-15 minutes (includes Serena compilation)"
    fi
fi
echo ""

# ... (after sourcing dependencies)

# Step 0: Environment Setup
log_info "Step 0/2: Checking environment configuration..."
if [ ! -f "${PROJECT_DIR}/.env" ]; then
    if [ -f "${SCIAGENT_SCRIPTS}/../templates/.env.template" ]; then
        log_info "Creating .env from template..."
        cp "${SCIAGENT_SCRIPTS}/../templates/.env.template" "${PROJECT_DIR}/.env"
        log_warn "Created .env file. PLEASE EDIT IT to add your API keys!"
    elif [ -f "${PROJECT_DIR}/.devcontainer/.env" ]; then
        log_info "Using existing .devcontainer/.env"
    else
        log_warn "No .env file or template found. You may need to create one manually for API keys."
    fi
else
    log_info "Found existing .env file."
fi

# Step 1: Run the main setup script
log_info "Step 1/2: Installing and configuring AI tools..."

# Ensure PATH includes potential install locations for re-runs
export PATH="$HOME/.local/bin:$HOME/.npm-global/bin:$PATH"

# The setup script needs to know where to create tooluniverse-env
# We run it from the project directory so paths are relative to project
if [ -f "${SCIAGENT_SCRIPTS}/setup_mcp_infrastructure.sh" ]; then
    bash "${SCIAGENT_SCRIPTS}/setup_mcp_infrastructure.sh" $SETUP_ARGS || {
        log_error "Setup failed. Check the output above for errors."
        exit 1
    }
else
    log_error "setup_mcp_infrastructure.sh not found"
    exit 1
fi

echo ""
echo -e "${GREEN}╔══════════════════════════════════════════╗${NC}"
echo -e "${GREEN}║       AI Setup Complete!                 ║${NC}"
echo -e "${GREEN}╚══════════════════════════════════════════╝${NC}"
echo ""

# Summary
log_info "Files created:"
if [ -f "${PROJECT_DIR}/.mcp.json" ]; then
    log_ok "  .mcp.json (MCP server configuration)"
fi
if [ -d "${PROJECT_DIR}/tooluniverse-env" ]; then
    log_ok "  tooluniverse-env/ (ToolUniverse Python environment)"
fi

echo ""
log_info "To start using AI tools:"
echo "  1. Run: claude"
echo "  2. Inside claude, check MCP servers: /mcp"
echo ""

if [ "$SKIP_SERENA" = true ]; then
    log_warn "Serena was skipped. To add it later, run:"
    echo "  uvx --from git+https://github.com/oraios/serena serena start-mcp-server --help"
    echo "  # Then re-run: ./setup-ai.sh --force"
fi

echo ""
