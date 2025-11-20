#!/usr/bin/env bash
#
# MCP Infrastructure Setup - Main Orchestrator
#
# This is the main entry point for setting up Claude Code, Codex CLI, and all MCP servers.
# It orchestrates the installation of individual components in the correct order.
#
# Usage:
#   ./setup_mcp_infrastructure.sh [OPTIONS]
#
# Options:
#   --claude-only          Install only Claude Code
#   --codex-only           Install only Codex CLI
#   --mcp-only             Install only MCP servers (assumes Claude/Codex already installed)
#   --skip-claude          Skip Claude Code installation
#   --skip-codex           Skip Codex CLI installation
#   --skip-pubmed          Skip PubMed MCP server
#   --skip-tooluniverse    Skip ToolUniverse MCP server
#   --help                 Show this help message
#
# Example:
#   # Full installation
#   ./setup_mcp_infrastructure.sh
#
#   # Install Claude Code and MCP servers only (skip Codex)
#   ./setup_mcp_infrastructure.sh --skip-codex
#
#   # Install only ToolUniverse MCP server
#   ./setup_mcp_infrastructure.sh --mcp-only --skip-pubmed

set -euo pipefail

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
CYAN='\033[0;36m'
MAGENTA='\033[0;35m'
NC='\033[0m' # No Color

# Helper functions
log_info()  { echo -e "${BLUE}[INFO]${NC} $*"; }
log_ok()    { echo -e "${GREEN}[OK]${NC} $*"; }
log_warn()  { echo -e "${YELLOW}[WARN]${NC} $*"; }
log_error() { echo -e "${RED}[ERROR]${NC} $*"; }
log_step()  { echo -e "${CYAN}[STEP]${NC} $*"; }
separator() { echo -e "${MAGENTA}======================================${NC}"; echo -e "${MAGENTA}==== $* ${NC}"; echo -e "${MAGENTA}======================================${NC}"; }

# Default flags
INSTALL_CLAUDE=true
INSTALL_CODEX=true
INSTALL_MCP=true
INSTALL_PUBMED=true
INSTALL_TOOLUNIVERSE=true

# Parse command-line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --claude-only)
            INSTALL_CODEX=false
            INSTALL_MCP=false
            shift
            ;;
        --codex-only)
            INSTALL_CLAUDE=false
            INSTALL_MCP=false
            shift
            ;;
        --mcp-only)
            INSTALL_CLAUDE=false
            INSTALL_CODEX=false
            shift
            ;;
        --skip-claude)
            INSTALL_CLAUDE=false
            shift
            ;;
        --skip-codex)
            INSTALL_CODEX=false
            shift
            ;;
        --skip-pubmed)
            INSTALL_PUBMED=false
            shift
            ;;
        --skip-tooluniverse)
            INSTALL_TOOLUNIVERSE=false
            shift
            ;;
        --help)
            grep '^#' "$0" | grep -v '#!/usr/bin/env' | sed 's/^# //' | sed 's/^#//'
            exit 0
            ;;
        *)
            log_error "Unknown option: $1"
            echo "Use --help for usage information"
            exit 1
            ;;
    esac
done

# Get script directory and project root
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"

# Detect OS
OS_TYPE="$(uname -s)"
case "$OS_TYPE" in
    Linux*)     PLATFORM=Linux;;
    Darwin*)    PLATFORM=macOS;;
    *)          PLATFORM=Unknown;;
esac

log_info "Platform: $PLATFORM"
log_info "Project directory: $PROJECT_DIR"

separator "MCP Infrastructure Setup"

# ============================================================================
# Installation Summary
# ============================================================================
log_info "Installation plan:"
echo "  Claude Code:      $([ "$INSTALL_CLAUDE" = true ] && echo "YES" || echo "SKIP")"
echo "  Codex CLI:        $([ "$INSTALL_CODEX" = true ] && echo "YES" || echo "SKIP")"
echo "  MCP Servers:"
echo "    - serena:       $([ "$INSTALL_MCP" = true ] && echo "YES" || echo "SKIP")"
echo "    - sequential:   $([ "$INSTALL_MCP" = true ] && echo "YES" || echo "SKIP")"
echo "    - PubMed:       $([ "$INSTALL_PUBMED" = true ] && echo "YES" || echo "SKIP")"
echo "    - ToolUniverse: $([ "$INSTALL_TOOLUNIVERSE" = true ] && echo "YES" || echo "SKIP")"
echo ""

# ============================================================================
# 1. Install Claude Code
# ============================================================================
if [ "$INSTALL_CLAUDE" = true ]; then
    separator "Installing Claude Code"
    if [ -f "${SCRIPT_DIR}/install_claude.sh" ]; then
        bash "${SCRIPT_DIR}/install_claude.sh" || {
            log_error "Claude Code installation failed"
            exit 1
        }
    else
        log_error "install_claude.sh not found in ${SCRIPT_DIR}"
        exit 1
    fi
else
    log_warn "Skipping Claude Code installation"
fi

# ============================================================================
# 2. Install Codex CLI
# ============================================================================
if [ "$INSTALL_CODEX" = true ]; then
    separator "Installing Codex CLI"
    if [ -f "${SCRIPT_DIR}/install_codex.sh" ]; then
        bash "${SCRIPT_DIR}/install_codex.sh" || {
            log_error "Codex CLI installation failed"
            exit 1
        }
    else
        log_error "install_codex.sh not found in ${SCRIPT_DIR}"
        exit 1
    fi
else
    log_warn "Skipping Codex CLI installation"
fi

# ============================================================================
# 3. Install MCP Servers
# ============================================================================
if [ "$INSTALL_MCP" = true ]; then
    separator "Installing MCP Servers"

    # Serena MCP
    log_step "Setting up Serena MCP server..."
    if [ -f "${SCRIPT_DIR}/mcp_servers/setup_serena.sh" ]; then
        bash "${SCRIPT_DIR}/mcp_servers/setup_serena.sh" || {
            log_error "Serena MCP setup failed"
            exit 1
        }
    else
        log_error "mcp_servers/setup_serena.sh not found"
        exit 1
    fi

    # Sequential Thinking MCP
    log_step "Setting up Sequential Thinking MCP server..."
    if [ -f "${SCRIPT_DIR}/mcp_servers/setup_sequential_thinking.sh" ]; then
        bash "${SCRIPT_DIR}/mcp_servers/setup_sequential_thinking.sh" || {
            log_error "Sequential Thinking MCP setup failed"
            exit 1
        }
    else
        log_error "mcp_servers/setup_sequential_thinking.sh not found"
        exit 1
    fi
else
    log_warn "Skipping base MCP servers (serena, sequential-thinking)"
fi

# ============================================================================
# 4. Install PubMed MCP
# ============================================================================
if [ "$INSTALL_PUBMED" = true ]; then
    separator "Installing PubMed MCP Server"
    if [ -f "${SCRIPT_DIR}/mcp_servers/setup_pubmed.sh" ]; then
        bash "${SCRIPT_DIR}/mcp_servers/setup_pubmed.sh" || {
            log_error "PubMed MCP setup failed"
            exit 1
        }
    else
        log_error "mcp_servers/setup_pubmed.sh not found"
        exit 1
    fi
else
    log_warn "Skipping PubMed MCP server"
fi

# ============================================================================
# 5. Install ToolUniverse MCP
# ============================================================================
if [ "$INSTALL_TOOLUNIVERSE" = true ]; then
    separator "Installing ToolUniverse MCP Server"
    if [ -f "${SCRIPT_DIR}/mcp_servers/setup_tooluniverse.sh" ]; then
        bash "${SCRIPT_DIR}/mcp_servers/setup_tooluniverse.sh" || {
            log_error "ToolUniverse MCP setup failed"
            exit 1
        }
    else
        log_error "mcp_servers/setup_tooluniverse.sh not found"
        exit 1
    fi
else
    log_warn "Skipping ToolUniverse MCP server"
fi

# ============================================================================
# 6. Final Configuration
# ============================================================================
separator "Finalizing Configuration"

# Merge all MCP configurations
log_info "Merging MCP configurations..."
# Note: Config merger not currently implemented
log_info "Using individual MCP server configurations"

# ============================================================================
# 7. Verification
# ============================================================================
separator "Verification"

CHECKS_PASSED=0
CHECKS_FAILED=0

# Check Claude Code
if [ "$INSTALL_CLAUDE" = true ]; then
    printf "Claude Code: "
    if command -v claude &>/dev/null && claude --version &>/dev/null; then
        CLAUDE_VER=$(claude --version 2>/dev/null | head -1 || echo "unknown")
        log_ok "${CLAUDE_VER}"
        CHECKS_PASSED=$((CHECKS_PASSED + 1))
    else
        log_error "NOT INSTALLED"
        CHECKS_FAILED=$((CHECKS_FAILED + 1))
    fi
fi

# Check Codex CLI
if [ "$INSTALL_CODEX" = true ]; then
    printf "Codex CLI: "
    if command -v codex &>/dev/null && codex --version &>/dev/null; then
        CODEX_VER=$(codex --version 2>/dev/null | head -1 || echo "unknown")
        log_ok "${CODEX_VER}"
        CHECKS_PASSED=$((CHECKS_PASSED + 1))
    else
        log_error "NOT INSTALLED"
        CHECKS_FAILED=$((CHECKS_FAILED + 1))
    fi
fi

# Check MCP configs
printf "MCP Configuration: "
if [ -f "${PROJECT_DIR}/.mcp.json" ]; then
    if python3 -m json.tool "${PROJECT_DIR}/.mcp.json" &>/dev/null; then
        log_ok "Valid"
        CHECKS_PASSED=$((CHECKS_PASSED + 1))
    else
        log_error "Invalid JSON"
        CHECKS_FAILED=$((CHECKS_FAILED + 1))
    fi
else
    log_warn "Not found"
fi

# ============================================================================
# Summary
# ============================================================================
separator "Setup Complete"

if [ $CHECKS_FAILED -eq 0 ]; then
    log_ok "All checks passed (${CHECKS_PASSED})"
    echo ""
    log_info "Next steps:"
    echo ""

    if [ "$INSTALL_CLAUDE" = true ]; then
        echo "  For Claude Code:"
        echo "    1. Run 'claude' to start Claude Code"
        echo "    2. Use /mcp to verify MCP servers are loaded"
        echo ""
    fi

    if [ "$INSTALL_CODEX" = true ]; then
        echo "  For Codex CLI:"
        echo "    1. Run 'codex' to start Codex CLI"
        echo "    2. Sign in with ChatGPT or use OPENAI_API_KEY"
        echo "    3. Use /mcp to verify MCP servers are loaded"
        echo ""
    fi

    echo "  Available MCP Servers:"
    [ "$INSTALL_MCP" = true ] && echo "    - serena: Semantic code search and editing"
    [ "$INSTALL_MCP" = true ] && echo "    - sequential-thinking: Structured reasoning"
    [ "$INSTALL_PUBMED" = true ] && echo "    - pubmed: Biomedical literature access"
    [ "$INSTALL_TOOLUNIVERSE" = true ] && echo "    - tooluniverse: 600+ scientific tools"
    echo ""

    log_info "Configuration files:"
    echo "    - Claude Code: ${PROJECT_DIR}/.mcp.json"
    [ "$INSTALL_CODEX" = true ] && echo "    - Codex CLI: ~/.codex/config.toml"
    echo ""
else
    log_error "Some checks failed (${CHECKS_PASSED} passed, ${CHECKS_FAILED} failed)"
    log_info "Review error messages above for troubleshooting"
    exit 1
fi
