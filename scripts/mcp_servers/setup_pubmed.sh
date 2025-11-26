#!/usr/bin/env bash
#
# PubMed MCP Server Setup Script
#
# Sets up the PubMed MCP server for Claude Code.
# PubMed provides access to 36+ million biomedical research articles.
#
# For Claude Code users:
#   Uses the plugin marketplace: /plugin marketplace add anthropics/life-sciences
#
# Note: PubMed integration is provided by Anthropic and does not require
# additional dependencies beyond Claude Code itself.

set -euo pipefail

# Get script directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Source common utilities
if [ -f "${SCRIPT_DIR}/../common.sh" ]; then
    source "${SCRIPT_DIR}/../common.sh"
else
    echo "Error: common.sh not found in ${SCRIPT_DIR}/.."
    exit 1
fi

log_info "Setting up PubMed MCP Server..."

# Check if Claude Code is installed
if ! command -v claude &>/dev/null; then
    log_error "Claude Code is not installed. Please install it first."
    exit 1
fi

# For Claude Code, PubMed is installed via the plugin marketplace
log_info "PubMed for Claude Code uses the plugin marketplace system"
log_info "To install PubMed in Claude Code:"
log_info "  1. Run 'claude' to start Claude Code"
log_info "  2. Use command: /plugin marketplace add anthropics/life-sciences"
log_info "  3. Use command: /plugin install pubmed@life-sciences"
log_info "  4. Restart Claude Code"
log_info "  5. Verify with: /mcp"
log_info ""

# Check if user wants to auto-configure (requires interactive session)
log_warn "PubMed MCP requires manual configuration via Claude Code CLI"
log_info "After installing, PubMed will be available as a connector"
log_info ""

# Create a marker file to track that PubMed setup instructions were shown
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
mkdir -p "${PROJECT_DIR}/.mcp_setup"
touch "${PROJECT_DIR}/.mcp_setup/pubmed_instructions_shown"

log_ok "PubMed setup instructions provided"
log_info ""
log_info "Features available with PubMed:"
log_info "  - Search 36+ million biomedical articles"
log_info "  - Access full-text articles from PubMed Central (PMC)"
log_info "  - Get article metadata (authors, abstracts, citations)"
log_info "  - Find related articles across NCBI databases"
log_info "  - Convert between ID formats (PMID, PMC ID, DOI)"
log_info ""
log_info "Example queries after setup:"
log_info "  - 'Find recent studies about immunotherapy for melanoma'"
log_info "  - 'What are the most cited papers on CRISPR gene editing?'"
log_info "  - 'Get the full text of PMID:12345678'"
