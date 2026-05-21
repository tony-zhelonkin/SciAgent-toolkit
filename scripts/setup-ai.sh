#!/usr/bin/env bash
#
# setup-ai.sh - AI context and role setup for scbio-docker projects
#
# Run this script inside the container to set up AI context files and activate roles.
#
# Usage:
#   ./setup-ai.sh [OPTIONS]
#
# Options:
#   --force           Force reinstall even if already configured
#   --help            Show this help message
#
# Creates:
#   - ~/.local/bin/claude (Claude Code CLI, if not present)
#   - CLAUDE.md, GEMINI.md, AGENTS.md, context.md (AI context files)
#   - .claude/agents/, .claude/skills/ (populated by role activation)
#   - 02_analysis/config/analysis_config.yaml (project parameters)
#
# Example:
#   ./setup-ai.sh
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

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
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
if [ -f "${PROJECT_DIR}/CLAUDE.md" ] && [ "$FORCE" = false ]; then
    log_ok "AI tools already configured (CLAUDE.md exists)"
    log_info "Run with --force to reinstall"

    # Quick status check
    if command -v claude &> /dev/null; then
        log_ok "Claude Code: $(claude --version 2>/dev/null || echo 'installed')"
    fi

    echo ""
    log_info "To start Claude Code: claude"
    exit 0
fi

echo ""

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

# Ensure PATH includes potential install locations
export PATH="$HOME/.local/bin:$HOME/.npm-global/bin:$PATH"

# Step 1: Create AI context files
log_info "Step 1/2: Creating AI context files..."

VENDOR_TEMPLATES="${SCIAGENT_SCRIPTS}/../templates/vendor"
PROJECT_NAME=$(basename "$PROJECT_DIR")

# Template substitution function
substitute_template() {
    local template="$1"
    local output="$2"

    if [ -f "$template" ]; then
        cp "$template" "$output"
        sed -i "s|{{PROJECT_ID}}|${PROJECT_NAME}|g" "$output"
        sed -i "s|{{PROJECT_TITLE}}|${PROJECT_NAME}|g" "$output"
        sed -i "s|{{DATE}}|$(date +%Y-%m-%d)|g" "$output"
        sed -i "s|{{SPECIES}}|Mus musculus|g" "$output"
        sed -i "s|{{EXPERIMENTAL_DESIGN}}|TBD|g" "$output"
        log_ok "  Created $(basename "$output")"
    else
        log_warn "  Template not found: $(basename "$template")"
    fi
}

# Only create if not exists (preserve user customizations)
if [ ! -f "${PROJECT_DIR}/CLAUDE.md" ]; then
    substitute_template "${VENDOR_TEMPLATES}/CLAUDE.md.template" "${PROJECT_DIR}/CLAUDE.md"
else
    log_info "  CLAUDE.md already exists (skipped)"
fi

if [ ! -f "${PROJECT_DIR}/GEMINI.md" ]; then
    substitute_template "${VENDOR_TEMPLATES}/GEMINI.md.template" "${PROJECT_DIR}/GEMINI.md"
else
    log_info "  GEMINI.md already exists (skipped)"
fi

if [ ! -f "${PROJECT_DIR}/AGENTS.md" ]; then
    substitute_template "${VENDOR_TEMPLATES}/AGENTS.md.template" "${PROJECT_DIR}/AGENTS.md"
else
    log_info "  AGENTS.md already exists (skipped)"
fi

if [ ! -f "${PROJECT_DIR}/context.md" ]; then
    substitute_template "${VENDOR_TEMPLATES}/context.md.template" "${PROJECT_DIR}/context.md"
else
    log_info "  context.md already exists (skipped)"
fi

# Create analysis_config.yaml in 02_analysis/config/ if directory exists
if [ -d "${PROJECT_DIR}/02_analysis/config" ]; then
    if [ ! -f "${PROJECT_DIR}/02_analysis/config/analysis_config.yaml" ]; then
        substitute_template "${VENDOR_TEMPLATES}/analysis_config.yaml.template" \
            "${PROJECT_DIR}/02_analysis/config/analysis_config.yaml"
    else
        log_info "  analysis_config.yaml already exists (skipped)"
    fi
fi

# Step 2: Activate base role (populate .claude/agents/, .claude/skills/)
log_info "Step 2/2: Activating base role..."

if [ -f "${SCIAGENT_SCRIPTS}/activate-role.sh" ]; then
    bash "${SCIAGENT_SCRIPTS}/activate-role.sh" base --project-dir "${PROJECT_DIR}" || {
        log_warn "Role activation failed. You can run it manually later:"
        log_warn "  ${SCIAGENT_SCRIPTS}/activate-role.sh base --project-dir ${PROJECT_DIR}"
    }
else
    log_warn "activate-role.sh not found - skipping role activation"
fi

echo ""
echo -e "${GREEN}╔══════════════════════════════════════════╗${NC}"
echo -e "${GREEN}║       AI Setup Complete!                 ║${NC}"
echo -e "${GREEN}╚══════════════════════════════════════════╝${NC}"
echo ""

# Summary
log_info "Files created:"
if [ -f "${PROJECT_DIR}/CLAUDE.md" ]; then
    log_ok "  CLAUDE.md, GEMINI.md, AGENTS.md (AI context files)"
fi
if [ -f "${PROJECT_DIR}/context.md" ]; then
    log_ok "  context.md (project scientific context)"
fi
if [ -d "${PROJECT_DIR}/.claude/agents" ]; then
    log_ok "  .claude/agents/, .claude/skills/ (role-activated)"
fi

echo ""
log_info "To start using AI tools:"
echo "  1. Run: claude"
echo ""

log_info "Don't forget to customize:"
echo "  - context.md: Add your scientific question and hypotheses"
echo "  - 02_analysis/config/analysis_config.yaml: Project parameters (if exists)"
echo ""

# --------------------------
# Security Validation
# --------------------------
echo ""
log_info "Running security validation..."
VALIDATE_SCRIPT="${SCIAGENT_SCRIPTS}/validate-secrets.sh"
if [ -f "${VALIDATE_SCRIPT}" ] && [ -x "${VALIDATE_SCRIPT}" ]; then
    if "${VALIDATE_SCRIPT}" "${PROJECT_DIR}"; then
        log_ok "Security validation passed"
    else
        log_warn "Security issues detected!"
        log_warn "Review the warnings above and fix .gitignore before committing."
        echo ""
        log_info "Ensure sensitive files (.env, .claude/) are in .gitignore"
    fi
else
    log_warn "validate-secrets.sh not found - skipping security validation"
    log_info "Ensure sensitive files (.env, .claude/) are in .gitignore"
fi

echo ""
