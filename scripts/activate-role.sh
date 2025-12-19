#!/usr/bin/env bash
#
# activate-role.sh - Activates a role by populating .claude/ directories
#
# Usage:
#   ./activate-role.sh [role-name] [mcp-profile] [--project-dir DIR]
#
# Reads role YAML definition, creates .claude/agents/ and .claude/skills/
# directories, symlinks specified agents/skills from the toolkit, and
# automatically switches to the recommended MCP profile.
#
# Arguments:
#   role-name    - Role to activate (default: base)
#   mcp-profile  - Optional MCP profile override (default: from role YAML)
#
# MCP Profile Behavior:
#   - If mcp-profile is specified, uses that instead of the role's default
#   - If not specified, uses the mcp_profile defined in the role YAML
#   - User can switch profiles independently later with switch-mcp-profile.sh
#
# Roles:
#   base       - Default bioinformatics analysis role (hybrid-research profile)
#   sc-atac    - Single-cell ATAC-seq role
#
# Examples:
#   ./activate-role.sh base                              # Use role's default MCP profile
#   ./activate-role.sh base coding                       # Override with 'coding' profile
#   ./activate-role.sh base --project-dir ~/project      # Specify project directory
#   ./activate-role.sh base research-lite --project-dir ~/project  # Both overrides

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TOOLKIT_DIR="$(dirname "$SCRIPT_DIR")"

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

# Show help first (before any other processing)
for arg in "$@"; do
    if [[ "$arg" == "--help" || "$arg" == "-h" ]]; then
        # Extract only the header comment block (lines 2-29)
        sed -n '2,29p' "$0" | sed 's/^# //' | sed 's/^#//'
        exit 0
    fi
done

# Default values
ROLE="base"
MCP_PROFILE_OVERRIDE=""
PROJECT_DIR="${PWD}"

# Parse positional arguments first (role and optional mcp-profile)
if [[ $# -gt 0 && ! "$1" =~ ^-- ]]; then
    ROLE="$1"
    shift
fi

# Check if next argument is an MCP profile override (not starting with --)
if [[ $# -gt 0 && ! "$1" =~ ^-- ]]; then
    MCP_PROFILE_OVERRIDE="$1"
    shift
fi

# Parse remaining options
while [[ $# -gt 0 ]]; do
    case $1 in
        --project-dir)
            PROJECT_DIR="$2"
            shift 2
            ;;
        *)
            log_error "Unknown option: $1"
            exit 1
            ;;
    esac
done

# Resolve paths
ROLE_FILE="${TOOLKIT_DIR}/roles/${ROLE}.yaml"
AGENTS_DIR="${TOOLKIT_DIR}/agents"
SKILLS_DIR="${TOOLKIT_DIR}/skills"
CLAUDE_DIR="${PROJECT_DIR}/.claude"

# Check if role file exists
if [ ! -f "${ROLE_FILE}" ]; then
    log_error "Role not found: ${ROLE}"
    log_info "Available roles:"
    ls -1 "${TOOLKIT_DIR}/roles/"*.yaml 2>/dev/null | xargs -I{} basename {} .yaml | sed 's/^/  - /' || echo "  (none)"
    exit 1
fi

# Parse YAML list items (simple sed/grep for portability, no yq dependency)
# Handles format:
#   agents:
#     - bioinf-librarian            # inline comments are stripped
#     - rnaseq-methods-writer
get_yaml_list() {
    local file="$1"
    local key="$2"
    # Extract section from key: until next top-level key (line starting with letter + colon)
    # Then extract list items (lines with "- item")
    # Strip inline comments (everything after #) and trailing whitespace
    sed -n "/^${key}:/,/^[a-z]*:/p" "$file" | \
        grep '^[[:space:]]*-' | \
        sed 's/^[[:space:]]*-[[:space:]]*//' | \
        sed 's/[[:space:]]*#.*//' | \
        sed 's/[[:space:]]*$//' | \
        tr -d '\r' | \
        grep -v '^$'
}

# Get role metadata
get_yaml_value() {
    local file="$1"
    local key="$2"
    grep "^${key}:" "$file" | sed "s/${key}:\s*//" | tr -d '\r"' | xargs
}

log_info "Activating role: ${ROLE}"

# Read role metadata
ROLE_NAME=$(get_yaml_value "$ROLE_FILE" "name")
ROLE_DESC=$(get_yaml_value "$ROLE_FILE" "description")
MCP_PROFILE=$(get_yaml_value "$ROLE_FILE" "mcp_profile")

if [ -n "$ROLE_DESC" ]; then
    log_info "Description: ${ROLE_DESC}"
fi

# Create directories
mkdir -p "${CLAUDE_DIR}/agents"
mkdir -p "${CLAUDE_DIR}/skills"

# Clear existing symlinks (for role switching)
find "${CLAUDE_DIR}/agents" -type l -delete 2>/dev/null || true
find "${CLAUDE_DIR}/skills" -type l -delete 2>/dev/null || true

# Symlink agents
AGENTS_COUNT=0
for agent in $(get_yaml_list "$ROLE_FILE" "agents"); do
    src="${AGENTS_DIR}/${agent}.md"
    if [ -f "$src" ]; then
        ln -sf "$src" "${CLAUDE_DIR}/agents/${agent}.md"
        log_ok "  Agent: ${agent}"
        AGENTS_COUNT=$((AGENTS_COUNT + 1))
    else
        log_warn "  Agent not found: ${agent} (expected at ${src})"
    fi
done

# Symlink skills
SKILLS_COUNT=0
for skill in $(get_yaml_list "$ROLE_FILE" "skills"); do
    src="${SKILLS_DIR}/${skill}.md"
    if [ -f "$src" ]; then
        ln -sf "$src" "${CLAUDE_DIR}/skills/${skill}.md"
        log_ok "  Skill: ${skill}"
        SKILLS_COUNT=$((SKILLS_COUNT + 1))
    else
        log_warn "  Skill not found: ${skill} (expected at ${src})"
    fi
done

echo ""
log_ok "Role '${ROLE}' activated"
log_info "  Agents: ${AGENTS_COUNT}"
log_info "  Skills: ${SKILLS_COUNT}"
log_info "  Location: ${CLAUDE_DIR}/"

# Determine which MCP profile to use
FINAL_MCP_PROFILE=""
if [ -n "$MCP_PROFILE_OVERRIDE" ]; then
    FINAL_MCP_PROFILE="$MCP_PROFILE_OVERRIDE"
    log_info "Using MCP profile override: ${FINAL_MCP_PROFILE}"
elif [ -n "$MCP_PROFILE" ]; then
    FINAL_MCP_PROFILE="$MCP_PROFILE"
    log_info "Using role's recommended MCP profile: ${FINAL_MCP_PROFILE}"
fi

# Automatically switch MCP profile
if [ -n "$FINAL_MCP_PROFILE" ]; then
    echo ""
    log_info "Switching MCP profile..."

    # Call switch-mcp-profile.sh
    SWITCH_SCRIPT="${SCRIPT_DIR}/switch-mcp-profile.sh"
    if [ -x "$SWITCH_SCRIPT" ]; then
        if "$SWITCH_SCRIPT" "$FINAL_MCP_PROFILE" --project-dir "$PROJECT_DIR"; then
            echo ""
            log_ok "MCP profile '${FINAL_MCP_PROFILE}' applied"
        else
            log_warn "Failed to switch MCP profile. You can try manually:"
            log_info "  ${SWITCH_SCRIPT} ${FINAL_MCP_PROFILE} --project-dir ${PROJECT_DIR}"
        fi
    else
        log_warn "switch-mcp-profile.sh not found at ${SWITCH_SCRIPT}"
        log_info "To apply MCP profile manually: ./switch-mcp-profile.sh ${FINAL_MCP_PROFILE}"
    fi

    echo ""
    log_info "To switch to a different MCP profile later:"
    log_info "  ${SCRIPT_DIR}/switch-mcp-profile.sh <profile-name>"
else
    log_warn "No MCP profile defined for role '${ROLE}'"

    # Run security validation even without MCP profile switch
    echo ""
    log_info "Running security validation..."
    VALIDATE_SCRIPT="${SCRIPT_DIR}/validate-secrets.sh"
    if [ -f "${VALIDATE_SCRIPT}" ] && [ -x "${VALIDATE_SCRIPT}" ]; then
        if "${VALIDATE_SCRIPT}" "${PROJECT_DIR}"; then
            log_ok "Security validation passed"
        else
            log_warn "Security issues detected - see warnings above"
        fi
    fi
fi
