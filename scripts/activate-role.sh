#!/usr/bin/env bash
#
# activate-role.sh - Activates a role by populating .claude/ directories
#
# Usage:
#   ./activate-role.sh [role-name] [--project-dir DIR]
#
# Reads role YAML definition, creates .claude/agents/ and .claude/skills/
# directories, and symlinks specified agents/skills from the toolkit.
#
# Roles:
#   base       - Default bioinformatics analysis role (bioinf-librarian, rnaseq-methods-writer)
#   coding     - Code-focused role (future)
#   research   - Research-focused role (future)
#
# Examples:
#   ./activate-role.sh base                          # Use current directory
#   ./activate-role.sh base --project-dir ~/project  # Specify project directory

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

# Default values
ROLE="${1:-base}"
PROJECT_DIR="${PWD}"

# Parse arguments
shift || true  # Shift past role name if provided
while [[ $# -gt 0 ]]; do
    case $1 in
        --project-dir)
            PROJECT_DIR="$2"
            shift 2
            ;;
        --help|-h)
            grep '^#' "$0" | grep -v '#!/usr/bin/env' | sed 's/^# //' | sed 's/^#//'
            exit 0
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
#     - bioinf-librarian
#     - rnaseq-methods-writer
get_yaml_list() {
    local file="$1"
    local key="$2"
    # Extract section from key: until next top-level key (line starting with letter + colon)
    # Then extract list items (lines with "- item")
    sed -n "/^${key}:/,/^[a-z]*:/p" "$file" | \
        grep '^[[:space:]]*-' | \
        sed 's/^[[:space:]]*-[[:space:]]*//' | \
        tr -d '\r'
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

# Suggest MCP profile if defined
if [ -n "$MCP_PROFILE" ]; then
    echo ""
    log_info "Recommended MCP profile: ${MCP_PROFILE}"
    log_info "  To apply: ./switch-mcp-profile.sh ${MCP_PROFILE}"
fi
