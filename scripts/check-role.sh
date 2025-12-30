#!/usr/bin/env bash
#
# check-role.sh - Check currently active role
#
# Usage:
#   ./check-role.sh [--project-dir DIR] [--quiet]
#
# Inspects the .claude/ directory to determine which role is currently
# active by comparing symlinked agents/skills against role definitions.
#
# Options:
#   --project-dir DIR  Project directory to check (default: current directory)
#   --quiet, -q        Only output role name (no decorations)
#
# Examples:
#   ./check-role.sh                              # Check current directory
#   ./check-role.sh --project-dir ~/project      # Check specific project
#   ./check-role.sh -q                           # Quiet mode for scripting

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TOOLKIT_DIR="$(dirname "$SCRIPT_DIR")"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
CYAN='\033[0;36m'
NC='\033[0m'

log_info()  { [[ "$QUIET" == "true" ]] || echo -e "${BLUE}[INFO]${NC} $*"; }
log_ok()    { [[ "$QUIET" == "true" ]] || echo -e "${GREEN}[OK]${NC} $*"; }
log_warn()  { [[ "$QUIET" == "true" ]] || echo -e "${YELLOW}[WARN]${NC} $*"; }
log_error() { [[ "$QUIET" == "true" ]] || echo -e "${RED}[ERROR]${NC} $*" >&2; }

# Show help
for arg in "$@"; do
    if [[ "$arg" == "--help" || "$arg" == "-h" ]]; then
        sed -n '2,17p' "$0" | sed 's/^# //' | sed 's/^#//'
        exit 0
    fi
done

# Default values
PROJECT_DIR="${PWD}"
QUIET="false"

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --project-dir)
            PROJECT_DIR="$2"
            shift 2
            ;;
        --quiet|-q)
            QUIET="true"
            shift
            ;;
        *)
            log_error "Unknown option: $1"
            exit 1
            ;;
    esac
done

# Resolve paths
CLAUDE_DIR="${PROJECT_DIR}/.claude"
ROLES_DIR="${TOOLKIT_DIR}/roles"

# Parse YAML list items (same helper as activate-role.sh)
get_yaml_list() {
    local file="$1"
    local key="$2"
    sed -n "/^${key}:/,/^[a-z]*:/p" "$file" | \
        grep '^[[:space:]]*-' | \
        sed 's/^[[:space:]]*-[[:space:]]*//' | \
        sed 's/[[:space:]]*#.*//' | \
        sed 's/[[:space:]]*$//' | \
        tr -d '\r' | \
        grep -v '^$' | \
        sort
}

# Get role metadata
get_yaml_value() {
    local file="$1"
    local key="$2"
    grep "^${key}:" "$file" 2>/dev/null | sed "s/${key}:\s*//" | tr -d '\r"' | xargs
}

# Get currently symlinked items (basename without .md)
get_current_items() {
    local dir="$1"
    if [[ -d "$dir" ]]; then
        find "$dir" -type l -name "*.md" 2>/dev/null | xargs -I{} basename {} .md | sort
    fi
}

# Check if .claude directory exists
if [[ ! -d "$CLAUDE_DIR" ]]; then
    if [[ "$QUIET" == "true" ]]; then
        echo "none"
    else
        log_warn "No .claude directory found at ${CLAUDE_DIR}"
        log_info "No role is currently active"
        log_info "To activate a role, run: activate-role.sh <role-name> --project-dir ${PROJECT_DIR}"
    fi
    exit 0
fi

# Get current agents and skills
CURRENT_AGENTS=$(get_current_items "${CLAUDE_DIR}/agents")
CURRENT_SKILLS=$(get_current_items "${CLAUDE_DIR}/skills")

# Check if any agents/skills are present
if [[ -z "$CURRENT_AGENTS" && -z "$CURRENT_SKILLS" ]]; then
    if [[ "$QUIET" == "true" ]]; then
        echo "none"
    else
        log_warn "No agents or skills symlinked in ${CLAUDE_DIR}"
        log_info "No role is currently active"
        log_info "To activate a role, run: activate-role.sh <role-name> --project-dir ${PROJECT_DIR}"
    fi
    exit 0
fi

# Try to match against each role
MATCHED_ROLE=""
MATCHED_ROLE_NAME=""
PARTIAL_MATCHES=()

for role_file in "${ROLES_DIR}"/*.yaml; do
    [[ -f "$role_file" ]] || continue

    role_id=$(basename "$role_file" .yaml)
    role_name=$(get_yaml_value "$role_file" "name")
    role_agents=$(get_yaml_list "$role_file" "agents")
    role_skills=$(get_yaml_list "$role_file" "skills")

    # Compare agents and skills
    if [[ "$CURRENT_AGENTS" == "$role_agents" && "$CURRENT_SKILLS" == "$role_skills" ]]; then
        MATCHED_ROLE="$role_id"
        MATCHED_ROLE_NAME="$role_name"
        break
    fi

    # Check for partial match (at least one agent matches)
    if [[ -n "$CURRENT_AGENTS" && -n "$role_agents" ]]; then
        # Check if any current agent is in role agents
        while IFS= read -r agent; do
            if echo "$role_agents" | grep -qx "$agent"; then
                PARTIAL_MATCHES+=("$role_id")
                break
            fi
        done <<< "$CURRENT_AGENTS"
    fi
done

# Output results
if [[ -n "$MATCHED_ROLE" ]]; then
    if [[ "$QUIET" == "true" ]]; then
        echo "$MATCHED_ROLE"
    else
        echo ""
        echo -e "${GREEN}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
        echo -e "${GREEN}  Active Role: ${CYAN}${MATCHED_ROLE}${NC}"
        [[ -n "$MATCHED_ROLE_NAME" ]] && echo -e "${GREEN}  Name: ${NC}${MATCHED_ROLE_NAME}"
        echo -e "${GREEN}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
        echo ""

        # Show agents
        if [[ -n "$CURRENT_AGENTS" ]]; then
            log_info "Agents:"
            while IFS= read -r agent; do
                echo "    • $agent"
            done <<< "$CURRENT_AGENTS"
        fi

        # Show skills
        if [[ -n "$CURRENT_SKILLS" ]]; then
            log_info "Skills:"
            while IFS= read -r skill; do
                echo "    • $skill"
            done <<< "$CURRENT_SKILLS"
        fi

        echo ""
        log_info "Location: ${CLAUDE_DIR}/"
    fi
else
    if [[ "$QUIET" == "true" ]]; then
        echo "custom"
    else
        echo ""
        echo -e "${YELLOW}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
        echo -e "${YELLOW}  Active Role: ${CYAN}custom (no exact match)${NC}"
        echo -e "${YELLOW}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
        echo ""

        # Show current configuration
        if [[ -n "$CURRENT_AGENTS" ]]; then
            log_info "Current agents:"
            while IFS= read -r agent; do
                echo "    • $agent"
            done <<< "$CURRENT_AGENTS"
        fi

        if [[ -n "$CURRENT_SKILLS" ]]; then
            log_info "Current skills:"
            while IFS= read -r skill; do
                echo "    • $skill"
            done <<< "$CURRENT_SKILLS"
        fi

        # Show partial matches if any
        if [[ ${#PARTIAL_MATCHES[@]} -gt 0 ]]; then
            echo ""
            log_info "Similar roles:"
            printf '    • %s\n' "${PARTIAL_MATCHES[@]}" | sort -u
        fi

        echo ""
        log_info "Available roles:"
        for role_file in "${ROLES_DIR}"/*.yaml; do
            [[ -f "$role_file" ]] || continue
            role_id=$(basename "$role_file" .yaml)
            echo "    • $role_id"
        done

        echo ""
        log_info "To activate a role: activate-role.sh <role-name> --project-dir ${PROJECT_DIR}"
    fi
fi
