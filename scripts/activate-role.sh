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
# Arguments:
#   role-name    - Role to activate (default: base)
#
# Roles:
#   base       - Default bioinformatics analysis role
#   sc-atac    - Single-cell ATAC-seq role
#
# Examples:
#   ./activate-role.sh base                              # Activate base role
#   ./activate-role.sh base --project-dir ~/project      # Specify project directory

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
        # Extract only the header comment block (lines 2-22)
        sed -n '2,22p' "$0" | sed 's/^# //' | sed 's/^#//'
        exit 0
    fi
done

# Default values
ROLE="base"
PROJECT_DIR="${PWD}"

# Parse positional arguments first (role name)
if [[ $# -gt 0 && ! "$1" =~ ^-- ]]; then
    ROLE="$1"
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
COMMANDS_DIR="${TOOLKIT_DIR}/commands"
SYSTEM_PROMPTS_DIR="${TOOLKIT_DIR}/system-prompts"
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
    grep "^${key}:" "$file" 2>/dev/null | sed "s/${key}:\s*//" | tr -d '\r"' | xargs || true
}

log_info "Activating role: ${ROLE}"

# Read role metadata
ROLE_NAME=$(get_yaml_value "$ROLE_FILE" "name")
ROLE_DESC=$(get_yaml_value "$ROLE_FILE" "description")
OUTPUT_STYLE=$(get_yaml_value "$ROLE_FILE" "output_style")

if [ -n "$ROLE_DESC" ]; then
    log_info "Description: ${ROLE_DESC}"
fi

# Create directories
mkdir -p "${CLAUDE_DIR}/agents"
mkdir -p "${CLAUDE_DIR}/skills"
mkdir -p "${CLAUDE_DIR}/commands"
mkdir -p "${CLAUDE_DIR}/output-styles"

# Clear existing symlinks (for role switching)
find "${CLAUDE_DIR}/agents" -type l -delete 2>/dev/null || true
find "${CLAUDE_DIR}/skills" -type l -delete 2>/dev/null || true
find "${CLAUDE_DIR}/commands" -type l -delete 2>/dev/null || true
find "${CLAUDE_DIR}/output-styles" -type l -delete 2>/dev/null || true

# Symlink every system-prompt into .claude/output-styles/ so Claude Code
# can discover them all. One is selected at a time via settings.json.
OUTPUT_STYLES_COUNT=0
if [ -d "$SYSTEM_PROMPTS_DIR" ]; then
    for sp in "${SYSTEM_PROMPTS_DIR}"/*.md; do
        [ -f "$sp" ] || continue
        ln -sf "$sp" "${CLAUDE_DIR}/output-styles/$(basename "$sp")"
        OUTPUT_STYLES_COUNT=$((OUTPUT_STYLES_COUNT + 1))
    done
fi

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
# Resolves skills in two formats (directory-first, then flat):
#   1. Canonical: skills/<name>/SKILL.md   -> symlinks the whole directory
#   2. Legacy:    skills/<name>.md         -> symlinks the single file
SKILLS_COUNT=0
for skill in $(get_yaml_list "$ROLE_FILE" "skills"); do
    dir_src="${SKILLS_DIR}/${skill}"
    file_src="${SKILLS_DIR}/${skill}.md"
    if [ -f "${dir_src}/SKILL.md" ]; then
        # -sfn: replace existing symlink without following it (critical for dirs)
        ln -sfn "$dir_src" "${CLAUDE_DIR}/skills/${skill}"
        log_ok "  Skill: ${skill} (dir)"
        SKILLS_COUNT=$((SKILLS_COUNT + 1))
    elif [ -f "$file_src" ]; then
        ln -sf "$file_src" "${CLAUDE_DIR}/skills/${skill}.md"
        log_ok "  Skill: ${skill} (flat)"
        SKILLS_COUNT=$((SKILLS_COUNT + 1))
    else
        log_warn "  Skill not found: ${skill} (expected ${dir_src}/SKILL.md or ${file_src})"
    fi
done

# Symlink commands (optional — role YAML may omit commands: list)
# Commands live at toolkit/commands/<name>.md and symlink into .claude/commands/
COMMANDS_COUNT=0
for cmd in $(get_yaml_list "$ROLE_FILE" "commands"); do
    src="${COMMANDS_DIR}/${cmd}.md"
    if [ -f "$src" ]; then
        ln -sf "$src" "${CLAUDE_DIR}/commands/${cmd}.md"
        log_ok "  Command: /${cmd}"
        COMMANDS_COUNT=$((COMMANDS_COUNT + 1))
    else
        log_warn "  Command not found: ${cmd} (expected at ${src})"
    fi
done

# Pin the output style in .claude/settings.json if the role specifies one.
# Merges into existing settings.json (preserves other keys). Requires python3.
if [ -n "$OUTPUT_STYLE" ]; then
    SETTINGS_FILE="${CLAUDE_DIR}/settings.json"
    if command -v python3 >/dev/null 2>&1; then
        python3 - "$SETTINGS_FILE" "$OUTPUT_STYLE" <<'PY'
import json, sys, os
path, style = sys.argv[1], sys.argv[2]
data = {}
if os.path.exists(path):
    try:
        with open(path) as f:
            data = json.load(f)
    except (json.JSONDecodeError, OSError):
        data = {}
data["outputStyle"] = style
os.makedirs(os.path.dirname(path), exist_ok=True)
with open(path, "w") as f:
    json.dump(data, f, indent=2)
    f.write("\n")
PY
        log_ok "  Output style: ${OUTPUT_STYLE} (pinned in settings.json)"
    else
        log_warn "  python3 not found — cannot pin outputStyle=${OUTPUT_STYLE}"
        log_info "  Set manually: add \"outputStyle\": \"${OUTPUT_STYLE}\" to ${SETTINGS_FILE}"
    fi
fi

echo ""
log_ok "Role '${ROLE}' activated"
log_info "  Agents:        ${AGENTS_COUNT}"
log_info "  Skills:        ${SKILLS_COUNT}"
log_info "  Commands:      ${COMMANDS_COUNT}"
log_info "  Output styles: ${OUTPUT_STYLES_COUNT} available${OUTPUT_STYLE:+, active=${OUTPUT_STYLE}}"
log_info "  Location:      ${CLAUDE_DIR}/"

# Run security validation
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
