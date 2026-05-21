#!/usr/bin/env bash
#
# inject-skills.sh - Symlink specific skills on top of the current role (no wipe)
#
# Unlike activate-role.sh, this appends instead of replacing.
# Use it to layer skills onto any role — including `min`.
#
# Usage:
#   ./inject-skills.sh <skill1> [skill2 ...] [--project-dir DIR]
#
# Examples:
#   ./inject-skills.sh anndata multimodal-anndata-mudata
#   ./inject-skills.sh anndata anndatar-seurat-scanpy-conversion cellranger-multi-to-anndata --project-dir ~/project
#   ./inject-skills.sh --list                    # Show all available skills
#   ./inject-skills.sh --list --project-dir ~/p  # List from specific project
#   ./inject-skills.sh --uninstall <skill>        # Remove a single injected skill
#   ./inject-skills.sh --uninstall-all            # Remove all injected skills only

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TOOLKIT_DIR="$(dirname "$SCRIPT_DIR")"
SKILLS_DIR="${TOOLKIT_DIR}/skills"
CLAUDE_SKILLS_DIR=""
PROJECT_DIR=""

# Colors
RED='\033[0;31m'; GREEN='\033[0;32m'; YELLOW='\033[1;33m'; BLUE='\033[0;34m'; NC='\033[0m'
log_info () { echo -e "${BLUE}[INFO]${NC} $*"; }
log_ok   () { echo -e "${GREEN}[OK]${NC} $*"; }
log_warn () { echo -e "${YELLOW}[WARN]${NC} $*"; }
log_err  () { echo -e "${RED}[ERROR]${NC} $*" >&2; }

# Parse args
ACTION="inject"
SKILLS=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        --project-dir) PROJECT_DIR="$2"; shift 2 ;;
        --list)        ACTION="list"; shift ;;
        --uninstall)   ACTION="uninstall"; SKILLS+=("$2"); shift 2 ;;
        --uninstall-all) ACTION="uninstall-all"; shift ;;
        --help|-h)
            sed -n '2,15p' "$0" | sed 's/^# *//'
            exit 0
            ;;
        *) SKILLS+=("$1"); shift ;;
    esac
done

# Resolve project dir
if [ -z "$PROJECT_DIR" ]; then
    PROJECT_DIR="$PWD"
fi
CLAUDE_SKILLS_DIR="${PROJECT_DIR}/.claude/skills"

# Resolve skills_dir if we're inside SciAgent-toolkit itself
if [ "${PROJECT_DIR}" = "${TOOLKIT_DIR}" ]; then
    SKILLS_DIR="${TOOLKIT_DIR}/skills"
fi

# ---- LIST ----
if [ "$ACTION" = "list" ]; then
    echo "Available skills in ${SKILLS_DIR}/"
    echo ""
    for sk in "${SKILLS_DIR}"/*/; do
        [ -d "$sk" ] || continue
        name=$(basename "$sk")
        [[ "$name" == "_"* ]] && continue
        entry_point="$sk/SKILL.md"
        has_refs=""
        if [ -d "${sk}references" ]; then has_refs=" (with refs)"; fi
        log_info "  $name$has_refs"
    done
    exit 0
fi

# ---- ERROR: not inside a project with .claude/ ----
if [ ! -d "${PROJECT_DIR}/.claude" ]; then
    log_err "No .claude/ found at ${PROJECT_DIR}"
    log_info "Use --project-dir to point to a project directory"
    exit 1
fi

# ---- UNINSTALL ----
if [ "$ACTION" = "uninstall-all" ]; then
    # Only remove skills that are symlinks AND come from our toolkit
    FOUND=0
    for sk_dir in "${CLAUDE_SKILLS_DIR}"/*/; do
        [ -L "$sk_dir" ] || continue
        target=$(readlink -f "$sk_dir" 2>/dev/null || true)
        # Only remove if it points to this toolkit's skills/
        if [[ "$target" == "${TOOLKIT_DIR}/skills/"* ]]; then
            rm -rf "$sk_dir"
            log_ok "  Removed: $(basename "$sk_dir")"
            FOUND=1
        fi
    done
    # Also remove flat-file symlinks
    for sk_file in "${CLAUDE_SKILLS_DIR}"/*.md; do
        [ -L "$sk_file" ] || continue
        target=$(readlink -f "$sk_file" 2>/dev/null || true)
        if [[ "$target" == "${TOOLKIT_DIR}/skills/"* ]]; then
            rm -f "$sk_file"
            log_ok "  Removed: $(basename "$sk_file")"
            FOUND=1
        fi
    done
    [ "$FOUND" -eq 1 ] && log_ok "All injected skills removed" || log_warn "No injected skills found"
    exit 0
fi

if [ "$ACTION" = "uninstall" ]; then
    for skill in "${SKILLS[@]}"; do
        dst_dir="${CLAUDE_SKILLS_DIR}/${skill}"
        dst_file="${CLAUDE_SKILLS_DIR}/${skill}.md"
        rm -rf "$dst_dir" 2>/dev/null || true
        rm -f  "$dst_file" 2>/dev/null || true
        log_ok "Uninstalled: ${skill}"
    done
    exit 0
fi

# ---- INJECT ----
# Resolve the actual skills/ directory (handles being called from within toolkit)
SKILLS_SRC="${SKILLS_DIR}"
if [ ! -d "$SKILLS_SRC" ]; then
    SKILLS_SRC="${TOOLKIT_DIR}/skills"
fi

if [ ! -d "$SKILLS_SRC" ]; then
    log_err "Skills directory not found: ${SKILLS_SRC}"
    exit 1
fi

mkdir -p "$CLAUDE_SKILLS_DIR"

for skill in "${SKILLS[@]}"; do
    dir_src="${SKILLS_SRC}/${skill}"
    flat_src="${SKILLS_SRC}/${skill}.md"

    FOUND=0
    if [ -f "${dir_src}/SKILL.md" ]; then
        # Directory-format skill: symlink the whole directory
        ln -sfn "$dir_src" "${CLAUDE_SKILLS_DIR}/${skill}"
        log_ok "  Skill: ${skill} (dir) -> installed"
        FOUND=1
    elif [ -f "$flat_src" ]; then
        # Flat-format skill: symlink the file
        ln -sf "$flat_src" "${CLAUDE_SKILLS_DIR}/${skill}.md"
        log_ok "  Skill: ${skill} (flat) -> installed"
        FOUND=1
    fi

    if [ "$FOUND" -eq 0 ]; then
        log_warn "  Skill '${skill}' not found (skipped)"
    fi
done

echo ""
# Show what's currently in .claude/skills/
echo "Current state of ${CLAUDE_SKILLS_DIR}/:"
COUNT=0
for entry in "${CLAUDE_SKILLS_DIR}"/*; do
    [ -e "$entry" ] || continue
    COUNT=$((COUNT+1))
    if [ -L "$entry" ]; then
        target=$(readlink "$entry" 2>/dev/null || echo "?")
        echo "  $(basename "$entry") -> $target"
    else
        echo "  $(basename "$entry") [file]"
    fi
done
[ "$COUNT" -eq 0 ] && echo "  (empty)"
echo ""
log_ok "Done — ${#SKILLS[@]} skill(s) processed"
