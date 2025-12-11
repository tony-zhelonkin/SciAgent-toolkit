#!/usr/bin/env bash
#
# MCP Profile Switcher
#
# Switches between pre-configured MCP profiles to manage context window usage.
# Each profile is optimized for different workflows with curated tool sets.
#
# Usage:
#   ./switch-mcp-profile.sh [profile-name] [--project-dir DIR]
#
# Profiles:
#   minimal       (~3k tokens)   - Sequential Thinking + Context7 only
#   coding        (~25k tokens)  - + PAL for multi-model collaboration
#   codebase      (~75k tokens)  - + PAL + Serena for large codebase exploration
#   research-lite (~30k tokens)  - + Targeted ToolUniverse (6 tools)
#   research-full (~50k tokens)  - + Broader ToolUniverse (14 tools) + SummarizationHook
#   full          (~100k tokens) - + PAL + ToolUniverse (14 tools) + Serena
#
# IMPORTANT: No profile loads ALL ToolUniverse tools (600+) as that would
# consume 500k+ tokens and exceed the 200k context window before conversation starts.

set -euo pipefail

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

log_info()  { echo -e "${BLUE}[INFO]${NC} $*"; }
log_ok()    { echo -e "${GREEN}[OK]${NC} $*"; }
log_warn()  { echo -e "${YELLOW}[WARN]${NC} $*"; }
log_error() { echo -e "${RED}[ERROR]${NC} $*"; }

# Defaults
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TOOLKIT_DIR="$(dirname "${SCRIPT_DIR}")"
PROFILES_DIR="${TOOLKIT_DIR}/templates/mcp-profiles"
PROJECT_DIR="${PWD}"
PROFILE=""

# Parse arguments
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
        -*)
            log_error "Unknown option: $1"
            exit 1
            ;;
        *)
            PROFILE="$1"
            shift
            ;;
    esac
done

# Show available profiles if none specified
if [ -z "${PROFILE}" ]; then
    echo "Available MCP profiles:"
    echo ""
    echo "  minimal        (~3k tokens)   Default coding - Sequential Thinking + Context7"
    echo "  coding         (~25k tokens)  + PAL multi-model collaboration"
    echo "  codebase       (~75k tokens)  + PAL + Serena for large codebase exploration"
    echo "  research-lite  (~30k tokens)  + ToolUniverse (6 curated tools)"
    echo "  research-full  (~50k tokens)  + ToolUniverse (14 tools) with summarization"
    echo "  full           (~100k tokens) + PAL + ToolUniverse + Serena (max useful config)"
    echo ""
    echo "Usage: $0 <profile-name> [--project-dir DIR]"
    echo ""
    echo "NOTE: Even 'full' profile uses curated ToolUniverse tools (14 tools)."
    echo "      Loading all 600 tools would exceed the 200k context window!"
    exit 0
fi

# Validate profile exists
PROFILE_FILE="${PROFILES_DIR}/${PROFILE}.mcp.json"
if [ ! -f "${PROFILE_FILE}" ]; then
    log_error "Profile not found: ${PROFILE}"
    echo "Available profiles:"
    ls -1 "${PROFILES_DIR}"/*.mcp.json 2>/dev/null | xargs -n1 basename | sed 's/.mcp.json//' | sed 's/^/  /'
    exit 1
fi

# Detect ToolUniverse environment path
TOOLUNIVERSE_ENV=""
for search_path in \
    "${PROJECT_DIR}/tooluniverse-env" \
    "${SCRIPT_DIR}/tooluniverse-env" \
    "${TOOLKIT_DIR}/scripts/tooluniverse-env"; do
    if [ -d "$search_path" ]; then
        TOOLUNIVERSE_ENV="$(cd "$search_path" && pwd)"
        break
    fi
done

# Copy and process profile
log_info "Switching to profile: ${PROFILE}"

# Replace ${TOOLUNIVERSE_ENV} placeholder with actual path
if [ -n "${TOOLUNIVERSE_ENV}" ]; then
    sed "s|\${TOOLUNIVERSE_ENV}|${TOOLUNIVERSE_ENV}|g" "${PROFILE_FILE}" > "${PROJECT_DIR}/.mcp.json"
else
    # If no ToolUniverse env found, copy as-is (may fail if profile needs it)
    cp "${PROFILE_FILE}" "${PROJECT_DIR}/.mcp.json"
    if grep -q 'TOOLUNIVERSE_ENV' "${PROJECT_DIR}/.mcp.json"; then
        log_warn "ToolUniverse environment not found - tooluniverse server may not work"
        log_info "Run setup-ai.sh first or set TOOLUNIVERSE_ENV manually"
    fi
fi

# Update .claude/settings.local.json with enabled servers
CLAUDE_DIR="${PROJECT_DIR}/.claude"
mkdir -p "${CLAUDE_DIR}"

# Extract server names from .mcp.json
if command -v jq &>/dev/null; then
    SERVERS=$(jq -r '.mcpServers | keys | @json' "${PROJECT_DIR}/.mcp.json")
    echo "{\"enabledMcpjsonServers\": ${SERVERS}}" > "${CLAUDE_DIR}/settings.local.json"
elif command -v python3 &>/dev/null; then
    python3 << PYEOF
import json
with open("${PROJECT_DIR}/.mcp.json") as f:
    config = json.load(f)
servers = list(config.get("mcpServers", {}).keys())
with open("${CLAUDE_DIR}/settings.local.json", "w") as f:
    json.dump({"enabledMcpjsonServers": servers}, f, indent=2)
PYEOF
else
    log_warn "Neither jq nor python3 found - settings.local.json not updated"
fi

# Show result
log_ok "Switched to profile: ${PROFILE}"
echo ""
echo "Context estimate:"
case $PROFILE in
    minimal)       echo "  ~3k tokens (1.5% of 200k)" ;;
    coding)        echo "  ~25k tokens (12.5% of 200k)" ;;
    codebase)      echo "  ~75k tokens (37.5% of 200k)" ;;
    research-lite) echo "  ~30k tokens (15% of 200k)" ;;
    research-full) echo "  ~50k tokens (25% of 200k)" ;;
    full)          echo "  ~100k tokens (50% of 200k)" ;;
    *)             echo "  Unknown" ;;
esac
echo ""
echo "Enabled servers:"
if command -v jq &>/dev/null; then
    jq -r '.mcpServers | keys[]' "${PROJECT_DIR}/.mcp.json" | sed 's/^/  - /'
else
    grep -o '"[^"]*":' "${PROJECT_DIR}/.mcp.json" | head -10 | tr -d '":' | sed 's/^/  - /'
fi
echo ""
log_info "Restart Claude Code to apply: exit && claude"
