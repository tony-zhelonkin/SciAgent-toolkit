#!/usr/bin/env bash
#
# MCP Server Configuration Script
#
# Automatically configures all installed MCP servers for Claude Code.
# Creates .mcp.json file with proper server definitions.
#
# Usage:
#   ./configure_mcp_servers.sh [OPTIONS]
#
# Options:
#   --project-dir DIR      Project directory for .mcp.json (default: current dir)
#   --force                Overwrite existing configuration
#   --dry-run              Show what would be configured without making changes
#   --help                 Show this help message

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
log_error() { echo -e "${RED}[ERROR]${NC} $*"; }

# Default settings
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# Default to current working directory, not script location
# This allows the script to be run from a submodule while creating config in user's project
PROJECT_DIR="${PWD}"
FORCE_OVERWRITE=false
DRY_RUN=false

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --project-dir)
            PROJECT_DIR="$2"
            shift 2
            ;;
        --force)
            FORCE_OVERWRITE=true
            shift
            ;;
        --dry-run)
            DRY_RUN=true
            shift
            ;;
        --help)
            grep '^#' "$0" | grep -v '#!/usr/bin/env' | sed 's/^# //' | sed 's/^#//'
            exit 0
            ;;
        *)
            log_error "Unknown option: $1"
            exit 1
            ;;
    esac
done

log_info "Configuring MCP servers..."
log_info "Project directory: ${PROJECT_DIR}"

# Configuration file path
MCP_CONFIG="${PROJECT_DIR}/.mcp.json"

# Check for existing configuration
if [ -f "${MCP_CONFIG}" ] && [ "${FORCE_OVERWRITE}" = false ]; then
    log_warn "Configuration file already exists: ${MCP_CONFIG}"
    log_info "Use --force to overwrite, or manually edit the file"
    exit 0
fi

# Detect installed MCP servers
log_info "Detecting installed MCP servers..."

SERVERS_TO_CONFIGURE=()

# Check for Sequential Thinking
if command -v npx &>/dev/null; then
    log_ok "Sequential Thinking available (npx found)"
    SERVERS_TO_CONFIGURE+=("sequential-thinking")
else
    log_warn "Sequential Thinking not available (npx not found)"
fi

# Check for ToolUniverse (check multiple possible locations)
# When running as submodule, the env may be in various places
TOOLUNIVERSE_ENV=""
for search_path in \
    "${PROJECT_DIR}/scripts/tooluniverse-env" \
    "${PROJECT_DIR}/tooluniverse-env" \
    "${SCRIPT_DIR}/tooluniverse-env" \
    "${SCRIPT_DIR}/../tooluniverse-env" \
    ; do
    if [ -d "$search_path" ]; then
        # Convert to absolute path
        TOOLUNIVERSE_ENV="$(cd "$search_path" && pwd)"
        break
    fi
done
# Also search common submodule paths
if [ -z "${TOOLUNIVERSE_ENV}" ]; then
    for subdir in "${PROJECT_DIR}/"*"/SciAgent-toolkit/scripts/tooluniverse-env" \
                  "${PROJECT_DIR}/"*"/"*"/SciAgent-toolkit/scripts/tooluniverse-env"; do
        if [ -d "$subdir" ] 2>/dev/null; then
            TOOLUNIVERSE_ENV="$(cd "$subdir" && pwd)"
            break
        fi
    done
fi
if [ -n "${TOOLUNIVERSE_ENV}" ] && [ -d "${TOOLUNIVERSE_ENV}" ]; then
    # Detect correct command name (use -h instead of --version)
    if uv --directory "${TOOLUNIVERSE_ENV}" run tooluniverse-mcp -h &>/dev/null; then
        TOOLUNIVERSE_CMD="tooluniverse-mcp"
    elif uv --directory "${TOOLUNIVERSE_ENV}" run tooluniverse-smcp-stdio -h &>/dev/null; then
        TOOLUNIVERSE_CMD="tooluniverse-smcp-stdio"
    else
        TOOLUNIVERSE_CMD=""
    fi

    if [ -n "${TOOLUNIVERSE_CMD}" ]; then
        log_ok "ToolUniverse available at ${TOOLUNIVERSE_ENV} (${TOOLUNIVERSE_CMD})"
        SERVERS_TO_CONFIGURE+=("tooluniverse")
    else
        log_warn "ToolUniverse installed but command not working"
    fi
else
    log_warn "ToolUniverse not available (not found in scripts/ or project root)"
fi

# Check for Serena (optional - may take long to build on first run)
if command -v uvx &>/dev/null; then
    log_info "Checking Serena availability (this may take a moment)..."
    # Quick test with short timeout - if it times out, skip Serena
    if timeout 30 uvx --from git+https://github.com/oraios/serena serena --version &>/dev/null 2>&1; then
        log_ok "Serena available (pre-cached)"
        SERVERS_TO_CONFIGURE+=("serena")
    else
        log_warn "Serena not pre-cached (would take 5-15 min to build)"
        log_info "Skipping Serena - install manually later if needed"
        log_info "  To install: uvx --from git+https://github.com/oraios/serena serena start-mcp-server --help"
    fi
else
    log_warn "Serena not available (uvx not found)"
fi

# Check for PAL
if command -v uvx &>/dev/null; then
    log_info "Checking PAL availability..."
    if timeout 30 uvx --from git+https://github.com/BeehiveInnovations/pal-mcp-server.git pal-mcp-server --help &>/dev/null 2>&1; then
        log_ok "PAL available"
        SERVERS_TO_CONFIGURE+=("pal")
    else
        log_warn "PAL not pre-cached or available"
        log_info "Skipping PAL"
    fi
fi

# Summary
log_info "Servers to configure: ${#SERVERS_TO_CONFIGURE[@]}"
for server in "${SERVERS_TO_CONFIGURE[@]}"; do
    echo "  - ${server}"
done

if [ ${#SERVERS_TO_CONFIGURE[@]} -eq 0 ]; then
    log_error "No MCP servers found to configure"
    log_info "Please run the installation script first:"
    log_info "  ./scripts/setup_mcp_infrastructure.sh"
    exit 1
elif [ ${#SERVERS_TO_CONFIGURE[@]} -lt 2 ]; then
    log_warn "Only ${#SERVERS_TO_CONFIGURE[@]} server(s) will be configured"
    log_info "This is functional but you may want to install more servers"
fi

# Generate .mcp.json configuration
log_info "Generating MCP configuration..."

if [ "${DRY_RUN}" = true ]; then
    log_warn "Dry run mode - not writing configuration file"
    MCP_CONFIG="/tmp/mcp-config-preview.json"
fi

# Build JSON using a more robust approach
# Collect server configs into an array, then join with commas
declare -a SERVER_CONFIGS=()

# Add Sequential Thinking if available
if [[ " ${SERVERS_TO_CONFIGURE[*]} " =~ " sequential-thinking " ]]; then
    SERVER_CONFIGS+=('    "sequential-thinking": {
      "command": "npx",
      "args": [
        "-y",
        "@modelcontextprotocol/server-sequential-thinking"
      ]
    }')
fi

# Add ToolUniverse if available
if [[ " ${SERVERS_TO_CONFIGURE[*]} " =~ " tooluniverse " ]]; then
    # Check for Azure OpenAI credentials
    if [ -n "${AZURE_OPENAI_API_KEY:-}" ] && [ -n "${AZURE_OPENAI_ENDPOINT:-}" ]; then
        SERVER_CONFIGS+=("    \"tooluniverse\": {
      \"command\": \"uv\",
      \"args\": [
        \"--directory\",
        \"${TOOLUNIVERSE_ENV}\",
        \"run\",
        \"${TOOLUNIVERSE_CMD}\",
        \"--exclude-tool-types\",
        \"PackageTool\",
        \"--hook-type\",
        \"SummarizationHook\"
      ],
      \"env\": {
        \"AZURE_OPENAI_API_KEY\": \"${AZURE_OPENAI_API_KEY}\",
        \"AZURE_OPENAI_ENDPOINT\": \"${AZURE_OPENAI_ENDPOINT}\"
      }
    }")
    else
        SERVER_CONFIGS+=("    \"tooluniverse\": {
      \"command\": \"uv\",
      \"args\": [
        \"--directory\",
        \"${TOOLUNIVERSE_ENV}\",
        \"run\",
        \"${TOOLUNIVERSE_CMD}\",
        \"--exclude-tool-types\",
        \"PackageTool\"
      ]
    }")
    fi
fi

# Add Serena if available
if [[ " ${SERVERS_TO_CONFIGURE[*]} " =~ " serena " ]]; then
    SERVER_CONFIGS+=('    "serena": {
      "command": "uvx",
      "args": [
        "--from",
        "git+https://github.com/oraios/serena",
        "serena",
        "start-mcp-server"
      ]
    }')
fi

# Add PAL if available
if [[ " ${SERVERS_TO_CONFIGURE[*]} " =~ " pal " ]]; then
    # Note: Using sh wrapper as per PAL documentation for robust path handling
    SERVER_CONFIGS+=('    "pal": {
      "command": "sh",
      "args": [
        "-c",
        "for p in $(which uvx 2>/dev/null) $HOME/.local/bin/uvx /opt/homebrew/bin/uvx /usr/local/bin/uvx uvx; do [ -x \"$p\" ] && exec \"$p\" --from git+https://github.com/BeehiveInnovations/pal-mcp-server.git pal-mcp-server; done; echo '"'"'uvx not found'"'"' >&2; exit 1"
      ],
      "env": {
        "PATH": "/usr/local/bin:/usr/bin:/bin:/opt/homebrew/bin:~/.local/bin"
      }
    }')
fi

# Write JSON file with proper comma handling
{
    echo '{'
    echo '  "mcpServers": {'

    # Join array elements with ",\n"
    first_config=true
    for config in "${SERVER_CONFIGS[@]}"; do
        if [ "$first_config" = true ]; then
            first_config=false
        else
            echo ','
        fi
        echo -n "$config"
    done
    echo ''

    echo '  }'
    echo '}'
} > "${MCP_CONFIG}"

# Validate JSON
if command -v python3 &>/dev/null; then
    if python3 -m json.tool "${MCP_CONFIG}" > /dev/null 2>&1; then
        log_ok "Configuration file is valid JSON"
    else
        log_error "Generated configuration has invalid JSON syntax"
        log_info "File location: ${MCP_CONFIG}"
        cat "${MCP_CONFIG}"
        exit 1
    fi
else
    log_warn "Python3 not available - skipping JSON validation"
fi

if [ "${DRY_RUN}" = true ]; then
    log_info "Dry run preview:"
    cat "${MCP_CONFIG}"
    rm -f "${MCP_CONFIG}"
    log_info "No changes made (dry run mode)"
else
    log_ok "MCP configuration created: ${MCP_CONFIG}"

    # Create .claude/settings.json to auto-enable MCP servers
    # Without this, Claude Code discovers servers but doesn't enable them
    CLAUDE_SETTINGS_DIR="${PROJECT_DIR}/.claude"
    CLAUDE_SETTINGS_FILE="${CLAUDE_SETTINGS_DIR}/settings.json"

    mkdir -p "${CLAUDE_SETTINGS_DIR}"

    # Build list of server names for enabledMcpjsonServers
    ENABLED_SERVERS_JSON="["
    first_server=true
    for server in "${SERVERS_TO_CONFIGURE[@]}"; do
        if [ "$first_server" = true ]; then
            first_server=false
        else
            ENABLED_SERVERS_JSON+=", "
        fi
        ENABLED_SERVERS_JSON+="\"${server}\""
    done
    ENABLED_SERVERS_JSON+="]"

    # Create or update settings.json
    if [ -f "${CLAUDE_SETTINGS_FILE}" ]; then
        # Check if it already has MCP settings
        if grep -q "enabledMcpjsonServers\|enableAllProjectMcpServers" "${CLAUDE_SETTINGS_FILE}" 2>/dev/null; then
            log_info "Claude settings already configured for MCP servers"
        else
            log_info "Updating existing Claude settings with MCP enablement..."
            # Use python to merge JSON if available
            if command -v python3 &>/dev/null; then
                python3 << EOF
import json
with open("${CLAUDE_SETTINGS_FILE}", "r") as f:
    settings = json.load(f)
settings["enabledMcpjsonServers"] = ${ENABLED_SERVERS_JSON}
with open("${CLAUDE_SETTINGS_FILE}", "w") as f:
    json.dump(settings, f, indent=2)
EOF
                log_ok "Updated Claude settings with MCP enablement"
            else
                log_warn "Cannot update settings.json - python3 not available"
            fi
        fi
    else
        # Create new settings file
        cat > "${CLAUDE_SETTINGS_FILE}" << EOF
{
  "enabledMcpjsonServers": ${ENABLED_SERVERS_JSON}
}
EOF
        log_ok "Created Claude settings: ${CLAUDE_SETTINGS_FILE}"
    fi

    log_info ""
    log_info "To verify the configuration:"
    log_info "  1. Start Claude Code: claude"
    log_info "  2. Check MCP servers: /mcp"
    log_info ""
    log_info "Available MCP servers:"
    for server in "${SERVERS_TO_CONFIGURE[@]}"; do
        echo "  - ${server}"
    done
fi

log_ok "MCP configuration complete"
