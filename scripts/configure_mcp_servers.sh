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
PROJECT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
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

# Check for ToolUniverse (check both possible locations)
if [ -d "${PROJECT_DIR}/scripts/tooluniverse-env" ]; then
    TOOLUNIVERSE_ENV="${PROJECT_DIR}/scripts/tooluniverse-env"
elif [ -d "${PROJECT_DIR}/tooluniverse-env" ]; then
    TOOLUNIVERSE_ENV="${PROJECT_DIR}/tooluniverse-env"
else
    TOOLUNIVERSE_ENV=""
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

# Start JSON configuration
cat > "${MCP_CONFIG}" << 'EOF_START'
{
  "mcpServers": {
EOF_START

# Track if we need a comma before the next entry
FIRST_ENTRY=true

# Add Sequential Thinking if available
if [[ " ${SERVERS_TO_CONFIGURE[*]} " =~ " sequential-thinking " ]]; then
    if [ "$FIRST_ENTRY" = false ]; then
        echo "," >> "${MCP_CONFIG}"
    fi
    FIRST_ENTRY=false

    cat >> "${MCP_CONFIG}" << 'EOF_SEQ'
    "sequential-thinking": {
      "command": "npx",
      "args": [
        "-y",
        "@modelcontextprotocol/server-sequential-thinking"
      ]
    }
EOF_SEQ
fi

# Add ToolUniverse if available
if [[ " ${SERVERS_TO_CONFIGURE[*]} " =~ " tooluniverse " ]]; then
    if [ "$FIRST_ENTRY" = false ]; then
        echo "," >> "${MCP_CONFIG}"
    fi
    FIRST_ENTRY=false

    # Check for Azure OpenAI credentials
    if [ -n "${AZURE_OPENAI_API_KEY:-}" ] && [ -n "${AZURE_OPENAI_ENDPOINT:-}" ]; then
        cat >> "${MCP_CONFIG}" << EOF_TOOL
    "tooluniverse": {
      "command": "uv",
      "args": [
        "--directory",
        "${TOOLUNIVERSE_ENV}",
        "run",
        "${TOOLUNIVERSE_CMD}",
        "--exclude-tool-types",
        "PackageTool",
        "--hook-type",
        "SummarizationHook"
      ],
      "env": {
        "AZURE_OPENAI_API_KEY": "${AZURE_OPENAI_API_KEY}",
        "AZURE_OPENAI_ENDPOINT": "${AZURE_OPENAI_ENDPOINT}"
      }
    }
EOF_TOOL
    else
        cat >> "${MCP_CONFIG}" << EOF_TOOL
    "tooluniverse": {
      "command": "uv",
      "args": [
        "--directory",
        "${TOOLUNIVERSE_ENV}",
        "run",
        "${TOOLUNIVERSE_CMD}",
        "--exclude-tool-types",
        "PackageTool"
      ]
    }
EOF_TOOL
    fi
fi

# Add Serena if available
if [[ " ${SERVERS_TO_CONFIGURE[*]} " =~ " serena " ]]; then
    if [ "$FIRST_ENTRY" = false ]; then
        echo "," >> "${MCP_CONFIG}"
    fi
    FIRST_ENTRY=false

    cat >> "${MCP_CONFIG}" << 'EOF_SERENA'
    "serena": {
      "command": "uvx",
      "args": [
        "--from",
        "git+https://github.com/oraios/serena",
        "serena",
        "start-mcp-server"
      ]
    }
EOF_SERENA
fi

# Add PAL if available
if [[ " ${SERVERS_TO_CONFIGURE[*]} " =~ " pal " ]]; then
    if [ "$FIRST_ENTRY" = false ]; then
        echo "," >> "${MCP_CONFIG}"
    fi
    FIRST_ENTRY=false

    # Note: Using sh wrapper as per PAL documentation for robust path handling
    cat >> "${MCP_CONFIG}" << 'EOF_PAL'
    "pal": {
      "command": "sh",
      "args": [
        "-c", 
        "for p in $(which uvx 2>/dev/null) $HOME/.local/bin/uvx /opt/homebrew/bin/uvx /usr/local/bin/uvx uvx; do [ -x \"$p\" ] && exec \"$p\" --from git+https://github.com/BeehiveInnovations/pal-mcp-server.git pal-mcp-server; done; echo 'uvx not found' >&2; exit 1"
      ],
      "env": {
        "PATH": "/usr/local/bin:/usr/bin:/bin:/opt/homebrew/bin:~/.local/bin"
      }
    }
EOF_PAL
fi

# Close JSON structure
cat >> "${MCP_CONFIG}" << 'EOF_END'
  }
}
EOF_END

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
