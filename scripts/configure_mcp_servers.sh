#!/usr/bin/env bash
#
# MCP Server Configuration Script
#
# Automatically configures all installed MCP servers for Claude Code.
# Uses `claude mcp add` CLI for robust configuration.
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

# Remove existing config if force mode
if [ -f "${MCP_CONFIG}" ] && [ "${FORCE_OVERWRITE}" = true ]; then
    log_info "Removing existing configuration (--force mode)"
    rm -f "${MCP_CONFIG}"
fi

# Load environment variables from .env files if present
if [ -f "${PROJECT_DIR}/.env" ]; then
    log_info "Loading environment variables from ${PROJECT_DIR}/.env"
    set -a
    source "${PROJECT_DIR}/.env"
    set +a
fi
if [ -f "${PROJECT_DIR}/.devcontainer/.env" ]; then
    log_info "Loading environment variables from ${PROJECT_DIR}/.devcontainer/.env"
    set -a
    source "${PROJECT_DIR}/.devcontainer/.env"
    set +a
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

# Check for Context7 (API key is OPTIONAL per official docs - works without it)
if command -v npx &>/dev/null; then
    log_ok "Context7 available (npx found)"
    SERVERS_TO_CONFIGURE+=("context7")
    if [ -z "${CONTEXT7_API_KEY:-}" ]; then
        log_info "Tip: Set CONTEXT7_API_KEY for higher rate limits (https://context7.com/dashboard)"
    fi
fi

# Check for ToolUniverse
TOOLUNIVERSE_ENV=""
TOOLUNIVERSE_CMD=""

# Define potential search paths for tooluniverse-env
POTENTIAL_PATHS=(
    "${PROJECT_DIR}/tooluniverse-env"
    "${SCRIPT_DIR}/tooluniverse-env"
    "${SCRIPT_DIR}/../tooluniverse-env"
    "${PROJECT_DIR}/scripts/tooluniverse-env"
    "/workspaces/12868-EH/01_scripts/SciAgent-toolkit/scripts/tooluniverse-env"
)

# Search for the environment
for search_path in "${POTENTIAL_PATHS[@]}"; do
    if [ -d "$search_path" ]; then
        TOOLUNIVERSE_ENV="$(cd "$search_path" && pwd)"
        break
    fi
done

# If not found, search submodule paths
if [ -z "${TOOLUNIVERSE_ENV}" ]; then
    for subdir in "${PROJECT_DIR}/"*"/SciAgent-toolkit/scripts/tooluniverse-env" \
                  "${PROJECT_DIR}/"*"/"*"/SciAgent-toolkit/scripts/tooluniverse-env"; do
        if [ -d "$subdir" ] 2>/dev/null; then
            TOOLUNIVERSE_ENV="$(cd "$subdir" && pwd)"
            break
        fi
    done
fi

if [ -n "${TOOLUNIVERSE_ENV}" ]; then
    # Detect command name - prefer stdio versions for Claude Code compatibility
    # tooluniverse-mcp uses HTTP by default (1.0.14+), tooluniverse-smcp-stdio uses stdio
    if uv --directory "${TOOLUNIVERSE_ENV}" run tooluniverse-smcp-stdio -h &>/dev/null 2>&1; then
        TOOLUNIVERSE_CMD="tooluniverse-smcp-stdio"
    elif uv --directory "${TOOLUNIVERSE_ENV}" run python -m tooluniverse.smcp_server --transport stdio -h &>/dev/null 2>&1; then
        TOOLUNIVERSE_CMD="python -m tooluniverse.smcp_server --transport stdio"
    elif uv --directory "${TOOLUNIVERSE_ENV}" run tooluniverse-mcp -h &>/dev/null 2>&1; then
        # Fallback to HTTP version (may not work with Claude Code)
        TOOLUNIVERSE_CMD="tooluniverse-mcp"
        log_warn "Only HTTP transport available - may not work with Claude Code"
    fi
    if [ -n "${TOOLUNIVERSE_CMD}" ]; then
        log_ok "ToolUniverse available at ${TOOLUNIVERSE_ENV}"
        SERVERS_TO_CONFIGURE+=("tooluniverse")
    else
        log_warn "ToolUniverse installed but command not working"
    fi
else
    log_warn "ToolUniverse not available (env directory not found)"
fi

# Check for Serena
if command -v uvx &>/dev/null; then
    log_info "uvx found, enabling Serena configuration..."
    SERVERS_TO_CONFIGURE+=("serena")
fi

# Check for PAL (always configure, warn if no keys)
if command -v uvx &>/dev/null; then
    log_info "uvx found, enabling PAL configuration..."
    SERVERS_TO_CONFIGURE+=("pal")

    # Check if any API key is set and warn if not
    if [ -z "${GEMINI_API_KEY:-}" ] && [ -z "${OPENAI_API_KEY:-}" ] && \
       [ -z "${XAI_API_KEY:-}" ] && [ -z "${OPENROUTER_API_KEY:-}" ]; then
        log_warn "PAL configured but NO API keys found!"
        log_info "PAL requires at least one AI provider key to function."
        log_info "Set one of these in .devcontainer/.env:"
        log_info "  - GEMINI_API_KEY  (https://aistudio.google.com/apikey)"
        log_info "  - OPENAI_API_KEY  (https://platform.openai.com/api-keys)"
        log_info "  - XAI_API_KEY     (https://console.x.ai/)"
        log_info "Then run: ./scripts/configure_mcp_servers.sh --force"
    fi
fi

# Summary
log_info "Servers to configure: ${#SERVERS_TO_CONFIGURE[@]}"
for server in "${SERVERS_TO_CONFIGURE[@]}"; do
    echo "  - ${server}"
done

if [ ${#SERVERS_TO_CONFIGURE[@]} -eq 0 ]; then
    log_error "No MCP servers found to configure"
    exit 1
fi

if [ "${DRY_RUN}" = true ]; then
    log_warn "Dry run mode - showing what would be configured"
    exit 0
fi

# Change to project directory for claude mcp commands
cd "${PROJECT_DIR}"

# Helper to build env args
build_env_args() {
    local args=""
    for var in "$@"; do
        if [ -n "${!var:-}" ]; then
            args="$args --env $var=${!var}"
        fi
    done
    echo "$args"
}

# Check if claude CLI is available
if command -v claude &>/dev/null; then
    log_info "Using claude CLI to configure MCP servers..."

    # Configure each server using claude mcp add
    for server in "${SERVERS_TO_CONFIGURE[@]}"; do
        case "$server" in
            sequential-thinking)
                log_info "Adding sequential-thinking server..."
                if claude mcp add --transport stdio --scope project sequential-thinking -- \
                    npx -y @modelcontextprotocol/server-sequential-thinking; then
                    log_ok "Added sequential-thinking"
                else
                    log_warn "Failed to add sequential-thinking via CLI, will add to JSON"
                fi
                ;;
            tooluniverse)
                log_info "Adding tooluniverse server..."
                # Note: --exclude-tool-types removed in ToolUniverse 1.0.14+
                # Use --compact-mode for minimal tool set if needed
                if claude mcp add --transport stdio --scope project tooluniverse -- \
                    uv --directory "${TOOLUNIVERSE_ENV}" run ${TOOLUNIVERSE_CMD}; then
                    log_ok "Added tooluniverse"
                else
                    log_warn "Failed to add tooluniverse via CLI, will add to JSON"
                fi
                ;;
            serena)
                log_info "Adding serena server..."
                if claude mcp add --transport stdio --scope project serena -- \
                    uvx --from git+https://github.com/oraios/serena serena start-mcp-server; then
                    log_ok "Added serena"
                else
                    log_warn "Failed to add serena via CLI"
                fi
                ;;
            pal)
                log_info "Adding pal server..."
                PAL_ENV_ARGS=$(build_env_args GEMINI_API_KEY OPENAI_API_KEY XAI_API_KEY OPENROUTER_API_KEY DIAL_API_KEY)
                # Note: variables must be expanded before passing to claude mcp add because it runs as a subprocess
                # Using eval to handle the dynamic args construction safely
                CMD="claude mcp add --transport stdio --scope project pal $PAL_ENV_ARGS -- uvx --from git+https://github.com/BeehiveInnovations/pal-mcp-server.git pal-mcp-server"
                if eval "$CMD"; then
                    log_ok "Added pal"
                else
                    log_warn "Failed to add pal via CLI"
                fi
                ;;
            context7)
                log_info "Adding context7 server..."
                if [ -n "${CONTEXT7_API_KEY:-}" ]; then
                    # With API key (higher rate limits)
                    if claude mcp add --transport stdio --scope project context7 \
                        --env CONTEXT7_API_KEY="${CONTEXT7_API_KEY}" -- \
                        npx -y @upstash/context7-mcp; then
                        log_ok "Added context7 (with API key)"
                    else
                        log_warn "Failed to add context7 via CLI"
                    fi
                else
                    # Without API key (basic rate limits - still works)
                    if claude mcp add --transport stdio --scope project context7 -- \
                        npx -y @upstash/context7-mcp; then
                        log_ok "Added context7 (no API key - basic rate limits)"
                    else
                        log_warn "Failed to add context7 via CLI"
                    fi
                fi
                ;;
        esac
    done

    log_ok "MCP servers configured via claude CLI"
    log_info "Verify with: claude mcp list"
else
    log_warn "claude CLI not found, generating .mcp.json manually"

    # Generate JSON manually as fallback
    # Using python for reliable JSON generation
    if command -v python3 &>/dev/null; then
        python3 << PYEOF
import json
import os

servers = {}
env_vars = os.environ

# Sequential Thinking
if "sequential-thinking" in """${SERVERS_TO_CONFIGURE[*]}""".split():
    servers["sequential-thinking"] = {
        "type": "stdio",
        "command": "npx",
        "args": ["-y", "@modelcontextprotocol/server-sequential-thinking"]
    }

# ToolUniverse (--exclude-tool-types removed in 1.0.14+)
if "tooluniverse" in """${SERVERS_TO_CONFIGURE[*]}""".split():
    env_path = """${TOOLUNIVERSE_ENV}"""
    cmd = """${TOOLUNIVERSE_CMD}"""
    servers["tooluniverse"] = {
        "type": "stdio",
        "command": "uv",
        "args": ["--directory", env_path, "run"] + cmd.split()
    }

# Serena
if "serena" in """${SERVERS_TO_CONFIGURE[*]}""".split():
    servers["serena"] = {
        "type": "stdio",
        "command": "uvx",
        "args": ["--from", "git+https://github.com/oraios/serena", "serena", "start-mcp-server"]
    }

# PAL
if "pal" in """${SERVERS_TO_CONFIGURE[*]}""".split():
    pal_env = {}
    for key in ["GEMINI_API_KEY", "OPENAI_API_KEY", "XAI_API_KEY", "OPENROUTER_API_KEY", "DIAL_API_KEY"]:
        # Only include non-empty values
        if key in env_vars and env_vars[key]:
            pal_env[key] = env_vars[key]

    pal_config = {
        "type": "stdio",
        "command": "uvx",
        "args": ["--from", "git+https://github.com/BeehiveInnovations/pal-mcp-server.git", "pal-mcp-server"]
    }
    # Only add env section if we have actual keys
    if pal_env:
        pal_config["env"] = pal_env
    servers["pal"] = pal_config

# Context7 (API key is optional)
if "context7" in """${SERVERS_TO_CONFIGURE[*]}""".split():
    context7_config = {
        "type": "stdio",
        "command": "npx",
        "args": ["-y", "@upstash/context7-mcp"]
    }
    # Only add env if API key is actually set
    api_key = env_vars.get("CONTEXT7_API_KEY", "")
    if api_key:
        context7_config["env"] = {"CONTEXT7_API_KEY": api_key}
    servers["context7"] = context7_config

config = {"mcpServers": servers}

with open("${MCP_CONFIG}", "w") as f:
    json.dump(config, f, indent=2)

print(f"Created {len(servers)} server(s) in .mcp.json")
PYEOF
        log_ok "Generated .mcp.json"
    else
        log_error "Neither claude CLI nor python3 available"
        exit 1
    fi
fi

# Create .claude/settings.json AND settings.local.json to auto-enable servers
CLAUDE_SETTINGS_DIR="${PROJECT_DIR}/.claude"
CLAUDE_SETTINGS_FILE="${CLAUDE_SETTINGS_DIR}/settings.json"
CLAUDE_LOCAL_SETTINGS_FILE="${CLAUDE_SETTINGS_DIR}/settings.local.json"

mkdir -p "${CLAUDE_SETTINGS_DIR}"

# Build enabled servers list
ENABLED_LIST=$(printf '"%s",' "${SERVERS_TO_CONFIGURE[@]}" | sed 's/,$//')

update_settings_file() {
    local file="$1"
    if [ -f "$file" ]; then
        if grep -q "enabledMcpjsonServers" "$file" 2>/dev/null; then
            log_info "Updating existing settings file: $file"
        else
            log_info "Adding enabledMcpjsonServers to: $file"
        fi
        
        python3 << PYEOF 2>/dev/null || true
import json
try:
    with open("$file", "r") as f:
        settings = json.load(f)
except:
    settings = {}

enabled_servers = [${ENABLED_LIST}]
settings["enabledMcpjsonServers"] = enabled_servers

# Remove these servers from disabled list if present
if "disabledMcpjsonServers" in settings:
    disabled = set(settings["disabledMcpjsonServers"])
    for server in enabled_servers:
        if server in disabled:
            disabled.remove(server)
    settings["disabledMcpjsonServers"] = list(disabled)

with open("$file", "w") as f:
    json.dump(settings, f, indent=2)
PYEOF
        log_ok "Updated $file"
    else
        cat > "$file" << EOF
{
  "enabledMcpjsonServers": [${ENABLED_LIST}]
}
EOF
        log_ok "Created settings: $file"
    fi
}

update_settings_file "${CLAUDE_SETTINGS_FILE}"
update_settings_file "${CLAUDE_LOCAL_SETTINGS_FILE}"

# Verify configuration
log_info ""
log_info "Configuration complete!"
log_info ""
log_info "To verify:"
log_info "  1. Run: claude mcp list"
log_info "  2. Start Claude Code: claude"
log_info "  3. Inside Claude: /mcp"
log_info ""
log_info "Configured servers:"
for server in "${SERVERS_TO_CONFIGURE[@]}"; do
    echo "  - ${server}"
done

log_ok "MCP configuration complete"
