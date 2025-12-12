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

# Determine script directory FIRST (before any variable references)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TOOLKIT_DIR="$(dirname "${SCRIPT_DIR}")"

# Source common utilities
if [ -f "${SCRIPT_DIR}/common.sh" ]; then
    source "${SCRIPT_DIR}/common.sh"
elif [ -f "${TOOLKIT_DIR}/scripts/common.sh" ]; then
    source "${TOOLKIT_DIR}/scripts/common.sh"
else
    # Fallback if common.sh not found (e.g. standalone script usage)
    echo "Warning: common.sh not found. Define log functions manually."
    log_info()  { echo -e "[INFO] $*"; }
    log_ok()    { echo -e "[OK] $*"; }
    log_warn()  { echo -e "[WARN] $*"; }
    log_error() { echo -e "[ERROR] $*"; }
    
    # Fallback load_env
    load_env() {
        local project_dir="${1:-$PWD}"
        if [ -f "${project_dir}/.env" ]; then set -a; source "${project_dir}/.env"; set +a; fi
        if [ -f "${project_dir}/.devcontainer/.env" ]; then set -a; source "${project_dir}/.devcontainer/.env"; set +a; fi
    }
fi

# Default values
PROJECT_DIR="${PWD}"
PROFILE="coding"  # Default profile

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --project-dir)
            PROJECT_DIR="$2"
            shift 2
            ;;
        --help|-h)
            echo "Usage: $0 [PROFILE] [--project-dir DIR]"
            echo ""
            echo "Profiles:"
            echo "  minimal       (~3k tokens)   - Sequential Thinking + Context7 only"
            echo "  coding        (~25k tokens)  - + PAL for multi-model collaboration"
            echo "  codebase      (~75k tokens)  - + PAL + Serena for large codebase exploration"
            echo "  hybrid        (~35k tokens)  - + PAL + Context7 + Research tools (Gemini only)"
            echo "  research-lite (~30k tokens)  - + Targeted ToolUniverse (6 tools)"
            echo "  research-full (~50k tokens)  - + Broader ToolUniverse (14 tools) + SummarizationHook"
            echo "  full          (~100k tokens) - + PAL + ToolUniverse (14 tools) + Serena"
            exit 0
            ;;
        -*)
            log_error "Unknown option: $1"
            exit 1
            ;;
        *)
            PROFILE="$1"
            # Alias hybrid to hybrid-research
            if [[ "$PROFILE" == "hybrid" ]]; then PROFILE="hybrid-research"; fi
            shift
            ;;
    esac
done

# Validate profile and find profile file
PROFILES_DIR="${TOOLKIT_DIR}/templates/mcp-profiles"
PROFILE_FILE="${PROFILES_DIR}/${PROFILE}.mcp.json"

if [ ! -f "${PROFILE_FILE}" ]; then
    log_error "Profile not found: ${PROFILE}"
    log_info "Available profiles:"
    ls -1 "${PROFILES_DIR}"/*.mcp.json 2>/dev/null | xargs -I{} basename {} .mcp.json | sed 's/^/  - /'
    exit 1
fi

# Auto-detect ToolUniverse environment location
TOOLUNIVERSE_ENV=""
if [ -d "${PROJECT_DIR}/tooluniverse-env" ]; then
    TOOLUNIVERSE_ENV="${PROJECT_DIR}/tooluniverse-env"
elif [ -d "${SCRIPT_DIR}/tooluniverse-env" ]; then
    TOOLUNIVERSE_ENV="${SCRIPT_DIR}/tooluniverse-env"
fi

# Load environment variables for injection
load_env "${PROJECT_DIR}"

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

# ... (previous code)

# --------------------------
# Gemini Configuration
# --------------------------
GEMINI_PROFILES_DIR="${TOOLKIT_DIR}/templates/gemini-profiles"
GEMINI_SETTINGS_DIR="${PROJECT_DIR}/.gemini"
mkdir -p "${GEMINI_SETTINGS_DIR}"

# Map profile name to Gemini profile
case $PROFILE in
    minimal|coding) GEMINI_PROFILE="coding" ;;
    codebase)       GEMINI_PROFILE="codebase" ;;
    research*)      GEMINI_PROFILE="research" ;;
    full)           GEMINI_PROFILE="research" ;; 
    hybrid-research) GEMINI_PROFILE="research" ;; # Claude=Coding, Gemini=Research
    *)              GEMINI_PROFILE="coding" ;;
esac

GEMINI_TEMPLATE="${GEMINI_PROFILES_DIR}/${GEMINI_PROFILE}.json"

if [ -f "${GEMINI_TEMPLATE}" ]; then
    log_info "Configuring Gemini with profile: ${GEMINI_PROFILE}"
    
    # Read template and substitute variables using Python for safety
    if command -v python3 &>/dev/null; then
        python3 -c "
import sys, os

# Read template
try:
    with open('${GEMINI_TEMPLATE}', 'r') as f:
        content = f.read()
except Exception as e:
    sys.stderr.write(f'Error reading template: {e}\\n')
    sys.exit(1)

# Variables to replace
replacements = {
    '\${TOOLUNIVERSE_ENV}': '${TOOLUNIVERSE_ENV}',
    '\${GEMINI_API_KEY}': os.environ.get('GEMINI_API_KEY', ''),
    '\${OPENAI_API_KEY}': os.environ.get('OPENAI_API_KEY', ''),
    '\${CONTEXT7_API_KEY}': os.environ.get('CONTEXT7_API_KEY', ''),
}

for key, value in replacements.items():
    content = content.replace(key, value)

# Write output
try:
    with open('${GEMINI_SETTINGS_DIR}/settings.json', 'w') as f:
        f.write(content)
except Exception as e:
    sys.stderr.write(f'Error writing settings: {e}\\n')
    sys.exit(1)
"
        log_ok "Updated .gemini/settings.json"
        log_info "Note: .gemini/settings.json contains API keys but is git-ignored (safe)."
    else
        log_warn "python3 not found - skipping Gemini configuration"
    fi
else
    log_warn "Gemini profile template not found: ${GEMINI_PROFILE}"
fi

# --------------------------
# Codex Configuration
# --------------------------
CODEX_CONFIG_FILE="${HOME}/.codex/config.toml"
if [ -f "${CODEX_CONFIG_FILE}" ] && command -v python3 &>/dev/null; then
    # Check if toml module is available
    if ! python3 -c "import toml" &>/dev/null; then
        log_warn "Python 'toml' module not installed - skipping Codex configuration"
        log_info "To enable: pip3 install toml (may require sudo in some environments)"
    else
        log_info "Updating Codex configuration..."

        # Determine if ToolUniverse should be enabled for Codex
        # For now, we align Codex with the *Claude* profile settings (from .mcp.json)
        # If the profile has ToolUniverse, we enable it in Codex

        python3 -c "
import sys, os, toml, json

codex_config_path = '${CODEX_CONFIG_FILE}'
mcp_json_path = '${PROJECT_DIR}/.mcp.json'
tooluniverse_env = '${TOOLUNIVERSE_ENV}'

try:
    # Load Codex Config
    with open(codex_config_path, 'r') as f:
        codex_config = toml.load(f)
    
    # Load MCP JSON to check if ToolUniverse is active
    with open(mcp_json_path, 'r') as f:
        mcp_config = json.load(f)
        
    has_tu = 'tooluniverse' in mcp_config.get('mcpServers', {})
    
    if 'mcp_servers' not in codex_config:
        codex_config['mcp_servers'] = {}
        
    if has_tu and tooluniverse_env:
        # Add ToolUniverse to Codex
        # We use the python module invocation for robustness
        codex_config['mcp_servers']['tooluniverse'] = {
            'command': 'uv',
            'args': [
                '--directory', tooluniverse_env,
                'run', 'python', '-m', 'tooluniverse.smcp_server',
                '--transport', 'stdio',
                '--compact-mode' # Default to compact for Codex to save tokens
            ]
        }
    elif 'tooluniverse' in codex_config['mcp_servers']:
        # Remove ToolUniverse if not in profile
        del codex_config['mcp_servers']['tooluniverse']
        
    # Write back
    with open(codex_config_path, 'w') as f:
        toml.dump(codex_config, f)
        
except Exception as e:
    # simplified error handling
    sys.stderr.write(f'Warning: Could not update Codex config: {e}\\n')
"
        log_ok "Updated ~/.codex/config.toml"
    fi
fi

# Update .claude/settings.local.json with enabled servers
CLAUDE_DIR="${PROJECT_DIR}/.claude"
mkdir -p "${CLAUDE_DIR}"

# Extract server names from .mcp.json
if command -v python3 &>/dev/null; then
    python3 << PYEOF
import json
import os

claude_dir = os.path.join(os.environ.get('PROJECT_DIR', ''), '.claude')
mcp_json_path = os.path.join(os.environ.get('PROJECT_DIR', ''), '.mcp.json')
settings_local_path = os.path.join(claude_dir, 'settings.local.json')

os.makedirs(claude_dir, exist_ok=True)

try:
    with open(mcp_json_path) as f:
        config = json.load(f)
    servers = list(config.get("mcpServers", {}).keys())
except FileNotFoundError:
    servers = []
except Exception as e:
    sys.stderr.write(f'Warning: Could not read .mcp.json: {e}\n')
    servers = []

# Define default permissions, including the corrected bash globbing
default_permissions = {
    "allow": [
        "Fs(read:*)",
        "Fs(write:*)",
        "Shell(git:*)",
        "Shell(npx:*)",
        "Shell(npm:*)",
        "Shell(uvx:*)",
        "Bash(for f in :*)",
        "Web(fetch:*)",
        "Web(google_web_search:*)"
    ],
    "deny": []
}

settings_content = {
    "enabledMcpjsonServers": servers,
    "permissions": default_permissions
}

with open(settings_local_path, "w") as f:
    json.dump(settings_content, f, indent=2)

print(f"Updated {settings_local_path}")
PYEOF
else
    log_warn "python3 not found - settings.local.json not updated"
fi

# Show result
log_ok "Switched to profile: ${PROFILE}"
echo ""
echo "Context estimate:"
case $PROFILE in
    minimal)       echo "  ~3k tokens (1.5% of 200k)" ;;
    coding)        echo "  ~25k tokens (12.5% of 200k)" ;;
    codebase)      echo "  ~75k tokens (37.5% of 200k)" ;;
    hybrid-research) echo "  ~35k tokens (17.5% of 200k)" ;;
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
