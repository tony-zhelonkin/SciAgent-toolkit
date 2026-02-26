#!/usr/bin/env bash
#
# MCP Addon Manager
#
# Manages standalone MCP server addons that persist across profile switches.
# Addons are layered on top of the active profile — they survive profile changes.
#
# Usage:
#   ./manage-addon.sh <command> [addon-name] [--project-dir DIR]
#
# Commands:
#   list                          Show available + active addons
#   enable <name> [--project-dir] Enable addon, merge into .mcp.json
#   disable <name> [--project-dir] Disable addon, remove from .mcp.json
#   status <name> [--project-dir] Check if addon is active + deps OK
#
# Examples:
#   ./manage-addon.sh list
#   ./manage-addon.sh enable jupyter --project-dir /workspaces/myproject
#   ./manage-addon.sh disable jupyter --project-dir /workspaces/myproject

set -euo pipefail

# Determine script directory and toolkit root
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TOOLKIT_DIR="$(dirname "${SCRIPT_DIR}")"
ADDONS_TEMPLATE_DIR="${TOOLKIT_DIR}/templates/mcp-addons"

# Source common utilities
if [ -f "${SCRIPT_DIR}/common.sh" ]; then
    source "${SCRIPT_DIR}/common.sh"
elif [ -f "${TOOLKIT_DIR}/scripts/common.sh" ]; then
    source "${TOOLKIT_DIR}/scripts/common.sh"
else
    echo "Warning: common.sh not found."
    log_info()  { echo -e "[INFO] $*"; }
    log_ok()    { echo -e "[OK] $*"; }
    log_warn()  { echo -e "[WARN] $*"; }
    log_error() { echo -e "[ERROR] $*"; }

    load_env() {
        local project_dir="${1:-$PWD}"
        if [ -f "${project_dir}/.env" ]; then set -a; source "${project_dir}/.env"; set +a; fi
        if [ -f "${project_dir}/.devcontainer/.env" ]; then set -a; source "${project_dir}/.devcontainer/.env"; set +a; fi
    }
fi

# ============================================================================
# Helper Functions
# ============================================================================

# Get path to active-addons.json for a project
get_addons_state_path() {
    local project_dir="$1"
    echo "${project_dir}/.claude/active-addons.json"
}

# Read active-addons.json, return empty object if missing
read_addons_state() {
    local state_path
    state_path="$(get_addons_state_path "$1")"
    if [ -f "${state_path}" ]; then
        cat "${state_path}"
    else
        echo '{"addons":{}}'
    fi
}

# Write active-addons.json
write_addons_state() {
    local project_dir="$1"
    local content="$2"
    local state_path
    state_path="$(get_addons_state_path "${project_dir}")"
    mkdir -p "$(dirname "${state_path}")"
    echo "${content}" > "${state_path}"
}

# Check if an addon is enabled in state
is_addon_enabled() {
    local project_dir="$1"
    local addon_name="$2"
    local state
    state="$(read_addons_state "${project_dir}")"
    python3 -c "
import json, sys
state = json.loads('''${state}''')
addon = state.get('addons', {}).get('${addon_name}', {})
sys.exit(0 if addon.get('enabled', False) else 1)
" 2>/dev/null
}

# List available addon templates
list_available_addons() {
    local addons=()
    if [ -d "${ADDONS_TEMPLATE_DIR}" ]; then
        for f in "${ADDONS_TEMPLATE_DIR}"/*.addon.json; do
            [ -f "$f" ] || continue
            addons+=("$(basename "$f" .addon.json)")
        done
    fi
    echo "${addons[@]}"
}

# Read addon template metadata
read_addon_meta() {
    local addon_name="$1"
    local template="${ADDONS_TEMPLATE_DIR}/${addon_name}.addon.json"
    if [ ! -f "${template}" ]; then
        return 1
    fi
    python3 -c "
import json
with open('${template}') as f:
    data = json.load(f)
meta = data.get('_meta', {})
print(json.dumps(meta))
"
}

# Merge active addons into .mcp.json
# This is the core merge function, also called by switch-mcp-profile.sh
merge_addons_into_mcp() {
    local project_dir="$1"
    local mcp_json="${project_dir}/.mcp.json"
    local state_path
    state_path="$(get_addons_state_path "${project_dir}")"

    # If no .mcp.json, nothing to merge into
    if [ ! -f "${mcp_json}" ]; then
        log_warn "No .mcp.json found at ${mcp_json} — skipping addon merge"
        return 0
    fi

    # If no active-addons.json, nothing to merge
    if [ ! -f "${state_path}" ]; then
        return 0
    fi

    # Load environment variables for substitution
    load_env "${project_dir}"

    python3 << 'PYEOF'
import json
import os
import sys

project_dir = os.environ.get('PROJECT_DIR', '')
toolkit_dir = os.environ.get('TOOLKIT_DIR', '')
mcp_json_path = os.path.join(project_dir, '.mcp.json')
state_path = os.path.join(project_dir, '.claude', 'active-addons.json')
addons_template_dir = os.path.join(toolkit_dir, 'templates', 'mcp-addons')

# Read current .mcp.json
try:
    with open(mcp_json_path) as f:
        mcp_config = json.load(f)
except (FileNotFoundError, json.JSONDecodeError) as e:
    print(f"Warning: Could not read .mcp.json: {e}", file=sys.stderr)
    sys.exit(0)

# Read active-addons.json
try:
    with open(state_path) as f:
        state = json.load(f)
except (FileNotFoundError, json.JSONDecodeError):
    # No addons configured — nothing to do
    sys.exit(0)

addons = state.get('addons', {})
merged_count = 0

for addon_name, addon_state in addons.items():
    if not addon_state.get('enabled', False):
        continue

    # Load addon template
    template_path = os.path.join(addons_template_dir, f'{addon_name}.addon.json')
    if not os.path.exists(template_path):
        print(f"Warning: Addon template not found: {template_path}", file=sys.stderr)
        continue

    try:
        with open(template_path) as f:
            addon_template = json.load(f)
    except json.JSONDecodeError as e:
        print(f"Warning: Invalid addon template {addon_name}: {e}", file=sys.stderr)
        continue

    # Merge mcpServers entries
    addon_servers = addon_template.get('mcpServers', {})
    if 'mcpServers' not in mcp_config:
        mcp_config['mcpServers'] = {}

    for server_name, server_config in addon_servers.items():
        mcp_config['mcpServers'][server_name] = server_config
        merged_count += 1

# Substitute environment variables in the merged result
content = json.dumps(mcp_config, indent=2)

# Replace ${VAR} patterns with environment variable values
import re
def env_replace(match):
    var_name = match.group(1)
    return os.environ.get(var_name, '')

content = re.sub(r'\$\{([^}]+)\}', env_replace, content)

# Write back
with open(mcp_json_path, 'w') as f:
    f.write(content)

if merged_count > 0:
    print(f"Merged {merged_count} addon server(s) into .mcp.json")
PYEOF
}

# Update settings.local.json to include addon servers and tool permissions
update_settings_local() {
    local project_dir="$1"
    local mcp_json="${project_dir}/.mcp.json"
    local settings_path="${project_dir}/.claude/settings.local.json"

    if [ ! -f "${mcp_json}" ]; then
        return 0
    fi

    if [ ! -f "${settings_path}" ]; then
        return 0
    fi

    python3 << 'PYEOF'
import json
import os
import sys

project_dir = os.environ.get('PROJECT_DIR', '')
toolkit_dir = os.environ.get('TOOLKIT_DIR', '')
mcp_json_path = os.path.join(project_dir, '.mcp.json')
settings_path = os.path.join(project_dir, '.claude', 'settings.local.json')
state_path = os.path.join(project_dir, '.claude', 'active-addons.json')
addons_template_dir = os.path.join(toolkit_dir, 'templates', 'mcp-addons')

# Read current server list from .mcp.json
try:
    with open(mcp_json_path) as f:
        mcp_config = json.load(f)
    servers = list(mcp_config.get('mcpServers', {}).keys())
except Exception:
    sys.exit(0)

# Read current settings.local.json
try:
    with open(settings_path) as f:
        settings = json.load(f)
except Exception:
    sys.exit(0)

# Update enabledMcpjsonServers
settings['enabledMcpjsonServers'] = servers

# Collect tool permissions from all enabled addons
addon_permissions = []
try:
    with open(state_path) as f:
        state = json.load(f)
    for addon_name, addon_state in state.get('addons', {}).items():
        if not addon_state.get('enabled', False):
            continue
        template_path = os.path.join(addons_template_dir, f'{addon_name}.addon.json')
        if not os.path.exists(template_path):
            continue
        with open(template_path) as f:
            addon_template = json.load(f)
        perms = addon_template.get('_meta', {}).get('tool_permissions', [])
        addon_permissions.extend(perms)
except (FileNotFoundError, json.JSONDecodeError):
    pass  # No active-addons.json or parse error — skip

# Merge addon tool permissions into permissions.allow (deduplicated)
if 'permissions' not in settings:
    settings['permissions'] = {'allow': [], 'deny': []}
if 'allow' not in settings['permissions']:
    settings['permissions']['allow'] = []

allow_list = settings['permissions']['allow']

# Remove any stale mcp__*__* entries that are NOT in the current addon_permissions set
# (This handles disable: old perms get cleaned, new ones get added)
non_addon_entries = [e for e in allow_list if not e.startswith('mcp__')]
# Keep mcp__ entries that belong to currently enabled addons
allow_list = non_addon_entries + addon_permissions

# Deduplicate while preserving order
seen = set()
deduped = []
for entry in allow_list:
    if entry not in seen:
        seen.add(entry)
        deduped.append(entry)

settings['permissions']['allow'] = deduped

# Write back
with open(settings_path, 'w') as f:
    json.dump(settings, f, indent=2)

if addon_permissions:
    print(f"Merged {len(addon_permissions)} tool permission(s) into settings.local.json")
PYEOF
}

# ============================================================================
# Commands
# ============================================================================

cmd_list() {
    local project_dir="${1:-$PWD}"

    echo "Available addons:"
    echo ""

    local available
    available="$(list_available_addons)"

    if [ -z "${available}" ]; then
        echo "  (none found in ${ADDONS_TEMPLATE_DIR})"
        return 0
    fi

    for addon_name in ${available}; do
        local status_marker="disabled"
        if is_addon_enabled "${project_dir}" "${addon_name}"; then
            status_marker="ENABLED"
        fi

        # Get description from metadata
        local description=""
        local meta
        if meta="$(read_addon_meta "${addon_name}" 2>/dev/null)"; then
            description="$(python3 -c "import json; print(json.loads('''${meta}''').get('description', ''))" 2>/dev/null || true)"
        fi

        if [ "${status_marker}" = "ENABLED" ]; then
            echo "  * ${addon_name} [${status_marker}]"
        else
            echo "    ${addon_name} [${status_marker}]"
        fi
        if [ -n "${description}" ]; then
            echo "      ${description}"
        fi
    done
    echo ""
}

cmd_enable() {
    local addon_name="$1"
    local project_dir="$2"

    # Validate addon template exists
    local template="${ADDONS_TEMPLATE_DIR}/${addon_name}.addon.json"
    if [ ! -f "${template}" ]; then
        log_error "Addon template not found: ${addon_name}"
        log_info "Available addons:"
        list_available_addons | tr ' ' '\n' | sed 's/^/  - /'
        exit 1
    fi

    # Check if requires_setup and warn
    local meta
    if meta="$(read_addon_meta "${addon_name}" 2>/dev/null)"; then
        local requires_setup
        requires_setup="$(python3 -c "import json; print(json.loads('''${meta}''').get('requires_setup', False))" 2>/dev/null || echo "False")"
        if [ "${requires_setup}" = "True" ]; then
            local setup_script
            setup_script="$(python3 -c "import json; print(json.loads('''${meta}''').get('setup_script', ''))" 2>/dev/null || echo "")"
            log_warn "Addon '${addon_name}' requires setup before use."
            if [ -n "${setup_script}" ]; then
                log_info "Run: ${SCRIPT_DIR}/mcp_servers/${setup_script}"
            fi

            # Check env vars
            local env_vars
            env_vars="$(python3 -c "
import json
meta = json.loads('''${meta}''')
vars = meta.get('env_vars', [])
print(' '.join(vars))
" 2>/dev/null || echo "")"

            if [ -n "${env_vars}" ]; then
                local missing_vars=()
                for var in ${env_vars}; do
                    if [ -z "${!var:-}" ]; then
                        missing_vars+=("${var}")
                    fi
                done
                if [ ${#missing_vars[@]} -gt 0 ]; then
                    log_warn "Missing environment variables: ${missing_vars[*]}"
                    log_info "Set these in .env or .devcontainer/.env"
                fi
            fi
        fi
    fi

    # Update active-addons.json
    local state
    state="$(read_addons_state "${project_dir}")"
    local new_state
    new_state="$(python3 -c "
import json
state = json.loads('''${state}''')
if 'addons' not in state:
    state['addons'] = {}
state['addons']['${addon_name}'] = {'enabled': True}
print(json.dumps(state, indent=2))
")"
    write_addons_state "${project_dir}" "${new_state}"

    log_ok "Addon '${addon_name}' enabled in active-addons.json"

    # Merge into .mcp.json
    export PROJECT_DIR="${project_dir}"
    export TOOLKIT_DIR="${TOOLKIT_DIR}"
    merge_addons_into_mcp "${project_dir}"

    # Update settings.local.json
    update_settings_local "${project_dir}"

    log_ok "Addon '${addon_name}' merged into .mcp.json"
    log_info "Restart Claude Code to apply: exit && claude"
}

cmd_disable() {
    local addon_name="$1"
    local project_dir="$2"

    # Update active-addons.json
    local state
    state="$(read_addons_state "${project_dir}")"
    local new_state
    new_state="$(python3 -c "
import json
state = json.loads('''${state}''')
if 'addons' not in state:
    state['addons'] = {}
state['addons']['${addon_name}'] = {'enabled': False}
print(json.dumps(state, indent=2))
")"
    write_addons_state "${project_dir}" "${new_state}"

    log_ok "Addon '${addon_name}' disabled in active-addons.json"

    # Remove addon servers from .mcp.json
    local mcp_json="${project_dir}/.mcp.json"
    local template="${ADDONS_TEMPLATE_DIR}/${addon_name}.addon.json"

    if [ -f "${mcp_json}" ] && [ -f "${template}" ]; then
        python3 -c "
import json

with open('${mcp_json}') as f:
    mcp_config = json.load(f)

with open('${template}') as f:
    addon = json.load(f)

addon_servers = addon.get('mcpServers', {})
for server_name in addon_servers:
    if server_name in mcp_config.get('mcpServers', {}):
        del mcp_config['mcpServers'][server_name]

with open('${mcp_json}', 'w') as f:
    json.dump(mcp_config, f, indent=2)
"
        log_ok "Removed addon servers from .mcp.json"
    fi

    # Update settings.local.json
    export PROJECT_DIR="${project_dir}"
    update_settings_local "${project_dir}"

    log_info "Restart Claude Code to apply: exit && claude"
}

cmd_status() {
    local addon_name="$1"
    local project_dir="$2"

    # Check template exists
    local template="${ADDONS_TEMPLATE_DIR}/${addon_name}.addon.json"
    if [ ! -f "${template}" ]; then
        log_error "Unknown addon: ${addon_name}"
        exit 1
    fi

    echo "Addon: ${addon_name}"
    echo ""

    # Enabled status
    if is_addon_enabled "${project_dir}" "${addon_name}"; then
        echo "  Status: ENABLED"
    else
        echo "  Status: disabled"
    fi

    # Metadata
    local meta
    if meta="$(read_addon_meta "${addon_name}" 2>/dev/null)"; then
        local description
        description="$(python3 -c "import json; print(json.loads('''${meta}''').get('description', 'N/A'))" 2>/dev/null || echo "N/A")"
        local tokens
        tokens="$(python3 -c "import json; print(json.loads('''${meta}''').get('context_tokens_estimate', 'unknown'))" 2>/dev/null || echo "unknown")"
        local runtime_dep
        runtime_dep="$(python3 -c "import json; print(json.loads('''${meta}''').get('runtime_dependency', 'none'))" 2>/dev/null || echo "none")"

        echo "  Description: ${description}"
        echo "  Context cost: ${tokens}"
        echo "  Runtime dependency: ${runtime_dep}"
    fi

    # Check if server appears in .mcp.json
    local mcp_json="${project_dir}/.mcp.json"
    if [ -f "${mcp_json}" ]; then
        local in_mcp
        in_mcp="$(python3 -c "
import json
with open('${mcp_json}') as f:
    config = json.load(f)
servers = list(config.get('mcpServers', {}).keys())

with open('${template}') as f:
    addon = json.load(f)
addon_servers = list(addon.get('mcpServers', {}).keys())

for s in addon_servers:
    if s in servers:
        print(f'  In .mcp.json: YES ({s})')
    else:
        print(f'  In .mcp.json: NO ({s})')
" 2>/dev/null || echo "  In .mcp.json: unknown")"
        echo "${in_mcp}"
    fi

    # Check env vars
    load_env "${project_dir}"
    local env_vars
    if meta="$(read_addon_meta "${addon_name}" 2>/dev/null)"; then
        env_vars="$(python3 -c "
import json
meta = json.loads('''${meta}''')
for var in meta.get('env_vars', []):
    import os
    val = os.environ.get(var, '')
    status = 'SET' if val else 'MISSING'
    print(f'  {var}: {status}')
" 2>/dev/null || true)"
        if [ -n "${env_vars}" ]; then
            echo ""
            echo "  Environment variables:"
            echo "${env_vars}"
        fi
    fi

    echo ""
}

# ============================================================================
# Main
# ============================================================================

COMMAND="${1:-}"
ADDON_NAME=""
PROJECT_DIR="${PWD}"

if [ -z "${COMMAND}" ]; then
    echo "Usage: $0 <command> [addon-name] [--project-dir DIR]"
    echo ""
    echo "Commands:"
    echo "  list                  Show available + active addons"
    echo "  enable <name>         Enable addon, merge into .mcp.json"
    echo "  disable <name>        Disable addon, remove from .mcp.json"
    echo "  status <name>         Check addon status and dependencies"
    echo ""
    echo "Options:"
    echo "  --project-dir DIR     Project directory (default: current directory)"
    exit 0
fi

shift

# Parse remaining arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --project-dir)
            PROJECT_DIR="$2"
            shift 2
            ;;
        -*)
            log_error "Unknown option: $1"
            exit 1
            ;;
        *)
            if [ -z "${ADDON_NAME}" ]; then
                ADDON_NAME="$1"
            fi
            shift
            ;;
    esac
done

# Dispatch command
case "${COMMAND}" in
    list)
        cmd_list "${PROJECT_DIR}"
        ;;
    enable)
        if [ -z "${ADDON_NAME}" ]; then
            log_error "Addon name required. Usage: $0 enable <name>"
            exit 1
        fi
        cmd_enable "${ADDON_NAME}" "${PROJECT_DIR}"
        ;;
    disable)
        if [ -z "${ADDON_NAME}" ]; then
            log_error "Addon name required. Usage: $0 disable <name>"
            exit 1
        fi
        cmd_disable "${ADDON_NAME}" "${PROJECT_DIR}"
        ;;
    status)
        if [ -z "${ADDON_NAME}" ]; then
            log_error "Addon name required. Usage: $0 status <name>"
            exit 1
        fi
        cmd_status "${ADDON_NAME}" "${PROJECT_DIR}"
        ;;
    *)
        log_error "Unknown command: ${COMMAND}"
        echo "Valid commands: list, enable, disable, status"
        exit 1
        ;;
esac
