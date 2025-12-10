#!/usr/bin/env bash
#
# toggle-mcp.sh - Simple toggle for MCP servers
#
# Usage:
#   ./toggle-mcp.sh [list|enable|disable] [server_name]
#
# Examples:
#   ./toggle-mcp.sh list
#   ./toggle-mcp.sh enable pal
#   ./toggle-mcp.sh disable tooluniverse
#

set -euo pipefail

# Colors
GREEN='\033[0;32m'
RED='\033[0;31m'
NC='\033[0m'

CMD="${1:-list}"
SERVER="${2:-}"

case "$CMD" in
    list)
        if command -v claude &>/dev/null; then
            echo "Current MCP Servers:"
            claude mcp list
        else
            echo -e "${RED}Error: claude CLI not found.${NC}"
            exit 1
        fi
        ;;
    enable)
        if [ -z "$SERVER" ]; then
            echo "Usage: $0 enable <server_name>"
            exit 1
        fi
        echo "Enabling $SERVER..."
        # Note: claude mcp enable currently might not exist in all versions, 
        # but if it does or if we map it to settings update:
        # Currently Claude Code CLI doesn't have a direct 'enable' command for toggling 
        # without removing/adding. 
        # Instead, we instruct user to use /mcp inside Claude, OR we modify settings.local.json.
        
        SETTINGS_FILE=".claude/settings.local.json"
        if [ -f "$SETTINGS_FILE" ]; then
             python3 -c "import json; f='$SETTINGS_FILE'; d=json.load(open(f)); s='$SERVER'; 
d['disabledMcpjsonServers'] = [x for x in d.get('disabledMcpjsonServers', []) if x != s];
if 'enabledMcpjsonServers' not in d: d['enabledMcpjsonServers'] = []
if s not in d['enabledMcpjsonServers']: d['enabledMcpjsonServers'].append(s);
json.dump(d, open(f,'w'), indent=2); print(f'Enabled {s} in {f}')"
        else
            echo -e "${RED}Settings file not found: $SETTINGS_FILE${NC}"
        fi
        ;;
    disable)
        if [ -z "$SERVER" ]; then
            echo "Usage: $0 disable <server_name>"
            exit 1
        fi
        echo "Disabling $SERVER..."
        SETTINGS_FILE=".claude/settings.local.json"
        if [ -f "$SETTINGS_FILE" ]; then
             python3 -c "import json; f='$SETTINGS_FILE'; d=json.load(open(f)); s='$SERVER'; 
if 'disabledMcpjsonServers' not in d: d['disabledMcpjsonServers'] = []
if s not in d['disabledMcpjsonServers']: d['disabledMcpjsonServers'].append(s);
d['enabledMcpjsonServers'] = [x for x in d.get('enabledMcpjsonServers', []) if x != s];
json.dump(d, open(f,'w'), indent=2); print(f'Disabled {s} in {f}')"
        else
            echo -e "${RED}Settings file not found: $SETTINGS_FILE${NC}"
        fi
        ;;
    *)
        echo "Unknown command: $CMD"
        echo "Usage: $0 [list|enable|disable] [server_name]"
        exit 1
        ;;
esac
