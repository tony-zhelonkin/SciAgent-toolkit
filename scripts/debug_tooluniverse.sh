#!/usr/bin/env bash
#
# Debug script for ToolUniverse MCP server issues
#
# Run this script to diagnose why ToolUniverse shows "connecting..." indefinitely
# in Claude Code. It tests the Python import, command detection, and server startup.
#
# Usage:
#   ./debug_tooluniverse.sh
#   ./debug_tooluniverse.sh --verbose
#

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

VERBOSE=false
if [[ "${1:-}" == "--verbose" ]]; then
    VERBOSE=true
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Search for tooluniverse-env in common locations
TOOLUNIVERSE_ENV=""
for search_path in \
    "${SCRIPT_DIR}/tooluniverse-env" \
    "${SCRIPT_DIR}/../tooluniverse-env" \
    "${PWD}/tooluniverse-env" \
    "${PWD}/01_scripts/SciAgent-toolkit/scripts/tooluniverse-env"; do
    if [ -d "$search_path" ]; then
        TOOLUNIVERSE_ENV="$(cd "$search_path" && pwd)"
        break
    fi
done

echo ""
echo "=============================================="
echo "   ToolUniverse MCP Server Debug Script"
echo "=============================================="
echo ""

# Check if tooluniverse-env found
if [ -z "${TOOLUNIVERSE_ENV}" ]; then
    log_error "tooluniverse-env directory not found!"
    log_info "Searched in:"
    log_info "  - ${SCRIPT_DIR}/tooluniverse-env"
    log_info "  - ${SCRIPT_DIR}/../tooluniverse-env"
    log_info "  - ${PWD}/tooluniverse-env"
    log_info ""
    log_info "Run setup_tooluniverse.sh first to install ToolUniverse."
    exit 1
fi

log_ok "Environment found: ${TOOLUNIVERSE_ENV}"
echo ""

# Check uv
echo "=== 1. Checking uv package manager ==="
if command -v uv &>/dev/null; then
    UV_VERSION=$(uv --version 2>/dev/null || echo "unknown")
    log_ok "uv installed: ${UV_VERSION}"
else
    log_error "uv not found! Install with: curl -LsSf https://astral.sh/uv/install.sh | sh"
    exit 1
fi
echo ""

# Test Python import
echo "=== 2. Testing Python import ==="
if uv --directory "${TOOLUNIVERSE_ENV}" run python -c "
import sys
print(f'Python: {sys.version}')
import tooluniverse
print(f'tooluniverse version: {getattr(tooluniverse, \"__version__\", \"unknown\")}')
print('Import successful!')
" 2>&1; then
    log_ok "Python import works"
else
    log_error "Python import failed!"
    log_info "Try reinstalling: uv --directory ${TOOLUNIVERSE_ENV} pip install --force-reinstall tooluniverse"
    exit 1
fi
echo ""

# Detect command name
echo "=== 3. Detecting MCP command ==="
TOOLUNIVERSE_CMD=""
for cmd in "tooluniverse-mcp" "python -m tooluniverse.smcp_server"; do
    if uv --directory "${TOOLUNIVERSE_ENV}" run ${cmd} --help &>/dev/null 2>&1; then
        TOOLUNIVERSE_CMD="${cmd}"
        log_ok "Command works: ${cmd}"
        break
    else
        log_warn "Command not available: ${cmd}"
    fi
done

if [ -z "${TOOLUNIVERSE_CMD}" ]; then
    log_error "No working MCP command found!"
    exit 1
fi
echo ""

# Test help output
echo "=== 4. Testing --help output ==="
if timeout 10 uv --directory "${TOOLUNIVERSE_ENV}" run ${TOOLUNIVERSE_CMD} --help 2>&1; then
    log_ok "--help works"
else
    log_error "--help failed or timed out"
fi
echo ""

# Test MCP startup with timeout
echo "=== 5. Testing MCP server startup (5s timeout) ==="
log_info "Starting server... (will timeout after 5s - this is expected)"
echo ""

if [ "${VERBOSE}" = true ]; then
    # Verbose mode - show all output
    # Note: --exclude-tool-types removed in ToolUniverse 1.0.14+
    timeout 5 uv --directory "${TOOLUNIVERSE_ENV}" run ${TOOLUNIVERSE_CMD} 2>&1 || true
else
    # Normal mode - capture and summarize
    OUTPUT=$(timeout 5 uv --directory "${TOOLUNIVERSE_ENV}" run ${TOOLUNIVERSE_CMD} 2>&1 || true)

    if echo "$OUTPUT" | grep -qi "error\|exception\|traceback"; then
        log_error "Server produced errors:"
        echo "$OUTPUT" | tail -20
    elif [ -z "$OUTPUT" ]; then
        log_warn "Server produced no output (may be waiting for stdin - normal for MCP)"
    else
        log_info "Server output (first 10 lines):"
        echo "$OUTPUT" | head -10
    fi
fi

echo ""
log_ok "Timeout reached (expected - MCP servers wait for stdio input)"
echo ""

# Check for common issues
echo "=== 6. Common Issue Checks ==="

# Check if there are multiple Python versions
PYTHON_COUNT=$(find "${TOOLUNIVERSE_ENV}" -name "python*" -type f 2>/dev/null | wc -l)
if [ "$PYTHON_COUNT" -gt 2 ]; then
    log_warn "Multiple Python executables found in venv (may cause conflicts)"
fi

# Check for missing dependencies
if uv --directory "${TOOLUNIVERSE_ENV}" run python -c "
try:
    import mcp
    print('mcp: OK')
except ImportError as e:
    print(f'mcp: MISSING - {e}')

try:
    import httpx
    print('httpx: OK')
except ImportError as e:
    print(f'httpx: MISSING - {e}')
" 2>&1; then
    log_ok "Core dependencies present"
fi
echo ""

# Summary
echo "=============================================="
echo "   Debug Summary"
echo "=============================================="
echo ""
log_info "Environment: ${TOOLUNIVERSE_ENV}"
log_info "Command: uv --directory ${TOOLUNIVERSE_ENV} run ${TOOLUNIVERSE_CMD}"
echo ""
log_info "If ToolUniverse still shows 'connecting...' in Claude Code:"
log_info "  1. Check Claude logs: ls -la ~/.cache/claude-cli-nodejs/"
log_info "  2. Run with PYTHONUNBUFFERED=1 for immediate output"
log_info "  3. Try: claude mcp remove tooluniverse && claude mcp add ..."
echo ""
log_info "To test MCP JSON-RPC handshake manually:"
log_info "  echo '{\"jsonrpc\":\"2.0\",\"method\":\"initialize\",\"id\":1}' | \\"
log_info "    uv --directory ${TOOLUNIVERSE_ENV} run ${TOOLUNIVERSE_CMD}"
echo ""
