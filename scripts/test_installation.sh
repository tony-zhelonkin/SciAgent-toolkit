#!/usr/bin/env bash
#
# Installation Testing Script
#
# Runs comprehensive tests on a completed installation
#

set -euo pipefail

# Get script directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Source common utilities
if [ -f "${SCRIPT_DIR}/common.sh" ]; then
    source "${SCRIPT_DIR}/common.sh"
else
    echo "Error: common.sh not found in ${SCRIPT_DIR}"
    exit 1
fi

PASSED=0
FAILED=0
WARNINGS=0

test_pass() {
    echo -e "${GREEN}✓${NC} $1"
    PASSED=$((PASSED + 1))
}

test_fail() {
    echo -e "${RED}✗${NC} $1"
    FAILED=$((FAILED + 1))
}

test_warn() {
    echo -e "${YELLOW}⚠${NC} $1"
    WARNINGS=$((WARNINGS + 1))
}

echo "=========================================="
echo "SciAgent Toolkit Installation Test"
echo "=========================================="
echo ""

# Test 1: MCP Configuration File
echo "Test 1: MCP Configuration"
if [ -f ".mcp.json" ]; then
    test_pass "MCP config file exists"

    if python3 -m json.tool .mcp.json > /dev/null 2>&1; then
        test_pass "MCP config is valid JSON"
    else
        test_fail "MCP config has invalid JSON"
    fi

    if command -v jq &>/dev/null; then
        SERVER_COUNT=$(jq '.mcpServers | length' .mcp.json)
        if [ "$SERVER_COUNT" -ge 3 ]; then
            test_pass "Expected number of servers configured ($SERVER_COUNT)"
        else
            test_warn "Only $SERVER_COUNT servers configured (expected 3+)"
        fi
    fi
else
    test_fail "MCP config file not found"
fi

# Test 2: ToolUniverse
echo ""
echo "Test 2: ToolUniverse"
if [ -d "tooluniverse-env" ]; then
    test_pass "ToolUniverse environment exists"

    if uv --directory ./tooluniverse-env run python -c "import tooluniverse" 2>/dev/null; then
        VERSION=$(uv --directory ./tooluniverse-env run python -c "import tooluniverse; print(tooluniverse.__version__)" 2>/dev/null)
        test_pass "ToolUniverse imports successfully (v$VERSION)"
    else
        test_fail "ToolUniverse import failed"
    fi

    if uv --directory ./tooluniverse-env run tooluniverse-mcp --help > /dev/null 2>&1; then
        test_pass "ToolUniverse MCP command works (tooluniverse-mcp)"
    elif uv --directory ./tooluniverse-env run tooluniverse-smcp-stdio --help > /dev/null 2>&1; then
        test_pass "ToolUniverse MCP command works (tooluniverse-smcp-stdio)"
    else
        test_fail "ToolUniverse MCP command not working"
    fi
else
    test_fail "ToolUniverse environment not found"
fi

# Test 3: Sequential Thinking
echo ""
echo "Test 3: Sequential Thinking"
if command -v npx &>/dev/null; then
    test_pass "npx is available"

    if timeout 10 npx -y @modelcontextprotocol/server-sequential-thinking --help > /dev/null 2>&1; then
        test_pass "Sequential Thinking server works"
    else
        test_warn "Sequential Thinking server test timeout (may need download)"
    fi
else
    test_fail "npx not found (Node.js not installed)"
fi

# Test 4: Serena
echo ""
echo "Test 4: Serena"
if command -v uvx &>/dev/null; then
    test_pass "uvx is available"

    if timeout 30 uvx --from git+https://github.com/oraios/serena serena --help > /dev/null 2>&1; then
        test_pass "Serena server works"
    else
        test_warn "Serena server test failed (may need first-time download)"
    fi
else
    test_fail "uvx not found (uv not installed)"
fi

# Test 5: PAL
echo ""
echo "Test 5: PAL"
if command -v uvx &>/dev/null; then
    if timeout 30 uvx --from git+https://github.com/BeehiveInnovations/pal-mcp-server.git pal-mcp-server --help > /dev/null 2>&1; then
        test_pass "PAL server works"
    else
        test_warn "PAL server test failed (may need first-time download)"
    fi
else
    test_warn "uvx not found, skipping PAL test"
fi

# Test 6: Claude Code
echo ""
echo "Test 6: Claude Code"
if command -v claude &>/dev/null; then
    test_pass "Claude Code is installed"

    if claude --version > /dev/null 2>&1; then
        VERSION=$(claude --version 2>/dev/null | head -1)
        test_pass "Claude Code version check: $VERSION"
    else
        test_warn "Claude Code version check failed"
    fi
else
    test_warn "Claude Code not installed (may be intentional)"
fi

# Summary
echo ""
echo "=========================================="
echo "Test Summary"
echo "=========================================="
echo -e "Passed:   ${GREEN}$PASSED${NC}"
echo -e "Failed:   ${RED}$FAILED${NC}"
echo -e "Warnings: ${YELLOW}$WARNINGS${NC}"
echo ""

if [ $FAILED -eq 0 ]; then
    echo -e "${GREEN}All critical tests passed!${NC}"
    echo ""
    echo "Next steps:"
    echo "  1. Start Claude Code: claude"
    echo "  2. Check MCP servers: /mcp"
    echo "  3. Try a query: 'What scientific tools are available?'"
    exit 0
else
    echo -e "${RED}Some tests failed. Please review errors above.${NC}"
    echo ""
    echo "For troubleshooting, see: docs/TROUBLESHOOTING.md"
    exit 1
fi
