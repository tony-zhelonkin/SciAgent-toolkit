#!/usr/bin/env bash
#
# Comprehensive testing script for ToolUniverse Docker integration
#
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"

# Colors
GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

log_info()  { echo -e "${YELLOW}[INFO]${NC} $*"; }
log_ok()    { echo -e "${GREEN}[OK]${NC} $*"; }
log_error() { echo -e "${RED}[ERROR]${NC} $*"; }
log_step()  { echo -e "${BLUE}[STEP]${NC} $*"; }

echo "==================================================================="
echo "           ToolUniverse Docker Integration Test Suite"
echo "==================================================================="
echo ""
echo "Project root: $PROJECT_ROOT"
echo "Test directory: $SCRIPT_DIR"
echo ""

cd "$PROJECT_ROOT"

TOTAL_TESTS=0
PASSED_TESTS=0
FAILED_TESTS=0

# Test 0: Architecture tests (role system, templates, profile switching)
echo ""
log_step "Test 0: Building architecture test image..."
TOTAL_TESTS=$((TOTAL_TESTS + 1))

if docker build -f docker/test/Dockerfile.architecture-test \
     -t architecture-test:latest . 2>&1 | tee /tmp/architecture-test-build.log; then
    log_ok "Architecture tests passed"
    PASSED_TESTS=$((PASSED_TESTS + 1))
else
    log_error "Architecture test build failed. See /tmp/architecture-test-build.log"
    FAILED_TESTS=$((FAILED_TESTS + 1))
    echo "==================================================================="
    echo "Architecture tests failed. Review the changes before proceeding."
    echo "==================================================================="
fi

# Test 1: Base ToolUniverse installation
echo ""
log_step "Test 1: Building base ToolUniverse test image..."
TOTAL_TESTS=$((TOTAL_TESTS + 1))

if docker build -f docker/test/Dockerfile.tooluniverse-test \
     -t tooluniverse-test:latest . 2>&1 | tee /tmp/tooluniverse-test-build.log; then
    log_ok "Base image built successfully"
    PASSED_TESTS=$((PASSED_TESTS + 1))
else
    log_error "Base image build failed. See /tmp/tooluniverse-test-build.log"
    FAILED_TESTS=$((FAILED_TESTS + 1))
    echo "==================================================================="
    echo "Build failed. Cannot proceed with further tests."
    echo "==================================================================="
    exit 1
fi

log_info "Testing ToolUniverse MCP server --help..."
TOTAL_TESTS=$((TOTAL_TESTS + 1))

if docker run --rm tooluniverse-test:latest \
     bash -c "uv --directory /opt/sciagent/scripts/tooluniverse-env run tooluniverse-smcp-stdio --help" > /dev/null 2>&1; then
    log_ok "ToolUniverse MCP server responds to --help"
    PASSED_TESTS=$((PASSED_TESTS + 1))
else
    log_error "ToolUniverse MCP server test failed"
    FAILED_TESTS=$((FAILED_TESTS + 1))
fi

# Test 2: Claude Code integration
echo ""
log_step "Test 2: Building Claude Code MCP test image..."
TOTAL_TESTS=$((TOTAL_TESTS + 1))

if docker build -f docker/test/Dockerfile.claude-test \
     -t claude-mcp-test:latest . 2>&1 | tee /tmp/claude-test-build.log; then
    log_ok "Claude test image built successfully"
    PASSED_TESTS=$((PASSED_TESTS + 1))
else
    log_error "Claude test image build failed. See /tmp/claude-test-build.log"
    FAILED_TESTS=$((FAILED_TESTS + 1))
fi

log_info "Testing Claude Code installation..."
TOTAL_TESTS=$((TOTAL_TESTS + 1))

if docker run --rm claude-mcp-test:latest claude --version > /dev/null 2>&1; then
    log_ok "Claude Code is installed"
    PASSED_TESTS=$((PASSED_TESTS + 1))
else
    log_error "Claude Code test failed"
    FAILED_TESTS=$((FAILED_TESTS + 1))
fi

log_info "Testing Claude MCP config validity..."
TOTAL_TESTS=$((TOTAL_TESTS + 1))

if docker run --rm claude-mcp-test:latest \
     python3 -m json.tool /root/.config/claude/mcp-config.json > /dev/null 2>&1; then
    log_ok "Claude MCP config is valid JSON"
    PASSED_TESTS=$((PASSED_TESTS + 1))
else
    log_error "Claude MCP config validation failed"
    FAILED_TESTS=$((FAILED_TESTS + 1))
fi

# Test 3: Codex CLI integration
echo ""
log_step "Test 3: Building Codex CLI MCP test image..."
TOTAL_TESTS=$((TOTAL_TESTS + 1))

if docker build -f docker/test/Dockerfile.codex-test \
     -t codex-mcp-test:latest . 2>&1 | tee /tmp/codex-test-build.log; then
    log_ok "Codex test image built successfully"
    PASSED_TESTS=$((PASSED_TESTS + 1))
else
    log_error "Codex test image build failed. See /tmp/codex-test-build.log"
    FAILED_TESTS=$((FAILED_TESTS + 1))
fi

log_info "Testing Codex CLI installation..."
TOTAL_TESTS=$((TOTAL_TESTS + 1))

if docker run --rm codex-mcp-test:latest codex --version > /dev/null 2>&1; then
    log_ok "Codex CLI is installed"
    PASSED_TESTS=$((PASSED_TESTS + 1))
else
    log_error "Codex CLI test failed"
    FAILED_TESTS=$((FAILED_TESTS + 1))
fi

log_info "Testing Codex config validity..."
TOTAL_TESTS=$((TOTAL_TESTS + 1))

if docker run --rm codex-mcp-test:latest \
     python3 -c "import tomli; tomli.load(open('/root/.codex/config.toml', 'rb'))" > /dev/null 2>&1; then
    log_ok "Codex config is valid TOML"
    PASSED_TESTS=$((PASSED_TESTS + 1))
else
    log_error "Codex config validation failed"
    FAILED_TESTS=$((FAILED_TESTS + 1))
fi

# Test 4: Gemini CLI integration
echo ""
log_step "Test 4: Building Gemini CLI MCP test image..."
TOTAL_TESTS=$((TOTAL_TESTS + 1))

if docker build -f docker/test/Dockerfile.gemini-test \
     -t gemini-mcp-test:latest . 2>&1 | tee /tmp/gemini-test-build.log; then
    log_ok "Gemini test image built successfully"
    PASSED_TESTS=$((PASSED_TESTS + 1))
else
    log_error "Gemini test image build failed. See /tmp/gemini-test-build.log"
    FAILED_TESTS=$((FAILED_TESTS + 1))
fi

log_info "Testing Gemini CLI installation..."
TOTAL_TESTS=$((TOTAL_TESTS + 1))

if docker run --rm gemini-mcp-test:latest gemini --version > /dev/null 2>&1; then
    log_ok "Gemini CLI is installed"
    PASSED_TESTS=$((PASSED_TESTS + 1))
else
    log_error "Gemini CLI test failed"
    FAILED_TESTS=$((FAILED_TESTS + 1))
fi

# Summary
echo ""
echo "==================================================================="
echo "                         Test Summary"
echo "==================================================================="
echo "Total tests:  $TOTAL_TESTS"
echo -e "Passed:       ${GREEN}$PASSED_TESTS${NC}"
if [ $FAILED_TESTS -gt 0 ]; then
    echo -e "Failed:       ${RED}$FAILED_TESTS${NC}"
else
    echo -e "Failed:       $FAILED_TESTS"
fi
echo "==================================================================="

if [ $FAILED_TESTS -eq 0 ]; then
    echo ""
    log_ok "All tests passed!"
    echo ""
    echo "Next steps:"
    echo "  1. Review test images:     docker images | grep -E 'tooluniverse|claude|codex|gemini'"
    echo "  2. Interactive testing:    docker run --rm -it claude-mcp-test:latest bash"
    echo "  3. Clean up test images:   docker rmi tooluniverse-test claude-mcp-test codex-mcp-test gemini-mcp-test"
    echo "  4. Proceed to Phase 2:     Integration into scbio-docker"
    echo ""
    exit 0
else
    echo ""
    log_error "Some tests failed. Review logs in /tmp/*-build.log"
    echo ""
    exit 1
fi
