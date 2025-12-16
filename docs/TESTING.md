# Testing Guide

**Last Updated:** 2025-12-16

This guide explains how to test SciAgent-toolkit installations and MCP server configurations.

---

## Testing Overview

SciAgent-toolkit provides three types of testing:

| Test Type | Location | Purpose |
|-----------|----------|---------|
| **Installation Tests** | `scripts/test_installation.sh` | Verify MCP setup in current environment |
| **Architecture Tests** | `docker/test/Dockerfile.architecture-test` | Verify role system, templates, profiles |
| **Docker Tests** | `docker/test/test-all.sh` | Full CI/CD testing with isolated containers (11 tests) |

---

## Quick Start

### Test Current Installation

After running `setup-ai.sh` or `configure_mcp_servers.sh`:

```bash
# From project directory with .mcp.json
./toolkits/SciAgent-toolkit/scripts/test_installation.sh
```

This verifies:
- `.mcp.json` exists and is valid JSON
- Claude Code CLI is installed
- MCP servers are configured
- ToolUniverse environment is functional

### Test with Docker (CI/CD)

For isolated testing or CI pipelines:

```bash
cd toolkits/SciAgent-toolkit/docker/test
./test-all.sh
```

This builds and tests:
- Base ToolUniverse image
- Claude Code integration
- Codex CLI integration
- Gemini CLI integration

---

## Installation Testing

### What `test_installation.sh` Checks

```
Test 1: MCP Configuration
  ✓ .mcp.json exists
  ✓ Valid JSON syntax
  ✓ Contains mcpServers key

Test 2: Claude Code CLI
  ✓ claude command available
  ✓ Version check passes

Test 3: ToolUniverse
  ✓ tooluniverse-env/ directory exists
  ✓ uv package manager available
  ✓ tooluniverse-smcp-stdio responds to --help

Test 4: MCP Server Connectivity
  ✓ Sequential Thinking (npx available)
  ✓ Context7 (npx available)
  ✓ PAL (uvx available, API keys optional)
  ✓ Serena (uvx available)
```

### Running Specific Tests

```bash
# Full test suite
./scripts/test_installation.sh

# Test only ToolUniverse
./scripts/debug_tooluniverse.sh --verbose

# Verify MCP configuration
claude mcp list
```

### Expected Output

```
==========================================
SciAgent Toolkit Installation Test
==========================================

Test 1: MCP Configuration
✓ MCP config file exists
✓ MCP config is valid JSON
✓ MCP config has mcpServers key

Test 2: Claude Code CLI
✓ Claude Code CLI installed
✓ Claude version: 2.0.64

Test 3: ToolUniverse
✓ ToolUniverse environment found
✓ tooluniverse-smcp-stdio responds

Test 4: MCP Servers
✓ Sequential Thinking available
✓ Context7 available
⚠ PAL available (no API keys configured)
✓ Serena available

==========================================
Results: 10 passed, 0 failed, 1 warning
==========================================
```

---

## Docker Testing

### Test Images

| Image | Dockerfile | Purpose |
|-------|------------|---------|
| `architecture-test` | `Dockerfile.architecture-test` | Role system, templates, profiles |
| `tooluniverse-test` | `Dockerfile.tooluniverse-test` | Base with ToolUniverse |
| `claude-mcp-test` | `Dockerfile.claude-test` | + Claude Code |
| `codex-mcp-test` | `Dockerfile.codex-test` | + Codex CLI |
| `gemini-mcp-test` | `Dockerfile.gemini-test` | + Gemini CLI |

### Running Docker Tests

```bash
cd toolkits/SciAgent-toolkit/docker/test

# Run all tests (11 tests)
./test-all.sh

# Run architecture tests only
docker build -f Dockerfile.architecture-test -t architecture-test:latest ../..

# Build individual images
docker build -f Dockerfile.tooluniverse-test -t tooluniverse-test:latest ../..
docker build -f Dockerfile.claude-test -t claude-mcp-test:latest ../..

# Interactive testing
docker run --rm -it claude-mcp-test:latest bash
```

### Test Output

```
===================================================================
           ToolUniverse Docker Integration Test Suite
===================================================================

[STEP] Test 0: Building architecture test image...
[OK] Architecture tests passed (role system, templates, profiles)

[STEP] Test 1: Building base ToolUniverse test image...
[OK] Base image built successfully
[OK] ToolUniverse MCP server responds to --help

[STEP] Test 2: Building Claude Code MCP test image...
[OK] Claude test image built successfully
[OK] Claude Code is installed
[OK] Claude MCP config is valid JSON

[STEP] Test 3: Building Codex CLI MCP test image...
[OK] Codex test image built successfully
[OK] Codex CLI is installed
[OK] Codex config is valid TOML

[STEP] Test 4: Building Gemini CLI MCP test image...
[OK] Gemini test image built successfully
[OK] Gemini CLI is installed

===================================================================
                         Test Summary
===================================================================
Total tests:  11
Passed:       11
Failed:       0
===================================================================
```

### Architecture Tests (Test 0)

The architecture test validates the modularization changes:

| Test | What It Checks |
|------|----------------|
| Test 1 | Role system files exist (`roles/base.yaml`, `activate-role.sh`, `agents/`, `skills/`) |
| Test 2 | Template files exist (`CLAUDE.md.template`, `GEMINI.md.template`, etc.) |
| Test 3 | Templates reference `01_modules` (not old `01_scripts`) |
| Test 4 | Role activation creates correct symlinks |
| Test 5 | Profile switcher has API key substitution code |
| Test 6 | Profile switching creates valid JSON |
| Test 7 | `research-full.mcp.json` has `--include-tools` |
| Test 8 | `setup-ai.sh` has template support functions |

### Cleanup Test Images

```bash
docker rmi tooluniverse-test claude-mcp-test codex-mcp-test gemini-mcp-test
```

---

## MCP Profile Testing

### Testing Profile Switcher

```bash
# Show available profiles
./scripts/switch-mcp-profile.sh

# Switch and verify
./scripts/switch-mcp-profile.sh minimal
cat .mcp.json | jq '.mcpServers | keys'

# Expected output for minimal:
# ["context7", "sequential-thinking"]
```

### Testing Context Usage

After switching profiles, verify context usage in Claude Code:

```bash
claude
# Inside Claude:
/context
```

Expected context by profile:

| Profile | MCP Tokens | Percentage |
|---------|------------|------------|
| minimal | ~3k | 1.5% |
| coding | ~25k | 12.5% |
| codebase | ~75k | 37.5% |
| research-lite | ~30k | 15% |
| research-full | ~50k | 25% |
| full | ~100k | 50% |

---

## Debugging Failed Tests

### MCP Configuration Issues

```bash
# Check .mcp.json syntax
python3 -m json.tool .mcp.json

# Verify server paths are absolute
cat .mcp.json | grep -E '"command"|"args"'

# Test MCP server directly
npx -y @modelcontextprotocol/server-sequential-thinking --help
```

### ToolUniverse Issues

```bash
# Run debug script
./scripts/debug_tooluniverse.sh --verbose

# Manual verification
uv --directory ./tooluniverse-env run tooluniverse-smcp-stdio --help

# Check if correct command is used (1.0.14+ uses stdio)
uv --directory ./tooluniverse-env run pip show tooluniverse
```

### PAL Connection Issues

PAL requires at least one API key:

```bash
# Check if keys are set
echo $GEMINI_API_KEY
echo $OPENAI_API_KEY

# Set keys in .env
cat >> .devcontainer/.env << EOF
GEMINI_API_KEY=your-key-here
EOF

# Reconfigure
./scripts/configure_mcp_servers.sh --force
```

### Claude Code Issues

```bash
# Run diagnostics
claude doctor

# Check MCP server status
claude mcp list

# Verify installation path
which claude
claude --version
```

---

## CI/CD Integration

### GitHub Actions Example

```yaml
name: Test SciAgent Toolkit

on: [push, pull_request]

jobs:
  test:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4
        with:
          submodules: recursive

      - name: Build and test Docker images
        run: |
          cd toolkits/SciAgent-toolkit/docker/test
          ./test-all.sh
```

### Local CI Simulation

```bash
# Clean environment test
docker run --rm -v $(pwd):/workspace -w /workspace ubuntu:22.04 bash -c "
  apt-get update && apt-get install -y curl git python3
  cd toolkits/SciAgent-toolkit
  ./scripts/setup_mcp_infrastructure.sh --skip-codex
  ./scripts/test_installation.sh
"
```

---

## Test Files Reference

```
toolkits/SciAgent-toolkit/
├── scripts/
│   ├── test_installation.sh      # Main installation test
│   ├── debug_tooluniverse.sh     # ToolUniverse diagnostics
│   └── common.sh                 # Shared test utilities
├── docker/test/
│   ├── test-all.sh               # Docker test orchestrator
│   ├── Dockerfile.tooluniverse-test
│   ├── Dockerfile.claude-test
│   ├── Dockerfile.codex-test
│   └── Dockerfile.gemini-test
└── docs/
    └── TESTING.md                # This document
```

---

## Adding New Tests

### Script-Based Tests

Add to `scripts/test_installation.sh`:

```bash
# Test N: Your New Test
echo "Test N: Description"
if your_test_command; then
    test_pass "Test description"
else
    test_fail "Test description"
fi
```

### Docker-Based Tests

1. Create `docker/test/Dockerfile.yourtest-test`
2. Add build/test steps to `docker/test/test-all.sh`
3. Update this documentation

---

## Troubleshooting

| Issue | Solution |
|-------|----------|
| "npx not found" | Install Node.js 18+ or run `ensure_nvm` from common.sh |
| "uv not found" | Run `curl -LsSf https://astral.sh/uv/install.sh \| sh` |
| "claude not found" | Run `./scripts/install_claude.sh` |
| "Invalid JSON" | Check .mcp.json with `python3 -m json.tool .mcp.json` |
| "PAL failed" | Add API keys to `.devcontainer/.env` |
| "Context overflow" | Switch to lighter MCP profile |
