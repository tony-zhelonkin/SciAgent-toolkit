# Docker Testing Environment for SciAgent-toolkit

This directory contains Docker-based testing infrastructure for validating the SciAgent-toolkit installation, architecture, and MCP server integration.

## Overview

The test suite consists of five Docker images built in sequence (11 total tests):

1. **architecture-test** - Validates role system, templates, and profile switching
2. **tooluniverse-test** - Base image with ToolUniverse and dependencies
3. **claude-mcp-test** - Adds Claude Code with MCP configuration
4. **codex-mcp-test** - Adds Codex CLI with MCP configuration
5. **gemini-mcp-test** - Adds Gemini CLI

## Quick Start

### Run All Tests

```bash
cd /path/to/SciAgent-toolkit
./docker/test/test-all.sh
```

This will:
- Build all three test images
- Verify installations
- Validate MCP configurations
- Report success/failure

### Build Individual Images

```bash
# Build base ToolUniverse test image
docker build -f docker/test/Dockerfile.tooluniverse-test \
  -t tooluniverse-test:latest .

# Build Claude Code test image (requires base image)
docker build -f docker/test/Dockerfile.claude-test \
  -t claude-mcp-test:latest .

# Build Codex CLI test image (requires base image)
docker build -f docker/test/Dockerfile.codex-test \
  -t codex-mcp-test:latest .
```

## Image Details

### architecture-test

**Base:** ubuntu:22.04

**Purpose:** Validates modularization changes (role system, templates, profiles)

**Tests Performed:**
1. Role system files exist (`roles/base.yaml`, `activate-role.sh`, `agents/`, `skills/`)
2. Template files exist (`CLAUDE.md.template`, `GEMINI.md.template`, etc.)
3. Templates reference `01_modules` (not old `01_scripts`)
4. Role activation creates correct symlinks
5. Profile switcher has API key substitution code
6. Profile switching creates valid JSON
7. `research-full.mcp.json` has `--include-tools`
8. `setup-ai.sh` has template support functions

**Size:** ~1.5GB

### tooluniverse-test

**Base:** ubuntu:22.04

**Installed:**
- Python 3.10
- uv package manager
- Node.js 20+
- ToolUniverse MCP server (`/opt/sciagent/scripts/tooluniverse-env`)

**Validation:**
- ToolUniverse responds to `--help`

**Size:** ~2GB

### claude-mcp-test

**Base:** tooluniverse-test:latest

**Additional:**
- Claude Code (latest from claude.ai/install.sh)
- MCP configuration (`/root/.config/claude/mcp-config.json`)
  - ToolUniverse server (PackageTool excluded)
  - Sequential Thinking server

**Validation:**
- Claude Code `--version` works
- MCP config is valid JSON

**Size:** ~2.5GB

### codex-mcp-test

**Base:** tooluniverse-test:latest

**Additional:**
- Node.js 22+ (required for Codex CLI)
- Codex CLI (`@openai/codex` npm package)
- MCP configuration (`/root/.codex/config.toml`)
  - ToolUniverse server (PackageTool excluded)
  - Sequential Thinking server

**Validation:**
- Codex CLI `--version` works
- MCP config is valid TOML

**Size:** ~2.7GB

## Interactive Testing

### Test Claude Code Integration

```bash
docker run --rm -it claude-mcp-test:latest bash

# Inside container:
claude --version
cat /root/.config/claude/mcp-config.json

# Test MCP server manually
uv --directory /opt/sciagent/tooluniverse-env run tooluniverse-smcp-stdio --help

# Note: Authentication requires interactive login
# claude mcp list  # Would show MCP servers if authenticated
```

### Test Codex CLI Integration

```bash
docker run --rm -it codex-mcp-test:latest bash

# Inside container:
codex --version
cat /root/.codex/config.toml

# Test MCP server manually
uv --directory /opt/sciagent/tooluniverse-env run tooluniverse-smcp-stdio --help

# Note: Authentication requires API key or ChatGPT login
# codex mcp list  # Would show MCP servers if authenticated
```

### Test ToolUniverse Functionality

```bash
docker run --rm -it tooluniverse-test:latest bash

# Inside container:
cd /opt/sciagent/tooluniverse-env

# Test with help
uv run tooluniverse-smcp-stdio --help

# Test with list tools (requires authentication for some features)
uv run python -c "from tooluniverse import ToolUniverse; tu = ToolUniverse(); tu.load_tools(); print(f'Loaded {len(tu.all_tools)} tools')"
```

## Troubleshooting

### Build Failures

**Issue:** Base image build fails

**Solution:**
- Check internet connection (downloads uv, Node.js)
- Review build logs: `/tmp/tooluniverse-test-build.log`
- Ensure Docker has sufficient disk space

**Issue:** Claude Code installation fails

**Solution:**
- Verify https://claude.ai/install.sh is accessible
- Check build logs: `/tmp/claude-test-build.log`

**Issue:** Codex CLI installation fails

**Solution:**
- Verify npm package `@openai/codex` exists
- Check Node.js version is 22+: `docker run --rm codex-mcp-test:latest node --version`

### Runtime Issues

**Issue:** MCP server won't start

**Solution:**
```bash
# Test manually
docker run --rm tooluniverse-test:latest \
  bash -c "uv --directory /opt/sciagent/tooluniverse-env run tooluniverse-smcp-stdio --help"

# Check logs
docker run --rm -it tooluniverse-test:latest bash
# Inside: Check /opt/sciagent/tooluniverse-env/ directory exists
```

**Issue:** Claude/Codex can't find MCP servers

**Solution:**
- Verify config file exists and is valid JSON/TOML
- Check command paths are absolute
- Ensure MCP server binary exists in specified path

## Cleanup

### Remove Test Images

```bash
docker rmi architecture-test:latest
docker rmi tooluniverse-test:latest
docker rmi claude-mcp-test:latest
docker rmi codex-mcp-test:latest
docker rmi gemini-mcp-test:latest
```

### Remove Build Cache

```bash
docker builder prune -f
```

## Next Steps

After successful tests:

1. **Proceed to scbio-docker integration** - See main README.md
2. **Add API keys** - For full ToolUniverse functionality (optional)
3. **Test authentication** - Sign in to Claude/Codex and verify MCP servers load

## Test Results

### Expected Output (test-all.sh)

```
===================================================================
           ToolUniverse Docker Integration Test Suite
===================================================================

[STEP] Test 0: Building architecture test image...
[OK] Architecture tests passed

[STEP] Test 1: Building base ToolUniverse test image...
[OK] Base image built successfully
[INFO] Testing ToolUniverse MCP server --help...
[OK] ToolUniverse MCP server responds to --help

[STEP] Test 2: Building Claude Code MCP test image...
[OK] Claude test image built successfully
[INFO] Testing Claude Code installation...
[OK] Claude Code is installed
[INFO] Testing Claude MCP config validity...
[OK] Claude MCP config is valid JSON

[STEP] Test 3: Building Codex CLI MCP test image...
[OK] Codex test image built successfully
[INFO] Testing Codex CLI installation...
[OK] Codex CLI is installed
[INFO] Testing Codex config validity...
[OK] Codex config is valid TOML

[STEP] Test 4: Building Gemini CLI MCP test image...
[OK] Gemini test image built successfully
[INFO] Testing Gemini CLI installation...
[OK] Gemini CLI is installed

===================================================================
                         Test Summary
===================================================================
Total tests:  11
Passed:       11
Failed:       0
===================================================================

[OK] All tests passed!
```

## CI/CD Integration (Future)

These tests can be integrated into GitHub Actions:

```yaml
name: Docker Test Suite

on: [push, pull_request]

jobs:
  test:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v3
      - name: Run Docker tests
        run: ./docker/test/test-all.sh
```

## Support

For issues:
- **Test failures:** Check build logs in `/tmp/*-build.log`
- **Installation issues:** See main [TROUBLESHOOTING.md](../docs/TROUBLESHOOTING.md)
- **GitHub Issues:** https://github.com/your-repo/SciAgent-toolkit/issues
