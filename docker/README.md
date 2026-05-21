# Docker Testing Environment for SciAgent-toolkit

This directory contains a Docker-based smoke test for the SciAgent-toolkit
role/skill/template architecture.

## Overview

Only one test image is currently maintained:

- **architecture-test** (`docker/test/Dockerfile.architecture-test`) — validates
  that the role system, templates, and `sciagent activate` work end-to-end in a
  clean Ubuntu container.

The previous suite (`tooluniverse-test`, `claude-mcp-test`, `codex-mcp-test`,
`gemini-mcp-test`, and the `test-all.sh` runner) was removed along with the
MCP infrastructure carve-out. See `CHANGELOG.md` for context.

## Quick Start

```bash
# From the repository root
docker build -f docker/test/Dockerfile.architecture-test -t architecture-test:latest .
```

A successful build prints `=== All Architecture Tests Passed ===` near the end
of the log.

## What It Tests

1. Role-system files exist (`roles/base.yaml`, `sciagent activate`,
   `agents/`, `skills/`).
2. Vendor template files exist (`CLAUDE.md.template`, `AGENTS.md.template`,
   `context.md.template`, `analysis_config.yaml.template`).
3. Templates reference `01_modules` (the canonical submodule path).
4. `sciagent activate base` populates `.claude/agents/` and `.claude/skills/`
   with symlinks.
5. Role activation does not emit spurious warnings from YAML comment parsing.
6. Template files (`AGENTS.md.template`, `CLAUDE.md.template`, `context.md.template`) are present and substitutable.

## Cleanup

```bash
docker rmi architecture-test:latest
docker builder prune -f
```
