#!/usr/bin/env bash
# Quick test of ToolUniverse MCP server
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TOOLUNIVERSE_ENV="${SCRIPT_DIR}/tooluniverse-env"

echo "Testing ToolUniverse MCP server..."
uv --directory "${TOOLUNIVERSE_ENV}" run python -m tooluniverse.smcp_server --help
