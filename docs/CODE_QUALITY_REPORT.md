# Code Quality Report

## Overview
This report documents the findings from inspecting the `SciAgent-toolkit` repository for code quality, syntax, linting, separation of concerns, and DRY (Don't Repeat Yourself) principles.

## Findings

### 1. DRY Principles (Don't Repeat Yourself)
**Status: Critical**

There is significant code duplication across the shell scripts, particularly in the `scripts/` and `scripts/mcp_servers/` directories.

*   **Color Definitions & Logging:**
    *   Found in almost every script: `scripts/setup_mcp_infrastructure.sh`, `scripts/install_claude.sh`, `scripts/install_codex.sh`, `scripts/install_gemini.sh`, `scripts/test_installation.sh`, and all `scripts/mcp_servers/*.sh` files.
    *   All these files define the same color variables (`RED`, `GREEN`, etc.) and logging functions (`log_info`, `log_ok`, etc.).

*   **OS Detection:**
    *   `setup_mcp_infrastructure.sh` and others implement their own OS detection logic (`uname -s`).

*   **Dependency Checks:**
    *   Node.js/npm installation logic is repeated.
    *   `uv`/`uvx` installation logic is present in multiple places.
    *   `curl`/`git` checks are sporadic.

**Recommendation:**
Create a `scripts/common.sh` library file containing:
*   Color definitions
*   Logging functions
*   OS detection logic
*   Common dependency installers/checkers (`ensure_node`, `ensure_uv`, `ensure_git`, `ensure_curl`)

Source these utilities in each script:
```bash
source "$(dirname "${BASH_SOURCE[0]}")/../common.sh"
# Adjust path based on script location
```

### 2. Separation of Concerns
**Status: Good**

The repository generally follows good separation of concerns:
*   **Main Orchestrator:** `setup_mcp_infrastructure.sh` handles high-level coordination.
*   **Component Installers:** Individual scripts (`install_claude.sh`, `setup_serena.sh`) handle specific components.
*   **Configuration:** `configure_mcp_servers.sh` is dedicated to generating the config file.
*   **Testing:** `test_installation.sh` and `docker/test/` separate testing concerns.

### 3. Syntax & Linting
**Status: Acceptable**

*   Scripts use `set -euo pipefail` for robust error handling.
*   Variables are generally quoted.
*   Usage of `command -v` for checking existence is consistent.
*   Shebangs are consistently `#!/usr/bin/env bash`.

### 4. Hardcoded Paths & Versions
**Status: Needs Attention**

*   Some scripts have hardcoded versions (e.g., `node@20`).
*   `scripts/test_installation.sh` has hardcoded checks that might drift from the actual installation logic.

### 5. Testing
**Status: Good**

*   `scripts/test_installation.sh` provides a quick local check.
*   `docker/test/` contains a robust suite of Dockerfiles to test different components in isolation and together.
*   `docker/test/test-all.sh` orchestrates the docker tests.

## Action Plan
1.  **Refactor:** Create `scripts/common.sh` to centralize shared logic (colors, logs, OS detection, dependency checks).
2.  **Update:** Modify all scripts in `scripts/` and `scripts/mcp_servers/` to source `common.sh` and remove duplicated code.
3.  **Verify:** Run `scripts/test_installation.sh` and `docker/test/test-all.sh` to ensure no regressions.
