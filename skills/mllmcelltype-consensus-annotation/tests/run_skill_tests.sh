#!/usr/bin/env bash
# Skill-local test entry point, discovered by the toolkit's tests/run-all.sh sweep.
# Skips GRACEFULLY (exit 0 with a notice) when uv or the lock are absent, so the toolkit's
# bash CI never blocks on an un-bootstrapped machine. When bootstrapped, runs the offline
# pytest suite inside the version-locked sandbox.
set -euo pipefail

SKILL="$(cd "$(dirname "$0")/.." && pwd)"

if ! command -v uv >/dev/null 2>&1; then
  echo "SKIP [mllmct]: uv not on PATH (skill not bootstrapped)"; exit 0
fi
if [ ! -f "$SKILL/uv.lock" ]; then
  echo "SKIP [mllmct]: $SKILL/uv.lock absent (run: uv sync --project \"$SKILL\")"; exit 0
fi

# --- Cache dir: make the outcome independent of whether $HOME is writable ----
# uv writes to ${XDG_CACHE_HOME:-$HOME/.cache}/uv and hard-errors if it cannot
# create it. Under a sandboxed/read-only home that turned this into a FAIL —
# so the toolkit suite scored 75/0 or 74/1 depending purely on the ambient
# environment. That non-determinism was the defect, not the failure itself.
#
# Three tiers, each for a distinct reason:
#   1. UV_CACHE_DIR already set  → honour it; the caller/CI knows best.
#   2. ambient default writable  → leave uv alone; a warm cache is fast and is
#                                  the only tier that can bootstrap a missing
#                                  .venv without network.
#   3. neither                   → private temp cache, removed on exit. Costs
#                                  ~28K when the locked sandbox is already
#                                  synced (nothing to download); if it is not,
#                                  uv reports a real resolution error rather
#                                  than a permission error about $HOME.
if [ -z "${UV_CACHE_DIR:-}" ]; then
  _ambient="${XDG_CACHE_HOME:-$HOME/.cache}/uv"
  if mkdir -p "$_ambient" 2>/dev/null && [ -w "$_ambient" ]; then
    : # tier 2 — uv's own default is usable; do not override it
  else
    UV_CACHE_DIR=$(mktemp -d)
    export UV_CACHE_DIR
    # Armed only on tier 3, so the trap body is never a no-op whose nonzero
    # status could leak into the exit code under `set -e`.
    # shellcheck disable=SC2064  # expand the path now, not at trap time
    trap "rm -rf '$UV_CACHE_DIR'" EXIT
    echo "[mllmct] \$HOME cache not writable — using a temp UV_CACHE_DIR"
  fi
fi

echo "[mllmct] smoke_check_versions + offline pytest in locked sandbox"
# Drop any inherited VIRTUAL_ENV so uv targets THIS skill's sandbox cleanly (no warning).
env -u VIRTUAL_ENV uv run --project "$SKILL" python "$SKILL/checks/smoke_check_versions.py"
env -u VIRTUAL_ENV uv run --project "$SKILL" --extra test pytest "$SKILL/tests" -q
