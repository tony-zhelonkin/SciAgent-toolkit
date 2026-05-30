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

echo "[mllmct] smoke_check_versions + offline pytest in locked sandbox"
# Drop any inherited VIRTUAL_ENV so uv targets THIS skill's sandbox cleanly (no warning).
env -u VIRTUAL_ENV uv run --project "$SKILL" python "$SKILL/checks/smoke_check_versions.py"
env -u VIRTUAL_ENV uv run --project "$SKILL" --extra test pytest "$SKILL/tests" -q
