#!/usr/bin/env bash
# tests/test_validate_collision_quiet.sh
# `sciagent validate --quiet` must suppress cross-namespace collision warnings
# the same way it suppresses the "all checks passed" success line. Scripted
# callers want silent-on-success across every soft-warn surface.
set -u
. "$(dirname "$0")/_lib.sh"

setup_tmpdir
FAKE="$TMPDIR_TEST/fake-toolkit"
build_fake_toolkit "$FAKE"

# Same fixture shape as test_validate_collision_warn.sh.
echo "cmd s_a (deliberate collision)" > "$FAKE/commands/s_a.md"

export SCIAGENT_TOOLKIT="$FAKE"
SCIAGENT="$FAKE/bin/sciagent"

set +e
out=$("$SCIAGENT" validate --quiet 2>&1)
rc=$?
set -e

if [[ "$rc" -ne 0 ]]; then
    echo "FAIL [$_TEST_NAME] validate --quiet exited $rc on collision (expected 0)" >&2
    echo "--- output ---" >&2
    printf '%s\n' "$out" >&2
    exit 1
fi

if [[ -n "$out" ]]; then
    echo "FAIL [$_TEST_NAME] validate --quiet produced output on collision; expected none" >&2
    echo "--- output ---" >&2
    printf '%s\n' "$out" >&2
    exit 1
fi

pass
