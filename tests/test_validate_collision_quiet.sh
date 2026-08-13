#!/usr/bin/env bash
# tests/test_validate_collision_quiet.sh
# `sciagent lint --check toolkit --quiet` suppresses collision warnings.
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
out=$("$SCIAGENT" lint --check toolkit --quiet 2>&1)
rc=$?
set -e

if [[ "$rc" -ne 0 ]]; then
    echo "FAIL [$_TEST_NAME] toolkit lint --quiet exited $rc on collision (expected 0)" >&2
    echo "--- output ---" >&2
    printf '%s\n' "$out" >&2
    exit 1
fi

if [[ -n "$out" ]]; then
    echo "FAIL [$_TEST_NAME] toolkit lint --quiet produced output on collision; expected none" >&2
    echo "--- output ---" >&2
    printf '%s\n' "$out" >&2
    exit 1
fi

pass
