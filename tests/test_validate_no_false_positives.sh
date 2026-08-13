#!/usr/bin/env bash
# tests/test_validate_no_false_positives.sh
# A toolkit with disjoint name sets across skills/agents/commands must
# produce zero collision warnings. Guards against the enumeration accidentally
# flagging same-name entries that don't actually overlap (e.g. a basename
# emitted twice from one namespace).
set -u
. "$(dirname "$0")/_lib.sh"

setup_tmpdir
FAKE="$TMPDIR_TEST/fake-toolkit"
build_fake_toolkit "$FAKE"
# build_fake_toolkit ships disjoint names; no further setup is needed.

export SCIAGENT_TOOLKIT="$FAKE"
SCIAGENT="$FAKE/bin/sciagent"

set +e
stderr_out=$("$SCIAGENT" lint --check toolkit 2>&1 >/dev/null)
rc=$?
set -e

if [[ "$rc" -ne 0 ]]; then
    echo "FAIL [$_TEST_NAME] toolkit lint exited $rc on clean fixture (expected 0)" >&2
    echo "--- stderr ---" >&2
    printf '%s\n' "$stderr_out" >&2
    exit 1
fi

# STDERR must be empty — no collision warnings, no other soft-warns.
if [[ -n "$stderr_out" ]]; then
    echo "FAIL [$_TEST_NAME] toolkit lint emitted unexpected STDERR on clean fixture" >&2
    echo "--- stderr ---" >&2
    printf '%s\n' "$stderr_out" >&2
    exit 1
fi

pass
