#!/usr/bin/env bash
# tests/test_validate_quiet_flag.sh — `sciagent validate --quiet` suppresses the
# "all checks passed" summary line on success.
#
# Verifies the --quiet contract from validate.sh: exit 0 on clean tree AND
# empty stdout (the "all checks passed" line is the only stdout emission when
# validate succeeds; --quiet must suppress it).
set -u
. "$(dirname "$0")/_lib.sh"

setup_tmpdir
FAKE="$TMPDIR_TEST/fake-toolkit"
build_fake_toolkit "$FAKE"
export SCIAGENT_TOOLKIT="$FAKE"
SCIAGENT="$FAKE/bin/sciagent"

# --- Without --quiet: the summary line must appear. ---
out_verbose=$("$SCIAGENT" validate 2>&1)
if ! printf '%s\n' "$out_verbose" | grep -q 'all checks passed'; then
    echo "FAIL [$_TEST_NAME] 'all checks passed' absent from validate output without --quiet" >&2
    echo "--- output ---" >&2
    printf '%s\n' "$out_verbose" >&2
    exit 1
fi

# --- With --quiet: stdout must be empty and exit code 0. ---
set +e
out_quiet=$("$SCIAGENT" validate --quiet 2>&1)
rc=$?
set -e

if [[ "$rc" -ne 0 ]]; then
    echo "FAIL [$_TEST_NAME] validate --quiet exited $rc (expected 0) on clean tree" >&2
    echo "--- output ---" >&2
    printf '%s\n' "$out_quiet" >&2
    exit 1
fi

if [[ -n "$out_quiet" ]]; then
    echo "FAIL [$_TEST_NAME] validate --quiet produced output on clean tree; expected none" >&2
    echo "--- output ---" >&2
    printf '%s\n' "$out_quiet" >&2
    exit 1
fi

pass
