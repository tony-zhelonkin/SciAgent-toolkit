#!/usr/bin/env bash
# Every dispatched verb loads the helpers its execution path references.

set -u
. "$(dirname "$0")/_lib.sh"

setup_tmpdir
FAKE="$TMPDIR_TEST/fake-toolkit"
build_fake_toolkit "$FAKE"
export SCIO_TOOLKIT="$FAKE"
SCIO="$FAKE/bin/scio"

mkdir project
cd project

run_verb() {
    local out
    out=$("$SCIO" "$@" 2>&1 || true)
    if printf '%s\n' "$out" | grep -qi 'command not found'; then
        echo "FAIL [$_TEST_NAME] verb '$*' hit a missing command:" >&2
        printf '%s\n' "$out" >&2
        exit 1
    fi
}

run_verb link
run_verb lint --check toolkit --quiet
run_verb craft

pass
