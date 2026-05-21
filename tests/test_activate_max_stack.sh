#!/usr/bin/env bash
# Three role args → error.
set -u
. "$(dirname "$0")/_lib.sh"

setup_tmpdir
FAKE="$TMPDIR_TEST/fake-toolkit"
build_fake_toolkit "$FAKE"
export SCIAGENT_TOOLKIT="$FAKE"
SCIAGENT="$FAKE/bin/sciagent"

mkdir project && cd project
out=$("$SCIAGENT" activate base reviewer alpha 2>&1) && rc=0 || rc=$?

if [[ $rc -eq 0 ]]; then
    echo "FAIL: expected non-zero exit on 3 roles"
    echo "$out"
    exit 1
fi

if ! echo "$out" | grep -q 'Maximum stack depth'; then
    echo "FAIL: expected error mentioning 'Maximum stack depth'"
    echo "$out"
    exit 1
fi

# No state should have been written.
[[ -e .sciagent ]] && { echo "FAIL: .sciagent exists after error"; exit 1; }
pass
