#!/usr/bin/env bash
# tests/test_inject_unknown_name.sh — `sciagent inject <name>` for a name that
# resolves nowhere must exit 1 with the "not found" message and leave the
# project state untouched.
set -u
. "$(dirname "$0")/_lib.sh"

setup_tmpdir
FAKE="$TMPDIR_TEST/fake-toolkit"
build_fake_toolkit "$FAKE"
export SCIAGENT_TOOLKIT="$FAKE"
SCIAGENT="$FAKE/bin/sciagent"

mkdir project && cd project
"$SCIAGENT" activate base >/dev/null

cp .sciagent/manifest.json /tmp/manifest_before.$$.json

set +e
out=$("$SCIAGENT" inject nonexistent-thing 2>&1)
rc=$?
set -e

if [[ "$rc" -ne 1 ]]; then
    echo "FAIL [$_TEST_NAME] expected exit 1 for unknown name, got $rc" >&2
    echo "output: $out" >&2
    exit 1
fi

printf '%s\n' "$out" | grep -qi 'not found' || {
    echo "FAIL [$_TEST_NAME] expected 'not found' message, got: $out" >&2
    exit 1
}
printf '%s\n' "$out" | grep -q "nonexistent-thing" || {
    echo "FAIL [$_TEST_NAME] expected the unknown name in the error, got: $out" >&2
    exit 1
}

# Manifest untouched.
if ! diff -q /tmp/manifest_before.$$.json .sciagent/manifest.json >/dev/null; then
    echo "FAIL [$_TEST_NAME] manifest mutated by a failed inject" >&2
    rm -f /tmp/manifest_before.$$.json
    exit 1
fi
rm -f /tmp/manifest_before.$$.json

pass
