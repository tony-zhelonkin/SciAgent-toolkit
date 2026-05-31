#!/usr/bin/env bash
# activate base; attempt to eject a skill that is stack-mounted (s_a is in base);
# assert exit 1 and error message mentions 'sciagent deactivate'.
set -u
. "$(dirname "$0")/_lib.sh"

setup_tmpdir
FAKE="$TMPDIR_TEST/fake-toolkit"
build_fake_toolkit "$FAKE"
export SCIAGENT_TOOLKIT="$FAKE"
SCIAGENT="$FAKE/bin/sciagent"

mkdir project && cd project
"$SCIAGENT" activate base >/dev/null

# s_a is in base role — ejecting it should fail with the deactivate pointer.
out=$("$SCIAGENT" eject s_a 2>&1)
rc=$?

if [[ "$rc" -eq 0 ]]; then
    echo "FAIL [$_TEST_NAME] expected exit 1 when ejecting stack-mounted skill, got 0" >&2
    echo "output: $out" >&2
    exit 1
fi

printf '%s\n' "$out" | grep -qi 'deactivate' || {
    echo "FAIL [$_TEST_NAME] error message should mention 'sciagent deactivate', got: $out" >&2
    exit 1
}

# s_b is in both base and reviewer roles — also stack-mounted.
out2=$("$SCIAGENT" eject s_b 2>&1)
rc2=$?
if [[ "$rc2" -eq 0 ]]; then
    echo "FAIL [$_TEST_NAME] expected exit 1 when ejecting stack-mounted skill s_b, got 0" >&2
    exit 1
fi

pass
