#!/usr/bin/env bash
# tests/test_status_no_notes_when_clean.sh
# Activate a role with no cross-namespace collisions. `sciagent status` must
# NOT emit the Notes section — the section is reserved for actionable signal
# and would only add noise on a clean stack.
set -u
. "$(dirname "$0")/_lib.sh"

setup_tmpdir
FAKE="$TMPDIR_TEST/fake-toolkit"
build_fake_toolkit "$FAKE"
# The fake toolkit ships disjoint names — no collisions exist anywhere.

export SCIAGENT_TOOLKIT="$FAKE"
SCIAGENT="$FAKE/bin/sciagent"

mkdir project && cd project
"$SCIAGENT" activate base >/dev/null 2>&1

out=$("$SCIAGENT" status 2>&1)

if printf '%s\n' "$out" | grep -q '^Notes:'; then
    echo "FAIL [$_TEST_NAME] status emitted Notes: section on a clean stack" >&2
    echo "--- status output ---" >&2
    printf '%s\n' "$out" >&2
    exit 1
fi

# Sanity: the Harness line must still be present (we didn't accidentally
# truncate the rest of the renderer).
if ! printf '%s\n' "$out" | grep -q '^Harness:'; then
    echo "FAIL [$_TEST_NAME] status missing Harness: line (renderer broke?)" >&2
    printf '%s\n' "$out" >&2
    exit 1
fi

pass
