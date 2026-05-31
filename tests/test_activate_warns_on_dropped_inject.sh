#!/usr/bin/env bash
# tests/test_activate_warns_on_dropped_inject.sh — re-activate over a stack
# with injected entries must emit a STDERR warning listing each dropped row.
# Verifies the clean-slate behavior documented at docs/architecture.md §5
# "Activation semantics" is explicit at the CLI surface and not silent.
set -u
. "$(dirname "$0")/_lib.sh"

setup_tmpdir
FAKE="$TMPDIR_TEST/fake-toolkit"
build_fake_toolkit "$FAKE"
export SCIAGENT_TOOLKIT="$FAKE"
SCIAGENT="$FAKE/bin/sciagent"

mkdir project && cd project
"$SCIAGENT" activate base >/dev/null

# Inject one of each kind so the warning lists all three.
"$SCIAGENT" inject s_c   >/dev/null   # skill
"$SCIAGENT" inject ag_b  >/dev/null   # agent
"$SCIAGENT" inject c_b   >/dev/null   # command

# Re-activate with a different stack. STDERR must enumerate the dropped rows.
err=$("$SCIAGENT" activate base alpha 2>&1 >/dev/null)

echo "$err" | grep -q "activate is a clean-slate operation; dropping injected entries:" || {
    echo "FAIL [$_TEST_NAME] missing dropped-injected header on STDERR" >&2
    printf '%s\n' "$err" >&2
    exit 1
}
echo "$err" | grep -Eq '^  - skill s_c$'   || { echo "FAIL [$_TEST_NAME] missing 'skill s_c' line" >&2;   exit 1; }
echo "$err" | grep -Eq '^  - agent ag_b$'  || { echo "FAIL [$_TEST_NAME] missing 'agent ag_b' line" >&2;  exit 1; }
echo "$err" | grep -Eq '^  - command c_b$' || { echo "FAIL [$_TEST_NAME] missing 'command c_b' line" >&2; exit 1; }
echo "$err" | grep -q "to preserve, run 'sciagent deactivate' first" || {
    echo "FAIL [$_TEST_NAME] missing preservation hint on STDERR" >&2
    exit 1
}

# Side effect: injected rows are actually gone after activation (the warning
# is descriptive of real teardown, not a dry-run).
if grep -q '"skill": "s_c"' .sciagent/manifest.json; then
    echo "FAIL [$_TEST_NAME] injected skill s_c still present in manifest after re-activate" >&2
    exit 1
fi

# Negative: a fresh activate with NO prior injects emits no such warning.
"$SCIAGENT" deactivate >/dev/null 2>&1 || true
"$SCIAGENT" activate base >/dev/null
err2=$("$SCIAGENT" activate base alpha 2>&1 >/dev/null)
if echo "$err2" | grep -q 'dropping injected entries'; then
    echo "FAIL [$_TEST_NAME] warning fired with no injected entries present" >&2
    printf '%s\n' "$err2" >&2
    exit 1
fi

pass
