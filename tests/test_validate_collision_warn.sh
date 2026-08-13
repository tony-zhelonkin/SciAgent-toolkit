#!/usr/bin/env bash
# tests/test_validate_collision_warn.sh
# A name appears as both a skill and a command. The toolkit lint check must:
#   - exit 0 (cross-namespace collisions are soft-warn only)
#   - emit one warning to STDERR naming both the name and the two kinds
#   - keep STDOUT clean
set -u
. "$(dirname "$0")/_lib.sh"

setup_tmpdir
FAKE="$TMPDIR_TEST/fake-toolkit"
build_fake_toolkit "$FAKE"

# Plant a command whose basename also matches an existing skill (s_a).
# The fake toolkit already ships skills/s_a/SKILL.md, so creating
# commands/s_a.md is sufficient to trigger the collision check.
echo "cmd s_a (deliberate collision)" > "$FAKE/commands/s_a.md"

export SCIAGENT_TOOLKIT="$FAKE"
SCIAGENT="$FAKE/bin/sciagent"

set +e
stdout_out=$("$SCIAGENT" lint --check toolkit 2>/tmp/_validate_err_$$)
rc=$?
stderr_out=$(cat /tmp/_validate_err_$$)
rm -f /tmp/_validate_err_$$
set -e

# Exit code must be 0 — soft-warn only.
if [[ "$rc" -ne 0 ]]; then
    echo "FAIL [$_TEST_NAME] toolkit lint exited $rc on cross-namespace collision (expected 0)" >&2
    echo "--- stderr ---" >&2
    printf '%s\n' "$stderr_out" >&2
    exit 1
fi

# STDERR must carry the warning line, naming the colliding name + both kinds.
if ! printf '%s\n' "$stderr_out" | grep -q "'s_a'"; then
    echo "FAIL [$_TEST_NAME] STDERR warning does not name the colliding name 's_a'" >&2
    printf '%s\n' "$stderr_out" >&2
    exit 1
fi
if ! printf '%s\n' "$stderr_out" | grep -Eq 'both skill and command'; then
    echo "FAIL [$_TEST_NAME] STDERR warning does not name the kinds 'skill and command'" >&2
    printf '%s\n' "$stderr_out" >&2
    exit 1
fi

# STDOUT remains empty; warnings use stderr.
if [[ -n "$stdout_out" ]]; then
    echo "FAIL [$_TEST_NAME] toolkit lint produced stdout" >&2
    printf '%s\n' "$stdout_out" >&2
    exit 1
fi

pass
