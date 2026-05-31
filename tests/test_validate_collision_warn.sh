#!/usr/bin/env bash
# tests/test_validate_collision_warn.sh
# A name appears as both a skill and a command. `sciagent validate` must:
#   - exit 0 (cross-namespace collisions are soft-warn only)
#   - emit one warning to STDERR naming both the name and the two kinds
#   - keep STDOUT clean (only the success summary on stdout)
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
stdout_out=$("$SCIAGENT" validate 2>/tmp/_validate_err_$$)
rc=$?
stderr_out=$(cat /tmp/_validate_err_$$)
rm -f /tmp/_validate_err_$$
set -e

# Exit code must be 0 — soft-warn only.
if [[ "$rc" -ne 0 ]]; then
    echo "FAIL [$_TEST_NAME] validate exited $rc on cross-namespace collision (expected 0)" >&2
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

# STDOUT must remain clean: only the success summary line on stdout.
if ! printf '%s\n' "$stdout_out" | grep -q 'all checks passed'; then
    echo "FAIL [$_TEST_NAME] STDOUT missing 'all checks passed' summary" >&2
    printf '%s\n' "$stdout_out" >&2
    exit 1
fi
if printf '%s\n' "$stdout_out" | grep -q 'warning'; then
    echo "FAIL [$_TEST_NAME] STDOUT leaked a warning line (must be STDERR-only)" >&2
    printf '%s\n' "$stdout_out" >&2
    exit 1
fi

pass
