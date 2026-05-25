#!/usr/bin/env bash
# tests/test_status_shows_collisions.sh
# Activate a role whose YAML mounts a colliding skill+command pair.
# `sciagent status` must include the Notes section listing the overlap and
# annotating it against the allowlist (intentional family overlap when
# present, UNEXPECTED otherwise).
set -u
. "$(dirname "$0")/_lib.sh"

setup_tmpdir
FAKE="$TMPDIR_TEST/fake-toolkit"
build_fake_toolkit "$FAKE"

# Plant the colliding command and a role that mounts both the skill and the
# command under the same name (s_a). The fake toolkit already has skills/s_a.
echo "cmd s_a (deliberate overlap)" > "$FAKE/commands/s_a.md"
cat > "$FAKE/roles/overlap.yaml" <<'EOF'
name: overlap
description: fixture role with skill+command name overlap
skills:
  - s_a
commands:
  - s_a
EOF

# Plant a fixture allowlist documenting the s_a skill+command overlap.
# The status renderer reads $SCIAGENT_TOOLKIT/tests/collision-allowlist.txt;
# the fake toolkit has no tests/ dir until we add one.
mkdir -p "$FAKE/tests"
cat > "$FAKE/tests/collision-allowlist.txt" <<'EOF'
# fixture allowlist
s_a skill,command # fixture: deliberate skill+command overlap
EOF

export SCIAGENT_TOOLKIT="$FAKE"
SCIAGENT="$FAKE/bin/sciagent"

mkdir project && cd project
"$SCIAGENT" activate overlap >/dev/null 2>&1

out=$("$SCIAGENT" status 2>&1)

# Notes header must appear.
if ! printf '%s\n' "$out" | grep -q '^Notes:'; then
    echo "FAIL [$_TEST_NAME] status missing Notes: section" >&2
    echo "--- status output ---" >&2
    printf '%s\n' "$out" >&2
    exit 1
fi

# The colliding name + both kinds must appear in the Notes body.
if ! printf '%s\n' "$out" | grep -q "'s_a'"; then
    echo "FAIL [$_TEST_NAME] Notes body does not name 's_a'" >&2
    printf '%s\n' "$out" >&2
    exit 1
fi
if ! printf '%s\n' "$out" | grep -Eq 'both skill and command'; then
    echo "FAIL [$_TEST_NAME] Notes body does not name kinds 'skill and command'" >&2
    printf '%s\n' "$out" >&2
    exit 1
fi

# Allowlisted: annotation must read "intentional family overlap ...".
if ! printf '%s\n' "$out" | grep -q 'intentional family overlap per tests/collision-allowlist.txt'; then
    echo "FAIL [$_TEST_NAME] Notes annotation missing 'intentional family overlap' label" >&2
    printf '%s\n' "$out" >&2
    exit 1
fi
if printf '%s\n' "$out" | grep -q 'UNEXPECTED'; then
    echo "FAIL [$_TEST_NAME] allowlisted overlap was tagged UNEXPECTED" >&2
    printf '%s\n' "$out" >&2
    exit 1
fi

# Re-test with allowlist scrubbed: the same activation should now tag the
# overlap UNEXPECTED rather than intentional.
: > "$FAKE/tests/collision-allowlist.txt"
out_no_list=$("$SCIAGENT" status 2>&1)
if ! printf '%s\n' "$out_no_list" | grep -q 'UNEXPECTED' ; then
    echo "FAIL [$_TEST_NAME] Notes annotation did not flip to UNEXPECTED when allowlist scrubbed" >&2
    printf '%s\n' "$out_no_list" >&2
    exit 1
fi

pass
