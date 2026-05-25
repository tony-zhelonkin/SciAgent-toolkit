#!/usr/bin/env bash
# tests/test_collision_allowlist_blocks_unknown.sh
# The CI uniqueness test (test_no_duplicate_basenames.sh) hard-fails when the
# toolkit contains a cross-namespace collision that is NOT documented in
# tests/collision-allowlist.txt. Drives the test via env-var overrides so we
# don't have to mutate the real allowlist file.
set -u
. "$(dirname "$0")/_lib.sh"

setup_tmpdir
FAKE="$TMPDIR_TEST/fake-toolkit"
build_fake_toolkit "$FAKE"

# Plant a command whose basename collides with an existing skill (s_a),
# with NO matching entry in the empty allowlist below.
echo "cmd s_a (deliberate unlisted collision)" > "$FAKE/commands/s_a.md"

# Empty allowlist — every collision should hard-fail.
ALLOWLIST="$TMPDIR_TEST/allowlist.txt"
printf '# empty allowlist for test\n' > "$ALLOWLIST"

set +e
out=$(COLLISION_TEST_TOOLKIT="$FAKE" \
      COLLISION_TEST_ALLOWLIST="$ALLOWLIST" \
      bash "$TOOLKIT_ROOT/tests/test_no_duplicate_basenames.sh" 2>&1)
rc=$?
set -e

if [[ "$rc" -eq 0 ]]; then
    echo "FAIL [$_TEST_NAME] uniqueness test exited 0 on unlisted collision; expected non-zero" >&2
    echo "--- output ---" >&2
    printf '%s\n' "$out" >&2
    exit 1
fi

# Failure message must (a) name the colliding basename, (b) name both kinds,
# (c) point the contributor at the allowlist.
if ! printf '%s\n' "$out" | grep -q "'s_a'"; then
    echo "FAIL [$_TEST_NAME] failure message missing the colliding name 's_a'" >&2
    printf '%s\n' "$out" >&2
    exit 1
fi
if ! printf '%s\n' "$out" | grep -Eq 'across skill and command'; then
    echo "FAIL [$_TEST_NAME] failure message missing kinds 'skill and command'" >&2
    printf '%s\n' "$out" >&2
    exit 1
fi
if ! printf '%s\n' "$out" | grep -q 'tests/collision-allowlist.txt'; then
    echo "FAIL [$_TEST_NAME] failure message does not point at tests/collision-allowlist.txt" >&2
    printf '%s\n' "$out" >&2
    exit 1
fi

pass
