#!/usr/bin/env bash
# tests/test_collision_allowlist_passes_known.sh
# The CI uniqueness test (test_no_duplicate_basenames.sh) passes when the
# collision IS recorded in the allowlist, with a kinds csv that matches
# the actual overlap exactly.
set -u
. "$(dirname "$0")/_lib.sh"

setup_tmpdir
FAKE="$TMPDIR_TEST/fake-toolkit"
build_fake_toolkit "$FAKE"

# Same fixture as the blocks-unknown test: command s_a vs skill s_a.
echo "cmd s_a (deliberate listed collision)" > "$FAKE/commands/s_a.md"

# Allowlist that documents this exact overlap. Kinds csv must match the
# Catalog collision canonical order (skill,agent,command).
ALLOWLIST="$TMPDIR_TEST/allowlist.txt"
cat > "$ALLOWLIST" <<'EOF'
# fixture allowlist
s_a skill,command # test fixture: the s_a command intentionally shares a name with the s_a skill
EOF

set +e
out=$(COLLISION_TEST_TOOLKIT="$FAKE" \
      COLLISION_TEST_ALLOWLIST="$ALLOWLIST" \
      bash "$TOOLKIT_ROOT/tests/test_no_duplicate_basenames.sh" 2>&1)
rc=$?
set -e

if [[ "$rc" -ne 0 ]]; then
    echo "FAIL [$_TEST_NAME] uniqueness test exited $rc on listed collision (expected 0)" >&2
    echo "--- output ---" >&2
    printf '%s\n' "$out" >&2
    exit 1
fi

# Sanity: the final PASS line must be present (the test should run to
# completion, not exit early before its PASS print).
if ! printf '%s\n' "$out" | grep -q 'PASS \[test_no_duplicate_basenames.sh\]'; then
    echo "FAIL [$_TEST_NAME] uniqueness test did not emit its PASS line" >&2
    printf '%s\n' "$out" >&2
    exit 1
fi

pass
