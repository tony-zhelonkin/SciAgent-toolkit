#!/usr/bin/env bash
# tests/test_eject_pre_pr3_manifest.sh — `sciagent eject` tolerates manifest
# entries written by a pre-PR-3 version that lack the "via" field.
#
# Per symlinks.sh comment (2026-05-24 PR 3): entries without "via" are treated
# as named-skill injections (via=="") on read. This test verifies that `eject
# <skill>` successfully removes such an entry.
set -u
. "$(dirname "$0")/_lib.sh"

setup_tmpdir
FAKE="$TMPDIR_TEST/fake-toolkit"
build_fake_toolkit "$FAKE"
export SCIAGENT_TOOLKIT="$FAKE"
SCIAGENT="$FAKE/bin/sciagent"

mkdir project && cd project

# Bootstrap with a normal activate so .sciagent/ and AGENTS.md are present.
"$SCIAGENT" activate base >/dev/null

# Create the symlink pair that the pre-PR-3 inject would have left.
mkdir -p .claude/skills .agents/skills
ln -sfn "$FAKE/skills/s_c" .claude/skills/s_c
ln -sfn "$FAKE/skills/s_c" .agents/skills/s_c

# Rewrite the manifest to include a pre-PR-3 injected entry (no "via" field)
# and record the two new symlink paths.
# We read the current hash so the manifest stays internally consistent.
old_hash=$(grep '"block_hash"' .sciagent/manifest.json \
    | sed 's/.*"block_hash":[[:space:]]*"\([^"]*\)".*/\1/')

# Build a manifest that looks like pre-PR-3 output: the injected entry has
# no "via" field, and the stack already includes the synthetic _injected
# overlay (matching what the old `inject` would have written).
cat > .sciagent/manifest.json <<EOF
{
  "version": 1,
  "stack": ["base", "_injected"],
  "symlinks": [".claude/skills/s_a",".agents/skills/s_a",
               ".claude/skills/s_b",".agents/skills/s_b",
               ".claude/agents/ag_a.md",".agents/agents/ag_a.md",
               ".claude/commands/c_a.md",".agents/commands/c_a.md",
               ".claude/skills/s_c",".agents/skills/s_c"],
  "injected": [
    {"overlay": "_injected", "skill": "s_c"}
  ],
  "block_hash": "$old_hash"
}
EOF

# Verify the pre-condition: symlink is present.
assert_symlink .claude/skills/s_c "pre-condition: .claude/skills/s_c must be present"

# Eject by name — must succeed even though there is no "via" field in the entry.
set +e
out=$("$SCIAGENT" eject s_c 2>&1)
rc=$?
set -e

if [[ "$rc" -ne 0 ]]; then
    echo "FAIL [$_TEST_NAME] eject s_c exited $rc on pre-PR-3 manifest entry; expected 0" >&2
    echo "--- output ---" >&2
    printf '%s\n' "$out" >&2
    exit 1
fi

# Symlinks removed.
if [[ -L .claude/skills/s_c ]]; then
    echo "FAIL [$_TEST_NAME] .claude/skills/s_c still present after eject" >&2
    exit 1
fi
if [[ -L .agents/skills/s_c ]]; then
    echo "FAIL [$_TEST_NAME] .agents/skills/s_c still present after eject" >&2
    exit 1
fi

# Entry removed from manifest injected[].
if grep -q '"s_c"' .sciagent/manifest.json; then
    echo "FAIL [$_TEST_NAME] s_c still present in manifest injected[] after eject" >&2
    exit 1
fi

# Base-role symlinks still intact.
assert_symlink .claude/skills/s_a "s_a must survive eject of s_c"
assert_symlink .agents/skills/s_b "s_b must survive eject of s_c"

pass
