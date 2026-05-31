#!/usr/bin/env bash
# tests/test_eject_kind_forward_compat.sh — `sciagent eject` tolerates a
# manifest entry written before the "kind" field existed.
#
# Per symlinks.sh comment (PR-A): rows without a "kind" key are treated as
# "skill" on read — every pre-PR-A injection was a skill, so the legacy
# default carries the migration. This test plants a row that has neither a
# "via" nor a "kind" field (worst-case pre-PR-A shape) and verifies eject
# removes the skill, its symlinks, and the manifest row.
set -u
. "$(dirname "$0")/_lib.sh"

setup_tmpdir
FAKE="$TMPDIR_TEST/fake-toolkit"
build_fake_toolkit "$FAKE"
export SCIAGENT_TOOLKIT="$FAKE"
SCIAGENT="$FAKE/bin/sciagent"

mkdir project && cd project
"$SCIAGENT" activate base >/dev/null

# Create the symlink pair that the pre-PR-A inject would have left.
mkdir -p .claude/skills .agents/skills
ln -sfn "$FAKE/skills/s_c" .claude/skills/s_c
ln -sfn "$FAKE/skills/s_c" .agents/skills/s_c

# Rewrite the manifest with an injected entry that lacks BOTH "via" and "kind"
# fields — the absolute minimum a pre-PR-A version would have written.
old_hash=$(grep '"block_hash"' .sciagent/manifest.json \
    | sed 's/.*"block_hash":[[:space:]]*"\([^"]*\)".*/\1/')

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

assert_symlink .claude/skills/s_c "pre-condition: .claude/skills/s_c must be present"

# Eject by name must succeed even though the row has no "kind" field.
set +e
out=$("$SCIAGENT" eject s_c 2>&1)
rc=$?
set -e

if [[ "$rc" -ne 0 ]]; then
    echo "FAIL [$_TEST_NAME] eject s_c exited $rc on a row missing the 'kind' field; expected 0" >&2
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

# Manifest row removed.
if grep -q '"s_c"' .sciagent/manifest.json; then
    echo "FAIL [$_TEST_NAME] s_c still present in manifest after eject" >&2
    exit 1
fi

# Base symlinks intact.
assert_symlink .claude/skills/s_a "s_a base symlink must survive"

pass
