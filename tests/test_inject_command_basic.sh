#!/usr/bin/env bash
# tests/test_inject_command_basic.sh — inject a slash command by bare name.
# Verifies that auto-detection routes the name to the commands/ namespace,
# the manifest records kind="command", and the .claude/commands symlink lands
# in the right tree (including resolution through a nested subdir).
set -u
. "$(dirname "$0")/_lib.sh"

setup_tmpdir
FAKE="$TMPDIR_TEST/fake-toolkit"
build_fake_toolkit "$FAKE"
export SCIAGENT_TOOLKIT="$FAKE"
SCIAGENT="$FAKE/bin/sciagent"

# Plant a command in a nested subdir to mirror the real commands/architect/<x>.md
# pattern. Auto-detect must still find it via the recursive walk.
mkdir -p "$FAKE/commands/nested"
echo "cmd nested_cmd" > "$FAKE/commands/nested/nested_cmd.md"

mkdir project && cd project
"$SCIAGENT" activate base >/dev/null

# c_b is the overlay-role command in the fake toolkit — not in base.
"$SCIAGENT" inject c_b >/dev/null

# Symlinks land in the command tree, NOT the skill tree.
assert_symlink .claude/commands/c_b.md "command must mount under .claude/commands/"
assert_symlink .agents/commands/c_b.md "command must mount under .agents/commands/"
if [[ -L .claude/skills/c_b ]]; then
    echo "FAIL [$_TEST_NAME] command inject must not create a skill symlink" >&2
    exit 1
fi

# Manifest records kind=command.
assert_grep '"skill": "c_b"'     .sciagent/manifest.json "manifest names the entry"
assert_grep '"kind": "command"'  .sciagent/manifest.json "manifest carries kind=command"

# Nested command resolves and mounts to the flat .claude/commands/ path.
"$SCIAGENT" inject nested_cmd >/dev/null
assert_symlink .claude/commands/nested_cmd.md "nested command mounts flat"

# Idempotency: re-injecting is a no-op.
out=$("$SCIAGENT" inject c_b 2>&1)
echo "$out" | grep -q 'already injected' || {
    echo "FAIL [$_TEST_NAME] expected idempotency message, got: $out" >&2
    exit 1
}

pass
