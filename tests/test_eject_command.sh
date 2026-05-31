#!/usr/bin/env bash
# tests/test_eject_command.sh — round-trip: inject a command then eject it.
# Verifies the command-tree symlinks come back down, the manifest row drops,
# and the _injected overlay collapses when it was the only entry.
set -u
. "$(dirname "$0")/_lib.sh"

setup_tmpdir
FAKE="$TMPDIR_TEST/fake-toolkit"
build_fake_toolkit "$FAKE"
export SCIAGENT_TOOLKIT="$FAKE"
SCIAGENT="$FAKE/bin/sciagent"

mkdir project && cd project
"$SCIAGENT" activate base >/dev/null

# Pre-state: c_b not on disk.
if [[ -L .claude/commands/c_b.md ]]; then
    echo "FAIL [$_TEST_NAME] c_b should not be mounted before inject" >&2
    exit 1
fi

"$SCIAGENT" inject c_b >/dev/null

assert_symlink .claude/commands/c_b.md
assert_symlink .agents/commands/c_b.md
assert_grep '"kind": "command"' .sciagent/manifest.json "manifest carries kind=command before eject"

# Eject.
"$SCIAGENT" eject c_b >/dev/null

# Symlinks gone.
if [[ -L .claude/commands/c_b.md ]]; then
    echo "FAIL [$_TEST_NAME] .claude/commands/c_b.md should be removed after eject" >&2
    exit 1
fi
if [[ -L .agents/commands/c_b.md ]]; then
    echo "FAIL [$_TEST_NAME] .agents/commands/c_b.md should be removed after eject" >&2
    exit 1
fi

# Manifest row gone.
if grep -q '"c_b"' .sciagent/manifest.json; then
    echo "FAIL [$_TEST_NAME] c_b manifest row should be removed" >&2
    exit 1
fi

# _injected overlay collapsed back to [base].
stack_line=$(grep '"stack"' .sciagent/manifest.json)
if printf '%s\n' "$stack_line" | grep -q '"_injected"'; then
    echo "FAIL [$_TEST_NAME] _injected overlay should have collapsed" >&2
    exit 1
fi

pass
