#!/usr/bin/env bash
# tests/test_eject_agent.sh — round-trip: inject an agent then eject it.
# Verifies the agent-tree symlinks come back down, the manifest row drops,
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

# Pre-state: ag_b not on disk.
if [[ -L .claude/agents/ag_b.md ]]; then
    echo "FAIL [$_TEST_NAME] ag_b should not be mounted before inject" >&2
    exit 1
fi

"$SCIAGENT" inject ag_b >/dev/null

assert_symlink .claude/agents/ag_b.md
assert_symlink .agents/agents/ag_b.md
assert_grep '"kind": "agent"' .sciagent/manifest.json "manifest carries kind=agent before eject"

# Eject.
"$SCIAGENT" eject ag_b >/dev/null

# Symlinks gone.
if [[ -L .claude/agents/ag_b.md ]]; then
    echo "FAIL [$_TEST_NAME] .claude/agents/ag_b.md should be removed after eject" >&2
    exit 1
fi
if [[ -L .agents/agents/ag_b.md ]]; then
    echo "FAIL [$_TEST_NAME] .agents/agents/ag_b.md should be removed after eject" >&2
    exit 1
fi

# Manifest row gone.
if grep -q '"ag_b"' .sciagent/manifest.json; then
    echo "FAIL [$_TEST_NAME] ag_b manifest row should be removed" >&2
    exit 1
fi

# _injected overlay collapsed back to [base].
stack_line=$(grep '"stack"' .sciagent/manifest.json)
if printf '%s\n' "$stack_line" | grep -q '"_injected"'; then
    echo "FAIL [$_TEST_NAME] _injected overlay should have collapsed" >&2
    exit 1
fi

# Idempotency: eject again is a no-op.
out=$("$SCIAGENT" eject ag_b 2>&1)
echo "$out" | grep -qi 'not injected' || {
    echo "FAIL [$_TEST_NAME] expected 'not injected' on repeat eject, got: $out" >&2
    exit 1
}

pass
