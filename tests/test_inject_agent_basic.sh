#!/usr/bin/env bash
# tests/test_inject_agent_basic.sh — inject a sub-agent by bare name.
# Verifies that auto-detection routes the name to the agents/ namespace,
# the manifest records kind="agent", and the .claude/agents symlink lands
# in the right tree.
set -u
. "$(dirname "$0")/_lib.sh"

setup_tmpdir
FAKE="$TMPDIR_TEST/fake-toolkit"
build_fake_toolkit "$FAKE"
export SCIAGENT_TOOLKIT="$FAKE"
SCIAGENT="$FAKE/bin/sciagent"

mkdir project && cd project
"$SCIAGENT" activate base >/dev/null

# ag_b is the overlay-role agent in the fake toolkit — not in base.
"$SCIAGENT" inject ag_b >/dev/null

# Symlinks land in the agent tree, NOT the skill tree.
assert_symlink .claude/agents/ag_b.md "agent must mount under .claude/agents/"
assert_symlink .agents/agents/ag_b.md "agent must mount under .agents/agents/"
if [[ -L .claude/skills/ag_b ]]; then
    echo "FAIL [$_TEST_NAME] agent inject must not create a skill symlink" >&2
    exit 1
fi

# Manifest records kind=agent.
assert_grep '"skill": "ag_b"'   .sciagent/manifest.json "manifest names the entry"
assert_grep '"kind": "agent"'   .sciagent/manifest.json "manifest carries kind=agent"

# Stack synthesises the _injected overlay (no overlay was active before).
assert_grep '"_injected"' .sciagent/manifest.json "_injected overlay present"

# Idempotency: re-injecting is a no-op (exit 0).
out=$("$SCIAGENT" inject ag_b 2>&1)
echo "$out" | grep -q 'already injected' || {
    echo "FAIL [$_TEST_NAME] expected idempotency message, got: $out" >&2
    exit 1
}

pass
