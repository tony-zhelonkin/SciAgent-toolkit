#!/usr/bin/env bash
# Activate base+reviewer, deactivate the overlay → base remains solo.
set -u
. "$(dirname "$0")/_lib.sh"

setup_tmpdir
FAKE="$TMPDIR_TEST/fake-toolkit"
build_fake_toolkit "$FAKE"
export SCIAGENT_TOOLKIT="$FAKE"
SCIAGENT="$FAKE/bin/sciagent"

mkdir project && cd project
"$SCIAGENT" activate base reviewer >/dev/null
"$SCIAGENT" deactivate reviewer >/dev/null

# Stack should now be just `base`.
assert_grep '"base"' .sciagent/manifest.json "stack reduced to base"
# Ensure reviewer is NOT in stack after partial deactivate.
if grep -q '"reviewer"' .sciagent/manifest.json; then
    echo "FAIL [$_TEST_NAME] reviewer still in manifest after deactivate" >&2
    exit 1
fi

# NOTHING is role-scoped anymore (Phase 5d): the whole catalog — skills,
# agents, commands — mounts regardless of stack, so ag_b/c_b/s_c all survive
# deactivating the overlay that happened to name them. Roles supply
# provenance only, never visibility (see stack.sh:stack_walk).
assert_symlink .claude/agents/ag_b.md
assert_symlink .claude/commands/c_b.md
assert_symlink .claude/skills/s_c

# Base-only artifacts still present.
assert_symlink .claude/skills/s_a
assert_symlink .claude/agents/ag_a.md
assert_symlink .claude/commands/c_a.md
# s_b (in both) — base provides it now.
assert_symlink .claude/skills/s_b

# Removing the base instead implies full teardown.
"$SCIAGENT" deactivate base >/dev/null
[[ -e .sciagent/manifest.json ]] && { echo "FAIL: manifest survived base teardown"; exit 1; }
pass
