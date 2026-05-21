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
assert_grep '^STACK base$' .sciagent/manifest.json "stack reduced to base"

# Overlay-only artifacts gone.
[[ -e .claude/skills/s_c     ]] && { echo "FAIL: s_c (overlay-only) still present"; exit 1; }
[[ -e .claude/agents/ag_b.md ]] && { echo "FAIL: ag_b (overlay-only) still present"; exit 1; }
[[ -e .claude/commands/c_b.md ]] && { echo "FAIL: c_b (overlay-only) still present"; exit 1; }

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
