#!/usr/bin/env bash
# Pre-existing user-owned symlink survives activate + deactivate (we only
# touch symlinks recorded in our manifest).
set -u
. "$(dirname "$0")/_lib.sh"

setup_tmpdir
FAKE="$TMPDIR_TEST/fake-toolkit"
build_fake_toolkit "$FAKE"
export SCIAGENT_TOOLKIT="$FAKE"
SCIAGENT="$FAKE/bin/sciagent"

mkdir project && cd project

# Pre-existing user-owned symlink to an unrelated target.
mkdir -p .claude/skills
mkdir -p other
echo "user content" > other/user-skill.md
ln -s "$PWD/other/user-skill.md" .claude/skills/user-owned

# Activate something that DOESN'T collide on name.
"$SCIAGENT" activate alpha >/dev/null

assert_symlink .claude/skills/user-owned "user-owned symlink survives activate"
assert_eq "$(readlink .claude/skills/user-owned)" "$PWD/other/user-skill.md" "target unchanged"

"$SCIAGENT" deactivate >/dev/null

assert_symlink .claude/skills/user-owned "user-owned symlink survives deactivate"
pass
