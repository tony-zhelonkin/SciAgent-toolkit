#!/usr/bin/env bash
# Activate base + reviewer; verify last-wins, shadow annotation, all symlinks.
set -u
. "$(dirname "$0")/_lib.sh"

setup_tmpdir
FAKE="$TMPDIR_TEST/fake-toolkit"
build_fake_toolkit "$FAKE"
export SCIAGENT_TOOLKIT="$FAKE"
SCIAGENT="$FAKE/bin/sciagent"

mkdir project && cd project
"$SCIAGENT" activate base reviewer >/dev/null

assert_grep '"base"'     .sciagent/manifest.json "base in stack"
assert_grep '"reviewer"' .sciagent/manifest.json "reviewer in stack"
assert_grep 'overlay'                AGENTS.md             "overlay annotated"

# All effective skills: s_a (base), s_b (reviewer wins, shadows base), s_c (reviewer)
assert_symlink .claude/skills/s_a
assert_symlink .claude/skills/s_b
assert_symlink .claude/skills/s_c

# s_b symlink must point at reviewer's path (same file in fixture, but check the
# managed block records the shadow).
assert_grep 'shadows base' AGENTS.md "s_b shadow recorded"

# Agents/commands from both roles.
assert_symlink .claude/agents/ag_a.md
assert_symlink .claude/agents/ag_b.md
assert_symlink .claude/commands/c_a.md
assert_symlink .claude/commands/c_b.md

pass
