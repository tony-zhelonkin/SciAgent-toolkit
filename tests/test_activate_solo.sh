#!/usr/bin/env bash
# Activate a solo role; verify dual symlinks and AGENTS.md block.
set -u
. "$(dirname "$0")/_lib.sh"

setup_tmpdir
FAKE="$TMPDIR_TEST/fake-toolkit"
build_fake_toolkit "$FAKE"
export SCIAGENT_TOOLKIT="$FAKE"
SCIAGENT="$FAKE/bin/sciagent"

mkdir project && cd project
"$SCIAGENT" activate base >/dev/null

assert_file_exists .sciagent/manifest.json
assert_grep '^VERSION 1$' .sciagent/manifest.json
assert_grep '^STACK base$' .sciagent/manifest.json
assert_grep '^BLOCK_HASH ' .sciagent/manifest.json

# Dual symlinks created for all skills/agents/commands.
assert_symlink .claude/skills/s_a
assert_symlink .agents/skills/s_a
assert_symlink .claude/skills/s_b
assert_symlink .agents/skills/s_b
assert_symlink .claude/agents/ag_a.md
assert_symlink .agents/agents/ag_a.md
assert_symlink .claude/commands/c_a.md
assert_symlink .agents/commands/c_a.md

# Managed block content.
assert_grep 'Active roles' AGENTS.md
assert_grep '\*\*base\*\*'  AGENTS.md
assert_grep 's_a'            AGENTS.md
assert_grep 'ag_a'           AGENTS.md
assert_grep '/c_a'           AGENTS.md

pass
