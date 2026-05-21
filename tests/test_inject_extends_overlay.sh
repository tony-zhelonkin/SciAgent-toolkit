#!/usr/bin/env bash
# Inject a skill on top of base+overlay: extends overlay, no synthetic role.
set -u
. "$(dirname "$0")/_lib.sh"

setup_tmpdir
FAKE="$TMPDIR_TEST/fake-toolkit"
build_fake_toolkit "$FAKE"
export SCIAGENT_TOOLKIT="$FAKE"
SCIAGENT="$FAKE/bin/sciagent"

# Add a fixture skill not in any role, to inject into the real overlay.
mkdir -p "$FAKE/skills/s_inj"
echo "skill s_inj" > "$FAKE/skills/s_inj/SKILL.md"

mkdir project && cd project
"$SCIAGENT" activate base reviewer >/dev/null
"$SCIAGENT" inject s_inj >/dev/null

# Stack unchanged.
assert_grep '^STACK base reviewer$' .sciagent/manifest.json "stack unchanged"

# INJECTED record attributes the skill to the real overlay.
assert_grep '^INJECTED reviewer s_inj$' .sciagent/manifest.json "injected into overlay"

# Symlinks present in both trees.
assert_symlink .claude/skills/s_inj
assert_symlink .agents/skills/s_inj

# Block has Injected subsection and lists the skill.
assert_grep '## Injected (overlay)' AGENTS.md
assert_grep 's_inj' AGENTS.md

pass
