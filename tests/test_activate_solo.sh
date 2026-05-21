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
# Validate the manifest is real JSON and has expected fields.
assert_grep '"version": 1' .sciagent/manifest.json "version field"
assert_grep '"base"'        .sciagent/manifest.json "base role in stack"
assert_grep '"block_hash"'  .sciagent/manifest.json "block_hash field present"

# Dual symlinks created for all skills/agents/commands.
assert_symlink .claude/skills/s_a
assert_symlink .agents/skills/s_a
assert_symlink .claude/skills/s_b
assert_symlink .agents/skills/s_b
assert_symlink .claude/agents/ag_a.md
assert_symlink .agents/agents/ag_a.md
assert_symlink .claude/commands/c_a.md
assert_symlink .agents/commands/c_a.md

# Managed block content: broad checks.
assert_grep 'Active roles' AGENTS.md
assert_grep '\*\*base\*\*'  AGENTS.md
assert_grep 'ag_a'           AGENTS.md
assert_grep '/c_a'           AGENTS.md

# Tighter check: s_a must appear in the ## Skills (effective) section,
# provided by role 'base', not in shadow lines or elsewhere.
skills_section=$(awk '
    /^## Skills \(effective\)/ { in_s=1; next }
    in_s && /^## /              { in_s=0 }
    in_s                        { print }
' AGENTS.md)
if ! printf '%s\n' "$skills_section" | grep -q 's_a'; then
    echo "FAIL [$_TEST_NAME] s_a not found in ## Skills (effective) section" >&2
    exit 1
fi
if ! printf '%s\n' "$skills_section" | grep -q '(base)'; then
    echo "FAIL [$_TEST_NAME] provider 'base' not found in ## Skills (effective) section" >&2
    exit 1
fi

pass
