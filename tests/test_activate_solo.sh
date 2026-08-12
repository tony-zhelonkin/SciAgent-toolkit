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
assert_grep '"version": 2' .sciagent/manifest.json "version field (schema v2)"
assert_grep '"base"'        .sciagent/manifest.json "base role in stack"
# block_hash was schema v1's write-only field: activate wrote it, nothing read
# it. Assert its ABSENCE rather than deleting the assertion — the field is gone
# on purpose, and a silent reappearance (someone "restoring" it from an old
# manifest they found on a consumer) is exactly the regression worth catching.
if grep -q '"block_hash"' .sciagent/manifest.json; then
    echo "FAIL [$_TEST_NAME] manifest still carries the removed block_hash field" >&2
    cat .sciagent/manifest.json >&2
    exit 1
fi
# ...and the manifest must still be valid JSON after the key removal (the
# trailing comma on the preceding line has to go with it).
if command -v jq >/dev/null 2>&1; then
    jq -e . .sciagent/manifest.json >/dev/null 2>&1 || {
        echo "FAIL [$_TEST_NAME] manifest.json is not valid JSON" >&2
        cat .sciagent/manifest.json >&2
        exit 1
    }
fi

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
