#!/usr/bin/env bash
# Full round-trip: inject s_c; inject --tag pathway; eject s_c; eject --tag pathway
# ends at pre-inject state (manifest injected[] empty, _injected overlay collapsed,
# AGENTS.md has no Injected subsection, all injected symlinks removed).
set -u
. "$(dirname "$0")/_lib.sh"

setup_tmpdir
FAKE="$TMPDIR_TEST/fake-toolkit"
build_fake_toolkit "$FAKE"
export SCIAGENT_TOOLKIT="$FAKE"
SCIAGENT="$FAKE/bin/sciagent"

# Plant a pathway-tagged skill for the --tag leg.
tag_skill "$FAKE" "pw_rtrip" "pathway"

mkdir project && cd project
"$SCIAGENT" activate base >/dev/null

# Capture the manifest hash immediately after activate (pre-inject baseline).
baseline_hash=$(grep '"block_hash"' .sciagent/manifest.json | sed 's/.*"block_hash":[[:space:]]*"\([^"]*\)".*/\1/')

# --- Inject phase ---
"$SCIAGENT" inject s_c >/dev/null
"$SCIAGENT" inject --tag pathway >/dev/null

# Verify both are present.
assert_grep '"s_c"'      .sciagent/manifest.json "s_c injected"
assert_grep '"pw_rtrip"' .sciagent/manifest.json "pw_rtrip injected"
assert_grep '"tag:pathway"' .sciagent/manifest.json "via tag:pathway recorded"

assert_symlink .claude/skills/s_c
assert_symlink .agents/skills/s_c
assert_symlink .claude/skills/pw_rtrip
assert_symlink .agents/skills/pw_rtrip

# --- Eject phase ---
"$SCIAGENT" eject s_c >/dev/null
"$SCIAGENT" eject --tag pathway >/dev/null

# Manifest injected[] must be empty.
if grep -q '"s_c"' .sciagent/manifest.json; then
    echo "FAIL [$_TEST_NAME] s_c should be absent from manifest after eject" >&2
    exit 1
fi
if grep -q '"pw_rtrip"' .sciagent/manifest.json; then
    echo "FAIL [$_TEST_NAME] pw_rtrip should be absent from manifest after eject" >&2
    exit 1
fi

# Injected symlinks removed.
if [[ -L .claude/skills/s_c ]]; then
    echo "FAIL [$_TEST_NAME] .claude/skills/s_c should be removed" >&2
    exit 1
fi
if [[ -L .agents/skills/pw_rtrip ]]; then
    echo "FAIL [$_TEST_NAME] .agents/skills/pw_rtrip should be removed" >&2
    exit 1
fi

# _injected overlay collapsed: stack back to [base].
stack_line=$(grep '"stack"' .sciagent/manifest.json)
if printf '%s\n' "$stack_line" | grep -q '"_injected"'; then
    echo "FAIL [$_TEST_NAME] _injected overlay should be collapsed after all ejects" >&2
    exit 1
fi

# Manifest still valid and still references the base stack.
assert_grep '"base"' .sciagent/manifest.json "base still in stack"

# Base-role symlinks (s_a, s_b) must still be intact.
assert_symlink .claude/skills/s_a
assert_symlink .agents/skills/s_a
assert_symlink .claude/skills/s_b
assert_symlink .agents/skills/s_b

# AGENTS.md should not contain an Injected subsection.
if grep -q '## Injected (overlay)' AGENTS.md 2>/dev/null; then
    echo "FAIL [$_TEST_NAME] AGENTS.md Injected subsection should be gone after full eject" >&2
    exit 1
fi

pass
