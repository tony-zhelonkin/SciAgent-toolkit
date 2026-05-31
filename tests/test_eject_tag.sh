#!/usr/bin/env bash
# inject --tag pathway then eject --tag pathway:
# only pathway-mounted skills are removed; other injected entries remain.
set -u
. "$(dirname "$0")/_lib.sh"

setup_tmpdir
FAKE="$TMPDIR_TEST/fake-toolkit"
build_fake_toolkit "$FAKE"
export SCIAGENT_TOOLKIT="$FAKE"
SCIAGENT="$FAKE/bin/sciagent"

# Plant two pathway skills and one extra unrelated skill to inject manually.
tag_skill "$FAKE" "pw_alpha" "pathway"
tag_skill "$FAKE" "pw_beta"  "pathway"

mkdir project && cd project
"$SCIAGENT" activate base >/dev/null

# Inject by tag (two pathway skills) and one named skill (s_c — no tag).
"$SCIAGENT" inject --tag pathway >/dev/null
"$SCIAGENT" inject s_c >/dev/null

# Both pathway skills and s_c should be in the manifest.
assert_grep '"pw_alpha"' .sciagent/manifest.json "pw_alpha in manifest"
assert_grep '"pw_beta"'  .sciagent/manifest.json "pw_beta in manifest"
assert_grep '"s_c"'      .sciagent/manifest.json "s_c in manifest"

# Eject by tag.
"$SCIAGENT" eject --tag pathway >/dev/null

# Pathway skills gone.
if grep -q '"pw_alpha"' .sciagent/manifest.json; then
    echo "FAIL [$_TEST_NAME] pw_alpha should be ejected" >&2
    exit 1
fi
if grep -q '"pw_beta"' .sciagent/manifest.json; then
    echo "FAIL [$_TEST_NAME] pw_beta should be ejected" >&2
    exit 1
fi
if [[ -L .claude/skills/pw_alpha ]]; then
    echo "FAIL [$_TEST_NAME] .claude/skills/pw_alpha symlink should be removed" >&2
    exit 1
fi
if [[ -L .agents/skills/pw_beta ]]; then
    echo "FAIL [$_TEST_NAME] .agents/skills/pw_beta symlink should be removed" >&2
    exit 1
fi

# s_c (named-skill injection) must remain.
assert_grep '"s_c"' .sciagent/manifest.json "s_c still in manifest after tag-eject"
assert_symlink .claude/skills/s_c
assert_symlink .agents/skills/s_c

# _injected overlay must NOT collapse (s_c still there).
stack_line=$(grep '"stack"' .sciagent/manifest.json)
if ! printf '%s\n' "$stack_line" | grep -q '"_injected"'; then
    echo "FAIL [$_TEST_NAME] _injected overlay should persist while s_c remains" >&2
    exit 1
fi

pass
