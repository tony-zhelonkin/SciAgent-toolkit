#!/usr/bin/env bash
# inject --tag pathway on solo base: synthesizes _injected overlay and mounts
# all pathway-tagged skills; manifest entries carry via:"tag:pathway".
set -u
. "$(dirname "$0")/_lib.sh"

setup_tmpdir
FAKE="$TMPDIR_TEST/fake-toolkit"
build_fake_toolkit "$FAKE"
export SCIAGENT_TOOLKIT="$FAKE"
SCIAGENT="$FAKE/bin/sciagent"

# Plant two pathway-tagged skills and one untagged skill.
tag_skill "$FAKE" "pw_one" "pathway"
tag_skill "$FAKE" "pw_two" "pathway"
mkdir -p "$FAKE/skills/untagged"
echo "skill untagged" > "$FAKE/skills/untagged/SKILL.md"

mkdir project && cd project
"$SCIAGENT" activate base >/dev/null

# Inject by tag.
"$SCIAGENT" inject --tag pathway >/dev/null

# Stack rewritten to include synthetic _injected overlay.
assert_grep '"_injected"' .sciagent/manifest.json "synthetic overlay recorded"

# Both pathway skills injected.
assert_grep '"pw_one"'  .sciagent/manifest.json "pw_one present"
assert_grep '"pw_two"'  .sciagent/manifest.json "pw_two present"

# Manifest entries carry via: "tag:pathway".
assert_grep '"tag:pathway"' .sciagent/manifest.json "via tag:pathway recorded"

# Untagged skill NOT injected.
if grep -q '"untagged"' .sciagent/manifest.json; then
    echo "FAIL [$_TEST_NAME] untagged skill should not be in manifest" >&2
    exit 1
fi

# Dual symlinks created for each pathway skill.
assert_symlink .claude/skills/pw_one
assert_symlink .agents/skills/pw_one
assert_symlink .claude/skills/pw_two
assert_symlink .agents/skills/pw_two

# AGENTS.md Injected subsection lists both.
assert_grep '## Injected (overlay)' AGENTS.md "Injected subsection present"
assert_grep 'pw_one' AGENTS.md "pw_one listed in AGENTS.md"
assert_grep 'pw_two' AGENTS.md "pw_two listed in AGENTS.md"

# Idempotency: re-running inject --tag pathway is a no-op.
out=$("$SCIAGENT" inject --tag pathway 2>&1)
if printf '%s\n' "$out" | grep -qv 'already injected'; then
    # At least some lines should say "already injected".
    if ! printf '%s\n' "$out" | grep -q 'already injected'; then
        echo "FAIL [$_TEST_NAME] expected idempotency messages, got: $out" >&2
        exit 1
    fi
fi

# Empty-tag warn: a valid tag with no matching skills → STDERR warn, exit 0.
# 'tooling' exists in the fake tags.yaml but no skills carry it.
out=$("$SCIAGENT" inject --tag tooling 2>&1)
echo "$out" | grep -qi 'nothing to inject' || {
    echo "FAIL [$_TEST_NAME] expected empty-tag warn, got: $out" >&2
    exit 1
}

pass
