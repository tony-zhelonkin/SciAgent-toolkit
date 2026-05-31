#!/usr/bin/env bash
# tests/test_inject_ambiguous_hard_fail.sh — bare-name inject of a name that
# resolves in two namespaces must exit 1 with the ambiguity message and leave
# no symlink, no manifest mutation behind.
set -u
. "$(dirname "$0")/_lib.sh"

setup_tmpdir
FAKE="$TMPDIR_TEST/fake-toolkit"
build_fake_toolkit "$FAKE"
export SCIAGENT_TOOLKIT="$FAKE"
SCIAGENT="$FAKE/bin/sciagent"

# Plant a same-named skill + command to create the ambiguity.
mkdir -p "$FAKE/skills/twin"
cat > "$FAKE/skills/twin/SKILL.md" <<'EOF'
---
metadata:
  scope: implementation
  requires: []
  complementary-skills: []
  contraindications: []
  tags: []
---
twin as a skill
EOF
echo "twin as a command" > "$FAKE/commands/twin.md"

mkdir project && cd project
"$SCIAGENT" activate base >/dev/null

# Snapshot the manifest before to verify no mutation occurs on the failed inject.
cp .sciagent/manifest.json /tmp/manifest_before.$$.json

set +e
out=$("$SCIAGENT" inject twin 2>&1)
rc=$?
set -e

if [[ "$rc" -ne 1 ]]; then
    echo "FAIL [$_TEST_NAME] expected exit 1 on ambiguous bare-name inject, got $rc" >&2
    echo "output: $out" >&2
    exit 1
fi

printf '%s\n' "$out" | grep -qi 'ambiguous' || {
    echo "FAIL [$_TEST_NAME] expected message to name 'ambiguous', got: $out" >&2
    exit 1
}
printf '%s\n' "$out" | grep -q 'skill' || {
    echo "FAIL [$_TEST_NAME] expected message to mention 'skill', got: $out" >&2
    exit 1
}
printf '%s\n' "$out" | grep -q 'command' || {
    echo "FAIL [$_TEST_NAME] expected message to mention 'command', got: $out" >&2
    exit 1
}
printf '%s\n' "$out" | grep -q -- '--skill' || {
    echo "FAIL [$_TEST_NAME] expected message to advertise --skill flag, got: $out" >&2
    exit 1
}

# Nothing mounted.
if [[ -L .claude/skills/twin ]]; then
    echo "FAIL [$_TEST_NAME] ambiguous inject must not create the skill symlink" >&2
    exit 1
fi
if [[ -L .claude/commands/twin.md ]]; then
    echo "FAIL [$_TEST_NAME] ambiguous inject must not create the command symlink" >&2
    exit 1
fi

# Manifest unchanged.
if ! diff -q /tmp/manifest_before.$$.json .sciagent/manifest.json >/dev/null; then
    echo "FAIL [$_TEST_NAME] manifest was mutated by a failed inject" >&2
    diff -u /tmp/manifest_before.$$.json .sciagent/manifest.json >&2
    rm -f /tmp/manifest_before.$$.json
    exit 1
fi
rm -f /tmp/manifest_before.$$.json

pass
