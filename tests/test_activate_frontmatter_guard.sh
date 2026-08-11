#!/usr/bin/env bash
# tests/test_activate_frontmatter_guard.sh
# Pins the hardness boundary defects #1/#2 depended on: `activate`'s
# `cmd_validate --quiet` pre-flight (activate.sh) must still hard-block on
# GENUINE toolkit breakage (a malformed skill frontmatter — the toolkit-wide
# walk that never moved), while a PROJECT-layout finding (docs-layout, now in
# lint.sh) must NOT block it at all (see test_validate_docs_layout.sh test 7
# for that half of the pin). Together these two tests prove the fix drew the
# line in the right place rather than just disabling the pre-flight
# entirely.
set -u
. "$(dirname "$0")/_lib.sh"

setup_tmpdir
FAKE="$TMPDIR_TEST/fake-toolkit"
build_fake_toolkit "$FAKE"

# Break s_a's frontmatter: no `name:` field at all.
cat > "$FAKE/skills/s_a/SKILL.md" <<'EOF'
---
description: A skill with no name field — genuine toolkit breakage.
---

body
EOF

export SCIAGENT_TOOLKIT="$FAKE"
SCIAGENT="$FAKE/bin/sciagent"

mkdir "$TMPDIR_TEST/project" && cd "$TMPDIR_TEST/project"

set +e
out=$("$SCIAGENT" activate base 2>&1)
rc=$?
set -e

if [[ "$rc" -eq 0 ]]; then
    echo "FAIL [$_TEST_NAME] activate must still hard-fail on malformed skill frontmatter" >&2
    printf '%s\n' "$out" >&2
    exit 1
fi
if ! printf '%s\n' "$out" | grep -q 'aborting activation'; then
    echo "FAIL [$_TEST_NAME] expected the 'aborting activation' pre-flight message" >&2
    printf '%s\n' "$out" >&2
    exit 1
fi
if ! printf '%s\n' "$out" | grep -q 's_a.*no name:'; then
    echo "FAIL [$_TEST_NAME] expected validate's own message naming the offending skill" >&2
    printf '%s\n' "$out" >&2
    exit 1
fi

# Nothing must have been mounted — the pre-flight ran before any mutation.
if [[ -e .claude/skills || -e .sciagent/manifest.json ]]; then
    echo "FAIL [$_TEST_NAME] activate aborted but still wrote mounts/manifest" >&2
    exit 1
fi

pass
