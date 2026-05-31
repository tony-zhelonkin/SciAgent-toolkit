#!/usr/bin/env bash
# A skill has a complementary-skills entry that doesn't exist in skills/.
# `sciagent activate` must:
#   - exit 0 (soft-warn, not hard-fail)
#   - emit the STDERR warning summary block naming both skill and missing entry
#   - successfully create symlinks (activation completed)
set -u
. "$(dirname "$0")/_lib.sh"

setup_tmpdir
FAKE="$TMPDIR_TEST/fake-toolkit"
build_fake_toolkit "$FAKE"

# Plant a skill whose complementary-skills references a nonexistent skill.
cat > "$FAKE/skills/s_a/SKILL.md" <<'EOF'
---
name: s_a
description: complementary-miss fixture
metadata:
  requires: []
  complementary-skills:
  - ghost-skill-xyz
  contraindications: []
  tags: []
  scope: implementation
---
EOF

export SCIAGENT_TOOLKIT="$FAKE"
SCIAGENT="$FAKE/bin/sciagent"

mkdir project && cd project

set +e
stderr_out=$("$SCIAGENT" activate base 2>&1 >/dev/null)
rc=$?
set -e

# Must exit 0 — complementary-skills miss is soft-warn only.
if [[ "$rc" -ne 0 ]]; then
    echo "FAIL [$_TEST_NAME] activate exited $rc (expected 0) on missing complementary-skill" >&2
    echo "--- stderr ---" >&2
    printf '%s\n' "$stderr_out" >&2
    exit 1
fi

# Activation must have succeeded: symlink for s_a must exist.
assert_symlink .claude/skills/s_a "symlink absent despite soft-warn activation"

# STDERR summary block must appear.
if ! printf '%s\n' "$stderr_out" | grep -q 'complementary-skills\|warnings'; then
    echo "FAIL [$_TEST_NAME] STDERR summary block absent" >&2
    echo "--- stderr ---" >&2
    printf '%s\n' "$stderr_out" >&2
    exit 1
fi

# Summary must name the missing entry.
if ! printf '%s\n' "$stderr_out" | grep -q 'ghost-skill-xyz'; then
    echo "FAIL [$_TEST_NAME] STDERR summary does not name the missing entry 'ghost-skill-xyz'" >&2
    echo "--- stderr ---" >&2
    printf '%s\n' "$stderr_out" >&2
    exit 1
fi

# Summary must name the skill that carries the reference.
if ! printf '%s\n' "$stderr_out" | grep -q 's_a'; then
    echo "FAIL [$_TEST_NAME] STDERR summary does not name the skill 's_a'" >&2
    echo "--- stderr ---" >&2
    printf '%s\n' "$stderr_out" >&2
    exit 1
fi

pass
