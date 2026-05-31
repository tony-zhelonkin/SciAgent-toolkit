#!/usr/bin/env bash
# A skill requires a nonexistent skill. `sciagent validate` must exit 1.
set -u
. "$(dirname "$0")/_lib.sh"

setup_tmpdir
FAKE="$TMPDIR_TEST/fake-toolkit"
build_fake_toolkit "$FAKE"

cat > "$FAKE/skills/s_a/SKILL.md" <<'EOF'
---
name: s_a
description: missing-dep fixture
metadata:
  requires: [does-not-exist]
  complementary-skills: []
  contraindications: []
  tags: []
  scope: implementation
---
EOF

export SCIAGENT_TOOLKIT="$FAKE"
SCIAGENT="$FAKE/bin/sciagent"

set +e
err=$("$SCIAGENT" validate 2>&1)
rc=$?
set -e

if [[ "$rc" -eq 0 ]]; then
    echo "FAIL [$_TEST_NAME] validate returned 0 with missing requires target; expected 1" >&2
    exit 1
fi

# Stderr must mention the missing skill name.
if ! printf '%s\n' "$err" | grep -q 'does-not-exist\|s_a'; then
    echo "FAIL [$_TEST_NAME] error output does not name the missing skill or its dependent" >&2
    echo "--- stderr ---" >&2
    printf '%s\n' "$err" >&2
    exit 1
fi

pass
