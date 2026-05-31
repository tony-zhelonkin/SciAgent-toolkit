#!/usr/bin/env bash
# A skill carries a tag not present in tags.yaml. `sciagent validate` must
# exit 1 and name both the offending tag and the offending skill.
set -u
. "$(dirname "$0")/_lib.sh"

setup_tmpdir
FAKE="$TMPDIR_TEST/fake-toolkit"
build_fake_toolkit "$FAKE"

# Plant a skill with an unknown tag. The fake toolkit's tags.yaml only
# declares "tooling" (added by build_fake_toolkit in _lib.sh).
cat > "$FAKE/skills/s_a/SKILL.md" <<'EOF'
---
name: s_a
description: unknown-tag fixture
metadata:
  requires: []
  complementary-skills: []
  contraindications: []
  tags:
  - bogus-tag-xyz
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
    echo "FAIL [$_TEST_NAME] validate returned 0 with unknown tag; expected 1" >&2
    exit 1
fi

# Stderr must name the offending tag.
if ! printf '%s\n' "$err" | grep -q 'bogus-tag-xyz'; then
    echo "FAIL [$_TEST_NAME] error output does not name the offending tag 'bogus-tag-xyz'" >&2
    echo "--- stderr ---" >&2
    printf '%s\n' "$err" >&2
    exit 1
fi

# Stderr must name the offending skill.
if ! printf '%s\n' "$err" | grep -q 's_a'; then
    echo "FAIL [$_TEST_NAME] error output does not name the offending skill 's_a'" >&2
    echo "--- stderr ---" >&2
    printf '%s\n' "$err" >&2
    exit 1
fi

pass
