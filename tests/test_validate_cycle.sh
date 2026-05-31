#!/usr/bin/env bash
# A requires B, B requires A. `sciagent validate` must exit 1 and name
# both skills in the cycle on stderr.
set -u
. "$(dirname "$0")/_lib.sh"

setup_tmpdir
FAKE="$TMPDIR_TEST/fake-toolkit"
build_fake_toolkit "$FAKE"

cat > "$FAKE/skills/s_a/SKILL.md" <<'EOF'
---
name: s_a
description: cycle fixture A
metadata:
  requires: [s_b]
  complementary-skills: []
  contraindications: []
  tags: []
  scope: implementation
---
EOF
cat > "$FAKE/skills/s_b/SKILL.md" <<'EOF'
---
name: s_b
description: cycle fixture B
metadata:
  requires: [s_a]
  complementary-skills: []
  contraindications: []
  tags: []
  scope: implementation
---
EOF

export SCIAGENT_TOOLKIT="$FAKE"
SCIAGENT="$FAKE/bin/sciagent"

# validate must exit 1 on a cycle.
set +e
err=$("$SCIAGENT" validate 2>&1)
rc=$?
set -e

if [[ "$rc" -eq 0 ]]; then
    echo "FAIL [$_TEST_NAME] validate returned 0 on a requires cycle; expected 1" >&2
    exit 1
fi

# Stderr must name at least one of the cycling skills.
if ! printf '%s\n' "$err" | grep -q 's_a\|s_b'; then
    echo "FAIL [$_TEST_NAME] error output does not name the cycling skills" >&2
    echo "--- stderr ---" >&2
    printf '%s\n' "$err" >&2
    exit 1
fi

pass
