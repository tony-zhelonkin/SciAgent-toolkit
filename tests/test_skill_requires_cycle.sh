#!/usr/bin/env bash
# a → b → a. Activation must abort with exit 1, and crucially must NOT have
# created any symlinks or written .sciagent/manifest.json (pre-mutation fail).
set -u
. "$(dirname "$0")/_lib.sh"

setup_tmpdir
FAKE="$TMPDIR_TEST/fake-toolkit"
build_fake_toolkit "$FAKE"

cat > "$FAKE/skills/s_a/SKILL.md" <<'EOF'
---
name: s_a
description: cyclic A
metadata:
  requires: [s_b]
---
EOF
cat > "$FAKE/skills/s_b/SKILL.md" <<'EOF'
---
name: s_b
description: cyclic B
metadata:
  requires: [s_a]
---
EOF
cat > "$FAKE/roles/base.yaml" <<EOF
name: base
description: cycle fixture
skills:
  - s_a
agents:
  - ag_a
commands:
  - c_a
EOF

export SCIAGENT_TOOLKIT="$FAKE"
SCIAGENT="$FAKE/bin/sciagent"

mkdir project && cd project
if "$SCIAGENT" activate base >/dev/null 2>&1; then
    echo "FAIL [$_TEST_NAME] activate succeeded on cyclic requires" >&2
    exit 1
fi

# Pre-mutation contract: nothing landed.
if [[ -e .sciagent/manifest.json ]]; then
    echo "FAIL [$_TEST_NAME] manifest.json present after cycle-failed activate" >&2
    exit 1
fi
if [[ -L .claude/skills/s_a || -L .claude/skills/s_b ]]; then
    echo "FAIL [$_TEST_NAME] symlinks present after cycle-failed activate" >&2
    exit 1
fi

pass
