#!/usr/bin/env bash
# a requires a nonexistent skill. Activation must abort with exit 1 and
# leave no filesystem residue.
set -u
. "$(dirname "$0")/_lib.sh"

setup_tmpdir
FAKE="$TMPDIR_TEST/fake-toolkit"
build_fake_toolkit "$FAKE"

cat > "$FAKE/skills/s_a/SKILL.md" <<'EOF'
---
name: s_a
description: dangling-dep fixture
metadata:
  requires: [s_does_not_exist]
---
EOF
cat > "$FAKE/roles/base.yaml" <<EOF
name: base
description: missing-target fixture
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
    echo "FAIL [$_TEST_NAME] activate succeeded with missing requires target" >&2
    exit 1
fi

if [[ -e .sciagent/manifest.json ]]; then
    echo "FAIL [$_TEST_NAME] manifest.json present after missing-dep failed activate" >&2
    exit 1
fi
if [[ -L .claude/skills/s_a ]]; then
    echo "FAIL [$_TEST_NAME] symlink created despite missing requires target" >&2
    exit 1
fi

pass
