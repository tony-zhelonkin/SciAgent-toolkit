#!/usr/bin/env bash
# Deep chain a → b → c → d. Activating role with only `a` must pull all
# four into the symlink set.
set -u
. "$(dirname "$0")/_lib.sh"

setup_tmpdir
FAKE="$TMPDIR_TEST/fake-toolkit"
build_fake_toolkit "$FAKE"

# Need an extra fixture skill s_d (build_fake_toolkit only ships a/b/c).
mkdir -p "$FAKE/skills/s_d"

cat > "$FAKE/skills/s_a/SKILL.md" <<'EOF'
---
name: s_a
description: deep chain head
metadata:
  requires: [s_b]
---
EOF
cat > "$FAKE/skills/s_b/SKILL.md" <<'EOF'
---
name: s_b
description: deep chain mid 1
metadata:
  requires: [s_c]
---
EOF
cat > "$FAKE/skills/s_c/SKILL.md" <<'EOF'
---
name: s_c
description: deep chain mid 2
metadata:
  requires: [s_d]
---
EOF
cat > "$FAKE/skills/s_d/SKILL.md" <<'EOF'
---
name: s_d
description: deep chain leaf
metadata:
  requires: []
---
EOF
cat > "$FAKE/roles/base.yaml" <<EOF
name: base
description: deep chain fixture
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
"$SCIAGENT" activate base >/dev/null

assert_symlink .claude/skills/s_a
assert_symlink .claude/skills/s_b
assert_symlink .claude/skills/s_c
assert_symlink .claude/skills/s_d

# Verify the inherited block records all three transitive entries.
assert_grep 's_b' AGENTS.md
assert_grep 's_c' AGENTS.md
assert_grep 's_d' AGENTS.md

pass
