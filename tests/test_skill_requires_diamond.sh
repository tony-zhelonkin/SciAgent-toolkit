#!/usr/bin/env bash
# Diamond: a → {b, c}; b → d; c → d. Activating role with `a` only must
# symlink d exactly once (no duplicate entries).
set -u
. "$(dirname "$0")/_lib.sh"

setup_tmpdir
FAKE="$TMPDIR_TEST/fake-toolkit"
build_fake_toolkit "$FAKE"

mkdir -p "$FAKE/skills/s_d"

cat > "$FAKE/skills/s_a/SKILL.md" <<'EOF'
---
name: s_a
description: diamond top
metadata:
  requires: [s_b, s_c]
---
EOF
cat > "$FAKE/skills/s_b/SKILL.md" <<'EOF'
---
name: s_b
description: diamond left
metadata:
  requires: [s_d]
---
EOF
cat > "$FAKE/skills/s_c/SKILL.md" <<'EOF'
---
name: s_c
description: diamond right
metadata:
  requires: [s_d]
---
EOF
cat > "$FAKE/skills/s_d/SKILL.md" <<'EOF'
---
name: s_d
description: diamond bottom
metadata:
  requires: []
---
EOF
cat > "$FAKE/roles/base.yaml" <<EOF
name: base
description: diamond fixture
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

# Count occurrences of s_d as a symlink target — must be exactly one
# per mount-point (claude + agents).
claude_d_count=$(find .claude/skills -maxdepth 1 -name 's_d' | wc -l)
agents_d_count=$(find .agents/skills -maxdepth 1 -name 's_d' | wc -l)
assert_eq "$claude_d_count" "1" "s_d should appear exactly once under .claude/skills"
assert_eq "$agents_d_count" "1" "s_d should appear exactly once under .agents/skills"

# Manifest should list s_d entries (claude + agents = 2 total, not 4).
sd_manifest_count=$(grep -c '"\.claude/skills/s_d"\|"\.agents/skills/s_d"' .sciagent/manifest.json || true)
assert_eq "$sd_manifest_count" "2" "s_d should appear twice in manifest (claude + agents), not duplicated"

pass
