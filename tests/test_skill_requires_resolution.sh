#!/usr/bin/env bash
# Activate a role whose only skill is `a`; `a` requires `b`. Both must be
# symlinked into .claude/skills and .agents/skills.
set -u
. "$(dirname "$0")/_lib.sh"

setup_tmpdir
FAKE="$TMPDIR_TEST/fake-toolkit"
build_fake_toolkit "$FAKE"

# Overlay the default fixture skills with ones that have requires-bearing
# frontmatter. The fixture roles/base.yaml references s_a + s_b — we point
# s_a at s_b via `requires:` and verify both end up symlinked even though
# the role only mentions s_a (we patch base.yaml below to be a→b only).
cat > "$FAKE/skills/s_a/SKILL.md" <<'EOF'
---
name: s_a
description: requires-bearing fixture skill A
metadata:
  requires:
    - s_b
---
body
EOF
cat > "$FAKE/skills/s_b/SKILL.md" <<'EOF'
---
name: s_b
description: leaf fixture skill B
metadata:
  requires: []
---
body
EOF
# Trim base role to only mention s_a; s_b must be pulled in transitively.
cat > "$FAKE/roles/base.yaml" <<EOF
name: base
description: requires-resolution fixture
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
assert_symlink .agents/skills/s_a
assert_symlink .claude/skills/s_b "transitive requires should symlink s_b"
assert_symlink .agents/skills/s_b "transitive requires should symlink s_b under .agents"

# Block should record s_b under inherited subsection.
assert_grep 'inherited via requires' AGENTS.md
assert_grep 's_b.*via.*s_a' AGENTS.md

pass
