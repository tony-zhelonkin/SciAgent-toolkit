#!/usr/bin/env bash
# tests/test_status_shows_inherited.sh — verify Agent 2's inherited-skill
# rendering in AGENTS.md (per ADR-0002 §4.4).
#
# Fixture: orchestrator skill `orch` requires two leaves `leaf_x`, `leaf_y`.
# A role lists only `orch`. After activation, AGENTS.md must:
#   - List `orch` under a heading that does NOT contain the word "inherited"
#   - List `leaf_x` and `leaf_y` under a heading that contains "inherited via requires"
#   - Annotate each inherited leaf with its parent (`(via orch)`)

set -u
. "$(dirname "$0")/_lib.sh"

setup_tmpdir
FAKE="$TMPDIR_TEST/fake-toolkit"
build_fake_toolkit "$FAKE"

# Rewrite fixture skills: orch (orchestrator with requires), leaf_x, leaf_y.
mkdir -p "$FAKE/skills/orch" "$FAKE/skills/leaf_x" "$FAKE/skills/leaf_y"
cat > "$FAKE/skills/orch/SKILL.md" <<'EOF'
---
name: orch
description: orchestrator fixture
metadata:
  scope: orchestrator
  requires:
    - leaf_x
    - leaf_y
---
body
EOF
cat > "$FAKE/skills/leaf_x/SKILL.md" <<'EOF'
---
name: leaf_x
description: leaf x
metadata:
  scope: atomic
  requires: []
---
body
EOF
cat > "$FAKE/skills/leaf_y/SKILL.md" <<'EOF'
---
name: leaf_y
description: leaf y
metadata:
  scope: atomic
  requires: []
---
body
EOF

# Role lists only the orchestrator.
cat > "$FAKE/roles/base.yaml" <<EOF
name: base
description: inherited-render fixture
skills:
  - orch
agents:
  - ag_a
commands:
  - c_a
EOF

export SCIAGENT_TOOLKIT="$FAKE"
SCIAGENT="$FAKE/bin/sciagent"

mkdir project && cd project
"$SCIAGENT" activate base >/dev/null

assert_file_exists AGENTS.md

# 1. Inherited heading present.
assert_grep '## Skills (inherited via requires:)' AGENTS.md \
    "expected '## Skills (inherited via requires:)' heading"

# 2. Direct heading present and does NOT mention 'inherited'.
if ! grep -q '^## Skills (effective)' AGENTS.md; then
    echo "FAIL [$_TEST_NAME] missing '## Skills (effective)' heading" >&2
    cat AGENTS.md >&2
    exit 1
fi

# 3. Direct skill (orch) appears under the effective heading, not inherited.
#    Walk the file: orch must appear BEFORE the inherited section.
orch_line=$(grep -n '`orch`' AGENTS.md | head -1 | cut -d: -f1 || echo 0)
inh_line=$(grep -n 'inherited via requires' AGENTS.md | head -1 | cut -d: -f1 || echo 0)
if [[ "$orch_line" -eq 0 || "$inh_line" -eq 0 || "$orch_line" -ge "$inh_line" ]]; then
    echo "FAIL [$_TEST_NAME] orch should appear under direct heading (before inherited)" >&2
    echo "  orch_line=$orch_line  inh_line=$inh_line" >&2
    cat AGENTS.md >&2
    exit 1
fi

# 4. Each leaf appears under inherited with `(via orch)` annotation.
assert_grep 'leaf_x.*via.*orch' AGENTS.md
assert_grep 'leaf_y.*via.*orch' AGENTS.md

# 5. The inherited leaves must NOT appear in the direct (effective) section.
#    Extract lines strictly between '## Skills (effective)' and the next '## '.
effective_block=$(awk '
    /^## Skills \(effective\)/ { in_block=1; next }
    in_block && /^## / { exit }
    in_block { print }
' AGENTS.md)
if printf '%s\n' "$effective_block" | grep -q '`leaf_x`'; then
    echo "FAIL [$_TEST_NAME] leaf_x leaked into direct/effective section" >&2
    echo "--- effective block ---" >&2
    printf '%s\n' "$effective_block" >&2
    exit 1
fi
if printf '%s\n' "$effective_block" | grep -q '`leaf_y`'; then
    echo "FAIL [$_TEST_NAME] leaf_y leaked into direct/effective section" >&2
    exit 1
fi

pass
