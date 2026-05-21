#!/usr/bin/env bash
# tests/test_list_deps_verbs.sh — cover `sciagent list deps` and
# `sciagent list dependents` against a synthetic skill graph.
#
# Graph:
#   top requires: mid_a, leaf_z
#   mid_a requires: leaf_x, leaf_y
#   leaf_x, leaf_y, leaf_z: no deps
#   sibling: requires leaf_x (used to exercise `dependents`)

set -u
. "$(dirname "$0")/_lib.sh"

setup_tmpdir
FAKE="$TMPDIR_TEST/fake-toolkit"
build_fake_toolkit "$FAKE"

mkdir -p "$FAKE/skills/top" "$FAKE/skills/mid_a" "$FAKE/skills/leaf_x" \
         "$FAKE/skills/leaf_y" "$FAKE/skills/leaf_z" "$FAKE/skills/sibling"

cat > "$FAKE/skills/top/SKILL.md" <<'EOF'
---
name: top
description: top orchestrator
metadata:
  scope: orchestrator
  requires:
    - mid_a
    - leaf_z
---
body
EOF
cat > "$FAKE/skills/mid_a/SKILL.md" <<'EOF'
---
name: mid_a
description: middle node
metadata:
  scope: orchestrator
  requires:
    - leaf_x
    - leaf_y
---
body
EOF
for leaf in leaf_x leaf_y leaf_z; do
    cat > "$FAKE/skills/$leaf/SKILL.md" <<EOF
---
name: $leaf
description: $leaf
metadata:
  scope: atomic
  requires: []
---
body
EOF
done
cat > "$FAKE/skills/sibling/SKILL.md" <<'EOF'
---
name: sibling
description: another consumer of leaf_x
metadata:
  scope: atomic
  requires:
    - leaf_x
---
body
EOF

export SCIAGENT_TOOLKIT="$FAKE"
SCIAGENT="$FAKE/bin/sciagent"

# --- list deps top -----------------------------------------------------------
# Expected (post-order, root excluded): leaf_x, leaf_y, mid_a, leaf_z
deps_out=$("$SCIAGENT" list deps top)
expected_deps=$'leaf_x\nleaf_y\nmid_a\nleaf_z'
assert_eq "$deps_out" "$expected_deps" "list deps top"

# --- list deps mid_a ---------------------------------------------------------
deps_mid=$("$SCIAGENT" list deps mid_a)
expected_mid=$'leaf_x\nleaf_y'
assert_eq "$deps_mid" "$expected_mid" "list deps mid_a"

# --- list deps on a leaf prints nothing -------------------------------------
deps_leaf=$("$SCIAGENT" list deps leaf_x)
assert_eq "$deps_leaf" "" "list deps leaf_x is empty"

# --- list dependents leaf_x --------------------------------------------------
# Direct dependents of leaf_x: mid_a, sibling (NOT top — top is two hops away)
dep_out=$("$SCIAGENT" list dependents leaf_x)
expected_dep=$'mid_a\nsibling'
assert_eq "$dep_out" "$expected_dep" "list dependents leaf_x"

# --- list dependents mid_a -------------------------------------------------
dep_mid=$("$SCIAGENT" list dependents mid_a)
assert_eq "$dep_mid" "top" "list dependents mid_a"

# --- list dependents on a leaf with no consumers ----------------------------
dep_none=$("$SCIAGENT" list dependents leaf_y)
assert_eq "$dep_none" "mid_a" "list dependents leaf_y"

# --- Missing argument errors ------------------------------------------------
if "$SCIAGENT" list deps >/dev/null 2>&1; then
    echo "FAIL [$_TEST_NAME] 'list deps' without arg should error" >&2
    exit 1
fi
if "$SCIAGENT" list dependents >/dev/null 2>&1; then
    echo "FAIL [$_TEST_NAME] 'list dependents' without arg should error" >&2
    exit 1
fi

pass
