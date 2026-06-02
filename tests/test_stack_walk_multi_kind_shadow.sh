#!/usr/bin/env bash
# tests/test_stack_walk_multi_kind_shadow.sh — locks the _sw_record (5.1)
# behavior: an overlay that re-declares a base skill, agent, AND command
# simultaneously must produce the correct <shadowed-roles-csv> column for all
# three kinds in the stack_walk TSV.
set -u
. "$(dirname "$0")/_lib.sh"

setup_tmpdir
FAKE="$TMPDIR_TEST/fake-toolkit"
build_fake_toolkit "$FAKE"

# Overlay role that shadows base's s_a (skill), ag_a (agent), c_a (command).
cat > "$FAKE/roles/shadower.yaml" <<EOF
name: shadower
description: re-declares base entries to force a shadow across three kinds
skills:
  - s_a
agents:
  - ag_a
commands:
  - c_a
EOF

export SCIAGENT_TOOLKIT="$FAKE"

# Drive stack_walk directly: source the modules the way bin/sciagent does.
LIB="$FAKE/lib/sciagent"
# shellcheck source=/dev/null
. "$LIB/roles.sh"
. "$LIB/block.sh"
. "$LIB/stack.sh"

OUT=$(stack_walk base shadower)

# s_a: provider=shadower, shadows base.
line=$(printf '%s\n' "$OUT" | grep -P '^SKILL\ts_a\t')
assert_eq "$line" "$(printf 'SKILL\ts_a\tshadower\tbase')" "s_a shadow csv"

# ag_a: provider=shadower, shadows base.
line=$(printf '%s\n' "$OUT" | grep -P '^AGENT\tag_a\t')
assert_eq "$line" "$(printf 'AGENT\tag_a\tshadower\tbase')" "ag_a shadow csv"

# c_a: provider=shadower, shadows base.
line=$(printf '%s\n' "$OUT" | grep -P '^COMMAND\tc_a\t')
assert_eq "$line" "$(printf 'COMMAND\tc_a\tshadower\tbase')" "c_a shadow csv"

pass
