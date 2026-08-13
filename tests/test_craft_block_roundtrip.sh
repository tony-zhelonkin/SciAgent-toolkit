#!/usr/bin/env bash
# CRAFT writes preserve a stale ROLES block, which remains readable/removable.
set -u
. "$(dirname "$0")/_lib.sh"
. "$TOOLKIT_ROOT/lib/scio/block.sh"

setup_tmpdir

cat > AGENTS.md <<'EOF'
# Project

<!-- BEGIN SCIAGENT:ROLES v1 hash=0ec4c1c79b42a1256dd537137412717ca94ede07 -->
# Active roles
stack: base
<!-- END SCIAGENT:ROLES -->
EOF

assert_eq "$(block_read AGENTS.md ROLES)" $'# Active roles\nstack: base' \
    "stale ROLES body is readable"

craft_body=$'## Figures\n- fig1.png\n'
block_write AGENTS.md "$craft_body" CRAFT
assert_exit 0 block_hash_check AGENTS.md CRAFT
assert_grep 'BEGIN SCIAGENT:ROLES' AGENTS.md "CRAFT write preserves legacy ROLES"

block_remove AGENTS.md CRAFT
assert_exit 1 block_read AGENTS.md CRAFT
block_read AGENTS.md ROLES >/dev/null || {
    echo "FAIL [$_TEST_NAME] stale ROLES block became unreadable" >&2
    exit 1
}

block_remove AGENTS.md ROLES
assert_exit 1 block_read AGENTS.md ROLES
assert_grep '^# Project$' AGENTS.md "project prose survives ROLES removal"

cat > current-roles.md <<'EOF'
<!-- BEGIN SCIO:ROLES v1 hash=0ec4c1c79b42a1256dd537137412717ca94ede07 -->
# Active roles
stack: base
<!-- END SCIO:ROLES -->
EOF
assert_eq "$(block_read current-roles.md ROLES)" $'# Active roles\nstack: base' \
    "SCIO ROLES body is readable"
block_remove current-roles.md ROLES
assert_exit 1 block_read current-roles.md ROLES

pass
