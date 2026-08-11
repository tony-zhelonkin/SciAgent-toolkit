#!/usr/bin/env bash
# Write a block to a file with surrounding content, read it back, hash-check.
set -u
. "$(dirname "$0")/_lib.sh"
. "$TOOLKIT_ROOT/lib/sciagent/block.sh"

setup_tmpdir
cat > AGENTS.md <<'EOF'
# Project AGENTS.md

Some user-owned content above.
EOF

body=$'# Active roles\nstack: base\n'
block_write AGENTS.md "$body" ROLES

read_back=$(block_read AGENTS.md ROLES)
assert_eq "$read_back" "${body%$'\n'}" "block_read returned body"

assert_exit 0 block_hash_check AGENTS.md ROLES

# Surrounding content preserved.
assert_grep '^# Project AGENTS.md' AGENTS.md "header preserved"
assert_grep '^Some user-owned content above\.' AGENTS.md "user content preserved"

pass
