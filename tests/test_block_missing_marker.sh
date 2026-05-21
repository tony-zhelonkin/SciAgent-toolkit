#!/usr/bin/env bash
# Only one marker present → block_read exits 2, block_hash_check exits 2.
set -u
. "$(dirname "$0")/_lib.sh"
. "$TOOLKIT_ROOT/lib/sciagent/block.sh"

setup_tmpdir
cat > AGENTS.md <<'EOF'
# AGENTS
<!-- BEGIN SCIAGENT:ROLES v1 hash=0000000000000000000000000000000000000000 -->
some body
EOF

assert_exit 2 block_read AGENTS.md
assert_exit 2 block_hash_check AGENTS.md
pass
