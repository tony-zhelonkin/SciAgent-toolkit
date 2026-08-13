#!/usr/bin/env bash
# Either lone marker reports a corrupted block.
set -u
. "$(dirname "$0")/_lib.sh"
. "$TOOLKIT_ROOT/lib/sciagent/block.sh"

setup_tmpdir
cat > AGENTS.md <<'EOF'
# AGENTS
<!-- BEGIN SCIAGENT:CRAFT v1 hash=0000000000000000000000000000000000000000 -->
some body
EOF

assert_exit 2 block_read AGENTS.md CRAFT
assert_exit 2 block_hash_check AGENTS.md CRAFT

cat > AGENTS.md <<'EOF'
# AGENTS
some body
<!-- END SCIAGENT:CRAFT -->
EOF

assert_exit 2 block_read AGENTS.md CRAFT
assert_exit 2 block_hash_check AGENTS.md CRAFT
pass
