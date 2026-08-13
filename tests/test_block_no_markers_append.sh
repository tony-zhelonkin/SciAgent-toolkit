#!/usr/bin/env bash
# No markers in file → block_write appends with a single blank-line separator.
# block_remove later restores exact original bytes.
set -u
. "$(dirname "$0")/_lib.sh"
. "$TOOLKIT_ROOT/lib/scio/block.sh"

setup_tmpdir
cat > AGENTS.md <<'EOF'
# Project

Existing line A.
Existing line B.
EOF
cp AGENTS.md AGENTS.md.orig

assert_exit 1 block_read AGENTS.md CRAFT

block_write AGENTS.md $'craft body\n' CRAFT

# Now contains markers, hash matches.
assert_exit 0 block_hash_check AGENTS.md CRAFT

# Removing restores original bytes.
block_remove AGENTS.md CRAFT
assert_file_eq AGENTS.md AGENTS.md.orig "block_remove restores original bytes"

pass
