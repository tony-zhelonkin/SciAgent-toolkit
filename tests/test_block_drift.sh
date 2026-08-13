#!/usr/bin/env bash
# Edit the body in place; hash_check must return 3.
set -u
. "$(dirname "$0")/_lib.sh"
. "$TOOLKIT_ROOT/lib/scio/block.sh"

setup_tmpdir
printf '# AGENTS\n' > AGENTS.md
block_write AGENTS.md $'line one\nline two\n' CRAFT

assert_exit 0 block_hash_check AGENTS.md CRAFT

# Mutate body inside markers.
sed -i 's/line one/LINE ONE/' AGENTS.md

assert_exit 3 block_hash_check AGENTS.md CRAFT
pass
