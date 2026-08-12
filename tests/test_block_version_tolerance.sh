#!/usr/bin/env bash
# Readers locate future marker versions while the writer remains on v1.
set -u
. "$(dirname "$0")/_lib.sh"
. "$TOOLKIT_ROOT/lib/sciagent/block.sh"

setup_tmpdir

body=$'versioned body\n'
printf '# Project\n' > AGENTS.md
block_write AGENTS.md "$body" CRAFT
sed -i 's/SCIAGENT:CRAFT v1 hash=/SCIAGENT:CRAFT v2 hash=/' AGENTS.md
cp AGENTS.md v2.orig

assert_eq "$(block_read AGENTS.md CRAFT)" "${body%$'\n'}" "v2 body located"
assert_eq "$(block_stored_hash AGENTS.md CRAFT)" \
    "$(printf '%s' "$body" | sciagent_sha1_stream)" "v2 stored hash extracted"
assert_eq "$(block_line_range AGENTS.md CRAFT)" "3 5" "v2 marker range located"
assert_exit 0 block_hash_check AGENTS.md CRAFT

assert_exit 1 block_write AGENTS.md $'replacement\n' CRAFT
assert_file_eq AGENTS.md v2.orig "unsupported version preserved on rewrite refusal"

sed -i 's/versioned body/versioned drift/' AGENTS.md
assert_exit 3 block_hash_check AGENTS.md CRAFT

block_remove AGENTS.md CRAFT
printf '# Project\n' > without-block.expected
assert_file_eq AGENTS.md without-block.expected "v2 block removed"

: > v1.actual
block_write v1.actual $'some body\n' ROLES
cat > v1.expected <<'EOF'
<!-- BEGIN SCIAGENT:ROLES v1 hash=41c9105ed1d3fe88123bd360d947c909fb85fd07 -->
some body
<!-- END SCIAGENT:ROLES -->
EOF
assert_file_eq v1.actual v1.expected "v1 write bytes unchanged"

pass
