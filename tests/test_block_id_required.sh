#!/usr/bin/env bash
# tests/test_block_id_required.sh
# block.sh requires an explicit id. A missing/empty id
# must be a DISTINCT return code (4) from every function's own "no block"/
# "no markers" result code (1) — otherwise a caller bug (forgot to pass the
# id) is indistinguishable from a legitimate "there is no block here yet",
# and code that branches on rc=1 to mean "safe to write" would silently take
# that branch on a bug instead. This test fixture is built to actually catch
# a collision: it calls the SAME function against the SAME file twice — once
# genuinely id-less (file has no block at all) and once with a missing id
# argument — and asserts the two return codes differ.
set -u
. "$(dirname "$0")/_lib.sh"
. "$TOOLKIT_ROOT/lib/sciagent/block.sh"

setup_tmpdir
printf '# AGENTS\n\nno block here\n' > AGENTS.md

# --- block_read: rc=1 ("no markers") must differ from rc=4 (missing id) ---
assert_exit 1 block_read AGENTS.md CRAFT     # legitimate: no markers present
assert_exit 4 block_read AGENTS.md ""        # caller bug: empty id
assert_exit 4 block_read AGENTS.md           # caller bug: id omitted entirely

# --- block_hash_check: rc=1 ("no block") must differ from rc=4 (missing id) ---
assert_exit 1 block_hash_check AGENTS.md CRAFT
assert_exit 4 block_hash_check AGENTS.md ""
assert_exit 4 block_hash_check AGENTS.md

# --- Every other public function also rejects a missing id with rc=4, and
# does not touch the filesystem when it does (AGENTS.md unchanged). ---
cp AGENTS.md AGENTS.md.orig

assert_exit 4 block_write AGENTS.md "some body"
assert_exit 4 block_remove AGENTS.md
assert_exit 4 block_line_range AGENTS.md
out=$(block_stored_hash AGENTS.md 2>/dev/null); rc=$?
assert_eq "$rc" "4" "block_stored_hash missing id"
assert_eq "$out" "" "block_stored_hash missing id prints nothing"

assert_file_eq AGENTS.md AGENTS.md.orig "no public function touches the file when id is missing"

# The writer accepts only the live CRAFT id.
assert_exit 4 block_write AGENTS.md "retired body" ROLES
assert_exit 0 block_write AGENTS.md "real body" CRAFT
assert_exit 0 block_hash_check AGENTS.md CRAFT

pass
