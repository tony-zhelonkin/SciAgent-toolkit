#!/usr/bin/env bash
# tests/test_block_write_io_failure.sh — block_write must not claim success
# when the write failed.
#
# The create branch (file does not exist yet) used to `return 0`
# unconditionally, so a failed redirect printed bash's own "Permission denied"
# to stderr and then reported success: `scio craft --project-dir <read-only
# dir>` said "added SCIO:CRAFT block to: <path>" and exited 0 with no file
# at that path. The append branch in the same function propagated its failure
# all along — only because its printf happens to be the last command — so the
# two write paths disagreed about whether an I/O error is an error.
#
# This is the surviving half of a "craft/provision can leave a 0-byte
# AGENTS.md" report. The 0-byte claim does not reproduce: every craft path
# either writes the complete block or writes nothing at all. What it can do
# — could do — is write nothing and call it a success.
#
# Tests:
#   1. happy path still returns 0 and writes a complete block
#   2. unwritable directory  -> block_write returns 1
#   3. a directory where the file should be -> block_write returns 1
#   4. the same two conditions through `scio craft` -> exit 1, no false
#      "added ... block" line
#   5. an existing but unwritable FILE (the append branch) still returns 1
set -u
. "$(dirname "$0")/_lib.sh"

setup_tmpdir
# shellcheck source=/dev/null
. "$TOOLKIT_ROOT/lib/scio/block.sh"

fail() { echo "FAIL [$_TEST_NAME] $1" >&2; exit 1; }

if [[ "$(id -u)" -eq 0 ]]; then
    echo "SKIP [$_TEST_NAME] running as root: permission bits are not enforced"
    pass
    exit 0
fi

# --- 1. happy path ---------------------------------------------------------
mkdir ok && cd ok
block_write AGENTS.md "hello" CRAFT || fail "happy-path block_write returned non-zero"
assert_exit 0 block_hash_check AGENTS.md CRAFT
[[ -s AGENTS.md ]] || fail "happy-path block_write left an empty AGENTS.md"
cd ..

# --- 2. unwritable directory (create branch) -------------------------------
mkdir ro && chmod 555 ro
set +e
( cd ro && block_write AGENTS.md "hello" CRAFT ) 2>/dev/null
rc=$?
set -e
[[ "$rc" -ne 0 ]] || fail "block_write returned 0 into an unwritable directory"
[[ -e ro/AGENTS.md ]] && fail "a file appeared in an unwritable directory"
chmod 755 ro

# --- 3. a directory sits where the file should be --------------------------
mkdir -p dir/AGENTS.md
set +e
( cd dir && block_write AGENTS.md "hello" CRAFT ) 2>/dev/null
rc=$?
set -e
[[ "$rc" -ne 0 ]] || fail "block_write returned 0 with a directory in the target's place"

# --- 4. through the craft verb --------------------------------------------
SCIO="$TOOLKIT_ROOT/bin/scio"
mkdir ro2 && chmod 555 ro2
set +e
out=$("$SCIO" craft --project-dir "$PWD/ro2" 2>&1); rc=$?
set -e
chmod 755 ro2
[[ "$rc" -eq 1 ]] || fail "scio craft exited $rc against an unwritable dir (expected 1); output: $out"
printf '%s\n' "$out" | grep -q 'added SCIO:CRAFT block' \
    && fail "scio craft reported success after a failed write: $out"
[[ -e ro2/AGENTS.md ]] && fail "craft created a file in an unwritable directory"

# --- 5. append branch: existing, unwritable file ---------------------------
mkdir ap && : > ap/AGENTS.md && chmod 444 ap/AGENTS.md
set +e
( cd ap && block_write AGENTS.md "hello" CRAFT ) 2>/dev/null
rc=$?
set -e
chmod 644 ap/AGENTS.md
[[ "$rc" -ne 0 ]] || fail "block_write returned 0 writing to a read-only file"

pass
