#!/usr/bin/env bash
# tests/test_block_preserves_mode.sh — a managed-block rewrite must not change
# who can read AGENTS.md.
#
# block_write's in-place rewrite path and block_remove both wrote to mktemp
# (0600) and `mv`d it over the target, carrying 0600 onto a file that was 0644.
# That is not theoretical: it already happened across eleven analysis repos,
# including a shared lab tree where AGENTS.md ended up -rw------- while its own
# CLAUDE.md and README.md siblings stayed -rw-r--r--, so colleagues could no
# longer read it. AGENTS.md exists to be read by other people and by other
# tools; rendering a block into it has no business narrowing its permissions.
#
# The append path (no markers yet) always preserved mode — it appends with >>
# rather than replacing the inode — so this test pins the two paths that did
# not: rewrite-in-place, and remove.
set -u
. "$(dirname "$0")/_lib.sh"

setup_tmpdir
. "$TOOLKIT_ROOT/lib/scio/block.sh"

# ---------------------------------------------------------------------------
# 1. block_write, append path (no existing block) — preserves mode.
# ---------------------------------------------------------------------------
printf '# Project\n\nprose\n' > a.md
chmod 664 a.md
block_write a.md "body line
" CRAFT
assert_eq "$(stat -c '%a' a.md)" "664" "append path preserves 0664"

# ---------------------------------------------------------------------------
# 2. block_write, in-place rewrite path (block already present) — the
#    regression. A second write with different content takes the awk+mv
#    branch, which is where 0600 leaked in.
# ---------------------------------------------------------------------------
block_write a.md "different body
" CRAFT
assert_eq "$(stat -c '%a' a.md)" "664" "rewrite path preserves 0664"

# Non-default, non-0600 mode so the assertion cannot pass by coincidence with
# either the original mode or mktemp's.
chmod 640 a.md
block_write a.md "third body
" CRAFT
assert_eq "$(stat -c '%a' a.md)" "640" "rewrite path preserves 0640"

# Group-writable + setgid-adjacent bits survive too (shared lab trees use them).
chmod 664 a.md
block_write a.md "fourth body
" CRAFT
assert_eq "$(stat -c '%a' a.md)" "664" "rewrite path preserves group write"

# ---------------------------------------------------------------------------
# 3. block_remove — same mv, same flaw.
# ---------------------------------------------------------------------------
chmod 644 a.md
block_remove a.md CRAFT
assert_eq "$(stat -c '%a' a.md)" "644" "remove path preserves 0644"
grep -q 'BEGIN SCIO:CRAFT' a.md \
    && { echo "FAIL [$_TEST_NAME] block_remove did not remove the block" >&2; exit 1; }

# ---------------------------------------------------------------------------
# 4. Content is still correct — a mode-preserving mv must not corrupt the
#    rewrite it is preserving the mode of.
# ---------------------------------------------------------------------------
printf '# Project\n\nprose\n' > b.md
chmod 664 b.md
block_write b.md "first
" CRAFT
block_write b.md "second
" CRAFT
assert_eq "$(block_read b.md CRAFT)" "second" "rewrite still replaces the body"
block_hash_check b.md CRAFT
assert_eq "$?" "0" "stored hash still matches after a mode-preserving rewrite"
grep -q '^prose$' b.md || { echo "FAIL [$_TEST_NAME] prose outside markers lost" >&2; exit 1; }

# ---------------------------------------------------------------------------
# 5. A missing target must not break the write (mode capture is best-effort).
# ---------------------------------------------------------------------------
rm -f c.md
block_write c.md "fresh
" CRAFT
assert_file_exists c.md "block_write still creates a file that did not exist"
assert_eq "$(block_read c.md CRAFT)" "fresh" "new-file body correct"

pass
