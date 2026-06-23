#!/usr/bin/env bash
# tests/test_craft_block_roundtrip.sh
# Verifies that CRAFT blocks are fully independent from ROLES blocks:
#   - write, read, hash-check, drift-detect, and remove all scope to the
#     requested id and never disturb the other block.
set -u
. "$(dirname "$0")/_lib.sh"
. "$TOOLKIT_ROOT/lib/sciagent/block.sh"

setup_tmpdir

# ── 1. Build a file that already has a ROLES block. ──────────────────────────
roles_body=$'# Active roles\nstack: base\n'
block_write AGENTS.md "$roles_body"

assert_exit 0 block_hash_check AGENTS.md          # ROLES hash OK before CRAFT

# ── 2. Append a CRAFT block (second id). ─────────────────────────────────────
craft_body=$'## Figures\n- fig1.png\n'
block_write AGENTS.md "$craft_body" CRAFT

# Both blocks must now be present.
assert_grep 'BEGIN SCIAGENT:ROLES'  AGENTS.md "ROLES BEGIN marker present"
assert_grep 'END SCIAGENT:ROLES'    AGENTS.md "ROLES END marker present"
assert_grep 'BEGIN SCIAGENT:CRAFT'  AGENTS.md "CRAFT BEGIN marker present"
assert_grep 'END SCIAGENT:CRAFT'    AGENTS.md "CRAFT END marker present"

# ── 3. Read-back returns the correct body for each id. ───────────────────────
roles_read=$(block_read AGENTS.md)
craft_read=$(block_read AGENTS.md CRAFT)

assert_eq "$roles_read" "${roles_body%$'\n'}" "block_read (ROLES) returns roles body"
assert_eq "$craft_read" "${craft_body%$'\n'}" "block_read CRAFT returns craft body"

# ── 4. Hash checks pass for both blocks. ─────────────────────────────────────
assert_exit 0 block_hash_check AGENTS.md        # ROLES
assert_exit 0 block_hash_check AGENTS.md CRAFT  # CRAFT

# ── 5. Mutate CRAFT body; CRAFT drifts but ROLES stays clean. ────────────────
sed -i 's/fig1\.png/fig1_MUTATED.png/' AGENTS.md

assert_exit 3 block_hash_check AGENTS.md CRAFT  # drift
assert_exit 0 block_hash_check AGENTS.md        # ROLES unaffected

# ── 6. Rewrite CRAFT to fix drift; then mutate ROLES; ROLES drifts, not CRAFT. ──
block_write AGENTS.md "$craft_body" CRAFT        # restore CRAFT
assert_exit 0 block_hash_check AGENTS.md CRAFT

sed -i 's/stack: base/stack: BASE_MUTATED/' AGENTS.md
assert_exit 3 block_hash_check AGENTS.md        # ROLES drifts
assert_exit 0 block_hash_check AGENTS.md CRAFT  # CRAFT unaffected

# ── 7. block_remove CRAFT leaves ROLES intact. ───────────────────────────────
block_write AGENTS.md "$roles_body"              # restore ROLES before remove
block_remove AGENTS.md CRAFT

assert_exit 1 block_read AGENTS.md CRAFT        # CRAFT gone
assert_exit 0 block_hash_check AGENTS.md        # ROLES still valid

# Confirm CRAFT markers are absent.
if grep -qF 'SCIAGENT:CRAFT' AGENTS.md; then
    echo "FAIL [$_TEST_NAME] CRAFT markers still present after block_remove" >&2
    exit 1
fi

pass
