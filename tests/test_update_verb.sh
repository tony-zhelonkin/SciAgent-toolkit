#!/usr/bin/env bash
# tests/test_update_verb.sh
# Exercises `sciagent update`:
#   1. --no-pin path: re-activates the current stack; ROLES + CRAFT blocks
#      are present and hash-valid; manifest still records the same stack.
#   2. No active stack: returns non-zero with the guidance message.
set -u
. "$(dirname "$0")/_lib.sh"

setup_tmpdir
FAKE="$TMPDIR_TEST/fake-toolkit"
build_fake_toolkit "$FAKE"

# Provision a craft.yaml so the CRAFT block is rendered during activate.
cat > "$FAKE/craft.yaml" <<'EOF'
version: 3

floors:
  figure_base_size: 16

body: |
  # Craft standards
  - Figures: base >= {{figure_base_size}}pt.
EOF

export SCIAGENT_TOOLKIT="$FAKE"
SCIAGENT="$FAKE/bin/sciagent"

# -------------------------------------------------------------------------
# Test A: --no-pin path re-activates and produces valid blocks.
# -------------------------------------------------------------------------
mkdir project && cd project

# Activate a base stack first.
"$SCIAGENT" activate base >/dev/null

# Capture the stack recorded in the manifest before update.
. "$FAKE/lib/sciagent/symlinks.sh"
stack_before=$(manifest_stack)

# Run update --no-pin; should succeed.
update_out=$("$SCIAGENT" update --no-pin 2>&1)
update_rc=$?

if [[ "$update_rc" -ne 0 ]]; then
    echo "FAIL [$_TEST_NAME] sciagent update --no-pin exited $update_rc" >&2
    echo "output: $update_out" >&2
    exit 1
fi

# ROLES block must be present and hash-valid.
. "$FAKE/lib/sciagent/block.sh"
if ! block_hash_check AGENTS.md; then
    echo "FAIL [$_TEST_NAME] ROLES block hash invalid after update --no-pin" >&2
    exit 1
fi

# CRAFT block must be present and hash-valid.
if ! block_hash_check AGENTS.md CRAFT; then
    echo "FAIL [$_TEST_NAME] CRAFT block hash invalid after update --no-pin" >&2
    exit 1
fi

# Manifest must still record the same stack.
stack_after=$(manifest_stack)
if [[ "$stack_before" != "$stack_after" ]]; then
    echo "FAIL [$_TEST_NAME] manifest stack changed: before='$stack_before' after='$stack_after'" >&2
    exit 1
fi

# Summary output should mention the stack and CRAFT version.
if ! printf '%s\n' "$update_out" | grep -q "stack:"; then
    echo "FAIL [$_TEST_NAME] update summary missing 'stack:' line" >&2
    echo "output: $update_out" >&2
    exit 1
fi
if ! printf '%s\n' "$update_out" | grep -q "CRAFT ver:"; then
    echo "FAIL [$_TEST_NAME] update summary missing 'CRAFT ver:' line" >&2
    echo "output: $update_out" >&2
    exit 1
fi

# Submodule line should indicate --no-pin skipped.
if ! printf '%s\n' "$update_out" | grep -q "skipped"; then
    echo "FAIL [$_TEST_NAME] update --no-pin did not print 'skipped' for submodule step" >&2
    echo "output: $update_out" >&2
    exit 1
fi

# -------------------------------------------------------------------------
# Test B: --no-pin with two-role stack (base + reviewer).
# -------------------------------------------------------------------------
"$SCIAGENT" activate base reviewer >/dev/null
stack_before2=$(manifest_stack)

update_out2=$("$SCIAGENT" update --no-pin 2>&1)
update_rc2=$?

if [[ "$update_rc2" -ne 0 ]]; then
    echo "FAIL [$_TEST_NAME] sciagent update --no-pin (stack) exited $update_rc2" >&2
    echo "output: $update_out2" >&2
    exit 1
fi

stack_after2=$(manifest_stack)
if [[ "$stack_before2" != "$stack_after2" ]]; then
    echo "FAIL [$_TEST_NAME] two-role stack changed: before='$stack_before2' after='$stack_after2'" >&2
    exit 1
fi

if ! block_hash_check AGENTS.md; then
    echo "FAIL [$_TEST_NAME] ROLES block hash invalid after update --no-pin (two-role stack)" >&2
    exit 1
fi
if ! block_hash_check AGENTS.md CRAFT; then
    echo "FAIL [$_TEST_NAME] CRAFT block hash invalid after update --no-pin (two-role stack)" >&2
    exit 1
fi

# -------------------------------------------------------------------------
# Test C: no active stack → non-zero exit with guidance message.
# -------------------------------------------------------------------------
cd "$TMPDIR_TEST"
mkdir no-stack-project && cd no-stack-project

no_stack_out=$("$SCIAGENT" update --no-pin 2>&1)
no_stack_rc=$?

if [[ "$no_stack_rc" -eq 0 ]]; then
    echo "FAIL [$_TEST_NAME] expected non-zero exit when no stack is active" >&2
    exit 1
fi

if ! printf '%s\n' "$no_stack_out" | grep -qi "no active stack"; then
    echo "FAIL [$_TEST_NAME] expected guidance message about no active stack" >&2
    echo "output: $no_stack_out" >&2
    exit 1
fi

# -------------------------------------------------------------------------
# Test D: unknown verb argument is rejected.
# -------------------------------------------------------------------------
cd "$TMPDIR_TEST/project"
bad_out=$("$SCIAGENT" update --no-pin --unknown-flag 2>&1)
bad_rc=$?
if [[ "$bad_rc" -eq 0 ]]; then
    echo "FAIL [$_TEST_NAME] expected non-zero exit for unknown flag" >&2
    exit 1
fi

pass
