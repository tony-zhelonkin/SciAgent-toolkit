#!/usr/bin/env bash
# tests/test_eject_rejects_stack_mounted_command.sh — `sciagent eject --command
# <name>` (and bare-name eject) must refuse when the command is listed under
# `commands:` in any role YAML on the stack, with the same deactivate pointer
# as the existing skill-side guard. Symmetric extension of F-1.
set -u
. "$(dirname "$0")/_lib.sh"

setup_tmpdir
FAKE="$TMPDIR_TEST/fake-toolkit"
build_fake_toolkit "$FAKE"
export SCIAGENT_TOOLKIT="$FAKE"
SCIAGENT="$FAKE/bin/sciagent"

# Mirror the F-1 regression pattern: real role YAMLs carry trailing comments
# on `  - <name>    # description` lines. Synthesise that shape for the
# command section to confirm role_load strips the comment cleanly.
cat > "$FAKE/roles/base.yaml" <<'EOF'
name: base
description: fixture base role with commented command lines
skills:
  - s_a
  - s_b
agents:
  - ag_a
commands:
  - c_a   # short description
EOF

mkdir project && cd project
"$SCIAGENT" activate base >/dev/null

# c_a is in roles/base.yaml under `commands:` — eject --command must refuse.
out=$("$SCIAGENT" eject --command c_a 2>&1)
rc=$?
if [[ "$rc" -eq 0 ]]; then
    echo "FAIL [$_TEST_NAME] expected exit 1 ejecting stack-mounted command c_a, got 0" >&2
    echo "output: $out" >&2
    exit 1
fi
printf '%s\n' "$out" | grep -qi 'deactivate' || {
    echo "FAIL [$_TEST_NAME] error message should mention 'sciagent deactivate', got: $out" >&2
    exit 1
}

# Bare-name eject (no kind flag) also refuses with the same pointer.
out2=$("$SCIAGENT" eject c_a 2>&1)
rc2=$?
if [[ "$rc2" -eq 0 ]]; then
    echo "FAIL [$_TEST_NAME] bare-name eject of stack-mounted command should fail" >&2
    exit 1
fi
printf '%s\n' "$out2" | grep -qi 'deactivate' || {
    echo "FAIL [$_TEST_NAME] bare-name eject error should mention 'sciagent deactivate', got: $out2" >&2
    exit 1
}

# Stack symlink for c_a survives both refused ejects.
assert_symlink .claude/commands/c_a.md "command symlink damaged by refused eject"

pass
