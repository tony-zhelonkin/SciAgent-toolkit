#!/usr/bin/env bash
# tests/test_inject_command_with_companion_skill.sh — when an explicit
# --command inject targets a name that also exists as a skill, a one-line
# stderr note advertises the companion skill, but the skill is NOT mounted.
set -u
. "$(dirname "$0")/_lib.sh"

setup_tmpdir
FAKE="$TMPDIR_TEST/fake-toolkit"
build_fake_toolkit "$FAKE"
export SCIAGENT_TOOLKIT="$FAKE"
SCIAGENT="$FAKE/bin/sciagent"

# Plant same-named skill + command (the companion-skill fixture).
mkdir -p "$FAKE/skills/duo"
cat > "$FAKE/skills/duo/SKILL.md" <<'EOF'
---
metadata:
  scope: implementation
  requires: []
  complementary-skills: []
  contraindications: []
  tags: []
---
duo as a skill
EOF
echo "duo as a command" > "$FAKE/commands/duo.md"

mkdir project && cd project
"$SCIAGENT" activate base >/dev/null

# Inject the command explicitly. Capture stderr.
out=$("$SCIAGENT" inject --command duo 2>&1 1>/dev/null)
rc=$?
if [[ "$rc" -ne 0 ]]; then
    echo "FAIL [$_TEST_NAME] --command inject failed (rc=$rc): $out" >&2
    exit 1
fi

# The companion-skill note must appear on stderr.
printf '%s\n' "$out" | grep -q "companion skill 'duo'" || {
    echo "FAIL [$_TEST_NAME] expected companion-skill note on stderr, got: $out" >&2
    exit 1
}
printf '%s\n' "$out" | grep -q 'inject --skill duo' || {
    echo "FAIL [$_TEST_NAME] companion note must advertise the explicit --skill form, got: $out" >&2
    exit 1
}

# Command was mounted.
assert_symlink .claude/commands/duo.md "command duo mounts"
# Skill was NOT mounted.
if [[ -L .claude/skills/duo ]]; then
    echo "FAIL [$_TEST_NAME] companion skill must NOT be auto-mounted" >&2
    exit 1
fi

# Manifest has the command row, no skill row for duo.
duo_skill_rows=$(grep '"skill": "duo"' .sciagent/manifest.json | grep -c '"kind": "skill"' || true)
if [[ "$duo_skill_rows" -ne 0 ]]; then
    echo "FAIL [$_TEST_NAME] manifest should have no kind=skill row for duo, got $duo_skill_rows" >&2
    cat .sciagent/manifest.json >&2
    exit 1
fi
duo_cmd_rows=$(grep '"skill": "duo"' .sciagent/manifest.json | grep -c '"kind": "command"' || true)
if [[ "$duo_cmd_rows" -ne 1 ]]; then
    echo "FAIL [$_TEST_NAME] manifest should have exactly 1 kind=command row for duo, got $duo_cmd_rows" >&2
    exit 1
fi

pass
