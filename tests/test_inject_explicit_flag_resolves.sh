#!/usr/bin/env bash
# tests/test_inject_explicit_flag_resolves.sh — same ambiguous fixture as
# test_inject_ambiguous_hard_fail.sh, but injecting with an explicit
# --command or --skill flag must succeed and mount only that kind.
set -u
. "$(dirname "$0")/_lib.sh"

setup_tmpdir
FAKE="$TMPDIR_TEST/fake-toolkit"
build_fake_toolkit "$FAKE"
export SCIAGENT_TOOLKIT="$FAKE"
SCIAGENT="$FAKE/bin/sciagent"

# Plant same-named skill + command (ambiguous on bare-name).
mkdir -p "$FAKE/skills/twin"
cat > "$FAKE/skills/twin/SKILL.md" <<'EOF'
---
metadata:
  scope: implementation
  requires: []
  complementary-skills: []
  contraindications: []
  tags: []
---
twin as a skill
EOF
echo "twin as a command" > "$FAKE/commands/twin.md"

mkdir project && cd project
"$SCIAGENT" activate base >/dev/null

# Explicit --command must succeed.
"$SCIAGENT" inject --command twin >/dev/null

assert_symlink .claude/commands/twin.md "command twin mounts under .claude/commands/"
assert_grep '"kind": "command"' .sciagent/manifest.json "manifest records kind=command"
if [[ -L .claude/skills/twin ]]; then
    echo "FAIL [$_TEST_NAME] --command must not mount the skill" >&2
    exit 1
fi

# Explicit --skill must also succeed (independently — different kind, same name).
"$SCIAGENT" inject --skill twin >/dev/null

assert_symlink .claude/skills/twin "skill twin mounts under .claude/skills/"
# Both rows now in manifest.
local_skill_count=$(grep -c '"kind": "skill"' .sciagent/manifest.json || true)
local_cmd_count=$(grep -c '"kind": "command"' .sciagent/manifest.json || true)
if [[ "$local_skill_count" -lt 1 ]]; then
    echo "FAIL [$_TEST_NAME] expected at least one kind=skill row, got $local_skill_count" >&2
    exit 1
fi
if [[ "$local_cmd_count" -lt 1 ]]; then
    echo "FAIL [$_TEST_NAME] expected at least one kind=command row, got $local_cmd_count" >&2
    exit 1
fi

# Two rows name "twin" — one per kind.
twin_rows=$(grep -c '"skill": "twin"' .sciagent/manifest.json || true)
if [[ "$twin_rows" -ne 2 ]]; then
    echo "FAIL [$_TEST_NAME] expected 2 manifest rows for twin (one per kind), got $twin_rows" >&2
    cat .sciagent/manifest.json >&2
    exit 1
fi

pass
