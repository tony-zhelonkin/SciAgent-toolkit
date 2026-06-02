#!/usr/bin/env bash
# tests/test_inject_collision_warns.sh — inline cross-namespace collision check
# at inject time (5.2 / ADR-5.1 C). Injecting a name that also exists in another
# namespace must:
#   - warn on stderr by default (but still mount: exit 0)
#   - stay silent under --force
#   - hard-fail under SCIAGENT_STRICT_COLLISIONS=1
set -u
. "$(dirname "$0")/_lib.sh"

setup_tmpdir
FAKE="$TMPDIR_TEST/fake-toolkit"
build_fake_toolkit "$FAKE"

# Make "clash" exist as both a skill and a command → cross-namespace collision.
mkdir -p "$FAKE/skills/clash"
echo "skill clash" > "$FAKE/skills/clash/SKILL.md"
echo "cmd clash"   > "$FAKE/commands/clash.md"

export SCIAGENT_TOOLKIT="$FAKE"
SCIAGENT="$FAKE/bin/sciagent"

mkdir project && cd project
"$SCIAGENT" activate base >/dev/null 2>&1

# Default: warn on stderr, but inject succeeds (exit 0) and mounts the skill.
err=$("$SCIAGENT" inject --skill clash 2>&1 1>/dev/null)
rc=$?
[[ "$rc" -eq 0 ]] || {
    echo "FAIL [$_TEST_NAME] default collision inject should exit 0, got $rc" >&2
    exit 1
}
printf '%s\n' "$err" | grep -qi 'collides across' || {
    echo "FAIL [$_TEST_NAME] expected collision warning, got: $err" >&2
    exit 1
}
assert_symlink .claude/skills/clash "collision inject still mounts the skill"

# --force: silent, still mounts. Eject first to re-test from a clean slate.
"$SCIAGENT" eject --skill clash >/dev/null 2>&1
err=$("$SCIAGENT" inject --skill clash --force 2>&1 1>/dev/null)
if printf '%s\n' "$err" | grep -qi 'collides across'; then
    echo "FAIL [$_TEST_NAME] --force must silence the collision warning, got: $err" >&2
    exit 1
fi
assert_symlink .claude/skills/clash "--force inject still mounts the skill"

# Strict mode: hard-fail. Eject first.
"$SCIAGENT" eject --skill clash >/dev/null 2>&1
SCIAGENT_STRICT_COLLISIONS=1 "$SCIAGENT" inject --skill clash >/dev/null 2>&1
rc=$?
[[ "$rc" -ne 0 ]] || {
    echo "FAIL [$_TEST_NAME] SCIAGENT_STRICT_COLLISIONS=1 must hard-fail the inject" >&2
    exit 1
}
[[ ! -L .claude/skills/clash ]] || {
    echo "FAIL [$_TEST_NAME] strict-mode refusal must not mount the skill" >&2
    exit 1
}

pass
