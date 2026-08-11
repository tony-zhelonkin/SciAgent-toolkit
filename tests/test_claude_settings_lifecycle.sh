#!/usr/bin/env bash
# claude_settings: activate must write outputStyle into
# .claude/settings.local.json; deactivate must revert based on the saved
# state tag (delete the file if we created it; strip the key if we added
# it to a pre-existing file).
#
# Uses the real architect role with an explicit `--output-style
# architect-mentor` flag (Phase 5d: output_style is no longer role-scoped)
# rather than the fake-toolkit fixture, because the fixture roles don't
# ship a matching system-prompts/ fixture.
set -u
. "$(dirname "$0")/_lib.sh"

setup_tmpdir
SCIAGENT="$TOOLKIT_ROOT/bin/sciagent"
export SCIAGENT_TOOLKIT="$TOOLKIT_ROOT"

# ----- Case 1: file does not exist -----------------------------------------
mkdir case1 && cd case1
echo "# stub" > AGENTS.md

"$SCIAGENT" activate architect --output-style architect-mentor >/dev/null

assert_file_exists .claude/settings.local.json "case1: file created"
assert_file_exists .sciagent/claude_settings.state "case1: state file written"
state=$(cat .sciagent/claude_settings.state)
assert_eq "$state" "created" "case1: state tag is 'created'"
got=$(grep -o '"outputStyle"[[:space:]]*:[[:space:]]*"architect-mentor"' .claude/settings.local.json || true)
[[ -n "$got" ]] || { echo "FAIL [$_TEST_NAME] case1: outputStyle key not set"; cat .claude/settings.local.json >&2; exit 1; }

"$SCIAGENT" deactivate >/dev/null

[[ -e .claude/settings.local.json ]] && { echo "FAIL [$_TEST_NAME] case1: file not removed on deactivate"; exit 1; }
[[ -e .sciagent/claude_settings.state ]] && { echo "FAIL [$_TEST_NAME] case1: state file not cleaned"; exit 1; }

cd ..

# ----- Case 2: file exists, no outputStyle key -----------------------------
mkdir case2 && cd case2
echo "# stub" > AGENTS.md
mkdir -p .claude
cat > .claude/settings.local.json <<'EOF'
{
  "permissions": {
    "allow": ["Bash(git:*)"],
    "deny": []
  }
}
EOF
cp .claude/settings.local.json /tmp/.case2-orig.json

"$SCIAGENT" activate architect --output-style architect-mentor >/dev/null

state=$(cat .sciagent/claude_settings.state)
assert_eq "$state" "existed-no-style" "case2: state tag is 'existed-no-style'"
got=$(grep -o '"outputStyle"[[:space:]]*:[[:space:]]*"architect-mentor"' .claude/settings.local.json || true)
[[ -n "$got" ]] || { echo "FAIL [$_TEST_NAME] case2: outputStyle not added"; cat .claude/settings.local.json >&2; exit 1; }
# Pre-existing key must still be there.
grep -q '"permissions"' .claude/settings.local.json || { echo "FAIL [$_TEST_NAME] case2: pre-existing 'permissions' key lost"; cat .claude/settings.local.json >&2; exit 1; }

"$SCIAGENT" deactivate >/dev/null

assert_file_exists .claude/settings.local.json "case2: file preserved on deactivate"
[[ -e .sciagent/claude_settings.state ]] && { echo "FAIL [$_TEST_NAME] case2: state file not cleaned"; exit 1; }
# outputStyle stripped, permissions preserved.
if grep -q '"outputStyle"' .claude/settings.local.json; then
    echo "FAIL [$_TEST_NAME] case2: outputStyle not stripped on deactivate" >&2
    cat .claude/settings.local.json >&2
    exit 1
fi
grep -q '"permissions"' .claude/settings.local.json || { echo "FAIL [$_TEST_NAME] case2: 'permissions' key lost on deactivate"; cat .claude/settings.local.json >&2; exit 1; }
grep -q 'Bash(git:\*)' .claude/settings.local.json || { echo "FAIL [$_TEST_NAME] case2: 'permissions.allow' content lost on deactivate"; cat .claude/settings.local.json >&2; exit 1; }

cd ..

# ----- Case 3: file exists WITH outputStyle (toolkit overwrites + warns) ---
mkdir case3 && cd case3
echo "# stub" > AGENTS.md
mkdir -p .claude
cat > .claude/settings.local.json <<'EOF'
{
  "outputStyle": "my-custom-style"
}
EOF

out=$("$SCIAGENT" activate architect --output-style architect-mentor 2>&1)

state=$(cat .sciagent/claude_settings.state)
assert_eq "$state" "existed-with-style" "case3: state tag is 'existed-with-style'"
# Warning must mention the prior value so the user can record it.
printf '%s\n' "$out" | grep -q "my-custom-style" || { echo "FAIL [$_TEST_NAME] case3: overwrite warning did not mention prior value"; printf '%s\n' "$out" >&2; exit 1; }
got=$(grep -o '"outputStyle"[[:space:]]*:[[:space:]]*"architect-mentor"' .claude/settings.local.json || true)
[[ -n "$got" ]] || { echo "FAIL [$_TEST_NAME] case3: outputStyle not overwritten"; cat .claude/settings.local.json >&2; exit 1; }

"$SCIAGENT" deactivate >/dev/null

# v1 behavior: prior value NOT restored (documented limitation).
if [[ -f .claude/settings.local.json ]] && grep -q '"outputStyle"' .claude/settings.local.json; then
    echo "FAIL [$_TEST_NAME] case3: outputStyle still present after deactivate" >&2
    cat .claude/settings.local.json >&2
    exit 1
fi

pass
