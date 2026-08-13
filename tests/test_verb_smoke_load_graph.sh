#!/usr/bin/env bash
# tests/test_verb_smoke_load_graph.sh — guards the declarative load graph (5.4):
# every dispatched verb must source its full transitive closure, so no verb
# fails with a "command not found" for a helper it relies on. Runs each verb in
# a scratch project against the fixture toolkit and asserts the dispatcher never
# emits a missing-command error.
set -u
. "$(dirname "$0")/_lib.sh"

setup_tmpdir
FAKE="$TMPDIR_TEST/fake-toolkit"
build_fake_toolkit "$FAKE"
mkdir -p "$FAKE/templates"
cp -R "$TOOLKIT_ROOT/templates/skill" "$FAKE/templates/skill"
export SCIAGENT_TOOLKIT="$FAKE"
SCIAGENT="$FAKE/bin/sciagent"

mkdir project && cd project

# Establish a stack first so inject/eject/status/deactivate have something to
# operate on. activate exercises the heaviest closure.
"$SCIAGENT" activate base >/dev/null 2>&1 || {
    echo "FAIL [$_TEST_NAME] activate base failed unexpectedly" >&2
    exit 1
}

# Run each verb and capture combined output; assert no "command not found".
run_verb() {
    local out
    out=$("$SCIAGENT" "$@" 2>&1 || true)
    if printf '%s\n' "$out" | grep -qi 'command not found'; then
        echo "FAIL [$_TEST_NAME] verb '$*' hit a missing command (under-loaded arm):" >&2
        printf '%s\n' "$out" >&2
        exit 1
    fi
}

run_verb status
run_verb status --json
run_verb list
run_verb validate --quiet
run_verb inject s_c
run_verb eject s_c
run_verb craft
run_verb deactivate
run_verb new role smoke-role
"$SCIAGENT" new skill smoke-skill >/dev/null 2>&1 || {
    echo "FAIL [$_TEST_NAME] new skill failed with templates/skill" >&2
    exit 1
}
assert_file_exists "$FAKE/skills/smoke-skill/SKILL.md"
assert_file_eq "$FAKE/skills/smoke-skill/SKILL.md" \
    "$FAKE/templates/skill/SKILL.md" "new skill copied templates/skill"

pass
