#!/usr/bin/env bash
# Scaffold a project, activate base, run `sciagent roster`; assert the table and
# JSON list the known analysis-base agents and are non-empty.
set -u
. "$(dirname "$0")/_lib.sh"

setup_tmpdir
SCIAGENT="$TOOLKIT_ROOT/bin/sciagent"

# Scaffold an analysis project and activate the real base role inside it.
"$SCIAGENT" new project proj --type analysis >/dev/null
cd proj
"$SCIAGENT" activate base >/dev/null

# Table form: non-empty, has the header and a known agent.
table=$("$SCIAGENT" roster)
if [[ -z "$table" ]]; then
    echo "FAIL [$_TEST_NAME] roster table is empty" >&2
    exit 1
fi
printf '%s\n' "$table" | grep -q '^AGENT' \
    || { echo "FAIL [$_TEST_NAME] roster table missing header" >&2; exit 1; }
printf '%s\n' "$table" | grep -q 'handoff' \
    || { echo "FAIL [$_TEST_NAME] handoff agent missing from roster table" >&2;
         printf '%s\n' "$table" >&2; exit 1; }
printf '%s\n' "$table" | grep -q 'bio-interpreter' \
    || { echo "FAIL [$_TEST_NAME] bio-interpreter missing from roster table" >&2;
         printf '%s\n' "$table" >&2; exit 1; }

# Domain column populated from frontmatter for handoff.
printf '%s\n' "$table" | grep handoff | grep -q 'session-management' \
    || { echo "FAIL [$_TEST_NAME] handoff domain not shown" >&2;
         printf '%s\n' "$table" >&2; exit 1; }

# JSON form: non-empty, valid-looking, lists the agent and the generated date.
json=$("$SCIAGENT" roster --json)
printf '%s\n' "$json" | grep -q '"generated_date"' \
    || { echo "FAIL [$_TEST_NAME] roster --json missing generated_date" >&2;
         printf '%s\n' "$json" >&2; exit 1; }
printf '%s\n' "$json" | grep -q '"name": "handoff"' \
    || { echo "FAIL [$_TEST_NAME] roster --json missing handoff agent" >&2;
         printf '%s\n' "$json" >&2; exit 1; }
printf '%s\n' "$json" | grep -q '"role_stack"' \
    || { echo "FAIL [$_TEST_NAME] roster --json missing role_stack" >&2; exit 1; }

# Graceful failure when no role is active.
cd "$TMPDIR_TEST"
mkdir bare && cd bare
assert_exit 1 "$SCIAGENT" roster

pass
