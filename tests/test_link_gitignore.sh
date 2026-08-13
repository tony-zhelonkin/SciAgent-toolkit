#!/usr/bin/env bash
# Link refreshes the legacy gitignore block and converges silently.
set -u
. "$(dirname "$0")/_lib.sh"

setup_tmpdir
SCIAGENT="$TOOLKIT_ROOT/bin/sciagent"
mkdir project
cat > project/.gitignore <<'EOF'
user-entry/

# BEGIN SCIAGENT:GITIGNORE
stale-entry/
# END SCIAGENT:GITIGNORE
EOF
chmod 640 project/.gitignore

first=$("$SCIAGENT" link --project-dir "$TMPDIR_TEST/project" 2>&1) || {
    echo "FAIL [$_TEST_NAME] initial link failed" >&2
    printf '%s\n' "$first" >&2
    exit 1
}

assert_grep '^user-entry/$' project/.gitignore "user entry preserved"
assert_grep '^02_analysis/helpers/figure-style$' project/.gitignore \
    "analysis helper ignore missing"
assert_eq "$(stat -c '%a' project/.gitignore)" "640" "gitignore mode preserved"
[[ $(grep -c '^# BEGIN SCIAGENT:GITIGNORE$' project/.gitignore) -eq 1 ]] || {
    echo "FAIL [$_TEST_NAME] expected one managed block" >&2
    exit 1
}

before=$(cat project/.gitignore)
second=$("$SCIAGENT" link --project-dir "$TMPDIR_TEST/project" 2>&1) || {
    echo "FAIL [$_TEST_NAME] idempotent link failed" >&2
    printf '%s\n' "$second" >&2
    exit 1
}
assert_eq "$second" "" "second link is silent"
assert_eq "$(cat project/.gitignore)" "$before" "gitignore bytes converge"

pass
