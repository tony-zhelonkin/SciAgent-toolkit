#!/usr/bin/env bash
# Link refreshes the legacy gitignore block and converges silently.
set -u
. "$(dirname "$0")/_lib.sh"

setup_tmpdir
SCIO="$TOOLKIT_ROOT/bin/scio"
mkdir project
cat > project/.gitignore <<'EOF'
user-entry/

# BEGIN SCIAGENT:GITIGNORE
stale-entry/
# END SCIAGENT:GITIGNORE
EOF
chmod 640 project/.gitignore

first=$("$SCIO" link --project-dir "$TMPDIR_TEST/project" 2>&1) || {
    echo "FAIL [$_TEST_NAME] initial link failed" >&2
    printf '%s\n' "$first" >&2
    exit 1
}

assert_grep '^user-entry/$' project/.gitignore "user entry preserved"
assert_grep '^02_analysis/helpers/figure-style$' project/.gitignore \
    "analysis helper ignore missing"
# Every helper mount link creates must be ignored. interactive-style was absent
# from the block while link created the symlink, so a generated mount showed up
# as untracked project content in each consumer.
assert_grep '^02_analysis/helpers/interactive-style$' project/.gitignore \
    "interactive-style mount is not ignored"
for _mount in figure-style interactive-style; do
    [[ -L "project/02_analysis/helpers/$_mount" ]] || continue
    git -C project check-ignore -q "02_analysis/helpers/$_mount" 2>/dev/null || {
        [[ -d project/.git ]] || break
        echo "FAIL [$_TEST_NAME] link created $_mount but git does not ignore it" >&2
        exit 1
    }
done
assert_eq "$(stat -c '%a' project/.gitignore)" "640" "gitignore mode preserved"
[[ $(grep -c '^# BEGIN SCIO:GITIGNORE$' project/.gitignore) -eq 1 ]] || {
    echo "FAIL [$_TEST_NAME] expected one managed block" >&2
    exit 1
}
if grep -q 'SCIAGENT:GITIGNORE' project/.gitignore; then
    echo "FAIL [$_TEST_NAME] legacy gitignore marker survived migration" >&2
    exit 1
fi

before=$(cat project/.gitignore)
second=$("$SCIO" link --project-dir "$TMPDIR_TEST/project" 2>&1) || {
    echo "FAIL [$_TEST_NAME] idempotent link failed" >&2
    printf '%s\n' "$second" >&2
    exit 1
}
assert_eq "$second" "" "second link is silent"
assert_eq "$(cat project/.gitignore)" "$before" "gitignore bytes converge"

pass
