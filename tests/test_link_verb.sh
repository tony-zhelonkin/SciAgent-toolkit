#!/usr/bin/env bash
# The link verb creates six whole-tree links, ensures hooks, and converges quietly.

set -u
. "$(dirname "$0")/_lib.sh"

setup_tmpdir
SCIO="$TOOLKIT_ROOT/bin/scio"

mkdir project
first=$("$SCIO" link --project-dir "$TMPDIR_TEST/project" 2>&1) || {
    echo "FAIL [$_TEST_NAME] initial link failed" >&2
    printf '%s\n' "$first" >&2
    exit 1
}

for harness in .claude .agents; do
    for category in skills agents commands; do
        path="$TMPDIR_TEST/project/$harness/$category"
        assert_symlink "$path"
        assert_eq "$(readlink "$path")" "$TOOLKIT_ROOT/$category" "$path target"
    done
done

assert_file_exists "$TMPDIR_TEST/project/.claude/hooks/no_ephemeral.sh"
assert_file_exists "$TMPDIR_TEST/project/.claude/hooks/caption_sweep.sh"
assert_file_exists "$TMPDIR_TEST/project/.claude/settings.json"
assert_file_exists "$TMPDIR_TEST/project/.gitignore"
assert_grep 'BEGIN SCIO:GITIGNORE' "$TMPDIR_TEST/project/.gitignore"

second=$("$SCIO" link --project-dir "$TMPDIR_TEST/project" 2>&1) || {
    echo "FAIL [$_TEST_NAME] idempotent link failed" >&2
    printf '%s\n' "$second" >&2
    exit 1
}
assert_eq "$second" "" "second link is a silent no-op"

mkdir -p replacement/.claude elsewhere
ln -s "$TMPDIR_TEST/elsewhere" replacement/.claude/skills
mkdir -p replacement/.scio
printf '{"version":2}\n' > replacement/.scio/manifest.json
changed=$("$SCIO" link --project-dir "$TMPDIR_TEST/replacement" 2>&1) || {
    echo "FAIL [$_TEST_NAME] replacement link failed" >&2
    printf '%s\n' "$changed" >&2
    exit 1
}
assert_eq "$(readlink replacement/.claude/skills)" "$TOOLKIT_ROOT/skills" "wrong link replaced"
case "$changed" in
    *"replaced link: .claude/skills"*) : ;;
    *) echo "FAIL [$_TEST_NAME] replacement was not reported" >&2; printf '%s\n' "$changed" >&2; exit 1 ;;
esac
[[ ! -e replacement/.scio/manifest.json ]] || {
    echo "FAIL [$_TEST_NAME] legacy manifest survived link" >&2
    exit 1
}

help=$("$SCIO" link --help 2>&1) || {
    echo "FAIL [$_TEST_NAME] link --help failed" >&2
    exit 1
}
case "$help" in *"--project-dir D"*) : ;; *) echo "FAIL [$_TEST_NAME] help omits --project-dir" >&2; exit 1 ;; esac

pass
