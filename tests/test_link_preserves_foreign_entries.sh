#!/usr/bin/env bash
# Entries the toolkit does not own survive a bind, whatever shape they take:
# a plain file, and a symlink pointing outside the toolkit. A stale mount that
# names a catalog entry the toolkit no longer carries is still swept.

set -u
. "$(dirname "$0")/_lib.sh"

setup_tmpdir
SCIO="$TOOLKIT_ROOT/bin/scio"

mkdir -p project/.claude/skills project/private-target
printf 'private skill\n' > project/.claude/skills/private.md
ln -s "$TMPDIR_TEST/project/private-target" project/.claude/skills/outside-link
ln -s "$TOOLKIT_ROOT/skills/figure-style" project/.claude/skills/old-toolkit-mount

set +e
out=$("$SCIO" link --project-dir "$TMPDIR_TEST/project" 2>&1)
rc=$?
set -e

assert_eq "$rc" "0" "a project-owned category directory binds rather than refusing"
assert_file_exists project/.claude/skills/private.md "private file preserved"
assert_symlink project/.claude/skills/outside-link "outside link preserved"
[[ ! -e project/.claude/skills/old-toolkit-mount ]] || {
    echo "FAIL [$_TEST_NAME] stale toolkit mount was not swept" >&2
    printf '%s\n' "$out" >&2
    exit 1
}
[[ -d project/.claude/skills && ! -L project/.claude/skills ]] || {
    echo "FAIL [$_TEST_NAME] the directory holding foreign entries was replaced" >&2
    exit 1
}

assert_symlink project/.claude/skills/figure-style "catalog entry bound beside the foreign ones"

for path in .claude/agents .claude/commands .agents/skills .agents/agents .agents/commands; do
    assert_symlink "project/$path" "category with nothing foreign stays one symlink: $path"
done

pass
