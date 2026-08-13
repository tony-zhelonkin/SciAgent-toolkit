#!/usr/bin/env bash
# A real populated category directory is preserved and refused loudly.

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

assert_eq "$rc" "1" "populated directory refusal exits nonzero"
assert_file_exists project/.claude/skills/private.md "private file preserved"
assert_symlink project/.claude/skills/outside-link "outside link preserved"
[[ ! -e project/.claude/skills/old-toolkit-mount ]] || {
    echo "FAIL [$_TEST_NAME] legacy toolkit mount was not swept" >&2
    exit 1
}
[[ -d project/.claude/skills && ! -L project/.claude/skills ]] || {
    echo "FAIL [$_TEST_NAME] populated directory was replaced" >&2
    exit 1
}

for expected in ".claude/skills" "private.md" "outside-link" "Move these entries outside"; do
    case "$out" in
        *"$expected"*) : ;;
        *) echo "FAIL [$_TEST_NAME] refusal omits '$expected'" >&2; printf '%s\n' "$out" >&2; exit 1 ;;
    esac
done

for path in .claude/agents .claude/commands .agents/skills .agents/agents .agents/commands; do
    assert_symlink "project/$path" "unblocked category still converged: $path"
done

pass
