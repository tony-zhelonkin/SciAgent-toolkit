#!/usr/bin/env bash
# Legacy toolkit links are swept from a scratch fixture; user paths survive.

set -u
. "$(dirname "$0")/_lib.sh"

setup_tmpdir
SCIAGENT="$TOOLKIT_ROOT/bin/sciagent"

assert_absent() {
    [[ ! -e "$1" && ! -L "$1" ]] || {
        echo "FAIL [$_TEST_NAME] legacy link survived: $1" >&2
        exit 1
    }
}

mkdir -p project/.claude/output-styles project/02_analysis/helpers outside
for harness in .claude .agents; do
    for category in skills agents commands; do
        mkdir -p "project/$harness/$category"
        case "$category" in
            skills) source="$TOOLKIT_ROOT/skills/figure-style" ;;
            agents) source="$TOOLKIT_ROOT/agents/architect.md" ;;
            commands) source="$TOOLKIT_ROOT/commands/commit.md" ;;
        esac
        ln -s "$source" "project/$harness/$category/current-owned"
    done
done
ln -s /workspaces/demo/01_modules/SciAgent-toolkit/skills/removed project/.claude/skills/container-dangling

ln -s "$TOOLKIT_ROOT/output-styles/retired" project/.claude/output-styles/current-retired
ln -s /workspaces/demo/01_modules/SciAgent-toolkit/output-styles/retired project/.claude/output-styles/container-retired
ln -s "$TMPDIR_TEST/outside" project/.claude/output-styles/outside-link
printf 'user style\n' > project/.claude/output-styles/user-style.md

ln -s "$TOOLKIT_ROOT/lib/figure-style" project/02_analysis/helpers/figure-style
ln -s /workspaces/demo/01_modules/SciAgent-toolkit/lib/interactive-style project/02_analysis/helpers/interactive-style
ln -s "$TOOLKIT_ROOT/lib/figure-style" project/02_analysis/helpers/retired-helper
printf 'user shim\n' > project/02_analysis/helpers/user_helper.py

out=$("$SCIAGENT" link --project-dir "$TMPDIR_TEST/project" 2>&1) || {
    echo "FAIL [$_TEST_NAME] link failed" >&2
    printf '%s\n' "$out" >&2
    exit 1
}

for harness in .claude .agents; do
    for category in skills agents commands; do
        assert_symlink "project/$harness/$category" "legacy mount directory collapsed"
        assert_eq "$(readlink "project/$harness/$category")" "$TOOLKIT_ROOT/$category"
    done
done

assert_absent project/.claude/output-styles/current-retired
assert_absent project/.claude/output-styles/container-retired
assert_symlink project/.claude/output-styles/outside-link "outside output-style link preserved"
assert_file_exists project/.claude/output-styles/user-style.md

assert_symlink project/02_analysis/helpers/figure-style
assert_symlink project/02_analysis/helpers/interactive-style
assert_eq "$(realpath project/02_analysis/helpers/figure-style)" "$TOOLKIT_ROOT/lib/figure-style"
assert_eq "$(realpath project/02_analysis/helpers/interactive-style)" "$TOOLKIT_ROOT/lib/interactive-style"
assert_eq "$(readlink project/02_analysis/helpers/figure-style)" "$TOOLKIT_ROOT/lib/figure-style" \
    "live helper mount was rewritten"
assert_absent project/02_analysis/helpers/retired-helper
assert_file_exists project/02_analysis/helpers/user_helper.py

case "$out" in
    *"container-dangling"*"container-retired"*"interactive-style"*) : ;;
    *) echo "FAIL [$_TEST_NAME] sweep did not report every legacy link" >&2; printf '%s\n' "$out" >&2; exit 1 ;;
esac

pass
