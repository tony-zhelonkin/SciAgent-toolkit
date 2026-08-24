#!/usr/bin/env bash
# Both live toolkit directory names participate in discovery and lint checks.

set -u
. "$(dirname "$0")/_lib.sh"

setup_tmpdir
SCIO="$TOOLKIT_ROOT/bin/scio"

# _build_vendored_project <project> <toolkit-dir-name>
_build_vendored_project() {
    local project="$1" toolkit_dir="$2"
    mkdir -p "$project/01_modules"
    build_fake_toolkit "$project/01_modules/$toolkit_dir"
}

# _link_and_check_mounts <project> <toolkit-dir-name>
_link_and_check_mounts() {
    local project="$1" toolkit_dir="$2" out harness category mount
    out=$(SCIO_TOOLKIT="$project/01_modules/$toolkit_dir" \
        "$SCIO" link --project-dir "$project" 2>&1) || {
        echo "FAIL [$_TEST_NAME] link failed for $toolkit_dir" >&2
        printf '%s\n' "$out" >&2
        exit 1
    }
    for harness in .claude .agents; do
        for category in skills agents commands; do
            mount="$project/$harness/$category"
            assert_symlink "$mount"
            assert_eq "$(realpath "$mount")" \
                "$(realpath "$project/01_modules/$toolkit_dir/$category")" \
                "$mount resolves through $toolkit_dir"
        done
    done
}

# Both directory names remain valid project-local toolkit locations.
NEW_PROJECT="$TMPDIR_TEST/new-project"
_build_vendored_project "$NEW_PROJECT" scio
_link_and_check_mounts "$NEW_PROJECT" scio

LEGACY_PROJECT="$TMPDIR_TEST/legacy-project"
_build_vendored_project "$LEGACY_PROJECT" SciAgent-toolkit
_link_and_check_mounts "$LEGACY_PROJECT" SciAgent-toolkit

# The new name wins when both candidates are present.
BOTH_PROJECT="$TMPDIR_TEST/both-project"
_build_vendored_project "$BOTH_PROJECT" scio
build_fake_toolkit "$BOTH_PROJECT/01_modules/SciAgent-toolkit"
_link_and_check_mounts "$BOTH_PROJECT" scio

# _assert_vendor_path_finding <toolkit-dir-name> <case-name>
_assert_vendor_path_finding() {
    local toolkit_dir="$1" case_name="$2" project out rc
    project="$TMPDIR_TEST/$case_name"
    mkdir -p "$project"
    git -C "$project" init -q
    printf 'Read 01_modules/%s/skills/example/SKILL.md.\n' "$toolkit_dir" \
        > "$project/AGENTS.md"
    git -C "$project" add AGENTS.md

    set +e
    out=$(SCIO_TOOLKIT="$TOOLKIT_ROOT" "$SCIO" lint --check harness-links \
        --project-dir "$project" 2>&1)
    rc=$?
    set -e
    if [[ "$rc" -ne 0 ]] \
       || ! printf '%s\n' "$out" | grep -q 'AGENTS.md reaches a skill by vendor path'; then
        echo "FAIL [$_TEST_NAME] harness-links missed $toolkit_dir (rc=$rc)" >&2
        printf '%s\n' "$out" >&2
        exit 1
    fi
}

_assert_vendor_path_finding scio lint-new-name
_assert_vendor_path_finding SciAgent-toolkit lint-legacy-name

# A declared submodule may place the toolkit below any parent directory.
SUBMODULE_PROJECT="$TMPDIR_TEST/submodule-project"
build_fake_toolkit "$SUBMODULE_PROJECT/vendor/scio"
printf '%s\n' \
    '[submodule "toolkit"]' \
    '    path = vendor/scio' \
    '    url = ../toolkit.git' \
    > "$SUBMODULE_PROJECT/.gitmodules"
out=$(SCIO_TOOLKIT="$SUBMODULE_PROJECT/vendor/scio" \
    "$SCIO" link --project-dir "$SUBMODULE_PROJECT" 2>&1) || {
    echo "FAIL [$_TEST_NAME] .gitmodules discovery failed for scio" >&2
    printf '%s\n' "$out" >&2
    exit 1
}
assert_symlink "$SUBMODULE_PROJECT/.claude/skills"
assert_eq "$(realpath "$SUBMODULE_PROJECT/.claude/skills")" \
    "$(realpath "$SUBMODULE_PROJECT/vendor/scio/skills")" \
    ".gitmodules mount resolves through scio"

pass
