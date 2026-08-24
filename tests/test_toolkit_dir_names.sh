#!/usr/bin/env bash
# tests/test_toolkit_dir_names.sh — discovery accepts both live toolkit
# directory names, and evidence strength outranks the name.
#
# ADR-D9 renames the consumer vendor directory from SciAgent-toolkit to scio in
# one fleet pass, so both names are live meanwhile. Six places once hardcoded the
# old one; a project that moved would have lost link and freshness silently.
#
# These assert _link_in_repo_toolkit DIRECTLY. An earlier version of this file
# ran `scio link` with SCIO_TOOLKIT pointing at the expected copy, which proved
# nothing: link treats "no project toolkit discovered" as success, so four of six
# cases passed with discovery replaced by `return 1`.
#
# Tests:
#   1. 01_modules/scio is discovered
#   2. 01_modules/SciAgent-toolkit is discovered
#   3. both present -> scio wins, deterministically
#   4. a declared .gitmodules path outranks a same-named directory elsewhere
#   5. a DECLARED legacy submodule beats a stray new-name directory
#   6. a directory that is not a checkout does not shadow a real one
#   7. no toolkit anywhere -> nonzero, no output
#   8. harness-links reports the vendor path under either name
#   9. freshness sees a behind submodule under either name
set -u
. "$(dirname "$0")/_lib.sh"

setup_tmpdir
SCIO="$TOOLKIT_ROOT/bin/scio"

# _discover <project> — print what link.sh would resolve, from inside <project>.
# Sources the same modules bin/scio does for the link verb, in that order.
_discover() {
    ( cd "$1" && SCIO_TOOLKIT="$TOOLKIT_ROOT" bash -c '
        . "$1/lib/scio/common.sh"
        . "$1/lib/scio/link.sh"
        _link_in_repo_toolkit
      ' _ "$TOOLKIT_ROOT" 2>/dev/null )
}

# _expect <project> <expected-path> <label>
_expect() {
    local got
    got=$(_discover "$1")
    assert_eq "$got" "$2" "$3"
}

# --- 1+2. either name, in the conventional location -------------------------
P=$TMPDIR_TEST/p_new
build_fake_toolkit "$P/01_modules/scio"
_expect "$P" "01_modules/scio" "the new name is discovered"

P=$TMPDIR_TEST/p_legacy
build_fake_toolkit "$P/01_modules/SciAgent-toolkit"
_expect "$P" "01_modules/SciAgent-toolkit" "the legacy name is discovered"

# --- 3. both present: the array order decides, not the filesystem's ---------
P=$TMPDIR_TEST/p_both
build_fake_toolkit "$P/01_modules/SciAgent-toolkit"
build_fake_toolkit "$P/01_modules/scio"
_expect "$P" "01_modules/scio" "the new name wins when both are present"

# --- 4. a declared path is found wherever it sits --------------------------
P=$TMPDIR_TEST/p_declared
build_fake_toolkit "$P/vendor/scio"
printf '%s\n' '[submodule "toolkit"]' '    path = vendor/scio' \
    '    url = ../toolkit.git' > "$P/.gitmodules"
_expect "$P" "vendor/scio" "a declared path is discovered"

# --- 5. evidence outranks the name ----------------------------------------
# The project DECLARED a legacy submodule. An unrelated directory named scio
# must not win: the locality guard would then refuse the project's own pinned
# copy as external, which is a broken consumer rather than a preference.
P=$TMPDIR_TEST/p_declared_legacy
build_fake_toolkit "$P/vendor/SciAgent-toolkit"
build_fake_toolkit "$P/misc/scio"
printf '%s\n' '[submodule "toolkit"]' '    path = vendor/SciAgent-toolkit' \
    '    url = ../toolkit.git' > "$P/.gitmodules"
_expect "$P" "vendor/SciAgent-toolkit" \
    "a declared legacy submodule outranks a stray new-name directory"

# --- 6. a same-named directory is not a checkout --------------------------
P=$TMPDIR_TEST/p_shadow
mkdir -p "$P/01_modules/scio/docs"          # no bin/scio, no craft.yaml
build_fake_toolkit "$P/01_modules/SciAgent-toolkit"
_expect "$P" "01_modules/SciAgent-toolkit" \
    "an empty same-named directory does not shadow a real checkout"

# --- 7. nothing to find --------------------------------------------------
P=$TMPDIR_TEST/p_none
mkdir -p "$P/01_modules"
set +e
out=$(_discover "$P"); rc=$?
set -e
if [[ "$rc" -eq 0 ]] || [[ -n "$out" ]]; then
    echo "FAIL [$_TEST_NAME] discovery invented a toolkit (rc=$rc, out='$out')" >&2
    exit 1
fi

# --- 8. the vendor-path finding covers both names ------------------------
_assert_vendor_finding() {
    local toolkit_dir="$1" project out rc
    project="$TMPDIR_TEST/lint_$2"
    mkdir -p "$project"
    git -C "$project" init -q
    printf 'Read 01_modules/%s/skills/example/SKILL.md.\n' "$toolkit_dir" > "$project/AGENTS.md"
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
_assert_vendor_finding scio new
_assert_vendor_finding SciAgent-toolkit legacy

# --- 9. freshness reads a behind submodule under either name -------------
# A real clone parked one commit back, so the ancestor test the check performs
# has something true to find. Without this case, narrowing the freshness search
# to one name goes unnoticed.
_assert_freshness_behind() {
    local toolkit_dir="$1" project out
    project="$TMPDIR_TEST/fresh_$2"
    mkdir -p "$project/01_modules"
    git clone -q --local --no-hardlinks "$TOOLKIT_ROOT" \
        "$project/01_modules/$toolkit_dir" 2>/dev/null || {
        echo "SKIP [$_TEST_NAME] case9 ($toolkit_dir): clone unavailable" >&2
        return 0
    }
    git -C "$project/01_modules/$toolkit_dir" checkout -q HEAD~1 2>/dev/null || {
        echo "SKIP [$_TEST_NAME] case9 ($toolkit_dir): no parent commit" >&2
        return 0
    }
    set +e
    out=$(SCIO_TOOLKIT="$TOOLKIT_ROOT" "$SCIO" lint --check freshness \
        --project-dir "$project" 2>&1)
    set -e
    printf '%s\n' "$out" | grep -q "01_modules/$toolkit_dir is behind toolkit HEAD" || {
        echo "FAIL [$_TEST_NAME] case9: freshness missed a behind $toolkit_dir" >&2
        printf '%s\n' "$out" >&2
        exit 1
    }
}
_assert_freshness_behind scio new
_assert_freshness_behind SciAgent-toolkit legacy

pass
