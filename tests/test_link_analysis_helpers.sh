#!/usr/bin/env bash
# Analysis projects receive the shared helper libraries and owned import shims.

set -u
. "$(dirname "$0")/_lib.sh"

setup_tmpdir
SCIAGENT="$TOOLKIT_ROOT/bin/sciagent"
SHIM_TEMPLATES="$TOOLKIT_ROOT/templates/project/analysis/02_analysis/helpers"
STATE_DIR=".sciagent/helper_shim_state"

mkdir -p analysis/02_analysis
first=$("$SCIAGENT" link --project-dir "$TMPDIR_TEST/analysis" 2>&1) || {
    echo "FAIL [$_TEST_NAME] analysis link failed" >&2
    printf '%s\n' "$first" >&2
    exit 1
}

for libdir in figure-style interactive-style; do
    path="$TMPDIR_TEST/analysis/02_analysis/helpers/$libdir"
    assert_symlink "$path" "$libdir mount missing"
    assert_eq "$(realpath "$path")" "$TOOLKIT_ROOT/lib/$libdir" "$libdir target"
done

n_shims=0
for template in "$SHIM_TEMPLATES"/*.template; do
    [[ -f "$template" ]] || continue
    base=$(basename "${template%.template}")
    shim="$TMPDIR_TEST/analysis/02_analysis/helpers/$base"
    assert_file_eq "$shim" "$template" "$base materialized verbatim"
    assert_file_exists "$TMPDIR_TEST/analysis/$STATE_DIR/$base.sha1" "$base ownership state missing"
    [[ ! -x "$shim" ]] || {
        echo "FAIL [$_TEST_NAME] imported shim is executable: $base" >&2
        exit 1
    }
    n_shims=$((n_shims + 1))
done
assert_eq "$n_shims" 3 "expected three helper shims"

second=$("$SCIAGENT" link --project-dir "$TMPDIR_TEST/analysis" 2>&1) || {
    echo "FAIL [$_TEST_NAME] idempotent analysis link failed" >&2
    printf '%s\n' "$second" >&2
    exit 1
}
assert_eq "$second" "" "second analysis link is a silent no-op"

shim="$TMPDIR_TEST/analysis/02_analysis/helpers/figure_style.py"
printf '\n# project theme extension\n' >> "$shim"
edited_hash=$(sha1sum "$shim" | cut -d' ' -f1)
ceded=$("$SCIAGENT" link --project-dir "$TMPDIR_TEST/analysis" 2>&1) || {
    echo "FAIL [$_TEST_NAME] link failed while ceding a user shim" >&2
    printf '%s\n' "$ceded" >&2
    exit 1
}
assert_eq "$(sha1sum "$shim" | cut -d' ' -f1)" "$edited_hash" "user-edited shim was clobbered"
case "$ceded" in
    *"leaving it as yours"*) : ;;
    *) echo "FAIL [$_TEST_NAME] user-edited shim was not ceded" >&2; printf '%s\n' "$ceded" >&2; exit 1 ;;
esac
assert_file_exists "$TMPDIR_TEST/analysis/$STATE_DIR/figure_style.py.ceded" "cede marker missing"
[[ ! -f "$TMPDIR_TEST/analysis/$STATE_DIR/figure_style.py.sha1" ]] || {
    echo "FAIL [$_TEST_NAME] ceded shim retained toolkit ownership state" >&2
    exit 1
}

mkdir coordination
"$SCIAGENT" link --project-dir "$TMPDIR_TEST/coordination" >/dev/null 2>&1 || {
    echo "FAIL [$_TEST_NAME] non-analysis link failed" >&2
    exit 1
}
[[ ! -e "$TMPDIR_TEST/coordination/02_analysis" ]] || {
    echo "FAIL [$_TEST_NAME] link created an analysis layout in a non-analysis repo" >&2
    exit 1
}
[[ ! -e "$TMPDIR_TEST/coordination/$STATE_DIR" ]] || {
    echo "FAIL [$_TEST_NAME] link wrote helper state in a non-analysis repo" >&2
    exit 1
}

mkdir -p legacy/02_analysis/helpers
ln -s "$TOOLKIT_ROOT/lib/figure-style" legacy/02_analysis/helpers/figure-style
ln -s /workspaces/demo/01_modules/SciAgent-toolkit/lib/interactive-style \
    legacy/02_analysis/helpers/interactive-style
"$SCIAGENT" link --project-dir "$TMPDIR_TEST/legacy" >/dev/null 2>&1 || {
    echo "FAIL [$_TEST_NAME] legacy-link convergence failed" >&2
    exit 1
}
for libdir in figure-style interactive-style; do
    path="$TMPDIR_TEST/legacy/02_analysis/helpers/$libdir"
    assert_symlink "$path" "legacy $libdir mount missing"
    assert_eq "$(realpath "$path")" "$TOOLKIT_ROOT/lib/$libdir" "legacy $libdir target"
done

pass
