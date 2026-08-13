#!/usr/bin/env bash
# Every dispatched verb answers help without mutating the project.

set -u
. "$(dirname "$0")/_lib.sh"

setup_tmpdir
FAKE="$TMPDIR_TEST/fake-toolkit"
build_fake_toolkit "$FAKE"
export SCIO_TOOLKIT="$FAKE"
SCIO="$FAKE/bin/scio"

[[ ! -e "$TOOLKIT_ROOT/bin/sciagent" ]] || {
    echo "FAIL [$_TEST_NAME] legacy command still exists; the rebrand has no compatibility shim" >&2
    exit 1
}

fail() {
    echo "FAIL [$_TEST_NAME] $1" >&2
    shift
    [[ $# -gt 0 ]] && { echo "--- output ---" >&2; printf '%s\n' "$@" >&2; }
    exit 1
}

WORK="$TMPDIR_TEST/work"
mkdir -p "$WORK"

for verb in craft link lint; do
    for flag in -h --help; do
        set +e
        out=$(cd "$WORK" && "$SCIO" "$verb" "$flag" 2>&1)
        rc=$?
        set -e
        [[ "$rc" -eq 0 ]] || fail "scio $verb $flag exited $rc" "$out"
        printf '%s\n' "$out" | grep -q "scio $verb" \
            || fail "scio $verb $flag does not name the verb" "$out"
    done
done

[[ -z "$(ls -A "$WORK")" ]] \
    || fail "help created project files: $(ls -A "$WORK" | tr '\n' ' ')"

for retired in activate deactivate status list new update validate gitignore; do
    set +e
    out=$("$SCIO" "$retired" --help 2>&1)
    rc=$?
    set -e
    [[ "$rc" -ne 0 ]] || fail "retired verb '$retired' is still dispatched" "$out"
    printf '%s\n' "$out" | grep -q "unknown verb '$retired'" \
        || fail "retired verb '$retired' lacks an unknown-verb error" "$out"
done

top=$("$SCIO" --help 2>&1) || fail "top-level help failed"
printf '%s\n' "$top" | grep -q 'project bindings, CRAFT rendering' \
    || fail "top-level help has the old dispatcher description" "$top"

set +e
"$SCIO" help >/dev/null 2>&1; help_rc=$?
"$SCIO" >/dev/null 2>&1; bare_rc=$?
set -e
[[ "$help_rc" -eq 0 ]] || fail "scio help exited $help_rc"
[[ "$bare_rc" -eq 1 ]] || fail "bare scio exited $bare_rc"

pass
