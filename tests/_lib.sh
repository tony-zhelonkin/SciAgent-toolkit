# tests/_lib.sh — shared assert helpers and tmpdir scaffold.
# Each test sources this, calls setup_tmpdir, runs assertions.
# On failure: print actual + expected + context, exit non-zero.

# shellcheck shell=bash

set -u

TOOLKIT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
export SCIAGENT_TOOLKIT="$TOOLKIT_ROOT"

_TEST_NAME="${0##*/}"

assert_eq() {
    local actual="$1" expected="$2" msg="${3:-}"
    if [[ "$actual" != "$expected" ]]; then
        echo "FAIL [$_TEST_NAME] ${msg:-assert_eq}" >&2
        echo "  expected: $(printf '%q' "$expected")" >&2
        echo "  actual:   $(printf '%q' "$actual")" >&2
        exit 1
    fi
}

assert_file_eq() {
    local a="$1" b="$2" msg="${3:-}"
    if ! diff -u "$a" "$b" >/dev/null; then
        echo "FAIL [$_TEST_NAME] ${msg:-assert_file_eq}" >&2
        echo "--- expected ($b)" >&2
        echo "+++ actual   ($a)" >&2
        diff -u "$b" "$a" >&2 || true
        exit 1
    fi
}

assert_exit() {
    local expected="$1"; shift
    "$@"
    local rc=$?
    if [[ "$rc" -ne "$expected" ]]; then
        echo "FAIL [$_TEST_NAME] expected exit $expected got $rc: $*" >&2
        exit 1
    fi
}

assert_file_exists() {
    local f="$1" msg="${2:-}"
    if [[ ! -e "$f" ]]; then
        echo "FAIL [$_TEST_NAME] ${msg:-file missing: $f}" >&2
        exit 1
    fi
}

assert_symlink() {
    local f="$1" msg="${2:-}"
    if [[ ! -L "$f" ]]; then
        echo "FAIL [$_TEST_NAME] ${msg:-not a symlink: $f}" >&2
        ls -la "$(dirname "$f")" >&2 || true
        exit 1
    fi
}

assert_grep() {
    local pattern="$1" file="$2" msg="${3:-}"
    if ! grep -q "$pattern" "$file"; then
        echo "FAIL [$_TEST_NAME] ${msg:-pattern not found: $pattern in $file}" >&2
        echo "--- $file ---" >&2
        cat "$file" >&2
        exit 1
    fi
}

setup_tmpdir() {
    TMPDIR_TEST=$(mktemp -d)
    trap 'cleanup_tmpdir' EXIT
    cd "$TMPDIR_TEST"
}

cleanup_tmpdir() {
    if [[ -n "${TMPDIR_TEST:-}" && -d "$TMPDIR_TEST" ]]; then
        rm -rf "$TMPDIR_TEST"
    fi
}

pass() {
    echo "PASS [$_TEST_NAME]"
}
