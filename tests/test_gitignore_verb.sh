#!/usr/bin/env bash
# tests/test_gitignore_verb.sh — `sciagent gitignore`.
#
# Defect (arch review #3): `sciagent gitignore --help` had no -h|--help case,
# so `--help` was parsed as the target path — `touch "$target"` created a
# file literally named `--help`, and a second run against an existing block
# hit `mv "$target.tmp" "$target"` == `mv "--help.tmp" "--help"`, which fails
# because a leading dash is parsed as an option, stranding the `.tmp` file. A
# fossil of exactly that (`./--help.tmp`, the SCIAGENT:GITIGNORE block body)
# was found sitting in the toolkit checkout.
#
# Tests:
#   1. `gitignore --help` in an empty dir creates NO file at all, prints
#      usage to stdout, exits 0.
#   2. `gitignore -h` — same.
#   3. Regression: normal invocation (no args) still appends the block to
#      ./.gitignore, and a second run updates it in place (no `--help`-shaped
#      fallout from the option-parsing change).
#   4. A literal `--`-prefixed target path is only honored after `--`, never
#      swallowed as an unknown option.
set -u
. "$(dirname "$0")/_lib.sh"

setup_tmpdir
FAKE="$TMPDIR_TEST/fake-toolkit"
build_fake_toolkit "$FAKE"
export SCIAGENT_TOOLKIT="$FAKE"
SCIAGENT="$FAKE/bin/sciagent"

# ---------------------------------------------------------------------------
# Test 1: --help creates no file.
# ---------------------------------------------------------------------------
mkdir "$TMPDIR_TEST/t1" && cd "$TMPDIR_TEST/t1"

before=$(find . -mindepth 1 | sort)
set +e
out1=$("$SCIAGENT" gitignore --help 2>&1)
rc1=$?
set -e
after=$(find . -mindepth 1 | sort)

if [[ "$rc1" -ne 0 ]]; then
    echo "FAIL [$_TEST_NAME] test1: --help should exit 0, got $rc1" >&2
    printf '%s\n' "$out1" >&2
    exit 1
fi
if [[ "$before" != "$after" ]]; then
    echo "FAIL [$_TEST_NAME] test1: --help must not create/modify any file" >&2
    echo "before: $before" >&2
    echo "after:  $after" >&2
    exit 1
fi
if [[ -e "--help" || -e "--help.tmp" ]]; then
    echo "FAIL [$_TEST_NAME] test1: --help must not create a file literally named --help(.tmp)" >&2
    exit 1
fi
if ! printf '%s\n' "$out1" | grep -q 'sciagent gitignore'; then
    echo "FAIL [$_TEST_NAME] test1: expected usage text on stdout" >&2
    printf '%s\n' "$out1" >&2
    exit 1
fi

# ---------------------------------------------------------------------------
# Test 2: -h — same guarantee.
# ---------------------------------------------------------------------------
mkdir "$TMPDIR_TEST/t2" && cd "$TMPDIR_TEST/t2"
before2=$(find . -mindepth 1 | sort)
"$SCIAGENT" gitignore -h >/dev/null
after2=$(find . -mindepth 1 | sort)
if [[ "$before2" != "$after2" ]]; then
    echo "FAIL [$_TEST_NAME] test2: -h must not create/modify any file" >&2
    exit 1
fi

# ---------------------------------------------------------------------------
# Test 3: regression — normal invocation still works (append then update).
# ---------------------------------------------------------------------------
mkdir "$TMPDIR_TEST/t3" && cd "$TMPDIR_TEST/t3"
"$SCIAGENT" gitignore >/dev/null
assert_file_exists .gitignore "gitignore file created on first run"
assert_grep 'BEGIN SCIAGENT:GITIGNORE' .gitignore
first_block=$(cat .gitignore)

# Second run: block already present → update path (the one that used to do
# `mv "$target.tmp" "$target"` and would break if $target ever carried a
# leading dash). Must still succeed and leave exactly one block.
"$SCIAGENT" gitignore >/dev/null
occurrences=$(grep -c 'BEGIN SCIAGENT:GITIGNORE' .gitignore)
assert_eq "$occurrences" "1" "exactly one SCIAGENT:GITIGNORE block after two runs"
assert_eq "$(cat .gitignore)" "$first_block" "idempotent re-run produces the same block"

# ---------------------------------------------------------------------------
# Test 4: `--` forces a dash-prefixed path to be treated as the target, not
# an unknown option — proves the hardening didn't just move the bug.
# ---------------------------------------------------------------------------
mkdir "$TMPDIR_TEST/t4" && cd "$TMPDIR_TEST/t4"
"$SCIAGENT" gitignore -- ./-oddname >/dev/null
assert_file_exists "./-oddname" "-- forces a dash-led path to be honored as the target"
assert_grep 'BEGIN SCIAGENT:GITIGNORE' ./-oddname

pass
