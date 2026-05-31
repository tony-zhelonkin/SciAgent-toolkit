#!/usr/bin/env bash
# tests/test_validate_docs_layout.sh — _validate_docs_layout() called via
# `sciagent validate --project-dir <dir>`.
#
# Tests:
#   1. Project dir with no docs/ → exit 0, output contains "WARN docs: no docs/"
#   2. docs/_internal/ NOT gitignored in a git repo → exit nonzero (hard fail)
#   3. docs/_internal/ gitignored → check C passes (exit 0)
#   4. .md file in 03_results/ → exit 0, output contains "WARN docs: report in results"
#   5. Non-standard handoff filename → exit 0, output contains "WARN docs: non-standard handoff"
set -u
. "$(dirname "$0")/_lib.sh"

setup_tmpdir
FAKE="$TMPDIR_TEST/fake-toolkit"
build_fake_toolkit "$FAKE"
export SCIAGENT_TOOLKIT="$FAKE"
SCIAGENT="$FAKE/bin/sciagent"

# ---------------------------------------------------------------------------
# Test 1: project dir with no docs/ directory → warn, exit 0
# ---------------------------------------------------------------------------
PROJ1="$TMPDIR_TEST/proj1"
mkdir -p "$PROJ1"

set +e
out1=$("$SCIAGENT" validate --project-dir "$PROJ1" 2>&1)
rc1=$?
set -e

if [[ "$rc1" -ne 0 ]]; then
    echo "FAIL [$_TEST_NAME] test1: expected exit 0 when no docs/ present, got $rc1" >&2
    printf '%s\n' "$out1" >&2
    exit 1
fi

if ! printf '%s\n' "$out1" | grep -q 'WARN docs: no docs/'; then
    echo "FAIL [$_TEST_NAME] test1: expected 'WARN docs: no docs/' in output" >&2
    printf '%s\n' "$out1" >&2
    exit 1
fi

# ---------------------------------------------------------------------------
# Test 2: docs/_internal/ exists in a git repo but NOT gitignored → hard fail
# ---------------------------------------------------------------------------
PROJ2="$TMPDIR_TEST/proj2"
mkdir -p "$PROJ2/docs/_internal"
git -C "$PROJ2" init -q
git -C "$PROJ2" config user.email "test@example.com"
git -C "$PROJ2" config user.name "Test"
# No .gitignore at all — docs/_internal is not ignored.

set +e
out2=$("$SCIAGENT" validate --project-dir "$PROJ2" 2>&1)
rc2=$?
set -e

if [[ "$rc2" -eq 0 ]]; then
    echo "FAIL [$_TEST_NAME] test2: expected non-zero exit when _internal not gitignored, got 0" >&2
    printf '%s\n' "$out2" >&2
    exit 1
fi

# ---------------------------------------------------------------------------
# Test 3: docs/_internal/ gitignored → check C passes, exit 0
# ---------------------------------------------------------------------------
PROJ3="$TMPDIR_TEST/proj3"
mkdir -p "$PROJ3/docs/_internal"
git -C "$PROJ3" init -q
git -C "$PROJ3" config user.email "test@example.com"
git -C "$PROJ3" config user.name "Test"
printf 'docs/_internal/\n' > "$PROJ3/.gitignore"

set +e
out3=$("$SCIAGENT" validate --project-dir "$PROJ3" 2>&1)
rc3=$?
set -e

if [[ "$rc3" -ne 0 ]]; then
    echo "FAIL [$_TEST_NAME] test3: expected exit 0 when _internal is gitignored, got $rc3" >&2
    printf '%s\n' "$out3" >&2
    exit 1
fi

if printf '%s\n' "$out3" | grep -q 'ERROR docs: docs/_internal/ is NOT gitignored'; then
    echo "FAIL [$_TEST_NAME] test3: should NOT emit ERROR when _internal is gitignored" >&2
    printf '%s\n' "$out3" >&2
    exit 1
fi

# ---------------------------------------------------------------------------
# Test 4: .md file in 03_results/ → warn, exit 0
# ---------------------------------------------------------------------------
PROJ4="$TMPDIR_TEST/proj4"
mkdir -p "$PROJ4/03_results"
touch "$PROJ4/03_results/summary_report.md"

set +e
out4=$("$SCIAGENT" validate --project-dir "$PROJ4" 2>&1)
rc4=$?
set -e

if [[ "$rc4" -ne 0 ]]; then
    echo "FAIL [$_TEST_NAME] test4: expected exit 0 with md in 03_results/, got $rc4" >&2
    printf '%s\n' "$out4" >&2
    exit 1
fi

if ! printf '%s\n' "$out4" | grep -q 'WARN docs: report in results'; then
    echo "FAIL [$_TEST_NAME] test4: expected 'WARN docs: report in results' in output" >&2
    printf '%s\n' "$out4" >&2
    exit 1
fi

# ---------------------------------------------------------------------------
# Test 5: non-standard handoff filename → warn, exit 0
# ---------------------------------------------------------------------------
PROJ5="$TMPDIR_TEST/proj5"
mkdir -p "$PROJ5"
touch "$PROJ5/handoff_notes.md"
touch "$PROJ5/handoff_20260101_120000.md"  # This one is valid — should NOT warn.

set +e
out5=$("$SCIAGENT" validate --project-dir "$PROJ5" 2>&1)
rc5=$?
set -e

if [[ "$rc5" -ne 0 ]]; then
    echo "FAIL [$_TEST_NAME] test5: expected exit 0 with non-standard handoff filename, got $rc5" >&2
    printf '%s\n' "$out5" >&2
    exit 1
fi

if ! printf '%s\n' "$out5" | grep -q 'WARN docs: non-standard handoff filename'; then
    echo "FAIL [$_TEST_NAME] test5: expected 'WARN docs: non-standard handoff filename' in output" >&2
    printf '%s\n' "$out5" >&2
    exit 1
fi

# The valid filename should NOT appear in warnings.
if printf '%s\n' "$out5" | grep -q 'handoff_20260101_120000.md'; then
    echo "FAIL [$_TEST_NAME] test5: valid handoff filename incorrectly flagged as non-standard" >&2
    printf '%s\n' "$out5" >&2
    exit 1
fi

pass
