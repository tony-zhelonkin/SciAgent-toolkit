#!/usr/bin/env bash
# tests/test_validate_docs_layout.sh — docs-layout project check, now
# `scio lint --check docs-layout --project-dir <dir>` (moved out of
# the toolkit-subject catalog path.
#
# Tests:
#   1. Project dir with no docs/ at all → CLEAN no-op: exit 0, NO output
#      whatsoever (same absent-subject pattern as figure-style/results-layout/
#      captions/provenance — docs-layout audits the STRUCTURE of an existing
#      docs/ tree, it does not mandate one exist; thresholds here are
#      aspirational, not baseline). This also covers the false positive
#      caught by test_validate_figure_style.sh/test_validate_provenance.sh/
#      test_validate_results_layout.sh's `--check all --strict` "conformant"
#      fixtures, which have 02_analysis/03_results but no docs/ at all.
#   1b. Same, but explicitly through `--check all --strict` on a minimal
#       analysis-shaped project with no docs/ — pins that docs-layout being a
#       full member of `all` again does not reintroduce the false positive.
#   2. docs/_internal/ NOT gitignored, no --strict → soft WARN only, exit 0
#      (severity is now consistent with every other opt-in lint check: WARN
#      unless --strict, not an unconditional hard-fail like the old
#      project-check behavior).
#   2b. Same fixture with --strict → hard ERROR, exit nonzero.
#   3. docs/_internal/ gitignored → no ERROR/WARN about it, exit 0 (even
#      under --strict).
#   4. .md file in 03_results/ → exit 0, output contains "WARN docs-layout: report in results"
#   5. Non-standard handoff filename → exit 0, output contains "WARN docs-layout: non-standard handoff"
#   6. `lint --check toolkit` against the same not-gitignored fixture remains
#      silent about docs-layout and exits 0.
set -u
. "$(dirname "$0")/_lib.sh"

setup_tmpdir
FAKE="$TMPDIR_TEST/fake-toolkit"
build_fake_toolkit "$FAKE"
export SCIO_TOOLKIT="$FAKE"
SCIO="$FAKE/bin/scio"

# ---------------------------------------------------------------------------
# Test 1: project dir with no docs/ directory at all → CLEAN no-op.
# ---------------------------------------------------------------------------
PROJ1="$TMPDIR_TEST/proj1"
mkdir -p "$PROJ1"

set +e
out1=$("$SCIO" lint --check docs-layout --project-dir "$PROJ1" 2>&1)
rc1=$?
set -e

if [[ "$rc1" -ne 0 ]]; then
    echo "FAIL [$_TEST_NAME] test1: expected exit 0 when no docs/ present, got $rc1" >&2
    printf '%s\n' "$out1" >&2
    exit 1
fi
if [[ -n "$out1" ]]; then
    echo "FAIL [$_TEST_NAME] test1: expected NO output at all when docs/ is entirely absent" >&2
    printf '%s\n' "$out1" >&2
    exit 1
fi

# Even --strict must stay clean — an absent subject is not a finding.
set +e
out1s=$("$SCIO" lint --check docs-layout --strict --project-dir "$PROJ1" 2>&1)
rc1s=$?
set -e
if [[ "$rc1s" -ne 0 || -n "$out1s" ]]; then
    echo "FAIL [$_TEST_NAME] test1: --strict must also no-op when docs/ is absent (rc=$rc1s)" >&2
    printf '%s\n' "$out1s" >&2
    exit 1
fi

# ---------------------------------------------------------------------------
# Test 1b: the exact shape of the false positive flagged in review —
# `lint --check all --strict` on a minimal analysis-shaped project (has
# 02_analysis/ + 03_results/, per the OTHER --check-all fixtures) but no
# docs/ tree at all — must stay silent about docs-layout and exit 0.
# ---------------------------------------------------------------------------
PROJ1B="$TMPDIR_TEST/proj1b"
mkdir -p "$PROJ1B/02_analysis" "$PROJ1B/03_results"

set +e
out1b=$("$SCIO" lint --check all --strict --project-dir "$PROJ1B" 2>&1)
rc1b=$?
set -e
if [[ "$rc1b" -ne 0 ]]; then
    echo "FAIL [$_TEST_NAME] test1b: --check all --strict on a docs-less analysis project should stay clean, got exit $rc1b" >&2
    printf '%s\n' "$out1b" >&2
    exit 1
fi
if printf '%s\n' "$out1b" | grep -qi 'docs-layout'; then
    echo "FAIL [$_TEST_NAME] test1b: docs-layout produced a finding on a project with no docs/ at all" >&2
    printf '%s\n' "$out1b" >&2
    exit 1
fi

# ---------------------------------------------------------------------------
# Test 2 / 2b: docs/_internal/ exists in a git repo but NOT gitignored.
# Without --strict: soft WARN, exit 0. With --strict: hard ERROR, exit nonzero.
# ---------------------------------------------------------------------------
PROJ2="$TMPDIR_TEST/proj2"
mkdir -p "$PROJ2/docs/_internal"
git -C "$PROJ2" init -q
git -C "$PROJ2" config user.email "test@example.com"
git -C "$PROJ2" config user.name "Test"
# No .gitignore at all — docs/_internal is not ignored.

set +e
out2=$("$SCIO" lint --check docs-layout --project-dir "$PROJ2" 2>&1)
rc2=$?
set -e

if [[ "$rc2" -ne 0 ]]; then
    echo "FAIL [$_TEST_NAME] test2: expected exit 0 (soft-warn, no --strict), got $rc2" >&2
    printf '%s\n' "$out2" >&2
    exit 1
fi
if ! printf '%s\n' "$out2" | grep -q 'WARN docs-layout: docs/_internal/ is NOT gitignored'; then
    echo "FAIL [$_TEST_NAME] test2: expected soft WARN about ungitignored docs/_internal/" >&2
    printf '%s\n' "$out2" >&2
    exit 1
fi

set +e
out2b=$("$SCIO" lint --check docs-layout --project-dir "$PROJ2" --strict 2>&1)
rc2b=$?
set -e

if [[ "$rc2b" -eq 0 ]]; then
    echo "FAIL [$_TEST_NAME] test2b: expected non-zero exit under --strict when _internal not gitignored, got 0" >&2
    printf '%s\n' "$out2b" >&2
    exit 1
fi
if ! printf '%s\n' "$out2b" | grep -q 'ERROR docs-layout: docs/_internal/ is NOT gitignored'; then
    echo "FAIL [$_TEST_NAME] test2b: expected hard ERROR about ungitignored docs/_internal/ under --strict" >&2
    printf '%s\n' "$out2b" >&2
    exit 1
fi

# ---------------------------------------------------------------------------
# Test 3: docs/_internal/ gitignored → the gitignore rule passes, exit 0,
# even under --strict.
# ---------------------------------------------------------------------------
PROJ3="$TMPDIR_TEST/proj3"
mkdir -p "$PROJ3/docs/_internal"
git -C "$PROJ3" init -q
git -C "$PROJ3" config user.email "test@example.com"
git -C "$PROJ3" config user.name "Test"
printf 'docs/_internal/\n' > "$PROJ3/.gitignore"

set +e
out3=$("$SCIO" lint --check docs-layout --project-dir "$PROJ3" --strict 2>&1)
rc3=$?
set -e

if [[ "$rc3" -ne 0 ]]; then
    echo "FAIL [$_TEST_NAME] test3: expected exit 0 when _internal is gitignored, got $rc3" >&2
    printf '%s\n' "$out3" >&2
    exit 1
fi

if printf '%s\n' "$out3" | grep -q 'docs/_internal/ is NOT gitignored'; then
    echo "FAIL [$_TEST_NAME] test3: should NOT flag gitignored docs/_internal/" >&2
    printf '%s\n' "$out3" >&2
    exit 1
fi

# ---------------------------------------------------------------------------
# Test 4: .md file in 03_results/ → warn, exit 0. Needs a docs/ dir present
# (even empty) — the absent-subject guard (test 1) means this check no-ops
# entirely without one.
# ---------------------------------------------------------------------------
PROJ4="$TMPDIR_TEST/proj4"
mkdir -p "$PROJ4/03_results" "$PROJ4/docs"
touch "$PROJ4/03_results/summary_report.md"

set +e
out4=$("$SCIO" lint --check docs-layout --project-dir "$PROJ4" 2>&1)
rc4=$?
set -e

if [[ "$rc4" -ne 0 ]]; then
    echo "FAIL [$_TEST_NAME] test4: expected exit 0 with md in 03_results/, got $rc4" >&2
    printf '%s\n' "$out4" >&2
    exit 1
fi

if ! printf '%s\n' "$out4" | grep -q 'WARN docs-layout: report in results'; then
    echo "FAIL [$_TEST_NAME] test4: expected 'WARN docs-layout: report in results' in output" >&2
    printf '%s\n' "$out4" >&2
    exit 1
fi

# ---------------------------------------------------------------------------
# Test 5: non-standard handoff filename → warn, exit 0. Needs a docs/ dir
# present for the same reason as test 4.
# ---------------------------------------------------------------------------
PROJ5="$TMPDIR_TEST/proj5"
mkdir -p "$PROJ5/docs"
touch "$PROJ5/handoff_notes.md"
touch "$PROJ5/handoff_20260101_120000.md"  # This one is valid — should NOT warn.

set +e
out5=$("$SCIO" lint --check docs-layout --project-dir "$PROJ5" 2>&1)
rc5=$?
set -e

if [[ "$rc5" -ne 0 ]]; then
    echo "FAIL [$_TEST_NAME] test5: expected exit 0 with non-standard handoff filename, got $rc5" >&2
    printf '%s\n' "$out5" >&2
    exit 1
fi

if ! printf '%s\n' "$out5" | grep -q 'WARN docs-layout: non-standard handoff filename'; then
    echo "FAIL [$_TEST_NAME] test5: expected 'WARN docs-layout: non-standard handoff filename' in output" >&2
    printf '%s\n' "$out5" >&2
    exit 1
fi

# The valid filename should NOT appear in warnings.
if printf '%s\n' "$out5" | grep -q 'handoff_20260101_120000.md'; then
    echo "FAIL [$_TEST_NAME] test5: valid handoff filename incorrectly flagged as non-standard" >&2
    printf '%s\n' "$out5" >&2
    exit 1
fi

# ---------------------------------------------------------------------------
# Test 6: the toolkit-subject check is silent about project docs layout.
# ---------------------------------------------------------------------------
set +e
out6=$("$SCIO" lint --check toolkit --project-dir "$PROJ2" 2>&1)
rc6=$?
set -e

if [[ "$rc6" -ne 0 ]]; then
    echo "FAIL [$_TEST_NAME] test6: expected exit 0 from toolkit check, got $rc6" >&2
    printf '%s\n' "$out6" >&2
    exit 1
fi
if printf '%s\n' "$out6" | grep -qi 'docs-layout\|docs/_internal'; then
    echo "FAIL [$_TEST_NAME] test6: toolkit check mentioned docs-layout" >&2
    printf '%s\n' "$out6" >&2
    exit 1
fi

pass
