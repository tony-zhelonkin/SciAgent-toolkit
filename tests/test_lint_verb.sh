#!/usr/bin/env bash
# tests/test_lint_verb.sh — `scio lint` dispatch and hardness behavior.
#
# Tests:
#   1. `lint --help` exits 0.
#   2. unknown `--check` name exits 1 and names the valid set.
#   3. no `--check` given → runs `all` (a planted figure-style finding shows up).
#   4. `--strict` turns a planted finding into exit 1; default is exit 0.
set -u
. "$(dirname "$0")/_lib.sh"

setup_tmpdir
FAKE="$TMPDIR_TEST/fake-toolkit"
build_fake_toolkit "$FAKE"
export SCIO_TOOLKIT="$FAKE"
SCIO="$FAKE/bin/scio"

fail() {
    echo "FAIL [$_TEST_NAME] $1" >&2
    shift
    [[ $# -gt 0 ]] && { echo "--- output ---" >&2; printf '%s\n' "$@" >&2; }
    exit 1
}

# ---------------------------------------------------------------------------
# Test 1: `lint --help` exits 0 and documents the checks.
# ---------------------------------------------------------------------------
set +e
out=$("$SCIO" lint --help 2>&1); rc=$?
set -e
[[ "$rc" -eq 0 ]] || fail "lint --help exited $rc (expected 0)" "$out"
printf '%s\n' "$out" | grep -q 'hooks' \
    || fail "lint --help omits 'hooks' from the valid --check names" "$out"
printf '%s\n' "$out" | grep -q 'toolkit' \
    || fail "lint --help omits 'toolkit' from the valid --check names" "$out"

# ---------------------------------------------------------------------------
# Test 2: unknown --check name exits 1 and names the valid set.
# ---------------------------------------------------------------------------
set +e
out=$("$SCIO" lint --check bogus-name 2>&1); rc=$?
set -e
[[ "$rc" -eq 1 ]] || fail "lint --check bogus-name exited $rc (expected 1)" "$out"
printf '%s\n' "$out" | grep -q 'valid:.*figure-style.*results-layout.*captions.*provenance.*freshness.*hooks' \
    || fail "unknown --check error does not name the valid set" "$out"

# ---------------------------------------------------------------------------
# Fixture: a project with a planted figure-style finding (raw hex color).
# ---------------------------------------------------------------------------
PROJ="$TMPDIR_TEST/proj"
mkdir -p "$PROJ/02_analysis/stages"
cat > "$PROJ/02_analysis/stages/10_qc_viz.R" <<'R'
library(ggplot2)
project_theme()
p <- ggplot(df) + scale_color_manual(values = c("#1a2b3c"))
ggsave("x.png", p)
R

# ---------------------------------------------------------------------------
# Test 3: no --check given → default is `all` (planted finding surfaces).
# ---------------------------------------------------------------------------
set +e
out=$("$SCIO" lint --project-dir "$PROJ" 2>&1); rc=$?
set -e
[[ "$rc" -eq 0 ]] || fail "lint with no --check exited $rc by default (expected soft 0)" "$out"
printf '%s\n' "$out" | grep -q 'WARN figure-style:.*raw hex color literal' \
    || fail "lint with no --check did not run figure-style (no 'all' default)" "$out"

# ---------------------------------------------------------------------------
# Test 4: --strict promotes the planted finding to exit 1; default is 0.
# ---------------------------------------------------------------------------
set +e
out=$("$SCIO" lint --check figure-style --project-dir "$PROJ" 2>&1); rc=$?
set -e
[[ "$rc" -eq 0 ]] || fail "lint --check figure-style (default) exited $rc (expected 0)" "$out"

set +e
out=$("$SCIO" lint --check figure-style --strict --project-dir "$PROJ" 2>&1); rc=$?
set -e
[[ "$rc" -eq 1 ]] || fail "lint --check figure-style --strict exited $rc (expected 1)" "$out"
printf '%s\n' "$out" | grep -q 'ERROR figure-style:.*raw hex color literal' \
    || fail "--strict finding is not an ERROR" "$out"

pass
