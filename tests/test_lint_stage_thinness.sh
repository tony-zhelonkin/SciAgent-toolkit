#!/usr/bin/env bash
# tests/test_lint_stage_thinness.sh — `scio lint --check stage-thinness`.
#
# Covers docs 09 §3.1: def-count (>2), summed function-body lines (>60), total
# stage LOC (>500), the `# stage-detail: <reason>` escape hatch (and its
# non-empty-reason requirement), the `helpers/` exemption, both
# `02_analysis/stages` and `02_analysis/scripts` spellings, `--strict`
# promotion, and `_scratch/` exemption.
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

run_check() {
    # run_check <projdir> [extra args...]
    local projdir="$1"; shift
    "$SCIO" lint --check stage-thinness --project-dir "$projdir" "$@" 2>&1
}

# ---------------------------------------------------------------------------
# Test 1: def-count rule fires (>2 defs), clean stage does not.
# ---------------------------------------------------------------------------
PROJ1="$TMPDIR_TEST/proj1"
mkdir -p "$PROJ1/02_analysis/stages"
cat > "$PROJ1/02_analysis/stages/01_three_defs.R" <<'R'
a <- function(x) {
  x + 1
}
b <- function(x) {
  x + 2
}
c <- function(x) {
  x + 3
}
R
cat > "$PROJ1/02_analysis/stages/02_clean.R" <<'R'
d <- function(x) {
  x + 1
}
R

set +e
out=$(run_check "$PROJ1"); rc=$?
set -e
[[ "$rc" -eq 0 ]] || fail "def-count finding exited $rc (expected soft 0)" "$out"
printf '%s\n' "$out" | grep -q '01_three_defs.R: 3 function definitions (> 2)' \
    || fail "def-count rule did not fire on 3-def stage" "$out"
printf '%s\n' "$out" | grep -q '02_clean.R' \
    && fail "clean 1-def stage incorrectly flagged" "$out"

# ---------------------------------------------------------------------------
# Test 2: helpers/ exemption — same 3-def content under helpers/ is silent.
# ---------------------------------------------------------------------------
PROJ2="$TMPDIR_TEST/proj2"
mkdir -p "$PROJ2/02_analysis/helpers"
cp "$PROJ1/02_analysis/stages/01_three_defs.R" "$PROJ2/02_analysis/helpers/util.R"
set +e
out=$(run_check "$PROJ2"); rc=$?
set -e
[[ "$rc" -eq 0 ]] || fail "helpers-only project exited $rc (expected 0)" "$out"
[[ -z "$out" ]] || fail "helpers/ is not exempt — got a finding" "$out"

# ---------------------------------------------------------------------------
# Test 3: escape hatch. Same 3-def file, but one def carries a non-empty
# `# stage-detail: <reason>` on the line above it -> only 2 non-exempt defs
# remain -> def-count rule must NOT fire. A sibling file with an EMPTY reason
# on the same escape-hatch comment must still fire (reason is required).
# ---------------------------------------------------------------------------
PROJ3="$TMPDIR_TEST/proj3"
mkdir -p "$PROJ3/02_analysis/stages"
cat > "$PROJ3/02_analysis/stages/03_hatch.R" <<'R'
# stage-detail: training spec is the narrative here, keep it inline
a <- function(x) {
  x + 1
}
b <- function(x) {
  x + 2
}
c <- function(x) {
  x + 3
}
R
cat > "$PROJ3/02_analysis/stages/04_empty_reason.R" <<'R'
# stage-detail:
a <- function(x) { x + 1 }
b <- function(x) { x + 1 }
c <- function(x) { x + 1 }
R

set +e
out=$(run_check "$PROJ3"); rc=$?
set -e
[[ "$rc" -eq 0 ]] || fail "escape-hatch project exited $rc (expected 0)" "$out"
printf '%s\n' "$out" | grep -q '03_hatch.R' \
    && fail "non-empty stage-detail reason did not exempt its definition" "$out"
printf '%s\n' "$out" | grep -q '04_empty_reason.R: 3 function definitions (> 2)' \
    || fail "an EMPTY stage-detail reason wrongly exempted (must require non-empty reason)" "$out"

# ---------------------------------------------------------------------------
# Test 4: summed function-body-lines rule (>60), independent of def-count.
# ---------------------------------------------------------------------------
PROJ4="$TMPDIR_TEST/proj4"
mkdir -p "$PROJ4/02_analysis/stages"
{
    echo 'big <- function(x) {'
    for i in $(seq 1 65); do echo "  y$i <- x"; done
    echo '}'
} > "$PROJ4/02_analysis/stages/05_bigbody.R"

set +e
out=$(run_check "$PROJ4"); rc=$?
set -e
[[ "$rc" -eq 0 ]] || fail "body-lines finding exited $rc (expected soft 0)" "$out"
printf '%s\n' "$out" | grep -q '05_bigbody.R: 66 lines inside function bodies (> 60)' \
    || fail "function-body-lines rule did not fire" "$out"
printf '%s\n' "$out" | grep -q 'function definitions (> 2)' \
    && fail "single-def file wrongly triggered the def-count rule too" "$out"

# ---------------------------------------------------------------------------
# Test 5: total-stage-LOC rule (>500), independent of defs/bodies.
# ---------------------------------------------------------------------------
PROJ5="$TMPDIR_TEST/proj5"
mkdir -p "$PROJ5/02_analysis/stages"
for i in $(seq 1 510); do echo "# line $i"; done > "$PROJ5/02_analysis/stages/06_biglines.R"

set +e
out=$(run_check "$PROJ5"); rc=$?
set -e
[[ "$rc" -eq 0 ]] || fail "LOC finding exited $rc (expected soft 0)" "$out"
printf '%s\n' "$out" | grep -q '06_biglines.R: 510 lines (> 500 LOC)' \
    || fail "total-LOC rule did not fire" "$out"

# ---------------------------------------------------------------------------
# Test 6: both stage-dir spellings fire — 02_analysis/scripts (migration
# window) as well as 02_analysis/stages (canonical).
# ---------------------------------------------------------------------------
PROJ6="$TMPDIR_TEST/proj6"
mkdir -p "$PROJ6/02_analysis/scripts"
cp "$PROJ1/02_analysis/stages/01_three_defs.R" "$PROJ6/02_analysis/scripts/01_three_defs.R"

set +e
out=$(run_check "$PROJ6"); rc=$?
set -e
[[ "$rc" -eq 0 ]] || fail "scripts/ spelling exited $rc (expected 0)" "$out"
printf '%s\n' "$out" | grep -q '02_analysis/scripts/01_three_defs.R: 3 function definitions (> 2)' \
    || fail "02_analysis/scripts/ (pre-rename spelling) was not scanned" "$out"

# ---------------------------------------------------------------------------
# Test 7: --strict promotes a finding to a hard failure (exit 1, ERROR).
# ---------------------------------------------------------------------------
set +e
out=$(run_check "$PROJ1" --strict); rc=$?
set -e
[[ "$rc" -eq 1 ]] || fail "--strict did not promote to exit 1" "$out"
printf '%s\n' "$out" | grep -q 'ERROR stage-thinness:.*01_three_defs.R: 3 function definitions (> 2)' \
    || fail "--strict finding is not an ERROR" "$out"

# ---------------------------------------------------------------------------
# Test 8: _scratch/ exemption.
# ---------------------------------------------------------------------------
PROJ8="$TMPDIR_TEST/proj8"
mkdir -p "$PROJ8/02_analysis/stages/_scratch"
cp "$PROJ1/02_analysis/stages/01_three_defs.R" "$PROJ8/02_analysis/stages/_scratch/01_three_defs.R"

set +e
out=$(run_check "$PROJ8"); rc=$?
set -e
[[ "$rc" -eq 0 ]] || fail "_scratch/ project exited $rc (expected 0)" "$out"
[[ -z "$out" ]] || fail "_scratch/ is not exempt — got a finding" "$out"

# ---------------------------------------------------------------------------
# Boundary-exact tests. These are the primary correctness evidence in place
# of empirical calibration against a real corpus (mouse_anchor/human_treg_
# arthritis are not checked out on this filesystem — see the task report).
# Each rule's threshold is a strict ">" — exactly-at-threshold must NOT fire,
# one-past-threshold MUST fire.
# ---------------------------------------------------------------------------

# --- Test 9: def-count boundary — exactly 2 defs clean, 3 defs fires. ------
PROJ9="$TMPDIR_TEST/proj9"
mkdir -p "$PROJ9/02_analysis/stages"
cat > "$PROJ9/02_analysis/stages/09a_at_two.R" <<'R'
a <- function(x) {
  x + 1
}
b <- function(x) {
  x + 2
}
R
cat > "$PROJ9/02_analysis/stages/09b_at_three.R" <<'R'
a <- function(x) {
  x + 1
}
b <- function(x) {
  x + 2
}
c <- function(x) {
  x + 3
}
R

set +e
out=$(run_check "$PROJ9"); rc=$?
set -e
[[ "$rc" -eq 0 ]] || fail "def-count boundary project exited $rc (expected 0)" "$out"
printf '%s\n' "$out" | grep -q '09a_at_two.R' \
    && fail "exactly 2 defs must NOT fire the def-count rule (threshold is > 2)" "$out"
printf '%s\n' "$out" | grep -q '09b_at_three.R: 3 function definitions (> 2)' \
    || fail "exactly 3 defs (one past threshold) must fire the def-count rule" "$out"

# --- Test 10: body-lines boundary — a function whose brace-balance body is
# exactly 60 lines is clean; 61 lines fires. Bodies constructed as N
# assignment statements between `{` (on the def line) and `}`: body-line
# count = N + 1 (see _vcheck_stage_scan_defs's brace-balance comment), so
# N=59 -> 60 (at threshold, clean) and N=60 -> 61 (fires).
PROJ10="$TMPDIR_TEST/proj10"
mkdir -p "$PROJ10/02_analysis/stages"
{
    echo 'f <- function(x) {'
    for i in $(seq 1 59); do echo "  y$i <- x"; done
    echo '}'
} > "$PROJ10/02_analysis/stages/10a_at_sixty.R"
{
    echo 'f <- function(x) {'
    for i in $(seq 1 60); do echo "  y$i <- x"; done
    echo '}'
} > "$PROJ10/02_analysis/stages/10b_at_sixtyone.R"

set +e
out=$(run_check "$PROJ10"); rc=$?
set -e
[[ "$rc" -eq 0 ]] || fail "body-lines boundary project exited $rc (expected 0)" "$out"
printf '%s\n' "$out" | grep -q '10a_at_sixty.R' \
    && fail "exactly 60 body lines must NOT fire the body-lines rule (threshold is > 60)" "$out"
printf '%s\n' "$out" | grep -q '10b_at_sixtyone.R: 61 lines inside function bodies (> 60)' \
    || fail "exactly 61 body lines (one past threshold) must fire the body-lines rule" "$out"

# --- Test 11: total-LOC boundary — a 500-line stage is clean; 501 fires. ---
PROJ11="$TMPDIR_TEST/proj11"
mkdir -p "$PROJ11/02_analysis/stages"
for i in $(seq 1 500); do echo "# line $i"; done > "$PROJ11/02_analysis/stages/11a_at_500.R"
for i in $(seq 1 501); do echo "# line $i"; done > "$PROJ11/02_analysis/stages/11b_at_501.R"

set +e
out=$(run_check "$PROJ11"); rc=$?
set -e
[[ "$rc" -eq 0 ]] || fail "LOC boundary project exited $rc (expected 0)" "$out"
printf '%s\n' "$out" | grep -q '11a_at_500.R' \
    && fail "exactly 500 LOC must NOT fire the total-LOC rule (threshold is > 500)" "$out"
printf '%s\n' "$out" | grep -q '11b_at_501.R: 501 lines (> 500 LOC)' \
    || fail "exactly 501 LOC (one past threshold) must fire the total-LOC rule" "$out"

pass
