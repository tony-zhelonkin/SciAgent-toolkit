#!/usr/bin/env bash
# tests/test_lint_comment_intent.sh — `sciagent lint --check comment-intent`
# (docs 09 §3.2). Two rules only: banner/separator comments, and mid-file
# comment runs >= 6 lines starting after line 25.
#
# Tests:
#   1. banner comment fires for '=', '-', '*', '#' fills at >= 6 chars.
#   2. a 5-char fill does NOT fire (threshold is >= 6).
#   3. a 6-line mid-file comment run (starting after line 25) fires.
#   4. a 5-line run does NOT fire.
#   5. a header block within the first 25 lines is exempt.
#   6. a clean stage produces no findings.
#   7. both `stages/` and `scripts/` spellings are scanned.
#   8. `--strict` promotes a finding to exit 1; default is exit 0 (WARN).
#   9. `_scratch/` is exempt.
set -u
. "$(dirname "$0")/_lib.sh"

setup_tmpdir
PROJ="$TMPDIR_TEST/proj"
mkdir -p "$PROJ/02_analysis/stages" "$PROJ/02_analysis/scripts"

SCIAGENT="$TOOLKIT_ROOT/bin/sciagent"

fail() {
    echo "FAIL [$_TEST_NAME] $1" >&2
    shift
    [[ $# -gt 0 ]] && { echo "--- output ---" >&2; printf '%s\n' "$@" >&2; }
    exit 1
}

run_lint() {
    local dir="$1"; shift
    set +e
    _LINT_OUT=$("$SCIAGENT" lint --check comment-intent --project-dir "$dir" "$@" 2>&1)
    _LINT_RC=$?
    set -e
}

# Helper: N numbered filler lines, e.g. "x1 <- 1" ... to push a run past
# line 25 without itself being a comment.
filler_lines() {
    local n="$1" i
    for ((i = 1; i <= n; i++)); do
        printf 'x%d <- %d\n' "$i" "$i"
    done
}

# ---------------------------------------------------------------------------
# Test 1: banner comments fire for '=', '-', '*', '#' at >= 6 chars.
# ---------------------------------------------------------------------------
cat > "$PROJ/02_analysis/stages/01_banners.R" <<'R'
x <- 1
# ======
y <- 1
# ------
z <- 1
# ******
w <- 1
# ######
v <- 1
R
run_lint "$PROJ/02_analysis/stages" --check comment-intent 2>/dev/null || true
run_lint "$PROJ"
printf '%s\n' "$_LINT_OUT" | grep -q '01_banners.R:2: banner' || fail "'=' fill (6) did not fire" "$_LINT_OUT"
printf '%s\n' "$_LINT_OUT" | grep -q '01_banners.R:4: banner' || fail "'-' fill (6) did not fire" "$_LINT_OUT"
printf '%s\n' "$_LINT_OUT" | grep -q '01_banners.R:6: banner' || fail "'*' fill (6) did not fire" "$_LINT_OUT"
printf '%s\n' "$_LINT_OUT" | grep -q '01_banners.R:8: banner' || fail "'#' fill (6) did not fire" "$_LINT_OUT"
rm "$PROJ/02_analysis/stages/01_banners.R"
pass_banner_all_fills=1

# ---------------------------------------------------------------------------
# Test 2: a 5-char fill does NOT fire (threshold is >= 6, not >= 5).
# ---------------------------------------------------------------------------
cat > "$PROJ/02_analysis/stages/02_short_fill.R" <<'R'
x <- 1
# =====
y <- 1
R
run_lint "$PROJ"
printf '%s\n' "$_LINT_OUT" | grep -q '02_short_fill.R' && fail "5-char '=' fill incorrectly fired" "$_LINT_OUT"
rm "$PROJ/02_analysis/stages/02_short_fill.R"

# ---------------------------------------------------------------------------
# Test 3 + 5: a 6-line run starting after line 25 fires; a header block in
# the first 25 lines (2 comment lines) is exempt.
# ---------------------------------------------------------------------------
{
    echo "#!/usr/bin/env Rscript"
    echo "# Stage 03: does a thing, writes 03_results/03/tables/out.csv"
    filler_lines 24
    for i in 1 2 3 4 5 6; do echo "# comment $i"; done
    echo "z <- 1"
} > "$PROJ/02_analysis/stages/03_run.R"
run_lint "$PROJ"
printf '%s\n' "$_LINT_OUT" | grep -q '03_run.R:27: mid-file comment run of 6 lines' \
    || fail "6-line run starting after line 25 did not fire at the expected line" "$_LINT_OUT"
printf '%s\n' "$_LINT_OUT" | grep -qE '03_run\.R:2:' && fail "2-line header block incorrectly flagged" "$_LINT_OUT"
rm "$PROJ/02_analysis/stages/03_run.R"

# ---------------------------------------------------------------------------
# Test 5: a >= 6-line header block starting AT OR BEFORE line 25 is exempt in
# full (not just its first-25-lines portion). This is the case that actually
# exercises the line-25 exemption — a 2-line header (test 3) is too short to
# ever hit the run-length threshold regardless of the exemption.
# ---------------------------------------------------------------------------
{
    echo "#!/usr/bin/env Rscript"
    for i in 1 2 3 4 5 6 7; do echo "# header line $i"; done
    filler_lines 20
} > "$PROJ/02_analysis/stages/05_header_run.R"
run_lint "$PROJ"
printf '%s\n' "$_LINT_OUT" | grep -q '05_header_run.R' \
    && fail "header block of 7 comment lines starting at line 2 (<=25) incorrectly flagged" "$_LINT_OUT"
rm "$PROJ/02_analysis/stages/05_header_run.R"

# ---------------------------------------------------------------------------
# Test 4b: a 6-line run that STARTS at line 20 and OVERLAPS line 25 does NOT
# fire — §3.2 says "starting after line 25", literally, not "overlapping
# line 25". Proves the check reads run_start, not any line touched by the run.
# ---------------------------------------------------------------------------
{
    filler_lines 19
    for i in 1 2 3 4 5 6; do echo "# comment $i"; done   # lines 20-25
    filler_lines 5
} > "$PROJ/02_analysis/stages/04b_overlap25.R"
run_lint "$PROJ"
printf '%s\n' "$_LINT_OUT" | grep -q '04b_overlap25.R' \
    && fail "6-line run starting at line 20 (overlapping 25) incorrectly fired" "$_LINT_OUT"
rm "$PROJ/02_analysis/stages/04b_overlap25.R"

# ---------------------------------------------------------------------------
# Test 4: a 5-line run (same position) does NOT fire.
# ---------------------------------------------------------------------------
{
    echo "#!/usr/bin/env Rscript"
    echo "# header"
    filler_lines 24
    for i in 1 2 3 4 5; do echo "# comment $i"; done
    echo "z <- 1"
} > "$PROJ/02_analysis/stages/04_run5.R"
run_lint "$PROJ"
printf '%s\n' "$_LINT_OUT" | grep -q '04_run5.R' && fail "5-line run incorrectly fired" "$_LINT_OUT"
rm "$PROJ/02_analysis/stages/04_run5.R"

# ---------------------------------------------------------------------------
# Test 6: a clean stage produces no findings.
# ---------------------------------------------------------------------------
cat > "$PROJ/02_analysis/stages/05_clean.R" <<'R'
#!/usr/bin/env Rscript
# Stage 05: clean input, writes 03_results/05_clean/tables/clean.csv
library(dplyr)
x <- 1
y <- 2
R
run_lint "$PROJ"
printf '%s\n' "$_LINT_OUT" | grep -q '05_clean.R' && fail "clean stage produced a finding" "$_LINT_OUT"
[[ "$_LINT_RC" -eq 0 ]] || fail "clean-only run exited $_LINT_RC (expected 0)" "$_LINT_OUT"
rm "$PROJ/02_analysis/stages/05_clean.R"

# ---------------------------------------------------------------------------
# Test 7: both `stages/` and `scripts/` spellings are scanned.
# ---------------------------------------------------------------------------
cat > "$PROJ/02_analysis/scripts/06_banner.py" <<'PY'
x = 1
# ======
y = 1
PY
run_lint "$PROJ"
printf '%s\n' "$_LINT_OUT" | grep -q '02_analysis/scripts/06_banner.py:2: banner' \
    || fail "scripts/ spelling was not scanned" "$_LINT_OUT"
rm "$PROJ/02_analysis/scripts/06_banner.py"

# ---------------------------------------------------------------------------
# Test 8: --strict promotes a finding to exit 1; default is exit 0.
# ---------------------------------------------------------------------------
cat > "$PROJ/02_analysis/stages/07_banner.sh" <<'SH'
#!/usr/bin/env bash
# ------
echo hi
SH
run_lint "$PROJ"
[[ "$_LINT_RC" -eq 0 ]] || fail "default (non-strict) run exited $_LINT_RC (expected 0)" "$_LINT_OUT"
printf '%s\n' "$_LINT_OUT" | grep -q 'WARN comment-intent:.*07_banner.sh:2: banner' \
    || fail "expected WARN for 07_banner.sh banner" "$_LINT_OUT"
run_lint "$PROJ" --strict
[[ "$_LINT_RC" -eq 1 ]] || fail "--strict run exited $_LINT_RC (expected 1)" "$_LINT_OUT"
printf '%s\n' "$_LINT_OUT" | grep -q 'ERROR comment-intent:.*07_banner.sh:2: banner' \
    || fail "--strict finding is not an ERROR" "$_LINT_OUT"
rm "$PROJ/02_analysis/stages/07_banner.sh"

# ---------------------------------------------------------------------------
# Test 9: `_scratch/` is exempt.
# ---------------------------------------------------------------------------
mkdir -p "$PROJ/02_analysis/stages/_scratch"
cat > "$PROJ/02_analysis/stages/_scratch/08_banner.R" <<'R'
x <- 1
# ======
y <- 1
R
run_lint "$PROJ"
printf '%s\n' "$_LINT_OUT" | grep -q '_scratch' && fail "_scratch/ file was not exempt" "$_LINT_OUT"
[[ "$_LINT_RC" -eq 0 ]] || fail "_scratch-only run exited $_LINT_RC (expected 0)" "$_LINT_OUT"
rm -rf "$PROJ/02_analysis/stages/_scratch"

pass
