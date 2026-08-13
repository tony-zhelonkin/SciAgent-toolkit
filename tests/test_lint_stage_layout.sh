#!/usr/bin/env bash
# tests/test_lint_stage_layout.sh — opt-in `--check stage-layout` guardrail.
# Doc 09 §3.3: cheap structural companion to the stage/viz-twin/decade-number
# contract. Six rules; see lib/scio/lint.sh::_lint_check_stage_layout.
#
# Tests:
#   1. Clean layout (paired stage/viz, two-digit numbers, no strays) -> exit 0
#      under --strict, both `stages/` and `scripts/` spellings.
#   2. Rule 1: *_viz stage writing tables -> WARN default, exit 1 --strict.
#   3. Rule 2: non-viz stage writing figures -> WARN default, exit 1 --strict.
#   4. Rule 3: orphan viz (near-miss: same number different stem) -> finding.
#   5. Rule 3: near-miss (same stem different number) still orphan -> finding.
#   6. Rule 3: real pair (same number, same stem) -> NOT flagged.
#   7. Rule 4: stray file directly under 02_analysis/ -> finding; known entry
#      points (README.md, run_all.sh) do NOT fire.
#   8. Rule 5: letter-suffixed stage number (03d) -> finding.
#   9. Rule 6: single-digit stage number (4_...) -> finding.
#  10. `_scratch/` exemption: a violating file under _scratch/ is silent.
set -u
. "$(dirname "$0")/_lib.sh"

setup_tmpdir
FAKE="$TMPDIR_TEST/fake-toolkit"
build_fake_toolkit "$FAKE"
export SCIO_TOOLKIT="$FAKE"
SCIO="$FAKE/bin/scio"

# ---------------------------------------------------------------------------
# Test 1: clean layout, both stage-dir spellings, exit 0 under --strict.
# ---------------------------------------------------------------------------
for SPELLING in stages scripts; do
    P="$TMPDIR_TEST/clean-$SPELLING"
    mkdir -p "$P/02_analysis/$SPELLING" "$P/02_analysis/config" "$P/02_analysis/helpers" "$P/02_analysis/notebooks"
    cat > "$P/02_analysis/$SPELLING/10_qc.R" <<'R'
x <- 1
write.csv(x, "out.csv")
R
    cat > "$P/02_analysis/$SPELLING/10_qc_viz.R" <<'R'
library(ggplot2)
ggsave("plot.png")
R
    touch "$P/README.md"

    out=$("$SCIO" lint --check stage-layout --strict --project-dir "$P" 2>&1); rc=$?
    [[ "$rc" -eq 0 ]] || { echo "FAIL [$_TEST_NAME] test1 ($SPELLING): expected exit 0, got $rc"; printf '%s\n' "$out" >&2; exit 1; }
done

# ---------------------------------------------------------------------------
# Test 2: rule 1 — *_viz stage writes tables.
# ---------------------------------------------------------------------------
P="$TMPDIR_TEST/rule1"
mkdir -p "$P/02_analysis/stages"
cat > "$P/02_analysis/stages/10_qc.R" <<'R'
x <- 1
R
cat > "$P/02_analysis/stages/10_qc_viz.R" <<'R'
library(ggplot2)
ggsave("plot.png")
write.csv(data.frame(x=1), "sneaky_table.csv")
R
out=$("$SCIO" lint --check stage-layout --project-dir "$P" 2>&1); rc=$?
[[ "$rc" -eq 0 ]] || { echo "FAIL [$_TEST_NAME] test2: default (non-strict) must exit 0"; printf '%s\n' "$out" >&2; exit 1; }
printf '%s\n' "$out" | grep -q 'WARN stage-layout:.*10_qc_viz.R.*writes tables' \
    || { echo "FAIL [$_TEST_NAME] test2: expected WARN for viz-writes-tables"; printf '%s\n' "$out" >&2; exit 1; }
out=$("$SCIO" lint --check stage-layout --strict --project-dir "$P" 2>&1); rc=$?
[[ "$rc" -eq 1 ]] || { echo "FAIL [$_TEST_NAME] test2: --strict must exit 1"; printf '%s\n' "$out" >&2; exit 1; }
printf '%s\n' "$out" | grep -q 'ERROR stage-layout:.*10_qc_viz.R.*writes tables' \
    || { echo "FAIL [$_TEST_NAME] test2: expected ERROR under --strict"; printf '%s\n' "$out" >&2; exit 1; }

# ---------------------------------------------------------------------------
# Test 3: rule 2 — non-viz stage writes figures.
# ---------------------------------------------------------------------------
P="$TMPDIR_TEST/rule2"
mkdir -p "$P/02_analysis/stages"
cat > "$P/02_analysis/stages/10_qc.R" <<'R'
library(ggplot2)
x <- 1
ggsave("sneaky_plot.png")
R
cat > "$P/02_analysis/stages/10_qc_viz.R" <<'R'
library(ggplot2)
ggsave("plot.png")
R
out=$("$SCIO" lint --check stage-layout --project-dir "$P" 2>&1); rc=$?
[[ "$rc" -eq 0 ]] || { echo "FAIL [$_TEST_NAME] test3: default exit 0"; printf '%s\n' "$out" >&2; exit 1; }
printf '%s\n' "$out" | grep -q 'WARN stage-layout:.*10_qc\.R.*writes figures' \
    || { echo "FAIL [$_TEST_NAME] test3: expected WARN for compute-writes-figures"; printf '%s\n' "$out" >&2; exit 1; }
out=$("$SCIO" lint --check stage-layout --strict --project-dir "$P" 2>&1); rc=$?
[[ "$rc" -eq 1 ]] || { echo "FAIL [$_TEST_NAME] test3: --strict exit 1"; printf '%s\n' "$out" >&2; exit 1; }

# ---------------------------------------------------------------------------
# Test 4/5/6: rule 3 — orphan viz near-miss logic.
# ---------------------------------------------------------------------------
# 4: same number, DIFFERENT stem -> orphan (near-miss).
P="$TMPDIR_TEST/rule3-num-match"
mkdir -p "$P/02_analysis/stages"
cat > "$P/02_analysis/stages/10_qc.R" <<'R'
x <- 1
R
cat > "$P/02_analysis/stages/10_other_viz.R" <<'R'
library(ggplot2)
ggsave("plot.png")
R
out=$("$SCIO" lint --check stage-layout --project-dir "$P" 2>&1)
printf '%s\n' "$out" | grep -q 'WARN stage-layout:.*10_other_viz.R.*orphan viz' \
    || { echo "FAIL [$_TEST_NAME] test4: expected orphan-viz WARN (same number, different stem)"; printf '%s\n' "$out" >&2; exit 1; }

# 5: same stem, DIFFERENT number -> still orphan (the exact number+stem pair is absent).
P="$TMPDIR_TEST/rule3-stem-match"
mkdir -p "$P/02_analysis/stages"
cat > "$P/02_analysis/stages/11_qc.R" <<'R'
x <- 1
R
cat > "$P/02_analysis/stages/10_qc_viz.R" <<'R'
library(ggplot2)
ggsave("plot.png")
R
out=$("$SCIO" lint --check stage-layout --project-dir "$P" 2>&1)
printf '%s\n' "$out" | grep -q 'WARN stage-layout:.*10_qc_viz.R.*orphan viz' \
    || { echo "FAIL [$_TEST_NAME] test5: expected orphan-viz WARN (same stem, different number)"; printf '%s\n' "$out" >&2; exit 1; }

# 6: real pair (same number AND same stem) -> no orphan finding.
P="$TMPDIR_TEST/rule3-real-pair"
mkdir -p "$P/02_analysis/stages"
cat > "$P/02_analysis/stages/10_qc.R" <<'R'
x <- 1
R
cat > "$P/02_analysis/stages/10_qc_viz.R" <<'R'
library(ggplot2)
ggsave("plot.png")
R
out=$("$SCIO" lint --check stage-layout --strict --project-dir "$P" 2>&1); rc=$?
[[ "$rc" -eq 0 ]] || { echo "FAIL [$_TEST_NAME] test6: real pair must not warn"; printf '%s\n' "$out" >&2; exit 1; }
printf '%s\n' "$out" | grep -q 'orphan viz' \
    && { echo "FAIL [$_TEST_NAME] test6: real pair falsely flagged as orphan"; printf '%s\n' "$out" >&2; exit 1; }
true

# ---------------------------------------------------------------------------
# Test 7: rule 4 — stray file at 02_analysis/ root; known entry points exempt.
# ---------------------------------------------------------------------------
P="$TMPDIR_TEST/rule4"
mkdir -p "$P/02_analysis/stages"
echo "x = 1" > "$P/02_analysis/config.py"
echo "#!/bin/sh" > "$P/02_analysis/run_all.sh"
echo "# analysis" > "$P/02_analysis/README.md"
out=$("$SCIO" lint --check stage-layout --project-dir "$P" 2>&1)
printf '%s\n' "$out" | grep -q 'WARN stage-layout:.*config\.py.*stray file' \
    || { echo "FAIL [$_TEST_NAME] test7: expected stray-file WARN for config.py"; printf '%s\n' "$out" >&2; exit 1; }
printf '%s\n' "$out" | grep -q 'stray file.*run_all\.sh\|run_all\.sh.*stray file' \
    && { echo "FAIL [$_TEST_NAME] test7: run_all.sh is a known entry point, must not warn"; printf '%s\n' "$out" >&2; exit 1; }
printf '%s\n' "$out" | grep -q '02_analysis/README\.md: stray file' \
    && { echo "FAIL [$_TEST_NAME] test7: README.md is a known entry point, must not warn"; printf '%s\n' "$out" >&2; exit 1; }
true

# ---------------------------------------------------------------------------
# Test 8: rule 5 — letter-suffixed stage number.
# ---------------------------------------------------------------------------
P="$TMPDIR_TEST/rule5"
mkdir -p "$P/02_analysis/stages"
cat > "$P/02_analysis/stages/03d_interim.R" <<'R'
x <- 1
R
out=$("$SCIO" lint --check stage-layout --project-dir "$P" 2>&1)
printf '%s\n' "$out" | grep -q "WARN stage-layout:.*03d_interim.R.*letter suffix" \
    || { echo "FAIL [$_TEST_NAME] test8: expected letter-suffix WARN"; printf '%s\n' "$out" >&2; exit 1; }

# ---------------------------------------------------------------------------
# Test 9: rule 6 — single-digit stage number.
# ---------------------------------------------------------------------------
P="$TMPDIR_TEST/rule6"
mkdir -p "$P/02_analysis/stages"
cat > "$P/02_analysis/stages/4_gsea_set_prep.R" <<'R'
x <- 1
R
out=$("$SCIO" lint --check stage-layout --project-dir "$P" 2>&1)
printf '%s\n' "$out" | grep -q "WARN stage-layout:.*4_gsea_set_prep.R.*single-digit" \
    || { echo "FAIL [$_TEST_NAME] test9: expected single-digit WARN"; printf '%s\n' "$out" >&2; exit 1; }

# ---------------------------------------------------------------------------
# Test 10: _scratch/ exemption — a violating file under _scratch/ is silent.
# ---------------------------------------------------------------------------
P="$TMPDIR_TEST/scratch-exempt"
mkdir -p "$P/02_analysis/stages/_scratch"
cat > "$P/02_analysis/stages/_scratch/4_bad_viz.R" <<'R'
library(ggplot2)
ggsave("plot.png")
write.csv(data.frame(x=1), "sneaky.csv")
R
out=$("$SCIO" lint --check stage-layout --strict --project-dir "$P" 2>&1); rc=$?
[[ "$rc" -eq 0 ]] || { echo "FAIL [$_TEST_NAME] test10: _scratch/ must be exempt"; printf '%s\n' "$out" >&2; exit 1; }

pass
