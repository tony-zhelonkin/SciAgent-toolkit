#!/usr/bin/env bash
# tests/test_lint_verb.sh — `sciagent lint`, the extracted (c) GUARDRAIL verb.
#
# validate.sh used to own both the toolkit-wide skill-frontmatter walk AND the
# opt-in per-project guardrail checks (--check). This test covers the new
# standalone `sciagent lint` surface + the `sciagent validate --check`
# backward-compat delegation (deprecation note + still-functional checks).
#
# Tests:
#   1. `lint --help` exits 0.
#   2. unknown `--check` name exits 1 and names the valid set.
#   3. no `--check` given → runs `all` (a planted figure-style finding shows up).
#   4. `--strict` turns a planted finding into exit 1; default is exit 0.
#   5. `validate --check` still works and emits the deprecation note.
#   6. `--quiet` suppresses the deprecation note.
set -u
. "$(dirname "$0")/_lib.sh"

setup_tmpdir
FAKE="$TMPDIR_TEST/fake-toolkit"
build_fake_toolkit "$FAKE"
export SCIAGENT_TOOLKIT="$FAKE"
SCIAGENT="$FAKE/bin/sciagent"

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
out=$("$SCIAGENT" lint --help 2>&1); rc=$?
set -e
[[ "$rc" -eq 0 ]] || fail "lint --help exited $rc (expected 0)" "$out"
printf '%s\n' "$out" | grep -q 'hooks' \
    || fail "lint --help omits 'hooks' from the valid --check names" "$out"

# ---------------------------------------------------------------------------
# Test 2: unknown --check name exits 1 and names the valid set.
# ---------------------------------------------------------------------------
set +e
out=$("$SCIAGENT" lint --check bogus-name 2>&1); rc=$?
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
out=$("$SCIAGENT" lint --project-dir "$PROJ" 2>&1); rc=$?
set -e
[[ "$rc" -eq 0 ]] || fail "lint with no --check exited $rc by default (expected soft 0)" "$out"
printf '%s\n' "$out" | grep -q 'WARN figure-style:.*raw hex color literal' \
    || fail "lint with no --check did not run figure-style (no 'all' default)" "$out"

# ---------------------------------------------------------------------------
# Test 4: --strict promotes the planted finding to exit 1; default is 0.
# ---------------------------------------------------------------------------
set +e
out=$("$SCIAGENT" lint --check figure-style --project-dir "$PROJ" 2>&1); rc=$?
set -e
[[ "$rc" -eq 0 ]] || fail "lint --check figure-style (default) exited $rc (expected 0)" "$out"

set +e
out=$("$SCIAGENT" lint --check figure-style --strict --project-dir "$PROJ" 2>&1); rc=$?
set -e
[[ "$rc" -eq 1 ]] || fail "lint --check figure-style --strict exited $rc (expected 1)" "$out"
printf '%s\n' "$out" | grep -q 'ERROR figure-style:.*raw hex color literal' \
    || fail "--strict finding is not an ERROR" "$out"

# ---------------------------------------------------------------------------
# Test 5: `validate --check` still works and emits the deprecation note.
# ---------------------------------------------------------------------------
set +e
out=$("$SCIAGENT" validate --check figure-style --project-dir "$PROJ" 2>&1); rc=$?
set -e
[[ "$rc" -eq 0 ]] || fail "validate --check figure-style exited $rc (expected 0)" "$out"
printf '%s\n' "$out" | grep -q 'WARN figure-style:.*raw hex color literal' \
    || fail "validate --check delegation did not run the check" "$out"
printf '%s\n' "$out" | grep -q 'sciagent validate --check is deprecated; use: sciagent lint --check' \
    || fail "validate --check did not emit the deprecation note" "$out"

# ---------------------------------------------------------------------------
# Test 6: --quiet suppresses the deprecation note (but not a --strict ERROR).
# ---------------------------------------------------------------------------
set +e
out=$("$SCIAGENT" validate --check figure-style --quiet --project-dir "$PROJ" 2>&1); rc=$?
set -e
printf '%s\n' "$out" | grep -q 'deprecated' \
    && fail "--quiet did not suppress the deprecation note" "$out"

set +e
out=$("$SCIAGENT" validate --check figure-style --strict --quiet --project-dir "$PROJ" 2>&1); rc=$?
set -e
[[ "$rc" -eq 1 ]] || fail "validate --check --strict --quiet exited $rc (expected 1)" "$out"
printf '%s\n' "$out" | grep -q 'ERROR figure-style:' \
    || fail "--quiet suppressed a --strict ERROR (must never happen)" "$out"

pass
