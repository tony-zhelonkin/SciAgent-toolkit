#!/usr/bin/env bash
# tests/test_validate_figure_style.sh — opt-in `--check figure-style` guardrail.
#
# Tests:
#   1. CONFORMANT project (base_size 16, viz uses set_paper_style + save_overview,
#      no inline styling) → exit 0 default AND exit 0 under --strict.
#   2. VIOLATION project (base_size 12, inline theme()/ggsave(width=)/raw hex,
#      no theme call) → default exit 0 with WARN on stderr; --strict exit 1.
#   3. No-false-positive: `validate --check all` on the conformant project → 0.
#   4. Software / empty project (no analysis layout) → exit 0 (clean no-op).
#   5. Activate-internal path is unaffected: plain `validate` ignores --check.
set -u
. "$(dirname "$0")/_lib.sh"

setup_tmpdir
FAKE="$TMPDIR_TEST/fake-toolkit"
build_fake_toolkit "$FAKE"
export SCIAGENT_TOOLKIT="$FAKE"
SCIAGENT="$FAKE/bin/sciagent"

# ---------------------------------------------------------------------------
# Helper: build a CONFORMANT project under $1.
# ---------------------------------------------------------------------------
make_conformant() {
    local p="$1"
    mkdir -p "$p/02_analysis/config" "$p/02_analysis/scripts"
    mkdir -p "$p/03_results/01_qc/figures/_overview" "$p/03_results/01_qc/tables/_overview"
    cat > "$p/02_analysis/config/analysis_config.yaml" <<'YAML'
figures:
  base_size: 16
  top_n: 20
  by_contrast_dir: "by_contrast"
  overview_dir: "_overview"
stages:
  - id: "01_qc"
    title: "QC"
paths:
  results: "03_results/"
YAML
    cat > "$p/02_analysis/scripts/10_qc_viz.py" <<'PY'
from figure_helpers import set_paper_style, save_overview, load_figure_config
cfg = load_figure_config()
set_paper_style(config=cfg)
save_overview(fig, "01_qc", "qc_counts", table=rows, config=cfg)
PY
    touch "$p/03_results/01_qc/figures/_overview/qc_counts.png"
    touch "$p/03_results/01_qc/tables/_overview/qc_counts.csv"
    cat > "$p/03_results/01_qc/README.md" <<'MD'
# 01_qc — artifact captions

## figures/_overview/qc_counts.png

Counts per cell.

**How to read:** higher = more reads.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/10_qc_viz.py` | `save_overview` | `figures.top_n = 20` | `03_results/objects/qc.h5ad` |

## tables/_overview/qc_counts.csv

Source table.

**How to read:** raw counts.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/10_qc_viz.py` | `save_overview` | `figures.top_n = 20` | `03_results/objects/qc.h5ad` |
MD
}

# ---------------------------------------------------------------------------
# Test 1: conformant project passes default + strict.
# ---------------------------------------------------------------------------
CONF="$TMPDIR_TEST/conf"
make_conformant "$CONF"

set +e
out=$("$SCIAGENT" validate --check figure-style --project-dir "$CONF" 2>&1); rc=$?
set -e
[[ "$rc" -eq 0 ]] || { echo "FAIL [$_TEST_NAME] test1: conformant default exit $rc"; printf '%s\n' "$out" >&2; exit 1; }

set +e
out=$("$SCIAGENT" validate --check figure-style --strict --project-dir "$CONF" 2>&1); rc=$?
set -e
[[ "$rc" -eq 0 ]] || { echo "FAIL [$_TEST_NAME] test1: conformant strict exit $rc"; printf '%s\n' "$out" >&2; exit 1; }

# ---------------------------------------------------------------------------
# Test 2: violation project — default warns (exit 0), strict fails (exit 1).
# ---------------------------------------------------------------------------
VIOL="$TMPDIR_TEST/viol"
mkdir -p "$VIOL/02_analysis/config" "$VIOL/02_analysis/scripts"
cat > "$VIOL/02_analysis/config/analysis_config.yaml" <<'YAML'
figures:
  base_size: 12
stages:
  - id: "01_qc"
YAML
cat > "$VIOL/02_analysis/scripts/10_qc_viz.R" <<'R'
library(ggplot2)
p <- ggplot(df) + theme(text = element_text(size = 8))
ggsave("x.png", p, width = 7, height = 5)
col <- "#FF00AA"
R

set +e
out=$("$SCIAGENT" validate --check figure-style --project-dir "$VIOL" 2>&1); rc=$?
set -e
[[ "$rc" -eq 0 ]] || { echo "FAIL [$_TEST_NAME] test2: violation default expected exit 0 got $rc"; printf '%s\n' "$out" >&2; exit 1; }
printf '%s\n' "$out" | grep -q 'WARN figure-style:.*base_size = 12' \
    || { echo "FAIL [$_TEST_NAME] test2: expected base_size WARN"; printf '%s\n' "$out" >&2; exit 1; }
printf '%s\n' "$out" | grep -q 'WARN figure-style:.*ggsave(width' \
    || { echo "FAIL [$_TEST_NAME] test2: expected ggsave(width=) WARN"; printf '%s\n' "$out" >&2; exit 1; }
printf '%s\n' "$out" | grep -q 'WARN figure-style:.*hex color' \
    || { echo "FAIL [$_TEST_NAME] test2: expected raw hex WARN"; printf '%s\n' "$out" >&2; exit 1; }
printf '%s\n' "$out" | grep -q 'WARN figure-style:.*without calling project_theme' \
    || { echo "FAIL [$_TEST_NAME] test2: expected missing-theme WARN"; printf '%s\n' "$out" >&2; exit 1; }

set +e
out=$("$SCIAGENT" validate --check figure-style --strict --project-dir "$VIOL" 2>&1); rc=$?
set -e
[[ "$rc" -eq 1 ]] || { echo "FAIL [$_TEST_NAME] test2: violation strict expected exit 1 got $rc"; printf '%s\n' "$out" >&2; exit 1; }
printf '%s\n' "$out" | grep -q 'ERROR figure-style:' \
    || { echo "FAIL [$_TEST_NAME] test2: expected ERROR under strict"; printf '%s\n' "$out" >&2; exit 1; }

# Default WARNs go to STDERR (so a clean stdout stays usable by callers).
set +e
serr=$("$SCIAGENT" validate --check figure-style --project-dir "$VIOL" 2>&1 1>/dev/null)
set -e
printf '%s\n' "$serr" | grep -q 'WARN figure-style:' \
    || { echo "FAIL [$_TEST_NAME] test2: WARN must be on stderr"; printf '%s\n' "$serr" >&2; exit 1; }

# ---------------------------------------------------------------------------
# Test 3: no-false-positive — `--check all` on conformant exits 0.
# ---------------------------------------------------------------------------
set +e
out=$("$SCIAGENT" validate --check all --project-dir "$CONF" 2>&1); rc=$?
set -e
[[ "$rc" -eq 0 ]] || { echo "FAIL [$_TEST_NAME] test3: --check all conformant exit $rc"; printf '%s\n' "$out" >&2; exit 1; }
set +e
out=$("$SCIAGENT" validate --check all --strict --project-dir "$CONF" 2>&1); rc=$?
set -e
[[ "$rc" -eq 0 ]] || { echo "FAIL [$_TEST_NAME] test3: --check all --strict conformant exit $rc"; printf '%s\n' "$out" >&2; exit 1; }

# ---------------------------------------------------------------------------
# Test 4: software / empty project produces no findings (clean no-op).
# ---------------------------------------------------------------------------
SW="$TMPDIR_TEST/sw"
mkdir -p "$SW/src"; echo 'int main(){return 0;}' > "$SW/src/main.c"
set +e
out=$("$SCIAGENT" validate --check all --strict --project-dir "$SW" 2>&1); rc=$?
set -e
[[ "$rc" -eq 0 ]] || { echo "FAIL [$_TEST_NAME] test4: software project exit $rc"; printf '%s\n' "$out" >&2; exit 1; }

# ---------------------------------------------------------------------------
# Test 5: plain `validate` (no --check) ignores project checks entirely.
# Even with a flagrant figure-style violation present, plain validate must
# pass (it never runs the opt-in checks) — this is the activate-internal path.
# ---------------------------------------------------------------------------
set +e
out=$("$SCIAGENT" validate --project-dir "$VIOL" 2>&1); rc=$?
set -e
[[ "$rc" -eq 0 ]] || { echo "FAIL [$_TEST_NAME] test5: plain validate over violation expected exit 0 got $rc"; printf '%s\n' "$out" >&2; exit 1; }
printf '%s\n' "$out" | grep -q 'figure-style' \
    && { echo "FAIL [$_TEST_NAME] test5: plain validate must NOT run figure-style"; printf '%s\n' "$out" >&2; exit 1; }

# ---------------------------------------------------------------------------
# Test 6: migration window — viz scripts are scanned under the canonical
# 02_analysis/stages/ as well as the legacy 02_analysis/scripts/, including
# when both dirs coexist mid-rename.
# ---------------------------------------------------------------------------
STG="$TMPDIR_TEST/stg"
mkdir -p "$STG/02_analysis/stages"
cat > "$STG/02_analysis/stages/10_qc_viz.R" <<'R'
library(ggplot2)
p <- ggplot(df) + theme(text = element_text(size = 8))
ggsave("x.png", p, width = 7, height = 5)
R
set +e
out=$("$SCIAGENT" validate --check figure-style --project-dir "$STG" 2>&1); rc=$?
set -e
printf '%s\n' "$out" | grep -q 'WARN figure-style:.*02_analysis/stages/10_qc_viz.R' \
    || { echo "FAIL [$_TEST_NAME] test6: stages/ viz script not scanned"; printf '%s\n' "$out" >&2; exit 1; }

# Both dirs present: findings from each are reported.
mkdir -p "$STG/02_analysis/scripts"
cp "$STG/02_analysis/stages/10_qc_viz.R" "$STG/02_analysis/scripts/09_old_viz.R"
set +e
out=$("$SCIAGENT" validate --check figure-style --project-dir "$STG" 2>&1); rc=$?
set -e
printf '%s\n' "$out" | grep -q '02_analysis/scripts/09_old_viz.R' \
    || { echo "FAIL [$_TEST_NAME] test6: legacy scripts/ dir no longer scanned"; printf '%s\n' "$out" >&2; exit 1; }
printf '%s\n' "$out" | grep -q '02_analysis/stages/10_qc_viz.R' \
    || { echo "FAIL [$_TEST_NAME] test6: stages/ dir dropped when both present"; printf '%s\n' "$out" >&2; exit 1; }

pass
