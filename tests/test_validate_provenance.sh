#!/usr/bin/env bash
# tests/test_validate_provenance.sh — opt-in `--check provenance` (+ captions +
# freshness, which share this fixture file).
#
# Tests:
#   1. CONFORMANT git project (committed script cited by every caption;
#      path-qualified captions present) → provenance + captions exit 0 default
#      AND under --strict.
#   2. provenance VIOLATION: caption cites a missing script → default exit 0
#      WARN; --strict exit 1. And an untracked (uncommitted) script in a git
#      repo is also a finding.
#   3. captions VIOLATION: artifact with no path-qualified `## <rel>` heading
#      in the sibling stage README → default exit 0 WARN; --strict exit 1.
#   4. freshness: a CRAFT block stale vs craft.yaml content (hash mismatch)
#      warns (exit 0) and hard-fails under --strict; an up-to-date block is clean.
#   5. No-false-positive: `--check all` on the conformant project → exit 0.
set -u
. "$(dirname "$0")/_lib.sh"

setup_tmpdir
FAKE="$TMPDIR_TEST/fake-toolkit"
build_fake_toolkit "$FAKE"
export SCIO_TOOLKIT="$FAKE"
SCIO="$FAKE/bin/scio"

# The fake toolkit has no craft.yaml; give it one (version 1) so freshness has
# a toolkit baseline to compare against.
printf 'version: 1\nfloors: {}\nbody: |\n  craft\n' > "$FAKE/craft.yaml"

# ---------------------------------------------------------------------------
# Build a CONFORMANT git project.
# ---------------------------------------------------------------------------
CONF="$TMPDIR_TEST/conf"
mkdir -p "$CONF/02_analysis/config" "$CONF/02_analysis/scripts"
mkdir -p "$CONF/03_results/01_qc/figures/_overview" "$CONF/03_results/01_qc/tables/_overview"
cat > "$CONF/02_analysis/config/analysis_config.yaml" <<'YAML'
figures:
  base_size: 16
  overview_dir: "_overview"
stages:
  - id: "01_qc"
paths:
  results: "03_results/"
YAML
echo "print('viz')" > "$CONF/02_analysis/scripts/10_qc_viz.py"
# Compute sibling: doc 09 §3.3 rule 3 (orphan viz) requires a same-number,
# same-stem compute stage. Without it this fixture is itself non-conformant.
echo "print('qc')" > "$CONF/02_analysis/scripts/10_qc.py"
touch "$CONF/03_results/01_qc/figures/_overview/qc.screen.png"
touch "$CONF/03_results/01_qc/tables/_overview/qc.csv"
cat > "$CONF/03_results/01_qc/README.md" <<'MD'
# 01_qc — artifact captions

## figures/_overview/qc.screen.png

A figure.

**How to read:** glyphs.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/10_qc_viz.py` | `save_overview` | `figures.top_n = 20` | `03_results/objects/qc.h5ad` |

## tables/_overview/qc.csv

Source table.

**How to read:** raw.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/10_qc_viz.py` | `save_overview` | `figures.top_n = 20` | `03_results/objects/qc.h5ad` |
MD
# Render a REAL CRAFT block from the fake craft.yaml so its stored hash matches
# what freshness recomputes — the conformant case is up to date.
printf '# CONF project\n' > "$CONF/AGENTS.md"
( . "$FAKE/lib/scio/block.sh"; . "$FAKE/lib/scio/craft.sh"
  SCIO_TOOLKIT="$FAKE" craft_render_and_write "$CONF/AGENTS.md" )
git -C "$CONF" init -q
git -C "$CONF" config user.email t@e.com
git -C "$CONF" config user.name T
git -C "$CONF" add -A
git -C "$CONF" commit -qm init

# ---------------------------------------------------------------------------
# Test 1: conformant — provenance + captions pass default + strict.
# ---------------------------------------------------------------------------
for chk in provenance captions; do
    set +e
    out=$("$SCIO" lint --check "$chk" --project-dir "$CONF" 2>&1); rc=$?
    set -e
    [[ "$rc" -eq 0 ]] || { echo "FAIL [$_TEST_NAME] test1: conformant $chk default exit $rc"; printf '%s\n' "$out" >&2; exit 1; }
    set +e
    out=$("$SCIO" lint --check "$chk" --strict --project-dir "$CONF" 2>&1); rc=$?
    set -e
    [[ "$rc" -eq 0 ]] || { echo "FAIL [$_TEST_NAME] test1: conformant $chk strict exit $rc"; printf '%s\n' "$out" >&2; exit 1; }
done

# ---------------------------------------------------------------------------
# Test 2: provenance violation — caption cites a missing script.
# ---------------------------------------------------------------------------
PV="$TMPDIR_TEST/pv"
mkdir -p "$PV/02_analysis/config" "$PV/02_analysis/scripts"
mkdir -p "$PV/03_results/01_qc/figures/_overview" "$PV/03_results/01_qc/tables/_overview"
cat > "$PV/02_analysis/config/analysis_config.yaml" <<'YAML'
stages:
  - id: "01_qc"
YAML
touch "$PV/03_results/01_qc/figures/_overview/qc.screen.png"
touch "$PV/03_results/01_qc/tables/_overview/qc.csv"
cat > "$PV/03_results/01_qc/README.md" <<'MD'
# 01_qc

## figures/_overview/qc.screen.png

F.

**How to read:** x.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/does_not_exist.py` | `f` | `c` | `i` |

## tables/_overview/qc.csv

T.

**How to read:** x.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/does_not_exist.py` | `f` | `c` | `i` |
MD

set +e
out=$("$SCIO" lint --check provenance --project-dir "$PV" 2>&1); rc=$?
set -e
[[ "$rc" -eq 0 ]] || { echo "FAIL [$_TEST_NAME] test2: provenance default expected 0 got $rc"; printf '%s\n' "$out" >&2; exit 1; }
printf '%s\n' "$out" | grep -q 'WARN provenance:.*missing script: 02_analysis/scripts/does_not_exist.py' \
    || { echo "FAIL [$_TEST_NAME] test2: expected missing-script WARN"; printf '%s\n' "$out" >&2; exit 1; }

set +e
out=$("$SCIO" lint --check provenance --strict --project-dir "$PV" 2>&1); rc=$?
set -e
[[ "$rc" -eq 1 ]] || { echo "FAIL [$_TEST_NAME] test2: provenance strict expected 1 got $rc"; printf '%s\n' "$out" >&2; exit 1; }

# Untracked (uncommitted) script in a git repo is a provenance finding too.
PVG="$TMPDIR_TEST/pvg"
mkdir -p "$PVG/02_analysis/config" "$PVG/02_analysis/scripts"
mkdir -p "$PVG/03_results/01_qc/figures/_overview" "$PVG/03_results/01_qc/tables/_overview"
cat > "$PVG/02_analysis/config/analysis_config.yaml" <<'YAML'
stages:
  - id: "01_qc"
YAML
echo "print('viz')" > "$PVG/02_analysis/scripts/10_qc_viz.py"
# Compute sibling: doc 09 §3.3 rule 3 (orphan viz) requires a same-number,
# same-stem compute stage. Without it this fixture is itself non-conformant.
echo "print('qc')" > "$PVG/02_analysis/scripts/10_qc.py"
touch "$PVG/03_results/01_qc/figures/_overview/qc.screen.png"
touch "$PVG/03_results/01_qc/tables/_overview/qc.csv"
cat > "$PVG/03_results/01_qc/README.md" <<'MD'
# 01_qc

## figures/_overview/qc.screen.png

F.

**How to read:** x.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/10_qc_viz.py` | `f` | `c` | `i` |

## tables/_overview/qc.csv

T.

**How to read:** x.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/10_qc_viz.py` | `f` | `c` | `i` |
MD
git -C "$PVG" init -q
git -C "$PVG" config user.email t@e.com
git -C "$PVG" config user.name T
# Intentionally DO NOT commit the script (it stays untracked).
set +e
out=$("$SCIO" lint --check provenance --project-dir "$PVG" 2>&1); rc=$?
set -e
[[ "$rc" -eq 0 ]] || { echo "FAIL [$_TEST_NAME] test2b: untracked default expected 0 got $rc"; printf '%s\n' "$out" >&2; exit 1; }
printf '%s\n' "$out" | grep -q 'WARN provenance:.*untracked script' \
    || { echo "FAIL [$_TEST_NAME] test2b: expected untracked-script WARN"; printf '%s\n' "$out" >&2; exit 1; }
# Once committed, the finding clears.
git -C "$PVG" add -A
git -C "$PVG" commit -qm init
set +e
out=$("$SCIO" lint --check provenance --strict --project-dir "$PVG" 2>&1); rc=$?
set -e
[[ "$rc" -eq 0 ]] || { echo "FAIL [$_TEST_NAME] test2b: committed strict expected 0 got $rc"; printf '%s\n' "$out" >&2; exit 1; }

# ---------------------------------------------------------------------------
# Test 3: captions violation — artifact with no caption heading.
# ---------------------------------------------------------------------------
CV="$TMPDIR_TEST/cv"
mkdir -p "$CV/02_analysis/config"
mkdir -p "$CV/03_results/01_qc/figures/_overview" "$CV/03_results/01_qc/tables/_overview"
cat > "$CV/02_analysis/config/analysis_config.yaml" <<'YAML'
stages:
  - id: "01_qc"
YAML
touch "$CV/03_results/01_qc/figures/_overview/uncaptioned.screen.png"
touch "$CV/03_results/01_qc/tables/_overview/uncaptioned.csv"
cat > "$CV/03_results/01_qc/README.md" <<'MD'
# 01_qc

## figures/_overview/something_else.screen.png

Wrong heading — does not match the artifact.
MD

set +e
out=$("$SCIO" lint --check captions --project-dir "$CV" 2>&1); rc=$?
set -e
[[ "$rc" -eq 0 ]] || { echo "FAIL [$_TEST_NAME] test3: captions default expected 0 got $rc"; printf '%s\n' "$out" >&2; exit 1; }
printf '%s\n' "$out" | grep -q 'WARN captions:.*no caption section.*uncaptioned' \
    || { echo "FAIL [$_TEST_NAME] test3: expected uncaptioned WARN"; printf '%s\n' "$out" >&2; exit 1; }

set +e
out=$("$SCIO" lint --check captions --strict --project-dir "$CV" 2>&1); rc=$?
set -e
[[ "$rc" -eq 1 ]] || { echo "FAIL [$_TEST_NAME] test3: captions strict expected 1 got $rc"; printf '%s\n' "$out" >&2; exit 1; }

# ---------------------------------------------------------------------------
# Test 4: freshness — CRAFT block stale vs current craft.yaml content (hash).
# ---------------------------------------------------------------------------
# Change the toolkit craft.yaml body so the conformant project's rendered CRAFT
# block hash no longer matches what a fresh render would produce.
printf 'version: 1\nfloors: {}\nbody: |\n  craft CHANGED\n' > "$FAKE/craft.yaml"
set +e
out=$("$SCIO" lint --check freshness --project-dir "$CONF" 2>&1); rc=$?
set -e
[[ "$rc" -eq 0 ]] || { echo "FAIL [$_TEST_NAME] test4: freshness default expected 0 got $rc"; printf '%s\n' "$out" >&2; exit 1; }
printf '%s\n' "$out" | grep -q 'WARN freshness:.*CRAFT block is stale.*scio craft' \
    || { echo "FAIL [$_TEST_NAME] test4: expected stale-CRAFT WARN with craft hint"; printf '%s\n' "$out" >&2; exit 1; }

set +e
out=$("$SCIO" lint --check freshness --strict --project-dir "$CONF" 2>&1); rc=$?
set -e
[[ "$rc" -eq 1 ]] || { echo "FAIL [$_TEST_NAME] test4: freshness strict expected 1 got $rc"; printf '%s\n' "$out" >&2; exit 1; }

# Restore craft.yaml body → conformant block hash matches again → clean.
printf 'version: 1\nfloors: {}\nbody: |\n  craft\n' > "$FAKE/craft.yaml"
set +e
out=$("$SCIO" lint --check freshness --strict --project-dir "$CONF" 2>&1); rc=$?
set -e
[[ "$rc" -eq 0 ]] || { echo "FAIL [$_TEST_NAME] test4: up-to-date freshness strict expected 0 got $rc"; printf '%s\n' "$out" >&2; exit 1; }

# ---------------------------------------------------------------------------
# Test 5: no-false-positive — `--check all` on conformant exits 0.
# ---------------------------------------------------------------------------
set +e
out=$("$SCIO" lint --check all --project-dir "$CONF" 2>&1); rc=$?
set -e
[[ "$rc" -eq 0 ]] || { echo "FAIL [$_TEST_NAME] test5: --check all conformant default exit $rc"; printf '%s\n' "$out" >&2; exit 1; }
set +e
out=$("$SCIO" lint --check all --strict --project-dir "$CONF" 2>&1); rc=$?
set -e
[[ "$rc" -eq 0 ]] || { echo "FAIL [$_TEST_NAME] test5: --check all conformant strict exit $rc"; printf '%s\n' "$out" >&2; exit 1; }

# ---------------------------------------------------------------------------
# Test 6: migration window — 02_analysis/stages/ (canonical) is accepted just
# like the legacy 02_analysis/scripts/, and a path under neither still warns.
# ---------------------------------------------------------------------------
ST="$TMPDIR_TEST/stages"
mkdir -p "$ST/02_analysis/config" "$ST/02_analysis/stages" "$ST/02_analysis/helpers"
mkdir -p "$ST/03_results/01_qc/figures/_overview" "$ST/03_results/01_qc/tables/_overview"
cat > "$ST/02_analysis/config/analysis_config.yaml" <<'YAML'
stages:
  - id: "01_qc"
YAML
echo "1" > "$ST/02_analysis/stages/10_qc_viz.R"
# Compute sibling: doc 09 §3.3 rule 3 (orphan viz) requires a same-number,
# same-stem compute stage. Without it this fixture is itself non-conformant.
echo "1" > "$ST/02_analysis/stages/10_qc.R"
echo "1" > "$ST/02_analysis/helpers/plot_utils.R"
touch "$ST/03_results/01_qc/figures/_overview/qc.screen.png"
touch "$ST/03_results/01_qc/tables/_overview/qc.csv"
cat > "$ST/03_results/01_qc/README.md" <<'MD'
# 01_qc

## figures/_overview/qc.screen.png

F.

**How to read:** x.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/stages/10_qc_viz.R` | `f` | `c` | `i` |

## tables/_overview/qc.csv

T.

**How to read:** x.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/stages/10_qc_viz.R` | `f` | `c` | `i` |
MD

set +e
out=$("$SCIO" lint --check provenance --strict --project-dir "$ST" 2>&1); rc=$?
set -e
[[ "$rc" -eq 0 ]] || { echo "FAIL [$_TEST_NAME] test6: stages/ provenance strict expected 0 got $rc"; printf '%s\n' "$out" >&2; exit 1; }

# A script outside both accepted dirs is still a finding, hinted at stages/.
sed -i 's#02_analysis/stages/10_qc_viz.R#02_analysis/helpers/plot_utils.R#g' \
    "$ST/03_results/01_qc/README.md"
set +e
out=$("$SCIO" lint --check provenance --project-dir "$ST" 2>&1); rc=$?
set -e
[[ "$rc" -eq 0 ]] || { echo "FAIL [$_TEST_NAME] test6: off-stage default expected 0 got $rc"; printf '%s\n' "$out" >&2; exit 1; }
printf '%s\n' "$out" | grep -q 'WARN provenance:.*not under 02_analysis/stages/' \
    || { echo "FAIL [$_TEST_NAME] test6: expected stages/-worded off-stage WARN"; printf '%s\n' "$out" >&2; exit 1; }

# ---------------------------------------------------------------------------
# Test 7: the new `scio lint --check` surface covers the same ground.
# ---------------------------------------------------------------------------
set +e
out=$("$SCIO" lint --check provenance --project-dir "$ST" 2>&1); rc=$?
set -e
[[ "$rc" -eq 0 ]] || { echo "FAIL [$_TEST_NAME] test7: lint --check provenance default expected 0 got $rc"; printf '%s\n' "$out" >&2; exit 1; }
printf '%s\n' "$out" | grep -q 'WARN provenance:.*not under 02_analysis/stages/' \
    || { echo "FAIL [$_TEST_NAME] test7: expected off-stage WARN via lint --check"; printf '%s\n' "$out" >&2; exit 1; }

pass
