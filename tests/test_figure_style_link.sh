#!/usr/bin/env bash
# tests/test_figure_style_link.sh
# P06: figure-style helper-lib symlink via `sciagent activate`.
#
# Covers:
#   (a) analysis-type project (has 02_analysis/): activate creates
#       02_analysis/helpers/figure-style -> toolkit lib; symlink is recorded
#       in manifest; files are readable through it.
#   (b) deactivate removes the symlink (teardown works).
#   (c) non-analysis project (no 02_analysis/): activate succeeds (rc 0)
#       and does NOT create 02_analysis/helpers/figure-style.
set -u
. "$(dirname "$0")/_lib.sh"

setup_tmpdir
FAKE="$TMPDIR_TEST/fake-toolkit"
build_fake_toolkit "$FAKE"

# Give the fake toolkit a real lib/figure-style/ dir with stub helper files,
# so there is a valid target for the symlink.
mkdir -p "$FAKE/lib/figure-style"
echo "# stub R helper" > "$FAKE/lib/figure-style/figure_helpers.R"
echo "# stub py helper" > "$FAKE/lib/figure-style/figure_helpers.py"

export SCIAGENT_TOOLKIT="$FAKE"
SCIAGENT="$FAKE/bin/sciagent"

# ---------------------------------------------------------------------------
# (a) analysis-type project: symlink created + recorded in manifest
# ---------------------------------------------------------------------------
mkdir -p "$TMPDIR_TEST/ana_proj/02_analysis"
cd "$TMPDIR_TEST/ana_proj"

"$SCIAGENT" activate base >/dev/null

# Symlink must exist at the expected path.
assert_symlink "02_analysis/helpers/figure-style" \
    "figure-style symlink not created in 02_analysis/helpers/"

# Files must be readable through the symlink (resolves to the lib).
assert_file_exists "02_analysis/helpers/figure-style/figure_helpers.R" \
    "figure_helpers.R not readable through symlink"
assert_file_exists "02_analysis/helpers/figure-style/figure_helpers.py" \
    "figure_helpers.py not readable through symlink"

# Symlink path must appear in the manifest so teardown tracks it.
assert_grep '02_analysis/helpers/figure-style' .sciagent/manifest.json \
    "symlink not recorded in manifest.json"

# ---------------------------------------------------------------------------
# (b) deactivate removes the symlink
# ---------------------------------------------------------------------------
"$SCIAGENT" deactivate >/dev/null

if [[ -L "02_analysis/helpers/figure-style" || -e "02_analysis/helpers/figure-style" ]]; then
    echo "FAIL [$_TEST_NAME] figure-style symlink survived deactivate" >&2
    exit 1
fi

# ---------------------------------------------------------------------------
# (c) non-analysis project (no 02_analysis/): activate succeeds, no symlink
# ---------------------------------------------------------------------------
mkdir -p "$TMPDIR_TEST/sw_proj"
cd "$TMPDIR_TEST/sw_proj"

"$SCIAGENT" activate base >/dev/null   # must exit 0

if [[ -L "02_analysis/helpers/figure-style" || -e "02_analysis/helpers/figure-style" ]]; then
    echo "FAIL [$_TEST_NAME] figure-style symlink created in non-analysis project" >&2
    exit 1
fi
if [[ -d "02_analysis" ]]; then
    echo "FAIL [$_TEST_NAME] 02_analysis/ created in non-analysis project" >&2
    exit 1
fi

pass
