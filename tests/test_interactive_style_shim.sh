#!/usr/bin/env bash
# tests/test_interactive_style_shim.sh — the per-project interactive-explorer shim.
#
# The mount `symlink_create_helper_lib` creates is HYPHENATED
# (02_analysis/helpers/interactive-style/), which is not a legal Python package
# name. Projects therefore import the UNDERSCORED shim next to it, exactly as
# they do for figure-style. This test pins that shim end to end:
#
#   (a) the template exists in the analysis project-template tree;
#   (b) `sciagent new project --type analysis` MATERIALIZES it — asserted at the
#       path derived from the template's own location under templates/, so the
#       test cannot drift from where _render_tree actually writes;
#   (c) with the lib symlink present, `from helpers.interactive_style import ...`
#       resolves and re-exports the whole interactive_helpers contract;
#   (d) with the lib symlink ABSENT the import FAILS LOUDLY (ImportError) rather
#       than degrading to silent stubs — a deliberate divergence from
#       figure_style.py, whose stubs are safe because a mis-styled figure is
#       still a figure, while a stubbed explorer persists nothing.
#
# (c)/(d) skip gracefully when pyyaml is absent (interactive_helpers' only hard dep).
set -u
. "$(dirname "$0")/_lib.sh"

TPL_ROOT="$TOOLKIT_ROOT/templates/project/analysis"
TPL="$TPL_ROOT/02_analysis/helpers/interactive_style.py.template"

# ---------------------------------------------------------------------------
# (a) the template exists
# ---------------------------------------------------------------------------
assert_file_exists "$TPL" \
    "interactive_style.py.template missing — the hyphenated interactive-style mount has no importable shim"

# The rendered path is DERIVED from the template's own location (strip the
# templates/project/analysis/ prefix and the .template suffix), so this test
# tracks _render_tree instead of duplicating its convention.
REL="${TPL#"$TPL_ROOT"/}"
REL="${REL%.template}"
assert_eq "$REL" "02_analysis/helpers/interactive_style.py" "derived render path"

# ---------------------------------------------------------------------------
# (b) `sciagent new project --type analysis` materializes it
# ---------------------------------------------------------------------------
setup_tmpdir
PROJ="$TMPDIR_TEST/ana"
"$TOOLKIT_ROOT/bin/sciagent" new project "$PROJ" --type analysis >/dev/null 2>&1 \
    || { echo "FAIL [$_TEST_NAME] scaffold --type analysis failed" >&2; exit 1; }

assert_file_exists "$PROJ/$REL" \
    "new project did not materialize $REL (template present but never rendered)"
assert_file_eq "$PROJ/$REL" "$TPL" "rendered shim differs from its template"

# ---------------------------------------------------------------------------
# (c)/(d) import behaviour — needs pyyaml (interactive_helpers' only hard dep)
# ---------------------------------------------------------------------------
if ! python3 -c "import yaml" 2>/dev/null; then
    echo "SKIP: pyyaml absent — shim import behaviour not exercised"
    pass
    exit 0
fi

# (d) FIRST, while the lib symlink does NOT yet exist: import must raise.
set +e
OUT=$(cd "$PROJ" && python3 - <<'PY' 2>&1
import sys
sys.path.insert(0, "02_analysis")
try:
    import helpers.interactive_style  # noqa: F401
except ImportError as exc:
    assert "interactive-style" in str(exc), f"unhelpful message: {exc}"
    assert "sciagent activate" in str(exc), f"message must name the remedy: {exc}"
    print("NO_LIB_RAISES_OK")
else:
    print("FAIL: shim imported with no toolkit lib — silent stubs are not acceptable here")
    sys.exit(1)
PY
)
RC=$?
set -e
if [[ "$RC" -ne 0 ]] || ! printf '%s\n' "$OUT" | grep -q 'NO_LIB_RAISES_OK'; then
    echo "FAIL [$_TEST_NAME] shim did not fail loudly with the lib absent:" >&2
    printf '%s\n' "$OUT" >&2
    exit 1
fi

# (c) now link the real lib the way `sciagent activate` does, and import for real.
mkdir -p "$PROJ/02_analysis/helpers"
ln -sfn "$TOOLKIT_ROOT/lib/interactive-style" "$PROJ/02_analysis/helpers/interactive-style"

set +e
OUT=$(cd "$PROJ" && python3 - <<'PY' 2>&1
import sys
sys.path.insert(0, "02_analysis")
# Every name the explorer .qmd and SKILL.md import off the shim.
from helpers.interactive_style import (  # noqa: F401
    find_root, load_interactive_config, load_explorer, color_key,
    grid, close_panels, live_panels, first_selection, save_selection, snapshot)
assert live_panels() == [] and close_panels() == 0, "lifecycle registry not reachable"
print("SHIM_IMPORT_OK")
PY
)
RC=$?
set -e
if [[ "$RC" -ne 0 ]] || ! printf '%s\n' "$OUT" | grep -q 'SHIM_IMPORT_OK'; then
    echo "FAIL [$_TEST_NAME] shim import failed with the lib linked:" >&2
    printf '%s\n' "$OUT" >&2
    exit 1
fi

pass
