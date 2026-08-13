#!/usr/bin/env bash
# tests/test_interactive_helpers_contract.sh — the interactive-style contract (lib/interactive-style/).
#
# Validates the interactive-explorer helper lib WITHOUT jscatter / pandas / numpy / matplotlib:
#   (a) LAZY-IMPORT DISCIPLINE: with jscatter FORCED ABSENT (a sitecustomize that raises
#       ImportError on `import jscatter`), python3 still imports interactive_helpers cleanly —
#       proving no heavy top-level import. This is the load-bearing guarantee: the module must
#       import on a bare analysis box that has no jscatter.
#   (b) NO-BACKEND EXECUTION: in a tmpdir with a minimal analysis_config.yaml, the stdlib+yaml
#       functions run — load_interactive_config reads the `interactive:` block, the config
#       accessors return the defaults backfilled with project values, and the OOM lifecycle
#       registry is sane on an empty session (live_panels() == [], close_panels() == 0).
#   (c) GRACEFUL SKIP: if pyyaml is absent, execution is not exercised and the test still passes.
set -u
. "$(dirname "$0")/_lib.sh"

LIB="$TOOLKIT_ROOT/lib/interactive-style"
PY_FILE="$LIB/interactive_helpers.py"

assert_file_exists "$PY_FILE" "interactive_helpers.py missing"

# ---------------------------------------------------------------------------
# Contract function names must all be present (static, no toolchain).
# ---------------------------------------------------------------------------
CONTRACT_NAMES=(
    find_root load_interactive_config load_explorer color_key
    grid close_panels live_panels first_selection save_selection snapshot
)
for name in "${CONTRACT_NAMES[@]}"; do
    if ! grep -Eq "^[[:space:]]*def[[:space:]]+${name}[[:space:]]*\(" "$PY_FILE"; then
        echo "FAIL [$_TEST_NAME] '$name' not defined in interactive_helpers.py" >&2
        exit 1
    fi
done
echo "  names: all ${#CONTRACT_NAMES[@]} contract functions present"

# Skip execution gracefully if pyyaml (the only hard dep) is absent.
if ! python3 -c "import yaml" 2>/dev/null; then
    echo "SKIP: pyyaml absent — interactive-style Python execution not exercised"
    pass
    exit 0
fi

setup_tmpdir
cp "$PY_FILE" "$TMPDIR_TEST/interactive_helpers.py"

# A fake sitecustomize that makes `import jscatter` raise ImportError, simulating a bare box.
# Put it FIRST on sys.path so the import machinery sees the shim before any real jscatter.
mkdir -p "$TMPDIR_TEST/fakemods"
cat > "$TMPDIR_TEST/fakemods/jscatter.py" <<'EOF'
raise ImportError("jscatter is FORCED ABSENT by the interactive-helpers contract test")
EOF

cat > "$TMPDIR_TEST/analysis_config.yaml" <<'EOF'
interactive:
  categorical_obs: [leiden, batch]
  grid_rows: 3
  save_selection_cols: [x, y, leiden]
paths:
  results: "03_results/"
  interactive: "03_results/interactive/"
EOF

set +e
PY_OUT=$(cd "$TMPDIR_TEST" && python3 - <<'PYEOF' 2>&1
import sys
sys.path.insert(0, "fakemods")   # `import jscatter` now raises ImportError
sys.path.insert(0, ".")

# (a) MUST import with jscatter absent (proves no heavy top-level import).
try:
    import jscatter  # noqa: F401
    print("FAIL: jscatter should have been forced absent"); sys.exit(1)
except ImportError:
    pass
import interactive_helpers as ih   # this line is the actual contract: import must NOT need jscatter

# (b) no-backend execution: config read + accessors + lifecycle registry.
cfg = ih.load_interactive_config("analysis_config.yaml")
assert ih._int_get(cfg, "grid_rows") == 3, "project value not read"
assert ih._int_get(cfg, "grid_height") == 340, "default not backfilled"
assert ih._int_get(cfg, "categorical_obs") == ["leiden", "batch"], "list value not read"
assert ih._int_get(cfg, "save_selection_cols") == ["x", "y", "leiden"], "cols not read"
# a config WITHOUT an interactive block falls back entirely to defaults
assert ih._int_get({}, "save_selection_cols") == ["x", "y"], "empty-config default wrong"

# OOM lifecycle registry sane on an empty session (no jscatter needed for these).
assert ih.live_panels() == [], "live_panels should start empty"
assert ih.close_panels() == 0, "close_panels on empty registry should free 0"

print("PY_CONTRACT_OK")
PYEOF
)
PY_RC=$?
set -e

if [[ "$PY_RC" -ne 0 ]] || ! printf '%s\n' "$PY_OUT" | grep -q 'PY_CONTRACT_OK'; then
    echo "FAIL [$_TEST_NAME] python lazy-import / no-backend execution failed:" >&2
    printf '%s\n' "$PY_OUT" >&2
    exit 1
fi
echo "  python: imports with jscatter ABSENT + config/accessor/lifecycle execution OK"

pass
