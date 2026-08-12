#!/usr/bin/env bash
# tests/test_root_walk_guard.sh — a sentinel walk-up must terminate at "/".
#
# `Path("/").parent == Path("/")`, so a walk-up loop with no explicit
# fixed-point guard SPINS FOREVER when the sentinel is missing — the worst
# failure mode there is, because it looks like a slow notebook rather than a
# misplaced launch directory. The R twin
# (skills/decision-gate-notebook/assets/notebook.qmd) has always had the guard
# and a message naming the sentinel, the start dir, and the remedy; the Python
# side (skills/interactive-breakpoint-explorer/assets/explorer.qmd) did not.
#
#   (a) STATIC, repo-wide: every Python-side sentinel walk-up carries a
#       fixed-point guard. This is the part that catches a NEW copy of the
#       idiom, which is how the bug spread in the first place.
#   (b) DYNAMIC: the explorer .qmd's real setup chunk, executed verbatim from a
#       directory with no sentinel above it, EXITS with an error. Run under
#       `timeout`, and a timeout kill is an explicit failure — a hanging test is
#       worse than a failing one.
set -u
. "$(dirname "$0")/_lib.sh"

# ---------------------------------------------------------------------------
# (a) static scan: `while ... .exists()` + `d = d.parent` implies a guard.
# ---------------------------------------------------------------------------
scanned=0
walks=0
while IFS= read -r f; do
    scanned=$((scanned + 1))
    report=$(awk '
        # Enter a candidate walk when a while-condition tests .exists().
        /while[^#]*\.exists\(\)/ { inwalk = 1; wline = FNR; asc = 0; guard = 0; next }
        inwalk {
            if ($0 ~ /=[[:space:]]*[A-Za-z_][A-Za-z0-9_]*\.parent[[:space:]]*$/) asc = 1
            if ($0 ~ /\.parent[[:space:]]*==/ || $0 ~ /==[[:space:]]*[A-Za-z_][A-Za-z0-9_]*\.parent/) guard = 1
            # A blank line or a dedented statement ends the loop body.
            if ($0 !~ /^[[:space:]]/ || $0 ~ /^[[:space:]]*$/) {
                if (asc && !guard) print wline
                inwalk = 0
            }
        }
        END { if (inwalk && asc && !guard) print wline }
    ' "$f")
    if [[ -n "$report" ]]; then
        echo "FAIL [$_TEST_NAME] unguarded sentinel walk-up in ${f#"$TOOLKIT_ROOT"/} (line $report)" >&2
        echo "  Path('/').parent == Path('/') — without an explicit fixed-point check this" >&2
        echo "  loop never terminates when the sentinel is absent. Add:" >&2
        echo "      if d.parent == d: raise FileNotFoundError(...)" >&2
        exit 1
    fi
    if grep -qE 'while[^#]*\.exists\(\)' "$f"; then walks=$((walks + 1)); fi
done < <(grep -rlE 'while[^#]*\.exists\(\)' \
             --include='*.py' --include='*.qmd' --include='*.template' --include='*.ipynb' \
             "$TOOLKIT_ROOT" 2>/dev/null | grep -v '/\.git/' | grep -v '/__pycache__/')

if (( walks == 0 )); then
    echo "FAIL [$_TEST_NAME] scan found no walk-up loops at all — the grep is broken" >&2
    exit 1
fi
echo "  static: $walks Python walk-up loop(s) across $scanned file(s), all guarded"

# ---------------------------------------------------------------------------
# (b) dynamic: run the explorer .qmd's real setup chunk with no sentinel above.
# ---------------------------------------------------------------------------
QMD="$TOOLKIT_ROOT/skills/interactive-breakpoint-explorer/assets/explorer.qmd"
assert_file_exists "$QMD" "explorer.qmd missing"

setup_tmpdir
# Extract the FIRST python chunk verbatim, then cut at os.chdir so we exercise
# only the root resolution (not the helper import, which needs an activated project).
awk '/^```\{python\}/ { inchunk = 1; next } inchunk && /^```/ { exit } inchunk { print }' \
    "$QMD" | sed '/os\.chdir/,$d' > "$TMPDIR_TEST/walk.py"

if ! grep -q 'Path.cwd()' "$TMPDIR_TEST/walk.py"; then
    echo "FAIL [$_TEST_NAME] could not extract the setup chunk from explorer.qmd" >&2
    cat "$TMPDIR_TEST/walk.py" >&2
    exit 1
fi

# A directory guaranteed to have no analysis_config.yaml anywhere above it:
# mkdtemp under /tmp, and the walk ends at "/".
mkdir -p "$TMPDIR_TEST/nowhere/deep"
set +e
OUT=$(cd "$TMPDIR_TEST/nowhere/deep" && timeout 20 python3 "$TMPDIR_TEST/walk.py" 2>&1)
RC=$?
set -e

if [[ "$RC" -eq 124 ]]; then
    echo "FAIL [$_TEST_NAME] explorer.qmd root walk HUNG (killed at 20s) with the sentinel absent" >&2
    exit 1
fi
if [[ "$RC" -eq 0 ]]; then
    echo "FAIL [$_TEST_NAME] explorer.qmd root walk succeeded with no sentinel — impossible" >&2
    exit 1
fi
# Message quality: must name the sentinel, the start dir, and the remedy —
# matching the R twin in decision-gate-notebook/assets/notebook.qmd.
for needle in '02_analysis/config/analysis_config.yaml' "$TMPDIR_TEST/nowhere/deep" 'compartment'; do
    printf '%s\n' "$OUT" | grep -qF "$needle" || {
        echo "FAIL [$_TEST_NAME] error message does not mention '$needle':" >&2
        printf '%s\n' "$OUT" >&2
        exit 1
    }
done
echo "  dynamic: walk errors out (rc=$RC) instead of hanging, message names sentinel + start + remedy"

pass
