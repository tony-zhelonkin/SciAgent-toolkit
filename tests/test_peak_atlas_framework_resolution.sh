#!/usr/bin/env bash
# tests/test_peak_atlas_framework_resolution.sh
# peak-atlas-multiome sources its sibling skill's primitives — the default path
# must be SCRIPT-relative, and a miss must be fatal.
#
# Skills are mounted as sibling directories, so from
# <skills>/peak-atlas-multiome/scripts/ the framework is always at
# ../../peak-atlas-framework/scripts. R resolves a relative path against getwd(),
# though, so the old cwd-relative default was correct ONLY when the user happened
# to have cd'd into this skill's own scripts/ dir. And a miss merely warn()ed, so
# the run died much later on "could not find function clusterGRanges" — an error
# that points nowhere near the cause.
#
#   (a) STATIC: the default is derived from the script's own location (no bare
#       cwd-relative literal), the env override still exists, and a miss stop()s.
#   (b) DYNAMIC (skipped when Rscript is absent): the real resolution block,
#       placed in a faithful sibling-skill layout, resolves the framework from a
#       cwd that is NOT the script's directory — under both `Rscript file.R` and
#       `source(file.R)` — and ERRORS when the framework is missing.
set -u
. "$(dirname "$0")/_lib.sh"

SRC="$TOOLKIT_ROOT/skills/peak-atlas-multiome/scripts/call_peaks_multistrategy.R"
assert_file_exists "$SRC" "call_peaks_multistrategy.R missing"

# ---------------------------------------------------------------------------
# (a) static
# ---------------------------------------------------------------------------
if grep -qE 'unset[[:space:]]*=[[:space:]]*"\.\./\.\./peak-atlas-framework' "$SRC"; then
    echo "FAIL [$_TEST_NAME] framework default is a bare cwd-relative literal" >&2
    echo "  R resolves it against getwd(), so it only works if you cd into scripts/." >&2
    exit 1
fi
assert_grep 'PEAK_ATLAS_FRAMEWORK_SCRIPTS' "$SRC" "env override must remain supported"
assert_grep 'commandArgs' "$SRC" "script-relative default must handle the Rscript (--file=) case"
assert_grep 'ofile'       "$SRC" "script-relative default must handle the source() case"
if grep -qE 'warning\("Framework script not found' "$SRC"; then
    echo "FAIL [$_TEST_NAME] a missing framework script still only warn()s" >&2
    exit 1
fi
assert_grep 'stop("peak-atlas-framework script not found' "$SRC" \
    "a missing framework script must be fatal"

# Every framework file it sources must actually exist in the sibling skill.
FRAMEWORK="$TOOLKIT_ROOT/skills/peak-atlas-framework/scripts"
for f in iterative_overlap.R support_voting.R normalize_width.R blacklist.R; do
    assert_file_exists "$FRAMEWORK/$f" "sibling skill is missing $f"
done

# No OTHER peak-atlas script may reintroduce the cwd-relative idiom.
if grep -rnE '(source|sys\.source)\([[:space:]]*"\.\./' "$TOOLKIT_ROOT"/skills/peak-atlas-*/scripts/ 2>/dev/null; then
    echo "FAIL [$_TEST_NAME] a peak-atlas script sources a cwd-relative '../' path" >&2
    exit 1
fi

if ! command -v Rscript >/dev/null 2>&1; then
    echo "SKIP: Rscript absent — framework resolution not executed"
    pass
    exit 0
fi

# ---------------------------------------------------------------------------
# (b) dynamic — faithful sibling layout, foreign cwd
# ---------------------------------------------------------------------------
setup_tmpdir
SK="$TMPDIR_TEST/skills"
mkdir -p "$SK/peak-atlas-multiome/scripts" "$SK/peak-atlas-framework/scripts"

# The REAL resolution block, lifted verbatim out of the shipped script.
awk '/^# ---- source the framework primitives/ { on = 1 }
     on && /^# =====/ { exit }
     on { print }' "$SRC" > "$SK/peak-atlas-multiome/scripts/call_peaks_multistrategy.R"
cat >> "$SK/peak-atlas-multiome/scripts/call_peaks_multistrategy.R" <<'EOF'
# Printed IFF control flow got past the sourcing block. Case (4) below asserts
# this marker is ABSENT, which is what separates a fatal stop() from a
# warn-and-continue: a downstream failure is NOT evidence the miss was fatal.
cat("SOURCING_SURVIVED\n")
stopifnot(exists("clusterGRanges", mode = "function"),
          exists("calculate_strategy_support", mode = "function"),
          exists("normalize_to_501bp", mode = "function"),
          exists("load_blacklist", mode = "function"))
cat("FRAMEWORK_RESOLVED_OK\n")
EOF
if ! grep -q 'FRAMEWORK_SCRIPTS' "$SK/peak-atlas-multiome/scripts/call_peaks_multistrategy.R"; then
    echo "FAIL [$_TEST_NAME] could not extract the resolution block from the real script" >&2
    exit 1
fi

# Stub framework primitives, one per real file, with the real function names.
cat > "$SK/peak-atlas-framework/scripts/iterative_overlap.R" <<'EOF'
clusterGRanges <- function(...) NULL
convergeClusterGRanges <- function(...) NULL
EOF
cat > "$SK/peak-atlas-framework/scripts/support_voting.R" <<'EOF'
calculate_strategy_support <- function(...) NULL
add_adjusted_score <- function(...) NULL
EOF
echo 'normalize_to_501bp <- function(...) NULL' > "$SK/peak-atlas-framework/scripts/normalize_width.R"
cat > "$SK/peak-atlas-framework/scripts/blacklist.R" <<'EOF'
load_blacklist <- function(...) NULL
remove_blacklist_peaks <- function(...) NULL
EOF

TARGET="$SK/peak-atlas-multiome/scripts/call_peaks_multistrategy.R"
FOREIGN="$TMPDIR_TEST/foreign"; mkdir -p "$FOREIGN"   # NOT the script's dir

run_from_foreign() {   # $1 = how to invoke
    ( cd "$FOREIGN" && env -u PEAK_ATLAS_FRAMEWORK_SCRIPTS timeout 60 bash -c "$1" 2>&1 )
}

# 1. Rscript, foreign cwd.
set +e
OUT=$(run_from_foreign "Rscript '$TARGET'"); RC=$?
set -e
if [[ "$RC" -ne 0 ]] || ! printf '%s\n' "$OUT" | grep -q 'FRAMEWORK_RESOLVED_OK'; then
    echo "FAIL [$_TEST_NAME] Rscript from a foreign cwd did not resolve the framework (rc=$RC):" >&2
    printf '%s\n' "$OUT" >&2
    exit 1
fi

# 2. source()d, foreign cwd.
set +e
OUT=$(run_from_foreign "Rscript -e 'source(\"$TARGET\")'"); RC=$?
set -e
if [[ "$RC" -ne 0 ]] || ! printf '%s\n' "$OUT" | grep -q 'FRAMEWORK_RESOLVED_OK'; then
    echo "FAIL [$_TEST_NAME] source() from a foreign cwd did not resolve the framework (rc=$RC):" >&2
    printf '%s\n' "$OUT" >&2
    exit 1
fi

# 3. the env override still wins.
ALT="$TMPDIR_TEST/alt"; mkdir -p "$ALT"
cp "$SK/peak-atlas-framework/scripts/"*.R "$ALT/"
set +e
OUT=$( cd "$FOREIGN" && PEAK_ATLAS_FRAMEWORK_SCRIPTS="$ALT" timeout 60 Rscript "$TARGET" 2>&1 ); RC=$?
set -e
if [[ "$RC" -ne 0 ]] || ! printf '%s\n' "$OUT" | grep -q 'FRAMEWORK_RESOLVED_OK'; then
    echo "FAIL [$_TEST_NAME] PEAK_ATLAS_FRAMEWORK_SCRIPTS override stopped working (rc=$RC):" >&2
    printf '%s\n' "$OUT" >&2
    exit 1
fi

# 4. a miss is FATAL, and the message points at the fix.
rm "$SK/peak-atlas-framework/scripts/support_voting.R"
set +e
OUT=$(run_from_foreign "Rscript '$TARGET'"); RC=$?
set -e
if [[ "$RC" -eq 0 ]]; then
    echo "FAIL [$_TEST_NAME] a missing framework script did not abort the run:" >&2
    printf '%s\n' "$OUT" >&2
    exit 1
fi
if printf '%s\n' "$OUT" | grep -q 'SOURCING_SURVIVED'; then
    echo "FAIL [$_TEST_NAME] a missing framework script only warned — execution continued past" >&2
    echo "  the sourcing block and died later on an unrelated error, which is exactly the" >&2
    echo "  diagnostic disconnect this fix removes." >&2
    printf '%s\n' "$OUT" >&2
    exit 1
fi
printf '%s\n' "$OUT" | grep -q 'PEAK_ATLAS_FRAMEWORK_SCRIPTS' || {
    echo "FAIL [$_TEST_NAME] the fatal message does not name the env override:" >&2
    printf '%s\n' "$OUT" >&2
    exit 1
}
echo "  dynamic: resolved from a foreign cwd under Rscript + source(); miss is fatal"

pass
