#!/usr/bin/env bash
# tests/test_block_marker_boundary.sh — enforces the 5.3 invariant that
# block.sh is the sole authority for managed-block framing + hashing. The
# marker constants and the sha1 canonicalisation must not appear in any other
# library module; everything routes through block.sh accessors.
set -u
. "$(dirname "$0")/_lib.sh"

LIB="$TOOLKIT_ROOT/lib/sciagent"

# Patterns that may only appear in block.sh.
patterns=(
    '_BLOCK_BEGIN_PREFIX'
    '_BLOCK_END'
    'BEGIN SCIAGENT:ROLES'
    'sha1sum'
)

violations=()
for f in "$LIB"/*.sh; do
    [[ "$(basename "$f")" == "block.sh" ]] && continue
    for p in "${patterns[@]}"; do
        while IFS= read -r match; do
            violations+=("$(basename "$f"): $match")
        done < <(grep -nF -- "$p" "$f" 2>/dev/null)
    done
done

if (( ${#violations[@]} > 0 )); then
    echo "FAIL [$_TEST_NAME] managed-block internals leaked outside block.sh:" >&2
    printf '  %s\n' "${violations[@]}" >&2
    exit 1
fi

pass
