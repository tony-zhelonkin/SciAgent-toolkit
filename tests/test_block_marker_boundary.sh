#!/usr/bin/env bash
# tests/test_block_marker_boundary.sh — enforces the 5.3 invariant that
# block.sh is the sole authority for managed-block framing + hashing. The
# marker constants and the sha1 canonicalisation must not appear in any other
# library module; everything routes through block.sh accessors.
set -u
. "$(dirname "$0")/_lib.sh"

LIB="$TOOLKIT_ROOT/lib/scio"

# Patterns that may only appear in block.sh.
# NOTE: _BLOCK_BEGIN_PREFIX / _BLOCK_END were replaced by _block_begin_prefix /
# _block_end_marker helpers (parameterised by block-id). We protect the new
# helper names and the HTML-comment format substring '<!-- BEGIN SCIO:'
# (which subsumes the old literal 'BEGIN SCIO:ROLES' and prevents any new
# caller from hard-coding the HTML managed-block marker format). The gitignore
# module uses a separate '# BEGIN SCIO:GITIGNORE' comment-style marker
# that is intentionally NOT covered by block.sh and is not matched by the
# more specific '<!-- BEGIN SCIO:' pattern. sha1sum is still guarded.
patterns=(
    '_block_begin_prefix'
    '_block_end_marker'
    '<!-- BEGIN SCIO:'
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
