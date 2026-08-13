#!/usr/bin/env bash
# tests/test_validate_frontmatter_shape.sh — the toolkit-wide walk in
# `scio lint --check toolkit` hard-fails on malformed skill frontmatter.
#
# This is the check that replaced the requires-graph and tag-vocab walks when
# the metadata: block was retired. `name` and `description` are the entire
# surface a harness preloads, so they are the two fields worth failing on:
#   1. clean fixture toolkit passes
#   2. missing name:            -> exit 1, message names the skill
#   3. name: != directory       -> exit 1
#   4. missing description:     -> exit 1
#   5. description over the cap -> exit 1, message reports the length
#   6. description at the cap   -> exit 0 (boundary is inclusive)
set -u
. "$(dirname "$0")/_lib.sh"

setup_tmpdir

# _plant <frontmatter-body> — rewrite skills/s_a in a fresh fake toolkit.
# Echoes the toolkit root.
_plant() {
    local fm="$1"
    local root="$TMPDIR_TEST/tk_$RANDOM$RANDOM"
    build_fake_toolkit "$root"
    {
        printf -- '---\n'
        printf '%s\n' "$fm"
        printf -- '---\n\nbody\n'
    } > "$root/skills/s_a/SKILL.md"
    printf '%s\n' "$root"
}

# _run <toolkit-root> — run the catalog check and capture output and status.
_run() {
    set +e
    out=$(SCIO_TOOLKIT="$1" "$1/bin/scio" lint --check toolkit 2>&1)
    rc=$?
    set -e
}

# --- 1. clean fixture passes -------------------------------------------------
CLEAN=$(_plant "name: s_a
description: A perfectly ordinary fixture skill.")
_run "$CLEAN"
if [[ "$rc" -ne 0 ]]; then
    echo "FAIL [$_TEST_NAME] case1: clean fixture failed catalog lint (rc=$rc)" >&2
    printf '%s\n' "$out" >&2
    exit 1
fi

# --- 2. missing name: --------------------------------------------------------
NONAME=$(_plant "description: A fixture with no name field.")
_run "$NONAME"
if [[ "$rc" -eq 0 ]]; then
    echo "FAIL [$_TEST_NAME] case2: missing name passed catalog lint" >&2
    exit 1
fi
if ! printf '%s\n' "$out" | grep -q 's_a.*no name:'; then
    echo "FAIL [$_TEST_NAME] case2: message does not name the offending skill" >&2
    printf '%s\n' "$out" >&2
    exit 1
fi

# --- 3. name: does not match the directory -----------------------------------
MISMATCH=$(_plant "name: not_s_a
description: A fixture whose name disagrees with its directory.")
_run "$MISMATCH"
if [[ "$rc" -eq 0 ]]; then
    echo "FAIL [$_TEST_NAME] case3: mismatched name passed catalog lint" >&2
    exit 1
fi
if ! printf '%s\n' "$out" | grep -q "does not match its directory"; then
    echo "FAIL [$_TEST_NAME] case3: message does not explain the mismatch" >&2
    printf '%s\n' "$out" >&2
    exit 1
fi

# --- 4. missing description: -------------------------------------------------
NODESC=$(_plant "name: s_a")
_run "$NODESC"
if [[ "$rc" -eq 0 ]]; then
    echo "FAIL [$_TEST_NAME] case4: missing description passed catalog lint" >&2
    exit 1
fi
if ! printf '%s\n' "$out" | grep -q 's_a.*no description:'; then
    echo "FAIL [$_TEST_NAME] case4: message does not name the offending skill" >&2
    printf '%s\n' "$out" >&2
    exit 1
fi

# --- 5. description over the cap ---------------------------------------------
LONG=$(head -c 351 < /dev/zero | tr '\0' 'x')
OVER=$(_plant "name: s_a
description: $LONG")
_run "$OVER"
if [[ "$rc" -eq 0 ]]; then
    echo "FAIL [$_TEST_NAME] case5: 351-char description passed catalog lint" >&2
    exit 1
fi
if ! printf '%s\n' "$out" | grep -q 'description is 351 chars (max 350)'; then
    echo "FAIL [$_TEST_NAME] case5: message does not report the measured length" >&2
    printf '%s\n' "$out" >&2
    exit 1
fi

# --- 6. description exactly at the cap passes (inclusive boundary) -----------
EXACT=$(head -c 350 < /dev/zero | tr '\0' 'x')
AT=$(_plant "name: s_a
description: $EXACT")
_run "$AT"
if [[ "$rc" -ne 0 ]]; then
    echo "FAIL [$_TEST_NAME] case6: 350-char description rejected (cap must be inclusive)" >&2
    printf '%s\n' "$out" >&2
    exit 1
fi

pass
