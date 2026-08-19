#!/usr/bin/env bash
# tests/test_craft_block_budget.sh — `scio lint --check toolkit` enforces the
# rendered CRAFT body against craft.yaml `max_lines:`.
#
# The CRAFT block is the toolkit's only always-on text: it is preloaded in every
# consumer session whether or not anything routes to it, so its budget is a cap
# the toolkit owes rather than advice. The cap is declared beside the body it
# governs, so the check must read it rather than carry its own number:
#   1. body under the declared cap        -> exit 0
#   2. body over the declared cap         -> exit 1, message reports both numbers
#   3. body exactly at the cap            -> exit 0 (boundary is inclusive)
#   4. a non-default cap is honoured      -> the cap comes from craft.yaml
#   5. no craft.yaml at all               -> exit 0, silent
#   6. the shipped craft.yaml is in budget
set -u
. "$(dirname "$0")/_lib.sh"

setup_tmpdir

# _plant <max-lines> <body-lines> — a fresh fake toolkit whose craft.yaml
# declares <max-lines> and renders <body-lines> non-blank lines.
# Echoes the toolkit root.
_plant() {
    local max="$1" lines="$2"
    local root="$TMPDIR_TEST/tk_$RANDOM$RANDOM"
    build_fake_toolkit "$root"
    {
        printf 'version: 1\n'
        printf 'max_lines: %s\n' "$max"
        printf 'floors:\n  figure_base_size: 14\n'
        printf 'body: |\n'
        local i
        for (( i=1; i<=lines; i++ )); do
            printf '  - bullet %s at {{figure_base_size}}pt\n' "$i"
        done
    } > "$root/craft.yaml"
    printf '%s\n' "$root"
}

# _run <toolkit-root> — run the catalog check and capture output and status.
_run() {
    set +e
    out=$(SCIO_TOOLKIT="$1" "$1/bin/scio" lint --check toolkit 2>&1)
    rc=$?
    set -e
}

# --- 1. under the cap --------------------------------------------------------
UNDER=$(_plant 25 17)
_run "$UNDER"
if [[ "$rc" -ne 0 ]]; then
    echo "FAIL [$_TEST_NAME] case1: a 17-line body under a 25-line cap failed lint" >&2
    printf '%s\n' "$out" >&2
    exit 1
fi

# --- 2. over the cap ---------------------------------------------------------
OVER=$(_plant 25 26)
_run "$OVER"
if [[ "$rc" -eq 0 ]]; then
    echo "FAIL [$_TEST_NAME] case2: a 26-line body passed a 25-line cap" >&2
    exit 1
fi
if ! printf '%s\n' "$out" | grep -q 'CRAFT body is 26 lines (max 25)'; then
    echo "FAIL [$_TEST_NAME] case2: message does not report the measured length and the cap" >&2
    printf '%s\n' "$out" >&2
    exit 1
fi

# --- 3. exactly at the cap ---------------------------------------------------
AT=$(_plant 25 25)
_run "$AT"
if [[ "$rc" -ne 0 ]]; then
    echo "FAIL [$_TEST_NAME] case3: a 25-line body rejected (cap must be inclusive)" >&2
    printf '%s\n' "$out" >&2
    exit 1
fi

# --- 4. the cap is read from craft.yaml, not hardcoded -----------------------
# A body of 5 lines is far under the documented default and must still fail
# against a craft.yaml that declares 4.
TIGHT=$(_plant 4 5)
_run "$TIGHT"
if [[ "$rc" -eq 0 ]]; then
    echo "FAIL [$_TEST_NAME] case4: a declared cap of 4 was ignored" >&2
    exit 1
fi
if ! printf '%s\n' "$out" | grep -q 'CRAFT body is 5 lines (max 4)'; then
    echo "FAIL [$_TEST_NAME] case4: the cap does not come from craft.yaml" >&2
    printf '%s\n' "$out" >&2
    exit 1
fi

# --- 5. no craft.yaml is silent, not a finding -------------------------------
# A toolkit without a craft SSOT renders no block, so there is nothing to bound.
BARE="$TMPDIR_TEST/tk_bare"
build_fake_toolkit "$BARE"
_run "$BARE"
if [[ "$rc" -ne 0 ]]; then
    echo "FAIL [$_TEST_NAME] case5: a toolkit with no craft.yaml failed lint" >&2
    printf '%s\n' "$out" >&2
    exit 1
fi
if printf '%s\n' "$out" | grep -qi 'craft'; then
    echo "FAIL [$_TEST_NAME] case5: absent craft.yaml produced output" >&2
    printf '%s\n' "$out" >&2
    exit 1
fi

# --- 6. the shipped craft.yaml is within its own budget ----------------------
_run "$TOOLKIT_ROOT"
if [[ "$rc" -ne 0 ]]; then
    echo "FAIL [$_TEST_NAME] case6: the shipped toolkit fails its own catalog check" >&2
    printf '%s\n' "$out" >&2
    exit 1
fi

pass
