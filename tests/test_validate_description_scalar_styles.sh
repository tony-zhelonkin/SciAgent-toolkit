#!/usr/bin/env bash
# tests/test_validate_description_scalar_styles.sh — the `description:` length
# check in `sciagent validate` must measure the value a YAML parser actually
# produces, across every scalar style a skill author might reach for: plain,
# single-quoted, double-quoted, and folded/literal block scalars (`>-`, `|-`).
#
# This is a regression test for a real bug: the original awk measured only
# the FIRST PHYSICAL LINE after `description:`. For a folded scalar written
# as
#     description: >-
#       <content on continuation lines>
# that first line is literally the 2-char style indicator `>-`, so the
# measured length was 2 regardless of how long the folded value actually was
# — `sciagent validate` exited 0 no matter how far over SCIAGENT_DESC_MAX the
# real description ran. Five in-repo skills were over cap and undetected
# before the fix (see docs/changelog.md or the commit that added this test).
#
# Each style below is checked at exactly the cap (350, must pass) and one
# char over (351, must fail with the correctly measured length in the
# message) — the folded case is the one that catches the original bug: on
# the pre-fix parser, ALL of the folded sub-cases here pass (wrongly),
# because the measured length is always ~2.
set -u
. "$(dirname "$0")/_lib.sh"

setup_tmpdir

# _plant <frontmatter-body> — rewrite skills/s_a in a fresh fake toolkit.
# Echoes the toolkit root. (Same helper as test_validate_frontmatter_shape.sh;
# duplicated locally so this file has no cross-file dependency.)
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

# _run <toolkit-root> — run validate, capture combined output into $out, rc into $rc.
_run() {
    set +e
    out=$(SCIAGENT_TOOLKIT="$1" "$1/bin/sciagent" validate 2>&1)
    rc=$?
    set -e
}

# _fill <n> — n 'x' characters, no trailing newline.
_fill() { head -c "$1" < /dev/zero | tr '\0' 'x'; }

# ---------------------------------------------------------------------------
# Single-quoted scalar: at cap passes, one over fails.
# ---------------------------------------------------------------------------
SQ_AT=$(_plant "name: s_a
description: '$(_fill 350)'")
_run "$SQ_AT"
if [[ "$rc" -ne 0 ]]; then
    echo "FAIL [$_TEST_NAME] single-quoted@350: expected exit 0, got $rc" >&2
    printf '%s\n' "$out" >&2
    exit 1
fi

SQ_OVER=$(_plant "name: s_a
description: '$(_fill 351)'")
_run "$SQ_OVER"
if [[ "$rc" -eq 0 ]]; then
    echo "FAIL [$_TEST_NAME] single-quoted@351: expected nonzero exit, got 0" >&2
    exit 1
fi
if ! printf '%s\n' "$out" | grep -q 'description is 351 chars (max 350)'; then
    echo "FAIL [$_TEST_NAME] single-quoted@351: message does not report 351 chars" >&2
    printf '%s\n' "$out" >&2
    exit 1
fi

# ---------------------------------------------------------------------------
# Double-quoted scalar: at cap passes, one over fails.
# ---------------------------------------------------------------------------
DQ_AT=$(_plant "name: s_a
description: \"$(_fill 350)\"")
_run "$DQ_AT"
if [[ "$rc" -ne 0 ]]; then
    echo "FAIL [$_TEST_NAME] double-quoted@350: expected exit 0, got $rc" >&2
    printf '%s\n' "$out" >&2
    exit 1
fi

DQ_OVER=$(_plant "name: s_a
description: \"$(_fill 351)\"")
_run "$DQ_OVER"
if [[ "$rc" -eq 0 ]]; then
    echo "FAIL [$_TEST_NAME] double-quoted@351: expected nonzero exit, got 0" >&2
    exit 1
fi
if ! printf '%s\n' "$out" | grep -q 'description is 351 chars (max 350)'; then
    echo "FAIL [$_TEST_NAME] double-quoted@351: message does not report 351 chars" >&2
    printf '%s\n' "$out" >&2
    exit 1
fi

# ---------------------------------------------------------------------------
# Folded block scalar (>-): content lives on continuation lines below the
# `description:` line, split across two lines so the fix's line-accumulation
# (not just its per-line stripping) is exercised. Folding joins the two
# lines with a single space, so len(line1) + 1 + len(line2) is the measured
# length under `>-` (strip chomping, no trailing newline).
#
# THIS is the exact shape of the original bug: the naive old parser read
# only the physical line right after `description:` (literally the text
# ">-"), measured length 2, and passed regardless of the real folded length.
# ---------------------------------------------------------------------------
FOLD_AT=$(_plant "name: s_a
description: >-
  $(_fill 150)
  $(_fill 199)")
_run "$FOLD_AT"
if [[ "$rc" -ne 0 ]]; then
    echo "FAIL [$_TEST_NAME] folded>-@350: expected exit 0, got $rc" >&2
    printf '%s\n' "$out" >&2
    exit 1
fi

FOLD_OVER=$(_plant "name: s_a
description: >-
  $(_fill 150)
  $(_fill 200)")
_run "$FOLD_OVER"
if [[ "$rc" -eq 0 ]]; then
    echo "FAIL [$_TEST_NAME] folded>-@351: expected nonzero exit, got 0 -- this is the exact bug this test guards against (old parser measured the '>-' indicator line, length 2, and always passed)" >&2
    exit 1
fi
if ! printf '%s\n' "$out" | grep -q 'description is 351 chars (max 350)'; then
    echo "FAIL [$_TEST_NAME] folded>-@351: message does not report 351 chars" >&2
    printf '%s\n' "$out" >&2
    exit 1
fi

# ---------------------------------------------------------------------------
# Literal block scalar (|-): continuation lines are joined with a literal
# newline (not folded to a space), so length is len(line1) + 1 + len(line2)
# exactly as with the folded case above -- exercises the `|` branch of
# finalize() specifically.
# ---------------------------------------------------------------------------
LIT_AT=$(_plant "name: s_a
description: |-
  $(_fill 150)
  $(_fill 199)")
_run "$LIT_AT"
if [[ "$rc" -ne 0 ]]; then
    echo "FAIL [$_TEST_NAME] literal|-@350: expected exit 0, got $rc" >&2
    printf '%s\n' "$out" >&2
    exit 1
fi

LIT_OVER=$(_plant "name: s_a
description: |-
  $(_fill 150)
  $(_fill 200)")
_run "$LIT_OVER"
if [[ "$rc" -eq 0 ]]; then
    echo "FAIL [$_TEST_NAME] literal|-@351: expected nonzero exit, got 0" >&2
    exit 1
fi
if ! printf '%s\n' "$out" | grep -q 'description is 351 chars (max 350)'; then
    echo "FAIL [$_TEST_NAME] literal|-@351: message does not report 351 chars" >&2
    printf '%s\n' "$out" >&2
    exit 1
fi

pass
