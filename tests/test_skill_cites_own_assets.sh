#!/usr/bin/env bash
# tests/test_skill_cites_own_assets.sh — a skill citing its own asset has it.
#
# skill-creator told a reader to open `assets/eval_review.html`. The file is
# `assets/eval-review.html`. One character, and the instruction could not be
# followed — invisible because nothing read a skill's own citations.
#
# Ownership is what makes this checkable. `assets/x` inside a skill that HAS an
# assets/ directory is a claim about its own tree. `scripts/y` inside a skill
# with no scripts/ directory names somebody else's tree, normally a vendored
# toolkit: star-te-preprocessing cites the TE toolkit's
# scripts/runFeatureCounts_TE_and_genes.sh and owns no scripts/ at all. Flagging
# that would make the check useless, so the directory's existence is the gate.
#
# Tests:
#   1. cited file present in an owned directory     -> exit 0
#   2. cited file ABSENT from an owned directory    -> exit 1, names skill + path
#   3. the directory is not owned at all            -> silent (external citation)
#   4. a nested path inside an owned directory      -> checked
#   5. a cross-skill citation carrying its owner    -> silent
#   6. the shipped catalog passes
set -u
. "$(dirname "$0")/_lib.sh"

setup_tmpdir

_run() {
    set +e
    out=$(SCIO_TOOLKIT="$1" "$1/bin/scio" lint --check toolkit 2>&1)
    rc=$?
    set -e
}

_fresh() {
    local root="$TMPDIR_TEST/tk_$RANDOM$RANDOM"
    build_fake_toolkit "$root"
    printf '%s\n' "$root"
}

# _skill <root> <name> <body>
_skill() {
    mkdir -p "$1/skills/$2"
    cat > "$1/skills/$2/SKILL.md" <<EOF
---
name: $2
description: A fixture skill.
---

$3
EOF
}

# --- 1. present in an owned directory --------------------------------------
TK=$(_fresh)
_skill "$TK" demo 'Open `assets/report.html` to review.'
mkdir -p "$TK/skills/demo/assets"
echo '<html></html>' > "$TK/skills/demo/assets/report.html"
_run "$TK"
if [[ "$rc" -ne 0 ]]; then
    echo "FAIL [$_TEST_NAME] case1: a present asset was reported" >&2
    printf '%s\n' "$out" >&2
    exit 1
fi

# --- 2. absent from an owned directory ------------------------------------
TK=$(_fresh)
_skill "$TK" demo 'Open `assets/report.html` to review.'
mkdir -p "$TK/skills/demo/assets"
echo 'x' > "$TK/skills/demo/assets/other.html"
_run "$TK"
if [[ "$rc" -eq 0 ]]; then
    echo "FAIL [$_TEST_NAME] case2: a missing asset passed" >&2
    exit 1
fi
if ! printf '%s\n' "$out" | grep -q 'demo: cites assets/report.html'; then
    echo "FAIL [$_TEST_NAME] case2: the finding does not name skill and path" >&2
    printf '%s\n' "$out" >&2
    exit 1
fi

# --- 3. an unowned directory is somebody else's tree ---------------------
TK=$(_fresh)
_skill "$TK" demo 'The vendored driver is `scripts/runFeatureCounts.sh`.'
_run "$TK"
if [[ "$rc" -ne 0 ]]; then
    echo "FAIL [$_TEST_NAME] case3: an external citation was flagged" >&2
    printf '%s\n' "$out" >&2
    exit 1
fi

# --- 4. nesting is followed ---------------------------------------------
TK=$(_fresh)
_skill "$TK" demo 'See `references/deep/topic.md`.'
mkdir -p "$TK/skills/demo/references/deep"
_run "$TK"
if [[ "$rc" -eq 0 ]]; then
    echo "FAIL [$_TEST_NAME] case4: a missing nested reference passed" >&2
    exit 1
fi
echo 'note' > "$TK/skills/demo/references/deep/topic.md"
_run "$TK"
if [[ "$rc" -ne 0 ]]; then
    echo "FAIL [$_TEST_NAME] case4: a present nested reference was reported" >&2
    printf '%s\n' "$out" >&2
    exit 1
fi

# --- 5. a citation that names its owner is not a self-claim -------------
# This is how a skill should point at a neighbour's file, and it is the form
# te-reference-saf-build now uses for star-te-preprocessing's reference.
TK=$(_fresh)
_skill "$TK" demo 'See `other-skill/references/topic.md` for the contract.'
mkdir -p "$TK/skills/demo/references"
echo 'own' > "$TK/skills/demo/references/mine.md"
_run "$TK"
if [[ "$rc" -ne 0 ]]; then
    echo "FAIL [$_TEST_NAME] case5: a qualified cross-skill citation was flagged" >&2
    printf '%s\n' "$out" >&2
    exit 1
fi

# --- 6. the shipped catalog ---------------------------------------------
_run "$TOOLKIT_ROOT"
if [[ "$rc" -ne 0 ]]; then
    echo "FAIL [$_TEST_NAME] case6: the shipped catalog fails the check" >&2
    printf '%s\n' "$out" >&2
    exit 1
fi

pass
