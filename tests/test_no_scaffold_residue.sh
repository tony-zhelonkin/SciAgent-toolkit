#!/usr/bin/env bash
# tests/test_no_scaffold_residue.sh — `scio lint --check toolkit` refuses a
# directory that claims content it does not have.
#
# A .gitkeep is a claim with nothing behind it, and it propagates: templates/skill/
# shipped four of them, so every skill copied from the template inherited four
# directories asserting emptiness. Nineteen reached five real skills that way, and
# every one of those directories already held real files — so not one of the
# .gitkeep files was load-bearing even on its own terms.
#
# Tests:
#   1. a .gitkeep under skills/            -> exit 1, message names the file
#   2. a .gitkeep under templates/         -> exit 1 (the propagation source)
#   3. an empty directory under skills/    -> exit 1
#   4. a directory with a real file        -> exit 0
#   5. a .gitkeep under docs/_internal/    -> ignored here (internal-memory owns it)
#   6. the shipped tree passes
set -u
. "$(dirname "$0")/_lib.sh"

setup_tmpdir

# _run <toolkit-root> — run the catalog check, capturing output and status.
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

# _skill <root> <name> — a minimal valid skill, so findings come from the residue.
_skill() {
    mkdir -p "$1/skills/$2"
    cat > "$1/skills/$2/SKILL.md" <<EOF
---
name: $2
description: A fixture skill.
---
EOF
}

# --- 1. a .gitkeep under skills/ ---------------------------------------------
TK=$(_fresh)
_skill "$TK" demo
mkdir -p "$TK/skills/demo/references"
touch "$TK/skills/demo/references/.gitkeep"
_run "$TK"
if [[ "$rc" -eq 0 ]]; then
    echo "FAIL [$_TEST_NAME] case1: a .gitkeep under skills/ passed lint" >&2
    exit 1
fi
if ! printf '%s\n' "$out" | grep -q 'skills/demo/references/.gitkeep'; then
    echo "FAIL [$_TEST_NAME] case1: the finding does not name the file" >&2
    printf '%s\n' "$out" >&2
    exit 1
fi

# --- 2. a .gitkeep under templates/ is the propagation source ----------------
TK=$(_fresh)
_skill "$TK" demo
mkdir -p "$TK/templates/skill/scripts"
touch "$TK/templates/skill/scripts/.gitkeep"
_run "$TK"
if [[ "$rc" -eq 0 ]]; then
    echo "FAIL [$_TEST_NAME] case2: a .gitkeep under templates/ passed lint" >&2
    exit 1
fi
if ! printf '%s\n' "$out" | grep -q 'templates/skill/scripts/.gitkeep'; then
    echo "FAIL [$_TEST_NAME] case2: the finding does not name the file" >&2
    printf '%s\n' "$out" >&2
    exit 1
fi

# --- 3. an empty directory is the same claim without the file ---------------
TK=$(_fresh)
_skill "$TK" demo
mkdir -p "$TK/skills/demo/assets"
_run "$TK"
if [[ "$rc" -eq 0 ]]; then
    echo "FAIL [$_TEST_NAME] case3: an empty skill directory passed lint" >&2
    exit 1
fi
if ! printf '%s\n' "$out" | grep -q 'skills/demo/assets'; then
    echo "FAIL [$_TEST_NAME] case3: the finding does not name the directory" >&2
    printf '%s\n' "$out" >&2
    exit 1
fi

# --- 4. a directory that earned itself is silent ----------------------------
TK=$(_fresh)
_skill "$TK" demo
mkdir -p "$TK/skills/demo/references"
echo '# A real reference' > "$TK/skills/demo/references/topic.md"
_run "$TK"
if [[ "$rc" -ne 0 ]]; then
    echo "FAIL [$_TEST_NAME] case4: a populated directory failed lint" >&2
    printf '%s\n' "$out" >&2
    exit 1
fi

# --- 5. the memory tree belongs to internal-memory, not this check ----------
TK=$(_fresh)
_skill "$TK" demo
mkdir -p "$TK/docs/_internal/_project"
touch "$TK/docs/_internal/_project/.gitkeep"
_run "$TK"
if printf '%s\n' "$out" | grep -q '_internal'; then
    echo "FAIL [$_TEST_NAME] case5: the toolkit check reached into docs/_internal/" >&2
    printf '%s\n' "$out" >&2
    exit 1
fi

# --- 6. the shipped tree passes its own check -------------------------------
_run "$TOOLKIT_ROOT"
if [[ "$rc" -ne 0 ]]; then
    echo "FAIL [$_TEST_NAME] case6: the shipped tree fails the check" >&2
    printf '%s\n' "$out" >&2
    exit 1
fi

pass
