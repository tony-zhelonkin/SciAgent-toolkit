#!/usr/bin/env bash
# tests/test_attic_holds_skills.sh — `scio lint --check toolkit` enforces what
# _attic/README.md and docs/skill-lifecycle.md already claim: the attic holds
# retired skills.
#
# The attic sits outside every category mount, so nothing there is reachable by
# a harness and nothing there is read by lint. A non-skill parked in it becomes
# a second home for a claim, ages out of sync with the skill that supersedes it,
# and keeps asserting its own authority while unreadable. That happened: a
# guidelines/ tree called itself the single source of truth for coding
# conventions for months after the router skill took the job.
#
# Tests:
#   1. an attic entry WITH a SKILL.md          -> exit 0
#   2. an attic entry WITHOUT a SKILL.md       -> exit 1, message names it
#   3. a loose FILE at the attic root          -> ignored (README.md lives there)
#   4. no _attic/ directory at all             -> exit 0, silent
#   5. mixed: one valid, one not               -> exit 1
#   6. the shipped _attic/ passes
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

# --- 1. a retired skill is what the attic is for -----------------------------
TK=$(_fresh)
mkdir -p "$TK/_attic/retired-skill"
cat > "$TK/_attic/retired-skill/SKILL.md" <<'EOF'
---
name: retired-skill
description: A retired fixture skill.
---
EOF
_run "$TK"
if [[ "$rc" -ne 0 ]]; then
    echo "FAIL [$_TEST_NAME] case1: a retired skill in the attic failed lint" >&2
    printf '%s\n' "$out" >&2
    exit 1
fi

# --- 2. a non-skill directory is the drift this catches ----------------------
TK=$(_fresh)
mkdir -p "$TK/_attic/guidelines"
echo '# I am the single source of truth' > "$TK/_attic/guidelines/README.md"
_run "$TK"
if [[ "$rc" -eq 0 ]]; then
    echo "FAIL [$_TEST_NAME] case2: a non-skill directory in the attic passed lint" >&2
    exit 1
fi
if ! printf '%s\n' "$out" | grep -q '_attic/guidelines'; then
    echo "FAIL [$_TEST_NAME] case2: the finding does not name the offending entry" >&2
    printf '%s\n' "$out" >&2
    exit 1
fi

# --- 3. a loose file at the attic root is fine (README.md is one) ------------
TK=$(_fresh)
mkdir -p "$TK/_attic"
echo '# Skill attic' > "$TK/_attic/README.md"
_run "$TK"
if [[ "$rc" -ne 0 ]]; then
    echo "FAIL [$_TEST_NAME] case3: a file at the attic root was treated as an entry" >&2
    printf '%s\n' "$out" >&2
    exit 1
fi

# --- 4. no attic at all is silence, not a finding ---------------------------
TK=$(_fresh)
_run "$TK"
if [[ "$rc" -ne 0 ]]; then
    echo "FAIL [$_TEST_NAME] case4: a toolkit with no _attic/ failed lint" >&2
    printf '%s\n' "$out" >&2
    exit 1
fi
if printf '%s\n' "$out" | grep -q '_attic'; then
    echo "FAIL [$_TEST_NAME] case4: an absent attic produced output" >&2
    printf '%s\n' "$out" >&2
    exit 1
fi

# --- 5. one valid entry does not excuse an invalid one ----------------------
TK=$(_fresh)
mkdir -p "$TK/_attic/good" "$TK/_attic/bad"
cat > "$TK/_attic/good/SKILL.md" <<'EOF'
---
name: good
description: A retired fixture skill.
---
EOF
echo 'notes' > "$TK/_attic/bad/notes.md"
_run "$TK"
if [[ "$rc" -eq 0 ]]; then
    echo "FAIL [$_TEST_NAME] case5: a mixed attic passed lint" >&2
    exit 1
fi
if ! printf '%s\n' "$out" | grep -q '_attic/bad'; then
    echo "FAIL [$_TEST_NAME] case5: the finding does not name the offending entry" >&2
    printf '%s\n' "$out" >&2
    exit 1
fi

# --- 6. the shipped attic passes its own check ------------------------------
_run "$TOOLKIT_ROOT"
if [[ "$rc" -ne 0 ]]; then
    echo "FAIL [$_TEST_NAME] case6: the shipped _attic/ fails the check" >&2
    printf '%s\n' "$out" >&2
    exit 1
fi

pass
