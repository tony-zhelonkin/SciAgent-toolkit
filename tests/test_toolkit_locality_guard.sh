#!/usr/bin/env bash
# tests/test_toolkit_locality_guard.sh
# Reproducibility fix 2: the toolkit-locality guard refuses to run a mutating
# verb against an EXTERNAL toolkit when the project ships its own — and the
# override escape hatches bypass it.
#
# Covers:
#   (a) external toolkit + in-repo toolkit present → activate REFUSES (rc 1),
#       nothing written, message names the in-repo path.
#   (b) --allow-external-toolkit flag bypasses the guard (activate succeeds).
#   (c) SCIAGENT_ALLOW_EXTERNAL_TOOLKIT=1 env var bypasses the guard.
#   (d) no in-repo toolkit → guard is a no-op (external toolkit activates fine).
#   (e) NON-STANDARD container dir (01_Modules, 01_scripts) → guard still fires,
#       and the remedy message names the real path. Added 2026-08-11.
#   (f) no .gitmodules → discovery falls back to the container scan.
#   (g) `new` is guarded too — it writes INTO the toolkit, the mirror of the
#       case (a) protects.
#
# Cases (e)-(g) are regressions, not new features: the guard hardcoded
# `./01_modules/SciAgent-toolkit` and `new` was missing from MUTATING_VERB,
# leaving 5 live consumer projects unprotected — two of them in a shared lab
# tree. Every pre-existing fixture used `01_modules`, so the suite stayed green
# the whole time. The blind spot was the fixture, not the logic.
set -u
. "$(dirname "$0")/_lib.sh"

setup_tmpdir

# In-repo toolkit that the project ships.
PROJ="$TMPDIR_TEST/proj"
IN_REPO="$PROJ/01_modules/SciAgent-toolkit"
mkdir -p "$IN_REPO"
build_fake_toolkit "$IN_REPO"

# A separate EXTERNAL/global toolkit.
EXT_TK="$TMPDIR_TEST/ext-toolkit"
build_fake_toolkit "$EXT_TK"
EXT_SCIAGENT="$EXT_TK/bin/sciagent"

cd "$PROJ"

# ---------------------------------------------------------------------------
# (a) external toolkit while in-repo exists → REFUSE
# ---------------------------------------------------------------------------
export SCIAGENT_TOOLKIT="$EXT_TK"
if err="$("$EXT_SCIAGENT" activate base 2>&1 >/dev/null)"; then
    echo "FAIL [$_TEST_NAME] activate against external toolkit should have failed" >&2
    echo "$err" >&2
    exit 1
fi
case "$err" in
    *"refusing to activate against an external toolkit"*) : ;;
    *) echo "FAIL [$_TEST_NAME] guard message missing/unexpected:" >&2; echo "$err" >&2; exit 1 ;;
esac
# Nothing should have been written.
if [[ -e ".claude/skills/s_a" || -e ".sciagent/manifest.json" ]]; then
    echo "FAIL [$_TEST_NAME] guard fired but files were still written" >&2
    exit 1
fi

# ---------------------------------------------------------------------------
# (b) --allow-external-toolkit bypasses the guard
# ---------------------------------------------------------------------------
"$EXT_SCIAGENT" activate base --allow-external-toolkit >/dev/null
assert_symlink ".claude/skills/s_a" "override flag should have allowed activation"
"$EXT_SCIAGENT" deactivate --allow-external-toolkit >/dev/null

# ---------------------------------------------------------------------------
# (c) SCIAGENT_ALLOW_EXTERNAL_TOOLKIT=1 bypasses the guard
# ---------------------------------------------------------------------------
SCIAGENT_ALLOW_EXTERNAL_TOOLKIT=1 "$EXT_SCIAGENT" activate base >/dev/null
assert_symlink ".claude/skills/s_a" "override env var should have allowed activation"
"$EXT_SCIAGENT" deactivate >/dev/null

# ---------------------------------------------------------------------------
# (d) no in-repo toolkit → guard is a no-op
# ---------------------------------------------------------------------------
NOREPO="$TMPDIR_TEST/norepo_proj"
mkdir -p "$NOREPO"
cd "$NOREPO"
export SCIAGENT_TOOLKIT="$EXT_TK"
"$EXT_SCIAGENT" activate base >/dev/null   # must succeed (rc 0), no guard
assert_symlink ".claude/skills/s_a" "activation without an in-repo toolkit should succeed"

# ---------------------------------------------------------------------------
# (e) NON-STANDARD CONTAINER DIRECTORY → guard must still fire.
#
# Regression test for a defect found 2026-08-11: the guard hardcoded the
# literal `./01_modules/SciAgent-toolkit` and matched it CASE-SENSITIVELY,
# while the real fleet uses four different container directories — measured
# across 23 checkouts: 17 `01_modules`, 2 `01_Scripts`, 2 `01_scripts`,
# 1 `01_Modules`. In the six non-conforming projects the `-d` test failed, the
# guard concluded "no in-repo toolkit, nothing to protect", and mutation from
# ANY external checkout was permitted — silently, and including two projects in
# a shared multi-user lab tree.
#
# The whole suite passed throughout, because every fixture until now used
# `01_modules`. That is why this case exists: the absent-fixture blind spot,
# not the guard logic, is what let five live consumers go unprotected.
#
# `01_Modules` (capital M) is used deliberately — it is the case that proves
# the match is not case-sensitive, and it is a real fleet path.
# ---------------------------------------------------------------------------
_assert_guard_fires_in() {
    local container="$1"
    local proj="$TMPDIR_TEST/proj_${container}"
    local in_repo="$proj/$container/SciAgent-toolkit"
    mkdir -p "$in_repo"
    build_fake_toolkit "$in_repo"

    # git's own record of the submodule location — the primary resolution
    # source, and the one verified correct in all four real consumers.
    cat > "$proj/.gitmodules" <<EOF
[submodule "$container/SciAgent-toolkit"]
	path = $container/SciAgent-toolkit
	url = https://example.invalid/SciAgent-toolkit
EOF

    cd "$proj"
    export SCIAGENT_TOOLKIT="$EXT_TK"
    local err
    if err="$("$EXT_SCIAGENT" activate base 2>&1 >/dev/null)"; then
        echo "FAIL [$_TEST_NAME] (e) guard did NOT fire for container '$container' — external activate succeeded" >&2
        exit 1
    fi
    case "$err" in
        *"refusing to activate against an external toolkit"*) : ;;
        *) echo "FAIL [$_TEST_NAME] (e) wrong error for '$container':" >&2; echo "$err" >&2; exit 1 ;;
    esac
    # The message must name the REAL path, not a hardcoded guess — otherwise the
    # remedy it suggests sends the user to a directory that does not exist.
    case "$err" in
        *"$container/SciAgent-toolkit"*) : ;;
        *) echo "FAIL [$_TEST_NAME] (e) message does not name '$container/SciAgent-toolkit':" >&2
           echo "$err" >&2; exit 1 ;;
    esac
    if [[ -e ".claude/skills/s_a" || -e ".sciagent/manifest.json" ]]; then
        echo "FAIL [$_TEST_NAME] (e) guard fired for '$container' but files were written" >&2
        exit 1
    fi
}

_assert_guard_fires_in "01_Modules"
_assert_guard_fires_in "01_scripts"

# ---------------------------------------------------------------------------
# (f) Discovery must work with NO .gitmodules at all — a vendored copy, or a
# scaffold that is not yet a git repo. Falls through to the container scan.
# ---------------------------------------------------------------------------
NOGM="$TMPDIR_TEST/proj_nogitmodules"
mkdir -p "$NOGM/01_Scripts/SciAgent-toolkit"
build_fake_toolkit "$NOGM/01_Scripts/SciAgent-toolkit"
cd "$NOGM"
export SCIAGENT_TOOLKIT="$EXT_TK"
if err="$("$EXT_SCIAGENT" activate base 2>&1 >/dev/null)"; then
    echo "FAIL [$_TEST_NAME] (f) guard did not fire without .gitmodules" >&2
    exit 1
fi
case "$err" in
    *"refusing to activate against an external toolkit"*) : ;;
    *) echo "FAIL [$_TEST_NAME] (f) wrong error:" >&2; echo "$err" >&2; exit 1 ;;
esac

# ---------------------------------------------------------------------------
# (g) `new` is a MUTATING verb — it writes INTO $SCIAGENT_TOOLKIT (roles/,
# skills/, agents/), so an external toolkit must be refused for it too. It was
# absent from MUTATING_VERB until 2026-08-11, so
# `SCIAGENT_TOOLKIT=<external> sciagent new skill foo` silently authored into
# another checkout, leaving it differing from its recorded commit.
# ---------------------------------------------------------------------------
# build_fake_toolkit ships no templates/skill, so `new skill` would fail on a
# missing template whether or not the guard fired — which would make this case
# pass for the wrong reason. Give the EXTERNAL toolkit a real template so the
# UNGUARDED path genuinely succeeds and writes; then a refusal is proof the
# guard stopped an actual write, not proof that something else broke first.
mkdir -p "$EXT_TK/templates/skill"
cat > "$EXT_TK/templates/skill/SKILL.md" <<'EOF'
---
name: SKILL_IDENTIFIER
description: Fixture template for the new-verb guard case.
---

template body
EOF

cd "$PROJ"
export SCIAGENT_TOOLKIT="$EXT_TK"
if err="$("$EXT_SCIAGENT" new skill guard_probe_skill 2>&1 >/dev/null)"; then
    echo "FAIL [$_TEST_NAME] (g) 'new' against an external toolkit should have been refused" >&2
    exit 1
fi
case "$err" in
    *"refusing to new against an external toolkit"*) : ;;
    *) echo "FAIL [$_TEST_NAME] (g) wrong error for new:" >&2; echo "$err" >&2; exit 1 ;;
esac
# And it must not have authored anything into the external toolkit.
if [[ -e "$EXT_TK/skills/guard_probe_skill" ]]; then
    echo "FAIL [$_TEST_NAME] (g) guard fired but 'new' still wrote into the external toolkit" >&2
    exit 1
fi

pass
