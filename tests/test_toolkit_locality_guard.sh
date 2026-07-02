#!/usr/bin/env bash
# tests/test_toolkit_locality_guard.sh
# Reproducibility fix 2: the toolkit-locality guard refuses to run a mutating
# verb (activate/inject/eject) against an EXTERNAL toolkit when the project
# ships its own at ./01_modules/SciAgent-toolkit — and the override escape
# hatches bypass it.
#
# Covers:
#   (a) external toolkit + in-repo toolkit present → activate REFUSES (rc 1),
#       nothing written, message names the in-repo path.
#   (b) --allow-external-toolkit flag bypasses the guard (activate succeeds).
#   (c) SCIAGENT_ALLOW_EXTERNAL_TOOLKIT=1 env var bypasses the guard.
#   (d) no in-repo toolkit → guard is a no-op (external toolkit activates fine).
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

pass
