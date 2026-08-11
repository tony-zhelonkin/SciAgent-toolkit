#!/usr/bin/env bash
# tests/test_update_toolkit_guard.sh
# Defect (arch review #4): bin/sciagent's MUTATING_VERB map guards `activate`
# against running with an EXTERNAL toolkit while the project ships its own at
# ./01_modules/SciAgent-toolkit (see test_toolkit_locality_guard.sh), but
# `update` was missing from that map even though update.sh calls
# `cmd_activate` in-process (both the --no-pin path and the freshly re-exec'd
# process on the re-pin path) — bypassing the dispatcher-level guard entirely.
# `SCIAGENT_TOOLKIT=<external> sciagent update --no-pin` therefore silently
# re-mounted the project against the wrong toolkit, escaping the submodule
# pin the guard exists to protect, even though direct `activate` correctly
# refused the same mismatch.
#
# Covers:
#   (a) external toolkit + in-repo toolkit present → `update --no-pin`
#       REFUSES (rc 1), nothing re-mounted, message names the in-repo path.
#   (b) SCIAGENT_ALLOW_EXTERNAL_TOOLKIT=1 bypasses the guard for `update` too.
#   (c) no in-repo toolkit → guard is a no-op; `update --no-pin` against an
#       external-only setup still works (no regression for PATH/dev-symlink
#       installs).
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
IN_REPO_SCIAGENT="$IN_REPO/bin/sciagent"

cd "$PROJ"

# Establish an active stack the normal (in-repo) way first. _lib.sh already
# exported SCIAGENT_TOOLKIT=$TOOLKIT_ROOT (the real toolkit) for its own
# purposes, so it must be overridden explicitly here — bin/sciagent honors an
# already-exported SCIAGENT_TOOLKIT over its own self-location.
export SCIAGENT_TOOLKIT="$IN_REPO"
"$IN_REPO_SCIAGENT" activate base >/dev/null
assert_symlink .claude/skills/s_a "sanity: activated against the in-repo toolkit"
before_target=$(readlink .claude/skills/s_a)

# ---------------------------------------------------------------------------
# (a) external toolkit while in-repo exists → `update --no-pin` REFUSES
# ---------------------------------------------------------------------------
export SCIAGENT_TOOLKIT="$EXT_TK"
set +e
err=$("$EXT_SCIAGENT" update --no-pin 2>&1 >/dev/null)
rc=$?
set -e

if [[ "$rc" -eq 0 ]]; then
    echo "FAIL [$_TEST_NAME] update --no-pin against external toolkit should have failed" >&2
    echo "$err" >&2
    exit 1
fi
case "$err" in
    *"refusing to update against an external toolkit"*) : ;;
    *) echo "FAIL [$_TEST_NAME] guard message missing/unexpected:" >&2; echo "$err" >&2; exit 1 ;;
esac

# The mount must be UNCHANGED — still pointing at the in-repo toolkit, not
# silently re-pointed at the external one.
after_target=$(readlink .claude/skills/s_a)
assert_eq "$after_target" "$before_target" "guard fired but the mount target changed anyway"

# ---------------------------------------------------------------------------
# (b) SCIAGENT_ALLOW_EXTERNAL_TOOLKIT=1 bypasses the guard for `update` too.
# ---------------------------------------------------------------------------
out=$(SCIAGENT_ALLOW_EXTERNAL_TOOLKIT=1 "$EXT_SCIAGENT" update --no-pin 2>&1)
rc=$?
if [[ "$rc" -ne 0 ]]; then
    echo "FAIL [$_TEST_NAME] override env var should have allowed update --no-pin" >&2
    printf '%s\n' "$out" >&2
    exit 1
fi
assert_symlink .claude/skills/s_a "override still mounted after update --no-pin"

# ---------------------------------------------------------------------------
# (c) no in-repo toolkit → guard is a no-op for `update` (no regression for
# external-only installs).
# ---------------------------------------------------------------------------
unset SCIAGENT_TOOLKIT
cd "$TMPDIR_TEST"
NOREPO="$TMPDIR_TEST/norepo_proj"
mkdir -p "$NOREPO" && cd "$NOREPO"
export SCIAGENT_TOOLKIT="$EXT_TK"
"$EXT_SCIAGENT" activate base >/dev/null
out=$("$EXT_SCIAGENT" update --no-pin 2>&1)
rc=$?
if [[ "$rc" -ne 0 ]]; then
    echo "FAIL [$_TEST_NAME] update --no-pin without an in-repo toolkit should succeed" >&2
    printf '%s\n' "$out" >&2
    exit 1
fi

pass
