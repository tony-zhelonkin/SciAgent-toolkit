#!/usr/bin/env bash
# tests/test_deactivate_overlay_guard.sh
# Bonus finding from the arch review's defect #4 ("check whether any other
# verb re-enters a mutating code path the same way"): `deactivate <overlay>`
# calls `cmd_activate "$base"` IN-PROCESS to re-activate solo-base
# (deactivate.sh) — the exact same shape as update.sh's bypass of the
# dispatcher-level guard. `deactivate` is deliberately NOT in MUTATING_VERB
# (its delete-only paths are safe against any toolkit — see
# test_deactivate_external_toolkit.sh), so bin/sciagent's guard never fires
# for this sub-path either. Without an in-process check,
# `SCIAGENT_TOOLKIT=<external> sciagent deactivate <overlay>` would silently
# re-mount solo-base against the wrong toolkit.
#
# Covers:
#   (a) external toolkit + in-repo toolkit present → `deactivate <overlay>`
#       REFUSES (rc 1); stack/mounts unchanged (still base+overlay, still
#       pointing at the in-repo toolkit).
#   (b) SCIAGENT_ALLOW_EXTERNAL_TOOLKIT=1 bypasses the check.
#   (c) removing the BASE (not the overlay) is a delete-only path and must
#       still succeed against an external toolkit, unaffected by this guard.
set -u
. "$(dirname "$0")/_lib.sh"

setup_tmpdir

PROJ="$TMPDIR_TEST/proj"
IN_REPO="$PROJ/01_modules/SciAgent-toolkit"
mkdir -p "$IN_REPO"
build_fake_toolkit "$IN_REPO"

EXT_TK="$TMPDIR_TEST/ext-toolkit"
build_fake_toolkit "$EXT_TK"

cd "$PROJ"
export SCIAGENT_TOOLKIT="$IN_REPO"
"$IN_REPO/bin/sciagent" activate base reviewer >/dev/null
assert_symlink .claude/skills/s_a "sanity: base+reviewer activated in-repo"
before_target=$(readlink .claude/skills/s_a)

# ---------------------------------------------------------------------------
# (a) external toolkit + in-repo present → deactivate <overlay> REFUSES.
# ---------------------------------------------------------------------------
export SCIAGENT_TOOLKIT="$EXT_TK"
set +e
err=$("$EXT_TK/bin/sciagent" deactivate reviewer 2>&1 >/dev/null)
rc=$?
set -e

if [[ "$rc" -eq 0 ]]; then
    echo "FAIL [$_TEST_NAME] deactivate <overlay> against external toolkit should have failed" >&2
    echo "$err" >&2
    exit 1
fi
case "$err" in
    *"refusing to deactivate"*"external toolkit"*) : ;;
    *) echo "FAIL [$_TEST_NAME] guard message missing/unexpected:" >&2; echo "$err" >&2; exit 1 ;;
esac

# Stack must be UNCHANGED — still base+reviewer, still pointing in-repo.
export SCIAGENT_TOOLKIT="$IN_REPO"
. "$IN_REPO/lib/sciagent/symlinks.sh"
stack_after=$(manifest_stack)
assert_eq "$stack_after" "base reviewer" "guard fired but the stack changed anyway"
after_target=$(readlink .claude/skills/s_a)
assert_eq "$after_target" "$before_target" "guard fired but the mount target changed anyway"

# ---------------------------------------------------------------------------
# (b) SCIAGENT_ALLOW_EXTERNAL_TOOLKIT=1 bypasses the check.
# ---------------------------------------------------------------------------
export SCIAGENT_TOOLKIT="$EXT_TK"
out=$(SCIAGENT_ALLOW_EXTERNAL_TOOLKIT=1 "$EXT_TK/bin/sciagent" deactivate reviewer 2>&1)
rc=$?
if [[ "$rc" -ne 0 ]]; then
    echo "FAIL [$_TEST_NAME] override env var should have allowed deactivate <overlay>" >&2
    printf '%s\n' "$out" >&2
    exit 1
fi
assert_symlink .claude/skills/s_a "override re-activated solo-base"

# ---------------------------------------------------------------------------
# (c) removing the BASE is delete-only and must still work against an
# external toolkit — this guard must not have over-broadened.
# ---------------------------------------------------------------------------
export SCIAGENT_TOOLKIT="$IN_REPO"
"$IN_REPO/bin/sciagent" activate base reviewer >/dev/null

export SCIAGENT_TOOLKIT="$EXT_TK"
out=$("$EXT_TK/bin/sciagent" deactivate base 2>&1)
rc=$?
if [[ "$rc" -ne 0 ]]; then
    echo "FAIL [$_TEST_NAME] deactivate <base> (delete-only) against external toolkit should succeed" >&2
    printf '%s\n' "$out" >&2
    exit 1
fi
if [[ -e .claude/skills/s_a ]]; then
    echo "FAIL [$_TEST_NAME] deactivate <base> should have torn everything down" >&2
    exit 1
fi

pass
