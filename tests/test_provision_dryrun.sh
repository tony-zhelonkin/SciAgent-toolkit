#!/usr/bin/env bash
# tests/test_provision_dryrun.sh — `sciagent provision --dry-run` writes nothing,
# and `--harness codex` (absent) warns + exits without writing.
#
#   - --dry-run: snapshot the scratch HOME before/after; assert identical, while
#     the run still prints the intended actions.
#   - --harness codex (uninstalled): warns "requested but not detected", ends
#     with "nothing to do", exit 0, and writes nothing.
#
# Hermetic: scratch HOME under the tmpdir, PATH pruned so detection sees only the
# fake config dirs we create.
set -u
. "$(dirname "$0")/_lib.sh"
export PATH=/usr/bin:/bin

SCIAGENT="$TOOLKIT_ROOT/bin/sciagent"

# snapshot — recursive manifest (paths + per-file content hashes) of $HOME.
snapshot() {
    ( cd "$HOME" && find . | LC_ALL=C sort && find . -type f -exec sha1sum {} \; | LC_ALL=C sort )
}

setup_tmpdir
unset CLAUDE_CONFIG_DIR CODEX_HOME
export HOME="$TMPDIR_TEST/home"
mkdir -p "$HOME/.claude" "$HOME/.pi"

# ----- --dry-run writes nothing but prints intended actions ----------------
before=$(snapshot)
out=$("$SCIAGENT" provision --dry-run 2>&1) \
    || { echo "FAIL [$_TEST_NAME] provision --dry-run exited non-zero" >&2; printf '%s\n' "$out" >&2; exit 1; }
after=$(snapshot)

assert_eq "$after" "$before" "dry-run must not modify the scratch HOME"

printf '%s\n' "$out" | grep -q '\[dry-run\]' \
    || { echo "FAIL [$_TEST_NAME] dry-run output missing '[dry-run]' marker" >&2; printf '%s\n' "$out" >&2; exit 1; }
# It should still describe intended context + settings actions for detected harnesses.
printf '%s\n' "$out" | grep -q 'would stamp' \
    || { echo "FAIL [$_TEST_NAME] dry-run did not announce the context stamp" >&2; printf '%s\n' "$out" >&2; exit 1; }
printf '%s\n' "$out" | grep -q 'would ensure' \
    || { echo "FAIL [$_TEST_NAME] dry-run did not announce the settings ensure" >&2; printf '%s\n' "$out" >&2; exit 1; }

# ----- --harness codex (absent) warns + exits without writing --------------
before2=$(snapshot)
out2=$("$SCIAGENT" provision --harness codex 2>&1)
rc=$?
after2=$(snapshot)

assert_eq "$rc" "0" "provision --harness codex (absent) exits 0"
assert_eq "$after2" "$before2" "absent-harness request must not write anything"
printf '%s\n' "$out2" | grep -q "requested but not detected" \
    || { echo "FAIL [$_TEST_NAME] absent harness: missing 'requested but not detected' warning" >&2; printf '%s\n' "$out2" >&2; exit 1; }
printf '%s\n' "$out2" | grep -q "nothing to do" \
    || { echo "FAIL [$_TEST_NAME] absent harness: missing 'nothing to do' notice" >&2; printf '%s\n' "$out2" >&2; exit 1; }

pass
