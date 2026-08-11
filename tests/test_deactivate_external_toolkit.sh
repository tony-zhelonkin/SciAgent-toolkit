#!/usr/bin/env bash
# tests/test_deactivate_external_toolkit.sh
# Defect (arch review #5): symlink_teardown_all resolved ownership purely
# against the CURRENTLY ACTIVE $SCIAGENT_TOOLKIT. A project mounted with
# toolkit A (e.g. its in-repo checkout), then torn down while $SCIAGENT_TOOLKIT
# points at a DIFFERENT toolkit B, found none of its own mounts under B — the
# target-based sweep silently removed nothing — yet the caller still deleted
# the manifest and the ROLES/CRAFT managed blocks and printed "deactivated",
# rc 0. Net effect: every mount stranded, bookkeeping erased, and `status`
# afterward lies ("no active stack") even though the mounts are still there.
#
# Fix: symlink_teardown_all now ALSO sweeps every path recorded in a
# still-present manifest (first-party bookkeeping — safe to trust
# unconditionally, unlike inferring ownership from a target path) as a UNION
# with the target-based sweep, so teardown succeeds regardless of which
# toolkit happens to be active at deactivate time.
#
# Reproduction, mirroring the team lead's repro command
# (`SCIAGENT_TOOLKIT=<other checkout> sciagent deactivate`): activate with
# toolkit A, then deactivate with $SCIAGENT_TOOLKIT pointed at a completely
# separate toolkit B.
set -u
. "$(dirname "$0")/_lib.sh"

setup_tmpdir
TOOLKIT_A="$TMPDIR_TEST/toolkit-a"
TOOLKIT_B="$TMPDIR_TEST/toolkit-b"
build_fake_toolkit "$TOOLKIT_A"
build_fake_toolkit "$TOOLKIT_B"

mkdir "$TMPDIR_TEST/project" && cd "$TMPDIR_TEST/project"

# Activate against toolkit A (the "normal" mount).
SCIAGENT_TOOLKIT="$TOOLKIT_A" "$TOOLKIT_A/bin/sciagent" activate base >/dev/null
assert_symlink .claude/skills/s_a   "sanity: mounted against toolkit A"
assert_symlink .claude/agents/ag_a.md "sanity: agent mounted against toolkit A"
assert_file_exists .sciagent/manifest.json "sanity: manifest present after activate"
assert_grep 'BEGIN SCIAGENT:ROLES' AGENTS.md "sanity: ROLES block present after activate"

mounted_before=$(find .claude .agents -mindepth 2 -maxdepth 2 -type l 2>/dev/null | wc -l)
[[ "$mounted_before" -gt 0 ]] || { echo "FAIL [$_TEST_NAME] fixture: nothing mounted" >&2; exit 1; }

# Deactivate with SCIAGENT_TOOLKIT pointed at the UNRELATED toolkit B — the
# exact repro shape. Use toolkit A's own dispatcher binary (as a user
# normally would) but override the env var, since bin/sciagent honors an
# already-exported SCIAGENT_TOOLKIT over its own self-location.
out=$(SCIAGENT_TOOLKIT="$TOOLKIT_B" "$TOOLKIT_A/bin/sciagent" deactivate 2>&1)
rc=$?

if [[ "$rc" -ne 0 ]]; then
    echo "FAIL [$_TEST_NAME] deactivate against external toolkit should still succeed (rc=$rc)" >&2
    printf '%s\n' "$out" >&2
    exit 1
fi

# The core regression: no mount may be left behind.
left=$(find .claude .agents -mindepth 2 -maxdepth 2 -type l 2>/dev/null | wc -l)
if [[ "$left" -ne 0 ]]; then
    echo "FAIL [$_TEST_NAME] deactivate against external toolkit stranded $left of $mounted_before mounts" >&2
    printf '%s\n' "$out" >&2
    find .claude .agents -mindepth 2 -maxdepth 2 -type l 2>/dev/null >&2
    exit 1
fi

# Bookkeeping must be gone too (this already worked pre-fix — pin it stays true).
if [[ -e .sciagent/manifest.json ]]; then
    echo "FAIL [$_TEST_NAME] manifest still present after deactivate" >&2
    exit 1
fi
if grep -q 'BEGIN SCIAGENT:ROLES' AGENTS.md 2>/dev/null; then
    echo "FAIL [$_TEST_NAME] ROLES block still present after deactivate" >&2
    exit 1
fi
if grep -q 'BEGIN SCIAGENT:CRAFT' AGENTS.md 2>/dev/null; then
    echo "FAIL [$_TEST_NAME] CRAFT block still present after deactivate" >&2
    exit 1
fi

# `status` (against the normal toolkit A) must now honestly report nothing
# active — not a lie, because there genuinely is nothing left.
status_out=$(SCIAGENT_TOOLKIT="$TOOLKIT_A" "$TOOLKIT_A/bin/sciagent" status 2>&1)
if ! printf '%s\n' "$status_out" | grep -qi 'no active\|not active\|no role'; then
    echo "FAIL [$_TEST_NAME] status after cross-toolkit deactivate should report no active stack" >&2
    printf '%s\n' "$status_out" >&2
    exit 1
fi

pass
