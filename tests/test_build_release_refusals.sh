#!/usr/bin/env bash
# tests/test_build_release_refusals.sh
#
# Everything build-release.sh must REFUSE, and the invariant that a refusal
# never leaves a file behind — "a release that fails its own suite must not
# exist as a file" (10_packaging_contracts.md §1).
#
#   A. no ref                → refused (an implicit "whatever is checked out"
#                              is exactly what the contract forbids)
#   B. unresolvable ref      → refused
#   C. dirty working tree    → refused, and the offending path is named
#   D. `sciagent validate` fails on the exported tree → no artifact
#   E. `tests/run-all.sh` fails on the exported tree  → no artifact
#
# D and E are why the release tests build FIXTURE repos: the fixture's gate
# stubs carry a committed exit code, so both directions of the gate are
# exercised in milliseconds and this script never re-enters the real suite.
set -u
. "$(dirname "$0")/_lib.sh"
. "$(dirname "$0")/_release_lib.sh"

setup_tmpdir
release_git_env

# assert_no_artifacts <dir>
assert_no_artifacts() {
    local d="$1" msg="$2" f
    for f in "$d"/*.tar.gz "$d"/*.sha256 "$d"/*.metadata.json; do
        if [[ -e "$f" ]]; then
            echo "FAIL [$_TEST_NAME] $msg — but $f exists" >&2
            ls -la "$d" >&2
            exit 1
        fi
    done
}

# ---------------------------------------------------------------- A + B + C
REPO="$TMPDIR_TEST/ok"
mkdir -p "$REPO"
fixture_repo "$REPO"

# A. no ref at all
out=$( cd "$REPO" && ./scripts/build-release.sh --out "$TMPDIR_TEST/outA" 2>&1 )
rc=$?
[[ $rc -ne 0 ]] || { echo "FAIL [$_TEST_NAME] building with no ref succeeded" >&2; exit 1; }
printf '%s\n' "$out" | grep -qi "ref is required" \
    || { echo "FAIL [$_TEST_NAME] no-ref refusal does not explain itself" >&2; printf '%s\n' "$out" >&2; exit 1; }
[[ -d "$TMPDIR_TEST/outA" ]] && assert_no_artifacts "$TMPDIR_TEST/outA" "no-ref build produced files"

# B. unresolvable ref
out=$( cd "$REPO" && ./scripts/build-release.sh no-such-ref --out "$TMPDIR_TEST/outB" 2>&1 )
rc=$?
[[ $rc -ne 0 ]] || { echo "FAIL [$_TEST_NAME] building an unresolvable ref succeeded" >&2; exit 1; }
[[ -d "$TMPDIR_TEST/outB" ]] && assert_no_artifacts "$TMPDIR_TEST/outB" "bad-ref build produced files"

# C. dirty tree — an untracked, non-ignored file is enough
echo "scratch" > "$REPO/UNTRACKED.md"
out=$( cd "$REPO" && ./scripts/build-release.sh HEAD --out "$TMPDIR_TEST/outC" 2>&1 )
rc=$?
[[ $rc -ne 0 ]] || { echo "FAIL [$_TEST_NAME] building from a dirty tree succeeded" >&2; exit 1; }
printf '%s\n' "$out" | grep -qi "dirty working tree" \
    || { echo "FAIL [$_TEST_NAME] dirty-tree refusal does not say so" >&2; printf '%s\n' "$out" >&2; exit 1; }
printf '%s\n' "$out" | grep -q "UNTRACKED.md" \
    || { echo "FAIL [$_TEST_NAME] dirty-tree refusal does not name the offending path" >&2; printf '%s\n' "$out" >&2; exit 1; }
[[ -d "$TMPDIR_TEST/outC" ]] && assert_no_artifacts "$TMPDIR_TEST/outC" "dirty-tree build produced files"
rm -f "$REPO/UNTRACKED.md"

# Sanity: with the same repo clean, the build succeeds. Without this, C could
# pass for the wrong reason (e.g. the fixture never builds at all).
release_build "$REPO" "$TMPDIR_TEST/outOK" HEAD >/dev/null 2>&1 \
    || { echo "FAIL [$_TEST_NAME] clean fixture failed to build (control case)" >&2; exit 1; }
artifact_in "$TMPDIR_TEST/outOK" >/dev/null \
    || { echo "FAIL [$_TEST_NAME] control build produced no artifact" >&2; exit 1; }

# ---------------------------------------------------------------------- D
REPO_V="$TMPDIR_TEST/badvalidate"
mkdir -p "$REPO_V"
fixture_repo "$REPO_V" 1 0            # validate exits 1, tests exit 0
out=$( cd "$REPO_V" && ./scripts/build-release.sh HEAD --out "$TMPDIR_TEST/outD" 2>&1 )
rc=$?
[[ $rc -ne 0 ]] || { echo "FAIL [$_TEST_NAME] build succeeded despite failing validate" >&2; exit 1; }
printf '%s\n' "$out" | grep -qi "validate failed" \
    || { echo "FAIL [$_TEST_NAME] validate-gate failure not reported" >&2; printf '%s\n' "$out" >&2; exit 1; }
[[ -d "$TMPDIR_TEST/outD" ]] && assert_no_artifacts "$TMPDIR_TEST/outD" "failing validate still produced a release"

# ---------------------------------------------------------------------- E
REPO_T="$TMPDIR_TEST/badtests"
mkdir -p "$REPO_T"
fixture_repo "$REPO_T" 0 1            # validate exits 0, tests exit 1
out=$( cd "$REPO_T" && ./scripts/build-release.sh HEAD --out "$TMPDIR_TEST/outE" 2>&1 )
rc=$?
[[ $rc -ne 0 ]] || { echo "FAIL [$_TEST_NAME] build succeeded despite a failing test suite" >&2; exit 1; }
printf '%s\n' "$out" | grep -qi "test suite failed" \
    || { echo "FAIL [$_TEST_NAME] test-gate failure not reported" >&2; printf '%s\n' "$out" >&2; exit 1; }
[[ -d "$TMPDIR_TEST/outE" ]] && assert_no_artifacts "$TMPDIR_TEST/outE" "failing suite still produced a release"

pass
