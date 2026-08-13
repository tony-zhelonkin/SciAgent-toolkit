#!/usr/bin/env bash
# tests/test_install_verify_and_dryrun.sh
#
# Two contract clauses that share a fixture:
#
#   1. "Verifies the checksum before extracting anything." Tested at its
#      strongest: after a mismatch the PREFIX ITSELF must not exist. Asserting
#      "no version directory" would pass even if the installer had created the
#      tree and then cleaned up; asserting "nothing at all" cannot.
#   2. "Supports --dry-run" — and a dry run must write literally nothing, not
#      "nothing important". Both the empty-prefix case and the
#      already-populated-prefix case are snapshotted (paths + content hashes).
#
# Also covers the no-URL rule: a URL-shaped --archive is refused with an
# explanation, not fetched. There is no network in this test because there is
# no network in the installer.
set -u
. "$(dirname "$0")/_lib.sh"
. "$(dirname "$0")/_release_lib.sh"

INSTALL="$TOOLKIT_ROOT/install.sh"

setup_tmpdir
release_git_env

REPO="$TMPDIR_TEST/repo"
mkdir -p "$REPO"
fixture_repo "$REPO"
release_build "$REPO" "$TMPDIR_TEST/out" HEAD >/dev/null 2>&1 \
    || { echo "FAIL [$_TEST_NAME] fixture build failed" >&2; exit 1; }
ART=$(artifact_in "$TMPDIR_TEST/out")
SUM="$ART.sha256"

# ---------------------------------------------------------------- corruption
cp "$ART" "$TMPDIR_TEST/corrupt.tar.gz"
printf 'corruption' >> "$TMPDIR_TEST/corrupt.tar.gz"

P_BAD="$TMPDIR_TEST/prefix-bad"
out=$("$INSTALL" --archive "$TMPDIR_TEST/corrupt.tar.gz" --checksum "$SUM" --prefix "$P_BAD" 2>&1)
rc=$?
[[ $rc -ne 0 ]] || { echo "FAIL [$_TEST_NAME] corrupt archive installed successfully" >&2; exit 1; }
printf '%s\n' "$out" | grep -qi "checksum mismatch" \
    || { echo "FAIL [$_TEST_NAME] corruption not reported as a checksum mismatch" >&2; printf '%s\n' "$out" >&2; exit 1; }
if [[ -e "$P_BAD" ]]; then
    echo "FAIL [$_TEST_NAME] verification must precede every write; $P_BAD was created" >&2
    find "$P_BAD" >&2
    exit 1
fi

# A truncated-but-self-consistent archive (checksum recomputed) must also be
# refused — verification passing does not mean the payload is installable.
gunzip -c "$ART" > "$TMPDIR_TEST/full.tar"
head -c 3000 "$TMPDIR_TEST/full.tar" > "$TMPDIR_TEST/trunc.tar"
gzip -n -c "$TMPDIR_TEST/trunc.tar" > "$TMPDIR_TEST/trunc.tar.gz"
( cd "$TMPDIR_TEST" && sha256sum trunc.tar.gz > trunc.tar.gz.sha256 )
P_TRUNC="$TMPDIR_TEST/prefix-trunc"
out=$("$INSTALL" --archive "$TMPDIR_TEST/trunc.tar.gz" --checksum "$TMPDIR_TEST/trunc.tar.gz.sha256" --prefix "$P_TRUNC" 2>&1)
rc=$?
[[ $rc -ne 0 ]] || { echo "FAIL [$_TEST_NAME] truncated archive installed successfully" >&2; exit 1; }
if [[ -d "$P_TRUNC/share/scio/versions" ]] && [[ -n "$(ls -A "$P_TRUNC/share/scio/versions" 2>/dev/null)" ]]; then
    echo "FAIL [$_TEST_NAME] truncated archive left a version tree" >&2
    find "$P_TRUNC" >&2
    exit 1
fi

# ------------------------------------------------------------------ no URLs
out=$("$INSTALL" --archive "https://example.invalid/scio.tar.gz" --checksum "$SUM" --prefix "$TMPDIR_TEST/p-url" 2>&1)
rc=$?
[[ $rc -ne 0 ]] || { echo "FAIL [$_TEST_NAME] a URL was accepted as --archive" >&2; exit 1; }
printf '%s\n' "$out" | grep -qi "no network code path" \
    || { echo "FAIL [$_TEST_NAME] URL refusal does not state the design rule" >&2; printf '%s\n' "$out" >&2; exit 1; }
[[ -e "$TMPDIR_TEST/p-url" ]] && { echo "FAIL [$_TEST_NAME] URL refusal created a prefix" >&2; exit 1; }

# ----------------------------------------------------------------- dry-run
# (a) into a prefix that does not exist yet
P_DRY="$TMPDIR_TEST/prefix-dry"
before=$(tree_snapshot "$P_DRY")
out=$("$INSTALL" --archive "$ART" --checksum "$SUM" --prefix "$P_DRY" --dry-run 2>&1)
rc=$?
after=$(tree_snapshot "$P_DRY")
assert_eq "$rc" "0" "dry-run must succeed"
assert_eq "$after" "$before" "dry-run must not create the prefix"
printf '%s\n' "$out" | grep -q '\[dry-run\]' \
    || { echo "FAIL [$_TEST_NAME] dry-run output has no [dry-run] marker" >&2; printf '%s\n' "$out" >&2; exit 1; }
printf '%s\n' "$out" | grep -q "would link" \
    || { echo "FAIL [$_TEST_NAME] dry-run does not announce the symlink it would create" >&2; printf '%s\n' "$out" >&2; exit 1; }
full_sha=$(git -C "$REPO" rev-parse HEAD)
printf '%s\n' "$out" | grep -q "$full_sha" \
    || { echo "FAIL [$_TEST_NAME] dry-run does not name the content-addressed target dir" >&2; printf '%s\n' "$out" >&2; exit 1; }

# (b) into a prefix that already has unrelated content — nothing may change
P_DRY2="$TMPDIR_TEST/prefix-dry2"
mkdir -p "$P_DRY2/bin" "$P_DRY2/share/scio/versions"
echo "someone else's program" > "$P_DRY2/bin/other-tool"
before2=$(tree_snapshot "$P_DRY2")
"$INSTALL" --archive "$ART" --checksum "$SUM" --prefix "$P_DRY2" --dry-run >/dev/null 2>&1 \
    || { echo "FAIL [$_TEST_NAME] dry-run into a populated prefix failed" >&2; exit 1; }
after2=$(tree_snapshot "$P_DRY2")
assert_eq "$after2" "$before2" "dry-run must not touch an existing prefix"

pass
