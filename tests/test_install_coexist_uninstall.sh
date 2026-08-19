#!/usr/bin/env bash
# tests/test_install_coexist_uninstall.sh
#
# "Two versions coexist by construction" and "writes an installation receipt
# sufficient for an exact uninstall" (10_packaging_contracts.md §2).
#
# Coexistence:
#   two commits → two archives → both installed → both trees and both receipts
#   present, and the executable link points at the one installed last.
#
# Exact uninstall — the ownership discipline the toolkit uses everywhere else:
#   - uninstalling the newer version removes ITS tree, ITS receipt and the link
#     it owns, and leaves the older version bit-for-bit untouched;
#   - uninstalling both leaves no files or symlinks under the prefix at all;
#   - a receipt that does not exist is refused rather than guessed at;
#   - a bin link that has been re-pointed elsewhere is NOT removed: it is
#     reported and left alone (exit 3), because the receipt can no longer prove
#     ownership of it.
set -u
. "$(dirname "$0")/_lib.sh"
. "$(dirname "$0")/_release_lib.sh"

INSTALL="$TOOLKIT_ROOT/install.sh"

setup_tmpdir
release_git_env

REPO="$TMPDIR_TEST/repo"
mkdir -p "$REPO"
fixture_repo "$REPO"
SHA1=$(git -C "$REPO" rev-parse HEAD)

echo "second release content" > "$REPO/docs/second.md"
fixture_commit "$REPO" "second release commit" "2026-02-02T00:00:00+0000"
SHA2=$(git -C "$REPO" rev-parse HEAD)
[[ "$SHA1" != "$SHA2" ]] || { echo "FAIL [$_TEST_NAME] fixture produced one commit, not two" >&2; exit 1; }

release_build "$REPO" "$TMPDIR_TEST/out1" "$SHA1" --version 0.1.0 >/dev/null 2>&1 \
    || { echo "FAIL [$_TEST_NAME] build of commit 1 failed" >&2; exit 1; }
release_build "$REPO" "$TMPDIR_TEST/out2" "$SHA2" --version 0.2.0 >/dev/null 2>&1 \
    || { echo "FAIL [$_TEST_NAME] build of commit 2 failed" >&2; exit 1; }
A1=$(artifact_in "$TMPDIR_TEST/out1"); A2=$(artifact_in "$TMPDIR_TEST/out2")

P="$TMPDIR_TEST/prefix"
"$INSTALL" --archive "$A1" --checksum "$A1.sha256" --prefix "$P" >/dev/null 2>&1 \
    || { echo "FAIL [$_TEST_NAME] install of version 1 failed" >&2; exit 1; }
"$INSTALL" --archive "$A2" --checksum "$A2.sha256" --prefix "$P" >/dev/null 2>&1 \
    || { echo "FAIL [$_TEST_NAME] install of version 2 failed" >&2; exit 1; }

V="$P/share/scio/versions"
R="$P/share/scio/receipts"

# --- coexistence -----------------------------------------------------------
assert_file_exists "$V/$SHA1/bin/scio" "version 1 tree missing after installing version 2"
assert_file_exists "$V/$SHA2/bin/scio" "version 2 tree missing"
assert_file_exists "$R/$SHA1.json" "version 1 receipt missing"
assert_file_exists "$R/$SHA2.json" "version 2 receipt missing"
# Content-addressing is real: the two trees differ.
[[ -e "$V/$SHA2/docs/second.md" && ! -e "$V/$SHA1/docs/second.md" ]] \
    || { echo "FAIL [$_TEST_NAME] the two version trees are not distinct" >&2; exit 1; }

assert_symlink "$P/bin/scio"
assert_eq "$(readlink "$P/bin/scio")" "../share/scio/versions/$SHA2/bin/scio" \
    "the executable must point at the version installed last"

# Only the executable escapes the version directory.
binlist=$(ls -A "$P/bin")
assert_eq "$binlist" "scio" "nothing but the executable may be linked into <prefix>/bin"

list_out=$("$INSTALL" --list --prefix "$P" 2>&1)
printf '%s\n' "$list_out" | grep -q "$SHA1" || { echo "FAIL [$_TEST_NAME] --list omits version 1" >&2; printf '%s\n' "$list_out" >&2; exit 1; }
printf '%s\n' "$list_out" | grep -q "$SHA2" || { echo "FAIL [$_TEST_NAME] --list omits version 2" >&2; printf '%s\n' "$list_out" >&2; exit 1; }

# --- pruning the version that does NOT own the link -----------------------
# The ordinary housekeeping case: two versions installed, remove the one not in
# use. The link belongs to the other version's receipt, so it stays and the exit
# status stays clean — a warning here would break `--uninstall <old> && ...` and
# would read as damage where there is none.
out=$("$INSTALL" --uninstall "$SHA1" --prefix "$P" 2>&1)
rc=$?
assert_eq "$rc" "0" "pruning a non-current version must exit 0, not report an anomaly"
[[ -e "$V/$SHA1" ]]      && { echo "FAIL [$_TEST_NAME] the pruned version's tree survived" >&2; exit 1; }
[[ -e "$R/$SHA1.json" ]] && { echo "FAIL [$_TEST_NAME] the pruned version's receipt survived" >&2; exit 1; }
assert_eq "$(readlink "$P/bin/scio")" "../share/scio/versions/$SHA2/bin/scio" \
    "pruning one version must not disturb the link the other owns"
assert_file_exists "$V/$SHA2/bin/scio" "the in-use version was collateral damage"

# Restore both, in the original order, so the link points at SHA2 again.
"$INSTALL" --archive "$A1" --checksum "$A1.sha256" --prefix "$P" >/dev/null 2>&1 \
    || { echo "FAIL [$_TEST_NAME] reinstall of version 1 failed" >&2; exit 1; }
"$INSTALL" --archive "$A2" --checksum "$A2.sha256" --prefix "$P" >/dev/null 2>&1 \
    || { echo "FAIL [$_TEST_NAME] reinstall of version 2 failed" >&2; exit 1; }
assert_eq "$(readlink "$P/bin/scio")" "../share/scio/versions/$SHA2/bin/scio" \
    "the restored fixture must match the state the next section assumes"

# --- receipt-driven exact uninstall ---------------------------------------
snap1_before=$(tree_snapshot "$V/$SHA1")

out=$("$INSTALL" --uninstall "$SHA2" --prefix "$P" 2>&1)
rc=$?
assert_eq "$rc" "0" "uninstall of the linked version should succeed cleanly"
[[ -e "$V/$SHA2" ]]      && { echo "FAIL [$_TEST_NAME] version 2 tree survived uninstall" >&2; exit 1; }
[[ -e "$R/$SHA2.json" ]] && { echo "FAIL [$_TEST_NAME] version 2 receipt survived uninstall" >&2; exit 1; }
[[ -L "$P/bin/scio" ]] && { echo "FAIL [$_TEST_NAME] the link this version owned survived uninstall" >&2; exit 1; }

snap1_after=$(tree_snapshot "$V/$SHA1")
assert_eq "$snap1_after" "$snap1_before" "uninstalling one version must not touch the other"
assert_file_exists "$R/$SHA1.json" "the other version's receipt was collateral damage"

# --- a missing receipt is refused, not guessed ----------------------------
out=$("$INSTALL" --uninstall "deadbeefdeadbeefdeadbeefdeadbeefdeadbeef" --prefix "$P" 2>&1)
rc=$?
[[ $rc -ne 0 ]] || { echo "FAIL [$_TEST_NAME] uninstalling an unknown sha reported success" >&2; exit 1; }
printf '%s\n' "$out" | grep -qi "no receipt" \
    || { echo "FAIL [$_TEST_NAME] missing-receipt refusal does not say why" >&2; printf '%s\n' "$out" >&2; exit 1; }

# --- a re-pointed link is left alone (cannot prove ownership) -------------
# Reinstate version 1's link, then repoint it at something else by hand.
"$INSTALL" --archive "$A1" --checksum "$A1.sha256" --prefix "$P" >/dev/null 2>&1
ln -sfn "/somewhere/else/scio" "$P/bin/scio"
out=$("$INSTALL" --uninstall "$SHA1" --prefix "$P" 2>&1)
rc=$?
assert_eq "$rc" "3" "a link that no longer matches the receipt must exit 3, not silently delete"
assert_symlink "$P/bin/scio" "a hand-repointed link must be left alone"
assert_eq "$(readlink "$P/bin/scio")" "/somewhere/else/scio" "the foreign link target was altered"
[[ -e "$V/$SHA1" ]] && { echo "FAIL [$_TEST_NAME] the version tree should still be removed" >&2; exit 1; }

# --- nothing of ours remains ----------------------------------------------
rm -f "$P/bin/scio"          # the foreign link we planted, not ours
ours=$(find "$V" "$R" \( -type f -o -type l \) 2>/dev/null)
if [[ -n "$ours" ]]; then
    echo "FAIL [$_TEST_NAME] uninstall was not exact; residue:" >&2
    printf '%s\n' "$ours" >&2
    exit 1
fi

pass
