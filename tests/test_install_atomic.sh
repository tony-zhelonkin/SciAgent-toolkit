#!/usr/bin/env bash
# tests/test_install_atomic.sh
#
# "An interrupted install leaves no half-tree" (10_packaging_contracts.md §2).
#
# The interruption is produced deterministically rather than by racing a signal:
# the install runs under `ulimit -f`, a per-process file-size cap set just above
# the small files and below one deliberately large tracked member. `tar` then
# dies PART-WAY THROUGH extraction — several files already written, one member
# killed mid-write — which is precisely the state the atomicity claim is about.
# A truncated archive would not test this: it fails before extraction starts.
#
# Asserted after the failure:
#   - the install exits non-zero and says extraction failed
#   - no version directory exists (not even an empty one)
#   - no receipt exists
#   - no bin/ symlink exists
#   - no staging directory is left behind — the partial tree is really gone,
#     not merely hidden somewhere under the prefix
#   - a subsequent NORMAL install of the same archive succeeds, i.e. the failed
#     attempt poisoned nothing
set -u
. "$(dirname "$0")/_lib.sh"
. "$(dirname "$0")/_release_lib.sh"

INSTALL="$TOOLKIT_ROOT/install.sh"

setup_tmpdir
release_git_env

REPO="$TMPDIR_TEST/repo"
mkdir -p "$REPO"
fixture_repo "$REPO"
# One tracked member far larger than the cap below. `zbig.bin` sorts last, so
# the small members extract first and the failure really is mid-extraction.
dd if=/dev/urandom of="$REPO/zbig.bin" bs=1024 count=1024 status=none
fixture_commit "$REPO" "add a large tracked member"

release_build "$REPO" "$TMPDIR_TEST/out" HEAD >/dev/null 2>&1 \
    || { echo "FAIL [$_TEST_NAME] fixture build failed" >&2; exit 1; }
ART=$(artifact_in "$TMPDIR_TEST/out")
SUM="$ART.sha256"
SHA=$(git -C "$REPO" rev-parse HEAD)

P="$TMPDIR_TEST/prefix"

# ulimit -f 8 → 4 KiB. Big enough for the metadata read and the stub files,
# far too small for the 1 MiB member. ulimit -c 0 keeps a SIGXFSZ from dropping
# a core file anywhere.
out=$( ulimit -c 0; ulimit -f 8; "$INSTALL" --archive "$ART" --checksum "$SUM" --prefix "$P" 2>&1 )
rc=$?

[[ $rc -ne 0 ]] || { echo "FAIL [$_TEST_NAME] interrupted install reported success" >&2; printf '%s\n' "$out" >&2; exit 1; }
printf '%s\n' "$out" | grep -qi "extraction failed" \
    || { echo "FAIL [$_TEST_NAME] interrupted install did not report a failed extraction" >&2
         printf '%s\n' "$out" >&2; exit 1; }

if [[ -e "$P/share/scio/versions/$SHA" ]]; then
    echo "FAIL [$_TEST_NAME] half-tree left behind at $P/share/scio/versions/$SHA" >&2
    find "$P/share/scio/versions/$SHA" | head -20 >&2
    exit 1
fi
if [[ -e "$P/share/scio/receipts/$SHA.json" ]]; then
    echo "FAIL [$_TEST_NAME] a receipt was written for an install that never completed" >&2
    exit 1
fi
if [[ -e "$P/bin/sciagent" || -L "$P/bin/sciagent" ]]; then
    echo "FAIL [$_TEST_NAME] the executable was linked despite a failed extraction" >&2
    exit 1
fi
leftovers=$(find "$P" -name 'install.*' -o -name '.staging' 2>/dev/null)
if [[ -n "$leftovers" ]]; then
    echo "FAIL [$_TEST_NAME] staging residue survived the failure:" >&2
    printf '%s\n' "$leftovers" >&2
    find "$P" >&2
    exit 1
fi
# Nothing but (possibly empty) skeleton directories may remain.
files=$(find "$P" \( -type f -o -type l \) 2>/dev/null)
if [[ -n "$files" ]]; then
    echo "FAIL [$_TEST_NAME] failed install left files behind:" >&2
    printf '%s\n' "$files" >&2
    exit 1
fi

# The prefix is not poisoned: a normal install now works.
"$INSTALL" --archive "$ART" --checksum "$SUM" --prefix "$P" >/dev/null 2>&1 \
    || { echo "FAIL [$_TEST_NAME] a normal install after the failed one did not succeed" >&2; exit 1; }
assert_file_exists "$P/share/scio/versions/$SHA/bin/sciagent" "retry did not install the tree"
assert_symlink "$P/bin/sciagent" "retry did not link the executable"

pass
