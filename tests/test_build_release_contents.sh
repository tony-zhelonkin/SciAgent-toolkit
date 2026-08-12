#!/usr/bin/env bash
# tests/test_build_release_contents.sh
#
# `git archive` is mandatory rather than stylistic (00_INDEX.md §5: the working
# checkout is 180 MB against 5.6 MB tracked, 171 MB of it one skill's .venv).
# This test asserts the property that mandate exists for, on a fixture that
# reproduces the shape: a 2 MiB gitignored `.venv` next to a few KB of tracked
# content.
#
# Asserted:
#   - no ignored path (.venv) appears in the artifact, at all
#   - the artifact is orders of magnitude smaller than the working tree, so a
#     future `tar -czf .`-style regression cannot pass by accident
#   - the archive is prefixed (no tarbomb) and carries `.scio-release.json`
#   - that embedded metadata carries the FULL 40-char sha, equal to the ref
#   - the filename carries the abbreviated sha (the §1-vs-ADR-D7 resolution)
#   - the sidecar metadata's sha256 field describes the actual artifact
set -u
. "$(dirname "$0")/_lib.sh"
. "$(dirname "$0")/_release_lib.sh"

setup_tmpdir
release_git_env

REPO="$TMPDIR_TEST/repo"
mkdir -p "$REPO"
fixture_repo "$REPO"

full_sha=$(git -C "$REPO" rev-parse HEAD)
short_sha=$(git -C "$REPO" rev-parse --short=12 HEAD)

release_build "$REPO" "$TMPDIR_TEST/out" HEAD >/dev/null 2>&1 \
    || { echo "FAIL [$_TEST_NAME] build failed" >&2; exit 1; }

art=$(artifact_in "$TMPDIR_TEST/out") || { echo "FAIL [$_TEST_NAME] no artifact" >&2; exit 1; }
members=$(tar tzf "$art")

# --- 1. no ignored content ------------------------------------------------
if printf '%s\n' "$members" | grep -q '\.venv'; then
    echo "FAIL [$_TEST_NAME] artifact contains ignored .venv content" >&2
    printf '%s\n' "$members" | grep '\.venv' >&2
    exit 1
fi

# --- 2. size is the right order of magnitude ------------------------------
tree_kb=$(du -sk "$REPO" | awk '{print $1}')
art_bytes=$(wc -c < "$art" | tr -d ' ')
if (( tree_kb < 2048 )); then
    echo "FAIL [$_TEST_NAME] fixture is not fat enough to make the size claim meaningful (${tree_kb}K)" >&2
    exit 1
fi
if (( art_bytes > 65536 )); then
    echo "FAIL [$_TEST_NAME] artifact is $art_bytes bytes; tracked fixture content is a few KB." >&2
    echo "  A tar of the working tree would look like this. git archive would not." >&2
    exit 1
fi

# --- 3. prefixed, self-describing -----------------------------------------
stem=$(basename "${art%.tar.gz}")
top=$(printf '%s\n' "$members" | head -1 | cut -d/ -f1)
assert_eq "$top" "$stem" "archive must be prefixed with its own stem (no tarbomb)"

printf '%s\n' "$members" | grep -qx "$stem/.scio-release.json" \
    || { echo "FAIL [$_TEST_NAME] archive does not embed .scio-release.json" >&2
         printf '%s\n' "$members" >&2; exit 1; }

embedded=$(tar xzOf "$art" "$stem/.scio-release.json")
emb_commit=$(printf '%s\n' "$embedded" | grep -o '"commit"[[:space:]]*:[[:space:]]*"[^"]*"' | head -1 | sed 's/.*"\([0-9a-f]*\)"$/\1/')
assert_eq "$emb_commit" "$full_sha" "embedded metadata must carry the FULL commit sha"
assert_eq "${#emb_commit}" "40" "embedded commit sha must be 40 chars, never abbreviated"

# --- 4. filename carries the abbreviation ---------------------------------
case "$stem" in
    scio-*-"$short_sha") ;;
    *) echo "FAIL [$_TEST_NAME] artifact stem '$stem' is not scio-<version>-<short-sha>" >&2; exit 1 ;;
esac

# --- 5. sidecar metadata describes this artifact --------------------------
meta="$TMPDIR_TEST/out/$stem.metadata.json"
assert_file_exists "$meta"
meta_sha=$(grep -o '"sha256"[[:space:]]*:[[:space:]]*"[^"]*"' "$meta" | sed 's/.*"\([0-9a-f]*\)"$/\1/')
assert_eq "$meta_sha" "$(sha256_hex "$art")" "metadata sha256 must match the artifact"
assert_grep "\"commit\": \"$full_sha\"" "$meta" "metadata must record the full sha"
assert_grep '"size_bytes"' "$meta" "metadata must record the artifact size"

pass
