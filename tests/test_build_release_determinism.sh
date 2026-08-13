#!/usr/bin/env bash
# tests/test_build_release_determinism.sh
#
# The contract's own stated acceptance test (10_packaging_contracts.md §1,
# "Deterministic: a test builds the same ref twice and asserts identical
# checksums"), plus the two corollaries that make the claim useful:
#
#   1. Two builds of the same ref produce a byte-identical tarball — including
#      across different output directories and separated in time, so a leaked
#      build mtime would show up.
#   2. The sidecar metadata is identical too. It deliberately carries no
#      wall-clock build time, so the whole triple is reproducible, not just the
#      archive.
#   3. The emitted .sha256 actually verifies under `sha256sum -c`, i.e. it is a
#      usable checksum file and not just a hex string in a file.
#
# Hermetic: throwaway git fixture under the test's own tmpdir; no network; no
# $HOME access; no dependence on this repo being a git checkout.
set -u
. "$(dirname "$0")/_lib.sh"
. "$(dirname "$0")/_release_lib.sh"

setup_tmpdir
release_git_env

REPO="$TMPDIR_TEST/repo"
mkdir -p "$REPO"
fixture_repo "$REPO"

out1=$(release_build "$REPO" "$TMPDIR_TEST/out1" HEAD 2>/dev/null) \
    || { echo "FAIL [$_TEST_NAME] first build failed" >&2; exit 1; }
# A second apart, into a different directory: any embedded timestamp diverges.
sleep 2
out2=$(release_build "$REPO" "$TMPDIR_TEST/out2" HEAD 2>/dev/null) \
    || { echo "FAIL [$_TEST_NAME] second build failed" >&2; exit 1; }

a1=$(artifact_in "$TMPDIR_TEST/out1") || { echo "FAIL [$_TEST_NAME] no artifact in out1" >&2; exit 1; }
a2=$(artifact_in "$TMPDIR_TEST/out2") || { echo "FAIL [$_TEST_NAME] no artifact in out2" >&2; exit 1; }

assert_eq "$(basename "$a1")" "$(basename "$a2")" "same ref must yield the same filename"

h1=$(sha256_hex "$a1")
h2=$(sha256_hex "$a2")
assert_eq "$h1" "$h2" "same ref built twice must be byte-identical (gzip -n)"

# Exactly three files per build, no more.
n1=$(find "$TMPDIR_TEST/out1" -maxdepth 1 -type f | wc -l)
assert_eq "$n1" "3" "a build produces exactly three files (tarball, .sha256, metadata)"

# Metadata is reproducible too (no build timestamp inside it).
m1="$TMPDIR_TEST/out1/$(basename "${a1%.tar.gz}").metadata.json"
m2="$TMPDIR_TEST/out2/$(basename "${a2%.tar.gz}").metadata.json"
assert_file_exists "$m1" "sidecar metadata missing for build 1"
assert_file_exists "$m2" "sidecar metadata missing for build 2"
assert_file_eq "$m2" "$m1" "release metadata must be reproducible (no build timestamp)"

# The .sha256 sidecar must verify with the standard tool, from its own dir.
( cd "$TMPDIR_TEST/out1" && sha256sum -c "$(basename "$a1").sha256" >/dev/null 2>&1 ) \
    || { echo "FAIL [$_TEST_NAME] emitted .sha256 does not verify under sha256sum -c" >&2
         cat "$a1.sha256" >&2; exit 1; }

# Rebuilding into the SAME directory is also stable (overwrite path).
release_build "$REPO" "$TMPDIR_TEST/out1" HEAD >/dev/null 2>&1 \
    || { echo "FAIL [$_TEST_NAME] rebuild into an existing out dir failed" >&2; exit 1; }
assert_eq "$(sha256_hex "$a1")" "$h1" "rebuild in place must not change the bytes"

pass
