#!/usr/bin/env bash
# Re-activating the same stack should leave state identical (idempotent).
set -u
. "$(dirname "$0")/_lib.sh"

setup_tmpdir
FAKE="$TMPDIR_TEST/fake-toolkit"
build_fake_toolkit "$FAKE"
export SCIAGENT_TOOLKIT="$FAKE"
SCIAGENT="$FAKE/bin/sciagent"

mkdir project && cd project
"$SCIAGENT" activate base reviewer >/dev/null
cp AGENTS.md AGENTS.md.first
cp .sciagent/manifest.json manifest.json.first

"$SCIAGENT" activate base reviewer >/dev/null
assert_file_eq AGENTS.md AGENTS.md.first "AGENTS.md unchanged on re-activate"
# Manifest block-hash should match; symlink list identical.
diff -u manifest.json.first .sciagent/manifest.json >/dev/null \
    || { echo "manifest changed unexpectedly:"; diff -u manifest.json.first .sciagent/manifest.json; exit 1; }
pass
