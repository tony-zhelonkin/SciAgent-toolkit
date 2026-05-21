#!/usr/bin/env bash
# Activate then deactivate (no arg) — manifest, symlinks, and block all gone.
set -u
. "$(dirname "$0")/_lib.sh"

setup_tmpdir
FAKE="$TMPDIR_TEST/fake-toolkit"
build_fake_toolkit "$FAKE"
export SCIAGENT_TOOLKIT="$FAKE"
SCIAGENT="$FAKE/bin/sciagent"

mkdir project && cd project
cat > AGENTS.md <<'EOF'
# Project AGENTS.md

User-owned content.
EOF
cp AGENTS.md AGENTS.md.orig

"$SCIAGENT" activate base reviewer >/dev/null
"$SCIAGENT" deactivate >/dev/null

[[ -e .sciagent/manifest.json ]] && { echo "FAIL: manifest still exists"; exit 1; }
[[ -e .claude/skills/s_a ]] && { echo "FAIL: claude symlink still exists"; exit 1; }
[[ -e .agents/skills/s_a ]] && { echo "FAIL: agents-mirror symlink still exists"; exit 1; }

assert_file_eq AGENTS.md AGENTS.md.orig "AGENTS.md restored byte-for-byte"
pass
