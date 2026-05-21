#!/usr/bin/env bash
# A role that requests an output_style with no matching frontmatter name
# in system-prompts/ must fail cleanly at activate-time: non-zero exit,
# actionable stderr, and no filesystem mutation (no manifest, no
# symlinks, no AGENTS.md block).
set -u
. "$(dirname "$0")/_lib.sh"

setup_tmpdir
FAKE="$TMPDIR_TEST/fake-toolkit"
build_fake_toolkit "$FAKE"
export SCIAGENT_TOOLKIT="$FAKE"
SCIAGENT="$FAKE/bin/sciagent"

# Add a role whose output_style does not exist anywhere in system-prompts/.
cat > "$FAKE/roles/with_ghost_style.yaml" <<EOF
name: with_ghost_style
description: fixture role requesting a non-existent style
output_style: ghost-style
skills:
  - s_a
EOF

mkdir project && cd project
cat > AGENTS.md <<'EOF'
# Project AGENTS.md

User-owned content.
EOF
cp AGENTS.md AGENTS.md.orig

# Activate should fail with non-zero exit.
set +e
out=$("$SCIAGENT" activate with_ghost_style 2>&1)
rc=$?
set -e

if [[ "$rc" -eq 0 ]]; then
    echo "FAIL [$_TEST_NAME] activate with drifted output_style returned 0; expected non-zero" >&2
    echo "--- stdout/stderr ---" >&2
    printf '%s\n' "$out" >&2
    exit 1
fi

# Stderr must name the requested style so the user can fix the role.
if ! printf '%s\n' "$out" | grep -q 'ghost-style'; then
    echo "FAIL [$_TEST_NAME] error message does not mention the requested style 'ghost-style'" >&2
    echo "--- stderr ---" >&2
    printf '%s\n' "$out" >&2
    exit 1
fi

# Stderr must list available styles so the user knows the options.
if ! printf '%s\n' "$out" | grep -q 'fixture-style'; then
    echo "FAIL [$_TEST_NAME] error message does not enumerate available styles (expected 'fixture-style')" >&2
    echo "--- stderr ---" >&2
    printf '%s\n' "$out" >&2
    exit 1
fi

# No filesystem mutation: no manifest, no symlinks, no block.
[[ -e .sciagent/manifest.json ]] && { echo "FAIL [$_TEST_NAME] manifest created despite validation failure"; exit 1; }
[[ -e .claude ]] && { echo "FAIL [$_TEST_NAME] .claude/ created despite validation failure"; exit 1; }
[[ -e .agents ]] && { echo "FAIL [$_TEST_NAME] .agents/ created despite validation failure"; exit 1; }

assert_file_eq AGENTS.md AGENTS.md.orig "AGENTS.md untouched on validation failure"
pass
