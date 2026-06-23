#!/usr/bin/env bash
# tests/test_activate_renders_craft.sh
# `sciagent activate` renders the SCIAGENT:CRAFT block from <toolkit>/craft.yaml
# alongside SCIAGENT:ROLES; re-activate is idempotent; deactivate removes both.
set -u
. "$(dirname "$0")/_lib.sh"

setup_tmpdir
FAKE="$TMPDIR_TEST/fake-toolkit"
build_fake_toolkit "$FAKE"

# Provision a craft.yaml in the fake toolkit (build_fake_toolkit ships none, so
# the existing activate tests exercise the no-craft.yaml no-op path).
cat > "$FAKE/craft.yaml" <<'EOF'
version: 1

floors:
  figure_base_size: 16

body: |
  # Craft standards
  - Figures: base >= {{figure_base_size}}pt; style via the project theme.
  - Results: under 03_results/<stage>/{figures,tables}/.
EOF

export SCIAGENT_TOOLKIT="$FAKE"
SCIAGENT="$FAKE/bin/sciagent"

mkdir project && cd project
cat > AGENTS.md <<'EOF'
# Project AGENTS.md

User-owned content.
EOF
cp AGENTS.md AGENTS.md.orig

"$SCIAGENT" activate base >/dev/null

# Both managed blocks coexist.
assert_grep 'BEGIN SCIAGENT:ROLES' AGENTS.md "ROLES block present"
assert_grep 'BEGIN SCIAGENT:CRAFT' AGENTS.md "CRAFT block present"
assert_grep 'Craft standards'      AGENTS.md "CRAFT body rendered"
# Token substitution happened (16, not the literal token).
assert_grep 'base >= 16pt'         AGENTS.md "floor token substituted"
if grep -qF '{{figure_base_size}}' AGENTS.md; then
    echo "FAIL [$_TEST_NAME] unsubstituted {{token}} left in CRAFT block" >&2
    exit 1
fi

# Idempotent: re-activate leaves AGENTS.md byte-identical.
cp AGENTS.md AGENTS.md.first
"$SCIAGENT" activate base >/dev/null
assert_file_eq AGENTS.md AGENTS.md.first "AGENTS.md unchanged on re-activate (CRAFT idempotent)"

# Drift detection scopes to the CRAFT block independently of ROLES.
. "$FAKE/lib/sciagent/block.sh"
assert_exit 0 block_hash_check AGENTS.md CRAFT
assert_exit 0 block_hash_check AGENTS.md
sed -i 's/Craft standards/Craft STANDARDS/' AGENTS.md
assert_exit 3 block_hash_check AGENTS.md CRAFT
assert_exit 0 block_hash_check AGENTS.md

# Deactivate removes both blocks and restores the file byte-for-byte.
"$SCIAGENT" activate base >/dev/null      # re-render to clear the manual drift
"$SCIAGENT" deactivate >/dev/null
if grep -qF 'SCIAGENT:CRAFT' AGENTS.md; then
    echo "FAIL [$_TEST_NAME] CRAFT block survived deactivate" >&2
    exit 1
fi
assert_file_eq AGENTS.md AGENTS.md.orig "AGENTS.md restored byte-for-byte after deactivate"

pass
