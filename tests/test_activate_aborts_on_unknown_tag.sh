#!/usr/bin/env bash
# tests/test_activate_aborts_on_unknown_tag.sh — a skill carrying an unknown
# tag causes `sciagent activate` to exit 1 before any filesystem mutation.
#
# Exercises the pre-mutation validate path (cmd_validate --quiet called at the
# start of cmd_activate). This is distinct from test_skill_requires_cycle.sh
# (which exercises the Phase-B resolver abort) and
# test_validate_output_style_drift.sh (which exercises the output_style check):
# unknown-tag failures are detected by cmd_validate before Phase B runs.
#
# Per kickoff.md §9 ADR-007 / activate pre-mutation contract (2026-05-24):
# tag-vocab failure is a hard error; no .claude/, .agents/, or manifest appear.
set -u
. "$(dirname "$0")/_lib.sh"

setup_tmpdir
FAKE="$TMPDIR_TEST/fake-toolkit"
build_fake_toolkit "$FAKE"
export SCIAGENT_TOOLKIT="$FAKE"
SCIAGENT="$FAKE/bin/sciagent"

# Plant a skill with a tag that does not exist in the fake toolkit's tags.yaml.
# build_fake_toolkit seeds only "tooling"; "not-a-real-tag" is unknown.
cat > "$FAKE/skills/s_a/SKILL.md" <<'EOF'
---
name: s_a
description: unknown-tag abort fixture
metadata:
  scope: implementation
  requires: []
  complementary-skills: []
  contraindications: []
  tags:
    - not-a-real-tag
---
skill s_a body
EOF

mkdir project && cd project

set +e
out=$("$SCIAGENT" activate base 2>&1)
rc=$?
set -e

if [[ "$rc" -eq 0 ]]; then
    echo "FAIL [$_TEST_NAME] activate returned 0 on unknown tag; expected 1" >&2
    echo "--- output ---" >&2
    printf '%s\n' "$out" >&2
    exit 1
fi

# Pre-mutation contract: nothing landed on the filesystem.
if [[ -e .sciagent/manifest.json ]]; then
    echo "FAIL [$_TEST_NAME] manifest.json written despite unknown-tag failure" >&2
    exit 1
fi
if [[ -e .claude ]]; then
    echo "FAIL [$_TEST_NAME] .claude/ created despite unknown-tag failure" >&2
    exit 1
fi
if [[ -e .agents ]]; then
    echo "FAIL [$_TEST_NAME] .agents/ created despite unknown-tag failure" >&2
    exit 1
fi

# Error output must name the offending tag and the aborting verb.
if ! printf '%s\n' "$out" | grep -q 'not-a-real-tag'; then
    echo "FAIL [$_TEST_NAME] error output does not name offending tag 'not-a-real-tag'" >&2
    echo "--- output ---" >&2
    printf '%s\n' "$out" >&2
    exit 1
fi

pass
