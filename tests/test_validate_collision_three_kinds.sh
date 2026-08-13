#!/usr/bin/env bash
# tests/test_validate_collision_three_kinds.sh
# A three-namespace overlap names skill, agent, and command in the warning.
set -u
. "$(dirname "$0")/_lib.sh"

setup_tmpdir
FAKE="$TMPDIR_TEST/fake-toolkit"
build_fake_toolkit "$FAKE"

# Plant a 3-namespace overlap on the name 'tri'. The fake toolkit ships
# skills/s_a, s_b, s_c — none named 'tri' — so we add all three from scratch.
mkdir -p "$FAKE/skills/tri"
cat > "$FAKE/skills/tri/SKILL.md" <<'EOF'
---
name: tri
description: Fixture skill colliding with an agent and a command of the same name.
---
tri as a skill
EOF
echo "agent tri"   > "$FAKE/agents/tri.md"
echo "command tri" > "$FAKE/commands/tri.md"

export SCIO_TOOLKIT="$FAKE"
SCIO="$FAKE/bin/scio"

set +e
stdout_out=$("$SCIO" lint --check toolkit 2>/tmp/_validate_err_$$)
rc=$?
stderr_out=$(cat /tmp/_validate_err_$$)
rm -f /tmp/_validate_err_$$
set -e

# Exit code must be 0 — cross-namespace collisions are soft-warn only.
if [[ "$rc" -ne 0 ]]; then
    echo "FAIL [$_TEST_NAME] toolkit lint exited $rc on a soft-warn collision (expected 0)" >&2
    printf '%s\n' "$stderr_out" >&2
    exit 1
fi

# Locate the warning line for 'tri'.
tri_warning=$(printf '%s\n' "$stderr_out" | grep "'tri'" || true)
if [[ -z "$tri_warning" ]]; then
    echo "FAIL [$_TEST_NAME] no warning for 'tri' in toolkit lint stderr" >&2
    printf '%s\n' "$stderr_out" >&2
    exit 1
fi

# All three kinds must appear in the warning. The natural-list join renders
# 3 kinds as "skill, agent, and command" in canonical order.
for kind in skill agent command; do
    if ! printf '%s\n' "$tri_warning" | grep -qw "$kind"; then
        echo "FAIL [$_TEST_NAME] kind '$kind' missing from 3-namespace collision warning" >&2
        echo "warning line: $tri_warning" >&2
        exit 1
    fi
done

phrase_segment="${tri_warning#*appears as }"
if [[ -z "$phrase_segment" ]]; then
    echo "FAIL [$_TEST_NAME] could not extract 'appears as ...' segment" >&2
    echo "warning line: $tri_warning" >&2
    exit 1
fi
if printf '%s' "$phrase_segment" | grep -q 'both '; then
    echo "FAIL [$_TEST_NAME] 3-kind warning still uses the 2-kind 'both ...' template" >&2
    echo "phrase segment: $phrase_segment" >&2
    exit 1
fi

pass
