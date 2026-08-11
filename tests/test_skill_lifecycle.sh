#!/usr/bin/env bash
# tests/test_skill_lifecycle.sh — lightweight skill lifecycle (_attic).
#
# Asserts the three lifecycle invariants, all against a synthetic toolkit:
#   (a) _attic skills are NOT in the active enumeration and NOT resolvable by
#       activate (no symlink is created, activate still succeeds);
#   (b) `list skills` lists _attic skills in a separate Attic section and does
#       not mix them in with the active ones;
#   (c) validate ignores _attic (no choke on the retired skill).
#
# Retirement is now the _attic move alone: the frontmatter `status:` field went
# with the rest of the metadata block in Phase 4.
#
# The lifecycle is a soft convention: no hooks, no fail-closed checks. See
# docs/skill-lifecycle.md.

set -u
. "$(dirname "$0")/_lib.sh"

setup_tmpdir
FAKE="$TMPDIR_TEST/fake-toolkit"
build_fake_toolkit "$FAKE"
export SCIAGENT_TOOLKIT="$FAKE"
SCIAGENT="$FAKE/bin/sciagent"

# --- Fixtures ----------------------------------------------------------------
# Plant a retired skill in the attic.
mkdir -p "$FAKE/skills/_attic/retired_skill"
cat > "$FAKE/skills/_attic/retired_skill/SKILL.md" <<'EOF'
---
name: WRONG_NAME_ON_PURPOSE
description: a skill that was retired to the attic
---
> **Deprecated.** Reference-only.
body
EOF
# A README in the attic (as the convention prescribes) must not be mistaken
# for a skill either.
echo "# attic" > "$FAKE/skills/_attic/README.md"

# =============================================================================
# (a) _attic excluded from active enumeration + not resolvable by activate
# =============================================================================

skills_out=$("$SCIAGENT" list skills)

# retired_skill must NOT appear as an active (bare or tagged) skill line. It may
# only appear under the Attic section, prefixed with `_attic/`.
active_lines=$(printf '%s\n' "$skills_out" | sed -n '/Attic (retired/q;p')
if printf '%s\n' "$active_lines" | grep -qw 'retired_skill'; then
    echo "FAIL [$_TEST_NAME] retired_skill leaked into the active skill list" >&2
    printf '%s\n' "$skills_out" >&2
    exit 1
fi

# It must NOT be resolvable by activate: build a role that references it and
# confirm activate fails to mount it (the resolver looks at skills/<name>/,
# never skills/_attic/<name>/).
cat > "$FAKE/roles/atticref.yaml" <<'EOF'
name: atticref
description: role that (wrongly) references an attic'd skill
skills:
  - retired_skill
EOF
mkdir -p "$TMPDIR_TEST/proj_attic"
# The resolver must HARD-FAIL on an unresolvable (attic'd) skill, not silently
# skip it — assert the non-zero exit (captured from the subshell), not just
# symlink-absence.
attic_rc=0
( cd "$TMPDIR_TEST/proj_attic" && "$SCIAGENT" activate atticref ) >/dev/null 2>&1 || attic_rc=$?
if [[ "$attic_rc" -eq 0 ]]; then
    echo "FAIL [$_TEST_NAME] activate succeeded on an unresolvable attic'd skill (expected hard-fail)" >&2
    exit 1
fi
if [[ -L "$TMPDIR_TEST/proj_attic/.claude/skills/retired_skill" ]]; then
    echo "FAIL [$_TEST_NAME] activate mounted an attic'd skill (symlink created)" >&2
    exit 1
fi

# A clean role (base) must still activate and mount its real skills, and must
# NOT create any retired_skill symlink.
mkdir -p "$TMPDIR_TEST/proj_base"
(
    cd "$TMPDIR_TEST/proj_base"
    "$SCIAGENT" activate base >/dev/null
)
assert_symlink "$TMPDIR_TEST/proj_base/.claude/skills/s_a" "base mounts s_a"
if [[ -e "$TMPDIR_TEST/proj_base/.claude/skills/retired_skill" ]]; then
    echo "FAIL [$_TEST_NAME] retired_skill symlink present after activating base" >&2
    exit 1
fi

# =============================================================================
# (b) list skills separates the Attic from the active catalog
# =============================================================================

# The Attic section must be present and list the retired skill.
if ! printf '%s\n' "$skills_out" | grep -q 'Attic (retired'; then
    echo "FAIL [$_TEST_NAME] 'list skills' missing the Attic section" >&2
    printf '%s\n' "$skills_out" >&2
    exit 1
fi
if ! printf '%s\n' "$skills_out" | grep -q '_attic/retired_skill'; then
    echo "FAIL [$_TEST_NAME] retired_skill not listed in the Attic section" >&2
    printf '%s\n' "$skills_out" >&2
    exit 1
fi

# =============================================================================
# (c) validate ignores _attic
# =============================================================================
# retired_skill's frontmatter `name:` deliberately disagrees with its directory,
# which hard-fails validate IF the attic is walked. It must not be.
if ! "$SCIAGENT" validate --quiet; then
    echo "FAIL [$_TEST_NAME] validate --quiet failed; _attic skill was likely walked" >&2
    "$SCIAGENT" validate >&2 || true
    exit 1
fi

pass
