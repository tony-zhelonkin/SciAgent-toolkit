#!/usr/bin/env bash
# A project keeping its own skills gets one link per catalog entry beside them,
# instead of a single category symlink that would swallow the directory.
#
#   1. A foreign entry turns the category into a fanout: local skill survives,
#      every catalog entry is reachable.
#   2. The fanout is idempotent — a second run rewrites nothing.
#   3. A local entry named like a catalog entry keeps the local one and says so.
#   4. Removing the local entry converges back to the single category symlink.
#   5. The other categories are untouched and stay single symlinks.

set -u
. "$(dirname "$0")/_lib.sh"

setup_tmpdir
SCIO="$TOOLKIT_ROOT/bin/scio"
PROJ="$TMPDIR_TEST/project"

mkdir -p "$PROJ/.claude/skills/my-local-skill"
cat >"$PROJ/.claude/skills/my-local-skill/SKILL.md" <<'EOF'
---
name: my-local-skill
description: A skill this project installed for itself.
---
EOF

# --- 1. fanout keeps the local skill and binds the catalog ------------------
out=$("$SCIO" link --project-dir "$PROJ" 2>&1) || {
    echo "FAIL [$_TEST_NAME] link failed with a local skill present" >&2
    printf '%s\n' "$out" >&2
    exit 1
}

assert_file_exists "$PROJ/.claude/skills/my-local-skill/SKILL.md"
if [[ -L "$PROJ/.claude/skills" ]]; then
    echo "FAIL [$_TEST_NAME] .claude/skills became a symlink, swallowing the local skill" >&2
    exit 1
fi

catalog_count=$(find "$TOOLKIT_ROOT/skills" -mindepth 1 -maxdepth 1 -type d | wc -l)
linked_count=$(find "$PROJ/.claude/skills" -mindepth 1 -maxdepth 1 -type l | wc -l)
assert_eq "$linked_count" "$catalog_count" "one link per catalog skill"

# Pick a catalog entry by inspection rather than by name, so a rename in the
# catalog cannot silently turn this into a test of nothing.
CATALOG_SKILL=$(find "$TOOLKIT_ROOT/skills" -mindepth 1 -maxdepth 1 -type d -printf '%f\n' \
    | LC_ALL=C sort | head -1)
assert_symlink "$PROJ/.claude/skills/$CATALOG_SKILL"
assert_file_exists "$PROJ/.claude/skills/$CATALOG_SKILL/SKILL.md"

# --- 2. the fanout converges quietly ----------------------------------------
second=$("$SCIO" link --project-dir "$PROJ" 2>&1) || {
    echo "FAIL [$_TEST_NAME] second link failed" >&2
    printf '%s\n' "$second" >&2
    exit 1
}
if printf '%s\n' "$second" | grep -qE 'replaced link|^linked:|removed stale'; then
    echo "FAIL [$_TEST_NAME] the fanout is not idempotent" >&2
    printf '%s\n' "$second" >&2
    exit 1
fi

# --- 3. a shadowing local entry is kept and reported ------------------------
rm "$PROJ/.claude/skills/$CATALOG_SKILL"
mkdir -p "$PROJ/.claude/skills/$CATALOG_SKILL"
echo "local override" >"$PROJ/.claude/skills/$CATALOG_SKILL/SKILL.md"

shadow=$("$SCIO" link --project-dir "$PROJ" 2>&1) || {
    echo "FAIL [$_TEST_NAME] link failed on a shadowing entry" >&2
    printf '%s\n' "$shadow" >&2
    exit 1
}
assert_eq "$(cat "$PROJ/.claude/skills/$CATALOG_SKILL/SKILL.md")" "local override" \
    "the project copy survives"
if ! printf '%s\n' "$shadow" | grep -q 'shadows the catalog entry'; then
    echo "FAIL [$_TEST_NAME] a shadowing entry was not reported" >&2
    printf '%s\n' "$shadow" >&2
    exit 1
fi

# --- 4. removing every local entry converges to the single symlink ----------
rm -rf "$PROJ/.claude/skills/my-local-skill" "$PROJ/.claude/skills/$CATALOG_SKILL"
converge=$("$SCIO" link --project-dir "$PROJ" 2>&1) || {
    echo "FAIL [$_TEST_NAME] link failed while converging" >&2
    printf '%s\n' "$converge" >&2
    exit 1
}
assert_symlink "$PROJ/.claude/skills"
assert_eq "$(readlink "$PROJ/.claude/skills")" "$TOOLKIT_ROOT/skills" \
    ".claude/skills target after converging"

# --- 5. the untouched categories never left single-symlink form -------------
for path in "$PROJ/.claude/agents" "$PROJ/.claude/commands" \
            "$PROJ/.agents/skills" "$PROJ/.agents/agents" "$PROJ/.agents/commands"; do
    assert_symlink "$path"
done

pass
