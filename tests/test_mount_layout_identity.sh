#!/usr/bin/env bash
# The flat catalog and both mounted trees expose the same category name sets.

set -u
. "$(dirname "$0")/_lib.sh"

setup_tmpdir

expected_skills="$TMPDIR_TEST/expected-skills"
expected_agents="$TMPDIR_TEST/expected-agents"
expected_commands="$TMPDIR_TEST/expected-commands"

for d in "$TOOLKIT_ROOT"/skills/*/; do
    [[ -f "$d/SKILL.md" ]] || {
        echo "FAIL [$_TEST_NAME] non-skill directory in skills/: $d" >&2
        exit 1
    }
    basename "$d"
done | LC_ALL=C sort > "$expected_skills"

for f in "$TOOLKIT_ROOT"/agents/*; do
    [[ -f "$f" && "$f" == *.md ]] || {
        echo "FAIL [$_TEST_NAME] non-agent entry in agents/: $f" >&2
        exit 1
    }
    basename "$f"
done | LC_ALL=C sort > "$expected_agents"

for f in "$TOOLKIT_ROOT"/commands/*; do
    [[ -f "$f" && "$f" == *.md ]] || {
        echo "FAIL [$_TEST_NAME] non-command entry in commands/: $f" >&2
        exit 1
    }
    basename "$f"
done | LC_ALL=C sort > "$expected_commands"

mkdir project
cd project
"$TOOLKIT_ROOT/bin/scio" link >/dev/null

for harness in .claude .agents; do
    assert_symlink "$harness/skills"
    assert_symlink "$harness/agents"
    assert_symlink "$harness/commands"

    find -H "$harness/skills" -mindepth 1 -maxdepth 1 -printf '%f\n' \
        | LC_ALL=C sort > "$TMPDIR_TEST/actual-skills"
    find -H "$harness/agents" -mindepth 1 -maxdepth 1 -printf '%f\n' \
        | LC_ALL=C sort > "$TMPDIR_TEST/actual-agents"
    find -H "$harness/commands" -mindepth 1 -maxdepth 1 -printf '%f\n' \
        | LC_ALL=C sort > "$TMPDIR_TEST/actual-commands"

    assert_file_eq "$TMPDIR_TEST/actual-skills" "$expected_skills" \
        "$harness skill names match skills/"
    assert_file_eq "$TMPDIR_TEST/actual-agents" "$expected_agents" \
        "$harness agent names match agents/"
    assert_file_eq "$TMPDIR_TEST/actual-commands" "$expected_commands" \
        "$harness command names match commands/"
done

grep -Fxq 'bio-interpreter.md' "$expected_agents" || {
    echo "FAIL [$_TEST_NAME] bio-interpreter.md lost its mounted name" >&2
    exit 1
}
grep -Fxq 'add-figure-variant.md' "$expected_commands" || {
    echo "FAIL [$_TEST_NAME] add-figure-variant.md lost its mounted name" >&2
    exit 1
}

pass
