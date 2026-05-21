#!/usr/bin/env bash
# tests/test_skill_frontmatter_valid.sh — every skill SKILL.md must have:
#   1. A YAML frontmatter block (lines 1..N delimited by ---)
#   2. A `name:` top-level field that matches the parent directory name
#   3. A `metadata:` top-level block
#
# The _TEMPLATE directory is exempt.

set -u
. "$(dirname "$0")/_lib.sh"

declare -i fail_count=0

for skill_dir in "$TOOLKIT_ROOT"/skills/*/; do
    name="$(basename "$skill_dir")"
    [[ "$name" == "_TEMPLATE" ]] && continue

    file="$skill_dir/SKILL.md"
    if [[ ! -f "$file" ]]; then
        echo "FAIL [$_TEST_NAME] $name: missing SKILL.md" >&2
        fail_count+=1
        continue
    fi

    # 1. Frontmatter delimited by --- ... ---
    first_line="$(head -n1 "$file")"
    if [[ "$first_line" != "---" ]]; then
        echo "FAIL [$_TEST_NAME] $name: no opening frontmatter '---' (got '$first_line')" >&2
        fail_count+=1
        continue
    fi

    # Find closing ---
    if ! awk 'NR==1 && /^---$/ {next} /^---$/ {found=1; exit} END{exit !found}' "$file"; then
        echo "FAIL [$_TEST_NAME] $name: no closing frontmatter '---'" >&2
        fail_count+=1
        continue
    fi

    # Extract frontmatter
    fm="$(awk 'BEGIN{c=0} /^---$/{c++; if(c==2) exit; next} c==1{print}' "$file")"

    # 2. name: matches dirname
    skill_name="$(printf '%s\n' "$fm" | awk '/^name:[ \t]/ {sub(/^name:[ \t]*/, ""); gsub(/^["'\'']|["'\'']$/, ""); print; exit}')"
    if [[ -z "$skill_name" ]]; then
        echo "FAIL [$_TEST_NAME] $name: missing top-level 'name:' in frontmatter" >&2
        fail_count+=1
        continue
    fi
    if [[ "$skill_name" != "$name" ]]; then
        echo "FAIL [$_TEST_NAME] $name: frontmatter name='$skill_name' does not match dirname" >&2
        fail_count+=1
        continue
    fi

    # 3. metadata: block present
    if ! printf '%s\n' "$fm" | grep -q '^metadata:[ \t]*$'; then
        echo "FAIL [$_TEST_NAME] $name: missing top-level 'metadata:' block" >&2
        fail_count+=1
        continue
    fi
done

if [[ "$fail_count" -gt 0 ]]; then
    echo "FAIL [$_TEST_NAME] $fail_count skill(s) with invalid frontmatter" >&2
    exit 1
fi

pass
