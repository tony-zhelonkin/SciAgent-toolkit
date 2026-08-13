#!/usr/bin/env bash
# tests/test_skill_frontmatter_valid.sh — every skill SKILL.md must have:
#   1. A YAML frontmatter block (lines 1..N delimited by ---)
#   2. A `name:` top-level field that matches the parent directory name
#   3. NO `metadata:` block — taxonomy lives in the body, not the frontmatter
#   4. A `description:` within the toolkit catalog length cap
#
set -u
. "$(dirname "$0")/_lib.sh"

declare -i fail_count=0
declare desc=""

for skill_dir in "$TOOLKIT_ROOT"/skills/*/; do
    name="$(basename "$skill_dir")"
    [[ "$name" == _* ]] && continue
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

    # 3. metadata: block ABSENT. Cross-references and caveats belong in the
    # body (## Prerequisites / ## When not to use / ## See also); a taxonomy
    # block is dead weight preloaded for every skill in the catalog.
    if printf '%s\n' "$fm" | grep -q '^metadata:[ \t]*$'; then
        echo "FAIL [$_TEST_NAME] $name: has a 'metadata:' block; move it into the body" >&2
        fail_count+=1
        continue
    fi

    # 4. description: present and within the cap.
    desc="$(printf '%s\n' "$fm" | awk '/^description:[ \t]/ {sub(/^description:[ \t]*/, ""); sub(/[ \t]+$/, ""); gsub(/^["'\'']|["'\'']$/, ""); print; exit}')"
    if [[ -z "$desc" ]]; then
        echo "FAIL [$_TEST_NAME] $name: missing top-level 'description:' in frontmatter" >&2
        fail_count+=1
        continue
    fi
    if (( ${#desc} > 350 )); then
        echo "FAIL [$_TEST_NAME] $name: description is ${#desc} chars (max 350)" >&2
        fail_count+=1
        continue
    fi
done

if [[ "$fail_count" -gt 0 ]]; then
    echo "FAIL [$_TEST_NAME] $fail_count skill(s) with invalid frontmatter" >&2
    exit 1
fi

pass
