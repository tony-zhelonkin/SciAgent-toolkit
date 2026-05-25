#!/usr/bin/env bash
# tests/test_tags_vocabulary.sh — every metadata.tags: value across all skills
# must exist as a `name:` entry in tags.yaml at the toolkit root.
#
# Per kickoff.md §9 (2026-05-24) ADR-001 tag vocabulary: tags.yaml is the
# single source of truth; unknown tags in skill frontmatter are a FAIL.
#
# The _TEMPLATE skill is exempt (it carries placeholder values that are not
# real tags).

set -u
. "$(dirname "$0")/_lib.sh"

TAGS_FILE="$TOOLKIT_ROOT/tags.yaml"

if [[ ! -f "$TAGS_FILE" ]]; then
    echo "FAIL [$_TEST_NAME] tags.yaml missing at toolkit root ($TAGS_FILE)" >&2
    exit 1
fi

# Build the set of known tag names from tags.yaml.
# Parses lines of the form: `  - name: <value>` under the `tags:` block.
known_tags="$(awk '
    /^tags:/ { intags=1; next }
    intags && /^  - name:/ {
        sub(/^  - name:[ \t]*/, "")
        sub(/[ \t]*#.*$/, "")
        sub(/[ \t]+$/, "")
        gsub(/^["'"'"']|["'"'"']$/, "")
        print
        next
    }
    intags && /^[^ ]/ { intags=0 }
' "$TAGS_FILE")"

if [[ -z "$known_tags" ]]; then
    echo "FAIL [$_TEST_NAME] tags.yaml contains no parseable tag names" >&2
    exit 1
fi

declare -i fail_count=0
declare -i pass_count=0

for skill_dir in "$TOOLKIT_ROOT"/skills/*/; do
    name="$(basename "$skill_dir")"
    [[ "$name" == "_TEMPLATE" ]] && continue

    file="$skill_dir/SKILL.md"
    [[ -f "$file" ]] || continue

    # Extract the tags list items from the metadata block.
    # Matches lines like `  - some-tag` that appear after `  tags:` inside
    # the metadata block and before any non-list continuation.
    skill_tags="$(awk '
        BEGIN { infm=0; closed=0; inmeta=0; intags=0 }
        NR==1 && /^---[ \t]*$/ { infm=1; next }
        infm && !closed && /^---[ \t]*$/ { closed=1; exit }
        !infm { next }
        /^metadata:[ \t]*$/ { inmeta=1; next }
        inmeta && /^  tags:/ { intags=1; next }
        intags && /^  - / {
            val=$0
            sub(/^  - /, "", val)
            sub(/[ \t]*#.*$/, "", val)
            sub(/[ \t]+$/, "", val)
            gsub(/^["'"'"']|["'"'"']$/, "", val)
            print val
            next
        }
        intags && /^  [^ ]/ { intags=0 }
        intags && /^[^ ]/ { intags=0 }
    ' "$file")"

    skill_ok=1
    while IFS= read -r tag; do
        [[ -z "$tag" ]] && continue
        if ! printf '%s\n' "$known_tags" | grep -qxF "$tag"; then
            echo "FAIL [$_TEST_NAME] $name: unknown tag '$tag' (not in tags.yaml)" >&2
            fail_count+=1
            skill_ok=0
        fi
    done <<< "$skill_tags"

    [[ "$skill_ok" -eq 1 ]] && pass_count+=1
done

if [[ "$fail_count" -gt 0 ]]; then
    echo "FAIL [$_TEST_NAME] $fail_count unknown tag reference(s) across skills; $pass_count skills ok" >&2
    exit 1
fi

pass
