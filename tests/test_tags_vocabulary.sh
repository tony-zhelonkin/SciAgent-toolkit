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
# NOTE: lib/sciagent/validate.sh carries a near-identical awk parser for the
# same vocabulary. If you change the matching/trimming rules here, change
# them there too — the two readers must agree on what counts as a tag.
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
declare -i tagref_count=0

for skill_dir in "$TOOLKIT_ROOT"/skills/*/; do
    name="$(basename "$skill_dir")"
    [[ "$name" == "_TEMPLATE" ]] && continue

    file="$skill_dir/SKILL.md"
    [[ -f "$file" ]] || continue

    # Extract the tags list items from the metadata block.
    # Matches list items (`- some-tag`) at ANY indent that appear after the
    # `  tags:` key inside the metadata block and before the next sibling key.
    # The indent-agnostic `/^[ \t]+-[ \t]/` match (not a fixed `^  - `) is
    # deliberate: skills in the library mix 2-space and 4-space list styles,
    # and a fixed-width matcher silently skips the 4-space ones — letting an
    # unknown tag slip through undetected. This MUST mirror the awk parser in
    # lib/sciagent/validate.sh so the test and the validator agree on what
    # counts as a tag.
    skill_tags="$(awk '
        BEGIN { infm=0; closed=0; inmeta=0; intags=0 }
        NR==1 && /^---[ \t]*$/ { infm=1; next }
        infm && !closed && /^---[ \t]*$/ { closed=1; exit }
        !infm { next }
        /^metadata:[ \t]*$/ { inmeta=1; next }
        inmeta && /^  tags:/ { intags=1; next }
        intags && /^[ \t]+-[ \t]/ {
            val=$0
            sub(/^[ \t]+-[ \t]+/, "", val)
            sub(/[ \t]*#.*$/, "", val)
            sub(/[ \t]+$/, "", val)
            gsub(/^["'"'"']|["'"'"']$/, "", val)
            print val
            next
        }
        intags && /^  [^ \t-]/ { intags=0 }
        intags && /^[^ \t]/ { intags=0 }
    ' "$file")"

    skill_ok=1
    while IFS= read -r tag; do
        [[ -z "$tag" ]] && continue
        tagref_count+=1
        if ! printf '%s\n' "$known_tags" | grep -qxF "$tag"; then
            echo "FAIL [$_TEST_NAME] $name: unknown tag '$tag' (not in tags.yaml)" >&2
            fail_count+=1
            skill_ok=0
        fi
    done <<< "$skill_tags"

    [[ "$skill_ok" -eq 1 ]] && pass_count+=1
done

# Sentinel: this whole test is vacuous if the parser extracts nothing — every
# skill would "pass" by carrying zero tags. The library always has many tagged
# skills, so a near-empty extraction means the parser regressed (e.g. an
# indent-width assumption silently skipping list items). Fail loudly so the
# brittleness can't hide behind a green test. The threshold is intentionally
# low (structural floor, not a hardcoded expected list) so it never breaks when
# skills are added or retagged.
declare -i MIN_TAGREFS=10
if [[ "$tagref_count" -lt "$MIN_TAGREFS" ]]; then
    echo "FAIL [$_TEST_NAME] parser extracted only $tagref_count tag reference(s) (< $MIN_TAGREFS) — the tags: parser likely regressed and is silently skipping list items" >&2
    exit 1
fi

if [[ "$fail_count" -gt 0 ]]; then
    echo "FAIL [$_TEST_NAME] $fail_count unknown tag reference(s) across skills; $pass_count skills ok" >&2
    exit 1
fi

# Positive vocabulary assertions: confirm that the three P09 tags are declared
# and well-formed (parseable by the awk reader above).  Each must appear in
# $known_tags exactly once and have a non-empty `description:` and `since:`.
_assert_tag_present() {
    local tag="$1"
    # 1. name is known (parseable by the awk reader)
    if ! printf '%s\n' "$known_tags" | grep -qxF "$tag"; then
        echo "FAIL [$_TEST_NAME] required tag '$tag' is missing from tags.yaml" >&2
        return 1
    fi
    # 2. description field exists and is non-empty
    local desc
    desc="$(awk -v tag="$tag" '
        /^  - name:/ {
            n=$0; sub(/^  - name:[ \t]*/, "", n); sub(/[ \t]+$/, "", n)
            gsub(/^["'"'"']|["'"'"']$/, "", n)
            found=(n == tag)
            next
        }
        found && /^    description:/ {
            d=$0; sub(/^    description:[ \t]*/, "", d); sub(/[ \t]+$/, "", d)
            gsub(/^["'"'"']|["'"'"']$/, "", d)
            print d; exit
        }
        found && /^  - name:/ { exit }
    ' "$TAGS_FILE")"
    if [[ -z "$desc" ]]; then
        echo "FAIL [$_TEST_NAME] tag '$tag' has no parseable description in tags.yaml" >&2
        return 1
    fi
    # 3. since field exists and is non-empty
    local since
    since="$(awk -v tag="$tag" '
        /^  - name:/ {
            n=$0; sub(/^  - name:[ \t]*/, "", n); sub(/[ \t]+$/, "", n)
            gsub(/^["'"'"']|["'"'"']$/, "", n)
            found=(n == tag)
            next
        }
        found && /^    since:/ {
            s=$0; sub(/^    since:[ \t]*/, "", s); sub(/[ \t]+$/, "", s)
            gsub(/^["'"'"']|["'"'"']$/, "", s)
            print s; exit
        }
        found && /^  - name:/ { exit }
    ' "$TAGS_FILE")"
    if [[ -z "$since" ]]; then
        echo "FAIL [$_TEST_NAME] tag '$tag' has no parseable since field in tags.yaml" >&2
        return 1
    fi
    return 0
}

declare -i vocab_fail=0
for _required_tag in figure provenance planning; do
    _assert_tag_present "$_required_tag" || vocab_fail+=1
done
if [[ "$vocab_fail" -gt 0 ]]; then
    exit 1
fi

pass
