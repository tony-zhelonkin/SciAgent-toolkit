# lib/sciagent/inject.sh — sciagent inject <skill> | --tag <name>
# Adds skill(s) on top of the current stack (per architecture §5).
# - Stack [base]            → synthesize implicit overlay "_injected".
# - Stack [base, overlay]   → record skill as injected-into-overlay.
# Idempotent: re-injecting a tracked skill is a no-op.
#
# Forms:
#   sciagent inject <skill>        — inject one named skill (via: "")
#   sciagent inject --tag <name>   — inject all skills tagged <name>
#                                    (via: "tag:<name>" per entry)
#
# Block-body rendering is delegated to stack.sh:render_block_body.

# shellcheck shell=bash

cmd_inject() {
    if [[ $# -eq 0 ]]; then
        echo "usage: sciagent inject <skill> | --tag <name>" >&2
        return 1
    fi

    if ! manifest_exists; then
        echo "no role active; run \`sciagent activate <role>\` first" >&2
        return 1
    fi

    # Dispatch on argument form.
    case "${1:-}" in
        --tag)
            local tag_name="${2:-}"
            if [[ -z "$tag_name" ]]; then
                echo "usage: sciagent inject --tag <name>" >&2
                return 1
            fi
            _inject_by_tag "$tag_name"
            return $?
            ;;
        --*)
            echo "sciagent inject: unknown option '$1'" >&2
            echo "usage: sciagent inject <skill> | --tag <name>" >&2
            return 1
            ;;
        *)
            if [[ $# -ne 1 ]]; then
                echo "usage: sciagent inject <skill> | --tag <name>" >&2
                return 1
            fi
            _inject_named_skill "$1" ""
            return $?
            ;;
    esac
}

# _inject_named_skill <skill> <via>
# Injects a single skill; records the given via value in the manifest.
# Returns 0 on success or idempotent no-op.
_inject_named_skill() {
    local skill="$1" via="$2"

    local skill_src
    skill_src=$(resolve_canonical skills "$skill") || return 1

    local stack base overlay
    stack=$(manifest_stack)
    base=$(printf '%s\n' "$stack" | awk '{print $1}')
    overlay=$(printf '%s\n' "$stack" | awk '{print $2}')

    # Determine target overlay.
    local target_overlay
    if [[ -z "$overlay" ]]; then
        target_overlay="_injected"
    else
        target_overlay="$overlay"
    fi

    # Idempotency: already injected?
    local inj_ov inj_sk _inj_via
    while read -r inj_ov inj_sk _inj_via; do
        if [[ "$inj_ov" == "$target_overlay" && "$inj_sk" == "$skill" ]]; then
            echo "already injected: $skill (into $target_overlay)"
            return 0
        fi
    done < <(manifest_injected)

    # Create symlinks.
    mkdir -p .claude/skills .agents/skills
    local claude_path=".claude/skills/$skill"
    local agents_path=".agents/skills/$skill"
    ln -sfn "$skill_src" "$claude_path"
    ln -sfn "$skill_src" "$agents_path"

    # Persist new symlinks and injected entry (with via).
    manifest_append_inject "$claude_path" "$agents_path" "$target_overlay" "$skill" "$via"

    # Update the stack line if synthesizing _injected.
    if [[ -z "$overlay" ]]; then
        manifest_update_stack "$base _injected"
    fi

    # Re-render the AGENTS.md managed block to reflect injection.
    _inject_rewrite_block

    echo "injected: $skill (into $target_overlay)"
    return 0
}

# _inject_by_tag <tag-name>
# Reads tags.yaml; asserts the tag exists in vocabulary; collects every skill
# whose metadata.tags: contains <tag-name>; injects each with via:"tag:<name>".
# Empty match (valid tag, no skills carry it): warn + exit 0.
_inject_by_tag() {
    local tag_name="$1"

    # Locate toolkit root.
    local tk_root="${SCIAGENT_TOOLKIT:-}"
    if [[ -z "$tk_root" ]]; then
        echo "sciagent inject --tag: SCIAGENT_TOOLKIT not set" >&2
        return 1
    fi

    local tags_file="$tk_root/tags.yaml"
    if [[ ! -f "$tags_file" ]]; then
        echo "sciagent inject --tag: tags.yaml not found at $tags_file" >&2
        return 1
    fi

    # Assert the tag exists in vocabulary.
    local known_tags
    known_tags="$(awk '
        /^tags:/ { intags=1; next }
        intags && /^  - name:/ {
            sub(/^  - name:[ \t]*/, "")
            sub(/[ \t]*#.*$/, "")
            sub(/[ \t]+$/, "")
            gsub(/^["'"'"']|["'"'"']$/, "")
            print
        }
        intags && /^[^ ]/ { intags=0 }
    ' "$tags_file")"

    if ! printf '%s\n' "$known_tags" | grep -qxF "$tag_name"; then
        echo "sciagent inject --tag: '$tag_name' is not in tags.yaml vocabulary" >&2
        echo "  Known tags:" >&2
        printf '%s\n' "$known_tags" | sed 's/^/    /' >&2
        return 1
    fi

    # Collect skills whose metadata.tags: contains tag_name.
    local skills_dir="$tk_root/skills"
    local -a matching_skills=()
    local skill_dir skill_name skill_file

    for skill_dir in "$skills_dir"/*/; do
        skill_name="$(basename "$skill_dir")"
        [[ "$skill_name" == "_TEMPLATE" ]] && continue
        skill_file="$skill_dir/SKILL.md"
        [[ -f "$skill_file" ]] || continue

        # Extract tags from frontmatter metadata.tags: block list.
        local skill_tags
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
        ' "$skill_file")"

        if printf '%s\n' "$skill_tags" | grep -qxF "$tag_name"; then
            matching_skills+=("$skill_name")
        fi
    done

    if [[ "${#matching_skills[@]}" -eq 0 ]]; then
        echo "sciagent inject --tag: no skills currently carry the '$tag_name' tag — nothing to inject" >&2
        return 0
    fi

    local via="tag:$tag_name"
    local injected_count=0
    local sk
    for sk in "${matching_skills[@]}"; do
        _inject_named_skill "$sk" "$via"
        local rc=$?
        if [[ "$rc" -eq 0 ]]; then
            (( injected_count++ )) || true
        fi
    done

    return 0
}

# _inject_rewrite_block — rebuild AGENTS.md block from current manifest + role YAMLs.
_inject_rewrite_block() {
    local stack base overlay
    stack=$(manifest_stack)
    base=$(printf '%s\n' "$stack" | awk '{print $1}')
    overlay=$(printf '%s\n' "$stack" | awk '{print $2}')

    # Collect injected skill names from manifest.
    local -a INJECTED_NAMES=()
    local ov nm _via
    while read -r ov nm _via; do
        [[ -n "$ov" ]] && INJECTED_NAMES+=("$nm")
    done < <(manifest_injected)

    local body
    body=$(render_block_body "$base" "$overlay" "${INJECTED_NAMES[@]+"${INJECTED_NAMES[@]}"}")
    block_write AGENTS.md "$body"

    # Refresh manifest BLOCK_HASH.
    manifest_update_block_hash "$(block_stored_hash AGENTS.md)"
}
