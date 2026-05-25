# lib/sciagent/eject.sh — sciagent eject <skill> | --tag <name>
# Symmetric counterpart to `inject` (per kickoff.md §9, 2026-05-24
# ADR-001 inject↔eject verb pair).
#
# Removes skills previously placed via `sciagent inject`. Skills mounted
# as part of a role-stack layer are NOT eligible — use `sciagent deactivate`
# for those. The two verbs cover orthogonal concepts (stack composition vs.
# ad-hoc injection).
#
# Forms:
#   sciagent eject <skill>        — remove one previously-injected skill
#   sciagent eject --tag <name>   — remove all entries with via=="tag:<name>"

# shellcheck shell=bash

cmd_eject() {
    if ! manifest_exists; then
        echo "no role active; nothing to eject" >&2
        return 1
    fi

    if [[ $# -eq 0 ]]; then
        echo "usage: sciagent eject <skill> | --tag <name>" >&2
        return 1
    fi

    # Parse: either `eject <skill>` or `eject --tag <name>`.
    local mode="" target=""
    case "${1:-}" in
        --tag)
            mode="tag"
            target="${2:-}"
            if [[ -z "$target" ]]; then
                echo "usage: sciagent eject --tag <name>" >&2
                return 1
            fi
            ;;
        --*)
            echo "sciagent eject: unknown option '$1'" >&2
            echo "usage: sciagent eject <skill> | --tag <name>" >&2
            return 1
            ;;
        *)
            mode="skill"
            target="$1"
            ;;
    esac

    case "$mode" in
        skill) _eject_named_skill "$target" ;;
        tag)   _eject_by_tag "$target" ;;
    esac
}

# _eject_named_skill <skill>
# Removes one previously-injected skill by name.
# Refuses to remove stack-mounted skills (pointer to deactivate).
# Idempotent: ejecting a skill not in the injected list is a no-op (exit 0).
_eject_named_skill() {
    local skill="$1"

    local stack base overlay
    stack=$(manifest_stack)
    base=$(printf '%s\n' "$stack" | awk '{print $1}')
    overlay=$(printf '%s\n' "$stack" | awk '{print $2}')

    # Guard: refuse if the skill is stack-mounted (not injected).
    # Stack-mounted means it appears in the role YAML, not in injected[].
    # Check: is the skill in injected[]?
    local found_injected=0
    local found_overlay="" found_via=""
    local inj_ov inj_sk inj_via
    while read -r inj_ov inj_sk inj_via; do
        if [[ "$inj_sk" == "$skill" ]]; then
            found_injected=1
            found_overlay="$inj_ov"
            found_via="$inj_via"
            break
        fi
    done < <(manifest_injected)

    if [[ "$found_injected" -eq 0 ]]; then
        # Not in injected[]. Check if it is stack-mounted.
        if _skill_is_stack_mounted "$skill" "$base" "$overlay"; then
            echo "sciagent eject: skill '$skill' is part of stack-mounted role '$base'${overlay:+/$overlay}; use 'sciagent deactivate' instead" >&2
            return 1
        fi
        # Not in injected and not stack-mounted: simply not present.
        echo "not injected: $skill — nothing to eject"
        return 0
    fi

    # Remove the dual symlinks for this skill.
    _eject_remove_symlinks "$skill"

    # Rebuild the manifest without this entry.
    _eject_drop_injected_entry "$skill" ""

    # Collapse _injected overlay if it is now empty and was synthetic.
    _eject_maybe_collapse_injected_overlay

    # Re-render AGENTS.md block.
    _eject_rewrite_block

    echo "ejected: $skill"
    return 0
}

# _eject_by_tag <tag-name>
# Removes all injected entries whose via value equals "tag:<tag-name>".
# Idempotent: if no entries match, prints a notice and exits 0.
_eject_by_tag() {
    local tag_name="$1"
    local via_key="tag:$tag_name"

    # Collect all injected entries matching this tag.
    local -a to_eject_skills=()
    local inj_ov inj_sk inj_via
    while read -r inj_ov inj_sk inj_via; do
        if [[ "$inj_via" == "$via_key" ]]; then
            to_eject_skills+=("$inj_sk")
        fi
    done < <(manifest_injected)

    if [[ "${#to_eject_skills[@]}" -eq 0 ]]; then
        echo "not injected via tag '$tag_name' — nothing to eject"
        return 0
    fi

    local sk
    for sk in "${to_eject_skills[@]}"; do
        # Remove symlinks.
        _eject_remove_symlinks "$sk"
        # Drop manifest entry (match skill + via).
        _eject_drop_injected_entry "$sk" "$via_key"
        echo "ejected: $sk (tag:$tag_name)"
    done

    # Collapse _injected overlay if now empty and synthetic.
    _eject_maybe_collapse_injected_overlay

    # Re-render AGENTS.md block.
    _eject_rewrite_block

    return 0
}

# _skill_is_stack_mounted <skill> <base> <overlay>
# Returns 0 if the skill appears in the roles for base or overlay.
_skill_is_stack_mounted() {
    local skill="$1" base="$2" overlay="$3"
    local tk_root="${SCIAGENT_TOOLKIT:-}"

    _role_contains_skill() {
        local role="$1" sk="$2"
        local role_file="$tk_root/roles/$role.yaml"
        [[ -f "$role_file" ]] || return 1
        grep -q "^  - $sk$" "$role_file"
    }

    if _role_contains_skill "$base" "$skill"; then
        return 0
    fi
    if [[ -n "$overlay" ]] && _role_contains_skill "$overlay" "$skill"; then
        return 0
    fi
    return 1
}

# _eject_remove_symlinks <skill>
# Removes the dual .claude/skills/<skill> and .agents/skills/<skill> symlinks.
_eject_remove_symlinks() {
    local skill="$1"
    local claude_path=".claude/skills/$skill"
    local agents_path=".agents/skills/$skill"

    if [[ -L "$claude_path" ]]; then
        rm "$claude_path"
    fi
    if [[ -L "$agents_path" ]]; then
        rm "$agents_path"
    fi
}

# _eject_drop_injected_entry <skill> <via>
# Rewrites the manifest, omitting the injected entry for <skill>.
# When <via> is non-empty, also matches via (used by --tag eject).
# When <via> is empty, matches by skill name alone.
_eject_drop_injected_entry() {
    local skill="$1" via_filter="$2"

    local old_stack old_hash
    old_stack=$(manifest_stack)
    old_hash=$(manifest_block_hash)

    _manifest_stack_val="$old_stack"
    _manifest_syms=()
    _manifest_injected_overlays=()
    _manifest_injected_skills=()
    _manifest_injected_vias=()

    # Re-read symlinks, dropping the two entries for this skill.
    local p
    while IFS= read -r p; do
        [[ -z "$p" ]] && continue
        # Drop the two symlink entries for this skill (both tracks).
        if [[ "$p" == ".claude/skills/$skill" || "$p" == ".agents/skills/$skill" ]]; then
            continue
        fi
        _manifest_syms+=("$p")
    done < <(manifest_symlinks)

    # Re-read injected entries, dropping the matched one.
    local inj_ov inj_sk inj_via
    local dropped=0
    while read -r inj_ov inj_sk inj_via; do
        [[ -z "$inj_ov" ]] && continue
        # Match: skill name matches AND (via_filter empty OR via matches).
        if [[ "$inj_sk" == "$skill" && "$dropped" -eq 0 ]]; then
            if [[ -z "$via_filter" || "${inj_via:-}" == "$via_filter" ]]; then
                dropped=1
                continue  # Skip this entry.
            fi
        fi
        _manifest_injected_overlays+=("$inj_ov")
        _manifest_injected_skills+=("$inj_sk")
        _manifest_injected_vias+=("${inj_via:-}")
    done < <(manifest_injected)

    _manifest_staging=$(mktemp)
    _manifest_write_json "$old_hash"
    mv "$_manifest_staging" "$_MANIFEST_PATH"
    _manifest_staging=""
}

# _eject_maybe_collapse_injected_overlay
# If the _injected synthetic overlay exists in the stack AND injected[] now
# contains no entries for _injected, collapse the stack back to [base].
_eject_maybe_collapse_injected_overlay() {
    local stack base overlay
    stack=$(manifest_stack)
    base=$(printf '%s\n' "$stack" | awk '{print $1}')
    overlay=$(printf '%s\n' "$stack" | awk '{print $2}')

    # Only act when the synthetic overlay is active.
    [[ "$overlay" == "_injected" ]] || return 0

    # Check if any injected entries still reference _injected.
    local inj_ov inj_sk _inj_via
    while read -r inj_ov inj_sk _inj_via; do
        if [[ "$inj_ov" == "_injected" ]]; then
            return 0  # Still has entries; don't collapse.
        fi
    done < <(manifest_injected)

    # No entries remain for _injected: collapse stack to [base].
    manifest_update_stack "$base"
}

# _eject_rewrite_block — rebuild AGENTS.md block from current manifest + role YAMLs.
_eject_rewrite_block() {
    local stack base overlay
    stack=$(manifest_stack)
    base=$(printf '%s\n' "$stack" | awk '{print $1}')
    overlay=$(printf '%s\n' "$stack" | awk '{print $2}')

    # Collect remaining injected skill names from manifest.
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
