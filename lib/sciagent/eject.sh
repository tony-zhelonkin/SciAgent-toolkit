# lib/sciagent/eject.sh — sciagent eject <name> [--skill|--agent|--command]
#                                          | --tag <name>
# Symmetric counterpart to `inject` (per kickoff.md §9, 2026-05-24 ADR-001
# inject↔eject verb pair; PR-A cross-kind extension).
#
# Removes entries previously placed via `sciagent inject`. Entries mounted as
# part of a role-stack layer are NOT eligible — use `sciagent deactivate` for
# those. The two verbs cover orthogonal concepts (stack composition vs.
# ad-hoc injection).
#
# Forms:
#   sciagent eject <name>               — find the matching injected row by
#                                          name, use its `kind` field to pick
#                                          the right symlink tree. Hard-fails
#                                          when more than one row matches by
#                                          name across kinds.
#   sciagent eject --skill <name>       — restrict to kind=skill rows
#   sciagent eject --agent <name>       — restrict to kind=agent rows
#   sciagent eject --command <name>     — restrict to kind=command rows
#   sciagent eject --tag <name>         — remove all rows with via=="tag:<name>"

# shellcheck shell=bash

cmd_eject() {
    if ! manifest_exists; then
        echo "no role active; nothing to eject" >&2
        return 1
    fi

    if [[ $# -eq 0 ]]; then
        echo "usage: sciagent eject <name> | --skill <name> | --agent <name> | --command <name> | --tag <name>" >&2
        return 1
    fi

    # Parse: `eject <name>`, `eject --<kind> <name>`, or `eject --tag <name>`.
    local mode="" target="" kind_filter=""
    case "${1:-}" in
        --tag)
            mode="tag"
            target="${2:-}"
            if [[ -z "$target" ]]; then
                echo "usage: sciagent eject --tag <name>" >&2
                return 1
            fi
            ;;
        --skill|--agent|--command)
            mode="kind"
            kind_filter="${1#--}"
            target="${2:-}"
            if [[ -z "$target" ]]; then
                echo "usage: sciagent eject --${kind_filter} <name>" >&2
                return 1
            fi
            ;;
        --*)
            echo "sciagent eject: unknown option '$1'" >&2
            echo "usage: sciagent eject <name> | --skill <name> | --agent <name> | --command <name> | --tag <name>" >&2
            return 1
            ;;
        '')
            echo "sciagent eject: empty name" >&2
            echo "usage: sciagent eject <name> | --skill <name> | --agent <name> | --command <name> | --tag <name>" >&2
            return 1
            ;;
        *)
            mode="name"
            target="$1"
            ;;
    esac

    case "$mode" in
        name) _eject_named_entry "$target" "" ;;
        kind) _eject_named_entry "$target" "$kind_filter" ;;
        tag)  _eject_by_tag "$target" ;;
    esac
}

# _eject_named_entry <name> <kind_filter>
# Removes one previously-injected entry by name. When kind_filter is non-empty
# (set by --skill/--agent/--command), only matches rows whose kind equals
# kind_filter — disambiguates the cross-kind case. When kind_filter is empty
# (bare-name eject), matches by name alone and hard-fails on ambiguity (two
# or more rows share the name across kinds).
#
# Refuses to remove stack-mounted entries (pointer to deactivate).
# Idempotent: ejecting a missing target is a no-op (exit 0).
_eject_named_entry() {
    local name="$1" kind_filter="$2"

    local stack base overlay
    stack=$(manifest_stack)
    base=$(printf '%s\n' "$stack" | awk '{print $1}')
    overlay=$(printf '%s\n' "$stack" | awk '{print $2}')

    # Collect every manifest row matching <name> (and kind_filter when set).
    local -a matched_overlays=() matched_vias=() matched_kinds=()
    local inj_ov inj_sk inj_via inj_kind
    while IFS='|' read -r inj_ov inj_sk inj_via inj_kind; do
        [[ -z "$inj_ov" ]] && continue
        # Forward-compat: pre-PR-A rows lack kind; default to "skill".
        local row_kind="${inj_kind:-skill}"
        if [[ "$inj_sk" != "$name" ]]; then
            continue
        fi
        if [[ -n "$kind_filter" && "$row_kind" != "$kind_filter" ]]; then
            continue
        fi
        matched_overlays+=("$inj_ov")
        matched_vias+=("$inj_via")
        matched_kinds+=("$row_kind")
    done < <(manifest_injected)

    local n=${#matched_overlays[@]}

    if (( n == 0 )); then
        # Not in injected[] (under the filter, if any). Check stack-mount for
        # a helpful error. The check uses kind_filter when set; otherwise it
        # matches any kind so the user gets the deactivate-pointer message
        # regardless of which YAML section the entry lives in.
        if _eject_entry_is_stack_mounted "$name" "$kind_filter" "$base" "$overlay"; then
            local kdesc=""
            [[ -n "$kind_filter" ]] && kdesc=" ($kind_filter)"
            echo "sciagent eject: '$name'$kdesc is part of stack-mounted role '$base'${overlay:+/$overlay}; use 'sciagent deactivate' instead" >&2
            return 1
        fi
        local notice="not injected: $name"
        [[ -n "$kind_filter" ]] && notice="not injected as $kind_filter: $name"
        echo "$notice — nothing to eject"
        return 0
    fi

    if (( n >= 2 )); then
        # Ambiguity: bare-name eject when two or more kinds carry the same
        # name. Tell the user to disambiguate with the same flag set that
        # inject accepts.
        local kinds_csv
        kinds_csv=$(IFS=,; printf '%s' "${matched_kinds[*]}")
        if (( n == 2 )); then
            echo "error: ambiguous — '$name' is injected as both ${matched_kinds[0]} and ${matched_kinds[1]}. use --skill <name>, --agent <name>, or --command <name>" >&2
        else
            echo "error: ambiguous — '$name' is injected as $kinds_csv. use --skill <name>, --agent <name>, or --command <name>" >&2
        fi
        return 1
    fi

    # Exactly one row matched.
    local target_kind="${matched_kinds[0]}"
    local target_via="${matched_vias[0]}"

    # Dependency guard: a row carrying via="requires:<root>" was auto-mounted to
    # satisfy another skill's closure, not injected on its own. Ejecting it by
    # name would silently break that closure (and the prune below would not
    # re-add it). Point the user at the root instead; orphaned deps whose root
    # is already gone are cleaned by the prune, never reached here.
    if [[ "$target_via" == requires:* ]]; then
        echo "sciagent eject: '$name' was auto-mounted as a dependency (${target_via}); eject that skill instead, or run 'sciagent deactivate'" >&2
        return 1
    fi

    # Remove the dual symlinks for this kind+name.
    _eject_remove_symlinks "$name" "$target_kind"

    # Rebuild the manifest without this entry (matched by skill + kind).
    _eject_drop_injected_entry "$name" "$target_kind" ""

    # Prune any requires-closure deps this skill pulled in that no real root
    # still needs (orphan-prevention counterpart to inject's F4 dep mounting).
    # Safe/idempotent for any kind, but only meaningful after a skill eject.
    [[ "$target_kind" == skill ]] && _eject_prune_orphaned_requires_deps

    # Collapse _injected overlay if it is now empty.
    _eject_maybe_collapse_injected_overlay

    # Re-render AGENTS.md block.
    _eject_rewrite_block || return 1

    echo "ejected: $name ($target_kind)"
    return 0
}

# _eject_by_tag <tag-name>
# Removes all injected rows whose via value equals "tag:<tag-name>".
# Idempotent: if no rows match, prints a notice and exits 0.
_eject_by_tag() {
    local tag_name="$1"
    local via_key="tag:$tag_name"

    # Collect all rows matching this tag.
    local -a to_eject_names=() to_eject_kinds=()
    local inj_ov inj_sk inj_via inj_kind
    while IFS='|' read -r inj_ov inj_sk inj_via inj_kind; do
        [[ -z "$inj_ov" ]] && continue
        if [[ "$inj_via" == "$via_key" ]]; then
            to_eject_names+=("$inj_sk")
            to_eject_kinds+=("${inj_kind:-skill}")
        fi
    done < <(manifest_injected)

    if [[ "${#to_eject_names[@]}" -eq 0 ]]; then
        echo "not injected via tag '$tag_name' — nothing to eject"
        return 0
    fi

    local i
    for (( i=0; i<${#to_eject_names[@]}; i++ )); do
        local nm="${to_eject_names[$i]}" kd="${to_eject_kinds[$i]}"
        _eject_remove_symlinks "$nm" "$kd"
        _eject_drop_injected_entry "$nm" "$kd" "$via_key"
        echo "ejected: $nm (tag:$tag_name)"
    done

    # A tag-injected skill may have pulled requires-closure deps; prune any now
    # orphaned (no real root still needs them).
    _eject_prune_orphaned_requires_deps

    # Collapse _injected overlay if now empty.
    _eject_maybe_collapse_injected_overlay

    # Re-render AGENTS.md block.
    _eject_rewrite_block || return 1

    return 0
}

# _eject_entry_is_stack_mounted <name> <kind_filter> <base> <overlay>
# Returns 0 if the entry appears in a role on the stack. When kind_filter is
# empty, matches ANY kind (skill OR agent OR command). When kind_filter is
# non-empty, restricts to the matching YAML section.
#
# NOTE: inject.sh carries a mirror of this check (_inject_entry_is_stack_mounted).
# The two helpers must stay in sync: both reuse role_load (which emits
# SKILL/AGENT/COMMAND tagged lines) so the YAML-section logic stays in one
# place. The F-1 fix is preserved because role_load strips trailing-comment
# YAML lines (`  - <name>  # description`) before tagging.
_eject_entry_is_stack_mounted() {
    local name="$1" kind_filter="$2" base="$3" overlay="$4"

    _eject_role_contains_entry() {
        local role="$1"
        if [[ -z "$kind_filter" ]]; then
            role_load "$role" 2>/dev/null \
                | grep -Eq "^(SKILL|AGENT|COMMAND) ${name}\$"
        else
            local want
            case "$kind_filter" in
                skill)   want="SKILL" ;;
                agent)   want="AGENT" ;;
                command) want="COMMAND" ;;
                *) return 1 ;;
            esac
            role_load "$role" 2>/dev/null \
                | grep -Eq "^${want} ${name}\$"
        fi
    }

    if _eject_role_contains_entry "$base"; then
        return 0
    fi
    if [[ -n "$overlay" && "$overlay" != "_injected" ]] \
        && _eject_role_contains_entry "$overlay"; then
        return 0
    fi
    return 1
}

# _eject_remove_symlinks <name> <kind>
# Removes the dual symlinks for this kind+name. The .claude/ and .agents/
# paths mirror the layout symlink_create_dual builds for stack-mounted
# entries (see symlinks.sh) and the layout _inject_named_entry creates.
_eject_remove_symlinks() {
    local name="$1" kind="$2"
    local claude_path agents_path
    case "$kind" in
        skill)
            claude_path=".claude/skills/$name"
            agents_path=".agents/skills/$name"
            ;;
        agent)
            claude_path=".claude/agents/${name}.md"
            agents_path=".agents/agents/${name}.md"
            ;;
        command)
            claude_path=".claude/commands/${name}.md"
            agents_path=".agents/commands/${name}.md"
            ;;
        *)
            echo "_eject_remove_symlinks: unknown kind '$kind'" >&2
            return 1 ;;
    esac

    if [[ -L "$claude_path" ]]; then
        rm "$claude_path"
    fi
    if [[ -L "$agents_path" ]]; then
        rm "$agents_path"
    fi
}

# _eject_drop_injected_entry <name> <kind> <via_filter>
# Rewrites the manifest, omitting the injected entry for <name> + <kind>.
# When <via_filter> is non-empty, also matches via (used by --tag eject so a
# named-injected and tag-injected entry of the same name don't collide).
_eject_drop_injected_entry() {
    local name="$1" target_kind="$2" via_filter="$3"

    local old_stack old_hash
    old_stack=$(manifest_stack)
    old_hash=$(manifest_block_hash)

    _manifest_stack_val="$old_stack"
    _manifest_syms=()
    _manifest_injected_overlays=()
    _manifest_injected_skills=()
    _manifest_injected_vias=()
    _manifest_injected_kinds=()

    # Compute the symlink paths to remove for this kind+name.
    local rm_claude rm_agents
    case "$target_kind" in
        skill)
            rm_claude=".claude/skills/$name"
            rm_agents=".agents/skills/$name"
            ;;
        agent)
            rm_claude=".claude/agents/${name}.md"
            rm_agents=".agents/agents/${name}.md"
            ;;
        command)
            rm_claude=".claude/commands/${name}.md"
            rm_agents=".agents/commands/${name}.md"
            ;;
        *)
            rm_claude=""
            rm_agents=""
            ;;
    esac

    # Re-read symlinks, dropping the two entries for this kind+name.
    local p
    while IFS= read -r p; do
        [[ -z "$p" ]] && continue
        if [[ -n "$rm_claude" && "$p" == "$rm_claude" ]]; then
            continue
        fi
        if [[ -n "$rm_agents" && "$p" == "$rm_agents" ]]; then
            continue
        fi
        _manifest_syms+=("$p")
    done < <(manifest_symlinks)

    # Re-read injected entries, dropping the matched one (first match only).
    local inj_ov inj_sk inj_via inj_kind
    local dropped=0
    while IFS='|' read -r inj_ov inj_sk inj_via inj_kind; do
        [[ -z "$inj_ov" ]] && continue
        local row_kind="${inj_kind:-skill}"
        if [[ "$inj_sk" == "$name" && "$row_kind" == "$target_kind" && "$dropped" -eq 0 ]]; then
            if [[ -z "$via_filter" || "${inj_via:-}" == "$via_filter" ]]; then
                dropped=1
                continue  # Skip this row.
            fi
        fi
        _manifest_injected_overlays+=("$inj_ov")
        _manifest_injected_skills+=("$inj_sk")
        _manifest_injected_vias+=("${inj_via:-}")
        _manifest_injected_kinds+=("$row_kind")
    done < <(manifest_injected)

    _manifest_staging=$(mktemp)
    _manifest_write_json "$old_hash"
    mv "$_manifest_staging" "$_MANIFEST_PATH"
    _manifest_staging=""
}

# _eject_prune_orphaned_requires_deps
# Orphan-prevention counterpart to inject's F4 dep mounting. Removes any
# injected skill row tagged via="requires:*" whose dep is no longer reachable
# from any "real root":
#   - every skill declared by the active stack roles (base + overlay), and
#   - every injected skill row whose via is NOT "requires:*" (user-injected roots).
# A requires:* row is itself never a root — it exists only because some root
# pulled it in. We recompute the live requires-closure over all real roots and
# drop any requires:* dep that falls outside it.
_eject_prune_orphaned_requires_deps() {
    # Gather the real roots.
    local -a _roots=()
    local stack base overlay role
    stack=$(manifest_stack)
    base=$(printf '%s\n' "$stack" | awk '{print $1}')
    overlay=$(printf '%s\n' "$stack" | awk '{print $2}')

    for role in "$base" "$overlay"; do
        [[ -z "$role" || "$role" == "_injected" ]] && continue
        local line skill_name
        while IFS= read -r line; do
            skill_name="${line#SKILL }"
            [[ -n "$skill_name" ]] && _roots+=("$skill_name")
        done < <(role_load "$role" 2>/dev/null | grep '^SKILL ')
    done

    # User-injected skill roots: injected rows whose via is not "requires:*".
    local inj_ov inj_sk inj_via inj_kind
    while IFS='|' read -r inj_ov inj_sk inj_via inj_kind; do
        [[ -z "$inj_ov" ]] && continue
        [[ "${inj_kind:-skill}" != "skill" ]] && continue
        [[ "$inj_via" == requires:* ]] && continue
        _roots+=("$inj_sk")
    done < <(manifest_injected)

    # Compute the live needed-set (requires-closure over all real roots). With
    # zero roots the closure is empty, so every requires:* dep is orphaned.
    declare -A _needed=()
    if [[ "${#_roots[@]}" -gt 0 ]]; then
        local n
        while IFS= read -r n; do
            [[ -n "$n" ]] && _needed["$n"]=1
        done < <(skill_resolve_transitive "${_roots[@]}")
    fi

    # Drop every requires:* dep that is no longer needed.
    while IFS='|' read -r inj_ov inj_sk inj_via inj_kind; do
        [[ -z "$inj_ov" ]] && continue
        [[ "${inj_kind:-skill}" != "skill" ]] && continue
        [[ "$inj_via" == requires:* ]] || continue
        if [[ -z "${_needed[$inj_sk]:-}" ]]; then
            _eject_remove_symlinks "$inj_sk" skill
            _eject_drop_injected_entry "$inj_sk" skill ""
        fi
    done < <(manifest_injected)
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
    local inj_ov inj_sk _inj_via _inj_kind
    while IFS='|' read -r inj_ov inj_sk _inj_via _inj_kind; do
        if [[ "$inj_ov" == "_injected" ]]; then
            return 0  # Still has entries; don't collapse.
        fi
    done < <(manifest_injected)

    # No entries remain for _injected: collapse stack to [base].
    manifest_update_stack "$base"
}

# _eject_rewrite_block — rebuild AGENTS.md block from current manifest + role YAMLs.
# Thin wrapper over the shared block_render_and_write recipe (stack.sh).
_eject_rewrite_block() {
    block_render_and_write AGENTS.md || return 1
    craft_render_and_write AGENTS.md
}
