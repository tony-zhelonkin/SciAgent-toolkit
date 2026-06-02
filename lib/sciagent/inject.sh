# lib/sciagent/inject.sh — sciagent inject <name> [--skill|--agent|--command]
#                                           | --tag <name>
# Adds one skill, sub-agent, or slash-command on top of the current stack
# For tag form, expands to every skill carrying the tag.
# - Stack [base]            → synthesize implicit overlay "_injected".
# - Stack [base, overlay]   → record entry as injected-into-overlay.
# Idempotent: re-injecting a tracked entry is a no-op.
#
# Forms:
#   sciagent inject <name>               — auto-detect kind across the three
#                                          canonical dirs (skills/, agents/,
#                                          commands/). Hard-fails if <name>
#                                          resolves in two or more namespaces;
#                                          use an explicit kind flag to
#                                          disambiguate. Hard-fails if <name>
#                                          resolves in none.
#   sciagent inject --skill <name>       — force skill resolution
#   sciagent inject --agent <name>       — force agent resolution
#   sciagent inject --command <name>     — force command resolution
#   sciagent inject --tag <name>         — inject all skills tagged <name>
#                                          (via: "tag:<name>" per entry).
#                                          Tags remain skill-only.
#
# When injecting a non-skill entry whose name also matches a skill, a one-line
# stderr note advertises the companion skill; the skill is NOT auto-mounted.
#
# Cross-namespace collision check (ADR-5.1 C): before any write, inject warns
# on stderr if <name> collides across ≥2 toolkit namespaces. `--force` silences
# the warning; SCIAGENT_STRICT_COLLISIONS=1 upgrades it to a hard fail;
# allowlisted overlaps are downgraded to a single informational note.
#
# Block-body rendering is delegated to stack.sh:render_block_body.

# shellcheck shell=bash

cmd_inject() {
    if [[ $# -eq 0 ]]; then
        echo "usage: sciagent inject <name> | --skill <name> | --agent <name> | --command <name> | --tag <name>" >&2
        return 1
    fi

    if ! manifest_exists; then
        echo "no role active; run \`sciagent activate <role>\` first" >&2
        return 1
    fi

    # Strip the (positional-independent) --force flag before form dispatch.
    # --force silences the inline cross-namespace collision warning. Exported
    # for _inject_named_entry / _inject_handle_collision via _INJECT_FORCE.
    _INJECT_FORCE=0
    local -a _args=()
    local a
    for a in "$@"; do
        if [[ "$a" == "--force" ]]; then
            _INJECT_FORCE=1
        else
            _args+=("$a")
        fi
    done
    set -- "${_args[@]+"${_args[@]}"}"

    if [[ $# -eq 0 ]]; then
        echo "usage: sciagent inject <name> | --skill <name> | --agent <name> | --command <name> | --tag <name>" >&2
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
        --skill|--agent|--command)
            local flag="$1" name="${2:-}"
            if [[ -z "$name" ]]; then
                echo "usage: sciagent inject $flag <name>" >&2
                return 1
            fi
            local kind="${flag#--}"
            _inject_named_entry "$name" "$kind" ""
            return $?
            ;;
        --*)
            echo "sciagent inject: unknown option '$1'" >&2
            echo "usage: sciagent inject <name> | --skill <name> | --agent <name> | --command <name> | --tag <name>" >&2
            return 1
            ;;
        *)
            if [[ $# -ne 1 ]]; then
                echo "usage: sciagent inject <name> | --skill <name> | --agent <name> | --command <name> | --tag <name>" >&2
                return 1
            fi
            # Auto-detect kind: search the three canonical dirs and hard-fail
            # on ambiguity (two or more matches) or zero matches.
            local detected_kind
            detected_kind=$(_inject_detect_kind "$1") || return 1
            _inject_named_entry "$1" "$detected_kind" ""
            return $?
            ;;
    esac
}

# _inject_detect_kind <name>
# Walks the three canonical dirs and prints the resolved kind on stdout.
# Hard-fails (exit 1, stderr message) when:
#   - <name> resolves in zero namespaces ("not found")
#   - <name> resolves in two or more namespaces ("ambiguous")
# Reuses the same enumeration logic status.sh uses for agents/commands
# (recursive walk so subfolders under agents/ and commands/ are transparent).
_inject_detect_kind() {
    local name="$1"
    local tk_root="${SCIAGENT_TOOLKIT:-}"
    if [[ -z "$tk_root" ]]; then
        echo "sciagent inject: SCIAGENT_TOOLKIT not set" >&2
        return 1
    fi

    local -a found_kinds=()

    # Skill: skills/<name>/SKILL.md (flat — skills are never nested).
    if [[ -d "$tk_root/skills/$name" && -f "$tk_root/skills/$name/SKILL.md" ]]; then
        found_kinds+=("skill")
    fi

    # Agent: agents/**/<name>.md (recursive search, hidden dirs skipped).
    if find "$tk_root/agents" -type f -name "${name}.md" -not -path '*/.*' \
        2>/dev/null | grep -q .; then
        found_kinds+=("agent")
    fi

    # Command: commands/**/<name>.md (recursive — commands often live under
    # nested subdirs like commands/architect/<name>.md).
    if find "$tk_root/commands" -type f -name "${name}.md" -not -path '*/.*' \
        2>/dev/null | grep -q .; then
        found_kinds+=("command")
    fi

    local n=${#found_kinds[@]}
    if (( n == 0 )); then
        echo "sciagent inject: '$name' not found as skill, agent, or command" >&2
        return 1
    fi
    if (( n >= 2 )); then
        # Ambiguity: hard-fail and point at the explicit flags.
        # Message names the two kinds for the most common 2-namespace case;
        # 3-namespace overlap is rare but still listed in full.
        local kind1="${found_kinds[0]}" kind2="${found_kinds[1]}"
        if (( n == 2 )); then
            echo "sciagent inject: ambiguous — '$name' exists as both $kind1 and $kind2. use --skill <name>, --agent <name>, or --command <name>" >&2
        else
            # 3+ namespaces. Render the list as "A, B, and C" rather than the
            # space-separated "${found_kinds[*]}" — the latter prints as
            # "skill agent command" and reads as a single garbled noun.
            local _last_idx=$(( n - 1 ))
            local _i _joined=""
            for (( _i=0; _i<_last_idx; _i++ )); do
                _joined+="${found_kinds[$_i]}, "
            done
            _joined+="and ${found_kinds[$_last_idx]}"
            echo "sciagent inject: ambiguous — '$name' exists as $_joined. use --skill <name>, --agent <name>, or --command <name>" >&2
        fi
        return 1
    fi
    printf '%s\n' "${found_kinds[0]}"
    return 0
}

# _inject_named_entry <name> <kind> <via>
# Injects a single skill, agent, or command. Records the kind + via in the
# manifest. Returns 0 on success or idempotent no-op.
_inject_named_entry() {
    local name="$1" kind="$2" via="$3"

    # Resolve canonical source for this kind.
    local src
    case "$kind" in
        skill)   src=$(resolve_canonical skills   "$name") || return 1 ;;
        agent)   src=$(resolve_canonical agents   "$name") || return 1 ;;
        command) src=$(resolve_canonical commands "$name") || return 1 ;;
        *)
            echo "sciagent inject: internal error — unknown kind '$kind'" >&2
            return 1 ;;
    esac

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

    # Stack-mount guard: if the entry is already mounted via the active role
    # stack (kind-aware — only matches entries declared under the matching
    # YAML key in the role file), the desired end-state is already satisfied.
    # Refusing-but-exit-0 keeps `inject` idempotent for scripts.
    if _inject_entry_is_stack_mounted "$name" "$kind" "$base" "$overlay"; then
        local where="$base"
        [[ -n "$overlay" && "$overlay" != "_injected" ]] && where="$base/$overlay"
        echo "already mounted via stack-role '$where': $name — nothing to inject" >&2
        return 0
    fi

    # Idempotency: already injected? Match on (overlay, name, kind).
    local inj_ov inj_sk _inj_via inj_kind
    while IFS='|' read -r inj_ov inj_sk _inj_via inj_kind; do
        if [[ "$inj_ov" == "$target_overlay" \
            && "$inj_sk" == "$name" \
            && "${inj_kind:-skill}" == "$kind" ]]; then
            echo "already injected: $name (into $target_overlay)"
            return 0
        fi
    done < <(manifest_injected)

    # Companion-skill note: when injecting a non-skill entry, advertise (but
    # do NOT auto-mount) any skill of the same name. Surfaces the "did you
    # also mean the skill?" discoverability cue without taking the action.
    if [[ "$kind" != "skill" ]]; then
        local tk_root="${SCIAGENT_TOOLKIT:-}"
        if [[ -n "$tk_root" \
            && -d "$tk_root/skills/$name" \
            && -f "$tk_root/skills/$name/SKILL.md" ]]; then
            echo "note: companion skill '$name' available — \`inject --skill $name\` to add" >&2
        fi
    fi

    # Inline cross-namespace collision check. Read-only; runs before the first
    # side-effecting line (mkdir), so a hard-fail policy leaves nothing written.
    local _coll_csv
    if _coll_csv=$(collisions_for_name "$name"); then
        _inject_handle_collision "$name" "$kind" "$_coll_csv" || return $?
    fi

    # Create symlinks per kind. The canonical layout mirrors what activate.sh
    # builds for stack-mounted entries (see symlinks.sh:symlink_create_dual).
    local claude_path agents_path claude_dir agents_dir
    case "$kind" in
        skill)
            claude_dir=.claude/skills;   agents_dir=.agents/skills
            claude_path=".claude/skills/$name"
            agents_path=".agents/skills/$name"
            ;;
        agent)
            claude_dir=.claude/agents;   agents_dir=.agents/agents
            claude_path=".claude/agents/${name}.md"
            agents_path=".agents/agents/${name}.md"
            ;;
        command)
            claude_dir=.claude/commands; agents_dir=.agents/commands
            claude_path=".claude/commands/${name}.md"
            agents_path=".agents/commands/${name}.md"
            ;;
    esac
    mkdir -p "$claude_dir" "$agents_dir" || {
        echo "sciagent inject: failed to create symlink dirs for '$name'" >&2
        return 1
    }
    ln -sfn "$src" "$claude_path" || {
        echo "sciagent inject: failed to symlink $claude_path" >&2
        return 1
    }
    ln -sfn "$src" "$agents_path" || {
        echo "sciagent inject: failed to symlink $agents_path" >&2
        return 1
    }

    # Persist new symlinks and injected entry (with via + kind).
    manifest_append_inject "$claude_path" "$agents_path" "$target_overlay" "$name" "$via" "$kind" || {
        echo "sciagent inject: failed to record '$name' in manifest" >&2
        return 1
    }

    # Update the stack line if synthesizing _injected.
    if [[ -z "$overlay" ]]; then
        manifest_update_stack "$base _injected" || {
            echo "sciagent inject: failed to update stack for '$name'" >&2
            return 1
        }
    fi

    # Re-render the AGENTS.md managed block to reflect injection.
    block_render_and_write AGENTS.md || return 1

    echo "injected: $name (into $target_overlay)"
    return 0
}

# _inject_handle_collision <name> <kind> <kinds-csv>
# Policy for a cross-namespace collision detected at inject time (ADR-5.1 C):
#   - allowlisted overlap  → single informational note (or silence under --force)
#   - SCIAGENT_STRICT_COLLISIONS=1 → hard-fail (return 1)
#   - --force (_INJECT_FORCE=1)     → silent, proceed
#   - default                       → warn on stderr, proceed (return 0)
_inject_handle_collision() {
    local name="$1" kind="$2" csv="$3"

    if collisions_is_allowlisted "$name" "$csv"; then
        [[ "${_INJECT_FORCE:-0}" == "1" ]] && return 0
        echo "note: '$name' is an approved family overlap across $csv (per tests/collision-allowlist.txt)" >&2
        return 0
    fi

    if [[ "${SCIAGENT_STRICT_COLLISIONS:-0}" == "1" ]]; then
        echo "sciagent inject: '$name' collides across $csv — disambiguation (/foo vs @foo vs foo) is load-bearing; refusing under SCIAGENT_STRICT_COLLISIONS=1" >&2
        return 1
    fi

    [[ "${_INJECT_FORCE:-0}" == "1" ]] && return 0

    echo "sciagent inject: '$name' now collides across $csv — disambiguation (/foo vs @foo vs foo) is load-bearing" >&2
    return 0
}

# _inject_by_tag <tag-name>
# Reads tags.yaml; asserts the tag exists in vocabulary; collects every skill
# whose metadata.tags: contains <tag-name>; injects each with via:"tag:<name>".
# Tags remain skill-only (per kickoff.md §9, ADR-001 mechanism).
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
            next
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
        _inject_named_entry "$sk" "skill" "$via"
        local rc=$?
        if [[ "$rc" -eq 0 ]]; then
            (( injected_count++ )) || true
        fi
    done

    return 0
}

# _inject_entry_is_stack_mounted <name> <kind> <base> <overlay>
# Returns 0 if the entry (under the matching YAML key for its kind) appears in
# any role on the stack. Kind-aware: a skill named "foo" does not collide with
# a command named "foo" when checking for stack-mount. Reuses role_load (which
# emits SKILL/AGENT/COMMAND tagged lines) so the YAML-section logic stays
# in one place. The `_injected` synthetic overlay has no role file; skipped.
_inject_entry_is_stack_mounted() {
    local name="$1" kind="$2" base="$3" overlay="$4"

    local want=""
    case "$kind" in
        skill)   want="SKILL" ;;
        agent)   want="AGENT" ;;
        command) want="COMMAND" ;;
        *) return 1 ;;
    esac

    _inject_role_contains_entry() {
        local role="$1"
        # role_load emits `<KIND> <name>` lines; match the exact kind+name.
        role_load "$role" 2>/dev/null \
            | grep -Eq "^$want ${name}\$"
    }

    if _inject_role_contains_entry "$base"; then
        return 0
    fi
    if [[ -n "$overlay" && "$overlay" != "_injected" ]] \
        && _inject_role_contains_entry "$overlay"; then
        return 0
    fi
    return 1
}
