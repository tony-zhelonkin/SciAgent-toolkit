# lib/sciagent/activate.sh — sciagent activate <base> [overlay]
# Computes effective merged stack (last-wins on name collisions), creates
# dual symlinks, rewrites the AGENTS.md managed block, writes manifest.

# shellcheck shell=bash

cmd_activate() {
    if [[ $# -lt 1 ]]; then
        echo "usage: sciagent activate <base> [overlay]" >&2
        return 1
    fi
    if [[ $# -gt 2 ]]; then
        echo "Maximum stack depth is 2 (base + overlay)" >&2
        return 1
    fi

    local base="$1"
    local overlay="${2:-}"

    if ! role_exists "$base"; then
        echo "role not found: $base ($(role_path "$base"))" >&2
        return 1
    fi
    if [[ -n "$overlay" ]] && ! role_exists "$overlay"; then
        echo "role not found: $overlay ($(role_path "$overlay"))" >&2
        return 1
    fi

    # Auto-deactivate if a stack is already active.
    if manifest_exists; then
        symlink_teardown_all
        block_remove AGENTS.md 2>/dev/null || true
    fi

    # Walk the stack: base first, then overlay (last-wins).
    # Use parallel arrays "name|source-role" per category.
    declare -A SKILLS=()        # name -> role
    declare -A AGENTS_M=()      # name -> role
    declare -A COMMANDS_M=()    # name -> role
    declare -A SKILL_SHADOWS=()
    declare -A AGENT_SHADOWS=()
    declare -A COMMAND_SHADOWS=()
    local OUTPUT_STYLE=""
    local OUTPUT_STYLE_ROLE=""

    # Preserve insertion order (bash assoc arrays don't).
    local -a SKILL_ORDER=() AGENT_ORDER=() COMMAND_ORDER=()

    local role
    for role in $base ${overlay:+$overlay}; do
        local line kind name
        while IFS= read -r line; do
            kind="${line%% *}"
            name="${line#* }"
            case "$kind" in
                SKILL)
                    if [[ -n "${SKILLS[$name]:-}" ]]; then
                        SKILL_SHADOWS[$name]="${SKILLS[$name]}"
                    else
                        SKILL_ORDER+=("$name")
                    fi
                    SKILLS[$name]="$role"
                    ;;
                AGENT)
                    if [[ -n "${AGENTS_M[$name]:-}" ]]; then
                        AGENT_SHADOWS[$name]="${AGENTS_M[$name]}"
                    else
                        AGENT_ORDER+=("$name")
                    fi
                    AGENTS_M[$name]="$role"
                    ;;
                COMMAND)
                    if [[ -n "${COMMANDS_M[$name]:-}" ]]; then
                        COMMAND_SHADOWS[$name]="${COMMANDS_M[$name]}"
                    else
                        COMMAND_ORDER+=("$name")
                    fi
                    COMMANDS_M[$name]="$role"
                    ;;
                OUTPUT_STYLE)
                    OUTPUT_STYLE="$name"
                    OUTPUT_STYLE_ROLE="$role"
                    ;;
            esac
        done < <(role_load "$role")
    done

    # Start a manifest staging file BEFORE creating any symlinks.
    local stack="$base"
    [[ -n "$overlay" ]] && stack="$base $overlay"
    manifest_begin "$stack"

    # Create symlinks.
    local n
    for n in "${SKILL_ORDER[@]:-}"; do
        [[ -z "$n" ]] && continue
        symlink_create_dual skills "$n" "$SCIAGENT_TOOLKIT/skills/$n"
    done
    for n in "${AGENT_ORDER[@]:-}"; do
        [[ -z "$n" ]] && continue
        symlink_create_dual agents "$n" "$SCIAGENT_TOOLKIT/agents/${n}.md"
    done
    for n in "${COMMAND_ORDER[@]:-}"; do
        [[ -z "$n" ]] && continue
        symlink_create_dual commands "$n" "$SCIAGENT_TOOLKIT/commands/${n}.md"
    done
    if [[ -n "$OUTPUT_STYLE" ]]; then
        # Output styles resolve by `name:` frontmatter, not filename. We
        # symlink every available style file; consumers pick by name.
        # For now, link the canonical cs101.md as the catch-all.
        local style_src="$SCIAGENT_TOOLKIT/output-styles/cs101.md"
        if [[ -f "$style_src" ]]; then
            symlink_create_dual output-styles "$OUTPUT_STYLE" "$style_src"
        fi
    fi

    # Render the managed block.
    local body
    body=$(_render_block_body "$base" "$overlay")

    block_write AGENTS.md "$body"
    local hash
    hash=$(printf '%s' "$body" | sha1sum | awk '{print $1}')
    manifest_finalize "$hash"

    # Summary.
    echo "Activated stack: $stack"
    echo "  skills:   ${#SKILL_ORDER[@]}"
    echo "  agents:   ${#AGENT_ORDER[@]}"
    echo "  commands: ${#COMMAND_ORDER[@]}"
    if [[ -n "$OUTPUT_STYLE" ]]; then
        echo "  output_style: $OUTPUT_STYLE ($OUTPUT_STYLE_ROLE)"
    fi
}

_render_block_body() {
    local base="$1" overlay="$2"
    {
        printf '# Active roles\n\n'
        printf 'Stack (in order, last-wins on name collisions):\n'
        local desc
        desc=$(role_description "$base")
        printf '1. **%s** — `roles/%s.yaml` — %s\n' "$base" "$base" "$desc"
        if [[ -n "$overlay" ]]; then
            desc=$(role_description "$overlay")
            printf '2. **%s** — `roles/%s.yaml` — %s  *(overlay)*\n' "$overlay" "$overlay" "$desc"
        fi
        printf '\n'

        local n shadow
        if [[ ${#SKILL_ORDER[@]} -gt 0 ]]; then
            printf '## Skills (effective)\n'
            for n in "${SKILL_ORDER[@]}"; do
                shadow=""
                [[ -n "${SKILL_SHADOWS[$n]:-}" ]] && shadow=" [shadows ${SKILL_SHADOWS[$n]}]"
                printf -- '- `%s` — (%s)%s\n' "$n" "${SKILLS[$n]}" "$shadow"
            done
            printf '\n'
        fi

        if [[ ${#AGENT_ORDER[@]} -gt 0 ]]; then
            printf '## Sub-agents (effective, Claude-only)\n'
            for n in "${AGENT_ORDER[@]}"; do
                shadow=""
                [[ -n "${AGENT_SHADOWS[$n]:-}" ]] && shadow=" [shadows ${AGENT_SHADOWS[$n]}]"
                printf -- '- `%s` — (%s)%s\n' "$n" "${AGENTS_M[$n]}" "$shadow"
            done
            printf '\n'
        fi

        if [[ ${#COMMAND_ORDER[@]} -gt 0 ]]; then
            printf '## Slash commands (effective, Claude-only)\n'
            for n in "${COMMAND_ORDER[@]}"; do
                shadow=""
                [[ -n "${COMMAND_SHADOWS[$n]:-}" ]] && shadow=" [shadows ${COMMAND_SHADOWS[$n]}]"
                printf -- '- `/%s` — (%s)%s\n' "$n" "${COMMANDS_M[$n]}" "$shadow"
            done
            printf '\n'
        fi

        if [[ -n "$OUTPUT_STYLE" ]]; then
            printf '## Output style\n'
            printf -- '- `%s` — (%s)\n' "$OUTPUT_STYLE" "$OUTPUT_STYLE_ROLE"
        fi
    }
}
