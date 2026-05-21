# lib/sciagent/activate.sh — sciagent activate <base> [overlay]
# Computes the effective merged stack (last-wins on name collisions), creates
# dual symlinks, rewrites the AGENTS.md managed block, writes manifest.
# Stack-walking and block-body rendering are delegated to stack.sh.

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

    # Gather resolved entries via stack_walk, preserving insertion order.
    local -a SKILL_ORDER=() AGENT_ORDER=() COMMAND_ORDER=()
    declare -A SKILLS=() AGENTS_M=() COMMANDS_M=()
    local OUTPUT_STYLE="" OUTPUT_STYLE_ROLE=""

    local kind name provider _shadow
    while IFS=$'\t' read -r kind name provider _shadow; do
        case "$kind" in
            SKILL)   SKILLS[$name]="$provider";   SKILL_ORDER+=("$name") ;;
            AGENT)   AGENTS_M[$name]="$provider";  AGENT_ORDER+=("$name") ;;
            COMMAND) COMMANDS_M[$name]="$provider"; COMMAND_ORDER+=("$name") ;;
            STYLE)   OUTPUT_STYLE="$name"; OUTPUT_STYLE_ROLE="$provider" ;;
        esac
    done < <(stack_walk "$base" "$overlay")

    # Start the manifest staging buffer before creating any symlinks.
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
        local style_src="$SCIAGENT_TOOLKIT/system-prompts/cs101.md"
        if [[ -f "$style_src" ]]; then
            symlink_create_dual output-styles "$OUTPUT_STYLE" "$style_src"
        fi
    fi

    # Render and write the managed block.
    local body
    body=$(render_block_body "$base" "$overlay")
    block_write AGENTS.md "$body"
    manifest_finalize "$(block_stored_hash AGENTS.md)"

    # Summary.
    echo "Activated stack: $stack"
    echo "  skills:   ${#SKILL_ORDER[@]}"
    echo "  agents:   ${#AGENT_ORDER[@]}"
    echo "  commands: ${#COMMAND_ORDER[@]}"
    if [[ -n "$OUTPUT_STYLE" ]]; then
        echo "  output_style: $OUTPUT_STYLE ($OUTPUT_STYLE_ROLE)"
    fi
}
