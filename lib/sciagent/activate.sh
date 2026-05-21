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

    # Auto-deactivate if a stack is already active. claude_settings_teardown
    # must run BEFORE symlink_teardown_all rmdirs .sciagent (the state file
    # lives there) and BEFORE the new symlinks/block land so that the
    # settings.local.json revert is observable for the duration of the
    # re-activation rather than racing against the new apply.
    if manifest_exists; then
        claude_settings_teardown
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

    # Resolve output_style → system-prompts/<file>.md by frontmatter name.
    # Validate before any filesystem mutation so a drifted role spec leaves
    # the project untouched (no half-state).
    local STYLE_SRC=""
    if [[ -n "$OUTPUT_STYLE" ]]; then
        STYLE_SRC=$(system_prompt_path "$OUTPUT_STYLE") || true
        if [[ -z "$STYLE_SRC" ]]; then
            echo "sciagent: role '$OUTPUT_STYLE_ROLE' requests output_style '$OUTPUT_STYLE'," >&2
            echo "  but no file in system-prompts/ has frontmatter 'name: $OUTPUT_STYLE'." >&2
            echo "  Available styles (frontmatter name → file):" >&2
            local pname pfile
            while IFS=$'\t' read -r pname pfile; do
                printf '    %s\t(%s)\n' "$pname" "$pfile" >&2
            done < <(system_prompt_inventory)
            return 1
        fi
    fi

    # Start the manifest staging buffer before creating any symlinks.
    local stack="$base"
    [[ -n "$overlay" ]] && stack="$base $overlay"
    manifest_begin "$stack"

    # Create symlinks. Canonical source paths are resolved via
    # resolve_canonical so that subfolders under agents/ and commands/ are
    # transparent to the consumer (symlinks stay flat in .claude/ and .agents/).
    local n src
    for n in "${SKILL_ORDER[@]:-}"; do
        [[ -z "$n" ]] && continue
        src=$(resolve_canonical skills "$n") || return 1
        symlink_create_dual skills "$n" "$src"
    done
    for n in "${AGENT_ORDER[@]:-}"; do
        [[ -z "$n" ]] && continue
        src=$(resolve_canonical agents "$n") || return 1
        symlink_create_dual agents "$n" "$src"
    done
    for n in "${COMMAND_ORDER[@]:-}"; do
        [[ -z "$n" ]] && continue
        src=$(resolve_canonical commands "$n") || return 1
        symlink_create_dual commands "$n" "$src"
    done
    local STYLE_APPLIED_TAG=""
    if [[ -n "$OUTPUT_STYLE" ]]; then
        symlink_create_dual output-styles "$OUTPUT_STYLE" "$STYLE_SRC"
        # Make the style Claude's active one by setting outputStyle in
        # .claude/settings.local.json. Symlinking alone makes the file
        # visible but does not select it as the active style.
        STYLE_APPLIED_TAG=$(claude_settings_apply "$OUTPUT_STYLE")
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
        echo "  output_style: $OUTPUT_STYLE ($OUTPUT_STYLE_ROLE) [settings.local.json: $STYLE_APPLIED_TAG]"
    fi
}
