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

    # Phase A: gather direct entries via stack_walk, preserving insertion order.
    # No mutation yet — all validation/resolution must succeed before we touch
    # the filesystem (per ADR-0002 §4.6).
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

    # Phase B: transitive `requires:` resolution. For each direct skill, fold
    # its closure into SKILL_ORDER; new entries get provider=":requires:<parent>"
    # so the manifest is self-describing (per ADR-0002 Decision 1).
    # Resolver failures (cycles, missing targets) abort before any mutation.
    local direct
    for direct in "${SKILL_ORDER[@]+"${SKILL_ORDER[@]}"}"; do
        # Skip skills with no SKILL.md frontmatter (e.g. test fixtures).
        if ! skill_frontmatter_path "$direct" >/dev/null 2>&1; then
            continue
        fi
        # Capture closure via command-substitution so the resolver's exit
        # status propagates (process substitution + `done` would swallow it).
        local closure_out
        if ! closure_out=$(skill_resolve_transitive "$direct"); then
            echo "sciagent: aborting activation; requires resolution failed for '$direct'" >&2
            return 1
        fi
        local dep
        while IFS= read -r dep; do
            [[ -z "$dep" ]] && continue
            # Skip self and skills already in the effective set.
            [[ "$dep" == "$direct" ]] && continue
            [[ -n "${SKILLS[$dep]:-}" ]] && continue
            SKILLS[$dep]=":requires:$direct"
            SKILL_ORDER+=("$dep")
        done <<< "$closure_out"
    done

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

    # Auto-deactivate if a stack is already active. claude_settings_teardown
    # must run BEFORE symlink_teardown_all rmdirs .sciagent (the state file
    # lives there) and BEFORE the new symlinks/block land so that the
    # settings.local.json revert is observable for the duration of the
    # re-activation rather than racing against the new apply.
    # Deferred until AFTER skill_resolve_transitive succeeds, so a failed
    # resolution leaves the previous stack intact.
    if manifest_exists; then
        claude_settings_teardown
        symlink_teardown_all
        block_remove AGENTS.md 2>/dev/null || true
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

    # Render and write the managed block. Inherited (`:requires:`) skills are
    # passed via the SCIAGENT_INHERITED env var so render_block_body can put
    # them under a dedicated subsection (per ADR-0002 §4.4).
    local SCIAGENT_INHERITED=""
    local sk
    for sk in "${SKILL_ORDER[@]+"${SKILL_ORDER[@]}"}"; do
        if [[ "${SKILLS[$sk]:-}" == :requires:* ]]; then
            SCIAGENT_INHERITED+="${sk}=${SKILLS[$sk]#:requires:};"
        fi
    done
    export SCIAGENT_INHERITED
    local body
    body=$(render_block_body "$base" "$overlay")
    unset SCIAGENT_INHERITED
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
