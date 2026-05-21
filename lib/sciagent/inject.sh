# lib/sciagent/inject.sh — sciagent inject <skill>
# Adds one skill on top of the current stack (per architecture §5).
# - Stack [base]            → synthesize implicit overlay "_injected".
# - Stack [base, overlay]   → record skill as injected-into-overlay.
# Idempotent: re-injecting a tracked skill is a no-op.

# shellcheck shell=bash

# Manifest is appended with INJECTED lines:
#   INJECTED <overlay-role> <skill-name>
# The corresponding symlinks are recorded with SYMLINK lines as normal.

cmd_inject() {
    if [[ $# -ne 1 ]]; then
        echo "usage: sciagent inject <skill>" >&2
        return 1
    fi
    local skill="$1"

    if ! manifest_exists; then
        echo "no role active; run \`sciagent activate <role>\` first" >&2
        return 1
    fi

    local skill_src="$SCIAGENT_TOOLKIT/skills/$skill"
    if [[ ! -d "$skill_src" ]]; then
        echo "skill not found: $skill ($skill_src)" >&2
        return 1
    fi

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
    if grep -qxF "INJECTED $target_overlay $skill" "$_MANIFEST_PATH" 2>/dev/null; then
        echo "already injected: $skill (into $target_overlay)"
        return 0
    fi

    # Create the symlinks and append manifest records in place.
    mkdir -p .claude/skills .agents/skills
    local claude_path=".claude/skills/$skill"
    local agents_path=".agents/skills/$skill"
    ln -sfn "$skill_src" "$claude_path"
    ln -sfn "$skill_src" "$agents_path"
    {
        printf 'SYMLINK %s\n' "$claude_path"
        printf 'SYMLINK %s\n' "$agents_path"
        printf 'INJECTED %s %s\n' "$target_overlay" "$skill"
    } >> "$_MANIFEST_PATH"

    # Update the stack line if synthesizing _injected.
    if [[ -z "$overlay" ]]; then
        local tmp
        tmp=$(mktemp)
        awk -v base="$base" '
            /^STACK / { print "STACK " base " _injected"; next }
            { print }
        ' "$_MANIFEST_PATH" > "$tmp"
        mv "$tmp" "$_MANIFEST_PATH"
    fi

    # Re-render the AGENTS.md managed block to reflect injection.
    _inject_rewrite_block

    echo "injected: $skill (into $target_overlay)"
}

# _inject_rewrite_block — rebuild the AGENTS.md block from current manifest +
# role YAMLs, adding an "## Injected (overlay)" subsection per architecture §5.
_inject_rewrite_block() {
    local stack base overlay
    stack=$(manifest_stack)
    base=$(printf '%s\n' "$stack" | awk '{print $1}')
    overlay=$(printf '%s\n' "$stack" | awk '{print $2}')

    # Collect injected skills per overlay from the manifest.
    local -a INJECTED_NAMES=()
    local kind ov name
    while read -r kind ov name; do
        [[ "$kind" == "INJECTED" ]] || continue
        INJECTED_NAMES+=("$name")
    done < "$_MANIFEST_PATH"

    # Walk roles to compute effective tables (mirrors activate.sh logic).
    declare -A SKILLS=() AGENTS_M=() COMMANDS_M=()
    declare -A SKILL_SHADOWS=() AGENT_SHADOWS=() COMMAND_SHADOWS=()
    local -a SKILL_ORDER=() AGENT_ORDER=() COMMAND_ORDER=()
    local OUTPUT_STYLE="" OUTPUT_STYLE_ROLE=""

    local role line kind2 nm
    for role in $base ${overlay:+$overlay}; do
        # _injected has no YAML file; its content is the INJECTED records.
        if [[ "$role" == "_injected" ]]; then
            continue
        fi
        while IFS= read -r line; do
            kind2="${line%% *}"
            nm="${line#* }"
            case "$kind2" in
                SKILL)
                    if [[ -n "${SKILLS[$nm]:-}" ]]; then
                        SKILL_SHADOWS[$nm]="${SKILLS[$nm]}"
                    else
                        SKILL_ORDER+=("$nm")
                    fi
                    SKILLS[$nm]="$role" ;;
                AGENT)
                    if [[ -n "${AGENTS_M[$nm]:-}" ]]; then
                        AGENT_SHADOWS[$nm]="${AGENTS_M[$nm]}"
                    else
                        AGENT_ORDER+=("$nm")
                    fi
                    AGENTS_M[$nm]="$role" ;;
                COMMAND)
                    if [[ -n "${COMMANDS_M[$nm]:-}" ]]; then
                        COMMAND_SHADOWS[$nm]="${COMMANDS_M[$nm]}"
                    else
                        COMMAND_ORDER+=("$nm")
                    fi
                    COMMANDS_M[$nm]="$role" ;;
                OUTPUT_STYLE)
                    OUTPUT_STYLE="$nm"; OUTPUT_STYLE_ROLE="$role" ;;
            esac
        done < <(role_load "$role")
    done

    local body
    body=$(_render_inject_body "$base" "$overlay" "${INJECTED_NAMES[@]+"${INJECTED_NAMES[@]}"}")
    block_write AGENTS.md "$body"

    # Refresh manifest BLOCK_HASH (match block_write's canonicalisation).
    [[ "${body: -1}" == $'\n' ]] || body="${body}"$'\n'
    local hash
    hash=$(printf '%s' "$body" | sha1sum | awk '{print $1}')
    local tmp
    tmp=$(mktemp)
    awk -v h="$hash" '
        /^BLOCK_HASH / { print "BLOCK_HASH " h; next }
        { print }
    ' "$_MANIFEST_PATH" > "$tmp"
    mv "$tmp" "$_MANIFEST_PATH"
}

# _render_inject_body <base> <overlay> [injected...]
# Renders the managed-block body including an Injected subsection.
_render_inject_body() {
    local base="$1"; shift
    local overlay="$1"; shift
    local -a injected=( "$@" )

    {
        printf '# Active roles\n\n'
        printf 'Stack (in order, last-wins on name collisions):\n'
        local desc
        desc=$(role_description "$base")
        printf '1. **%s** — `roles/%s.yaml` — %s\n' "$base" "$base" "$desc"
        if [[ -n "$overlay" ]]; then
            if [[ "$overlay" == "_injected" ]]; then
                printf '2. **_injected** — (synthetic overlay for injected skills)  *(overlay)*\n'
            else
                desc=$(role_description "$overlay")
                printf '2. **%s** — `roles/%s.yaml` — %s  *(overlay)*\n' "$overlay" "$overlay" "$desc"
            fi
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

        if [[ ${#injected[@]} -gt 0 ]]; then
            printf '## Injected (overlay)\n'
            for n in "${injected[@]}"; do
                printf -- '- `%s` — (injected)\n' "$n"
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
