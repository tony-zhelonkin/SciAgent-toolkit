# lib/sciagent/inject.sh — sciagent inject <skill>
# Adds one skill on top of the current stack (per architecture §5).
# - Stack [base]            → synthesize implicit overlay "_injected".
# - Stack [base, overlay]   → record skill as injected-into-overlay.
# Idempotent: re-injecting a tracked skill is a no-op.
#
# Block-body rendering is delegated to stack.sh:render_block_body.

# shellcheck shell=bash

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
    local inj_line
    while IFS=' ' read -r inj_line; do
        if [[ "$inj_line" == "$target_overlay $skill" ]]; then
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

    # Persist new symlinks and injected entry.
    manifest_append_inject "$claude_path" "$agents_path" "$target_overlay" "$skill"

    # Update the stack line if synthesizing _injected.
    if [[ -z "$overlay" ]]; then
        manifest_update_stack "$base _injected"
        overlay="_injected"
    fi

    # Re-render the AGENTS.md managed block to reflect injection.
    _inject_rewrite_block

    echo "injected: $skill (into $target_overlay)"
}

# _inject_rewrite_block — rebuild AGENTS.md block from current manifest + role YAMLs.
_inject_rewrite_block() {
    local stack base overlay
    stack=$(manifest_stack)
    base=$(printf '%s\n' "$stack" | awk '{print $1}')
    overlay=$(printf '%s\n' "$stack" | awk '{print $2}')

    # Collect injected skill names from manifest.
    local -a INJECTED_NAMES=()
    local ov nm
    while read -r ov nm; do
        [[ -n "$ov" ]] && INJECTED_NAMES+=("$nm")
    done < <(manifest_injected)

    local body
    body=$(render_block_body "$base" "$overlay" "${INJECTED_NAMES[@]+"${INJECTED_NAMES[@]}"}")
    block_write AGENTS.md "$body"

    # Refresh manifest BLOCK_HASH.
    manifest_update_block_hash "$(block_stored_hash AGENTS.md)"
}
