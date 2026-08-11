# lib/sciagent/deactivate.sh — sciagent deactivate [<role>]
# No arg: full teardown. One arg: partial.

# shellcheck shell=bash

cmd_deactivate() {
    if [[ $# -eq 0 ]]; then
        # Full teardown does NOT require a manifest. Ownership is derived from
        # each link's target, so a lost/deleted manifest.json must not strand
        # the mounts — recovering that state is the point of target-based
        # teardown. Only the partial path below genuinely needs the manifest,
        # because it has to know the stack to decide what to keep.
        local had_manifest=0
        manifest_exists && had_manifest=1

        claude_settings_teardown
        symlink_teardown_all
        # ROLES explicitly: full teardown removes the ROLES block sciagent
        # itself owns; CRAFT is a separate id, torn down via craft_remove.
        block_remove AGENTS.md ROLES 2>/dev/null || true
        craft_remove AGENTS.md 2>/dev/null || true
        if [[ "$had_manifest" -eq 1 ]]; then
            echo "deactivated"
        elif [[ "${_SCIAGENT_TEARDOWN_COUNT:-0}" -gt 0 ]]; then
            # Orphaned mounts with no manifest: say so, because the count is
            # evidence the manifest was lost rather than never written.
            echo "deactivated ($_SCIAGENT_TEARDOWN_COUNT orphaned mounts removed; no manifest present)"
        else
            echo "no active stack"
        fi
        return 0
    fi

    if [[ $# -gt 1 ]]; then
        echo "usage: sciagent deactivate [<role>]" >&2
        return 1
    fi

    local target="$1"
    if ! manifest_exists; then
        echo "no active stack" >&2
        return 1
    fi

    local stack
    stack=$(manifest_stack)
    local base overlay
    base=$(printf '%s\n' "$stack" | awk '{print $1}')
    overlay=$(printf '%s\n' "$stack" | awk '{print $2}')

    if [[ "$target" == "$base" ]]; then
        # Removing the base implies removing the overlay too.
        claude_settings_teardown
        symlink_teardown_all
        # ROLES explicitly: full teardown removes the ROLES block sciagent
        # itself owns; CRAFT is a separate id, torn down via craft_remove.
        block_remove AGENTS.md ROLES 2>/dev/null || true
        craft_remove AGENTS.md 2>/dev/null || true
        echo "deactivated (removed base implies overlay too)"
        return 0
    fi

    if [[ "$target" == "$overlay" ]]; then
        # Re-activate solo-base. cmd_activate handles teardown-first.
        cmd_activate "$base"
        return 0
    fi

    echo "role '$target' not in active stack ($stack)" >&2
    return 1
}
