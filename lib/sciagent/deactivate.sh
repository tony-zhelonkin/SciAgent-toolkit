# lib/sciagent/deactivate.sh — sciagent deactivate [<role>]
# No arg: full teardown. One arg: partial.

# shellcheck shell=bash

cmd_deactivate() {
    if [[ $# -eq 0 ]]; then
        if ! manifest_exists; then
            echo "no active stack"
            return 0
        fi
        claude_settings_teardown
        symlink_teardown_all
        block_remove AGENTS.md 2>/dev/null || true
        craft_remove AGENTS.md 2>/dev/null || true
        echo "deactivated"
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
        block_remove AGENTS.md 2>/dev/null || true
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
