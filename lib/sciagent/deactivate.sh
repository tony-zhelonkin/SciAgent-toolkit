# lib/sciagent/deactivate.sh — sciagent deactivate [<role>]
# No arg: full teardown. One arg: partial.

# shellcheck shell=bash

cmd_deactivate() {
    case "${1:-}" in
        -h|--help)
            cat <<'USAGE'
sciagent deactivate [<role>]
  No argument: full teardown of the active stack (all mounted symlinks,
  the ROLES/CRAFT managed blocks, and the manifest).
  <role>: partial teardown — removing the base implies removing the
  overlay too; removing the overlay re-activates solo-base.
USAGE
            return 0 ;;
    esac

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
        #
        # This calls cmd_activate IN-PROCESS, the same shape as update.sh's
        # bypass of the dispatcher-level guard: bin/sciagent's
        # _guard_toolkit_locality already ran for the top-level `deactivate`
        # verb, but `deactivate` is deliberately NOT in MUTATING_VERB (its
        # delete-only paths are safe against any toolkit), so that check said
        # nothing about the toolkit this re-activation is about to mount
        # against. Check it here, or `SCIAGENT_TOOLKIT=<external> sciagent
        # deactivate <overlay>` would silently re-mount solo-base against the
        # wrong checkout — the exact escape the guard exists to prevent.
        if [[ "${SCIAGENT_ALLOW_EXTERNAL_TOOLKIT:-0}" != "1" ]] && ! _sciagent_toolkit_locality_ok; then
            local _in_repo_real _active_real
            _in_repo_real="$(realpath "./01_modules/SciAgent-toolkit" 2>/dev/null || printf '%s' "./01_modules/SciAgent-toolkit")"
            _active_real="$(realpath "${SCIAGENT_TOOLKIT:-}" 2>/dev/null || printf '%s' "${SCIAGENT_TOOLKIT:-}")"
            echo "sciagent: refusing to deactivate '$overlay' against an external toolkit" >&2
            echo "  active toolkit : $_active_real" >&2
            echo "  in-repo toolkit: $_in_repo_real" >&2
            echo "  removing '$overlay' re-activates solo-'$base', which would mount it against" >&2
            echo "  the wrong toolkit and escape the submodule pin. Fix by one of:" >&2
            echo "    - run ./01_modules/SciAgent-toolkit/bin/sciagent deactivate $overlay" >&2
            echo "    - export SCIAGENT_TOOLKIT=$_in_repo_real" >&2
            echo "    - set SCIAGENT_ALLOW_EXTERNAL_TOOLKIT=1 to override" >&2
            return 1
        fi
        cmd_activate "$base"
        return 0
    fi

    echo "role '$target' not in active stack ($stack)" >&2
    return 1
}
