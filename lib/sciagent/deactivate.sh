# lib/sciagent/deactivate.sh — sciagent deactivate [<role>]
# No arg: full teardown. One arg: partial.

# shellcheck shell=bash

# _deactivate_prune_empty_agents_md
# Remove AGENTS.md iff removing our blocks left it completely empty.
#
# block_write's append path creates AGENTS.md when it is absent, so activating
# in a project that had no AGENTS.md leaves a 0-byte file behind once both
# blocks are torn down — residue of exactly the kind this teardown exists to
# eliminate. A user's own AGENTS.md always survives, because it still has their
# bytes in it after block_remove and so is not empty.
#
# No ownership record is needed for this one: at 0 bytes the "is it ours?"
# question has no consequence — an empty AGENTS.md carries no information, so
# neither keeping nor removing it can lose anything. That is the only reason a
# size test is sufficient here rather than a content hash.
_deactivate_prune_empty_agents_md() {
    [[ -f AGENTS.md && ! -s AGENTS.md ]] && rm -f AGENTS.md
    return 0
}

cmd_deactivate() {
    case "${1:-}" in
        -h|--help)
            cat <<'USAGE'
sciagent deactivate [<role>]
  No argument: full teardown of the active stack (all mounted symlinks,
  the ROLES/CRAFT managed blocks, and the manifest) AND of the
  project-level artifacts activate ensures unconditionally: .claude/
  statusline.sh, .claude/hooks/*.sh, the settings.json key-backfill, and
  the CLAUDE.md @AGENTS.md import. Each is reversed only if unchanged
  since sciagent created/modified it; anything edited or replaced since
  is left in place with a warning on stderr, never silently discarded.
  <role>: partial teardown — removing the base implies removing the
  overlay too (same full reversal as above); removing the overlay
  re-activates solo-base (stack-only; the project-level artifacts above
  are untouched, since a stack remains active).
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
        # Full teardown also reverses the project-level artifacts activate
        # ensures unconditionally (statusline.sh, the hook bodies, the
        # settings.json key-backfill, the CLAUDE.md shim) — these are not
        # stack-specific, so they are only reversed here (a full deactivate),
        # never on activate.sh's re-activation preamble. Each is a no-op if we
        # never created/modified it, or if it has already been reversed.
        claude_settings_teardown_project_artifacts
        claude_md_shim_teardown
        _deactivate_prune_empty_agents_md
        # Best-effort: .claude/ and .sciagent/ themselves, now that every
        # subpath we own has been removed above. rmdir (never rm -r) — a
        # non-empty directory means real content remains (ours or the
        # user's) and is left in place by construction.
        rmdir .claude 2>/dev/null || true
        rmdir .sciagent 2>/dev/null || true
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
        # See the no-arg branch above: this is also a full teardown (removing
        # the base implies the overlay too), so the project-level artifacts
        # get reversed here as well.
        claude_settings_teardown_project_artifacts
        claude_md_shim_teardown
        _deactivate_prune_empty_agents_md
        rmdir .claude 2>/dev/null || true
        rmdir .sciagent 2>/dev/null || true
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
