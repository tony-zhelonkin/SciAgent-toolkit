# lib/sciagent/update.sh — sciagent update [--to <ref>] [--no-pin] [--quiet]
#
# Propagation verb (architecture §5 item 3).
# Three steps:
#   1. Re-pin the 01_modules/SciAgent-toolkit submodule (unless --no-pin).
#   2. Re-activate the current stack via re-exec of the dispatcher so the
#      freshly-pinned lib is loaded (not the stale copy already sourced).
#   3. Print a concise summary: submodule delta, CRAFT version, next steps.
#
# Re-exec rationale: after a submodule pin the on-disk lib changes but THIS
# process already sourced the old code.  Rather than trying to re-source
# selectively (fragile: global variables, associative arrays, already-declared
# functions), we exec the dispatcher binary from the NEW toolkit path with
# `update --no-pin`.  The new process starts clean, sources the new lib, and
# runs cmd_update --no-pin, which does the correct re-activate
# in-process.  This is safe because exec replaces the process image; no output
# or state from the current process is lost (we print the update summary
# before exec'ing).
#
# Injected-stack handling: when the manifest stack contains the synthetic
# `_injected` token (written by the retired `sciagent inject`), passing
# `cmd_activate base _injected` fails because `_injected` is not a real role.
# Fix: strip `_injected` from the stack before calling cmd_activate. Legacy
# manifests can still carry the token, so the strip stays for compatibility.
#
# --no-pin mode: skip the submodule step and only re-activate.
# This is useful when the toolkit is not a submodule (PATH install, dev
# symlink, etc.).  In --no-pin mode the re-activate runs in THIS process
# (lib is already current) so no re-exec is needed.
#
# Depends on: block.sh craft.sh symlinks.sh roles.sh
#             collisions.sh validate.sh stack.sh claude_settings.sh activate.sh

# shellcheck shell=bash

# Default only. cmd_update re-derives this per invocation via
# _sciagent_in_repo_toolkit (symlinks.sh), because the container directory is
# NOT always `01_modules` — measured across the fleet: 01_modules, 01_Modules,
# 01_scripts, 01_Scripts. Where the literal did not match, the `-d` test below
# failed, `update` printed "not a submodule — skipped" and then RE-ACTIVATED
# AGAINST A STALE TOOLKIT without re-pinning. Silently, and in ~6 of 23
# consumers.
_UPDATE_SUBMODULE_PATH="01_modules/SciAgent-toolkit"

cmd_update() {
    local _ref="" _no_pin=false _quiet=false

    while (( $# )); do
        case "$1" in
            --to)
                if [[ -z "${2:-}" ]]; then
                    echo "sciagent update: --to requires an argument" >&2
                    return 1
                fi
                _ref="$2"; shift 2 ;;
            --no-pin)  _no_pin=true;  shift ;;
            --quiet)   _quiet=true;   shift ;;
            -h|--help)
                cat <<'EOF'
usage: sciagent update [--to <ref>] [--no-pin] [--quiet]

  --to <ref>   git ref/tag/branch to check out in the submodule
               (default: submodule update --remote to configured branch)
  --no-pin     skip submodule re-pin; only re-activate the current stack
  --quiet      suppress verbose output from the re-activation step
EOF
                return 0 ;;
            *)
                echo "sciagent update: unknown argument '$1'" >&2
                return 1 ;;
        esac
    done

    # Resolve where THIS project's toolkit actually lives before any step uses
    # the path. Falls back to the module-level default when the project ships
    # no in-repo toolkit at all, so the "not a submodule" branch below still
    # reports a sensible path.
    local _discovered
    if _discovered=$(_sciagent_in_repo_toolkit); then
        _UPDATE_SUBMODULE_PATH="${_discovered#./}"
    fi

    # -----------------------------------------------------------------------
    # Step 1: re-pin the submodule (unless --no-pin).
    # -----------------------------------------------------------------------
    local _sub_old="" _sub_new="" _sub_count=0 _sub_status=""

    if [[ "$_no_pin" == "false" ]]; then
        if [[ -d "$_UPDATE_SUBMODULE_PATH" && \
              ( -d "$_UPDATE_SUBMODULE_PATH/.git" || \
                -f "$_UPDATE_SUBMODULE_PATH/.git" ) ]]; then

            # Record the current commit (short SHA).
            _sub_old=$(git -C "$_UPDATE_SUBMODULE_PATH" rev-parse --short HEAD 2>/dev/null \
                       || echo "unknown")

            if [[ -n "$_ref" ]]; then
                # Explicit ref: fetch then checkout.
                if ! git -C "$_UPDATE_SUBMODULE_PATH" fetch --quiet 2>/dev/null; then
                    echo "sciagent update: warning — git fetch in $_UPDATE_SUBMODULE_PATH failed (continuing)" >&2
                fi
                if ! git -C "$_UPDATE_SUBMODULE_PATH" checkout --quiet "$_ref" 2>/dev/null; then
                    echo "sciagent update: error — could not checkout '$_ref' in $_UPDATE_SUBMODULE_PATH" >&2
                    return 1
                fi
            else
                # Default: update to the configured remote tracking branch.
                if ! git submodule update --remote --quiet "$_UPDATE_SUBMODULE_PATH" 2>/dev/null; then
                    echo "sciagent update: warning — 'git submodule update --remote' failed (continuing)" >&2
                fi
            fi

            _sub_new=$(git -C "$_UPDATE_SUBMODULE_PATH" rev-parse --short HEAD 2>/dev/null \
                       || echo "unknown")

            if [[ "$_sub_old" != "$_sub_new" ]]; then
                _sub_count=$(git -C "$_UPDATE_SUBMODULE_PATH" \
                    log --oneline "${_sub_old}..${_sub_new}" 2>/dev/null | wc -l \
                    | tr -d ' ')
                _sub_status="updated ${_sub_old}..${_sub_new} (${_sub_count} commit(s))"
            else
                _sub_status="already up-to-date (${_sub_old})"
            fi

            [[ "$_quiet" == "false" ]] && \
                echo "submodule: $_UPDATE_SUBMODULE_PATH — $_sub_status"

        else
            _sub_status="not a submodule / not found — skipped"
            [[ "$_quiet" == "false" ]] && \
                echo "submodule: $_UPDATE_SUBMODULE_PATH — $_sub_status"
        fi
    else
        _sub_status="skipped (--no-pin)"
        [[ "$_quiet" == "false" ]] && echo "submodule: $_sub_status"
    fi

    # -----------------------------------------------------------------------
    # Step 2: require an active stack.
    # -----------------------------------------------------------------------
    if ! manifest_exists; then
        echo "sciagent update: no active stack — run 'sciagent activate <role>' first" >&2
        return 1
    fi

    local _stack
    _stack=$(manifest_stack) || {
        echo "sciagent update: could not read active stack from manifest" >&2
        return 1
    }
    if [[ -z "$_stack" ]]; then
        echo "sciagent update: no active stack — run 'sciagent activate <role>' first" >&2
        return 1
    fi

    # -----------------------------------------------------------------------
    # Compute the REAL stack: drop the synthetic `_injected` token so that
    # cmd_activate is never called with a non-role argument.
    # e.g. "base _injected" → base="_base" real_overlay=""
    #      "base reviewer"  → base="base"  real_overlay="reviewer"
    #      "base"           → base="base"  real_overlay=""
    # -----------------------------------------------------------------------
    local _base _overlay_raw="" _real_overlay=""
    read -r _base _overlay_raw <<< "$_stack"
    if [[ "$_overlay_raw" != "_injected" ]]; then
        _real_overlay="$_overlay_raw"
    fi
    # The real stack (without _injected) for display purposes.
    local _real_stack="$_base"
    [[ -n "$_real_overlay" ]] && _real_stack="$_base $_real_overlay"

    # -----------------------------------------------------------------------
    # Step 3 (re-activate):
    #   --no-pin path: lib is already current; call cmd_activate in-process,
    #                  (the whole skill catalog is remounted by activate).
    #   re-pin path:   on-disk lib may have changed; exec `update --no-pin`
    #                  on the freshly-loaded dispatcher so the new process
    #                  sources the fresh lib and performs the combined
    #                  re-activation in one place.
    # -----------------------------------------------------------------------

    # Print the summary BEFORE exec so it is visible even when exec replaces us.
    [[ "$_quiet" == "false" ]] && echo "re-activating stack: $_stack"

    if [[ "$_no_pin" == "true" ]]; then
        # In-process re-activate: lib is already sourced and current.
        # Activate the REAL stack (never pass _injected to cmd_activate).
        local _activate_args=("$_base")
        [[ -n "$_real_overlay" ]] && _activate_args+=("$_real_overlay")
        if [[ "$_quiet" == "true" ]]; then
            cmd_activate "${_activate_args[@]}" >/dev/null
        else
            cmd_activate "${_activate_args[@]}"
        fi
        local _rc=$?
        if [[ $_rc -ne 0 ]]; then
            echo "sciagent update: re-activation failed (exit $_rc)" >&2
            return $_rc
        fi

    else
        # Re-exec path: the dispatcher binary in the POTENTIALLY-UPDATED toolkit.
        # The submodule path is relative to $PWD (the project root); resolve the
        # dispatcher from the updated toolkit location.
        local _new_toolkit="$_UPDATE_SUBMODULE_PATH"
        local _dispatcher

        # Prefer the new submodule's bin/sciagent if it exists; fall back to the
        # ambient SCIAGENT_TOOLKIT (useful when submodule was not found/updated).
        if [[ -x "$_new_toolkit/bin/sciagent" ]]; then
            _dispatcher="$_new_toolkit/bin/sciagent"
        else
            _dispatcher="$SCIAGENT_TOOLKIT/bin/sciagent"
        fi

        # Print step 4 summary now — exec will replace us.
        _update_print_summary "$_sub_status" "$_stack"

        # Exec `update --no-pin` on the fresh dispatcher: the new process sources
        # the updated lib and performs the re-activation in-process.
        # This correctly handles _injected stacks (no _injected passed to activate).
        local -a _reexec_args=("update" "--no-pin")
        [[ "$_quiet" == "true" ]] && _reexec_args+=("--quiet")

        exec "$_dispatcher" "${_reexec_args[@]}"
        # exec only returns on error.
        echo "sciagent update: exec failed — $_dispatcher update --no-pin" >&2
        return 1
    fi

    # -----------------------------------------------------------------------
    # Step 4 (report) — only reached in --no-pin path (re-exec path printed
    # summary before exec).
    # -----------------------------------------------------------------------
    _update_print_summary "$_sub_status" "$_stack"
}

# _update_print_summary <sub_status> <stack>
_update_print_summary() {
    local sub_status="$1" stack="$2"
    local cv
    cv=$(craft_version 2>/dev/null || echo "n/a")
    echo ""
    echo "sciagent update — summary"
    echo "  stack:      $stack"
    echo "  submodule:  $sub_status"
    echo "  CRAFT ver:  $cv"
    echo ""
    echo "  Run 'sciagent status' to see the effective table."
    echo "  Run 'sciagent validate' to check toolkit graph health."
    echo "  (sciagent validate --check freshness) will be available in a future release."
}
