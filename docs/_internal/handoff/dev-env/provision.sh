# lib/sciagent/provision.sh — sciagent provision [flags]
#
# User-level / global cross-harness provisioning — an OPT-IN PERSONAL BOOTSTRAP
# (ADR-D4), not a tier of normal activation. It is the one verb that writes
# outside a project directory, so nothing may call it implicitly: no project
# verb invokes it, and it is not part of what the toolkit distributes (ADR-D3 —
# packaging must not select or configure harnesses). See docs/architecture.md §1
# and docs/proposals/2026-08-11-offline-distribution/50_ADRs.md.
#
# Seeds the baseline every project inherits, idempotently and non-clobberingly:
#   --context   drop the shared global AGENTS.md-class context file into each
#               detected harness's global root (block.sh SCIAGENT:CONTEXT block).
#   --settings  seed power-user defaults. claude is implemented fully (folds in
#               claude_settings_ensure_user_defaults + _statusline); the other
#               harnesses' settings adapters are Tier 2 — logged honestly and
#               skipped, never faked.
#
# This is the verb a devcontainer runs once at create time. It does NOT activate
# a role (that is per-project) and never calls activate/validate.
#
# Depends: harness.sh (detection + path map), block.sh (managed-block writer),
# claude_settings.sh (user-scope seeders).

# shellcheck shell=bash

_provision_usage() {
    cat <<EOF
sciagent provision — seed user-level / global context + settings per harness

Usage:
  sciagent provision [--harness claude,pi,codex,agy,opencode|all]
                     [--user] [--context] [--settings] [--dry-run]

  --harness <csv|all>  Target harnesses (default: all detected). Requested-but-
                       absent harnesses are skipped with a warning.
  --user               Seed user/global scope (default; the only scope in Tier 1).
  --context            Drop the global AGENTS.md-class context file.
  --settings           Seed power-user defaults (claude settings.json + statusline).
  --dry-run            Print what would be written, touch nothing.

  With neither --context nor --settings, both are done.
EOF
}

cmd_provision() {
    local harness_arg="all"
    local do_user=1        # --user default on (only scope in Tier 1)
    local do_context=0 do_settings=0
    local dry_run=0

    while [[ $# -gt 0 ]]; do
        case "$1" in
            --harness)
                harness_arg="${2:-}"
                if [[ -z "$harness_arg" ]]; then
                    echo "usage: sciagent provision --harness <csv|all>" >&2
                    return 1
                fi
                shift 2 ;;
            --harness=*) harness_arg="${1#*=}"; shift ;;
            --user)      do_user=1; shift ;;
            --context)   do_context=1; shift ;;
            --settings)  do_settings=1; shift ;;
            --dry-run)   dry_run=1; shift ;;
            -h|--help)   _provision_usage; return 0 ;;
            *)
                echo "sciagent provision: unknown flag '$1'" >&2
                _provision_usage >&2
                return 1 ;;
        esac
    done

    # Default: neither surface named → do both.
    if (( do_context == 0 && do_settings == 0 )); then
        do_context=1; do_settings=1
    fi

    # Detected harnesses (pure probe).
    local -a detected=()
    local h
    while IFS= read -r h; do
        [[ -n "$h" ]] && detected+=("$h")
    done < <(harness_detect)

    # Resolve the requested target set.
    local -a requested=()
    if [[ "$harness_arg" == "all" ]]; then
        requested=("${detected[@]+"${detected[@]}"}")
    else
        local -a names=()
        local _ifs="$IFS"
        IFS=','; read -r -a names <<< "$harness_arg"; IFS="$_ifs"
        local n
        for n in "${names[@]+"${names[@]}"}"; do
            n="${n// /}"
            [[ -z "$n" ]] && continue
            if ! harness_bin "$n" >/dev/null 2>&1; then
                echo "sciagent provision: unknown harness '$n' (known: $_HARNESS_ALL)" >&2
                continue
            fi
            if harness_is_present "$n"; then
                requested+=("$n")
            else
                echo "sciagent provision: harness '$n' requested but not detected; skipping" >&2
            fi
        done
    fi

    if [[ "${#requested[@]}" -eq 0 ]]; then
        echo "sciagent provision: no target harnesses detected; nothing to do"
        return 0
    fi

    local prefix=""
    (( dry_run )) && prefix="[dry-run] "
    echo "${prefix}provisioning harnesses: ${requested[*]}"

    local rc=0
    for h in "${requested[@]}"; do
        (( do_settings )) && { _provision_settings "$h" "$dry_run" || rc=1; }
        (( do_context ))  && { _provision_context  "$h" "$dry_run" || rc=1; }
    done
    return $rc
}

# _provision_settings <harness> <dry_run>
# claude: seed user settings.json + statusline. Others: Tier 2 — skip honestly.
_provision_settings() {
    local h="$1" dry="$2"
    case "$h" in
        claude)
            local cfg="${CLAUDE_CONFIG_DIR:-$HOME/.claude}"
            if (( dry )); then
                echo "[dry-run] claude settings: would ensure $cfg/settings.json + $cfg/statusline.sh"
                return 0
            fi
            echo "claude settings:"
            claude_settings_ensure_user_defaults
            claude_settings_ensure_user_statusline
            ;;
        *)
            echo "$h settings: not yet implemented for $h (Tier 2); skipping"
            ;;
    esac
}

# _provision_context <harness> <dry_run>
# Stamp the shared global context template into the harness's global context
# path as a SCIAGENT:CONTEXT managed block (block.sh) — idempotent, and any
# user-authored bytes in that file are preserved outside the markers.
_provision_context() {
    local h="$1" dry="$2"
    local dst src
    dst=$(harness_global_context_path "$h") || return 1
    src="$SCIAGENT_TOOLKIT/templates/global/AGENTS.md.template"

    if [[ ! -f "$src" ]]; then
        echo "sciagent provision: context template missing ($src); skipping $h context" >&2
        return 0
    fi

    if (( dry )); then
        echo "[dry-run] $h context: would stamp SCIAGENT:CONTEXT block into $dst"
        return 0
    fi

    mkdir -p "$(dirname "$dst")" || {
        echo "sciagent provision: cannot create $(dirname "$dst")" >&2
        return 1
    }

    local before="" after body
    [[ -f "$dst" ]] && before=$(block_stored_hash "$dst" CONTEXT 2>/dev/null)
    body=$(cat "$src")
    block_write "$dst" "$body" CONTEXT || {
        echo "sciagent provision: failed to write $dst" >&2
        return 1
    }
    after=$(block_stored_hash "$dst" CONTEXT 2>/dev/null)

    if [[ -n "$before" && "$before" == "$after" ]]; then
        echo "$h context: up to date ($dst)"
    else
        echo "$h context: wrote SCIAGENT:CONTEXT block -> $dst"
    fi
}
