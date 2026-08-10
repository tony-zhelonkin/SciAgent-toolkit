# lib/sciagent/craft_verb.sh — sciagent craft [--project-dir D] [--force] [--quiet]
#
# Renders/refreshes the SCIAGENT:CRAFT block in D/AGENTS.md and does nothing
# else: no symlinks, no manifest, no settings, no role. CRAFT used to arrive
# only via `activate`, which is why coordination repos never got it — activating
# there would have mounted the whole bioinformatics catalog. This verb is the
# split: craft conventions travel on their own.
#
# Drift guard: a hand-edit inside the markers leaves the body disagreeing with
# the hash in the BEGIN marker. Refuse to overwrite that silently — the edit is
# someone's intent, and the toolkit is not the only author of the file it is a
# guest in. --force renders anyway, discarding the edit.
#
# Depends on block.sh (framing/hashing) and craft.sh (the renderer); both are
# reused unchanged.

# shellcheck shell=bash

_craft_usage() {
    cat <<'EOF'
Usage: sciagent craft [--project-dir D] [--force] [--quiet]

Render or refresh the SCIAGENT:CRAFT block in D/AGENTS.md (default: cwd).
Mounts nothing. Re-running against an up-to-date block is a no-op.

  --project-dir D   Project root holding AGENTS.md (default: current directory)
  --force           Overwrite a drifted block (hand-edited inside the markers)
  --quiet           Print nothing on success
EOF
}

cmd_craft() {
    local projdir="" force=0 quiet=0
    while [[ $# -gt 0 ]]; do
        case "$1" in
            --project-dir)
                projdir="${2:-}"
                [[ -n "$projdir" ]] || { echo "sciagent craft: --project-dir needs a path" >&2; return 1; }
                shift 2 ;;
            --project-dir=*) projdir="${1#*=}"; shift ;;
            --force)  force=1; shift ;;
            --quiet)  quiet=1; shift ;;
            -h|--help) _craft_usage; return 0 ;;
            *) echo "sciagent craft: unknown argument '$1'" >&2; _craft_usage >&2; return 1 ;;
        esac
    done

    projdir="${projdir:-$(pwd)}"
    if [[ ! -d "$projdir" ]]; then
        echo "sciagent craft: not a directory: $projdir" >&2
        return 1
    fi

    # The library renderer treats a missing craft.yaml as "nothing to do"; an
    # explicit verb must say so instead of exiting 0 having written nothing.
    local yaml
    yaml=$(_craft_yaml_path)
    if [[ ! -f "$yaml" ]]; then
        echo "sciagent craft: no craft.yaml in this toolkit ($yaml)" >&2
        return 1
    fi

    local target="$projdir/AGENTS.md"

    block_hash_check "$target" "$CRAFT_BLOCK_ID"
    local state=$?
    case "$state" in
        2)
            echo "sciagent craft: $target has one CRAFT marker but not the other" >&2
            echo "  repair the markers by hand, then re-run." >&2
            return 1 ;;
        3)
            if (( force == 0 )); then
                echo "sciagent craft: CRAFT block in $target has drifted (hand-edited)" >&2
                echo "  the block is toolkit-managed; edit craft.yaml in the toolkit instead." >&2
                echo "  re-run with --force to discard the local edit." >&2
                return 1
            fi
            (( quiet )) || echo "sciagent craft: discarding drifted block (--force)" ;;
    esac

    local before after
    before=$(block_stored_hash "$target" "$CRAFT_BLOCK_ID" 2>/dev/null)
    craft_render_and_write "$target" || return 1
    after=$(block_stored_hash "$target" "$CRAFT_BLOCK_ID" 2>/dev/null)

    (( quiet )) && return 0
    if [[ -z "$before" ]]; then
        echo "added SCIAGENT:CRAFT block to: $target"
    elif [[ "$before" == "$after" ]]; then
        echo "SCIAGENT:CRAFT block already current: $target"
    else
        echo "updated SCIAGENT:CRAFT block in: $target"
    fi
}
