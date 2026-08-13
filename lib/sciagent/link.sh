# lib/sciagent/link.sh — bind the toolkit catalog into one project.

# shellcheck shell=bash

_link_usage() {
    cat <<'EOF'
Usage: sciagent link [--project-dir D]

Link the toolkit's skills, agents, and commands directories into both
D/.claude/ and D/.agents/ (default: cwd). Materialize the two project
guardrail hooks and merge their registrations into .claude/settings.json.

Existing toolkit-owned mounts are swept first. A populated real category
directory is preserved and refused with relocation instructions.
EOF
}

_link_resolve() {
    local path="$1"
    if command -v realpath >/dev/null 2>&1; then
        realpath -m "$path" 2>/dev/null || printf '%s\n' "$path"
    elif [[ -d "$path" ]]; then
        (cd "$path" 2>/dev/null && pwd -P) || printf '%s\n' "$path"
    elif [[ -e "$path" ]]; then
        local dir base
        dir=$(dirname "$path")
        base=$(basename "$path")
        printf '%s/%s\n' "$(cd "$dir" 2>/dev/null && pwd -P)" "$base"
    else
        printf '%s\n' "$path"
    fi
}

_link_target_path() {
    local path="$1" target
    target=$(readlink "$path") || return 1
    if [[ "$target" == /* ]]; then
        _link_resolve "$target"
    else
        _link_resolve "$(dirname "$path")/$target"
    fi
}

# Print the project's pinned toolkit path when one is present.
_link_in_repo_toolkit() {
    local path
    if [[ -f .gitmodules ]] && command -v git >/dev/null 2>&1; then
        path=$(git config -f .gitmodules --get-regexp '^submodule\..*\.path$' 2>/dev/null \
            | awk '{print $2}' | grep -E '(^|/)SciAgent-toolkit$' | head -1)
        if [[ -n "$path" && -d "$path" ]]; then
            printf '%s\n' "$path"
            return 0
        fi
    fi
    if [[ -d 01_modules/SciAgent-toolkit ]]; then
        printf '%s\n' 01_modules/SciAgent-toolkit
        return 0
    fi
    for path in ./*/SciAgent-toolkit; do
        if [[ -d "$path" ]]; then
            printf '%s\n' "$path"
            return 0
        fi
    done
    return 1
}

_link_check_toolkit_locality() {
    local in_repo
    in_repo=$(_link_in_repo_toolkit) || return 0
    [[ "$(_link_resolve "$in_repo")" == "$(_link_resolve "$SCIAGENT_TOOLKIT")" ]] && return 0

    echo "sciagent: refusing to link against an external toolkit" >&2
    echo "  active toolkit : $(_link_resolve "$SCIAGENT_TOOLKIT")" >&2
    echo "  in-repo toolkit: $(_link_resolve "$in_repo")" >&2
    echo "  Run $in_repo/bin/sciagent link so this project stays bound to its pinned copy." >&2
    return 1
}

# True for a link into this toolkit, including legacy dangling container paths.
_link_owned_by_toolkit() {
    local path="$1" resolved root
    [[ -L "$path" ]] || return 1
    resolved=$(_link_target_path "$path") || return 1
    root=$(_link_resolve "$SCIAGENT_TOOLKIT")
    case "$resolved" in
        "$root"|"$root"/*) return 0 ;;
    esac

    # An existing outside target is user-owned. The suffix test is only for
    # dangling legacy mounts whose /workspaces project path is absent here.
    [[ -e "$path" ]] && return 1
    case "$resolved" in
        */SciAgent-toolkit/skills/*|*/SciAgent-toolkit/agents/*|*/SciAgent-toolkit/commands/*|*/SciAgent-toolkit/output-styles/*|*/SciAgent-toolkit/lib/figure-style|*/SciAgent-toolkit/lib/interactive-style)
            return 0 ;;
        *)  return 1 ;;
    esac
}

_link_sweep_dir() {
    local dir="$1" entry
    [[ -d "$dir" && ! -L "$dir" ]] || return 0
    while IFS= read -r -d '' entry; do
        if _link_owned_by_toolkit "$entry"; then
            rm "$entry" || { echo "sciagent link: failed to remove legacy link: $entry" >&2; return 1; }
            echo "removed stale toolkit link: $entry"
        fi
    done < <(find "$dir" -mindepth 1 -maxdepth 1 -type l -print0 2>/dev/null)
}

_link_refuse_directory() {
    local dst="$1"
    echo "sciagent link: refusing $dst; it is a real populated directory containing:" >&2
    find "$dst" -mindepth 1 -maxdepth 1 -printf '  - %f\n' 2>/dev/null | LC_ALL=C sort >&2
    echo "  Move these entries outside $dst, remove the empty directory, and run sciagent link again." >&2
    return 1
}

_link_category() {
    local harness="$1" category="$2"
    local dst="$harness/$category" src="$SCIAGENT_TOOLKIT/$category"

    [[ -d "$src" ]] || { echo "sciagent link: missing toolkit category: $src" >&2; return 1; }
    _link_sweep_dir "$dst" || return 1
    mkdir -p "$harness" || { echo "sciagent link: failed to create $harness" >&2; return 1; }

    if [[ -L "$dst" ]]; then
        if [[ "$(_link_target_path "$dst")" == "$(_link_resolve "$src")" ]]; then
            return 0
        fi
        local old
        old=$(readlink "$dst")
        ln -sfn "$src" "$dst" || { echo "sciagent link: failed to replace $dst" >&2; return 1; }
        echo "replaced link: $dst -> $src (was $old)"
        return 0
    fi

    if [[ -d "$dst" ]]; then
        if find "$dst" -mindepth 1 -maxdepth 1 -print -quit | grep -q .; then
            _link_refuse_directory "$dst"
            return 1
        fi
        rmdir "$dst" || { echo "sciagent link: failed to remove empty directory: $dst" >&2; return 1; }
    elif [[ -e "$dst" ]]; then
        echo "sciagent link: refusing $dst; move the existing non-symlink path aside and run sciagent link again." >&2
        return 1
    fi

    ln -s "$src" "$dst" || { echo "sciagent link: failed to create $dst" >&2; return 1; }
    echo "linked: $dst -> $src"
}

cmd_link() {
    local project_dir="" failed=0
    while [[ $# -gt 0 ]]; do
        case "$1" in
            --project-dir)
                project_dir="${2:-}"
                [[ -n "$project_dir" ]] || { echo "sciagent link: --project-dir needs a path" >&2; return 1; }
                shift 2 ;;
            --project-dir=*) project_dir="${1#*=}"; shift ;;
            -h|--help) _link_usage; return 0 ;;
            *) echo "sciagent link: unknown argument '$1'" >&2; _link_usage >&2; return 1 ;;
        esac
    done

    project_dir="${project_dir:-$(pwd)}"
    [[ -d "$project_dir" ]] || { echo "sciagent link: not a directory: $project_dir" >&2; return 1; }

    (
        cd "$project_dir" || return 1
        _link_check_toolkit_locality || return 1

        local harness category
        for harness in .claude .agents; do
            for category in skills agents commands; do
                _link_category "$harness" "$category" || failed=1
            done
        done

        _link_sweep_dir .claude/output-styles || failed=1
        # Best effort: a remaining entry means the retired directory is user-owned.
        if [[ -d .claude/output-styles ]]; then
            rmdir .claude/output-styles 2>/dev/null || true
        fi
        _link_sweep_dir 02_analysis/helpers || failed=1

        if [[ -f .sciagent/manifest.json ]]; then
            rm .sciagent/manifest.json || { echo "sciagent link: failed to remove legacy state" >&2; return 1; }
            echo "removed legacy state: .sciagent/manifest.json"
            rmdir .sciagent 2>/dev/null || true
        fi

        claude_settings_ensure_hooks || return 1
        return "$failed"
    )
}
