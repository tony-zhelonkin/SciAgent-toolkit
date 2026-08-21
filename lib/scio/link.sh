# lib/scio/link.sh — bind the toolkit catalog into one project.

# shellcheck shell=bash

_link_usage() {
    cat <<'EOF'
Usage: scio link [--project-dir D]

Link the toolkit's skills, agents, and commands directories into both
D/.claude/ and D/.agents/ (default: cwd). Materialize the two project
guardrail hooks, merge their registrations into .claude/settings.json, and
refresh the SCIO:GITIGNORE block in .gitignore.

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

# _link_symlink_target <link_dir> <src>
# Echo the target to write into a symlink placed in <link_dir>.
#
# Relative whenever the source lives inside this project, so a single `link` run
# resolves from the host and from inside a container: the same project tree is
# mounted at two different absolute paths, and an absolute target can only ever
# name one of them. Absolute when the toolkit is installed outside the project,
# where a relative chain would break as soon as the project directory moves.
_link_symlink_target() {
    local link_dir="$1" src="$2"
    local proj src_abs rel
    proj=$(pwd -P 2>/dev/null) || { printf '%s' "$src"; return 0; }
    src_abs=$(_link_resolve "$src")
    case "$src_abs/" in
        "$proj"/*) ;;
        *) printf '%s' "$src_abs"; return 0 ;;
    esac
    rel=$(realpath --relative-to="$link_dir" "$src_abs" 2>/dev/null) \
        || { printf '%s' "$src_abs"; return 0; }
    printf '%s' "$rel"
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
    [[ "$(_link_resolve "$in_repo")" == "$(_link_resolve "$SCIO_TOOLKIT")" ]] && return 0

    echo "scio: refusing to link against an external toolkit" >&2
    echo "  active toolkit : $(_link_resolve "$SCIO_TOOLKIT")" >&2
    echo "  in-repo toolkit: $(_link_resolve "$in_repo")" >&2
    echo "  Run $in_repo/bin/scio link so this project stays bound to its pinned copy." >&2
    return 1
}

# Carry ownership evidence forward before any materialized body is inspected.
# Remove this migration once all 9 consumer checkouts have re-pinned past the
# rebrand.
_link_migrate_state_dir() {
    [[ -e .sciagent || -L .sciagent ]] || return 0

    if [[ ! -d .sciagent || -L .sciagent ]]; then
        echo "scio link: refusing to migrate .sciagent; it is not a real directory" >&2
        return 1
    fi
    if [[ -e .scio || -L .scio ]]; then
        echo "scio link: refusing to migrate .sciagent because .scio already exists" >&2
        echo "  Reconcile the two state directories, then run scio link again." >&2
        return 1
    fi

    mv .sciagent .scio || { echo "scio link: failed to migrate .sciagent to .scio" >&2; return 1; }
    echo "migrated state: .sciagent -> .scio"
}

# True for a link into this toolkit, including legacy dangling container paths.
_link_owned_by_toolkit() {
    local path="$1" resolved root
    [[ -L "$path" ]] || return 1
    resolved=$(_link_target_path "$path") || return 1
    root=$(_link_resolve "$SCIO_TOOLKIT")
    case "$resolved" in
        "$root"|"$root"/*) return 0 ;;
    esac

    # An existing outside target is user-owned. The suffix test is only for
    # dangling legacy mounts whose /workspaces project path is absent here.
    [[ -e "$path" ]] && return 1
    case "$resolved" in
        */SciAgent-toolkit/skills/*|*/SciAgent-toolkit/agents/*|*/SciAgent-toolkit/commands/*|*/SciAgent-toolkit/system-prompts/*|*/SciAgent-toolkit/lib/figure-style|*/SciAgent-toolkit/lib/interactive-style)
            return 0 ;;
        *)  return 1 ;;
    esac
}

_link_sweep_dir() {
    local dir="$1" entry
    [[ -d "$dir" && ! -L "$dir" ]] || return 0
    while IFS= read -r -d '' entry; do
        if _link_owned_by_toolkit "$entry"; then
            rm "$entry" || { echo "scio link: failed to remove legacy link: $entry" >&2; return 1; }
            echo "removed stale toolkit link: $entry"
        fi
    done < <(find "$dir" -mindepth 1 -maxdepth 1 -type l -print0 2>/dev/null)
}

_link_refuse_directory() {
    local dst="$1"
    echo "scio link: refusing $dst; it is a real populated directory containing:" >&2
    find "$dst" -mindepth 1 -maxdepth 1 -printf '  - %f\n' 2>/dev/null | LC_ALL=C sort >&2
    echo "  Move these entries outside $dst, remove the empty directory, and run scio link again." >&2
    return 1
}

_link_category() {
    local harness="$1" category="$2"
    local dst="$harness/$category" src="$SCIO_TOOLKIT/$category"

    [[ -d "$src" ]] || { echo "scio link: missing toolkit category: $src" >&2; return 1; }
    _link_sweep_dir "$dst" || return 1
    mkdir -p "$harness" || { echo "scio link: failed to create $harness" >&2; return 1; }

    local target
    target=$(_link_symlink_target "$harness" "$src")

    # Compare the written target, not where it resolves to. An absolute link
    # resolves correctly on the side that wrote it, so a resolution test reads
    # as current and leaves every already-bound project unrepaired.
    if [[ -L "$dst" ]]; then
        if [[ "$(readlink "$dst")" == "$target" ]]; then
            return 0
        fi
        local old
        old=$(readlink "$dst")
        ln -sfn "$target" "$dst" || { echo "scio link: failed to replace $dst" >&2; return 1; }
        echo "replaced link: $dst -> $target (was $old)"
        return 0
    fi

    if [[ -d "$dst" ]]; then
        if find "$dst" -mindepth 1 -maxdepth 1 -print -quit | grep -q .; then
            _link_refuse_directory "$dst"
            return 1
        fi
        rmdir "$dst" || { echo "scio link: failed to remove empty directory: $dst" >&2; return 1; }
    elif [[ -e "$dst" ]]; then
        echo "scio link: refusing $dst; move the existing non-symlink path aside and run scio link again." >&2
        return 1
    fi

    ln -s "$target" "$dst" || { echo "scio link: failed to create $dst" >&2; return 1; }
    echo "linked: $dst -> $target"
}

_HELPER_SHIM_TPL_REL="project/analysis/02_analysis/helpers"
_HELPER_SHIM_STATE_DIR=".scio/helper_shim_state"

# Sweep retired helper mounts while preserving canonical live bindings.
_link_sweep_helper_links() {
    local dir="02_analysis/helpers" entry name src
    [[ -d "$dir" && ! -L "$dir" ]] || return 0

    while IFS= read -r -d '' entry; do
        name=$(basename "$entry")
        case "$name" in
            figure-style|interactive-style)
                src="$SCIO_TOOLKIT/lib/$name"
                if [[ -d "$src" && "$(_link_target_path "$entry")" == "$(_link_resolve "$src")" ]]; then
                    continue
                fi ;;
        esac
        if _link_owned_by_toolkit "$entry"; then
            rm "$entry" || { echo "scio link: failed to remove legacy link: $entry" >&2; return 1; }
            echo "removed stale toolkit link: $entry"
        fi
    done < <(find "$dir" -mindepth 1 -maxdepth 1 -type l -print0 2>/dev/null)
}

# Bind the shared libraries beside their importable project shims.
_link_helper_libs() {
    local link_dir="02_analysis/helpers" libdir src dst target old failed=0
    mkdir -p "$link_dir" || { echo "scio link: failed to create $link_dir" >&2; return 1; }

    for libdir in figure-style interactive-style; do
        src="$SCIO_TOOLKIT/lib/$libdir"
        [[ -d "$src" ]] || continue
        dst="$link_dir/$libdir"

        target=$(_link_symlink_target "$link_dir" "$src")

        if [[ -L "$dst" ]]; then
            if [[ "$(readlink "$dst")" == "$target" ]]; then
                continue
            fi
            if ! _link_owned_by_toolkit "$dst"; then
                echo "scio link: refusing $dst; its symlink target is user-owned" >&2
                failed=1
                continue
            fi
            old=$(readlink "$dst")
            if ! rm "$dst"; then
                echo "scio link: failed to remove legacy link: $dst" >&2
                failed=1
                continue
            fi
            echo "removed stale toolkit link: $dst -> $old"
        elif [[ -e "$dst" ]]; then
            echo "scio link: refusing $dst; move the existing non-symlink path aside and run scio link again." >&2
            failed=1
            continue
        fi

        ln -s "$target" "$dst" || { echo "scio link: failed to create $dst" >&2; return 1; }
        echo "linked: $dst -> $target"
    done
    return "$failed"
}

# Materialize import shims with content-provenance ownership.
_link_helper_shims() {
    local tpl_dir="$SCIO_TOOLKIT/templates/$_HELPER_SHIM_TPL_REL"
    [[ -d "$tpl_dir" ]] || return 0

    local src base
    for src in "$tpl_dir"/*.template; do
        [[ -f "$src" ]] || continue
        base=$(basename "${src%.template}")
        ownership_ensure_body \
            "$src" "02_analysis/helpers/$base" \
            "$_HELPER_SHIM_TPL_REL/$base.template" \
            "$_HELPER_SHIM_STATE_DIR/$base.sha1" \
            "$_HELPER_SHIM_STATE_DIR/$base.ceded" \
            plain || return 1
    done
}

# Materialize the analysis helper seam when the project declares that layout.
_link_analysis_helpers() {
    [[ -d 02_analysis ]] || return 0

    local failed=0
    _link_sweep_helper_links || failed=1
    _link_helper_libs || failed=1
    _link_helper_shims || failed=1
    return "$failed"
}

_LINK_GITIGNORE_BEGIN='# BEGIN SCIO:GITIGNORE'
_LINK_GITIGNORE_END='# END SCIO:GITIGNORE'
_LINK_GITIGNORE_LEGACY_BEGIN='# BEGIN SCIAGENT:GITIGNORE'
_LINK_GITIGNORE_LEGACY_END='# END SCIAGENT:GITIGNORE'

_link_gitignore_block() {
    cat <<'EOF'
# BEGIN SCIO:GITIGNORE
docs/_internal/
.claude/
.agents/
.gemini/
.scio/
02_analysis/helpers/figure-style
.mcp.json
.env
.env.*
# END SCIO:GITIGNORE
EOF
}

# Refresh the managed ignore block as part of project binding. Remove the
# SCIAGENT marker pair once all 9 consumer checkouts have re-pinned past the
# rebrand.
_link_ensure_gitignore() {
    local target=".gitignore" stripped desired target_mode="644"
    stripped=$(mktemp) || return 1
    desired=$(mktemp) || { rm -f "$stripped"; return 1; }

    if [[ -f "$target" ]]; then
        target_mode=$(stat -c '%a' "$target" 2>/dev/null || printf '644')
        awk -v begin="$_LINK_GITIGNORE_BEGIN" -v end="$_LINK_GITIGNORE_END" \
            -v legacy_begin="$_LINK_GITIGNORE_LEGACY_BEGIN" \
            -v legacy_end="$_LINK_GITIGNORE_LEGACY_END" '
            $0 == begin || $0 == legacy_begin { in_block=1; next }
            $0 == end || $0 == legacy_end { in_block=0; next }
            !in_block { lines[++n]=$0 }
            END {
                while (n > 0 && lines[n] == "") n--
                for (i=1; i<=n; i++) print lines[i]
            }
        ' "$target" > "$stripped"
    fi
    {
        cat "$stripped"
        [[ ! -s "$stripped" ]] || printf '\n'
        _link_gitignore_block
    } > "$desired"

    if [[ -f "$target" ]] && cmp -s "$desired" "$target"; then
        rm -f "$stripped" "$desired"
        return 0
    fi
    chmod "$target_mode" "$desired" || { rm -f "$stripped" "$desired"; return 1; }
    mv "$desired" "$target" || { rm -f "$stripped" "$desired"; return 1; }
    rm -f "$stripped"
    echo "refreshed: .gitignore SCIO:GITIGNORE block"
}

cmd_link() {
    local project_dir="" failed=0
    while [[ $# -gt 0 ]]; do
        case "$1" in
            --project-dir)
                project_dir="${2:-}"
                [[ -n "$project_dir" ]] || { echo "scio link: --project-dir needs a path" >&2; return 1; }
                shift 2 ;;
            --project-dir=*) project_dir="${1#*=}"; shift ;;
            -h|--help) _link_usage; return 0 ;;
            *) echo "scio link: unknown argument '$1'" >&2; _link_usage >&2; return 1 ;;
        esac
    done

    project_dir="${project_dir:-$(pwd)}"
    [[ -d "$project_dir" ]] || { echo "scio link: not a directory: $project_dir" >&2; return 1; }

    (
        cd "$project_dir" || return 1
        _link_check_toolkit_locality || return 1
        _link_migrate_state_dir || return 1

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
        _link_analysis_helpers || failed=1
        _link_ensure_gitignore || { echo "scio link: failed to refresh .gitignore" >&2; return 1; }

        if [[ -f .scio/manifest.json ]]; then
            rm .scio/manifest.json || { echo "scio link: failed to remove legacy state" >&2; return 1; }
            echo "removed legacy state: .scio/manifest.json"
            rmdir .scio 2>/dev/null || true
        fi

        claude_settings_ensure_hooks || return 1
        return "$failed"
    )
}
