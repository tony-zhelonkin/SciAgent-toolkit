# lib/sciagent/new.sh — sciagent new <kind> [args]
# Scaffolding per architecture §11.

# shellcheck shell=bash

cmd_new() {
    local kind="${1:-}"
    shift || true
    case "$kind" in
        project) _new_project "$@" ;;
        role)    _new_role    "$@" ;;
        skill)   _new_skill   "$@" ;;
        agent)   _new_agent   "$@" ;;
        ""|-h|--help)
            cat <<EOF
usage:
  sciagent new project <dir> [--type analysis|software-tool] [--species ...]
                             [--genome ...] [--title ...] [--git]
                             [--with-submodules] [--force]
      --type analysis        00_data/ 01_modules/ 02_analysis/ 03_results/ tree (default)
      --type software-tool   src/ tests/ docs/ examples/ packageable-library tree
  sciagent new role <name>       scaffold roles/<name>.yaml
  sciagent new skill <name>      copy skills/_TEMPLATE/ to skills/<name>/
  sciagent new agent <name>      scaffold agents/<name>.md
EOF
            return 0 ;;
        *)
            echo "sciagent new: unknown kind '$kind'" >&2
            return 1 ;;
    esac
}

# _derive_species_db <species> → echoes the short DB code (MM / HS / "").
_derive_species_db() {
    local species
    species=$(printf '%s' "${1:-}" | tr '[:upper:]' '[:lower:]')
    case "$species" in
        mouse|*mus*musculus*) echo "MM" ;;
        human|*homo*sapiens*) echo "HS" ;;
        *)                    echo "" ;;
    esac
}

# _subst <src> <dst> — renders one template, substituting the {{TOKENS}} held in the
# PROJECT_VARS associative array (populated by the caller). Errors out (returns) if the
# write fails so the caller can surface it.
_subst() {
    local src="$1" dst="$2"
    local -a sed_args=()
    local key
    for key in "${!PROJECT_VARS[@]}"; do
        sed_args+=(-e "s|{{$key}}|${PROJECT_VARS[$key]}|g")
    done
    sed "${sed_args[@]}" "$src" > "$dst" || {
        echo "render failed: $dst" >&2
        return 1
    }
}

# _mkdirs <base> <dir>... — mkdir -p each child of <base>, error-checked.
_mkdirs() {
    local base="$1"; shift
    local d
    for d in "$@"; do
        mkdir -p "$base/$d" || { echo "mkdir failed: $base/$d" >&2; return 1; }
    done
}

# _seed_gitkeep <root> — drop a .gitkeep into every empty directory under <root>.
_seed_gitkeep() {
    local root="$1" d
    while IFS= read -r d; do
        touch "$d/.gitkeep" || { echo "gitkeep failed: $d" >&2; return 1; }
    done < <(find "$root" -type d -empty)
}

# _render_tree <tpl_root> <dir> <force> — walk *.template under <tpl_root>, strip the
# suffix, re-root under <dir>, render via _subst. Non-template files (e.g. .gitkeep,
# .gitignore-seed) are handled by the caller, not here.
_render_tree() {
    local tpl_root="$1" dir="$2" force="$3"
    [[ -d "$tpl_root" ]] || return 0
    local src rel out
    while IFS= read -r src; do
        rel="${src#"$tpl_root"/}"
        out="$dir/${rel%.template}"
        if [[ -e "$out" && "$force" != "true" ]]; then
            echo "skip (exists): $out"
            continue
        fi
        mkdir -p "$(dirname "$out")" || return 1
        _subst "$src" "$out" || return 1
        echo "wrote: $out"
    done < <(find "$tpl_root" -type f -name '*.template')
}

# Directory lists per project type (relative to project root).
_dirs_for_type() {
    case "$1" in
        analysis)
            echo "00_data/raw 00_data/processed 00_data/references" \
                 "01_modules" \
                 "02_analysis/config 02_analysis/helpers 02_analysis/scripts 02_analysis/notebooks" \
                 "03_results/objects 03_results/master 03_results/interactive 03_results/_scratch" \
                 "03_results/01_qc/tables 03_results/01_qc/figures" \
                 "03_results/02_eda/tables 03_results/02_eda/figures" \
                 "docs/plan docs/_internal/reasoning docs/_internal/sessions docs/_internal/research" \
                 "logs"
            ;;
        software-tool)
            echo "src tests docs examples" \
                 "docs/_internal/reasoning docs/_internal/design docs/_internal/benchmarks"
            ;;
    esac
}

_new_project() {
    local dir="" type="analysis" species="Mus musculus" genome="mm10" title=""
    local do_git=false with_submodules=false force=false
    while (( $# )); do
        case "$1" in
            --type)             type="$2"; shift 2 ;;
            --species)          species="$2"; shift 2 ;;
            --genome)           genome="$2"; shift 2 ;;
            --title)            title="$2"; shift 2 ;;
            --git)              do_git=true; shift ;;
            --with-submodules)  with_submodules=true; shift ;;
            --force)            force=true; shift ;;
            -*)  echo "sciagent new project: unknown flag '$1'" >&2; return 1 ;;
            *)   if [[ -z "$dir" ]]; then dir="$1"; else
                     echo "sciagent new project: unexpected argument '$1'" >&2; return 1
                 fi; shift ;;
        esac
    done
    dir="${dir:-.}"

    case "$type" in
        analysis|software-tool) ;;
        *) echo "sciagent new project: --type must be 'analysis' or 'software-tool' (got '$type')" >&2
           return 1 ;;
    esac

    mkdir -p "$dir" || { echo "mkdir failed: $dir" >&2; return 1; }
    local pid abs
    abs=$(cd "$dir" && pwd) || return 1
    pid=$(basename "$abs")

    # Build the substitution vocabulary (ADR-2.2).
    local species_db
    species_db=$(_derive_species_db "$species")
    declare -gA PROJECT_VARS=(
        [PROJECT_ID]="$pid"
        [PROJECT_NAME]="$pid"                         # deprecated alias for PROJECT_ID
        [PROJECT_TITLE]="${title:-$pid Analysis}"
        [DATE]="$(date +%Y-%m-%d)"
        [PROJECT_TYPE]="$type"
        [SPECIES]="$species"
        [SPECIES_DB]="$species_db"
        [GENOME_BUILD]="$genome"
    )

    # 1. directory tree for this type.
    local -a dirlist
    read -r -a dirlist <<< "$(_dirs_for_type "$type")"
    _mkdirs "$dir" "${dirlist[@]}" || return 1

    # 2. render templates: _common/ first, then the type overlay.
    local proj_tpl="$SCIAGENT_TOOLKIT/templates/project"
    _render_tree "$proj_tpl/_common" "$dir" "$force" || return 1
    _render_tree "$proj_tpl/$type"   "$dir" "$force" || return 1

    # 3. seed the shared .gitignore (non-.template source) if absent.
    local gi_seed="$proj_tpl/_common/.gitignore-seed"
    if [[ -f "$gi_seed" && ( ! -e "$dir/.gitignore" || "$force" == "true" ) ]]; then
        cp "$gi_seed" "$dir/.gitignore" && echo "wrote: $dir/.gitignore" || return 1
    fi

    # 4. seed .gitkeep into empty dirs.
    _seed_gitkeep "$dir" || return 1

    # 5. post-render validation: warn on any unresolved {{TOKEN}}.
    _warn_unresolved_tokens "$dir"

    # 6. optional git init + submodule registration.
    if [[ "$do_git" == "true" || "$with_submodules" == "true" ]]; then
        _new_git "$abs" "$type" "$with_submodules" || return 1
    fi

    _new_project_next_steps "$dir" "$abs" "$type" "$with_submodules"
}

# _warn_unresolved_tokens <dir> — surface any {{TOKEN}} left in rendered output. Some
# templates intentionally carry no tokens; this only warns, never fails.
_warn_unresolved_tokens() {
    local dir="$1" hits
    hits=$(grep -rl '{{' "$dir" 2>/dev/null) || return 0
    [[ -n "$hits" ]] || return 0
    echo "warning: unresolved {{tokens}} in:" >&2
    printf '  %s\n' $hits >&2
}

# _default_role_for_type <type> — the role the "Next:" hint suggests activating.
#   analysis → base ; software-tool → software-tool (umbrella uses base, same as analysis).
_default_role_for_type() {
    case "$1" in
        software-tool) echo "software-tool" ;;
        *)             echo "base" ;;
    esac
}

# _note_child_software_tool <abs> — when a software-tool is scaffolded inside an existing
# project's 01_modules/, it gets its own independent activation root (ADR-6.3); flag that.
# Detection: the parent directory of <abs> is named "01_modules". Emits nothing otherwise.
_note_child_software_tool() {
    local abs="$1" parent grandparent
    parent=$(dirname "$abs")
    [[ "$(basename "$parent")" == "01_modules" ]] || return 0
    grandparent=$(dirname "$parent")
    echo
    echo "Note: scaffolded a child software-tool under $grandparent."
    echo "      It has its own activation root (independent of the parent, depth-2 cap intact)."
    echo "      Activate its context: cd $abs && sciagent activate software-tool"
}

# _new_project_next_steps <dir> <abs> <type> <with_submodules>
_new_project_next_steps() {
    local dir="$1" abs="$2" type="$3" with_submodules="$4"
    local role
    role=$(_default_role_for_type "$type")
    echo
    echo "Project scaffolded at $dir (type: $type)."
    if [[ "$with_submodules" != "true" ]]; then
        echo "Hint: re-run with --with-submodules to attach toolkits under 01_modules/."
    fi

    # Create the canonical docs/ tree.
    mkdir -p "$dir/docs/stages" "$dir/docs/reference"
    mkdir -p "$dir/docs/_internal/reasoning" "$dir/docs/_internal/research"
    mkdir -p "$dir/docs/_internal/plans" "$dir/docs/_internal/reports" "$dir/docs/_internal/handoffs"

    touch "$dir/docs/_internal/reasoning/.gitkeep"
    touch "$dir/docs/_internal/research/.gitkeep"
    touch "$dir/docs/_internal/plans/.gitkeep"
    touch "$dir/docs/_internal/reports/.gitkeep"
    touch "$dir/docs/_internal/handoffs/.gitkeep"

    if [[ ! -e "$dir/docs/README.md" ]]; then
        _subst <(cat <<'TMPL'
# docs — {{PROJECT_ID}}

Navigation hub for project documentation.

| Directory | Contents |
|-----------|----------|
| `stages/` | Durable per-phase methodology guides (committed) |
| `reference/` | External dataset READMEs and data lineage (committed) |
| `_internal/` | Gitignored: reasoning traces, research notes, plans, reports, handoffs |

Generated by `sciagent new project` on {{DATE}}.
TMPL
        ) "$dir/docs/README.md" "$pid"
        echo "wrote: $dir/docs/README.md"
    fi

    # Apply the managed gitignore block.
    . "$SCIAGENT_TOOLKIT/lib/sciagent/gitignore.sh"
    cmd_gitignore "${dir}/.gitignore"

    echo
    echo "Next steps:"
    echo "  1. cd $dir"
    echo "  2. Add container substrate (from scbio-docker):"
    echo "       <scbio-docker>/scripts/init-container.sh $dir --type $type"
    echo "  3. Open in VS Code → Reopen in Container"
    echo "  4. sciagent activate $role"

    if [[ "$type" == "software-tool" ]]; then
        _note_child_software_tool "$abs"
    fi
}

# _new_git <abs_dir> <type> <with_submodules> — git init (if absent) and, when requested,
# register the type's toolkit submodules under 01_modules/. SSH → gh-HTTPS fallback ported
# from init-project.sh.
_new_git() {
    local abs_dir="$1" type="$2" with_submodules="$3"
    local owner="tony-zhelonkin"

    if [[ ! -d "$abs_dir/.git" ]]; then
        git -C "$abs_dir" init -q || { echo "git init failed" >&2; return 1; }
        echo "git: initialized repository in $abs_dir"
    fi

    [[ "$with_submodules" == "true" ]] || return 0

    # Per-type toolkit set: "<repo> <branch> <target>".
    local -a subs=()
    case "$type" in
        analysis)
            subs+=("RNAseq-toolkit dev 01_modules/RNAseq-toolkit")
            subs+=("SciAgent-toolkit main 01_modules/SciAgent-toolkit")
            ;;
        software-tool)
            subs+=("SciAgent-toolkit main 01_modules/SciAgent-toolkit")
            ;;
    esac

    local spec repo branch target added=0
    for spec in "${subs[@]}"; do
        read -r repo branch target <<< "$spec"
        echo "git: adding $repo (branch: $branch)..."
        if _add_submodule_with_fallback "$abs_dir" "$owner" "$repo" "$branch" "$target"; then
            added=$((added + 1))
        fi
    done

    if (( added > 0 )); then
        git -C "$abs_dir" add .gitmodules >/dev/null 2>&1 || true
        git -C "$abs_dir" commit -q -m "add analysis toolkits as git submodules" \
            >/dev/null 2>&1 || true
    fi
}

# _add_submodule_with_fallback <abs_dir> <owner> <repo> <branch> <target>
# Tries SSH `git submodule add`, then falls back to a gh-CLI HTTPS clone that rewrites the
# recorded URL back to SSH. Returns non-zero (without aborting the scaffold) on failure.
_add_submodule_with_fallback() {
    local abs_dir="$1" owner="$2" repo="$3" branch="$4" target="$5"
    local ssh_url="git@github.com:${owner}/${repo}.git"

    if [[ -d "$abs_dir/$target" ]]; then
        echo "  $repo directory already exists, skipping"
        return 1
    fi

    if git -C "$abs_dir" submodule add -b "$branch" "$ssh_url" "$target" 2>/dev/null; then
        echo "  $repo added via git (SSH)"
        return 0
    fi

    echo "  SSH failed, trying gh CLI fallback (HTTPS)..."
    if command -v gh >/dev/null 2>&1 && gh auth status >/dev/null 2>&1; then
        if GIT_CONFIG_COUNT=1 \
           GIT_CONFIG_KEY_0="url.https://github.com/.insteadOf" \
           GIT_CONFIG_VALUE_0="git@github.com:" \
           gh repo clone "${owner}/${repo}" "$abs_dir/$target" -- -b "$branch" 2>/dev/null; then
            git -C "$abs_dir" config -f .gitmodules "submodule.${target}.path" "$target"
            git -C "$abs_dir" config -f .gitmodules "submodule.${target}.url" "$ssh_url"
            git -C "$abs_dir" config -f .gitmodules "submodule.${target}.branch" "$branch"
            git -C "$abs_dir" config "submodule.${target}.url" "$ssh_url"
            git -C "$abs_dir" config "submodule.${target}.active" "true"
            git -C "$abs_dir" add "$target" 2>/dev/null || true
            echo "  $repo added via gh CLI (HTTPS)"
            return 0
        fi
    fi

    echo "  failed to add $repo" >&2
    return 1
}

_new_role() {
    local name="${1:-}"
    if [[ -z "$name" ]]; then
        echo "usage: sciagent new role <name>" >&2
        return 1
    fi
    local out="$SCIAGENT_TOOLKIT/roles/$name.yaml"
    if [[ -e "$out" ]]; then
        echo "refuse to overwrite: $out" >&2
        return 1
    fi
    cat > "$out" <<EOF
name: $name
description: TODO — one-line description of $name role

skills: []
agents: []
commands: []
# output_style: architect-mentor   # optional, Claude-specific
EOF
    echo "wrote: $out"
}

_new_skill() {
    local name="${1:-}"
    if [[ -z "$name" ]]; then
        echo "usage: sciagent new skill <name>" >&2
        return 1
    fi
    local tpl="$SCIAGENT_TOOLKIT/skills/_TEMPLATE"
    if [[ ! -d "$tpl" ]]; then
        echo "skill template not found: $tpl" >&2
        return 1
    fi
    local dst="$SCIAGENT_TOOLKIT/skills/$name"
    if [[ -e "$dst" ]]; then
        echo "refuse to overwrite: $dst" >&2
        return 1
    fi
    cp -R "$tpl" "$dst"
    # Substitute {{SKILL_NAME}} in any text files (best-effort).
    local f
    while IFS= read -r f; do
        sed -i "s|{{SKILL_NAME}}|$name|g" "$f" 2>/dev/null || true
    done < <(find "$dst" -type f \( -name '*.md' -o -name '*.txt' -o -name '*.yaml' -o -name '*.yml' \))
    echo "wrote: $dst/"
}

_new_agent() {
    local name="${1:-}"
    if [[ -z "$name" ]]; then
        echo "usage: sciagent new agent <name>" >&2
        return 1
    fi
    local out="$SCIAGENT_TOOLKIT/agents/$name.md"
    if [[ -e "$out" ]]; then
        echo "refuse to overwrite: $out" >&2
        return 1
    fi
    cat > "$out" <<EOF
---
name: $name
description: TODO — when to invoke this agent (with example triggers)
model: sonnet
color: blue
---

# $name

TODO — describe the agent's role, methodology, and outputs.
EOF
    echo "wrote: $out"
}
