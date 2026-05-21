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
  sciagent new project [<dir>]   bootstrap a project (copy templates)
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

_subst() {
    # _subst <src> <dst> <project_id>
    local src="$1" dst="$2" pid="$3"
    local date_str
    date_str=$(date +%Y-%m-%d)
    sed -e "s|{{PROJECT_ID}}|$pid|g" \
        -e "s|{{DATE}}|$date_str|g" \
        -e "s|{{SKILL_NAME}}|$pid|g" \
        "$src" > "$dst"
}

_new_project() {
    local dir="${1:-.}"
    mkdir -p "$dir"
    local pid
    pid=$(basename "$(cd "$dir" && pwd)")

    local tpl_dir="$SCIAGENT_TOOLKIT/templates"
    local f base out
    local copied=0
    for f in AGENTS.md.template CLAUDE.md.template context.md.template; do
        [[ -f "$tpl_dir/$f" ]] || continue
        base="${f%.template}"
        out="$dir/$base"
        if [[ -e "$out" ]]; then
            echo "skip (exists): $out"
            continue
        fi
        _subst "$tpl_dir/$f" "$out" "$pid"
        echo "wrote: $out"
        copied=$((copied + 1))
    done
    if (( copied == 0 )); then
        echo "no new files written (all templates already present)"
    fi
    echo
    echo "Next: cd $dir && sciagent activate <role>"
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
