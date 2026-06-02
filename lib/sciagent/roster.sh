# lib/sciagent/roster.sh — enumerate the active agent set from `.claude/agents/`.
#
# Reads YAML frontmatter from each `.claude/agents/*.md` and prints either a
# human-readable table (default) or a JSON manifest (--json). Pure bash + awk;
# no yq, no jq. The optional `domain:` and `outputs:` frontmatter fields drive
# the richer columns — absent fields render blank, never an error.
#
# This module only reads the already-activated `.claude/agents/` directory; it
# does not touch the activation mechanism.

# shellcheck shell=bash

# _roster_frontmatter <agent.md>
# Print the lines between the first and second `---` fences (the frontmatter).
_roster_frontmatter() {
    awk '
        BEGIN { in_fm = 0 }
        /^---[ \t]*$/ { in_fm++; next }
        in_fm == 1 { print; next }
        in_fm > 1 { exit }
    ' "$1"
}

# _roster_scalar <frontmatter> <key>
# Print a top-level scalar value (key: value). Empty if absent.
_roster_scalar() {
    local fm="$1" key="$2"
    printf '%s\n' "$fm" | awk -v k="$key" '
        BEGIN { pat = "^" k ":[ \t]*" }
        $0 ~ pat {
            sub(pat, "")
            sub(/[ \t]*#.*$/, "")
            sub(/^["'\'']/, ""); sub(/["'\'']$/, "")
            sub(/[ \t]+$/, "")
            print
            exit
        }
    '
}

# _roster_list <frontmatter> <key>
# Print items of a top-level YAML list (`key:` then `  - item` lines) one per
# line. Stops at the next top-level key.
_roster_list() {
    local fm="$1" key="$2"
    printf '%s\n' "$fm" | awk -v k="$key" '
        BEGIN { in_block = 0 }
        {
            line = $0
            if (line ~ /^[a-zA-Z_][a-zA-Z0-9_]*:/) {
                in_block = (line ~ ("^" k ":[ \t]*(#.*)?$")) ? 1 : 0
                next
            }
            if (in_block && line ~ /^[ \t]*-[ \t]+/) {
                sub(/^[ \t]*-[ \t]+/, "", line)
                sub(/[ \t]*#.*$/, "", line)
                sub(/[ \t]+$/, "", line)
                if (line != "") print line
            }
        }
    '
}

# _roster_outputs_field <frontmatter> <field>
# Print a value nested one level under `outputs:` (e.g. default_path, kind).
# Empty if `outputs:` or the field is absent.
_roster_outputs_field() {
    local fm="$1" field="$2"
    printf '%s\n' "$fm" | awk -v f="$field" '
        BEGIN { in_outputs = 0 }
        {
            if ($0 ~ /^outputs:[ \t]*$/) { in_outputs = 1; next }
            if ($0 ~ /^[a-zA-Z_]/)       { in_outputs = 0 }
            if (in_outputs && $0 ~ ("^[ \t]+" f ":[ \t]*")) {
                sub(("^[ \t]+" f ":[ \t]*"), "")
                sub(/[ \t]*#.*$/, "")
                sub(/[ \t]+$/, "")
                print
                exit
            }
        }
    '
}

# _roster_description <agent.md>
# Print the description as a single line. Handles both inline `description: text`
# and block-scalar `description: |` forms; for the block form, the first
# non-blank content line is used as the brief.
_roster_description() {
    local file="$1" fm
    fm=$(_roster_frontmatter "$file")
    local inline
    inline=$(_roster_scalar "$fm" description)
    if [[ -n "$inline" && "$inline" != "|" && "$inline" != ">" ]]; then
        printf '%s\n' "$inline"
        return
    fi
    # Block scalar: take the first non-blank indented line after `description:`.
    printf '%s\n' "$fm" | awk '
        BEGIN { in_desc = 0 }
        /^description:[ \t]*[|>]/ { in_desc = 1; next }
        in_desc {
            if ($0 ~ /^[a-zA-Z_]/) exit          # next top-level key
            if ($0 ~ /^[ \t]*$/)   next           # skip blank lines
            line = $0
            sub(/^[ \t]+/, "", line)
            print line
            exit
        }
    '
}

# _roster_agent_files
# Print each active agent .md path, sorted by name. Empty if none.
_roster_agent_files() {
    local f
    for f in .claude/agents/*.md; do
        [[ -e "$f" ]] || continue
        printf '%s\n' "$f"
    done
}

# _roster_role_stack
# Print space-separated active role names from the manifest, or empty.
_roster_role_stack() {
    [[ -f .sciagent/manifest.json ]] || return 0
    grep '"stack"' .sciagent/manifest.json \
        | sed 's/.*"stack":[[:space:]]*\[//' \
        | sed 's/\].*//' \
        | tr ',' '\n' \
        | sed 's/[" ]//g' \
        | grep -v '^$' \
        | tr '\n' ' ' \
        | sed 's/ $//'
}

# _roster_print_table
# Human-readable table: AGENT | DOMAIN | MODEL | DESCRIPTION (first 72 chars).
_roster_print_table() {
    printf '%-18s %-22s %-8s %s\n' AGENT DOMAIN MODEL "DESCRIPTION (first 72 chars)"
    local file fm name domain model desc
    while IFS= read -r file; do
        [[ -n "$file" ]] || continue
        fm=$(_roster_frontmatter "$file")
        name=$(_roster_scalar "$fm" name)
        [[ -n "$name" ]] || name="${file##*/}"; name="${name%.md}"
        domain=$(_roster_list "$fm" domain | paste -sd, -)
        model=$(_roster_scalar "$fm" model)
        desc=$(_roster_description "$file")
        desc="${desc:0:72}"
        printf '%-18s %-22s %-8s %s\n' "$name" "${domain:--}" "${model:--}" "$desc"
    done < <(_roster_agent_files)
}

# _roster_json_escape <string>
# Minimal JSON string escaping (backslash, quote).
_roster_json_escape() {
    local s="$1"
    s="${s//\\/\\\\}"
    s="${s//\"/\\\"}"
    printf '%s' "$s"
}

# _roster_json_array <items...>
# Print a JSON array of strings from newline-separated stdin.
_roster_json_array() {
    local first=1 item
    printf '['
    while IFS= read -r item; do
        [[ -n "$item" ]] || continue
        [[ $first -eq 1 ]] || printf ', '
        printf '"%s"' "$(_roster_json_escape "$item")"
        first=0
    done
    printf ']'
}

# _roster_print_json
# JSON manifest: generated_date, role_stack, agents[].
_roster_print_json() {
    local gen_date role_stack
    gen_date=$(date +%Y-%m-%d)
    role_stack=$(_roster_role_stack)

    printf '{\n'
    printf '  "generated_date": "%s",\n' "$gen_date"
    printf '  "role_stack": '
    printf '%s\n' "$role_stack" | tr ' ' '\n' | _roster_json_array
    printf ',\n'
    printf '  "agents": [\n'

    local file fm name model desc kind default_path first=1
    while IFS= read -r file; do
        [[ -n "$file" ]] || continue
        fm=$(_roster_frontmatter "$file")
        name=$(_roster_scalar "$fm" name)
        [[ -n "$name" ]] || { name="${file##*/}"; name="${name%.md}"; }
        model=$(_roster_scalar "$fm" model)
        desc=$(_roster_description "$file")
        kind=$(_roster_outputs_field "$fm" kind)
        default_path=$(_roster_outputs_field "$fm" default_path)

        [[ $first -eq 1 ]] || printf ',\n'
        first=0
        printf '    {\n'
        printf '      "name": "%s",\n' "$(_roster_json_escape "$name")"
        printf '      "model": "%s",\n' "$(_roster_json_escape "$model")"
        printf '      "domain": '
        _roster_list "$fm" domain | _roster_json_array
        printf ',\n'
        printf '      "description_brief": "%s",\n' "$(_roster_json_escape "$desc")"
        if [[ -n "$kind" || -n "$default_path" ]]; then
            printf '      "outputs": { "kind": "%s", "default_path": "%s" }\n' \
                "$(_roster_json_escape "$kind")" "$(_roster_json_escape "$default_path")"
        else
            printf '      "outputs": null\n'
        fi
        printf '    }'
    done < <(_roster_agent_files)

    printf '\n  ]\n}\n'
}

# cmd_roster [--json]
# Entry point. Reads `.claude/agents/` in CWD; returns 1 gracefully if absent.
cmd_roster() {
    local as_json=0
    while [[ $# -gt 0 ]]; do
        case "$1" in
            --json) as_json=1; shift ;;
            -*) echo "sciagent roster: unknown flag '$1'" >&2; return 1 ;;
            *)  echo "sciagent roster: unexpected argument '$1'" >&2; return 1 ;;
        esac
    done

    if [[ ! -d .claude/agents ]]; then
        echo "sciagent roster: no active role — run 'sciagent activate <role>'" >&2
        return 1
    fi

    if [[ $as_json -eq 1 ]]; then
        _roster_print_json
    else
        _roster_print_table
    fi
    return 0
}
