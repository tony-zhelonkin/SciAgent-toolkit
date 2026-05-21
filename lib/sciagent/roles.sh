# lib/sciagent/roles.sh — minimal YAML role parser (bash-only).
#
# Schema (per architecture §6): name, description, optional skills[],
# agents[], commands[], output_style. Strict YAML edge cases (quoting,
# nesting, anchors) are out of scope.
#
# Resolution: roles live in <toolkit>/roles/<name>.yaml. The toolkit root
# is the SCIAGENT_TOOLKIT env var (set by bin/sciagent) or the parent of
# the dir holding this file.

# shellcheck shell=bash

_roles_toolkit_root() {
    if [[ -n "${SCIAGENT_TOOLKIT:-}" ]]; then
        printf '%s\n' "$SCIAGENT_TOOLKIT"
        return
    fi
    local self_dir
    self_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    # self_dir = <toolkit>/lib/sciagent
    (cd "$self_dir/../.." && pwd)
}

role_path() {
    local name="$1"
    printf '%s/roles/%s.yaml\n' "$(_roles_toolkit_root)" "$name"
}

role_exists() {
    local name="$1"
    [[ -f "$(role_path "$name")" ]]
}

# Print a single top-level scalar (key: value). Empty if absent.
role_scalar() {
    local file="$1" key="$2"
    awk -v k="$key" '
        BEGIN { pat = "^" k ":[ \t]*" }
        $0 ~ pat {
            sub(pat, "")
            sub(/[ \t]*#.*$/, "")    # strip trailing comment
            sub(/[ \t]+$/, "")
            print
            exit
        }
    ' "$file"
}

# Print items in a top-level YAML array `<key>:` (list of `  - item` lines).
# Comments after `#` are stripped. Stops at the next top-level key.
role_array() {
    local file="$1" key="$2"
    awk -v k="$key" '
        BEGIN { in_block=0 }
        {
            line = $0
            # Detect top-level key (no leading whitespace, ends with colon)
            if (line ~ /^[a-zA-Z_][a-zA-Z0-9_]*:/) {
                if (in_block) in_block = 0
                if (line ~ ("^" k ":[ \t]*(#.*)?$")) {
                    in_block = 1
                    next
                }
            }
            if (in_block) {
                # Match `  - item` or `- item`, possibly with trailing comment.
                if (line ~ /^[ \t]*-[ \t]+/) {
                    sub(/^[ \t]*-[ \t]+/, "", line)
                    sub(/[ \t]*#.*$/, "", line)
                    sub(/[ \t]+$/, "", line)
                    if (line != "") print line
                }
            }
        }
    ' "$file"
}

# role_load <name>
# Emits tagged lines on stdout:
#   SKILL <name>
#   AGENT <name>
#   COMMAND <name>
#   OUTPUT_STYLE <name>
role_load() {
    local name="$1"
    local file
    file=$(role_path "$name")
    if [[ ! -f "$file" ]]; then
        echo "role_load: $name not found at $file" >&2
        return 1
    fi
    local item
    while IFS= read -r item; do
        [[ -n "$item" ]] && printf 'SKILL %s\n' "$item"
    done < <(role_array "$file" skills)
    while IFS= read -r item; do
        [[ -n "$item" ]] && printf 'AGENT %s\n' "$item"
    done < <(role_array "$file" agents)
    while IFS= read -r item; do
        [[ -n "$item" ]] && printf 'COMMAND %s\n' "$item"
    done < <(role_array "$file" commands)
    local style
    style=$(role_scalar "$file" output_style)
    [[ -n "$style" ]] && printf 'OUTPUT_STYLE %s\n' "$style"
    return 0
}

role_description() {
    role_scalar "$(role_path "$1")" description
}

# ---------------------------------------------------------------------------
# system-prompts resolution
# ---------------------------------------------------------------------------
#
# A role's `output_style: <name>` field is a logical identifier resolved
# against the frontmatter `name:` field of files in `system-prompts/*.md`.
# This decouples the role spec from filenames — files can be renamed
# without breaking roles, and a single style can live in any filename.
#
# Validation lives here so that drift between a role's declared style and
# the canonical source files is caught at activation time, before any
# symlink or manifest mutation. Otherwise an unresolved style would
# silently fall back to a hard-coded default and confuse the user.

# _read_frontmatter_name <file>
# Print the value of the `name:` field inside the leading `---` ... `---`
# YAML frontmatter block. Empty if none.
_read_frontmatter_name() {
    awk '
        /^---[ \t]*$/ { fm++; if (fm == 2) exit; next }
        fm == 1 && /^name:[ \t]/ {
            sub(/^name:[ \t]*/, "")
            sub(/[ \t]*#.*$/, "")
            sub(/[ \t]+$/, "")
            print
            exit
        }
    ' "$1"
}

# system_prompt_path <name>
# Resolve a system-prompt source file by frontmatter name. Print the
# absolute path on stdout. Exit 1 if no file matches.
system_prompt_path() {
    local name="$1"
    local dir
    dir="$(_roles_toolkit_root)/system-prompts"
    [[ -d "$dir" ]] || return 1
    local f n
    for f in "$dir"/*.md; do
        [[ -f "$f" ]] || continue
        n=$(_read_frontmatter_name "$f")
        if [[ "$n" == "$name" ]]; then
            printf '%s\n' "$f"
            return 0
        fi
    done
    return 1
}

# system_prompt_inventory
# List every system-prompt source as a `<name>\t<basename>` line. Used by
# the validation error message and by `sciagent list system-prompts`.
system_prompt_inventory() {
    local dir
    dir="$(_roles_toolkit_root)/system-prompts"
    [[ -d "$dir" ]] || return 1
    local f n
    for f in "$dir"/*.md; do
        [[ -f "$f" ]] || continue
        n=$(_read_frontmatter_name "$f")
        printf '%s\t%s\n' "${n:-<missing-frontmatter-name>}" "${f##*/}"
    done
}
