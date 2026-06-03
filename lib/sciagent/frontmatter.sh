# lib/sciagent/frontmatter.sh — pure YAML frontmatter parser.
#
# Leaf module: no sciagent deps. Functions lifted and generalised from the old
# roster.sh so other modules (status.sh, list) can read agent/skill frontmatter
# without pulling in roster.
#
# Usage pattern (here-string or pipe):
#   fm=$(_fm_extract path/to/file.md)
#   _fm_scalar description <<< "$fm"
#   _fm_list   domain      <<< "$fm"
#   _fm_nested_scalar outputs default_path <<< "$fm"

# shellcheck shell=bash

# _fm_extract <file>
# Print the raw lines between the first and second `---` fences.
_fm_extract() {
    awk '
        BEGIN { in_fm = 0 }
        /^---[ \t]*$/ { in_fm++; next }
        in_fm == 1 { print; next }
        in_fm > 1  { exit }
    ' "$1"
}

# _fm_scalar <key>
# Read frontmatter from stdin; print the top-level scalar value for <key>.
# Empty if absent. Strips inline comments, quotes, and trailing whitespace.
_fm_scalar() {
    local key="$1"
    awk -v k="$key" '
        BEGIN { pat = "^" k ":[ \t]*" }
        $0 ~ pat {
            sub(pat, "")
            sub(/[ \t]*#.*$/, "")
            sub(/^["'"'"']/, ""); sub(/["'"'"']$/, "")
            sub(/[ \t]+$/, "")
            print
            exit
        }
    '
}

# _fm_list <key>
# Read frontmatter from stdin; print items of a top-level YAML list
# (`key:` then `  - item` lines) one per line. Stops at the next
# top-level key. Strips inline comments and trailing whitespace.
_fm_list() {
    local key="$1"
    awk -v k="$key" '
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

# _fm_nested_scalar <parent> <child>
# Read frontmatter from stdin; print the scalar value nested one level
# under <parent> (e.g. `outputs.default_path`). Empty if either level
# is absent.
_fm_nested_scalar() {
    local parent="$1" child="$2"
    awk -v p="$parent" -v f="$child" '
        BEGIN { in_parent = 0 }
        {
            if ($0 ~ ("^" p ":[ \t]*$")) { in_parent = 1; next }
            if ($0 ~ /^[a-zA-Z_]/)       { in_parent = 0 }
            if (in_parent && $0 ~ ("^[ \t]+" f ":[ \t]*")) {
                sub(("^[ \t]+" f ":[ \t]*"), "")
                sub(/[ \t]*#.*$/, "")
                sub(/[ \t]+$/, "")
                print
                exit
            }
        }
    '
}

# _fm_description <file>
# Print a one-line description from a markdown-with-frontmatter file.
# Handles both inline (`description: text`) and block-scalar
# (`description: |`) forms; for block, the first non-blank content line
# is used as the brief. Empty if absent.
_fm_description() {
    local file="$1" fm inline
    fm=$(_fm_extract "$file")
    inline=$(_fm_scalar description <<< "$fm")
    if [[ -n "$inline" && "$inline" != "|" && "$inline" != ">" ]]; then
        printf '%s\n' "$inline"
        return
    fi
    # Block scalar: take the first non-blank indented line after `description:`.
    awk '
        BEGIN { in_desc = 0 }
        /^description:[ \t]*[|>]/ { in_desc = 1; next }
        in_desc {
            if ($0 ~ /^[a-zA-Z_]/) exit
            if ($0 ~ /^[ \t]*$/)   next
            line = $0
            sub(/^[ \t]+/, "", line)
            print line
            exit
        }
    ' <<< "$fm"
}
