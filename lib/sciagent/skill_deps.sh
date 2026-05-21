# lib/sciagent/skill_deps.sh — skill frontmatter dependency parser & resolver.
#
# Public functions:
#   skill_frontmatter_path <name>
#       Print absolute path to skills/<name>/SKILL.md if it exists. Exit 1
#       with stderr message if the skill or its SKILL.md is missing.
#
#   skill_read_requires <name>
#       Print direct requires of <name>, one per line. Looks under
#       `metadata.requires:` (canonical, per ADR-0002 §3.1). If the field is
#       absent or empty, prints nothing and returns 0. The frontmatter parser
#       is intentionally narrow — it only handles the well-known schema used
#       by SciAgent-toolkit skills (`---\n...\n---` block at file head,
#       metadata as a nested mapping, requires as an inline `[a, b]` list
#       OR a block list of `  - item` lines).
#
#   skill_resolve_transitive <name1> [<name2> ...]
#       DFS from each input name. Emits the full closure on stdout, one skill
#       per line, in topological (post-order) order with duplicates removed.
#       Hard-fails (exit 1, no partial output) on:
#         - a cycle in the requires graph
#         - a `requires:` target with no matching skills/<target>/SKILL.md
#
# Design notes:
#   - Bash-only; no jq, no python in the hot path.
#   - No global state survives a call (associative arrays are scoped via
#     `declare -A` inside the resolver function).
#   - Pre-mutation contract: callers (activate.sh) MUST invoke
#     skill_resolve_transitive BEFORE any filesystem mutation (per ADR §4.6).

# shellcheck shell=bash

# ---------------------------------------------------------------------------
# skill_frontmatter_path <name>
# ---------------------------------------------------------------------------
skill_frontmatter_path() {
    local name="$1"
    local root
    if [[ -n "${SCIAGENT_TOOLKIT:-}" ]]; then
        root="$SCIAGENT_TOOLKIT"
    else
        local self_dir
        self_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
        root="$(cd "$self_dir/../.." && pwd)"
    fi
    local path="$root/skills/$name/SKILL.md"
    if [[ ! -f "$path" ]]; then
        echo "skill_deps: skill '$name' has no SKILL.md at $path" >&2
        return 1
    fi
    printf '%s\n' "$path"
}

# ---------------------------------------------------------------------------
# _skill_extract_frontmatter <file>
# Print just the YAML frontmatter (between the leading `---` markers).
# Empty output if there is no frontmatter.
# ---------------------------------------------------------------------------
_skill_extract_frontmatter() {
    awk '
        BEGIN { fm = 0 }
        /^---[ \t]*$/ {
            fm++
            if (fm == 1) next
            if (fm == 2) exit
        }
        fm == 1 { print }
    ' "$1"
}

# ---------------------------------------------------------------------------
# skill_read_requires <name>
# Emit direct `metadata.requires` entries, one per line.
# ---------------------------------------------------------------------------
skill_read_requires() {
    local name="$1"
    local file
    file=$(skill_frontmatter_path "$name") || return 1

    # Strategy: pull frontmatter, then walk it. We only honour
    # `metadata.requires:` (canonical, per ADR §3.1 / Decision 4).
    # Two formats supported:
    #   metadata:
    #     requires: [a, b, c]          # inline list
    #   metadata:
    #     requires:                     # block list
    #       - a
    #       - b
    _skill_extract_frontmatter "$file" | awk '
        BEGIN { in_meta = 0; in_req = 0 }
        # Top-level key resets context. A "top-level" key has no leading
        # whitespace and ends with a colon.
        /^[A-Za-z_][A-Za-z0-9_-]*:/ {
            in_meta = ($0 ~ /^metadata:[ \t]*(#.*)?$/)
            in_req = 0
            next
        }
        # Inside metadata: look for `  requires:` line (indented one level).
        in_meta && /^[ \t]+requires:/ {
            # Inline form?
            line = $0
            sub(/^[ \t]+requires:[ \t]*/, "", line)
            sub(/[ \t]*#.*$/, "", line)
            sub(/[ \t]+$/, "", line)
            if (line ~ /^\[.*\]$/) {
                # Inline list: strip [ ], split on commas.
                gsub(/^\[[ \t]*/, "", line)
                gsub(/[ \t]*\]$/, "", line)
                n = split(line, parts, /[ \t]*,[ \t]*/)
                for (i = 1; i <= n; i++) {
                    item = parts[i]
                    gsub(/^["'\'']|["'\'']$/, "", item)
                    if (item != "") print item
                }
                in_req = 0
                next
            }
            # Block form: subsequent `    - item` lines.
            in_req = 1
            next
        }
        # Another indented metadata key terminates the requires block.
        in_req && /^[ \t]+[A-Za-z_][A-Za-z0-9_-]*:/ {
            in_req = 0
            next
        }
        in_req && /^[ \t]+-[ \t]+/ {
            item = $0
            sub(/^[ \t]+-[ \t]+/, "", item)
            sub(/[ \t]*#.*$/, "", item)
            sub(/[ \t]+$/, "", item)
            gsub(/^["'\'']|["'\'']$/, "", item)
            if (item != "") print item
        }
    '
}

# ---------------------------------------------------------------------------
# skill_resolve_transitive <name>...
# DFS with grey/black coloring. Post-order emission = topological order
# (dependencies appear before dependants in the output).
# ---------------------------------------------------------------------------
skill_resolve_transitive() {
    local -a roots=( "$@" )
    declare -A _sd_state=()    # name -> "grey" | "black"
    local -a _sd_order=()

    local r
    for r in "${roots[@]}"; do
        _sd_visit "$r" || return 1
    done

    local n
    for n in "${_sd_order[@]+"${_sd_order[@]}"}"; do
        printf '%s\n' "$n"
    done
}

# Internal DFS visitor. Uses caller's _sd_state and _sd_order arrays.
_sd_visit() {
    local name="$1"
    local state="${_sd_state[$name]:-}"
    if [[ "$state" == "black" ]]; then
        return 0
    fi
    if [[ "$state" == "grey" ]]; then
        echo "skill_deps: cycle detected in requires graph at '$name'" >&2
        return 1
    fi

    # Verify the skill exists before we recurse — gives a precise error.
    if ! skill_frontmatter_path "$name" >/dev/null; then
        echo "skill_deps: required skill '$name' not found" >&2
        return 1
    fi

    _sd_state[$name]="grey"

    local dep
    local -a deps=()
    while IFS= read -r dep; do
        [[ -n "$dep" ]] && deps+=("$dep")
    done < <(skill_read_requires "$name")

    for dep in "${deps[@]+"${deps[@]}"}"; do
        if ! _sd_visit "$dep"; then
            echo "skill_deps:   (while resolving '$name' -> '$dep')" >&2
            return 1
        fi
    done

    _sd_state[$name]="black"
    _sd_order+=("$name")
    return 0
}
