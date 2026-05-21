# lib/sciagent/symlinks.sh — dual-track symlink + manifest helpers.
#
# Manifest format: real JSON (pretty-printed, 2-space indent). Schema v1:
#
#   {
#     "version": 1,
#     "stack": ["base", "reviewer"],
#     "symlinks": ["/abs/or/relative/path"],
#     "injected": [{"overlay": "_injected", "skill": "obsidian-vignette"}],
#     "block_hash": "abc123..."
#   }
#
# jq is used when available for reading; when absent, a targeted bash fallback
# uses grep+sed against this well-known flat structure. The fallback does NOT
# attempt general JSON parsing — it only supports this specific schema.
# Writing always uses the bash _json_escape helper (no jq dependency).

# shellcheck shell=bash

_MANIFEST_PATH=".sciagent/manifest.json"

# ---------------------------------------------------------------------------
# resolve_canonical <category> <name>
# Prints the absolute path to the canonical source for a skill/agent/command.
# Categories:
#   skills    — returns $SCIAGENT_TOOLKIT/skills/<name>  (must be a directory
#               containing SKILL.md; not searched recursively)
#   agents    — recursive search for <name>.md under $SCIAGENT_TOOLKIT/agents/
#   commands  — recursive search for <name>.md under $SCIAGENT_TOOLKIT/commands/
# Fails (rc=1) with a stderr message if zero matches or multiple matches.
# Uniqueness is enforced by tests/test_no_duplicate_basenames.sh.
# ---------------------------------------------------------------------------
resolve_canonical() {
    local category="$1"
    local name="$2"
    local root="$SCIAGENT_TOOLKIT/$category"

    case "$category" in
        skills)
            local p="$root/$name"
            if [[ -d "$p" && -f "$p/SKILL.md" ]]; then
                printf '%s\n' "$p"
                return 0
            fi
            echo "skill '$name' not found under skills/ ($p)" >&2
            return 1
            ;;
        agents|commands)
            local -a matches=()
            local f
            while IFS= read -r f; do
                matches+=("$f")
            done < <(find "$root" -type f -name "${name}.md" -not -path '*/.*' 2>/dev/null)
            local n=${#matches[@]}
            if (( n == 0 )); then
                echo "${category%s} '$name' not found under ${category}/" >&2
                return 1
            fi
            if (( n > 1 )); then
                echo "duplicate basename '$name' under ${category}/: ${matches[*]}" >&2
                return 1
            fi
            printf '%s\n' "${matches[0]}"
            return 0
            ;;
        *)
            echo "resolve_canonical: unknown category '$category'" >&2
            return 1
            ;;
    esac
}

# ---------------------------------------------------------------------------
# JSON helpers
# ---------------------------------------------------------------------------

# _json_escape <string>
# Escapes a string value for safe embedding between JSON double-quotes.
# Handles: backslash, double-quote, newline, carriage-return, tab.
_json_escape() {
    local s="$1"
    s="${s//\\/\\\\}"     # \ → \\
    s="${s//\"/\\\"}"     # " → \"
    s="${s//$'\n'/\\n}"   # LF → \n
    s="${s//$'\r'/\\r}"   # CR → \r
    s="${s//$'\t'/\\t}"   # TAB → \t
    printf '%s' "$s"
}

# ---------------------------------------------------------------------------
# Manifest write helpers (staging buffer)
# ---------------------------------------------------------------------------

_manifest_staging=""
_manifest_stack_val=""
declare -a _manifest_syms=()
declare -a _manifest_injected_overlays=()
declare -a _manifest_injected_skills=()

manifest_begin() {
    local stack="$1"   # space-separated role names
    _manifest_staging=$(mktemp)
    _manifest_stack_val="$stack"
    _manifest_syms=()
    _manifest_injected_overlays=()
    _manifest_injected_skills=()
}

_manifest_record_symlink() {
    [[ -n "$_manifest_staging" ]] || return 0
    _manifest_syms+=("$1")
}

_manifest_record_injected() {
    local overlay="$1" skill="$2"
    _manifest_injected_overlays+=("$overlay")
    _manifest_injected_skills+=("$skill")
}

# _manifest_write_json <block_hash>
# Serialise accumulated state into $_manifest_staging as pretty JSON.
_manifest_write_json() {
    local block_hash="$1"
    {
        printf '{\n'
        printf '  "version": 1,\n'

        # stack array
        printf '  "stack": ['
        local first=1 r
        for r in $_manifest_stack_val; do
            (( first )) && first=0 || printf ', '
            printf '"%s"' "$(_json_escape "$r")"
        done
        printf '],\n'

        # symlinks array
        printf '  "symlinks": ['
        first=1
        local s
        local _nsyms=${#_manifest_syms[@]}
        if (( _nsyms > 0 )); then
            for s in "${_manifest_syms[@]}"; do
                (( first )) && first=0 || printf ',\n              '
                printf '"%s"' "$(_json_escape "$s")"
            done
        fi
        printf '],\n'

        # injected array
        printf '  "injected": ['
        first=1
        local _ninj=${#_manifest_injected_overlays[@]}
        if (( _ninj > 0 )); then
            local i
            for (( i=0; i<_ninj; i++ )); do
                (( first )) && first=0 || printf ','
                printf '\n    {"overlay": "%s", "skill": "%s"}' \
                    "$(_json_escape "${_manifest_injected_overlays[$i]}")" \
                    "$(_json_escape "${_manifest_injected_skills[$i]}")"
            done
            printf '\n  '
        fi
        printf '],\n'

        printf '  "block_hash": "%s"\n' "$(_json_escape "$block_hash")"
        printf '}\n'
    } > "$_manifest_staging"
}

manifest_finalize() {
    local block_hash="$1"
    [[ -n "$_manifest_staging" ]] || return 1
    _manifest_write_json "$block_hash"
    mkdir -p "$(dirname "$_MANIFEST_PATH")"
    mv "$_manifest_staging" "$_MANIFEST_PATH"
    _manifest_staging=""
}

# ---------------------------------------------------------------------------
# symlink_create_dual
# ---------------------------------------------------------------------------

# symlink_create_dual <category> <name> <canonical_path>
# Categories: skills (dir symlink), agents/commands/output-styles (file symlink).
# Records both .claude/<cat>/<name> and .agents/<cat>/<name>.
symlink_create_dual() {
    local category="$1"
    local name="$2"
    local canonical="$3"

    local claude_path agents_path
    case "$category" in
        skills)
            claude_path=".claude/skills/$name"
            agents_path=".agents/skills/$name"
            ;;
        agents)
            claude_path=".claude/agents/${name}.md"
            agents_path=".agents/agents/${name}.md"
            ;;
        commands)
            claude_path=".claude/commands/${name}.md"
            agents_path=".agents/commands/${name}.md"
            ;;
        output-styles)
            # Claude-only. No .agents/ mirror.
            claude_path=".claude/output-styles/${name}.md"
            agents_path=""
            ;;
        *)
            echo "symlink_create_dual: unknown category '$category'" >&2
            return 1
            ;;
    esac

    mkdir -p "$(dirname "$claude_path")"
    ln -sfn "$canonical" "$claude_path"
    _manifest_record_symlink "$claude_path"

    if [[ -n "$agents_path" ]]; then
        mkdir -p "$(dirname "$agents_path")"
        ln -sfn "$canonical" "$agents_path"
        _manifest_record_symlink "$agents_path"
    fi
}

# ---------------------------------------------------------------------------
# symlink_teardown_all — remove only the symlinks recorded in our manifest.
# ---------------------------------------------------------------------------
symlink_teardown_all() {
    [[ -f "$_MANIFEST_PATH" ]] || return 0
    local -a syms
    readarray -t syms < <(manifest_symlinks)
    local path
    for path in "${syms[@]+"${syms[@]}"}"; do
        [[ -z "$path" ]] && continue
        if [[ ! -e "$path" && ! -L "$path" ]]; then
            echo "warning: $path already gone, skipping" >&2
            continue
        fi
        if [[ ! -L "$path" ]]; then
            echo "warning: $path is not a symlink, skipping" >&2
            continue
        fi
        rm "$path"
    done
    rm -f "$_MANIFEST_PATH"
    rmdir .sciagent 2>/dev/null || true
    # Best-effort: clean empty .claude/* and .agents/* dirs we created.
    for d in .claude/skills .claude/agents .claude/commands .claude/output-styles \
             .agents/skills .agents/agents .agents/commands; do
        [[ -d "$d" ]] && rmdir "$d" 2>/dev/null || true
    done
    rmdir .agents 2>/dev/null || true
}

# ---------------------------------------------------------------------------
# Manifest readers — jq when available, bash fallback otherwise.
#
# Fallback strategy: the schema has well-known flat keys at the top level.
# We grep for the key and extract the value with sed. This is intentionally
# narrow — it only handles this exact schema, not arbitrary JSON.
# ---------------------------------------------------------------------------

manifest_exists() {
    [[ -f "$_MANIFEST_PATH" ]]
}

# manifest_stack — print space-separated role names from "stack" array.
manifest_stack() {
    [[ -f "$_MANIFEST_PATH" ]] || return 1
    if command -v jq >/dev/null 2>&1; then
        jq -r '.stack | join(" ")' "$_MANIFEST_PATH"
    else
        # Fallback: "stack": ["base", "reviewer"]  — single-line as written.
        grep '"stack"' "$_MANIFEST_PATH" \
            | sed 's/.*"stack":[[:space:]]*\[//' \
            | sed 's/\].*//' \
            | tr -d '"' \
            | tr ',' ' ' \
            | tr -s ' ' \
            | sed 's/^[[:space:]]*//; s/[[:space:]]*$//'
    fi
}

# manifest_symlinks — print one symlink path per line from "symlinks" array.
manifest_symlinks() {
    [[ -f "$_MANIFEST_PATH" ]] || return 1
    if command -v jq >/dev/null 2>&1; then
        jq -r '.symlinks[]' "$_MANIFEST_PATH"
    else
        # Fallback: the writer formats the array as
        #     "symlinks": ["first",
        #                  "second",
        #                  ...
        #                  "last"],
        # so the first symlink shares a line with the `"symlinks":` key and the
        # last symlink shares a line with the `]` closer. An earlier version of
        # this awk used `next` on those framing lines and silently dropped both
        # the first and last symlink, leaving them un-tracked for teardown
        # (caused test_deactivate_full.sh to leave a `.claude/skills/s_a`
        # behind). Now we extract every quoted string on lines inside the
        # `"symlinks": [ ... ]` block — including the opener and closer.
        awk '
            /"symlinks"/ { in_s=1 }
            in_s {
                line = $0
                while (match(line, /"[^"]+"/)) {
                    v = substr(line, RSTART+1, RLENGTH-2)
                    if (v != "symlinks") print v
                    line = substr(line, RSTART+RLENGTH)
                }
            }
            in_s && /\]/ { in_s=0 }
        ' "$_MANIFEST_PATH"
    fi
}

# manifest_injected — print one "overlay skill" pair per line.
manifest_injected() {
    [[ -f "$_MANIFEST_PATH" ]] || return 1
    if command -v jq >/dev/null 2>&1; then
        jq -r '.injected[] | (.overlay + " " + .skill)' "$_MANIFEST_PATH"
    else
        # Fallback: each injected entry spans one line:
        #   {"overlay": "X", "skill": "Y"}
        grep '"overlay"' "$_MANIFEST_PATH" | while IFS= read -r line; do
            local ov sk
            ov=$(printf '%s' "$line" | sed 's/.*"overlay":[[:space:]]*"\([^"]*\)".*/\1/')
            sk=$(printf '%s' "$line" | sed 's/.*"skill":[[:space:]]*"\([^"]*\)".*/\1/')
            printf '%s %s\n' "$ov" "$sk"
        done
    fi
}

# manifest_block_hash — print stored block hash.
manifest_block_hash() {
    [[ -f "$_MANIFEST_PATH" ]] || return 1
    if command -v jq >/dev/null 2>&1; then
        jq -r '.block_hash' "$_MANIFEST_PATH"
    else
        grep '"block_hash"' "$_MANIFEST_PATH" \
            | sed 's/.*"block_hash":[[:space:]]*"\([^"]*\)".*/\1/'
    fi
}

# manifest_update_stack <new-stack>  — rewrite "stack" value in-place.
manifest_update_stack() {
    local new_stack="$1"   # space-separated
    [[ -f "$_MANIFEST_PATH" ]] || return 1

    # Re-read current manifest fields, rebuild JSON with updated stack.
    local old_syms old_inj old_hash
    old_syms=$(manifest_symlinks)
    old_inj=$(manifest_injected)
    old_hash=$(manifest_block_hash)

    # Reconstruct via staging.
    _manifest_stack_val="$new_stack"
    _manifest_syms=()
    _manifest_injected_overlays=()
    _manifest_injected_skills=()

    local p
    while IFS= read -r p; do
        [[ -n "$p" ]] && _manifest_syms+=("$p")
    done <<< "$old_syms"

    local ov sk
    while read -r ov sk; do
        [[ -n "$ov" ]] && _manifest_injected_overlays+=("$ov") && _manifest_injected_skills+=("$sk")
    done <<< "$old_inj"

    _manifest_staging=$(mktemp)
    _manifest_write_json "$old_hash"
    mv "$_manifest_staging" "$_MANIFEST_PATH"
    _manifest_staging=""
}

# manifest_update_block_hash <hash>  — rewrite "block_hash" in-place.
manifest_update_block_hash() {
    local new_hash="$1"
    [[ -f "$_MANIFEST_PATH" ]] || return 1

    local old_stack old_syms old_inj
    old_stack=$(manifest_stack)
    old_syms=$(manifest_symlinks)
    old_inj=$(manifest_injected)

    _manifest_stack_val="$old_stack"
    _manifest_syms=()
    _manifest_injected_overlays=()
    _manifest_injected_skills=()

    local p
    while IFS= read -r p; do
        [[ -n "$p" ]] && _manifest_syms+=("$p")
    done <<< "$old_syms"

    local ov sk
    while read -r ov sk; do
        [[ -n "$ov" ]] && _manifest_injected_overlays+=("$ov") && _manifest_injected_skills+=("$sk")
    done <<< "$old_inj"

    _manifest_staging=$(mktemp)
    _manifest_write_json "$new_hash"
    mv "$_manifest_staging" "$_MANIFEST_PATH"
    _manifest_staging=""
}

# manifest_append_symlinks_and_injected <sym1> <sym2> ... -- <overlay> <skill>
# Appends new symlinks and one injected entry to the existing manifest in-place.
# Call signature: manifest_append_inject <claude_sym> <agents_sym> <overlay> <skill>
manifest_append_inject() {
    local claude_sym="$1" agents_sym="$2" target_overlay="$3" skill="$4"
    [[ -f "$_MANIFEST_PATH" ]] || return 1

    local old_stack old_hash
    old_stack=$(manifest_stack)
    old_hash=$(manifest_block_hash)

    _manifest_stack_val="$old_stack"
    _manifest_syms=()
    _manifest_injected_overlays=()
    _manifest_injected_skills=()

    local p
    while IFS= read -r p; do
        [[ -n "$p" ]] && _manifest_syms+=("$p")
    done < <(manifest_symlinks)

    local ov sk
    while read -r ov sk; do
        [[ -n "$ov" ]] && _manifest_injected_overlays+=("$ov") && _manifest_injected_skills+=("$sk")
    done < <(manifest_injected)

    # Append new entries.
    _manifest_syms+=("$claude_sym" "$agents_sym")
    _manifest_injected_overlays+=("$target_overlay")
    _manifest_injected_skills+=("$skill")

    _manifest_staging=$(mktemp)
    _manifest_write_json "$old_hash"
    mv "$_manifest_staging" "$_MANIFEST_PATH"
    _manifest_staging=""
}
