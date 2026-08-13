# lib/sciagent/craft.sh — SCIAGENT:CRAFT managed-block renderer.
#
# The CRAFT block carries the owner's standing craft conventions (figure
# legibility, results placement, README adjacency, planning decomposition,
# reproducibility). Its single source of truth is <toolkit>/craft.yaml.
# Rendered into AGENTS.md by `sciagent craft` (see craft_verb.sh). SHA1
# drift-detected and idempotently re-rendered via block.sh (id=CRAFT).
#
# Depends on block.sh (block_write/block_remove with id=CRAFT). All managed-
# block framing + hashing stays in block.sh; this module only assembles body
# text, so the test_block_marker_boundary invariant is preserved.

# shellcheck shell=bash

CRAFT_BLOCK_ID="CRAFT"

# _craft_yaml_path — absolute path to the craft SSOT.
_craft_yaml_path() {
    printf '%s/craft.yaml' "${SCIAGENT_TOOLKIT:-.}"
}

# craft_version — emit the integer `version:` from craft.yaml (default 1).
craft_version() {
    local f
    f=$(_craft_yaml_path)
    [[ -f "$f" ]] || { echo 1; return; }
    local v
    v=$(awk '/^version:[[:space:]]*[0-9]+/ { gsub(/[^0-9]/, "", $2); print $2; exit }' "$f")
    echo "${v:-1}"
}

# _craft_render_body — assemble the CRAFT block body from craft.yaml:
#   * read the `floors:` map (key -> value)
#   * read the trailing `body: |` block scalar (must be the last top-level key)
#   * substitute {{key}} tokens with floor values
# Echoes the rendered body. Returns 1 when craft.yaml is absent.
_craft_render_body() {
    local f
    f=$(_craft_yaml_path)
    [[ -f "$f" ]] || return 1

    # Floors: "  key: value" lines under a top-level "floors:" key.
    local -A _floors=()
    local line k v
    while IFS='=' read -r k v; do
        [[ -n "$k" ]] && _floors["$k"]="$v"
    done < <(awk '
        /^[^[:space:]]/        { in_f=0 }
        /^floors:[[:space:]]*$/ { in_f=1; next }
        in_f && /^[[:space:]]+[A-Za-z0-9_]+:/ {
            entry=$0
            sub(/^[[:space:]]+/, "", entry)
            key=entry; sub(/:.*/, "", key)
            val=entry; sub(/^[^:]+:[[:space:]]*/, "", val)
            sub(/[[:space:]]*#.*/, "", val)     # strip inline comment
            sub(/[[:space:]]+$/, "", val)
            print key "=" val
        }
    ' "$f")

    # Body: the trailing "body: |" block scalar, de-indented by its own indent.
    local body
    body=$(awk '
        seen { raw[n++]=$0; next }
        /^body:[[:space:]]*\|[[:space:]]*$/ { seen=1 }
        END {
            indent=-1
            for (i=0; i<n; i++) {
                if (raw[i] ~ /[^[:space:]]/) {
                    match(raw[i], /^[[:space:]]*/)
                    indent=RLENGTH
                    break
                }
            }
            if (indent<0) indent=0
            # Drop trailing blank lines for a tight block.
            last=n-1
            while (last>=0 && raw[last] !~ /[^[:space:]]/) last--
            for (i=0; i<=last; i++) print substr(raw[i], indent+1)
        }
    ' "$f")

    # Substitute {{key}} tokens with floor values.
    for k in "${!_floors[@]}"; do
        body="${body//\{\{$k\}\}/${_floors[$k]}}"
    done
    printf '%s\n' "$body"
}

# craft_render_and_write [<agents-file>]
# Render the CRAFT block from craft.yaml and write it (id=CRAFT) via block.sh.
# No-op (return 0) when craft.yaml is absent, so toolkits/harnesses without a
# craft SSOT are unaffected. Returns block_write status otherwise.
craft_render_and_write() {
    local file="${1:-AGENTS.md}"
    local body
    body=$(_craft_render_body) || return 0     # absent craft.yaml -> skip
    [[ -n "$body" ]] || return 0
    block_write "$file" "$body" "$CRAFT_BLOCK_ID" || {
        echo "sciagent: failed to write CRAFT block to $file" >&2
        return 1
    }
}
