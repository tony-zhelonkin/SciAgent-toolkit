# lib/scio/claude_settings.sh — Claude Code project guardrail hooks.

# shellcheck shell=bash

_CLAUDE_PROJECT_SETTINGS=".claude/settings.json"
_CLAUDE_HOOKS_DIR=".claude/hooks"
_CLAUDE_HOOKS_STATE_DIR=".scio/hook_state"
_CLAUDE_HOOK_SETTINGS_STATE=".scio/hook_settings.state"

# Register the shipped hooks while preserving every unrelated project setting.
_claude_settings_ensure_hook_registration() {
    local dst="$_CLAUDE_PROJECT_SETTINGS"
    local src="$SCIO_TOOLKIT/templates/project/_common/.claude/settings.json.template"
    local state="$_CLAUDE_HOOK_SETTINGS_STATE"
    [[ -f "$src" ]] || return 0

    mkdir -p "$(dirname "$dst")" "$(dirname "$state")"

    if [[ ! -f "$dst" ]]; then
        cp "$src" "$dst" || { echo "scio: failed to write $dst" >&2; return 1; }
        if command -v jq >/dev/null 2>&1; then
            jq -n \
                --arg h "$(scio_sha1_file "$dst")" \
                --slurpfile template "$src" \
                '{tag:"created", hash:$h, added:$template[0].hooks}' > "$state"
        else
            printf '{"tag":"created","hash":"%s"}\n' "$(scio_sha1_file "$dst")" > "$state"
        fi
        echo "wrote: $dst"
        return 0
    fi

    if ! command -v jq >/dev/null 2>&1; then
        echo "scio: warning — jq not found; cannot register hooks in $dst" >&2
        return 0
    fi

    local added merged
    added=$(mktemp)
    merged=$(mktemp)

    if ! jq -n --slurpfile wanted "$src" --slurpfile current "$dst" '
        reduce (($wanted[0].hooks // {}) | to_entries[]) as $event ({};
          ($event.value | map(
            . as $entry
            | select(((($current[0].hooks[$event.key] // [])
              | any(. == $entry)) | not))
          )) as $missing
          | if ($missing | length) > 0
            then .[$event.key] = $missing
            else .
            end
        )
    ' > "$added" 2>/dev/null; then
        rm -f "$added" "$merged"
        echo "scio: warning — could not read hook settings from $dst (invalid JSON?)" >&2
        return 0
    fi

    if jq -e 'length == 0' "$added" >/dev/null 2>&1; then
        rm -f "$added" "$merged"
        return 0
    fi

    if ! jq -s '
        .[0].hooks as $wanted | .[1]
        | .hooks = reduce ($wanted | to_entries[]) as $event (.hooks // {};
            .[$event.key] = reduce $event.value[] as $entry (.[$event.key] // [];
              if any(.[]; . == $entry) then . else . + [$entry] end
            )
          )
    ' "$src" "$dst" > "$merged" 2>/dev/null || [[ ! -s "$merged" ]]; then
        rm -f "$added" "$merged"
        echo "scio: warning — could not register hooks in $dst (invalid JSON?)" >&2
        return 0
    fi

    if [[ ! -f "$state" ]]; then
        cp "$dst" "${state}.orig" || { rm -f "$added" "$merged"; return 1; }
        cp "$merged" "$dst" || { rm -f "$added" "$merged"; return 1; }
        jq -n \
            --arg h "$(scio_sha1_file "$dst")" \
            --slurpfile added "$added" \
            '{tag:"existed-registered", written_hash:$h, added:$added[0]}' > "$state"
    else
        cp "$merged" "$dst" || { rm -f "$added" "$merged"; return 1; }
        local next_state
        next_state=$(mktemp)
        if jq --slurpfile extra "$added" '
            .added = reduce (($extra[0] // {}) | to_entries[]) as $event (.added // {};
              .[$event.key] = reduce $event.value[] as $entry (.[$event.key] // [];
                if any(.[]; . == $entry) then . else . + [$entry] end
              )
            )
        ' "$state" > "$next_state" 2>/dev/null; then
            cp "$next_state" "$state"
        fi
        rm -f "$next_state"
    fi

    rm -f "$added" "$merged"
    echo "updated: $dst (registered project guardrail hooks)"
}

# Materialize both project guardrail bodies and register them in settings.json.
claude_settings_ensure_hooks() {
    local tpl_dir="$SCIO_TOOLKIT/templates/project/_common/$_CLAUDE_HOOKS_DIR"
    [[ -d "$tpl_dir" ]] || return 0

    local src dst base found=false
    for src in "$tpl_dir"/*.sh.template; do
        [[ -f "$src" ]] || continue
        found=true
        base="$(basename "${src%.template}")"
        dst="$_CLAUDE_HOOKS_DIR/$base"
        mkdir -p "$_CLAUDE_HOOKS_DIR"
        ownership_ensure_body \
            "$src" "$dst" \
            "project/_common/$_CLAUDE_HOOKS_DIR/$base.template" \
            "$_CLAUDE_HOOKS_STATE_DIR/$base.sha1" \
            "$_CLAUDE_HOOKS_STATE_DIR/$base.ceded" \
            exec || return 1
    done
    [[ "$found" == true ]] || return 0

    _claude_settings_ensure_hook_registration
}
