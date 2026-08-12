# lib/sciagent/claude_settings.sh — Claude Code project guardrail hooks.

# shellcheck shell=bash

_CLAUDE_PROJECT_SETTINGS=".claude/settings.json"
_CLAUDE_HOOKS_DIR=".claude/hooks"
_CLAUDE_HOOKS_STATE_DIR=".sciagent/hook_state"
_CLAUDE_HOOK_SETTINGS_STATE=".sciagent/hook_settings.state"

# _claude_settings_ensure_hook_registration
# Register the shipped hooks while preserving every unrelated project setting.
_claude_settings_ensure_hook_registration() {
    local dst="$_CLAUDE_PROJECT_SETTINGS"
    local src="$SCIAGENT_TOOLKIT/templates/project/_common/.claude/settings.json.template"
    local state="$_CLAUDE_HOOK_SETTINGS_STATE"
    [[ -f "$src" ]] || return 0

    mkdir -p "$(dirname "$dst")" "$(dirname "$state")"

    if [[ ! -f "$dst" ]]; then
        cp "$src" "$dst" || { echo "sciagent: failed to write $dst" >&2; return 1; }
        if command -v jq >/dev/null 2>&1; then
            jq -n \
                --arg h "$(sciagent_sha1_file "$dst")" \
                --slurpfile template "$src" \
                '{tag:"created", hash:$h, added:$template[0].hooks}' > "$state"
        else
            printf '{"tag":"created","hash":"%s"}\n' "$(sciagent_sha1_file "$dst")" > "$state"
        fi
        echo "wrote: $dst"
        return 0
    fi

    if ! command -v jq >/dev/null 2>&1; then
        echo "sciagent: warning — jq not found; cannot register hooks in $dst" >&2
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
        echo "sciagent: warning — could not read hook settings from $dst (invalid JSON?)" >&2
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
        echo "sciagent: warning — could not register hooks in $dst (invalid JSON?)" >&2
        return 0
    fi

    if [[ ! -f "$state" ]]; then
        cp "$dst" "${state}.orig" || { rm -f "$added" "$merged"; return 1; }
        cp "$merged" "$dst" || { rm -f "$added" "$merged"; return 1; }
        jq -n \
            --arg h "$(sciagent_sha1_file "$dst")" \
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

# _claude_settings_teardown_hook_registration
# Remove only registrations added by the toolkit and restore untouched files.
_claude_settings_teardown_hook_registration() {
    local dst="$_CLAUDE_PROJECT_SETTINGS"
    local state="$_CLAUDE_HOOK_SETTINGS_STATE"
    [[ -f "$state" ]] || return 0

    local tag stored current
    tag=$(sed -n 's/.*"tag"[[:space:]]*:[[:space:]]*"\([^"]*\)".*/\1/p' "$state" | head -1)
    stored=$(sed -n 's/.*"hash"[[:space:]]*:[[:space:]]*"\([0-9a-f]*\)".*/\1/p' "$state" | head -1)
    current=""
    [[ -f "$dst" ]] && current=$(sciagent_sha1_file "$dst")

    if [[ "$tag" == "created" && -n "$stored" && "$current" == "$stored" ]]; then
        rm -f "$dst" "$state" "${state}.orig"
        return 0
    fi

    if ! command -v jq >/dev/null 2>&1; then
        echo "sciagent: warning — jq not found; cannot safely unregister hooks from $dst" >&2
        return 0
    fi

    tag=$(jq -r '.tag // empty' "$state" 2>/dev/null)
    stored=$(jq -r 'if .tag == "created" then .hash else .written_hash end // empty' "$state" 2>/dev/null)

    if [[ "$tag" == "existed-registered" && -f "${state}.orig" \
          && -n "$stored" && "$current" == "$stored" ]]; then
        cp "${state}.orig" "$dst"
        rm -f "$state" "${state}.orig"
        return 0
    fi

    if [[ -f "$dst" ]] && jq -e '.added | type == "object"' "$state" >/dev/null 2>&1; then
        local original='{}' tmp
        if [[ -f "${state}.orig" ]]; then
            original=$(jq -c '.' "${state}.orig" 2>/dev/null || printf '{}')
        fi
        tmp=$(mktemp)
        if jq --slurpfile ownership "$state" --argjson original "$original" '
            ($ownership[0].added // {}) as $added
            | reduce ($added | to_entries[]) as $event (.;
                .hooks[$event.key] = [
                  (.hooks[$event.key] // [])[] as $entry
                  | select(($event.value | any(. == $entry)) | not)
                ]
                | if ((.hooks[$event.key] | length) == 0
                      and ((($original.hooks // {}) | has($event.key)) | not))
                  then del(.hooks[$event.key])
                  else .
                  end
              )
            | if (.hooks == {} and (($original | has("hooks")) | not))
              then del(.hooks)
              else .
              end
        ' "$dst" > "$tmp" 2>/dev/null; then
            cp "$tmp" "$dst"
        else
            echo "sciagent: warning — could not safely unregister hooks from $dst; leaving it in place" >&2
        fi
        rm -f "$tmp"
    fi

    rm -f "$state" "${state}.orig"
}

# claude_settings_ensure_hooks
# Materialize both project guardrail bodies and register them in settings.json.
claude_settings_ensure_hooks() {
    local tpl_dir="$SCIAGENT_TOOLKIT/templates/project/_common/$_CLAUDE_HOOKS_DIR"
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

# claude_settings_teardown_hooks
# Unregister and remove unchanged guardrail bodies; preserve user-owned files.
claude_settings_teardown_hooks() {
    _claude_settings_teardown_hook_registration

    if [[ -d "$_CLAUDE_HOOKS_STATE_DIR" ]]; then
        local sf base
        for sf in "$_CLAUDE_HOOKS_STATE_DIR"/*.sha1; do
            [[ -f "$sf" ]] || continue
            base="$(basename "$sf" .sha1)"
            ownership_teardown_body \
                "$_CLAUDE_HOOKS_DIR/$base" "$sf" "${sf%.sha1}.ceded"
        done
        rm -f "$_CLAUDE_HOOKS_STATE_DIR"/*.ceded 2>/dev/null || true
        rmdir "$_CLAUDE_HOOKS_STATE_DIR" 2>/dev/null || true
    fi
    [[ -d "$_CLAUDE_HOOKS_DIR" ]] && rmdir "$_CLAUDE_HOOKS_DIR" 2>/dev/null || true
    return 0
}
