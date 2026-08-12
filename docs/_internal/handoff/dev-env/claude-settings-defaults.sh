# Extracted from SciAgent-toolkit for adoption by dev-env.
# Requires ownership_ensure_body, ownership_teardown_body, and sciagent_sha1_file.

_CLAUDE_PROJECT_SETTINGS=".claude/settings.json"
_CLAUDE_STATUSLINE=".claude/statusline.sh"
_CLAUDE_STATUSLINE_STATE=".sciagent/statusline.sha1"
_CLAUDE_PROJECT_SETTINGS_STATE=".sciagent/project_settings.state"
_CLAUDE_STATUSLINE_CEDED=".sciagent/statusline.ceded"

# claude_settings_ensure_statusline
# Materialize .claude/statusline.sh from the toolkit template and keep it
# current, via the same ownership logic as the hooks
# (ownership.sh's ownership_ensure_body): create if absent, refresh a body whose bytes
# are a version this toolkit shipped, cede anything else. A user's customized
# status line is still never overwritten — that is what the cede branch is for.
# Executable bit re-asserted every run regardless.
claude_settings_ensure_statusline() {
    local dst="$_CLAUDE_STATUSLINE"
    local src="$SCIAGENT_TOOLKIT/templates/project/_common/.claude/statusline.sh.template"
    [[ -f "$src" ]] || return 0   # template missing (e.g. stripped install): silent no-op

    ownership_ensure_body \
        "$src" "$dst" \
        "project/_common/.claude/statusline.sh.template" \
        "$_CLAUDE_STATUSLINE_STATE" \
        "$_CLAUDE_STATUSLINE_CEDED" \
        exec
}

# claude_settings_teardown_statusline
# Reverse claude_settings_ensure_statusline using the saved content-hash
# snapshot: remove the file iff unchanged since we wrote it; otherwise warn
# and leave it (see ownership-record note above). No-op if we never created it.
claude_settings_teardown_statusline() {
    ownership_teardown_body \
        "$_CLAUDE_STATUSLINE" "$_CLAUDE_STATUSLINE_STATE" "$_CLAUDE_STATUSLINE_CEDED"
}

# claude_settings_ensure_project_defaults
# Materialize .claude/settings.json from the toolkit template if absent.
# If present, backfill any top-level keys missing from the user's file with
# the template's gold-standard defaults (statusLine, editorMode, effortLevel,
# alwaysThinkingEnabled, autoMemoryEnabled, hooks, …) without touching any
# key the user already set.
#
# Records what we did to $_CLAUDE_PROJECT_SETTINGS_STATE (JSON) so
# claude_settings_teardown_project_defaults can reverse exactly that and
# nothing else:
#   {"tag":"created","hash":"<sha1 of the file as we wrote it>"}
#   {"tag":"existed-backfilled","added":{"<key>":<value-as-written>, ...}}
# "added" is the set of top-level keys the reverse-merge introduced — the
# only keys we ever add to a pre-existing file — each paired with the exact
# value we wrote, so teardown can tell "still ours" from "user has since
# edited this key" per-key, not just for the file as a whole.
claude_settings_ensure_project_defaults() {
    local dst="$_CLAUDE_PROJECT_SETTINGS"
    local src="$SCIAGENT_TOOLKIT/templates/project/_common/.claude/settings.json.template"
    local state="$_CLAUDE_PROJECT_SETTINGS_STATE"
    [[ -f "$src" ]] || return 0   # template missing: silent no-op

    mkdir -p "$(dirname "$dst")" "$(dirname "$state")"

    if [[ ! -f "$dst" ]]; then
        cp "$src" "$dst" || { echo "sciagent: failed to write $dst" >&2; return 1; }
        echo "wrote: $dst"
        printf '{"tag":"created","hash":"%s"}\n' "$(sciagent_sha1_file "$dst")" > "$state"
        return 0
    fi

    if ! command -v jq >/dev/null 2>&1; then
        echo "sciagent: warning — jq not found; cannot backfill missing keys into $dst" >&2
        echo "  (statusLine + gold defaults not merged; install jq or add manually)" >&2
        return 0
    fi

    # Reverse merge: template first, existing file second — `*` is a
    # recursive merge where the RIGHT side wins on any key present on both
    # sides, so existing user values are always preserved; only keys absent
    # from the user's file are filled in from the template.
    local old_keys
    old_keys=$(jq -r 'keys[]' "$dst" 2>/dev/null | sort)

    local tmp
    tmp=$(mktemp)
    if jq -s '.[0] * .[1]' "$src" "$dst" > "$tmp" 2>/dev/null && [[ -s "$tmp" ]]; then
        if ! cmp -s "$tmp" "$dst"; then
            local new_keys added_keys
            new_keys=$(jq -r 'keys[]' "$tmp" 2>/dev/null | sort)
            added_keys=$(comm -13 <(printf '%s\n' "$old_keys") <(printf '%s\n' "$new_keys"))
            # Snapshot the user's bytes BEFORE the merge overwrites them —
            # after `cp "$tmp" "$dst"` the original is gone (see the
            # written_hash note below for why we keep it).
            cp "$dst" "${state}.orig" 2>/dev/null || true
            cp "$tmp" "$dst"
            rm -f "$tmp"
            echo "updated: $dst (backfilled missing default keys)"
            if [[ -n "$added_keys" ]]; then
                local state_json='{"tag":"existed-backfilled","added":{}}' k v
                while IFS= read -r k; do
                    [[ -z "$k" ]] && continue
                    v=$(jq -c --arg k "$k" '.[$k]' "$dst")
                    state_json=$(printf '%s' "$state_json" | jq --arg k "$k" --argjson v "$v" '.added[$k]=$v')
                done <<< "$added_keys"
                # Also keep the user's original bytes and the hash of what we
                # wrote. The `jq -s` merge above re-serialises the whole file,
                # so a hand-written single-line settings.json comes back
                # pretty-printed: semantically identical, byte-different.
                # Deleting our keys at teardown cannot undo that reformatting.
                # With the original bytes on hand we can restore them verbatim
                # — but ONLY when the file is still exactly what we wrote, or
                # we would discard whatever the user changed since. When it has
                # moved on, teardown falls back to per-key deletion, which
                # preserves their edit and accepts the reformatting.
                state_json=$(printf '%s' "$state_json" \
                    | jq --arg h "$(sciagent_sha1_file "$dst")" '.written_hash=$h')
                printf '%s\n' "$state_json" > "$state"
            else
                # No keys added (merge changed only formatting): nothing for
                # teardown to delete, so drop the snapshot rather than leave a
                # stale .orig lying around.
                rm -f "${state}.orig"
            fi
        else
            rm -f "$tmp"
        fi
    else
        rm -f "$tmp"
        echo "sciagent: warning — could not merge defaults into $dst (invalid JSON?)" >&2
    fi
}

# claude_settings_teardown_project_defaults
# Reverse claude_settings_ensure_project_defaults using the saved state:
#   created            → remove the file iff unchanged since we wrote it.
#   existed-backfilled → remove only the keys we added, and only the ones
#                         whose current value still matches what we wrote;
#                         a key the user has since edited is left in place
#                         with a warning (per-key, not all-or-nothing).
# Requires jq to read the state back (the same dependency the backfill path
# already has) — if jq is unavailable, warn and leave both the file and the
# state record untouched rather than guess.
claude_settings_teardown_project_defaults() {
    local dst="$_CLAUDE_PROJECT_SETTINGS"
    local state="$_CLAUDE_PROJECT_SETTINGS_STATE"
    [[ -f "$state" ]] || return 0

    if ! command -v jq >/dev/null 2>&1; then
        echo "sciagent: warning — jq not found; cannot safely revert $dst; leaving as-is" >&2
        return 0
    fi

    local tag
    tag=$(jq -r '.tag // empty' "$state" 2>/dev/null)

    case "$tag" in
        created)
            local stored cur
            stored=$(jq -r '.hash // empty' "$state" 2>/dev/null)
            if [[ -f "$dst" ]]; then
                cur=$(sciagent_sha1_file "$dst")
                if [[ "$cur" == "$stored" ]]; then
                    rm -f "$dst"
                else
                    echo "sciagent: warning — $dst was modified since sciagent created it; leaving it in place" >&2
                fi
            fi
            ;;
        existed-backfilled)
            # Fast path: the file is still byte-for-byte what we wrote, so the
            # user has changed nothing since — restore their original bytes
            # verbatim. This is the only way to undo the reformatting the
            # `jq -s` merge applied (a hand-written single-line settings.json
            # comes back pretty-printed otherwise). Guarded on the hash
            # precisely because restoring blindly would discard any edit the
            # user made after activation; if it does not match we fall through
            # to per-key deletion, which keeps their edit.
            local _orig="${state}.orig" _wh _cur
            _wh=$(jq -r '.written_hash // empty' "$state" 2>/dev/null)
            if [[ -f "$_orig" && -f "$dst" && -n "$_wh" ]]; then
                _cur=$(sciagent_sha1_file "$dst")
                if [[ "$_cur" == "$_wh" ]]; then
                    cp "$_orig" "$dst"
                    rm -f "$_orig" "$state"
                    return 0
                fi
            fi
            rm -f "$_orig"
            if [[ -f "$dst" ]]; then
                local keys k stored_val cur_val tmp
                keys=$(jq -r '.added | keys[]' "$state" 2>/dev/null)
                while IFS= read -r k; do
                    [[ -z "$k" ]] && continue
                    stored_val=$(jq -c --arg k "$k" '.added[$k]' "$state")
                    cur_val=$(jq -c --arg k "$k" 'if has($k) then .[$k] else null end' "$dst")
                    if [[ "$cur_val" == "$stored_val" ]]; then
                        tmp=$(mktemp)
                        jq --arg k "$k" 'del(.[$k])' "$dst" > "$tmp" && cp "$tmp" "$dst"
                        rm -f "$tmp"
                    else
                        echo "sciagent: warning — $dst key '$k' was modified since sciagent added it; leaving it in place" >&2
                    fi
                done <<< "$keys"
            fi
            ;;
    esac
    rm -f "$state"
}

# ----- User-level ~/.claude/settings.json + statusline.sh ------------------
#
# Root cause (see CHANGELOG / plans/2026-07-02-provider-agnostic): nothing ever
# seeds the container's USER-level ${CLAUDE_CONFIG_DIR:-$HOME/.claude}/settings.json.
# claude_settings_ensure_project_defaults only ever writes the PROJECT file
# (.claude/settings.json), and even that only lands on `activate` inside a project
# cwd. A fresh container home is ephemeral, so Claude's first-run writes only bare
# `{"theme":"dark"}`-class defaults and the gold power-user settings (editorMode
# vim, effortLevel, alwaysThinkingEnabled, autoMemoryEnabled false, …) are missing
# everywhere outside an activated project dir.
#
# Fix: mirror the project functions at USER level, from a hooks-free user template
# (templates/user/.claude/settings.json.template) whose statusLine.command is the
# user-level `~/.claude/statusline.sh` path (NOT $CLAUDE_PROJECT_DIR, which is
# undefined at user scope). Hooks stay a project-activate concern. Same
# non-clobbering guarantees as the project versions:
#   - settings.json: written verbatim if absent; if present, missing top-level keys
#     are backfilled via the same reverse `jq` merge that always lets the EXISTING
#     file's values win — never overwrites a user-set value. jq missing → warn +
#     skip (no lossy fallback).
#   - statusline.sh: written only if absent; executable bit (re-)asserted either way.

# claude_settings_ensure_user_statusline
# Materialize ${CLAUDE_CONFIG_DIR:-$HOME/.claude}/statusline.sh from the toolkit
# template if absent, and make sure it is executable either way. No-op (besides
# chmod) if the file already exists — a user's customized status line is never
# overwritten.
claude_settings_ensure_user_statusline() {
    local dst="${CLAUDE_CONFIG_DIR:-$HOME/.claude}/statusline.sh"
    local src="$SCIAGENT_TOOLKIT/templates/project/_common/.claude/statusline.sh.template"
    [[ -f "$src" ]] || return 0   # template missing (e.g. stripped install): silent no-op

    mkdir -p "$(dirname "$dst")"
    if [[ ! -f "$dst" ]]; then
        cp "$src" "$dst" || { echo "sciagent: failed to write $dst" >&2; return 1; }
        chmod +x "$dst"
        echo "wrote: $dst"
    else
        chmod +x "$dst" 2>/dev/null || true
    fi
}

# claude_settings_ensure_user_defaults
# Materialize ${CLAUDE_CONFIG_DIR:-$HOME/.claude}/settings.json from the hooks-free
# user template if absent. If present, backfill any top-level keys missing from the
# user's file with the template's gold-standard defaults (statusLine, editorMode,
# effortLevel, alwaysThinkingEnabled, autoMemoryEnabled, …) without touching any
# key the user already set.
claude_settings_ensure_user_defaults() {
    local dst="${CLAUDE_CONFIG_DIR:-$HOME/.claude}/settings.json"
    local src="$SCIAGENT_TOOLKIT/templates/user/.claude/settings.json.template"
    [[ -f "$src" ]] || return 0   # template missing: silent no-op

    mkdir -p "$(dirname "$dst")"
    if [[ ! -f "$dst" ]]; then
        cp "$src" "$dst" || { echo "sciagent: failed to write $dst" >&2; return 1; }
        echo "wrote: $dst"
        return 0
    fi

    if ! command -v jq >/dev/null 2>&1; then
        echo "sciagent: warning — jq not found; cannot backfill missing keys into $dst" >&2
        echo "  (statusLine + gold defaults not merged; install jq or add manually)" >&2
        return 0
    fi

    # Reverse merge: template first, existing file second — `*` is a
    # recursive merge where the RIGHT side wins on any key present on both
    # sides, so existing user values are always preserved; only keys absent
    # from the user's file are filled in from the template.
    local tmp
    tmp=$(mktemp)
    if jq -s '.[0] * .[1]' "$src" "$dst" > "$tmp" 2>/dev/null && [[ -s "$tmp" ]]; then
        if ! cmp -s "$tmp" "$dst"; then
            mv "$tmp" "$dst"
            echo "updated: $dst (backfilled missing default keys)"
        else
            rm -f "$tmp"
        fi
    else
        rm -f "$tmp"
        echo "sciagent: warning — could not merge defaults into $dst (invalid JSON?)" >&2
    fi
}
