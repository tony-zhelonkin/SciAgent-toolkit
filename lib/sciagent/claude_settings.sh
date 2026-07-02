# lib/sciagent/claude_settings.sh — Claude Code settings.local.json mgmt.
#
# Claude Code reads `outputStyle` from .claude/settings.local.json (and
# .claude/settings.json) to decide which output-style file under
# .claude/output-styles/ to load as its system prompt. Symlinking the
# style file into .claude/output-styles/ makes the file visible to
# Claude but does NOT make it the active style — that requires the
# `outputStyle` key to point at the symlink's basename.
#
# The toolkit therefore owns the `outputStyle` key in
# .claude/settings.local.json while a role is active. Other keys
# (permissions, enabledMcpjsonServers, …) are user-owned and untouched.
#
# Ownership tracking: a one-line state tag is written to
# .sciagent/claude_settings.state on apply, so revert knows whether to
# delete the file (we created it), strip the key (we added the key to a
# pre-existing file), or warn (we overwrote a pre-existing
# user-set value).
#
# Tag values:
#   created             — file did not exist; toolkit created it
#   existed-no-style    — file existed, no outputStyle key; toolkit added one
#   existed-with-style  — file existed with outputStyle; toolkit overwrote it
#                         (prior value is lost in v1 — overwrite warning emitted
#                          on apply so the user can record it themselves)
#
# JSON mutation: jq when available; sed/awk fallback for the canonical
# flat schema (one `"outputStyle": "value"` key, two-space indent).
# Non-canonical formats: install jq.

# shellcheck shell=bash

_CLAUDE_SETTINGS_LOCAL=".claude/settings.local.json"
_CLAUDE_SETTINGS_STATE=".sciagent/claude_settings.state"

# claude_settings_apply <output_style>
# Make `<output_style>` the active style by writing it into
# .claude/settings.local.json. Writes the tracking tag to the state
# file. Emits the tag on stdout as well so callers can echo a summary.
claude_settings_apply() {
    local name="$1"
    local f="$_CLAUDE_SETTINGS_LOCAL"
    mkdir -p "$(dirname "$f")" .sciagent

    local tag prior=""
    if [[ ! -f "$f" ]]; then
        printf '{\n  "outputStyle": "%s"\n}\n' "$name" > "$f"
        tag="created"
    else
        prior=$(claude_settings_get_output_style)
        if [[ -n "$prior" ]]; then
            _claude_settings_replace_value "$name" "$f"
            tag="existed-with-style"
        else
            _claude_settings_insert_key "$name" "$f"
            tag="existed-no-style"
        fi
    fi

    printf '%s' "$tag" > "$_CLAUDE_SETTINGS_STATE"

    if [[ "$tag" == "existed-with-style" && "$prior" != "$name" ]]; then
        echo "sciagent: warning — overwrote prior outputStyle '$prior' with '$name'" >&2
        echo "  (deactivate will not restore it; v1 limitation)" >&2
    fi

    printf '%s\n' "$tag"
}

# claude_settings_teardown
# Reverse claude_settings_apply using the saved state tag. Silent no-op
# if no state file exists.
claude_settings_teardown() {
    [[ -f "$_CLAUDE_SETTINGS_STATE" ]] || return 0
    local tag
    tag=$(cat "$_CLAUDE_SETTINGS_STATE")
    local f="$_CLAUDE_SETTINGS_LOCAL"
    case "$tag" in
        created)
            rm -f "$f"
            ;;
        existed-no-style|existed-with-style)
            [[ -f "$f" ]] && _claude_settings_remove_key "$f"
            ;;
        *)
            # Unknown tag — best-effort remove the key.
            [[ -f "$f" ]] && _claude_settings_remove_key "$f"
            ;;
    esac
    rm -f "$_CLAUDE_SETTINGS_STATE"
}

# claude_settings_get_output_style
# Print the current value of outputStyle from settings.local.json.
# Empty if file missing or key absent.
claude_settings_get_output_style() {
    local f="$_CLAUDE_SETTINGS_LOCAL"
    [[ -f "$f" ]] || return 1
    if command -v jq >/dev/null 2>&1; then
        jq -r '.outputStyle // empty' "$f" 2>/dev/null
    else
        sed -n 's/^[[:space:]]*"outputStyle"[[:space:]]*:[[:space:]]*"\([^"]*\)".*/\1/p' "$f" | head -1
    fi
}

# ----- Internal mutators -----------------------------------------------------

_claude_settings_replace_value() {
    local name="$1" f="$2"
    if command -v jq >/dev/null 2>&1; then
        local tmp; tmp=$(mktemp)
        jq --arg s "$name" '.outputStyle = $s' "$f" > "$tmp" && mv "$tmp" "$f"
    else
        local esc
        esc=$(printf '%s' "$name" | sed 's/[\&|]/\\&/g')
        sed -i 's|\("outputStyle"[[:space:]]*:[[:space:]]*\)"[^"]*"|\1"'"$esc"'"|' "$f"
    fi
}

_claude_settings_insert_key() {
    local name="$1" f="$2"
    if command -v jq >/dev/null 2>&1; then
        local tmp; tmp=$(mktemp)
        jq --arg s "$name" '. + {outputStyle: $s}' "$f" > "$tmp" && mv "$tmp" "$f"
    else
        # Insert `  "outputStyle": "name",` on the line after the first `{`.
        # Handles the empty-object case `{}` by rewriting the file.
        if grep -qE '^[[:space:]]*\{[[:space:]]*\}[[:space:]]*$' "$f"; then
            printf '{\n  "outputStyle": "%s"\n}\n' "$name" > "$f"
            return 0
        fi
        local tmp; tmp=$(mktemp)
        awk -v style="$name" '
            !done && /\{/ {
                print
                print "  \"outputStyle\": \"" style "\","
                done = 1
                next
            }
            { print }
        ' "$f" > "$tmp" && mv "$tmp" "$f"
    fi
}

# ----- Project-level .claude/settings.json + statusline.sh -----------------
#
# Root cause (see CHANGELOG): `sciagent new project` renders
# templates/project/_common/.claude/{settings.json,statusline.sh}.template
# via new.sh's generic _render_tree, but that only ever runs at project
# scaffold time. `sciagent activate` — the verb actually run against
# pre-existing / vendored-toolkit projects (the overwhelming majority in
# practice) — never touched .claude/settings.json or .claude/statusline.sh
# at all. Result: the status line + gold-standard defaults (editorMode vim,
# effortLevel xhigh, autoMemoryEnabled false, hooks, …) were materialized
# almost nowhere outside of freshly-`new project`-scaffolded repos.
#
# Fix: activate (and therefore update, which re-runs activate) now also
# ensures the status line script and settings.json exist/are current,
# using the same templates. Non-clobbering by construction:
#   - statusline.sh: written only if absent (a user's customized script is
#     left untouched); executable bit is (re-)asserted every time either way.
#   - settings.json: written verbatim if absent; if present, missing top-level
#     keys are backfilled from the template via a reverse `jq` merge that
#     always lets the EXISTING file's values win on any conflict — this only
#     ever adds keys the user's file lacks, never overwrites a user-set value.
#     Requires jq; a clear one-line warning is emitted (and materialization
#     skipped) if jq is unavailable, rather than attempting a lossy sed merge
#     of an arbitrary-shape JSON file.

_CLAUDE_PROJECT_SETTINGS=".claude/settings.json"
_CLAUDE_STATUSLINE=".claude/statusline.sh"

# claude_settings_ensure_statusline
# Materialize .claude/statusline.sh from the toolkit template if absent, and
# make sure it is executable either way. No-op (besides chmod) if the file
# already exists — a user's customized status line is never overwritten.
claude_settings_ensure_statusline() {
    local dst="$_CLAUDE_STATUSLINE"
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

# claude_settings_ensure_project_defaults
# Materialize .claude/settings.json from the toolkit template if absent.
# If present, backfill any top-level keys missing from the user's file with
# the template's gold-standard defaults (statusLine, editorMode, effortLevel,
# alwaysThinkingEnabled, autoMemoryEnabled, hooks, …) without touching any
# key the user already set.
claude_settings_ensure_project_defaults() {
    local dst="$_CLAUDE_PROJECT_SETTINGS"
    local src="$SCIAGENT_TOOLKIT/templates/project/_common/.claude/settings.json.template"
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

_claude_settings_remove_key() {
    local f="$1"
    if command -v jq >/dev/null 2>&1; then
        local tmp; tmp=$(mktemp)
        jq 'del(.outputStyle)' "$f" > "$tmp" && mv "$tmp" "$f"
        # If the result is `{}` and the file was previously non-trivial,
        # leave it — user owns the file; toolkit only owns the key.
    else
        local tmp; tmp=$(mktemp)
        # Drop the outputStyle line.
        sed '/^[[:space:]]*"outputStyle"[[:space:]]*:[[:space:]]*"[^"]*"[[:space:]]*,\?[[:space:]]*$/d' "$f" > "$tmp"
        # If the last remaining content line before the closing `}` ends
        # with a trailing comma, strip it.
        awk '
            { lines[++n] = $0 }
            END {
                for (i = 1; i <= n; i++) {
                    line = lines[i]
                    if (line ~ /,[[:space:]]*$/) {
                        j = i + 1
                        while (j <= n && lines[j] ~ /^[[:space:]]*$/) j++
                        if (j <= n && lines[j] ~ /^[[:space:]]*\}/) {
                            sub(/,[[:space:]]*$/, "", line)
                        }
                    }
                    print line
                }
            }
        ' "$tmp" > "$tmp.2" && mv "$tmp.2" "$f" && rm -f "$tmp"
    fi
}
