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
_CLAUDE_HOOKS_DIR=".claude/hooks"

# ----- Ownership records for the project-level artifacts -------------------
#
# statusline.sh, the hook bodies, and the settings.json key-backfill are
# materialized by `activate` but were never reversed by `deactivate` — the
# defect this section fixes. Each gets the same kind of ownership record the
# already-reversible artifacts have (symlinks carry their target, the ROLES/
# CRAFT blocks carry a hash-framed marker, outputStyle carries the tag above):
# a one-line/one-file SNAPSHOT taken at creation/modification time, consulted
# exactly once at teardown.
#
#   - State file ABSENT  → we made no change here; deactivate touches nothing.
#   - State file PRESENT, current content hash MATCHES the snapshot → our
#     content, untouched since → deactivate reverses it exactly.
#   - State file PRESENT, hash MISMATCH → the user edited or replaced it since
#     → deactivate leaves it alone and warns on stderr. It never re-checks:
#     the state file is removed either way, so the decision is made once and
#     the artifact is considered the user's from then on (this is what makes
#     `deactivate` idempotent — a second run finds no state file and is a
#     silent no-op, rather than re-warning forever).
#
# Directories are only ever rmdir'd (never rm -r) — a non-empty directory
# means real content remains and is left in place by construction.
_CLAUDE_STATUSLINE_STATE=".sciagent/statusline.sha1"
_CLAUDE_HOOKS_STATE_DIR=".sciagent/hook_state"
_CLAUDE_PROJECT_SETTINGS_STATE=".sciagent/project_settings.state"

# Content hashing for the ownership records above routes through block.sh's
# sciagent_sha1_file. block.sh is the single home for the hashing primitive in
# lib/ (the 5.3 invariant, enforced by tests/test_block_marker_boundary.sh,
# which greps for the tool name literally — so do not name it here either).
# Re-rolling the hash locally would duplicate the primitive and trip the guard.
_sciagent_sha1_file() {
    sciagent_sha1_file "$1"
}

# claude_settings_ensure_hooks
# Materialize .claude/hooks/*.sh from the toolkit templates, so the hooks that
# settings.json registers actually exist. Same non-clobbering contract as the
# status line: written only if absent, executable bit re-asserted every run.
# Each file we create gets a content-hash snapshot under
# $_CLAUDE_HOOKS_STATE_DIR so `deactivate` can reverse exactly the files we
# wrote (see claude_settings_teardown_hooks); a pre-existing hook (user's own,
# or ours from an earlier activate) gets no snapshot and is never touched.
#
# Why this exists: settings.json.template REGISTERS PreToolUse/Stop hooks by
# path, but the bodies were rendered only by `sciagent new project`. Any project
# retrofitted by `activate` therefore got hook registration with no hook files —
# a silently dead enforcement tier. Observed in Meta-Aging/14616-DM, which had
# both hooks registered and no .claude/hooks/ directory at all.
#
# Only *.sh.template is materialized. README.md.template carries a {{PROJECT_ID}}
# placeholder and belongs to `new project`'s substitution pass; rendering it here
# would emit a half-substituted file.
claude_settings_ensure_hooks() {
    local tpl_dir="$SCIAGENT_TOOLKIT/templates/project/_common/$_CLAUDE_HOOKS_DIR"
    [[ -d "$tpl_dir" ]] || return 0   # templates missing (e.g. stripped install): silent no-op

    local src dst found=false
    for src in "$tpl_dir"/*.sh.template; do
        [[ -f "$src" ]] || continue   # nullglob-safe: unmatched glob stays literal
        found=true
        dst="$_CLAUDE_HOOKS_DIR/$(basename "${src%.template}")"
        mkdir -p "$_CLAUDE_HOOKS_DIR"
        if [[ ! -f "$dst" ]]; then
            cp "$src" "$dst" || { echo "sciagent: failed to write $dst" >&2; return 1; }
            chmod +x "$dst"
            echo "wrote: $dst"
            mkdir -p "$_CLAUDE_HOOKS_STATE_DIR"
            _sciagent_sha1_file "$dst" > "$_CLAUDE_HOOKS_STATE_DIR/$(basename "$dst").sha1"
        else
            chmod +x "$dst" 2>/dev/null || true
        fi
    done
    [[ "$found" == true ]] || return 0
}

# claude_settings_teardown_hooks
# Reverse claude_settings_ensure_hooks: for every hook we have a snapshot for,
# remove the file iff its content is still exactly what we wrote; otherwise
# warn and leave it (see ownership-record note above). Removes the state dir
# and, if it ends up empty, .claude/hooks/ itself.
claude_settings_teardown_hooks() {
    [[ -d "$_CLAUDE_HOOKS_STATE_DIR" ]] || return 0
    local sf base dst stored cur
    for sf in "$_CLAUDE_HOOKS_STATE_DIR"/*.sha1; do
        [[ -f "$sf" ]] || continue
        base="$(basename "$sf" .sha1)"
        dst="$_CLAUDE_HOOKS_DIR/$base"
        stored=$(cat "$sf")
        if [[ -f "$dst" ]]; then
            cur=$(_sciagent_sha1_file "$dst")
            if [[ "$cur" == "$stored" ]]; then
                rm -f "$dst"
            else
                echo "sciagent: warning — $dst was modified since sciagent created it; leaving it in place" >&2
            fi
        fi
        rm -f "$sf"
    done
    rmdir "$_CLAUDE_HOOKS_STATE_DIR" 2>/dev/null || true
    [[ -d "$_CLAUDE_HOOKS_DIR" ]] && rmdir "$_CLAUDE_HOOKS_DIR" 2>/dev/null
    return 0
}

# claude_settings_ensure_statusline
# Materialize .claude/statusline.sh from the toolkit template if absent, and
# make sure it is executable either way. No-op (besides chmod) if the file
# already exists — a user's customized status line is never overwritten.
# Snapshots the content hash under $_CLAUDE_STATUSLINE_STATE only when WE
# create the file, so claude_settings_teardown_statusline knows to reverse it.
claude_settings_ensure_statusline() {
    local dst="$_CLAUDE_STATUSLINE"
    local src="$SCIAGENT_TOOLKIT/templates/project/_common/.claude/statusline.sh.template"
    [[ -f "$src" ]] || return 0   # template missing (e.g. stripped install): silent no-op

    mkdir -p "$(dirname "$dst")"
    if [[ ! -f "$dst" ]]; then
        cp "$src" "$dst" || { echo "sciagent: failed to write $dst" >&2; return 1; }
        chmod +x "$dst"
        echo "wrote: $dst"
        mkdir -p "$(dirname "$_CLAUDE_STATUSLINE_STATE")"
        _sciagent_sha1_file "$dst" > "$_CLAUDE_STATUSLINE_STATE"
    else
        chmod +x "$dst" 2>/dev/null || true
    fi
}

# claude_settings_teardown_statusline
# Reverse claude_settings_ensure_statusline using the saved content-hash
# snapshot: remove the file iff unchanged since we wrote it; otherwise warn
# and leave it (see ownership-record note above). No-op if we never created it.
claude_settings_teardown_statusline() {
    [[ -f "$_CLAUDE_STATUSLINE_STATE" ]] || return 0
    local stored cur
    stored=$(cat "$_CLAUDE_STATUSLINE_STATE")
    if [[ -f "$_CLAUDE_STATUSLINE" ]]; then
        cur=$(_sciagent_sha1_file "$_CLAUDE_STATUSLINE")
        if [[ "$cur" == "$stored" ]]; then
            rm -f "$_CLAUDE_STATUSLINE"
        else
            echo "sciagent: warning — $_CLAUDE_STATUSLINE was modified since sciagent created it; leaving it in place" >&2
        fi
    fi
    rm -f "$_CLAUDE_STATUSLINE_STATE"
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
        printf '{"tag":"created","hash":"%s"}\n' "$(_sciagent_sha1_file "$dst")" > "$state"
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
                printf '%s\n' "$state_json" > "$state"
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
                cur=$(_sciagent_sha1_file "$dst")
                if [[ "$cur" == "$stored" ]]; then
                    rm -f "$dst"
                else
                    echo "sciagent: warning — $dst was modified since sciagent created it; leaving it in place" >&2
                fi
            fi
            ;;
        existed-backfilled)
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

# claude_settings_teardown_project_artifacts
# Reverse everything claude_settings_ensure_{statusline,hooks,project_defaults}
# may have created/modified. Called only from a FULL deactivate (no stack
# remains active) — never from activate.sh's re-activation preamble, since
# these artifacts are not stack-specific and ensure_* is already idempotent
# without needing a teardown-then-recreate round trip.
claude_settings_teardown_project_artifacts() {
    claude_settings_teardown_statusline
    claude_settings_teardown_hooks
    claude_settings_teardown_project_defaults
}
