# lib/sciagent/symlinks.sh — dual-track symlink + manifest helpers, plus the
# 02_analysis/helpers seam (the hyphenated contract-lib mounts AND the
# underscored shim modules that make them importable — two halves of one
# design, kept together; see symlink_create_helper_lib / helper_shims_ensure).
#
# Manifest format: real JSON (pretty-printed, 2-space indent). Schema v2:
#
#   {
#     "version": 2,
#     "stack": ["base", "reviewer"],
#     "symlinks": ["/abs/or/relative/path"]
#   }
#
# What the manifest is FOR, now that teardown no longer depends on it:
#   - "stack" is the only machine-readable record of WHICH roles are active,
#     and its presence is the is-a-stack-active flag (manifest_exists).
#   - "symlinks" lets `status` notice a mount that has been DELETED. A
#     target-ownership scan cannot: it sees what is present, not what is
#     missing. Teardown reads it as a second source alongside link targets so
#     that a teardown run against a different toolkit checkout still cleans up.
#
# v1 -> v2 (2026-08-11): dropped "block_hash". activate wrote it; nothing ever
# read it back — the drift guard reads the hash from the AGENTS.md BEGIN marker
# via block_hash_check (status.sh, craft_verb.sh), never from here. Dead data in
# a state file is worse than no data: the next reader assumes it is current,
# because nothing about a stale value looks stale.
#
# Compatibility, both directions, with no migration step:
#   - An OLD on-disk manifest (v1, still carrying block_hash) reads fine. Both
#     readers are key-targeted — manifest_stack takes `.stack`, manifest_symlinks
#     takes `.symlinks[]` — and neither cares about extra keys. That is the
#     normal case across the fleet until each project next activates, which
#     rewrites the manifest wholesale (manifest_begin -> manifest_finalize).
#   - A NEW manifest read by an OLDER toolkit also reads fine, since the key it
#     no longer finds is one it never read.
#
# "version" is bumped rather than left at 1 even though NOTHING BRANCHES ON IT
# today (no reader in lib/, none in tests beyond an existence assertion). The
# reason is narrow and not aesthetic: `version: 1` is now the only way to tell
# an on-disk manifest that MAY carry block_hash from one that cannot. Leaving
# both shapes labelled 1 destroys that distinction permanently, and it is only
# recordable at the moment of the change. Any future reader must therefore
# treat >= 2 as "no block_hash" and 1 as "may carry an ignorable one" — it must
# not treat an unknown version as fatal.
#
# An "injected" array was removed in 2026-08 along with the readers that
# round-tripped it; it was residue of the retired `inject`/`eject` verbs and
# was always emitted empty.
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
#   agents    — finds <name>.md under $SCIAGENT_TOOLKIT/agents/
#   commands  — finds <name>.md under $SCIAGENT_TOOLKIT/commands/
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

manifest_begin() {
    local stack="$1"   # space-separated role names
    _manifest_staging=$(mktemp)
    _manifest_stack_val="$stack"
    _manifest_syms=()
}

_manifest_record_symlink() {
    [[ -n "$_manifest_staging" ]] || return 0
    _manifest_syms+=("$1")
}

# _manifest_write_json
# Serialise accumulated state into $_manifest_staging as pretty JSON.
_manifest_write_json() {
    {
        printf '{\n'
        printf '  "version": 2,\n'

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
                (( first )) && first=0 || printf ',\n'
                printf '"%s"' "$(_json_escape "$s")"
            done
        fi
        printf ']\n'

        printf '}\n'
    } > "$_manifest_staging"
}

# manifest_finalize
# Takes no arguments: the block hash it used to receive was written to a field
# nothing read (see the schema note in this file's header).
manifest_finalize() {
    [[ -n "$_manifest_staging" ]] || return 1
    _manifest_write_json
    mkdir -p "$(dirname "$_MANIFEST_PATH")"
    mv "$_manifest_staging" "$_MANIFEST_PATH"
    _manifest_staging=""
}

# ---------------------------------------------------------------------------
# symlink_target_for <link_path> <canonical>
# Compute the target string to store in a symlink at <link_path> pointing at
# <canonical> (an absolute path into the toolkit).
#
# Portability rule (mirrors symlink_create_helper_lib): when the toolkit lives
# INSIDE the project tree (the activation CWD — where .claude/, .agents/,
# AGENTS.md, .sciagent/ are written), emit a RELATIVE target computed relative
# to the link's OWN directory. A relative link then resolves identically on any
# host/container and travels with the committed submodule pin, instead of
# hard-coding whatever absolute path $SCIAGENT_TOOLKIT happened to be at
# activation time.
#
# When the toolkit is NOT within the project tree (an external/global checkout,
# e.g. under /data1 or a container /workspaces path), fall back to the ABSOLUTE
# path. A relative link in that case would be a long, fragile
# `../../../../data1/...` chain that breaks the moment the project or the
# external checkout moves — the absolute path is the more robust choice there.
#
# Prints the target string on stdout.
symlink_target_for() {
    local link_path="$1"
    local canonical="$2"
    local link_dir
    link_dir="$(dirname "$link_path")"

    # Project root is the activation CWD. The toolkit is "in-repo" when the
    # canonical target path lies at or below it. We test the canonical path AS
    # GIVEN (via $SCIAGENT_TOOLKIT) rather than its realpath-resolved form:
    # this mirrors symlink_create_helper_lib and, crucially, keeps links
    # relative even when the in-repo toolkit is itself reached through a symlink
    # (e.g. a submodule surfaced under 01_modules/). Relativizing against a
    # symlink-resolved path could instead point outside the project and defeat
    # portability. The relative target is computed with `realpath -ms
    # --relative-to`: -s (no-symlinks) keeps it a purely lexical path
    # computation so a symlink-surfaced toolkit (e.g. a submodule reached via a
    # symlink under 01_modules/) still yields an in-project `../../…` target
    # rather than escaping to the symlink's resolved location; -m tolerates
    # not-yet-existing path components.
    local proj_root
    proj_root="$(pwd)"

    if [[ "$canonical" == "$proj_root"/* ]] \
        && command -v realpath >/dev/null 2>&1 \
        && realpath -ms --relative-to="$link_dir" "$canonical" >/dev/null 2>&1; then
        realpath -ms --relative-to="$link_dir" "$canonical"
    else
        printf '%s' "$canonical"
    fi
}

# ---------------------------------------------------------------------------
# symlink_create_dual
# ---------------------------------------------------------------------------

# symlink_create_dual <category> <name> <canonical_path>
# Categories: skills (dir symlink), agents/commands (file symlink).
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
        *)
            echo "symlink_create_dual: unknown category '$category'" >&2
            return 1
            ;;
    esac
    # Adding a category here REQUIRES adding its directories to
    # _SCIAGENT_MOUNT_DIRS. Since teardown derives ownership from link targets,
    # that list is the only record of where to look — a mount in an unlisted
    # directory is never torn down and never reported, silently, forever.

    mkdir -p "$(dirname "$claude_path")"
    ln -sfn "$(symlink_target_for "$claude_path" "$canonical")" "$claude_path"
    _manifest_record_symlink "$claude_path"

    if [[ -n "$agents_path" ]]; then
        mkdir -p "$(dirname "$agents_path")"
        ln -sfn "$(symlink_target_for "$agents_path" "$canonical")" "$agents_path"
        _manifest_record_symlink "$agents_path"
    fi
}

# ---------------------------------------------------------------------------
# _sciagent_toolkit_locality_ok — pure predicate, no output, no mutation.
# True (rc 0) iff EITHER the project ships no in-repo toolkit at all (nothing
# to protect — see _sciagent_in_repo_toolkit below for how that is decided),
# OR the currently active $SCIAGENT_TOOLKIT resolves to it. False (rc 1) iff an
# in-repo toolkit exists and the active one is a different (external)
# checkout.
#
# This is the boolean core of bin/sciagent's `_guard_toolkit_locality` —
# factored out here (rather than duplicated) so that lib code which re-enters
# a mutating path IN-PROCESS, bypassing the dispatcher's own guard call, can
# still check locality before it mounts anything. `update` calling
# `cmd_activate` in-process (update.sh) is exactly this shape, and so is
# `deactivate <overlay>` re-activating solo-base (deactivate.sh) — the
# dispatcher guard only runs once, for the TOP-LEVEL verb, before any lib code
# is even sourced, so a verb that internally calls into another mutating verb
# needs its own check. bin/sciagent sources this file unconditionally and
# early (before its own guard call) so both call sites share one definition.
# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------
# _sciagent_in_repo_toolkit — print the project's own toolkit path, or nothing.
#
# Returns 0 and echoes a path iff this project ships an in-repo toolkit;
# returns 1 and echoes nothing otherwise.
#
# WHY THIS IS NOT JUST `./01_modules/SciAgent-toolkit`: that literal was
# hardcoded in three places and matched CASE-SENSITIVELY, while the real fleet
# uses four different container directories — measured 2026-08-11 across 23
# checkouts: 17 `01_modules`, 2 `01_Scripts`, 2 `01_scripts`, 1 `01_Modules`.
# In the six non-conforming projects the `-d` test failed, so the locality
# guard concluded "no in-repo toolkit, nothing to protect" and permitted
# mutation from ANY external checkout — precisely the escape it exists to
# prevent, and silently. Two of those projects live in a SHARED lab tree
# (/data2/users/JCRLab), which is exactly where a stray $SCIAGENT_TOOLKIT is
# most likely to come from.
#
# Resolution order, most authoritative first:
#   1. .gitmodules — git's own record of where the submodule lives. Verified
#      correct in all four non-conforming consumers.
#   2. the conventional path, for a project with no .gitmodules (vendored copy,
#      not-yet-a-repo scaffold).
#   3. any single-level container holding a SciAgent-toolkit directory.
# ---------------------------------------------------------------------------
_sciagent_in_repo_toolkit() {
    local p
    if [[ -f .gitmodules ]] && command -v git >/dev/null 2>&1; then
        p=$(git config -f .gitmodules --get-regexp '^submodule\..*\.path$' 2>/dev/null \
            | awk '{print $2}' | grep -E '(^|/)SciAgent-toolkit$' | head -1)
        if [[ -n "$p" && -d "$p" ]]; then
            printf '%s\n' "./$p"
            return 0
        fi
    fi
    if [[ -d "./01_modules/SciAgent-toolkit" ]]; then
        printf '%s\n' "./01_modules/SciAgent-toolkit"
        return 0
    fi
    # Unmatched globs stay literal, and `-d` is false for them, so this is safe
    # in a project with no container directory at all.
    for p in ./*/SciAgent-toolkit; do
        if [[ -d "$p" ]]; then
            printf '%s\n' "$p"
            return 0
        fi
    done
    return 1
}

_sciagent_toolkit_locality_ok() {
    local in_repo
    in_repo=$(_sciagent_in_repo_toolkit) || return 0

    local in_repo_real active_real
    if command -v realpath >/dev/null 2>&1; then
        in_repo_real="$(realpath "$in_repo" 2>/dev/null || printf '%s' "$in_repo")"
        active_real="$(realpath "${SCIAGENT_TOOLKIT:-}" 2>/dev/null || printf '%s' "${SCIAGENT_TOOLKIT:-}")"
    else
        in_repo_real="$(cd "$in_repo" 2>/dev/null && pwd -P || printf '%s' "$in_repo")"
        active_real="$(cd "${SCIAGENT_TOOLKIT:-}" 2>/dev/null && pwd -P || printf '%s' "${SCIAGENT_TOOLKIT:-}")"
    fi

    [[ "$active_real" == "$in_repo_real" ]]
}

# ---------------------------------------------------------------------------
# Toolkit-ownership check (Phase 5c follow-up, "Option D").
#
# The manifest's "symlinks" array used to be the sole record of what
# activate created, and teardown trusted it blindly. That has two failure
# modes: (a) anything the manifest lost track of (crashed activate,
# hand-edited/deleted manifest.json, a manual `ln -sfn` repair) is stranded
# forever, never cleaned up; (b) it offers no protection beyond "not in this
# list" — it doesn't independently verify a path is actually ours before
# deleting it.
#
# Ownership is now a property of the LINK ITSELF: a mount is toolkit-owned
# iff it is a symlink whose fully-resolved target lies inside the
# (canonicalised) toolkit checkout. Everything else — a real user file/dir,
# or a user's own symlink pointing OUTSIDE the toolkit — is left strictly
# alone, independent of whatever the manifest does or doesn't know.
# ---------------------------------------------------------------------------

# _toolkit_canonical_root — fully-resolved (all symlink components followed)
# absolute path to $SCIAGENT_TOOLKIT, or empty + rc=1 if it can't be
# resolved (env unset, `realpath` unavailable, or the toolkit checkout is
# itself gone). Callers MUST treat that as "ownership can't be determined,
# do nothing" — never fall back to a bare string comparison, which is
# exactly how "/toolkit-evil" would wrongly match a "/toolkit" prefix.
_toolkit_canonical_root() {
    [[ -n "${SCIAGENT_TOOLKIT:-}" ]] || return 1
    command -v realpath >/dev/null 2>&1 || return 1
    realpath -e "$SCIAGENT_TOOLKIT" 2>/dev/null
}

# _link_owned_by_toolkit <path> <toolkit_canonical_root>
# True iff <path> is a symlink whose normalized target is the toolkit root
# itself or a proper descendant of it. Path-boundary comparison (trailing
# "/" on the toolkit side), not a string prefix — "/toolkit-evil" must not
# match "/toolkit". `realpath -m` preserves ownership detection when multiple
# target components are gone, including retired system-prompts/ paths.
#
# Deliberate, accepted behavior (not an oversight — owner sign-off): a user
# symlink that happens to point INTO the toolkit checkout is indistinguishable
# from a toolkit-created mount by any means available here, and is therefore
# treated as toolkit-owned and removed on deactivate, same as any other
# toolkit-owned link.
_link_owned_by_toolkit() {
    local path="$1" toolkit_root="$2"
    [[ -L "$path" ]] || return 1
    [[ -n "$toolkit_root" ]] || return 1
    local target resolved
    target=$(readlink "$path" 2>/dev/null) || return 1
    if [[ "$target" == /* ]]; then
        resolved=$(realpath -m "$target" 2>/dev/null) || return 1
    else
        resolved=$(realpath -m "$(dirname "$path")/$target" 2>/dev/null) || return 1
    fi
    [[ -n "$resolved" ]] || return 1
    case "$resolved" in
        "$toolkit_root")   return 0 ;;
        "$toolkit_root"/*) return 0 ;;
        *)                  return 1 ;;
    esac
}

# Every directory sciagent mounts into, plus the retired output-style location
# retained as a teardown sweep for legacy toolkit-owned symlinks.
# maxdepth 1 when scanning these — resolve_canonical hides toolkit-side
# subfolders, so a mount is always a flat entry directly inside one of these;
# we must never recurse into a mounted skill DIRECTORY itself (that would
# walk into the toolkit checkout through the symlink).
_SCIAGENT_MOUNT_DIRS=(
    .claude/skills .claude/agents .claude/commands .claude/output-styles
    .agents/skills .agents/agents .agents/commands
    02_analysis/helpers
)

# ---------------------------------------------------------------------------
# symlink_teardown_all — remove exactly the symlinks sciagent mounted, via the
# UNION of two independent sources:
#
#   1. Target-based ownership: sweep every _SCIAGENT_MOUNT_DIRS entry and
#      remove any symlink whose fully-resolved target lies inside the
#      CURRENTLY ACTIVE $SCIAGENT_TOOLKIT (see block comment above). This is
#      the orphan-recovery path — it needs no manifest, so a lost/hand-edited/
#      deleted manifest.json never strands a mount.
#   2. Manifest-recorded paths: every path in a still-present manifest's
#      "symlinks" array, if it is still a symlink. Trusting these unconditionally
#      (regardless of what they currently resolve to) is safe because they were
#      never guessed — every entry was written by OUR OWN
#      symlink_create_dual/symlink_create_helper_lib at mount time, first-party
#      bookkeeping, not an inference from a target path.
#
# Source 2 exists because source 1 alone regresses when $SCIAGENT_TOOLKIT
# differs between mount time and teardown time — e.g. a project mounted
# against its in-repo toolkit, then `deactivate` invoked with
# SCIAGENT_TOOLKIT=<a different external checkout>. None of the mounts resolve
# under that external root, so a target-only sweep silently finds nothing,
# while the caller still deletes the manifest and the managed blocks — exactly
# the "intermediate state worse than doing nothing or doing everything" this
# union avoids. Source 1 alone remains necessary for the complementary case
# (manifest missing/stale) that source 2 alone can't cover. `rm` unlinks the
# symlink itself in every case below — never `rm -r`/`rm -rf`, and never
# anything that would dereference into the toolkit checkout the link points
# at.
#
# Sets _SCIAGENT_TEARDOWN_COUNT to the number of links removed, so a caller can
# report what actually happened rather than guessing from manifest presence.
# ---------------------------------------------------------------------------
symlink_teardown_all() {
    _SCIAGENT_TEARDOWN_COUNT=0
    # No manifest gate on purpose. Ownership comes from each link's target, so
    # a missing manifest.json is not a reason to leave mounts behind — it is
    # the exact orphan case this function exists to recover. Gating here made
    # `rm -rf .sciagent` strand every mount permanently: the sweep below never
    # ran, and nothing else ever removes them.
    local toolkit_root=""
    if ! toolkit_root=$(_toolkit_canonical_root); then
        echo "warning: could not resolve \$SCIAGENT_TOOLKIT ('${SCIAGENT_TOOLKIT:-<unset>}'); skipping toolkit-owned symlink cleanup this run" >&2
        toolkit_root=""
    fi

    local d entry
    for d in "${_SCIAGENT_MOUNT_DIRS[@]}"; do
        [[ -d "$d" ]] || continue
        while IFS= read -r entry; do
            [[ -z "$entry" ]] && continue
            if [[ -n "$toolkit_root" ]] && _link_owned_by_toolkit "$entry" "$toolkit_root"; then
                rm "$entry"
                _SCIAGENT_TEARDOWN_COUNT=$((_SCIAGENT_TEARDOWN_COUNT + 1))
            fi
        done < <(find "$d" -mindepth 1 -maxdepth 1 -type l 2>/dev/null)
    done

    # Source 2: manifest-recorded paths (see header comment above). Runs AFTER
    # the target-based sweep and BEFORE the manifest is deleted below. A path
    # already removed by source 1 simply fails the `-L` test here and is
    # skipped — no double-count, no error.
    if manifest_exists; then
        local mpath
        while IFS= read -r mpath; do
            [[ -n "$mpath" ]] || continue
            [[ -L "$mpath" ]] || continue
            rm "$mpath"
            _SCIAGENT_TEARDOWN_COUNT=$((_SCIAGENT_TEARDOWN_COUNT + 1))
        done < <(manifest_symlinks 2>/dev/null)
    fi

    rm -f "$_MANIFEST_PATH"
    rmdir .sciagent 2>/dev/null || true
    # Best-effort: clean empty dirs we may have mkdir -p'd at mount time. A
    # non-empty dir here means real content remains (ours or the user's) —
    # rmdir fails on a non-empty dir and we leave it, by construction.
    for d in "${_SCIAGENT_MOUNT_DIRS[@]}"; do
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

# Removed 2026-08: manifest_injected / manifest_block_hash /
# manifest_update_block_hash (~110 lines).
#
# All three were readers with no callers. manifest_injected served the retired
# `inject`/`eject` verbs (Phase 2); its only remaining caller was
# manifest_update_block_hash, which had no callers of its own, and
# manifest_block_hash had none either — the drift guard reads the hash from the
# AGENTS.md marker (block_hash_check), never from the manifest. A value that is
# only ever read back in order to be written out unchanged is dead, so the
# whole cluster went together with the "injected" JSON key it round-tripped.
#
# The `block_hash` FIELD those functions served outlived them by one release,
# because dropping it changed the on-disk schema. That change has now been made
# (schema v2, 2026-08-11): see this file's header for the v1->v2 note and the
# both-directions compatibility argument.

# ---------------------------------------------------------------------------
# symlink_create_helper_lib
# If the project has a 02_analysis/ directory (analysis-type repo), create one
# DIRECTORY symlink per shared helper-lib into 02_analysis/helpers/:
#   02_analysis/helpers/figure-style      →  <SCIAGENT_TOOLKIT>/lib/figure-style
#   02_analysis/helpers/interactive-style →  <SCIAGENT_TOOLKIT>/lib/interactive-style
# Each symlink is RELATIVE (portable across container-mount changes) when
# `realpath --relative-to` is available; otherwise falls back to the absolute
# path. Records each symlink path in the manifest so symlink_teardown_all
# removes it on deactivate. Idempotent (ln -sfn). A lib dir that does not exist
# in the toolkit is skipped. Silently skips entirely if 02_analysis/ is absent.
symlink_create_helper_lib() {
    [[ -d "02_analysis" ]] || return 0   # not an analysis-type repo — skip

    local link_dir="02_analysis/helpers"
    mkdir -p "$link_dir"

    local libdir
    for libdir in figure-style interactive-style; do
        local target_dir="$SCIAGENT_TOOLKIT/lib/$libdir"
        [[ -d "$target_dir" ]] || continue   # toolkit lacks this lib — skip it
        local link_path="$link_dir/$libdir"

        # Compute a relative symlink target for portability.
        local rel_target
        if realpath --relative-to="$link_dir" "$target_dir" >/dev/null 2>&1; then
            rel_target=$(realpath --relative-to="$link_dir" "$target_dir")
        else
            rel_target="$target_dir"   # absolute fallback
        fi

        ln -sfn "$rel_target" "$link_path"
        _manifest_record_symlink "$link_path"
    done
}


# ---------------------------------------------------------------------------
# The other half of the 02_analysis/helpers seam: the SHIM MODULES.
#
# symlink_create_helper_lib (above) mounts HYPHENATED contract-lib directories,
# which no Python `import` can name. The house design pairs each mount with an
# UNDERSCORED shim module next to it — 02_analysis/helpers/figure_style.{py,R},
# interactive_style.py — that projects actually import (see
# tests/test_helper_shim_coherence.sh, which enforces the pairing).
#
# THE DEFECT THIS FIXES. Those shims were materialized ONLY by
# `sciagent new project`'s generic template walk (new.sh's _render_tree), so an
# already-provisioned project never received one, and `activate` created the two
# mounts with no importable module beside them. Figure-style,
# interactive-breakpoint-explorer, and decision-gate-notebook were therefore
# satisfiable only in a freshly scaffolded repo. This is the same blind spot as
# the hook bodies and the status line: content written at scaffold time only,
# never reaching the field.
#
# Ownership is the shared discipline in ownership.sh (create if absent, adopt a
# copy that already matches, refresh a copy matching any version we ever
# shipped, cede anything else). Content-provenance applies because these
# templates carry no {{PLACEHOLDER}}, so `new project`'s sed pass emits bytes
# identical to a plain `cp` — enforced by tests/test_template_provenance.sh,
# which fails a managed template that grows one.
#
# Mode `plain`: unlike the hooks and the status line these are imported, never
# executed, so the execute bit is deliberately NOT asserted.
# ---------------------------------------------------------------------------
# Path RELATIVE TO templates/ — the same frame templates/PROVENANCE.sha1 records
# its paths in, so the value handed to ownership_ensure_body can be looked up in
# the manifest verbatim. Prefixing it with templates/ here would silently make
# every provenance lookup miss (and every stale shim look user-authored).
_HELPER_SHIM_TPL_REL="project/analysis/02_analysis/helpers"
_HELPER_SHIM_STATE_DIR=".sciagent/helper_shim_state"

# helper_shims_ensure
# Materialize (and keep current) every helper shim under 02_analysis/helpers/.
# Gated on 02_analysis/ exactly like symlink_create_helper_lib: a coordination
# or software repo has no analysis layout and must be left untouched, including
# not having a helpers/ directory created for it. Silent no-op when the toolkit
# ships no template tree (e.g. a stripped install or the test fixture toolkit).
helper_shims_ensure() {
    [[ -d "02_analysis" ]] || return 0   # not an analysis-type repo — skip

    local tpl_dir="$SCIAGENT_TOOLKIT/templates/$_HELPER_SHIM_TPL_REL"
    [[ -d "$tpl_dir" ]] || return 0

    local src base
    for src in "$tpl_dir"/*.template; do
        [[ -f "$src" ]] || continue   # nullglob-safe: unmatched glob stays literal
        base="$(basename "${src%.template}")"
        ownership_ensure_body \
            "$src" "02_analysis/helpers/$base" \
            "$_HELPER_SHIM_TPL_REL/$base.template" \
            "$_HELPER_SHIM_STATE_DIR/$base.sha1" \
            "$_HELPER_SHIM_STATE_DIR/$base.ceded" \
            plain || return 1
    done
}

# helper_shims_teardown
# Reverse helper_shims_ensure: for every shim we have a snapshot for, remove it
# iff its content is still exactly what we wrote; otherwise warn and leave it
# (ownership_teardown_body). Then the state dir, and 02_analysis/helpers/ itself
# if it ends up empty — rmdir only, never rm -r, so a helpers/ directory holding
# anything else (a project-owned helper, a surviving mount) is left in place.
# 02_analysis/ itself is the project's and is never touched.
helper_shims_teardown() {
    [[ -d "$_HELPER_SHIM_STATE_DIR" ]] || return 0
    local sf base
    for sf in "$_HELPER_SHIM_STATE_DIR"/*.sha1; do
        [[ -f "$sf" ]] || continue
        base="$(basename "$sf" .sha1)"
        ownership_teardown_body \
            "02_analysis/helpers/$base" "$sf" "${sf%.sha1}.ceded"
    done
    # Ceded markers name files that are explicitly NOT ours: drop the marker,
    # never the file. Also unblocks the rmdir below.
    rm -f "$_HELPER_SHIM_STATE_DIR"/*.ceded 2>/dev/null || true
    rmdir "$_HELPER_SHIM_STATE_DIR" 2>/dev/null || true
    [[ -d "02_analysis/helpers" ]] && rmdir "02_analysis/helpers" 2>/dev/null
    return 0
}
