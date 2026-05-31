#!/usr/bin/env bash
# tests/test_no_duplicate_basenames.sh
# Two invariants in one test:
#   (a) intra-namespace uniqueness — no two .md files under commands/ or
#       agents/ (recursively) share a basename, and no two skill directories
#       share a name. The recursive resolver in lib/sciagent/symlinks.sh
#       (resolve_canonical) depends on this; duplicate basenames make role
#       bindings ambiguous.
#   (b) cross-namespace non-collision (with allowlist) — a basename appearing
#       in >=2 of skills/agents/commands/roles is the runtime soft-warn
#       surface of `sciagent validate`. The merge-time mirror is here:
#       hard-fail unless the overlap is documented in
#       tests/collision-allowlist.txt with a rationale.

set -u
# shellcheck source=/dev/null
. "$(dirname "$0")/_lib.sh"

# Allow tests to override the toolkit + allowlist via env vars; default to
# the live tree. The collision-allowlist override is what
# test_collision_allowlist_blocks_unknown.sh / _passes_known.sh exploit
# without having to plant a fixture file inside the real allowlist.
TOOLKIT="${COLLISION_TEST_TOOLKIT:-$TOOLKIT_ROOT}"
ALLOWLIST="${COLLISION_TEST_ALLOWLIST:-$TOOLKIT_ROOT/tests/collision-allowlist.txt}"

check_unique_basenames() {
    local subtree="$1"
    [[ -d "$TOOLKIT/$subtree" ]] || return 0
    local dupes
    dupes=$(
        find "$TOOLKIT/$subtree" -type f -name '*.md' \
            ! -name 'README.md' \
            -printf '%f\t%p\n' \
        | sort \
        | awk -F'\t' '
            { count[$1]++; paths[$1] = paths[$1] " " $2 }
            END {
                for (n in count) if (count[n] > 1)
                    printf "%s:%s\n", n, paths[n]
            }
          '
    )
    if [[ -n "$dupes" ]]; then
        echo "FAIL [$_TEST_NAME] duplicate basenames under $subtree/:" >&2
        echo "$dupes" >&2
        exit 1
    fi
}

check_unique_skill_dirs() {
    [[ -d "$TOOLKIT/skills" ]] || return 0
    # Skills live at skills/<name>/SKILL.md — directory names must be unique
    # within skills/. Since we only look one level deep, this is implicit on
    # POSIX filesystems; the check exists to flag accidental nested skills/
    # subfolders.
    local nested
    nested=$(find "$TOOLKIT/skills" -mindepth 2 -type d -name '*' \
        | grep -v -E '/(_TEMPLATE|[^/]+/[^/]+)$' || true)
    # Simpler: assert no SKILL.md exists deeper than skills/<name>/SKILL.md.
    local deep
    # Exclude hidden directories (e.g. .deprecated/) — they are off the
    # resolver's path and exist for historical reference only.
    deep=$(find "$TOOLKIT/skills" -mindepth 3 -name 'SKILL.md' \
        -not -path '*/.*' 2>/dev/null || true)
    if [[ -n "$deep" ]]; then
        echo "FAIL [$_TEST_NAME] SKILL.md found deeper than skills/<name>/SKILL.md:" >&2
        echo "$deep" >&2
        exit 1
    fi
}

check_unique_basenames commands
check_unique_basenames agents
check_unique_skill_dirs

# ---------------------------------------------------------------------------
# Cross-namespace collision check, gated by tests/collision-allowlist.txt.
# Reuses lib/sciagent/collisions.sh so the CI test and the runtime check
# (validate.sh check 5) cannot drift.
# ---------------------------------------------------------------------------

check_cross_namespace_collisions() {
    local helper="$TOOLKIT_ROOT/lib/sciagent/collisions.sh"
    if [[ ! -f "$helper" ]]; then
        echo "FAIL [$_TEST_NAME] missing helper: $helper" >&2
        exit 1
    fi
    # shellcheck source=/dev/null
    . "$helper"
    # The helper reads SCIAGENT_TOOLKIT at call time (not source time), so
    # we point it at the target toolkit for the duration of this function.
    # Saving + restoring keeps the rest of the test untouched.
    local _saved_tk="${SCIAGENT_TOOLKIT:-}"
    export SCIAGENT_TOOLKIT="$TOOLKIT"

    # Build allowed set: "name<TAB>kind1,kind2,..." lines from the allowlist,
    # comments and blanks dropped, first two whitespace-separated fields kept.
    declare -A allowed=()
    if [[ -f "$ALLOWLIST" ]]; then
        local line name kinds
        while IFS= read -r line || [[ -n "$line" ]]; do
            # Strip trailing comment, then leading/trailing whitespace.
            line="${line%%#*}"
            line="${line#"${line%%[![:space:]]*}"}"
            line="${line%"${line##*[![:space:]]}"}"
            [[ -z "$line" ]] && continue
            name="${line%%[[:space:]]*}"
            kinds="${line#"$name"}"
            kinds="${kinds#"${kinds%%[![:space:]]*}"}"
            [[ -z "$name" || -z "$kinds" ]] && continue
            allowed["$name"]="$kinds"
        done < "$ALLOWLIST"
    fi

    local fail=0
    local col_name col_kinds
    while IFS=$'\t' read -r col_name col_kinds; do
        [[ -z "$col_name" ]] && continue
        if [[ -z "${allowed[$col_name]+_}" ]]; then
            local k1="${col_kinds%%,*}" k2="${col_kinds#*,}"
            k2="${k2%%,*}"
            echo "error: name '$col_name' collides across $k1 and $k2. either rename one, or add to tests/collision-allowlist.txt with rationale." >&2
            fail=1
            continue
        fi
        # Allowlisted: assert the kinds csv matches exactly. A drift here
        # means a new namespace acquired a same-named entry — must come
        # through PR review, not silently piggyback on the existing entry.
        if [[ "${allowed[$col_name]}" != "$col_kinds" ]]; then
            echo "error: name '$col_name' collides across $col_kinds but allowlist records '${allowed[$col_name]}'. update tests/collision-allowlist.txt to reflect the new kinds, or rename one." >&2
            fail=1
        fi
    done < <(collisions_enumerate)

    # Restore the outer SCIAGENT_TOOLKIT (if any).
    if [[ -n "$_saved_tk" ]]; then
        export SCIAGENT_TOOLKIT="$_saved_tk"
    else
        unset SCIAGENT_TOOLKIT
    fi

    if [[ "$fail" -ne 0 ]]; then
        exit 1
    fi
}

check_cross_namespace_collisions

pass
