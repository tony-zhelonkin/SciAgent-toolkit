#!/usr/bin/env bash
# tests/test_no_duplicate_basenames.sh
# Enforces the uniqueness invariant for the resolver: no two .md files under
# commands/ or agents/ (recursively) may share a basename, and no two skill
# directories at skills/<name> may share a name. The recursive resolver in
# lib/sciagent/symlinks.sh (resolve_canonical) depends on this — duplicate
# basenames would make role bindings ambiguous.

set -u
# shellcheck source=/dev/null
. "$(dirname "$0")/_lib.sh"

TOOLKIT="$TOOLKIT_ROOT"

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

pass
