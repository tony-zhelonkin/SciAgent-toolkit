#!/usr/bin/env bash
# tests/test_skill_scope_lint.sh — enforce per-scope body line caps on skills.
#
# ADR-0002 §3.3: count body lines (after frontmatter close to EOF),
# EXCLUDING fenced code-block content.
#
# Caps:
#   atomic       <= 300
#   orchestrator <= 250
#   foundation   <= 800
#
# Legacy cutoff: a skill whose `metadata.last-reviewed` is absent OR strictly
# before 2026-05-21 is exempt from hard-fail (warning only). Skills reviewed
# on or after the cutoff MUST pass the cap.

set -u
. "$(dirname "$0")/_lib.sh"

CUTOFF="2026-05-21"

# Compute body line count excluding fenced code blocks.
body_loc_no_fences() {
    awk '
        BEGIN { infm = 0; closed = 0; infence = 0; n = 0 }
        NR == 1 && /^---[ \t]*$/ { infm = 1; next }
        infm && !closed && /^---[ \t]*$/ { closed = 1; next }
        !closed { next }
        /^```/ { infence = !infence; next }
        !infence { n++ }
        END { print n }
    ' "$1"
}

# Extract a scalar value from frontmatter top-level OR metadata: block.
# Usage: fm_scalar <file> <metadata-key>
fm_scalar() {
    local file="$1" key="$2"
    awk -v key="$key" '
        BEGIN { infm=0; closed=0; inmeta=0 }
        NR==1 && /^---[ \t]*$/ { infm=1; next }
        infm && !closed && /^---[ \t]*$/ { closed=1; exit }
        !infm { next }
        /^metadata:[ \t]*$/ { inmeta=1; next }
        /^[A-Za-z_][A-Za-z0-9_-]*:/ && !/^metadata:/ { inmeta=0 }
        inmeta {
            line = $0
            # match "  key: value"
            re = "^[ \t]+" key ":[ \t]*"
            if (line ~ re) {
                sub(re, "", line)
                sub(/[ \t]*#.*$/, "", line)
                sub(/[ \t]+$/, "", line)
                gsub(/^["'\'']|["'\'']$/, "", line)
                print line
                exit
            }
        }
    ' "$file"
}

declare -i fail_count=0
declare -i warn_count=0
declare -i pass_count=0

# Skip the _TEMPLATE skill — it has no frontmatter and is intentionally
# excluded from linting (it's literally a copy-target).
for skill_dir in "$TOOLKIT_ROOT"/skills/*/; do
    name="$(basename "$skill_dir")"
    [[ "$name" == "_TEMPLATE" ]] && continue

    file="$skill_dir/SKILL.md"
    [[ -f "$file" ]] || continue

    scope="$(fm_scalar "$file" scope)"
    [[ -z "$scope" ]] && scope="atomic"

    last_reviewed="$(fm_scalar "$file" last-reviewed)"

    case "$scope" in
        atomic)       cap=300 ;;
        orchestrator) cap=250 ;;
        foundation)   cap=800 ;;
        *)
            echo "FAIL [$_TEST_NAME] $name: unknown scope '$scope'" >&2
            fail_count+=1
            continue
            ;;
    esac

    loc="$(body_loc_no_fences "$file")"

    if [[ "$loc" -le "$cap" ]]; then
        pass_count+=1
        continue
    fi

    # Over cap — apply legacy exemption.
    if [[ -z "$last_reviewed" ]] || [[ "$last_reviewed" < "$CUTOFF" ]]; then
        echo "WARN [$_TEST_NAME] $name: scope=$scope body=$loc lines (cap=$cap) — legacy exemption (last-reviewed=${last_reviewed:-<none>})" >&2
        warn_count+=1
    else
        echo "FAIL [$_TEST_NAME] $name: scope=$scope body=$loc lines exceeds cap=$cap (last-reviewed=$last_reviewed >= $CUTOFF)" >&2
        fail_count+=1
    fi
done

if [[ "$fail_count" -gt 0 ]]; then
    echo "FAIL [$_TEST_NAME] $fail_count skill(s) over cap (post-cutoff); $warn_count legacy warning(s); $pass_count ok" >&2
    exit 1
fi

pass
