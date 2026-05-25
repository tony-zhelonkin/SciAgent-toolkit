# lib/sciagent/validate.sh — sciagent validate [--quiet]
#
# Per kickoff.md §9 ADR-007 dissolution: thin callable verb, not a subsystem.
# Called internally by `activate` before mounting; standalone for debugging.
#
# Checks (per kickoff.md §9, ADR-007 resolution 2026-05-24):
#   1. requires resolution    — every named dep exists as a skill
#   2. cycle detection        — DFS over requires graph (skill_deps.sh helper)
#   3. tag-vocab compliance   — every metadata.tags: entry is in tags.yaml
#   4. (optional) skills-ref  — if installed, invoke per skill; surface exit
#                               code as warning, not error; silently skip
#                               when absent
#
# Hardness boundary (kickoff.md §9, PR 2 2026-05-24):
#   Hard-fail (exit 1): cycle or missing requires target; unknown tag.
#   Soft-warn:          skills-ref findings (when present).
#
# Exit code:
#   0 — all checks pass
#   1 — any hard check fails
#
# Usage: cmd_validate [--quiet]
#   --quiet  suppress "all checks passed" summary on success

# shellcheck shell=bash

cmd_validate() {
    local quiet=0
    [[ "${1:-}" == "--quiet" ]] && quiet=1

    local fail=0
    local -a failures=()

    # Locate the toolkit root. SCIAGENT_TOOLKIT is exported by bin/sciagent;
    # fall back to a relative path from this file's own location.
    local tk_root
    if [[ -n "${SCIAGENT_TOOLKIT:-}" ]]; then
        tk_root="$SCIAGENT_TOOLKIT"
    else
        local self_dir
        self_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
        tk_root="$(cd "$self_dir/../.." && pwd)"
    fi

    local skills_dir="$tk_root/skills"
    local tags_file="$tk_root/tags.yaml"

    # -----------------------------------------------------------------------
    # Check 3 prep: build known-tag vocabulary from tags.yaml.
    # Parse lines:  `  - name: <value>`
    # -----------------------------------------------------------------------
    local known_tags=""
    if [[ -f "$tags_file" ]]; then
        known_tags="$(awk '
            /^tags:/ { intags=1; next }
            intags && /^  - name:/ {
                sub(/^  - name:[ \t]*/, "")
                sub(/[ \t]*#.*$/, "")
                sub(/[ \t]+$/, "")
                gsub(/^["'"'"']|["'"'"']$/, "")
                print
                next
            }
            intags && /^[^ ]/ { intags=0 }
        ' "$tags_file")"
    else
        failures+=("tags.yaml: file missing at $tags_file")
        fail=1
    fi

    # -----------------------------------------------------------------------
    # Walk every skill directory.
    # -----------------------------------------------------------------------
    local skill_dir skill_name skill_file
    for skill_dir in "$skills_dir"/*/; do
        skill_name="$(basename "$skill_dir")"
        [[ "$skill_name" == "_TEMPLATE" ]] && continue
        skill_file="$skill_dir/SKILL.md"
        [[ -f "$skill_file" ]] || continue

        # Check 1+2: transitive requires resolution + cycle detection.
        # skill_resolve_transitive hard-fails on missing targets and cycles.
        local resolve_err
        if ! resolve_err=$(skill_resolve_transitive "$skill_name" 2>&1); then
            failures+=("$skill_name: requires resolution failed — $resolve_err")
            fail=1
        fi

        # Check 3: tag-vocab compliance.
        # Parse metadata.tags block list items from frontmatter.
        if [[ -n "$known_tags" ]]; then
            local skill_tags
            skill_tags="$(awk '
                BEGIN { infm=0; closed=0; inmeta=0; intags=0 }
                NR==1 && /^---[ \t]*$/ { infm=1; next }
                infm && !closed && /^---[ \t]*$/ { closed=1; exit }
                !infm { next }
                /^metadata:[ \t]*$/ { inmeta=1; next }
                inmeta && /^  tags:/ { intags=1; next }
                intags && /^[ \t]+-[ \t]/ {
                    val=$0
                    sub(/^[ \t]+-[ \t]+/, "", val)
                    sub(/[ \t]*#.*$/, "", val)
                    sub(/[ \t]+$/, "", val)
                    gsub(/^["'"'"']|["'"'"']$/, "", val)
                    print val
                    next
                }
                intags && /^  [^ \t-]/ { intags=0 }
                intags && /^[^ \t]/ { intags=0 }
            ' "$skill_file")"

            local tag
            while IFS= read -r tag; do
                [[ -z "$tag" ]] && continue
                if ! printf '%s\n' "$known_tags" | grep -qxF "$tag"; then
                    failures+=("$skill_name: unknown tag '$tag' (not in tags.yaml)")
                    fail=1
                fi
            done <<< "$skill_tags"
        fi

        # Check 4 (optional): skills-ref shell-out.
        if command -v skills-ref >/dev/null 2>&1; then
            if ! skills-ref "$skill_file" >/dev/null 2>&1; then
                # Warn-only; does not set fail.
                echo "sciagent validate: warning: skills-ref reported issues for '$skill_name'" >&2
            fi
        fi
    done

    # -----------------------------------------------------------------------
    # Report.
    # -----------------------------------------------------------------------
    if [[ "$fail" -ne 0 ]]; then
        echo "sciagent validate: checks failed:" >&2
        local f
        for f in "${failures[@]}"; do
            echo "  - $f" >&2
        done
        return 1
    fi

    if [[ "$quiet" -eq 0 ]]; then
        echo "sciagent validate: all checks passed"
    fi
    return 0
}
