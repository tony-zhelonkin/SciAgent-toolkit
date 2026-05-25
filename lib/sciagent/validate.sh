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
#   5. cross-namespace collision — names appearing in >=2 of
#                               skills/agents/commands/roles. Soft-warn only;
#                               most overlaps are intentional family overlaps
#                               (e.g. /architect → @architect → roles/architect).
#                               tests/collision-allowlist.txt is the CI mirror
#                               that turns accidental overlaps into a merge gate.
#
# Hardness boundary (kickoff.md §9, PR 2 2026-05-24):
#   Hard-fail (exit 1): cycle or missing requires target; unknown tag.
#   Soft-warn:          skills-ref findings (when present);
#                       cross-namespace name collisions.
#
# Exit code:
#   0 — all checks pass
#   1 — any hard check fails
#
# Output streams:
#   stdout — "all checks passed" summary line on success (suppressed by --quiet)
#   stderr — per-failure stanza on hard-fail; skills-ref warnings when present
#
# Usage: cmd_validate [--quiet]
#   --quiet  suppress "all checks passed" summary on success

# shellcheck shell=bash

cmd_validate() {
    local quiet=0
    while [[ $# -gt 0 ]]; do
        case "$1" in
            -h|--help)
                cat <<'USAGE'
sciagent validate [--quiet]
  Check requires graph + tag-vocab compliance over every skill in the
  toolkit. Exits 0 on success, 1 on any hard-fail (missing requires
  target, cycle, or unknown tag).

  --quiet   Suppress the "all checks passed" summary on success.
USAGE
                return 0 ;;
            --quiet)
                quiet=1; shift ;;
            *)
                echo "sciagent validate: unknown option '$1'" >&2
                echo "usage: sciagent validate [--quiet]" >&2
                return 1 ;;
        esac
    done

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
    # NOTE: tests/test_tags_vocabulary.sh carries a near-identical awk parser
    # for the same vocabulary. If you change the matching/trimming rules here,
    # change them there too — the two readers must agree on what counts as a tag.
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
        # skill_resolve_transitive hard-fails on missing targets and cycles
        # and writes the diagnostic to stderr (which we capture into
        # resolve_err and re-emit). The kickoff spec asks for "message names
        # the cycle" — the resolver's "cycle detected in requires graph at
        # '<name>'" line satisfies that transitively, so we relay it as-is
        # rather than re-format here.
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

    # Check 5: cross-namespace collisions (soft-warn).
    # Buffered then emitted after the "all checks passed" line so a clean
    # tree still produces zero stderr output. Quiet mode mutes the warnings
    # for the same reason it mutes the success line — scripted callers want
    # silent-on-success. Allowlist annotation lives in status.sh; validate
    # speaks namespace-level (it doesn't know which collisions are "blessed").
    if [[ "$quiet" -eq 0 ]]; then
        local -a _collision_warnings=()
        local _col_name _col_kinds _col_k1 _col_k2
        while IFS=$'\t' read -r _col_name _col_kinds; do
            [[ -z "$_col_name" ]] && continue
            _col_k1="${_col_kinds%%,*}"
            _col_k2="${_col_kinds#*,}"
            _col_k2="${_col_k2%%,*}"
            _collision_warnings+=("validate: warning — name '$_col_name' appears as both $_col_k1 and $_col_k2 (mounting both is supported; ensure the overlap is intentional)")
        done < <(collisions_enumerate)
        local w
        for w in "${_collision_warnings[@]+"${_collision_warnings[@]}"}"; do
            echo "$w" >&2
        done
    fi

    if [[ "$quiet" -eq 0 ]]; then
        echo "sciagent validate: all checks passed"
    fi
    return 0
}
