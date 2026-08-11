# lib/sciagent/validate.sh — sciagent validate [--quiet]
#
# Thin callable verb, not a subsystem.
# Called internally by `activate` before mounting; standalone for debugging.
#
# Checks:
#   1. frontmatter shape      — every skill has `name:` matching its directory
#                               and a `description:` of at most 350 chars
#   2. (optional) skills-ref  — if installed, invoke per skill; surface exit
#                               code as warning, not error; silently skip
#                               when absent
#   3. cross-namespace collision — names appearing in >=2 of
#                               skills/agents/commands/roles. Soft-warn only;
#                               most overlaps are intentional family overlaps
#                               (e.g. /architect → @architect → roles/architect).
#                               tests/collision-allowlist.txt is the CI mirror
#                               that turns accidental overlaps into a merge gate.
#
# Opt-in PROJECT guardrail checks (the (c) GUARDRAIL layer, figure-style /
# results-layout / captions / provenance / freshness / hooks / ...) have moved
# to lib/sciagent/lint.sh — `sciagent lint --check <name>`. `sciagent validate
# --check <name>` still works: it delegates to _lint_run_checks and prints a
# one-line deprecation note to stderr (suppressed by --quiet). The default
# path and the activate-internal `cmd_validate --quiet` call pass NO --check,
# so they behave EXACTLY as before (never blocked by a lint finding).
#
# Hardness boundary:
#   Hard-fail (exit 1): malformed skill frontmatter;
#                       any --check finding when --strict (delegated to lint.sh).
#   Soft-warn:          skills-ref findings (when present);
#                       cross-namespace name collisions;
#                       --check findings without --strict.
#
# Exit code:
#   0 — all checks pass
#   1 — any hard check fails

# Ceiling on a skill `description:`. Harnesses preload name+description for
# every skill in the catalog, so this is the one frontmatter field with a
# per-context cost; the body is read only on activation.
: "${SCIAGENT_DESC_MAX:=350}"
#
# Output streams:
#   stdout — "all checks passed" summary line on success (suppressed by --quiet)
#   stderr — per-failure stanza on hard-fail; skills-ref warnings when present
#
# Usage: cmd_validate [--quiet]
#   --quiet  suppress "all checks passed" summary on success

# shellcheck shell=bash

# _validate_join_kinds <csv>
# Render a comma-separated kinds list (as emitted by collisions_enumerate) as
# a natural-language phrase. Scales from 2 kinds ("both A and B") through 3+
# ("A, B, and C"). The 1-kind path is unreachable from the current caller
# (collisions_enumerate only emits names with >=2 namespace matches) but is
# handled defensively so this helper is callable from elsewhere later.
_validate_join_kinds() {
    local csv="$1"
    local -a parts=()
    local IFS=','
    read -r -a parts <<< "$csv"
    unset IFS
    case ${#parts[@]} in
        0) printf '' ;;
        1) printf '%s' "${parts[0]}" ;;
        2) printf 'both %s and %s' "${parts[0]}" "${parts[1]}" ;;
        *)
            local last_idx=$(( ${#parts[@]} - 1 ))
            local i out=""
            for (( i=0; i<last_idx; i++ )); do
                out+="${parts[$i]}, "
            done
            out+="and ${parts[$last_idx]}"
            printf '%s' "$out"
            ;;
    esac
}

# _validate_docs_layout <projdir> [--quiet]
# Lints the docs/ layout of a project directory for common structural problems.
# Hard failures (exit 1): docs/_internal/ exists, is inside a git repo, and is
#   NOT gitignored (would expose internal notes in a push).
# Soft warnings (exit 0): missing docs/, missing _internal/, .md in 03_results/,
#   mixed archive-naming conventions, non-standard handoff filenames.
# --quiet suppresses soft-warn stdout; hard-fail ERROR still goes to stderr.
_validate_docs_layout() {
    local projdir="${1:-.}"
    local _quiet=0
    if [[ "${2:-}" == "--quiet" ]]; then
        _quiet=1
    fi
    local _docs_fail=0

    # Check A (soft warn): docs/ directory does not exist.
    [[ "$_quiet" -eq 0 && ! -d "$projdir/docs" ]] && echo "WARN docs: no docs/ directory found"

    # Check B (soft warn): docs/_internal/ does not exist.
    [[ "$_quiet" -eq 0 && ! -d "$projdir/docs/_internal" ]] && echo "WARN docs: docs/_internal/ missing — run: sciagent gitignore"

    # Check C (HARD fail): docs/_internal/ exists but is NOT gitignored.
    if [[ -d "$projdir/docs/_internal" ]]; then
        local _in_git
        _in_git=$(git -C "$projdir" rev-parse --is-inside-work-tree 2>/dev/null)
        if [[ "$_in_git" == "true" ]]; then
            if ! git -C "$projdir" check-ignore -q docs/_internal 2>/dev/null; then
                echo "ERROR docs: docs/_internal/ is NOT gitignored — add 'docs/_internal/' to .gitignore" >&2
                _docs_fail=1
            fi
        fi
    fi

    # Check D (soft warn): any .md file directly in 03_results/ at maxdepth 1.
    if [[ "$_quiet" -eq 0 ]]; then
        find "$projdir/03_results" -maxdepth 1 -name "*.md" 2>/dev/null | while IFS= read -r f; do
            echo "WARN docs: report in results dir: $(basename "$f") — move to docs/_internal/reports/"
        done
    fi

    # Check E (soft warn): multiple archive-style naming conventions coexist.
    # Look for dirs named .archive, _deprecated, _legacy, .deprecated.
    if [[ "$_quiet" -eq 0 ]]; then
        local _dir
        for _dir in "$projdir/02_analysis" "$projdir/03_results"; do
            [[ -d "$_dir" ]] || continue
            local _conventions
            _conventions=$(find "$_dir" -maxdepth 3 -type d \
                \( -name ".archive" -o -name "_deprecated" -o -name "_legacy" -o -name ".deprecated" \) \
                2>/dev/null | awk -F'/' '{
                    # get the parent path (all but last component)
                    n=split($0,a,"/")
                    parent=""
                    for(i=1;i<n;i++) parent=parent (i>1?"/":"") a[i]
                    # record which convention names appear under each parent
                    basename=a[n]
                    key=parent SUBSEP basename
                    if(!seen[key]++) {
                        count[parent]++
                    }
                }
                END {
                    for(p in count) {
                        if(count[p]>=2) print p
                    }
                }')
            if [[ -n "$_conventions" ]]; then
                local _p
                while IFS= read -r _p; do
                    echo "WARN docs: mixed archive-naming conventions under $_p — pick one (.archive / _deprecated / _legacy / .deprecated)"
                done <<< "$_conventions"
            fi
        done
    fi

    # Check F (soft warn): handoff_*.md files at repo root that don't match
    # handoff_YYYYMMDD_HHMMSS.md
    if [[ "$_quiet" -eq 0 ]]; then
        local f
        for f in "$projdir"/handoff_*.md; do
            [[ -f "$f" ]] || continue
            if ! printf '%s\n' "$(basename "$f")" | grep -qE '^handoff_[0-9]{8}_[0-9]{6}\.md$'; then
                echo "WARN docs: non-standard handoff filename: $(basename "$f")"
            fi
        done
    fi

    return $_docs_fail
}

# _validate_env_hygiene [--quiet]
# Environment-hygiene check (soft warn only; never hard-fails by default).
# Warns when CLAUDE_CODE_SKIP_PROMPT_HISTORY is set/truthy in the environment:
# in current Claude Code builds this flag does not merely drop up-arrow input
# history — it forces the whole "session persistence is disabled" state, which
# in turn disables backgrounding and the agents/teams panel ("this conversation
# cannot be backgrounded"). This is exactly why a host session couldn't be
# backgrounded (a stray `export CLAUDE_CODE_SKIP_PROMPT_HISTORY=1` in ~/.bashrc).
# Use ~/.claude/settings.json `cleanupPeriodDays` for bounded retention instead.
# --quiet suppresses the soft-warn stdout (mirrors _validate_docs_layout).
_validate_env_hygiene() {
    local _quiet=0
    if [[ "${1:-}" == "--quiet" ]]; then
        _quiet=1
    fi

    # Truthy = set and not one of the empty/0/false/off spellings.
    local v="${CLAUDE_CODE_SKIP_PROMPT_HISTORY:-}"
    case "$v" in
        ""|0|false|FALSE|off|OFF|no|NO) ;;
        *)
            [[ "$_quiet" -eq 0 ]] && \
                echo "WARN session-persistence: CLAUDE_CODE_SKIP_PROMPT_HISTORY is set — this disables Claude Code session persistence, which also disables backgrounding + agents/teams (a session cannot be backgrounded). Unset it (use settings.json cleanupPeriodDays for bounded retention instead)."
            ;;
    esac

    return 0
}

# ===========================================================================
# PROJECT GUARDRAIL CHECKS (opt-in via --check) have moved to lib/sciagent/
# lint.sh — `sciagent lint`. cmd_validate's --check branch below delegates to
# _lint_run_checks for backward compatibility. See lint.sh for the (c)
# GUARDRAIL layer's checks, helpers, and the _VCHECK_STAGE_DIRS migration
# window constant.
# ===========================================================================

cmd_validate() {
    local quiet=0
    local _projdir="."
    local _strict=0
    local -a _checks=()
    while [[ $# -gt 0 ]]; do
        case "$1" in
            -h|--help)
                cat <<'USAGE'
sciagent validate [--quiet] [--project-dir <dir>] [--check <name>...] [--strict]
  Check the frontmatter shape of every skill in the toolkit. Exits 0 on
  success, 1 on any hard-fail (missing/mismatched name, missing
  description, or a description over 350 chars).

  --quiet            Suppress the "all checks passed" summary on success.
  --project-dir <d>  Project directory to lint (default: .).
  --check <name>     DEPRECATED — delegates to `sciagent lint --check <name>`.
                     name ∈ figure-style | results-layout | captions |
                            provenance | freshness | hooks | stage-thinness |
                            comment-intent | stage-layout | all. Repeatable.
                     When given, ONLY the named project check(s) run (the
                     toolkit-wide frontmatter walk is skipped).
                     Findings are soft WARN (exit 0) unless --strict.
  --strict           Project-check findings become HARD failures (exit 1).

  _scratch/ and $TMPDIR are always exempt from project checks.
USAGE
                return 0 ;;
            --quiet)
                quiet=1; shift ;;
            --strict)
                _strict=1; shift ;;
            --check)
                _checks+=("${2:-}"); shift 2 ;;
            --project-dir)
                _projdir="${2:-.}"; shift 2 ;;
            *)
                # Accept first positional non-flag arg as project dir.
                if [[ "$1" != -* ]]; then
                    _projdir="$1"; shift
                else
                    echo "sciagent validate: unknown option '$1'" >&2
                    echo "usage: sciagent validate [--quiet] [--project-dir <dir>] [--check <name>] [--strict]" >&2
                    return 1
                fi ;;
        esac
    done

    # -----------------------------------------------------------------------
    # Opt-in PROJECT guardrail checks (--check). These are the (c) GUARDRAIL
    # layer. CRITICAL: they run ONLY when --check is given. The default path
    # and the activate-internal `cmd_validate --quiet` call (which passes no
    # --check) MUST behave EXACTLY as before — never let a figure/caption/
    # layout finding block activation. So when --check is present we run the
    # selected project checks and return; we do NOT run the toolkit-wide walk.
    # -----------------------------------------------------------------------
    if [[ ${#_checks[@]} -gt 0 ]]; then
        [[ "$quiet" -eq 0 ]] && echo "sciagent validate --check is deprecated; use: sciagent lint --check ..." >&2
        _lint_run_checks "$_projdir" "$_strict" "$quiet" "${_checks[@]}"
        return $?
    fi

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

    # -----------------------------------------------------------------------
    # Walk every skill directory.
    # -----------------------------------------------------------------------
    local skill_dir skill_name skill_file
    for skill_dir in "$skills_dir"/*/; do
        skill_name="$(basename "$skill_dir")"
        # Skip non-skill scaffolding: _TEMPLATE (copy-target), _attic (retired,
        # reference-only — see docs/skill-lifecycle.md), and any other
        # underscore-prefixed holding dir (_archive backups, etc.).
        [[ "$skill_name" == _* ]] && continue
        skill_file="$skill_dir/SKILL.md"
        [[ -f "$skill_file" ]] || continue

        # Check 1: frontmatter shape. `name` and `description` are the whole
        # preloaded surface a harness sees — a missing or mismatched name makes
        # the skill unaddressable, and an oversized description spends context
        # on every skill in the catalog whether it is used or not.
        local fm_name fm_desc_len
        fm_name="$(awk '
            NR==1 && /^---[ \t]*$/ { infm=1; next }
            infm && /^---[ \t]*$/ { exit }
            infm && /^name:[ \t]*/ {
                sub(/^name:[ \t]*/, "")
                sub(/[ \t]+$/, "")
                gsub(/^["'"'"']|["'"'"']$/, "")
                print; exit
            }
        ' "$skill_file")"
        if [[ -z "$fm_name" ]]; then
            failures+=("$skill_name: SKILL.md frontmatter has no name:")
            fail=1
        elif [[ "$fm_name" != "$skill_name" ]]; then
            failures+=("$skill_name: frontmatter name: '$fm_name' does not match its directory")
            fail=1
        fi

        fm_desc_len="$(awk '
            NR==1 && /^---[ \t]*$/ { infm=1; next }
            infm && /^---[ \t]*$/ { exit }
            infm && /^description:[ \t]*/ {
                sub(/^description:[ \t]*/, "")
                sub(/[ \t]+$/, "")
                gsub(/^["'"'"']|["'"'"']$/, "")
                print length($0); exit
            }
        ' "$skill_file")"
        if [[ -z "$fm_desc_len" ]]; then
            failures+=("$skill_name: SKILL.md frontmatter has no description:")
            fail=1
        elif (( fm_desc_len > SCIAGENT_DESC_MAX )); then
            failures+=("$skill_name: description is $fm_desc_len chars (max $SCIAGENT_DESC_MAX)")
            fail=1
        fi

        # Check 2 (optional): skills-ref shell-out.
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
        local _col_name _col_kinds
        while IFS=$'\t' read -r _col_name _col_kinds; do
            [[ -z "$_col_name" ]] && continue
            # collisions_enumerate emits a kinds csv that can carry 2, 3, or 4
            # entries. Render it as a natural-language list so the warning
            # stays accurate beyond the two-kind case (the prior templated
            # "both A and B" form silently dropped any third or fourth kind).
            local _kinds_phrase
            _kinds_phrase=$(_validate_join_kinds "$_col_kinds")
            _collision_warnings+=("validate: warning — name '$_col_name' appears as $_kinds_phrase (mounting both is supported; ensure the overlap is intentional)")
        done < <(collisions_enumerate)
        local w
        for w in "${_collision_warnings[@]+"${_collision_warnings[@]}"}"; do
            echo "$w" >&2
        done
    fi

    if [[ "$quiet" -eq 0 ]]; then
        echo "sciagent validate: all checks passed"
    fi

    # Docs-layout linter (runs after all other checks; hard-fail from it
    # increments fail and causes a non-zero exit).
    # Pass --quiet through so soft-warn stdout is suppressed in quiet mode.
    if [[ "$quiet" -eq 1 ]]; then
        _validate_docs_layout "$_projdir" --quiet
    else
        _validate_docs_layout "$_projdir"
    fi
    local _docs_rc=$?
    (( _docs_rc != 0 )) && fail=$(( fail + 1 ))

    # Environment-hygiene check (soft warn only; never touches the exit code).
    # Warns on CLAUDE_CODE_SKIP_PROMPT_HISTORY (silently disables session
    # persistence + backgrounding/agents). Quiet-aware like the docs linter.
    if [[ "$quiet" -eq 1 ]]; then
        _validate_env_hygiene --quiet
    else
        _validate_env_hygiene
    fi

    [[ "$fail" -eq 0 ]] && return 0 || return 1
}
