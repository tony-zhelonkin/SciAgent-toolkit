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
# results-layout / captions / provenance / freshness / hooks / docs-layout /
# ...) have moved to lib/sciagent/lint.sh — `sciagent lint --check <name>`.
# `sciagent validate --check <name>` still works: it delegates to
# _lint_run_checks and prints a one-line deprecation note to stderr
# (suppressed by --quiet). The default path and the activate-internal
# `cmd_validate --quiet` call pass NO --check, so they behave EXACTLY as
# before (never blocked by a lint finding).
#
# docs-layout (docs/_internal/ gitignore-status, results-root .md files,
# mixed archive-naming, non-standard handoff names) used to live here as
# _validate_docs_layout, called unconditionally on validate's default path —
# that was itself the bug: a PROJECT finding (a consumer repo's docs/ layout)
# could hard-fail the default `validate` path, and since `activate` calls
# `cmd_validate --quiet` as pre-flight, a docs/_internal gitignore miss could
# hard-block activation entirely, contradicting this file's own contract
# above. It now lives in lint.sh as `--check docs-layout`, opt-in like its
# siblings, never reachable from activate's pre-flight.
#
# env-hygiene (warns when CLAUDE_CODE_SKIP_PROMPT_HISTORY silently disables
# session persistence/backgrounding) is not a PROJECT check at all — it reads
# the invoking shell's environment, not anything under --project-dir — so it
# does not belong in lint.sh's "against a consumer project" layer either. It
# stays here, but off the default path: pass --env-hygiene to run it. Always
# soft-warn; never affects the exit code.
#
# Hardness boundary:
#   Hard-fail (exit 1): malformed skill frontmatter;
#                       any --check finding when --strict (delegated to lint.sh).
#   Soft-warn:          skills-ref findings (when present);
#                       cross-namespace name collisions;
#                       --env-hygiene findings (always, regardless of --strict);
#                       --check findings without --strict.
#
# Exit code:
#   0 — all checks pass
#   1 — any hard check fails

# Ceiling on a skill `description:`. Harnesses preload name+description for
# every skill in the catalog, so this is the one frontmatter field with a
# per-context cost; the body is read only on activation.
: "${SCIAGENT_DESC_MAX:=350}"

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

# _validate_env_hygiene [--quiet]
# Environment-hygiene check (soft warn only; never hard-fails by default).
# Warns when CLAUDE_CODE_SKIP_PROMPT_HISTORY is set/truthy in the environment:
# in current Claude Code builds this flag does not merely drop up-arrow input
# history — it forces the whole "session persistence is disabled" state, which
# in turn disables backgrounding and the agents/teams panel ("this conversation
# cannot be backgrounded"). This is exactly why a host session couldn't be
# backgrounded (a stray `export CLAUDE_CODE_SKIP_PROMPT_HISTORY=1` in ~/.bashrc).
# Use ~/.claude/settings.json `cleanupPeriodDays` for bounded retention instead.
# Off the default path — only runs when `sciagent validate --env-hygiene` is
# given. --quiet suppresses the soft-warn stdout.
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
    local _env_hygiene=0
    local -a _checks=()
    while [[ $# -gt 0 ]]; do
        case "$1" in
            -h|--help)
                cat <<'USAGE'
sciagent validate [--quiet] [--project-dir <dir>] [--check <name>...] [--strict] [--env-hygiene]
  Check the frontmatter shape of every skill in the toolkit. Exits 0 on
  success, 1 on any hard-fail (missing/mismatched name, missing
  description, or a description over 350 chars).

  --quiet            Suppress the "all checks passed" summary on success.
  --project-dir <d>  Project directory for --check (default: .).
  --check <name>     DEPRECATED — delegates to `sciagent lint --check <name>`.
                     name ∈ figure-style | results-layout | captions |
                            provenance | freshness | hooks | docs-layout |
                            stage-thinness | comment-intent | stage-layout |
                            all. Repeatable.
                     When given, ONLY the named project check(s) run (the
                     toolkit-wide frontmatter walk is skipped).
                     Findings are soft WARN (exit 0) unless --strict.
  --strict           Project-check findings become HARD failures (exit 1).
  --env-hygiene      Also warn if CLAUDE_CODE_SKIP_PROMPT_HISTORY is set in
                     the invoking shell's environment (disables session
                     persistence/backgrounding). Off by default; not a
                     project check — not affected by --project-dir/--strict.

  _scratch/ and $TMPDIR are always exempt from project checks.
USAGE
                return 0 ;;
            --quiet)
                quiet=1; shift ;;
            --strict)
                _strict=1; shift ;;
            --env-hygiene)
                _env_hygiene=1; shift ;;
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
                    echo "usage: sciagent validate [--quiet] [--project-dir <dir>] [--check <name>] [--strict] [--env-hygiene]" >&2
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
        # Reserved underscore-prefixed directories stay off the active path.
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

        # `description:` length must be measured on the value a YAML parser
        # actually produces, not on the raw source line. A folded/literal
        # block scalar (`>`, `>-`, `>+`, `|`, `|-`, `|+`) puts NO content on
        # the `description:` line itself -- the value lives in the indented
        # continuation lines below it. Measuring only the first physical line
        # (as a naive `print length($0)` would) sees just the 2-char style
        # indicator and silently passes any length, no matter how long the
        # folded value actually is. So: extract the frontmatter block once,
        # confirm a `description:` key exists at all, then fold it the way
        # YAML does --
        #   - plain / single-quoted / double-quoted scalars: single physical
        #     line, quotes stripped (minimal unescaping of \" and '').
        #   - `>`/`>-`/`>+` (folded): continuation lines are accumulated,
        #     stripped of the common indentation set by the first
        #     continuation line, and joined with a space; a blank line
        #     between content lines folds to a literal newline instead
        #     (one newline per blank line, matching YAML line folding).
        #   - `|`/`|-`/`|+` (literal): continuation lines are joined with
        #     newlines verbatim (indentation-stripped only).
        #   - chomping indicator (`-` strip / `+` keep / default clip) is
        #     honored: clip appends exactly one trailing newline when the
        #     folded content is non-empty, keep appends the trailing blank
        #     lines actually present plus one, strip appends none.
        # A continuation line is only "in" the block while its indentation is
        # >= the first continuation line's indentation (the base); the block
        # ends at the first equal-or-lower-indented line, at end of
        # frontmatter, or immediately (empty value) if the very next line is
        # not indented past the `description:` key itself.
        # Verified byte-for-byte against a PyYAML parse of the whole skills/
        # corpus (83 skills, all styles present in-repo) plus synthetic cases
        # covering every chomp/style combination and blank-line folding.
        local skill_fm fm_desc fm_desc_len
        skill_fm="$(awk '
            NR==1 && /^---[ \t]*$/ { infm=1; next }
            infm && /^---[ \t]*$/ { exit }
            infm { print }
        ' "$skill_file")"
        if ! printf '%s\n' "$skill_fm" | grep -q '^description:'; then
            failures+=("$skill_name: SKILL.md frontmatter has no description:")
            fail=1
        else
            fm_desc="$(printf '%s\n' "$skill_fm" | awk '
                function leadspaces(s) {
                    n = 0
                    while (n < length(s) && substr(s, n+1, 1) == " ") n++
                    return n
                }
                function finalize() {
                    out = ""
                    last = nlines
                    while (last > 0 && blank[last]) last--
                    if (style == "|") {
                        for (i=1; i<=last; i++) out = out (i>1 ? "\n" : "") content[i]
                        if (chomp == "clip") { if (last > 0) out = out "\n" }
                        else if (chomp == "keep") {
                            for (i=last+1; i<=nlines; i++) out = out "\n"
                            if (last > 0 || nlines > 0) out = out "\n"
                        }
                    } else {
                        blankrun = 0
                        first = 1
                        for (i=1; i<=last; i++) {
                            if (blank[i]) { blankrun++; continue }
                            if (!first) {
                                if (blankrun > 0) { for (b=0; b<blankrun; b++) out = out "\n" }
                                else out = out " "
                            }
                            out = out content[i]
                            first = 0
                            blankrun = 0
                        }
                        if (chomp == "clip") { if (last > 0) out = out "\n" }
                        else if (chomp == "keep") {
                            for (i=last+1; i<=nlines; i++) out = out "\n"
                            if (last > 0 || nlines > 0) out = out "\n"
                        }
                    }
                    printf "%s", out
                    state = 2
                }
                BEGIN { state = 0 }
                state==0 && /^description:[ \t]*/ {
                    key_indent = leadspaces($0)
                    line = $0
                    sub(/^description:[ \t]*/, "", line)
                    sub(/[ \t]+$/, "", line)
                    if (line ~ /^[|>][+-]?[0-9]*[ \t]*(#.*)?$/) {
                        style = substr(line, 1, 1)
                        rest = substr(line, 2)
                        chomp = "clip"
                        if (rest ~ /^-/) chomp = "strip"
                        else if (rest ~ /^\+/) chomp = "keep"
                        state = 1
                        nlines = 0
                        have_base = 0
                        next
                    } else {
                        val = line
                        if (val ~ /^".*"$/ && length(val) >= 2) {
                            val = substr(val, 2, length(val)-2)
                            gsub(/\\"/, "\"", val)
                        } else if (val ~ /^'"'"'.*'"'"'$/ && length(val) >= 2) {
                            val = substr(val, 2, length(val)-2)
                            gsub(/'"'"''"'"'/, "'"'"'", val)
                        }
                        printf "%s", val
                        state = 2
                        exit
                    }
                }
                state==1 {
                    if ($0 ~ /^[ \t]*$/) {
                        nlines++
                        blank[nlines] = 1
                        content[nlines] = ""
                        next
                    }
                    ls = leadspaces($0)
                    if (!have_base) {
                        if (ls <= key_indent) { finalize(); exit }
                        base = ls; have_base = 1
                    }
                    if (ls < base) { finalize(); exit }
                    nlines++
                    blank[nlines] = 0
                    content[nlines] = substr($0, base+1)
                    next
                }
                END { if (state == 1) finalize() }
            '; printf 'X')"
            # The trailing sentinel defeats command-substitution's stripping
            # of trailing newlines, which would otherwise undercount a
            # clip-chomped (`>`/`|`, no `-`) folded value by one character.
            fm_desc="${fm_desc%X}"
            fm_desc_len=${#fm_desc}
            if (( fm_desc_len > SCIAGENT_DESC_MAX )); then
                failures+=("$skill_name: description is $fm_desc_len chars (max $SCIAGENT_DESC_MAX)")
                fail=1
            fi
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

    # Check 3: cross-namespace collisions (soft-warn).
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

    # Optional environment-hygiene check (--env-hygiene only; off by default).
    # Soft warn only — never touches $fail — but it still runs BEFORE the
    # "all checks passed" line below, on the same principle as every other
    # check above: nothing that prints to this invocation runs after the
    # success summary. (Historically this — and the now-relocated
    # docs-layout project check — ran AFTER the summary line, so a run could
    # print "all checks passed" and then still fail/warn afterward. Project
    # checks moved to lint.sh entirely; this ordering fix is what's left of
    # that bug for any check that stays in cmd_validate itself.)
    if [[ "$_env_hygiene" -eq 1 ]]; then
        if [[ "$quiet" -eq 1 ]]; then
            _validate_env_hygiene --quiet
        else
            _validate_env_hygiene
        fi
    fi

    if [[ "$quiet" -eq 0 ]]; then
        echo "sciagent validate: all checks passed"
    fi

    [[ "$fail" -eq 0 ]] && return 0 || return 1
}
