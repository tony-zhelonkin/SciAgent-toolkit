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
# Opt-in PROJECT guardrail checks (--check, the (c) GUARDRAIL layer):
#   figure-style    inline theme()/ggsave(width=)/figsize/raw-hex in viz
#                   scripts; missing project_theme()/set_paper_style(); config
#                   figures.base_size below the 14pt base font floor.
#   results-layout  artifacts must sit under <stage>/{figures,tables}/; stage
#                   must be a known stages: id; no artifact at 03_results/ root;
#                   each figure needs a same-stem table neighbor.
#   captions        every figures/tables artifact has a path-qualified
#                   `## <rel>` heading in the sibling stage README.md.
#   provenance      each caption's Script: cell resolves to an existing (and,
#                   in a git repo, tracked) 02_analysis/stages/ path.
#   freshness       project CRAFT block version / submodule commit vs toolkit.
#   hooks           every hook registered in .claude/settings.json exists and is
#                   executable (a registered-but-absent hook silently disables
#                   the (c) GUARDRAIL layer).
# These run ONLY when --check <name> (or --check all) is given. They are SOFT
# warnings (exit 0) by default and HARD failures (exit 1) under --strict.
# `_scratch/` and $TMPDIR are always exempt. The default + activate-internal
# call passes NO --check, so it behaves EXACTLY as before (never blocked by a
# figure/caption/layout finding).
#
# Hardness boundary:
#   Hard-fail (exit 1): malformed skill frontmatter;
#                       any --check finding when --strict.
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
# PROJECT GUARDRAIL CHECKS (opt-in via --check). The (c) GUARDRAIL layer.
#
# Each check is `_validate_check_<name> <projdir> <strict> <quiet>` and emits
# zero or more findings via `_vcheck_emit`. Findings are soft WARN (exit 0) by
# default; under <strict>=1 they become HARD failures (exit 1). The helper
# centralises the soft-warn vs hard-fail hardness boundary + stderr discipline.
#
# Invariant: an absent dir / a software (non-analysis) project / an empty
# project produces NO findings (clean exit 0). Every find/grep is guarded so a
# missing dir is a clean no-op. `_scratch/` and $TMPDIR are always exempt.
# ===========================================================================

# _vcheck_emit <strict> <quiet> <check> <message>
# Emit one finding. Returns 1 when strict (so the caller can OR it into its
# fail tally), 0 otherwise. Soft WARN and hard ERROR both go to stderr
# (matching _validate_docs_layout's ERROR stream + cmd_validate's failure
# stanza); --quiet suppresses the soft WARN but never a strict ERROR.
_vcheck_emit() {
    local strict="$1" quiet="$2" check="$3" msg="$4"
    if [[ "$strict" -eq 1 ]]; then
        echo "ERROR $check: $msg" >&2
        return 1
    fi
    [[ "$quiet" -eq 0 ]] && echo "WARN $check: $msg" >&2
    return 0
}

# _vcheck_is_exempt <relpath>
# True (0) when a path is in a sanctioned ephemeral zone (_scratch/ anywhere,
# or under $TMPDIR). Such paths never produce findings.
_vcheck_is_exempt() {
    local rel="$1"
    case "/$rel/" in
        */_scratch/*) return 0 ;;
    esac
    if [[ -n "${TMPDIR:-}" ]]; then
        case "$rel" in
            "$TMPDIR"/*|"${TMPDIR%/}"/*) return 0 ;;
        esac
    fi
    return 1
}

# Accepted analysis stage directories, canonical name first. `02_analysis/scripts`
# is the pre-rename spelling: dual acceptance is a temporary migration window so
# lint keeps working in repos that have not renamed yet. Drop the second entry
# once they have; every stage-path decision below reads this one list.
_VCHECK_STAGE_DIRS=(02_analysis/stages 02_analysis/scripts)

# _vcheck_stage_dirs <projdir>
# Echo the absolute path of each accepted stage dir that exists, canonical first.
_vcheck_stage_dirs() {
    local projdir="$1" d
    for d in "${_VCHECK_STAGE_DIRS[@]}"; do
        [[ -d "$projdir/$d" ]] && printf '%s\n' "$projdir/$d"
    done
    return 0
}

# _vcheck_is_stage_path <relpath>
# True (0) when a project-relative path lies under an accepted stage dir.
_vcheck_is_stage_path() {
    local p="$1" d
    for d in "${_VCHECK_STAGE_DIRS[@]}"; do
        [[ "$p" == "$d"/* ]] && return 0
    done
    return 1
}

# _vcheck_config_path <projdir>
# Echo the analysis_config.yaml path (empty string if absent).
_vcheck_config_path() {
    local projdir="$1"
    local cfg="$projdir/02_analysis/config/analysis_config.yaml"
    [[ -f "$cfg" ]] && printf '%s' "$cfg"
}

# _vcheck_config_figures_int <cfg> <key>
# Echo the integer value of `figures.<key>` from the config (empty if absent).
# Reuses the tags-style awk scan: enter the `figures:` block, read `  <key>:`.
_vcheck_config_figures_int() {
    local cfg="$1" key="$2"
    [[ -f "$cfg" ]] || return 0
    awk -v K="$key" '
        /^figures:/ { infig=1; next }
        infig && /^[^ \t]/ { infig=0 }
        infig {
            line=$0
            sub(/[ \t]*#.*$/, "", line)
            if (match(line, "^  " K ":[ \t]*")) {
                val=substr(line, RLENGTH+1)
                gsub(/[ \t]/, "", val)
                if (val ~ /^[0-9]+$/) { print val; exit }
            }
        }
    ' "$cfg"
}

# _vcheck_config_stages <cfg>
# Echo the list of stage ids (one per line) from the `stages:` block.
# Matches `  - id: "<id>"` / `  - id: <id>` (tags-style awk parse).
_vcheck_config_stages() {
    local cfg="$1"
    [[ -f "$cfg" ]] || return 0
    awk '
        /^stages:/ { instages=1; next }
        instages && /^[^ \t-]/ { instages=0 }
        instages && /^[ \t]*-[ \t]*id:/ {
            line=$0
            sub(/^[ \t]*-[ \t]*id:[ \t]*/, "", line)
            sub(/[ \t]*#.*$/, "", line)
            sub(/[ \t]+$/, "", line)
            gsub(/^["'"'"']|["'"'"']$/, "", line)
            print line
        }
    ' "$cfg"
}

# ---------------------------------------------------------------------------
# Check 1: figure-style
# ---------------------------------------------------------------------------
# Scan viz scripts for inline styling that bypasses the project theme entry
# point, plus the config base_size floor. A viz script is any
# 02_analysis/stages/*_viz.{R,py} (or *viz* / a script that calls a plotting
# primitive). Each violation is a WARN (hard-fail under strict).
_validate_check_figure_style() {
    local projdir="$1" strict="$2" quiet="$3"
    local rc=0
    local -a stage_dirs=()
    local _d
    while IFS= read -r _d; do stage_dirs+=("$_d"); done < <(_vcheck_stage_dirs "$projdir")

    # --- config floor: figures.base_size < 14 ------------------------------
    local cfg
    cfg=$(_vcheck_config_path "$projdir")
    if [[ -n "$cfg" ]]; then
        local bs
        bs=$(_vcheck_config_figures_int "$cfg" "base_size")
        if [[ -n "$bs" && "$bs" -lt 14 ]]; then
            _vcheck_emit "$strict" "$quiet" figure-style \
                "analysis_config.yaml:figures.base_size = $bs (< 14 base font floor) — raise to >= 14" || rc=1
        fi
    fi

    [[ ${#stage_dirs[@]} -gt 0 ]] || return $rc

    # Collect candidate viz scripts: *_viz.{R,py}, any *viz* file, or any
    # script that calls a plotting primitive (ggsave/ggplot/savefig/plt.).
    local -a vizscripts=()
    local f rel
    while IFS= read -r f; do
        [[ -f "$f" ]] || continue
        rel="${f#$projdir/}"
        _vcheck_is_exempt "$rel" && continue
        case "$(basename "$f")" in
            *_viz.R|*_viz.py|*viz*) vizscripts+=("$f"); continue ;;
        esac
        if grep -Eq 'ggsave\(|ggplot\(|\.savefig\(|plt\.(plot|figure|subplots)' "$f" 2>/dev/null; then
            vizscripts+=("$f")
        fi
    done < <(find "${stage_dirs[@]}" -maxdepth 3 -type f \( -name '*.R' -o -name '*.py' \) 2>/dev/null)

    for f in "${vizscripts[@]+"${vizscripts[@]}"}"; do
        rel="${f#$projdir/}"
        local hit

        # inline theme(text=...) / element_text(size=...)
        if hit=$(grep -nE 'theme\([^)]*text|element_text\([^)]*size[ \t]*=' "$f" 2>/dev/null | head -n1); then
            [[ -n "$hit" ]] && { _vcheck_emit "$strict" "$quiet" figure-style \
                "$rel: inline theme()/element_text(size=) — style via the project theme entry point (line ${hit%%:*})" || rc=1; }
        fi
        # inline ggsave(width=<number>)
        if hit=$(grep -nE 'ggsave\([^)]*width[ \t]*=[ \t]*[0-9.]+' "$f" 2>/dev/null | head -n1); then
            [[ -n "$hit" ]] && { _vcheck_emit "$strict" "$quiet" figure-style \
                "$rel: inline ggsave(width=<n>) — use save_figure()/save_overview() (line ${hit%%:*})" || rc=1; }
        fi
        # plt.savefig(... figsize ...) — raw figsize geometry
        if hit=$(grep -nE 'savefig\(.*figsize|figsize[ \t]*=' "$f" 2>/dev/null | head -n1); then
            [[ -n "$hit" ]] && { _vcheck_emit "$strict" "$quiet" figure-style \
                "$rel: inline figsize= — use save_figure()/save_overview() (line ${hit%%:*})" || rc=1; }
        fi
        # raw 6-digit hex color literal
        if hit=$(grep -nE '#[0-9a-fA-F]{6}\b' "$f" 2>/dev/null | head -n1); then
            [[ -n "$hit" ]] && { _vcheck_emit "$strict" "$quiet" figure-style \
                "$rel: raw hex color literal — read colors from analysis_config.yaml:colors (line ${hit%%:*})" || rc=1; }
        fi
        # saves a figure but never calls the theme entry point
        if grep -Eq 'ggsave\(|\.savefig\(|save_figure\(|save_overview\(' "$f" 2>/dev/null; then
            if ! grep -Eq 'project_theme\(|set_paper_style\(' "$f" 2>/dev/null; then
                _vcheck_emit "$strict" "$quiet" figure-style \
                    "$rel: saves a figure without calling project_theme()/set_paper_style()" || rc=1
            fi
        fi
    done

    return $rc
}

# ---------------------------------------------------------------------------
# Check 2: results-layout
# ---------------------------------------------------------------------------
# Every artifact under 03_results/ (except objects/, master/, interactive/,
# _scratch/) must sit under <stage>/{figures,tables}/...; <stage> must be a
# known stages: id; no artifact at the 03_results/ root; every figure
# (.png/.pdf) under a figures/ dir needs a same-stem table neighbor under the
# sibling tables/ dir.
_validate_check_results_layout() {
    local projdir="$1" strict="$2" quiet="$3"
    local rc=0
    local results="$projdir/03_results"
    [[ -d "$results" ]] || return 0

    local cfg stages
    cfg=$(_vcheck_config_path "$projdir")
    stages=$(_vcheck_config_stages "$cfg")

    # --- no artifact directly at 03_results/ root --------------------------
    local f base
    while IFS= read -r f; do
        [[ -f "$f" ]] || continue
        base=$(basename "$f")
        [[ "$base" == "README.md" ]] && continue
        _vcheck_emit "$strict" "$quiet" results-layout \
            "artifact at 03_results/ root: $base — move under <stage>/{figures,tables}/" || rc=1
    done < <(find "$results" -maxdepth 1 -type f 2>/dev/null)

    # --- every artifact must be under <stage>/{figures,tables}/ ------------
    # Walk all files; classify by their first two path components under results.
    local rel comp1 comp2 stage kind
    while IFS= read -r f; do
        [[ -f "$f" ]] || continue
        rel="${f#$results/}"
        # Skip sanctioned roots: objects/, master/, interactive/, _scratch/.
        case "$rel" in
            objects/*|master/*|interactive/*|_scratch/*) continue ;;
            */_scratch/*) continue ;;
        esac
        comp1="${rel%%/*}"
        # Files directly at root already handled above.
        [[ "$comp1" == "$rel" ]] && continue
        stage="$comp1"
        local rest="${rel#*/}"
        comp2="${rest%%/*}"
        kind="$comp2"
        # README.md directly under a stage dir is allowed (caption file).
        [[ "$rest" == "README.md" ]] && continue

        # stage must be a known stages: id.
        if [[ -n "$stages" ]] && ! printf '%s\n' "$stages" | grep -qxF "$stage"; then
            _vcheck_emit "$strict" "$quiet" results-layout \
                "unknown stage '$stage' (not in analysis_config.yaml:stages): $rel" || rc=1
        fi
        # kind must be figures/ or tables/.
        if [[ "$kind" != "figures" && "$kind" != "tables" ]]; then
            _vcheck_emit "$strict" "$quiet" results-layout \
                "artifact not under <stage>/{figures,tables}/: $rel" || rc=1
        fi
    done < <(find "$results" -mindepth 2 -type f 2>/dev/null)

    # --- every figure needs a same-stem table neighbor --------------------
    # For each .png/.pdf under a figures/ dir, the sibling tables/ dir (same
    # relative sub-layout) must hold a file with the same stem (figure stems
    # are <name>.<variant>.<ext>; the table neighbor is <name>.<*>).
    local figdir tabdir stem name
    while IFS= read -r f; do
        [[ -f "$f" ]] || continue
        rel="${f#$results/}"
        case "$rel" in */_scratch/*) continue ;; esac
        figdir=$(dirname "$f")
        # sibling tables dir: replace the LAST 'figures' segment with 'tables'.
        tabdir="${figdir/\/figures\//\/tables\/}"
        if [[ "$tabdir" == "$figdir" ]]; then
            # figures is the last segment (no trailing sub-layout).
            tabdir="${figdir%/figures}/tables"
        fi
        base=$(basename "$f")
        # strip variant + extension: foo.screen.png -> foo ; foo.pdf -> foo
        name="${base%.*}"          # drop ext
        name="${name%.screen}"     # drop common variant suffixes
        name="${name%.print}"
        stem="$name"
        # Look for any same-stem file under the sibling tables dir.
        if ! ls "$tabdir/$stem".* >/dev/null 2>&1; then
            _vcheck_emit "$strict" "$quiet" results-layout \
                "figure without same-stem table neighbor: $rel (expected ${tabdir#$projdir/}/$stem.*)" || rc=1
        fi
    done < <(find "$results" -type f \( -name '*.png' -o -name '*.pdf' \) 2>/dev/null | grep '/figures/' )

    return $rc
}

# ---------------------------------------------------------------------------
# Check 3: captions
# ---------------------------------------------------------------------------
# Every artifact file under 03_results/<stage>/{figures,tables}/ must have a
# path-qualified `## <rel-path>` heading in the sibling
# 03_results/<stage>/README.md (the caption format from figure_helpers).
_validate_check_captions() {
    local projdir="$1" strict="$2" quiet="$3"
    local rc=0
    local results="$projdir/03_results"
    [[ -d "$results" ]] || return 0

    local f rel stage rest within readme heading
    while IFS= read -r f; do
        [[ -f "$f" ]] || continue
        rel="${f#$results/}"
        case "$rel" in
            objects/*|master/*|interactive/*|_scratch/*) continue ;;
            */_scratch/*) continue ;;
        esac
        stage="${rel%%/*}"
        rest="${rel#*/}"
        # only figures/ + tables/ artifacts carry captions
        case "$rest" in
            figures/*|tables/*) ;;
            *) continue ;;
        esac
        readme="$results/$stage/README.md"
        # The caption heading is path-qualified RELATIVE to the stage dir.
        within="$rest"
        if [[ ! -f "$readme" ]]; then
            _vcheck_emit "$strict" "$quiet" captions \
                "no README.md in 03_results/$stage/ to caption: $rel" || rc=1
            continue
        fi
        heading="## $within"
        if ! grep -qxF "$heading" "$readme" 2>/dev/null; then
            _vcheck_emit "$strict" "$quiet" captions \
                "no caption section '$heading' in 03_results/$stage/README.md for $rel" || rc=1
        fi
    done < <(find "$results" -mindepth 3 -type f \( -name '*.png' -o -name '*.pdf' -o -name '*.csv' -o -name '*.tsv' \) 2>/dev/null)

    return $rc
}

# ---------------------------------------------------------------------------
# Check 4: provenance
# ---------------------------------------------------------------------------
# For each caption section in a stage README, the `Script:` cell (the first
# cell of the `| Script | Function | Config | Input |` table) must resolve to
# an existing path under 02_analysis/stages/; if the project is a git repo,
# the script must also be tracked.
_validate_check_provenance() {
    local projdir="$1" strict="$2" quiet="$3"
    local rc=0
    local results="$projdir/03_results"
    [[ -d "$results" ]] || return 0

    local in_git
    in_git=$(git -C "$projdir" rev-parse --is-inside-work-tree 2>/dev/null)

    local readme rel script
    while IFS= read -r readme; do
        [[ -f "$readme" ]] || continue
        rel="${readme#$projdir/}"
        case "$rel" in */_scratch/*) continue ;; esac
        # Extract every Script cell from the caption tables. Caption rows look
        # like: | `<script>` | `<fn>` | `<config>` | `<input>` |
        # The Script cell is the first data cell, backtick-wrapped.
        while IFS= read -r script; do
            [[ -z "$script" ]] && continue
            # Resolve relative to project root.
            local abs="$projdir/$script"
            if [[ ! -f "$abs" ]]; then
                _vcheck_emit "$strict" "$quiet" provenance \
                    "${rel}: caption cites missing script: $script" || rc=1
                continue
            fi
            if ! _vcheck_is_stage_path "$script"; then
                _vcheck_emit "$strict" "$quiet" provenance \
                    "${rel}: caption Script not under 02_analysis/stages/: $script" || rc=1
            fi
            # If a git repo, the script must be tracked.
            if [[ "$in_git" == "true" ]]; then
                if ! git -C "$projdir" ls-files --error-unmatch "$script" >/dev/null 2>&1; then
                    _vcheck_emit "$strict" "$quiet" provenance \
                        "${rel}: caption cites untracked script (not committed): $script" || rc=1
                fi
            fi
        done < <(_vcheck_extract_scripts "$readme")
    done < <(find "$results" -mindepth 2 -name 'README.md' -type f 2>/dev/null)

    return $rc
}

# _vcheck_extract_scripts <readme>
# Echo the Script cell (first data cell) of every caption provenance table row.
# Matches a markdown table data row whose first cell is a backtick-wrapped path
# ending in .R/.py, immediately following the `|---|...` separator OR any row
# under a `| Script | Function | ...` header. We simply pull the first
# backtick-wrapped token from rows whose first cell looks like a script path.
_vcheck_extract_scripts() {
    local readme="$1"
    awk '
        # Header row marks the table shape; data rows follow the separator.
        /^\|[ \t]*Script[ \t]*\|/ { inrow=1; next }
        inrow && /^\|[ \t]*-+[ \t]*\|/ { next }     # separator row
        inrow && /^\|/ {
            line=$0
            # first data cell is between the first and second pipe
            n=split(line, cells, "|")
            cell=cells[2]
            gsub(/`/, "", cell)
            gsub(/^[ \t]+|[ \t]+$/, "", cell)
            if (cell != "") print cell
            next
        }
        inrow && !/^\|/ { inrow=0 }
    ' "$readme"
}

# ---------------------------------------------------------------------------
# Check 5: freshness (informational; soft-ish, hard-fail only under --strict)
# ---------------------------------------------------------------------------
# Compare the project's stored CRAFT block hash to what the toolkit's current
# craft.yaml renders (any staleness, not just a version-integer bump). Also, if
# 01_modules/SciAgent-toolkit is a git submodule, compare its checked-out
# commit to the toolkit HEAD. WARN with a `sciagent update` hint when behind.
# Uses block.sh/craft.sh accessors so marker-format knowledge stays in block.sh.
_validate_check_freshness() {
    local projdir="$1" strict="$2" quiet="$3"
    local rc=0

    local tk_root
    if [[ -n "${SCIAGENT_TOOLKIT:-}" ]]; then
        tk_root="$SCIAGENT_TOOLKIT"
    else
        local self_dir
        self_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
        tk_root="$(cd "$self_dir/../.." && pwd)"
    fi

    # --- CRAFT block content vs current craft.yaml render -------------------
    # The repo's CRAFT block is stale when the hash stored in its marker differs
    # from the hash craft.yaml would render to now. block_stored_hash reads the
    # marker; _craft_render_body + block.sh canonicalisation reproduce the hash
    # a fresh render would store. Both functions are sourced for this verb.
    local agents="$projdir/AGENTS.md"
    if [[ -f "$agents" && -f "$tk_root/craft.yaml" ]] \
        && declare -F _craft_render_body >/dev/null && declare -F block_stored_hash >/dev/null; then
        local stored_hash expected_body expected_hash
        stored_hash=$(block_stored_hash "$agents" CRAFT 2>/dev/null)
        if [[ -n "$stored_hash" ]]; then
            expected_body=$(SCIAGENT_TOOLKIT="$tk_root" _craft_render_body 2>/dev/null)
            # Match block_write's trailing-newline canonicalisation before hashing.
            [[ "${expected_body: -1}" == $'\n' ]] || expected_body="${expected_body}"$'\n'
            expected_hash=$(printf '%s' "$expected_body" | _sha1)
            if [[ -n "$expected_hash" && "$stored_hash" != "$expected_hash" ]]; then
                _vcheck_emit "$strict" "$quiet" freshness \
                    "CRAFT block is stale vs toolkit craft.yaml — run: sciagent craft" || rc=1
            fi
        fi
    fi

    # --- submodule commit vs toolkit HEAD ---------------------------------
    local submod="$projdir/01_modules/SciAgent-toolkit"
    if [[ -d "$submod" ]]; then
        local sub_head tk_head
        sub_head=$(git -C "$submod" rev-parse HEAD 2>/dev/null)
        tk_head=$(git -C "$tk_root" rev-parse HEAD 2>/dev/null)
        if [[ -n "$sub_head" && -n "$tk_head" && "$sub_head" != "$tk_head" ]]; then
            # Only warn when the submodule is an ANCESTOR of (behind) HEAD.
            if git -C "$tk_root" merge-base --is-ancestor "$sub_head" "$tk_head" 2>/dev/null; then
                _vcheck_emit "$strict" "$quiet" freshness \
                    "01_modules/SciAgent-toolkit is behind toolkit HEAD (${sub_head:0:8} < ${tk_head:0:8}) — run: sciagent update" || rc=1
            fi
        fi
    fi

    return $rc
}

# ---------------------------------------------------------------------------
# Check 6: hooks
# ---------------------------------------------------------------------------
# Assert that every hook settings.json REGISTERS actually exists and is
# executable. The (c) GUARDRAIL layer is the one tier that binds regardless of
# model cooperation, so a registered-but-absent hook is a silent hole: Claude
# Code fails the hook call quietly and the convention stops being enforced with
# no signal anywhere.
#
# This check exists because that state shipped. Meta-Aging/14616-DM registered
# PreToolUse -> .claude/hooks/no_ephemeral.sh and Stop -> caption_sweep.sh with
# no .claude/hooks/ directory at all, because hook bodies were rendered only by
# `new project` while `activate` merged in the registering settings template.
_validate_check_hooks() {
    local projdir="$1" strict="$2" quiet="$3"
    local rc=0

    local settings="$projdir/.claude/settings.json"
    [[ -f "$settings" ]] || return 0   # no settings: nothing registered, nothing to verify

    if ! command -v jq >/dev/null 2>&1; then
        _vcheck_emit "$strict" "$quiet" hooks \
            "jq not found — cannot verify registered hooks in .claude/settings.json" || rc=1
        return $rc
    fi

    # Every registered hook command, one per line. Absent/!object `hooks` yields
    # empty output rather than a jq error.
    local cmds
    cmds=$(jq -r '(.hooks // {}) | to_entries[]? | .value[]? | .hooks[]? | .command // empty' \
        "$settings" 2>/dev/null)
    [[ -n "$cmds" ]] || return 0

    local cmd path
    while IFS= read -r cmd; do
        [[ -n "$cmd" ]] || continue
        # Pull the first .claude/hooks/<file> token out of the command line. Hook
        # commands are `bash "$CLAUDE_PROJECT_DIR/.claude/hooks/<name>.sh"`, so
        # match on the project-relative tail and ignore the interpreter/quoting.
        path=$(printf '%s\n' "$cmd" | grep -oE '\.claude/hooks/[A-Za-z0-9._-]+' | head -1)
        [[ -n "$path" ]] || continue   # not a project-local hook script; not ours to verify

        if [[ ! -f "$projdir/$path" ]]; then
            _vcheck_emit "$strict" "$quiet" hooks \
                "$path is registered in .claude/settings.json but does not exist — run: sciagent activate" || rc=1
        elif [[ ! -x "$projdir/$path" ]] && ! [[ "$cmd" =~ (^|[[:space:]/])(bash|sh|zsh|python3?|uv)([[:space:]]|$) ]]; then
            # The exec bit only matters when the hook is invoked DIRECTLY. The
            # shipped templates register `bash "<path>"`, where mode 0644 runs
            # fine — demanding +x there would warn on a non-problem in every
            # correctly-provisioned repo.
            _vcheck_emit "$strict" "$quiet" hooks \
                "$path is registered for direct execution but is not executable — run: chmod +x $path" || rc=1
        fi
    done <<< "$cmds"

    return $rc
}

# _validate_run_project_checks <projdir> <strict> <quiet> <check>...
# Dispatch the named opt-in project checks. Accepts `all` to run every check.
# Returns 1 if any check reports a (strict) hard failure, 0 otherwise.
_validate_run_project_checks() {
    local projdir="$1"; shift
    local strict="$1"; shift
    local quiet="$1"; shift
    local -a names=("$@")
    local rc=0

    # Expand `all`.
    local -a run=()
    local n
    for n in "${names[@]}"; do
        case "$n" in
            all) run=(figure-style results-layout captions provenance freshness hooks); break ;;
            figure-style|results-layout|captions|provenance|freshness|hooks) run+=("$n") ;;
            "") ;;
            *)
                echo "sciagent validate: unknown --check name '$n'" >&2
                echo "  valid: figure-style results-layout captions provenance freshness hooks all" >&2
                return 1 ;;
        esac
    done

    for n in "${run[@]+"${run[@]}"}"; do
        case "$n" in
            figure-style)   _validate_check_figure_style   "$projdir" "$strict" "$quiet" || rc=1 ;;
            results-layout) _validate_check_results_layout "$projdir" "$strict" "$quiet" || rc=1 ;;
            captions)       _validate_check_captions       "$projdir" "$strict" "$quiet" || rc=1 ;;
            provenance)     _validate_check_provenance     "$projdir" "$strict" "$quiet" || rc=1 ;;
            freshness)      _validate_check_freshness      "$projdir" "$strict" "$quiet" || rc=1 ;;
            hooks)          _validate_check_hooks          "$projdir" "$strict" "$quiet" || rc=1 ;;
        esac
    done

    return $rc
}

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
  --check <name>     Run an opt-in PROJECT guardrail check against --project-dir.
                     name ∈ figure-style | results-layout | captions |
                            provenance | freshness | all. Repeatable.
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
        _validate_run_project_checks "$_projdir" "$_strict" "$quiet" "${_checks[@]}"
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
