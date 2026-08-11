# lib/sciagent/lint.sh — sciagent lint [--project-dir <dir>] [--check <name>...] [--strict] [--quiet]
#
# The (c) GUARDRAIL layer: opt-in PROJECT guardrail checks. These are aimed at
# analysis-repo authors (vs. validate.sh's toolkit-maintainer audience: skill
# frontmatter shape, cross-namespace collisions). Different audience, different
# failure semantics — warn-only by default, hard-fail only under --strict.
#
# Checks:
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
#   stage-thinness  see docs 09 §3.N (stub; implemented in a parallel change).
#   comment-intent  see docs 09 §3.N (stub; implemented in a parallel change).
#   stage-layout    see docs 09 §3.N (stub; implemented in a parallel change).
# These run ONLY when --check <name> (or --check all, or no --check at all —
# `all` is the default) is given. They are SOFT warnings (exit 0) by default
# and HARD failures (exit 1) under --strict. `_scratch/` and $TMPDIR are always
# exempt.
#
# Hardness boundary:
#   Hard-fail (exit 1): any finding when --strict.
#   Soft-warn:          findings without --strict.
#
# Exit code:
#   0 — no --strict hard failures
#   1 — a --strict hard failure, or an unknown --check name
#
# Output streams:
#   stdout — nothing on success (this verb is quiet by design; there is no
#            toolkit-wide walk to summarize)
#   stderr — WARN/ERROR per finding
#
# Depends on block.sh + craft.sh (the freshness check's CRAFT-hash comparison
# reuses their accessors unchanged).

# shellcheck shell=bash

# ===========================================================================
# PROJECT GUARDRAIL CHECKS (opt-in via --check). The (c) GUARDRAIL layer.
#
# Each check is `_lint_check_<name> <projdir> <strict> <quiet>` and emits
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
# fail tally), 0 otherwise. Soft WARN and hard ERROR both go to stderr;
# --quiet suppresses the soft WARN but never a strict ERROR.
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
_lint_check_figure_style() {
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
_lint_check_results_layout() {
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
_lint_check_captions() {
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
_lint_check_provenance() {
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
_lint_check_freshness() {
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
_lint_check_hooks() {
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

# _lint_run_checks <projdir> <strict> <quiet> <check>...
# Dispatch the named opt-in project checks. Accepts `all` to run every check.
# Returns 1 if any check reports a (strict) hard failure, 0 otherwise.
_lint_run_checks() {
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
            all) run=(figure-style results-layout captions provenance freshness hooks stage-thinness comment-intent stage-layout); break ;;
            figure-style|results-layout|captions|provenance|freshness|hooks|stage-thinness|comment-intent|stage-layout) run+=("$n") ;;
            "") ;;
            *)
                echo "sciagent lint: unknown --check name '$n'" >&2
                echo "  valid: figure-style results-layout captions provenance freshness hooks stage-thinness comment-intent stage-layout all" >&2
                return 1 ;;
        esac
    done

    for n in "${run[@]+"${run[@]}"}"; do
        case "$n" in
            figure-style)   _lint_check_figure_style   "$projdir" "$strict" "$quiet" || rc=1 ;;
            results-layout) _lint_check_results_layout "$projdir" "$strict" "$quiet" || rc=1 ;;
            captions)       _lint_check_captions       "$projdir" "$strict" "$quiet" || rc=1 ;;
            provenance)     _lint_check_provenance     "$projdir" "$strict" "$quiet" || rc=1 ;;
            freshness)      _lint_check_freshness      "$projdir" "$strict" "$quiet" || rc=1 ;;
            hooks)          _lint_check_hooks          "$projdir" "$strict" "$quiet" || rc=1 ;;
            stage-thinness) _lint_check_stage_thinness "$projdir" "$strict" "$quiet" || rc=1 ;;
            comment-intent) _lint_check_comment_intent "$projdir" "$strict" "$quiet" || rc=1 ;;
            stage-layout)   _lint_check_stage_layout   "$projdir" "$strict" "$quiet" || rc=1 ;;
        esac
    done

    return $rc
}

# cmd_lint [--project-dir <dir>] [--check <name>...] [--strict] [--quiet]
cmd_lint() {
    local quiet=0
    local _projdir="."
    local _strict=0
    local -a _checks=()
    while [[ $# -gt 0 ]]; do
        case "$1" in
            -h|--help)
                cat <<'USAGE'
sciagent lint [--project-dir <dir>] [--check <name>...] [--strict] [--quiet]
  Run the opt-in PROJECT guardrail checks (the (c) GUARDRAIL layer) against
  --project-dir. Aimed at analysis-repo authors, not the toolkit itself — see
  `sciagent validate` for the toolkit-wide skill-frontmatter walk.

  --project-dir <d>  Project directory to lint (default: .).
  --check <name>     Run only the named check(s). Repeatable.
                     name ∈ figure-style | results-layout | captions |
                            provenance | freshness | hooks | stage-thinness |
                            comment-intent | stage-layout | all.
                     Default (no --check given): all.
  --strict           Findings become HARD failures (exit 1) instead of WARN.
  --quiet            Suppress soft-WARN output on success/no-strict findings.

  _scratch/ and $TMPDIR are always exempt from every check.

Exit code: 0 unless --strict promotes a finding, or an unknown --check name
was given (exit 1 either way).
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
                    echo "sciagent lint: unknown option '$1'" >&2
                    echo "usage: sciagent lint [--project-dir <dir>] [--check <name>] [--strict] [--quiet]" >&2
                    return 1
                fi ;;
        esac
    done

    # Default when no --check given: run `all`.
    [[ ${#_checks[@]} -gt 0 ]] || _checks=(all)

    _lint_run_checks "$_projdir" "$_strict" "$quiet" "${_checks[@]}"
    return $?
}

# ---------------------------------------------------------------------------
# Check: stage-thinness
# ---------------------------------------------------------------------------
# Implemented in Phase 5; see docs 09 §3.N.
_lint_check_stage_thinness() {
    local projdir="$1" strict="$2" quiet="$3"
    # Implemented in Phase 5; see docs 09 §3.N.
    return 0
}

# ---------------------------------------------------------------------------
# Check: comment-intent
# ---------------------------------------------------------------------------
# Implemented in Phase 5; see docs 09 §3.N.
_lint_check_comment_intent() {
    local projdir="$1" strict="$2" quiet="$3"
    # Implemented in Phase 5; see docs 09 §3.N.
    return 0
}

# ---------------------------------------------------------------------------
# Check: stage-layout
# ---------------------------------------------------------------------------
# Implemented in Phase 5; see docs 09 §3.N.
_lint_check_stage_layout() {
    local projdir="$1" strict="$2" quiet="$3"
    # Implemented in Phase 5; see docs 09 §3.N.
    return 0
}
