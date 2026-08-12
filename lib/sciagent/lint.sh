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
#   docs-layout     docs/_internal/ must be gitignored when the project is a
#                   git repo (leaks internal notes on push otherwise); plus
#                   softer structural warnings — missing docs/_internal/, a
#                   stray .md at the 03_results/ root, mixed archive-naming
#                   conventions, non-standard handoff filenames. No-ops
#                   entirely (same absent-subject pattern as figure-style/
#                   results-layout/captions/provenance) when the project has
#                   no docs/ tree at all — it audits the STRUCTURE of an
#                   existing docs/ tree, it does not mandate that one exist.
#                   Moved here from validate.sh, where it used to run
#                   unconditionally on the default path and could hard-block
#                   `activate`'s pre-flight — see validate.sh's header
#                   comment for the incident this fixed.
#   stage-thinness  02_analysis/stages|scripts function-def count (>2),
#                   summed function-body lines (>60), total LOC (>500) per
#                   file; see docs 09 §3.1. `# stage-detail: <reason>` exempts
#                   one definition.
#   comment-intent  see docs 09 §3.N (stub; implemented in a parallel change).
#   stage-layout    see docs 09 §3.N (stub; implemented in a parallel change).
#   skill-coupling  drift guard on skills' `compatibility:` declarations. THE
#                   ODD ONE OUT: its subject is the TOOLKIT checkout, not
#                   --project-dir, and it is warn-only even under --strict.
#                   Opt-in BY NAME ONLY — deliberately not a member of `all`.
#                   See its own header block for both decisions.
# These run ONLY when --check <name> (or --check all, or no --check at all —
# `all` is the default) is given. They are SOFT warnings (exit 0) by default
# and HARD failures (exit 1) under --strict (skill-coupling excepted — it can
# never hard-fail). `_scratch/` and $TMPDIR are always exempt.
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

# _vcheck_stage_files <projdir>
# Emit every stage SCRIPT file, one per line, across all spellings in
# _VCHECK_STAGE_DIRS.
#
# Single definition on purpose. The doc-09 §3 checks each grew their own `find`
# and ended up disagreeing about which files are stages: depth 1, 2 and 3 with
# two different extension sets. In a repo using the nested
# 02_analysis/scripts/{compute,figures}/ layout that meant one check saw a file
# and its sibling check did not, so a stage could pass stage-layout while
# stage-thinness had never looked at it. Divergence between checks on "what is
# a stage" is a defect, not a tuning knob.
#
# Depth 3 covers both the flat doc-09 layout (stages/NN_topic.R) and the one
# level of nesting the pre-migration fleet still uses. .sh is included: an
# ordered stage may legitimately be a shell script.
_vcheck_stage_files() {
    local projdir="$1"
    local -a dirs=()
    local d
    while IFS= read -r d; do dirs+=("$d"); done < <(_vcheck_stage_dirs "$projdir")
    [[ ${#dirs[@]} -gt 0 ]] || return 0
    find "${dirs[@]}" -maxdepth 3 -type f \
        \( -name '*.R' -o -name '*.py' -o -name '*.sh' \) 2>/dev/null | sort
    return 0
}

# _vcheck_toolkit_root
# Echo the toolkit checkout this lint run belongs to: $SCIAGENT_TOOLKIT when the
# dispatcher exported one (always, in practice — bin/sciagent sets it before any
# module loads), else the checkout this very file lives in. Two checks need it
# for different reasons — freshness compares a project against it,
# skill-coupling audits it directly — so the resolution lives in one place.
_vcheck_toolkit_root() {
    if [[ -n "${SCIAGENT_TOOLKIT:-}" ]]; then
        printf '%s' "$SCIAGENT_TOOLKIT"
        return 0
    fi
    local self_dir
    self_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    (cd "$self_dir/../.." && pwd)
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
    tk_root="$(_vcheck_toolkit_root)"

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
# Check: docs-layout
# ---------------------------------------------------------------------------
# Moved from validate.sh (formerly _validate_docs_layout, called
# unconditionally on validate's default path). Findings, all soft-warn unless
# --strict (including the not-gitignored rule below — previously an
# unconditional hard-fail regardless of --strict/--quiet; now consistent with
# every other lint check's hardness boundary, since this check is opt-in and
# no longer reachable from activate's pre-flight):
#   - docs/_internal/ does not exist (docs/ itself does — see the absent-
#     subject guard below)
#   - docs/_internal/ exists, is inside a git repo, and is NOT gitignored
#     (would expose internal notes on a push)
#   - a .md file directly in 03_results/ (maxdepth 1)
#   - multiple archive-naming conventions (.archive/_deprecated/_legacy/
#     .deprecated) coexisting under the same parent
#   - a handoff_*.md file at the project root not matching
#     handoff_YYYYMMDD_HHMMSS.md
_lint_check_docs_layout() {
    local projdir="$1" strict="$2" quiet="$3"
    local rc=0

    # Absent-subject guard, same shape as every sibling check in this file
    # (figure-style's stage_dirs check; results-layout's/captions'/
    # provenance's `[[ -d "$results" ]] || return 0`): no docs/ tree at all
    # means this project hasn't adopted the docs/ convention (or is a
    # software/early-stage repo that doesn't need one) — that is an absent
    # subject, not a finding, and must produce NO output, not a WARN. Per
    # this file's header (and the standing project position that these
    # thresholds are aspirational, a backlog rather than a baseline), a
    # project that simply hasn't grown a docs/ tree yet is not "wrong" by
    # default. docs-layout audits the STRUCTURE of an existing docs/ tree; it
    # does not mandate that one exist.
    [[ -d "$projdir/docs" ]] || return 0

    # docs/_internal/ does not exist (docs/ itself does, so the project has
    # opted into the convention but hasn't finished scaffolding it).
    [[ -d "$projdir/docs/_internal" ]] || \
        { _vcheck_emit "$strict" "$quiet" docs-layout "docs/_internal/ missing — run: sciagent gitignore" || rc=1; }

    # docs/_internal/ exists but is NOT gitignored (in a git repo).
    if [[ -d "$projdir/docs/_internal" ]]; then
        local _in_git
        _in_git=$(git -C "$projdir" rev-parse --is-inside-work-tree 2>/dev/null)
        if [[ "$_in_git" == "true" ]] && ! git -C "$projdir" check-ignore -q docs/_internal 2>/dev/null; then
            _vcheck_emit "$strict" "$quiet" docs-layout \
                "docs/_internal/ is NOT gitignored — add 'docs/_internal/' to .gitignore" || rc=1
        fi
    fi

    # Any .md file directly in 03_results/ at maxdepth 1.
    local f
    while IFS= read -r f; do
        [[ -f "$f" ]] || continue
        _vcheck_emit "$strict" "$quiet" docs-layout \
            "report in results dir: $(basename "$f") — move to docs/_internal/reports/" || rc=1
    done < <(find "$projdir/03_results" -maxdepth 1 -name "*.md" 2>/dev/null)

    # Multiple archive-style naming conventions coexisting under one parent.
    local _dir
    for _dir in "$projdir/02_analysis" "$projdir/03_results"; do
        [[ -d "$_dir" ]] || continue
        local _conventions
        _conventions=$(find "$_dir" -maxdepth 3 -type d \
            \( -name ".archive" -o -name "_deprecated" -o -name "_legacy" -o -name ".deprecated" \) \
            2>/dev/null | awk -F'/' '{
                n=split($0,a,"/")
                parent=""
                for(i=1;i<n;i++) parent=parent (i>1?"/":"") a[i]
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
                _vcheck_emit "$strict" "$quiet" docs-layout \
                    "mixed archive-naming conventions under $_p — pick one (.archive / _deprecated / _legacy / .deprecated)" || rc=1
            done <<< "$_conventions"
        fi
    done

    # handoff_*.md files at the project root that don't match
    # handoff_YYYYMMDD_HHMMSS.md.
    for f in "$projdir"/handoff_*.md; do
        [[ -f "$f" ]] || continue
        if ! printf '%s\n' "$(basename "$f")" | grep -qE '^handoff_[0-9]{8}_[0-9]{6}\.md$'; then
            _vcheck_emit "$strict" "$quiet" docs-layout \
                "non-standard handoff filename: $(basename "$f")" || rc=1
        fi
    done

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
# `docs-layout` upholds the same "absent subject -> no findings" invariant
# every other check here follows (see its own header comment above), so it is
# a full member of `all` like everything else.
#
# `skill-coupling` is the ONE check `all` does not expand to, and the only one
# that is runnable but unlisted there: its subject is the toolkit checkout
# rather than <projdir>, so it is opt-in by name only. See its header block.
#
# Returns 1 if any check reports a (strict) hard failure, 0 otherwise.
# `skill-coupling` can never contribute to that tally — it always returns 0.
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
            all) run=(figure-style results-layout captions provenance freshness hooks docs-layout stage-thinness comment-intent stage-layout); break ;;
            figure-style|results-layout|captions|provenance|freshness|hooks|docs-layout|stage-thinness|comment-intent|stage-layout|skill-coupling) run+=("$n") ;;
            "") ;;
            *)
                echo "sciagent lint: unknown --check name '$n'" >&2
                echo "  valid: figure-style results-layout captions provenance freshness hooks docs-layout stage-thinness comment-intent stage-layout skill-coupling all" >&2
                echo "  ('all' covers every name above except skill-coupling, whose subject is the toolkit, not a project)" >&2
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
            docs-layout)    _lint_check_docs_layout    "$projdir" "$strict" "$quiet" || rc=1 ;;
            stage-thinness) _lint_check_stage_thinness "$projdir" "$strict" "$quiet" || rc=1 ;;
            comment-intent) _lint_check_comment_intent "$projdir" "$strict" "$quiet" || rc=1 ;;
            stage-layout)   _lint_check_stage_layout   "$projdir" "$strict" "$quiet" || rc=1 ;;
            skill-coupling) _lint_check_skill_coupling "$projdir" "$strict" "$quiet" || rc=1 ;;
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
                            provenance | freshness | hooks | docs-layout |
                            stage-thinness | comment-intent | stage-layout |
                            skill-coupling | all.
                     Default (no --check given): all.
  --strict           Findings become HARD failures (exit 1) instead of WARN.
  --quiet            Suppress soft-WARN output on success/no-strict findings.

  skill-coupling is the one check `all` does not include, and the one that
  ignores --project-dir: it audits the TOOLKIT checkout's skills/ for drift
  between a skill's `compatibility:` declaration and its actual content
  (undeclared coupling, stale declarations, scaffold paths the toolkit does
  not scaffold). Heuristic by nature, so it is warn-only even under --strict
  and can never change this verb's exit code. Run it explicitly:
      sciagent lint --check skill-coupling

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
# See docs 09 §3.1. Three rules, evaluated per stage file under an accepted
# stage dir (_VCHECK_STAGE_DIRS — both `stages/` and `scripts/` spellings, per
# the migration window; see the comment above that array). `helpers/` is out
# of scope by construction: it is never one of _VCHECK_STAGE_DIRS, so a
# function living there never reaches this check at all — that is where
# definitions belong, per §1.2/§2.1.
#
#   1. function definitions in one stage file            > 2
#   2. total lines inside (non-exempt) function bodies    > 60
#      in one stage file (summed across its definitions)
#   3. total stage file length                            > 500 LOC
#
# Per §0, thresholds are aspirational (a backlog, not a baseline): this check
# ships warn-only like its siblings, hard-fail only under --strict.
#
# Escape hatch (§3.1): a `# stage-detail: <reason>` comment on the line
# immediately above a definition exempts THAT definition from both the
# def-count and the body-line tally. The reason must be non-empty — a bare
# `# stage-detail:` does not exempt, since the reason text is itself the
# required intent documentation.
#
# Detection, and its known blind spots (see _vcheck_stage_scan_defs):
#   R      `name <- function(...)`, `name = function(...)`, and the `\(x)`
#          lambda shorthand when bound to a name (`f <- \(x) ...`). Anonymous
#          `function(...)` passed inline (e.g. to `lapply`) is NOT counted —
#          only named top-level bindings read as "a definition" to a human.
#   Python `def name(...)` / `async def name(...)` at ANY indent level — this
#          deliberately also catches methods/nested defs inside a class, per
#          spec.
#   sh     `name() { ... }` / `function name { ... }`.
#
#   Blind spots, documented rather than silently under-reported:
#   - R/sh: multiple defs on one line are only matched once per line (the
#     regex anchors at line start) — vanishingly rare in the corpus.
#   - R/sh body-line counting is brace-balance from the first `{` after the
#     def line to its match; a brace-less one-liner (`f <- \(x) x + 1`) has
#     no braces to balance and is approximated as a 1-line body.
#   - Python body-line counting is indentation-based: every subsequent line
#     indented strictly deeper than the `def` line is body, ending at the
#     first line at or below that indent. Blank lines strictly between two
#     body statements are counted as body; but trailing blank lines that
#     separate a function from the NEXT construct at shallower indent get
#     attributed to the preceding function's body too (indentation alone
#     can't distinguish "blank line inside" from "blank line after") — a
#     documented over-count, not a silent one.
#   - A `def`/`function(` token that only *looks* like one inside a string
#     literal or docstring is not special-cased; ordinary text lines rarely
#     start with `def `/`name <- function(` so this is low-incidence.
_lint_check_stage_thinness() {
    local projdir="$1" strict="$2" quiet="$3"
    local rc=0
    local -a stage_dirs=()
    local _d
    while IFS= read -r _d; do stage_dirs+=("$_d"); done < <(_vcheck_stage_dirs "$projdir")
    [[ ${#stage_dirs[@]} -gt 0 ]] || return 0

    local f rel lang loc scan defcount bodytotal
    while IFS= read -r f; do
        [[ -f "$f" ]] || continue
        rel="${f#$projdir/}"
        _vcheck_is_exempt "$rel" && continue
        case "$f" in
            *.R) lang=R ;;
            *.py) lang=py ;;
            *.sh) lang=sh ;;
            *) continue ;;
        esac

        loc=$(wc -l < "$f" 2>/dev/null); loc=${loc:-0}

        scan=$(_vcheck_stage_scan_defs "$f" "$lang")
        defcount=$(printf '%s\n' "$scan" | awk -F: '$1=="DEFCOUNT"{print $2}')
        bodytotal=$(printf '%s\n' "$scan" | awk -F: '$1=="BODYTOTAL"{print $2}')
        defcount=${defcount:-0}
        bodytotal=${bodytotal:-0}

        if [[ "$defcount" -gt 2 ]]; then
            _vcheck_emit "$strict" "$quiet" stage-thinness \
                "$rel: $defcount function definitions (> 2) — carve machinery into helpers/" || rc=1
        fi
        if [[ "$bodytotal" -gt 60 ]]; then
            _vcheck_emit "$strict" "$quiet" stage-thinness \
                "$rel: $bodytotal lines inside function bodies (> 60) — carve machinery into helpers/" || rc=1
        fi
        if [[ "$loc" -gt 500 ]]; then
            _vcheck_emit "$strict" "$quiet" stage-thinness \
                "$rel: $loc lines (> 500 LOC) — split the stage or extract helpers" || rc=1
        fi
    done < <(_vcheck_stage_files "$projdir")

    return $rc
}

# _vcheck_stage_scan_defs <file> <lang>
# Print `DEFCOUNT:<n>` and `BODYTOTAL:<n>` — the count of non-exempt function
# definitions and the sum of their body-line counts — for <file>. <lang> is
# one of R|py|sh. See the stage-thinness header comment above for the
# per-language detection rules, the escape-hatch handling, and known blind
# spots.
_vcheck_stage_scan_defs() {
    local file="$1" lang="$2"
    awk -v lang="$lang" '
    {
        lines[NR] = $0
    }
    END {
        n = NR
        defcount = 0
        bodytotal = 0
        for (i = 1; i <= n; i++) {
            line = lines[i]
            isdef = 0
            if (lang == "R") {
                if (line ~ /^[ \t]*[A-Za-z_.][A-Za-z0-9_.]*[ \t]*(<-|=)[ \t]*function[ \t]*\(/) isdef = 1
                else if (line ~ /^[ \t]*[A-Za-z_.][A-Za-z0-9_.]*[ \t]*(<-|=)[ \t]*\\\(/) isdef = 1
            } else if (lang == "py") {
                if (line ~ /^[ \t]*(async[ \t]+)?def[ \t]+[A-Za-z_][A-Za-z0-9_]*[ \t]*\(/) isdef = 1
            } else if (lang == "sh") {
                if (line ~ /^[ \t]*(function[ \t]+)?[A-Za-z_][A-Za-z0-9_]*[ \t]*\(\)/) isdef = 1
                else if (line ~ /^[ \t]*function[ \t]+[A-Za-z_][A-Za-z0-9_]*[ \t]*\{?[ \t]*$/) isdef = 1
            }
            if (!isdef) continue

            exempt = 0
            if (i > 1) {
                prev = lines[i-1]
                if (match(prev, /^[ \t]*#[ \t]*stage-detail:[ \t]*[^ \t]/)) exempt = 1
            }

            blen = 0
            if (lang == "py") {
                defindent = match(line, /[^ \t]/)
                if (defindent == 0) defindent = 1
                j = i + 1
                bodyend = i
                while (j <= n) {
                    bl = lines[j]
                    if (bl ~ /^[ \t]*$/) { bodyend = j; j++; continue }
                    ind = match(bl, /[^ \t]/)
                    if (ind == 0) ind = 1
                    if (ind <= defindent) break
                    bodyend = j
                    j++
                }
                blen = bodyend - i
                if (blen < 0) blen = 0
            } else {
                depth = 0
                started = 0
                bodystart = i
                bodyend = i
                k = i
                done = 0
                while (k <= n && !done) {
                    cl = lines[k]
                    # Blank out simple quoted-string contents before counting
                    # braces: R/glue()/f-string-style interpolation routinely
                    # embeds literal { }  inside "..."/(quotes (e.g. glue(
                    # "id: {tag}")), which would otherwise desync the brace
                    # balance. This is a line-local, escape-naive strip (does
                    # not track escaped quotes across lines) — a documented
                    # approximation, not a full tokenizer.
                    gsub(/"[^"]*"/, "\"\"", cl)
                    gsub(/\047[^\047]*\047/, "\047\047", cl)
                    clen = length(cl)
                    for (cpos = 1; cpos <= clen; cpos++) {
                        ch = substr(cl, cpos, 1)
                        if (ch == "{") {
                            if (!started) { started = 1; bodystart = k }
                            depth++
                        } else if (ch == "}") {
                            depth--
                            if (started && depth == 0) { bodyend = k; done = 1; break }
                        }
                    }
                    k++
                }
                if (started) {
                    blen = bodyend - bodystart
                } else {
                    blen = 1
                }
            }

            if (!exempt) {
                defcount++
                bodytotal += blen
            }
        }
        print "DEFCOUNT:" defcount
        print "BODYTOTAL:" bodytotal
    }
    ' "$file"
}

# ---------------------------------------------------------------------------
# Check: comment-intent
# ---------------------------------------------------------------------------
# See docs 09 §3.2. Exactly two rules — the rejected candidates (history
# markers, docstring-beyond-contract, comment-restates-next-line) are recorded
# in §3.4 and deliberately NOT implemented here.
#
#   1. Banner / separator comment: `^#\s*[=\-*#]{6,}` (any line, any file).
#   2. Mid-file comment run of >= 6 CONSECUTIVE comment lines, where the run
#      STARTS after line 25 (§3.2's literal wording — "starting after line
#      25", not merely overlapping it). This is what keeps a file HEADER
#      legal: a stage should still open with a comment block stating what it
#      does and what it writes (§3.2, §1.4). A run that starts at or before
#      line 25 is exempt in full, even if it continues past line 25 — the
#      whole run is "the header", not just its first 25 lines.
#
# Scope decision (documented per the task): STAGE FILES ONLY, via
# _vcheck_stage_dirs/_vcheck_is_stage_path — the same accepted-dirs list
# stage-thinness uses (both `stages/` and `scripts/` spellings; see the
# comment on _VCHECK_STAGE_DIRS). helpers/ is NOT scanned. Justification:
# §3.1 exempts helpers/ from stage-thinness "by definition" because that is
# where machinery/definitions belong; the comment-intent rules exist to keep
# a STAGE's narrative readable top-to-bottom (§1.4, §2.1 "a stage is
# narrative... implementation lives in helpers/, not inline") — the CRAFT
# comment bullet (§2.2) is written against that same narrative-readability
# rationale ("if intent needs six lines mid-body, the code needs a helper
# with a name"), which is a statement about STAGES, not about helpers'
# already-unordered machinery. Scoping to stages/ keeps this check
# consistent with its siblings stage-thinness/stage-layout, both of which are
# stage-only per §3.1/§3.3.
#
# Blank-line handling (documented): a blank line BREAKS a comment run. §3.2
# says "consecutive comment lines"; a blank line is not a comment line, so it
# ends the run. This also naturally separates a header block (lines 1..~20)
# from an unrelated mid-file block below it, without extra bookkeeping.
#
# False-positive sources, documented rather than special-cased away:
#   - A line that is `#` inside a string/docstring literal (not a real
#     comment) is not distinguished from a real comment — same blind spot
#     _vcheck_stage_scan_defs already accepts for `def`/`function(` tokens.
#   - A Python shebang `#!/usr/bin/env python3` on line 1 matches neither
#     rule in practice (it isn't 6+ fill chars, and a lone line 1 can't reach
#     a 6-line run by itself), so no special-case is needed; documented here
#     because the task explicitly flagged it as a thing to consider.
#   - A legitimate `# -----` section divider inside a long helper is out of
#     scope entirely under the stages-only decision above, so it is not a
#     false positive here (helpers/ is never scanned).
_lint_check_comment_intent() {
    local projdir="$1" strict="$2" quiet="$3"
    local rc=0
    local -a stage_dirs=()
    local _d
    while IFS= read -r _d; do stage_dirs+=("$_d"); done < <(_vcheck_stage_dirs "$projdir")
    [[ ${#stage_dirs[@]} -gt 0 ]] || return 0

    local f rel
    while IFS= read -r f; do
        [[ -f "$f" ]] || continue
        rel="${f#$projdir/}"
        _vcheck_is_exempt "$rel" && continue
        case "$f" in
            *.R|*.py|*.sh) ;;
            *) continue ;;
        esac

        # --- rule 1: banner / separator comments -------------------------
        local hit
        while IFS= read -r hit; do
            [[ -n "$hit" ]] || continue
            _vcheck_emit "$strict" "$quiet" comment-intent \
                "$rel:${hit%%:*}: banner/separator comment — drop it, per craft.yaml's comments bullet" || rc=1
        done < <(grep -nE '^#\s*[=*#-]{6,}' "$f" 2>/dev/null)

        # --- rule 2: mid-file comment run >= 6 lines starting after line 25
        local run_start=0 run_len=0 lineno=0 line
        while IFS= read -r line || [[ -n "$line" ]]; do
            lineno=$((lineno + 1))
            if [[ "$line" =~ ^[[:space:]]*#.*$ && -n "${line//[[:space:]]/}" ]]; then
                [[ "$run_len" -eq 0 ]] && run_start=$lineno
                run_len=$((run_len + 1))
            else
                if [[ "$run_len" -ge 6 && "$run_start" -gt 25 ]]; then
                    _vcheck_emit "$strict" "$quiet" comment-intent \
                        "$rel:$run_start: mid-file comment run of $run_len lines — extract a named helper instead" || rc=1
                fi
                run_len=0
            fi
        done < "$f"
        # File-end flush: a run may run to EOF without a trailing non-comment line.
        if [[ "$run_len" -ge 6 && "$run_start" -gt 25 ]]; then
            _vcheck_emit "$strict" "$quiet" comment-intent \
                "$rel:$run_start: mid-file comment run of $run_len lines — extract a named helper instead" || rc=1
        fi
    done < <(_vcheck_stage_files "$projdir")

    return $rc
}

# ---------------------------------------------------------------------------
# Check: stage-layout
# ---------------------------------------------------------------------------
# Cheap structural companion to doc 09 §1 / §3.3. Six rules, all emitted via
# _vcheck_emit (which itself collapses "violation" vs "warn" — see the doc
# comment above _vcheck_emit / the header block: there is only soft-WARN vs
# --strict hard-ERROR, no third severity channel, so §3.3's "violation" and
# "warn" language both land as the same finding class here):
#   1. a *_viz stage that writes tables         -> finding
#   2. a non-_viz stage that writes figures     -> finding
#   3. a NN_<topic>_viz with no NN_<topic> sibling (same number, same stem)
#                                                -> finding (orphan viz)
#   4. a stray file directly under 02_analysis/ that is neither in a subdir
#      nor a known entry point                  -> finding
#   5. a stage number with a letter suffix (03d, 06b)      -> finding
#   6. a single-digit stage number (4_... not 04_...)      -> finding
#
# Rules 1/2 reuse the exact figure-write regex from figure-style
# (ggsave/.savefig/save_figure/save_overview) so the two checks never diverge
# on what "writes a figure" means; there is no existing "writes a table"
# notion anywhere in lint.sh (results-layout reasons from artifacts already on
# disk, not from script source), so that regex is new here.
_lint_check_stage_layout() {
    local projdir="$1" strict="$2" quiet="$3"
    local rc=0
    local -a stage_dirs=()
    local _d
    while IFS= read -r _d; do stage_dirs+=("$_d"); done < <(_vcheck_stage_dirs "$projdir")

    # --- rule 4: stray files directly under 02_analysis/ --------------------
    # "Known entry point" = docs (README.md/AGENTS.md), build/orchestration
    # (Makefile, run_*.sh), env/dependency manifests (renv.lock,
    # requirements.txt, pyproject.toml, environment.yml, .Rprofile), and
    # dotfiles (.gitkeep, .gitignore, ...). Anything else at the namespace
    # root — e.g. a stray config.py — is neither narrative nor machinery nor
    # config-by-convention, so it warns.
    local analysis_root="$projdir/02_analysis"
    if [[ -d "$analysis_root" ]]; then
        local f base rel
        while IFS= read -r f; do
            [[ -f "$f" ]] || continue
            base=$(basename "$f")
            rel="${f#$projdir/}"
            _vcheck_is_exempt "$rel" && continue
            case "$base" in
                .*) continue ;;
                README.md|AGENTS.md|Makefile) continue ;;
                run_*.sh) continue ;;
                renv.lock|requirements.txt|pyproject.toml|environment.yml|.Rprofile) continue ;;
            esac
            _vcheck_emit "$strict" "$quiet" stage-layout \
                "$rel: stray file directly under 02_analysis/ — not in a subdir and not a known entry point (README.md/AGENTS.md/Makefile/run_*.sh/env manifests)" || rc=1
        done < <(find "$analysis_root" -maxdepth 1 -type f 2>/dev/null)
    fi

    [[ ${#stage_dirs[@]} -gt 0 ]] || return $rc

    # --- rules 1, 2, 3, 5, 6: per-stage-file checks --------------------------
    # Enumeration goes through _vcheck_stage_files so this check agrees with
    # stage-thinness and comment-intent about which files are stages. Rule 3's
    # sibling lookup therefore resolves against each file's OWN directory
    # rather than a single flat stage dir, which is what makes the nested
    # 02_analysis/scripts/{compute,figures}/ layout work.
    local f base rel dir num suffix stem found_sib e
    while IFS= read -r f; do
            [[ -f "$f" ]] || continue
            base=$(basename "$f")
            dir=$(dirname "$f")
            rel="${f#$projdir/}"
            _vcheck_is_exempt "$rel" && continue
            case "$base" in *.R|*.py) ;; *) continue ;; esac

            # rules 5 & 6: stage-number shape (letter suffix / single digit).
            if [[ "$base" =~ ^([0-9]+)([a-zA-Z]?)_ ]]; then
                num="${BASH_REMATCH[1]}"
                suffix="${BASH_REMATCH[2]}"
                if [[ -n "$suffix" ]]; then
                    _vcheck_emit "$strict" "$quiet" stage-layout \
                        "$rel: stage number carries a letter suffix ('$num$suffix') — use decade numbering instead (doc 09 §1.5)" || rc=1
                elif [[ "${#num}" -eq 1 ]]; then
                    _vcheck_emit "$strict" "$quiet" stage-layout \
                        "$rel: single-digit stage number ('$num') — zero-pad to two digits so lexical order matches narrative order (doc 09 §1.5)" || rc=1
                fi
            fi

            # rules 1 & 2: viz/compute write-target contract.
            if [[ "$base" == *_viz.R || "$base" == *_viz.py ]]; then
                if grep -Eq 'write\.csv\(|write\.table\(|saveRDS\(|write_csv\(|write_tsv\(|fwrite\(|\.to_csv\(|\.to_parquet\(|write_parquet\(' "$f" 2>/dev/null; then
                    _vcheck_emit "$strict" "$quiet" stage-layout \
                        "$rel: _viz stage writes tables — viz never computes (doc 09 §1.4)" || rc=1
                fi
            else
                if grep -Eq 'ggsave\(|\.savefig\(|save_figure\(|save_overview\(' "$f" 2>/dev/null; then
                    _vcheck_emit "$strict" "$quiet" stage-layout \
                        "$rel: non-viz stage writes figures — compute never plots (doc 09 §1.4)" || rc=1
                fi
            fi

            # rule 3: orphan viz — no same-number, same-stem compute sibling.
            # Applied literally per doc 09 §3.3, with no minimum-file-count
            # gate: a viz stage with nothing to visualise is exactly the
            # severed pair the rule names, and a lone one is the clearest case
            # of it rather than an exception to it.
            if [[ "$base" =~ ^([0-9]+[a-zA-Z]?)_(.+)_viz\.(R|py)$ ]]; then
                num="${BASH_REMATCH[1]}"; stem="${BASH_REMATCH[2]}"
                found_sib=0
                for e in R py; do
                    [[ -f "$dir/${num}_${stem}.$e" ]] && { found_sib=1; break; }
                done
                if [[ "$found_sib" -eq 0 ]]; then
                    _vcheck_emit "$strict" "$quiet" stage-layout \
                        "$rel: orphan viz — no sibling ${num}_${stem}.{R,py} compute stage (doc 09 §1.4)" || rc=1
                fi
            fi
    done < <(_vcheck_stage_files "$projdir")

    return $rc
}

# ---------------------------------------------------------------------------
# Check: skill-coupling — drift guard on skills' `compatibility:` declarations
# ---------------------------------------------------------------------------
# ADR-D6 chose per-skill `compatibility:` declarations over a curated public
# subset because "declaring a requirement is self-maintaining; maintaining a
# hand-picked list is not". On its own that is not true: nothing forces a
# declaration to track the skill's actual content, so the declaration set
# becomes exactly the second corpus ADR-D1 refuses to keep in sync. `sciagent
# validate` checks the declaration's GRAMMAR (and that toolkit/sibling items
# resolve); it cannot tell whether the declaration is still TRUE. This check
# is the missing half.
#
# Three rules, all warn-only:
#   1. undeclared coupling  — a skill's CODE shows a dependency of some flavour
#                             and the skill declares no clause of that flavour.
#   2. stale declaration    — a declared item is not mentioned ANYWHERE in the
#                             skill directory (outside the declaration itself).
#   3. unscaffolded root    — a `sciagent-scaffold` item whose first path
#                             component is not a directory `sciagent new
#                             project --type analysis` actually creates.
#
# ---- Why it is warn-only even under --strict ------------------------------
# Rules 1 and 2 are heuristics over prose-and-code, and heuristics must not be
# able to block. Note the ONE place strictness could leak: `_lint_run_checks`
# ORs each check's return into its tally, and `_vcheck_emit <strict=1>` returns
# 1. This check therefore calls `_vcheck_emit` with a HARDCODED strict=0 and
# always `return 0` — its <strict> parameter is accepted and ignored on purpose.
#
# ---- Why it is not in `all` -----------------------------------------------
# Every sibling check's subject is --project-dir; this one's subject is the
# TOOLKIT CHECKOUT (`$SCIAGENT_TOOLKIT/skills/`), and it never touches the
# project at all. `all` is the analysis-author's project sweep and is what a
# bare `sciagent lint` runs, so folding a toolkit-corpus audit into it would
# print findings about files the person linting does not own and cannot fix
# from where they are standing. It is opt-in by name: `--check skill-coupling`.
# The <projdir> parameter is likewise accepted and ignored, so the check keeps
# the uniform `_lint_check_<name> <projdir> <strict> <quiet>` signature that
# `_lint_run_checks` dispatches on.
#
# ---- Why the evidence scan reads CODE FILES ONLY, minus comment lines -----
# This is the whole design, and it is calibrated against a known-wrong prior:
# the static scan in 00_INDEX.md §5 flagged 15 skills of which the semantic
# audit found 6 to be false positives — they MENTION the analysis-repo layout,
# they do not REQUIRE it. Reproducing those 6 shows two distinct causes, and
# the scan must defeat both:
#   * provenance citations in `#` comments — `peak-atlas-framework`'s
#     scripts/checks carry `#   /data2/.../02_analysis/config/...` headers
#     crediting the project the code came from; `peak-atlas-unpaired` the same.
#     Hence: a line whose first non-blank characters are `#` or `//` is never
#     evidence.
#   * documentation prose — `scrna-pipeline-conventions` is a skill ABOUT the
#     layout, `anndatar-seurat-scanpy-conversion` and `coresh-signature-search`
#     cite paths in reference notes. Hence: `.md` is never scanned.
# Both are needed; each alone silences only half the six. A middle tier —
# "fenced code blocks inside .md" — was tried and REJECTED: it re-flags
# `coresh-signature-search` (references/coresh-to-gsea-bridge.md:76, a
# `write_gmt(..., "03_results/...")` example) and `scrna-pipeline-conventions`
# (SKILL.md's illustrative `03_results/objects/*.h5ad` tree), i.e. it fails on
# exactly the corpus the exclusion exists for.
#
# The price is recall, stated plainly: a skill whose coupling lives only in
# prose (`figure-style`, `reasoning-trace`) is invisible to rule 1. That is the
# deliberate trade — a drift guard that cries wolf gets ignored, and prose is
# where the wolves were. Rule 1 catches the case that matters most anyway: a
# skill shipping RUNNABLE code against the scaffold without saying so.
#
# Rule 2 is deliberately asymmetric to rule 1: it acquits on a mention
# ANYWHERE in the skill (prose, comments, any file type), because the question
# it asks is the opposite one — not "is this required?" but "is there any trace
# of this at all?". A declaration naming something the skill never once
# mentions is dead text by any reading.
#
# Calibration over the shipped corpus (87 skill dirs, 15 declaring): rule 1 = 0
# findings, rule 2 = 0 findings, rule 3 = 1 finding (`iterative-peak-merging`
# declares `01_scripts/...`, which the toolkit does not scaffold). Zero on
# rules 1/2 is the CORRECT reading of a corpus audited days earlier, not a
# broken check — the mutation tests in tests/test_lint_skill_coupling.sh pin
# both directions.

# Path roots `sciagent new project --type analysis` actually creates: the
# `00_data/... 02_analysis/... 03_results/...` list in new.sh's _new_project
# plus the `docs/` tree it scaffolds alongside. `01_modules/` is created too
# but is deliberately ABSENT here — docs/packaged-skills.md §6 forbids filing a
# sibling submodule under `sciagent-scaffold` ("SciAgent does not ship it, and
# the declaration would be a lie"), so a `01_modules/` item earns its own
# use-external-module message below rather than passing as a scaffold root.
_VCHECK_SCAFFOLD_ROOTS=(00_data 02_analysis 03_results docs)

# Evidence regexes, one per flavour. Anchored with a leading
# `(^|[^A-Za-z0-9_./-])` so an ABSOLUTE foreign path (`/data2/users/.../
# 02_analysis/config/`) does not read as a repo-root-relative requirement —
# those appear in the corpus and are provenance, not coupling.
_VCHECK_COUPLE_RE_SCAFFOLD='(^|[^A-Za-z0-9_./-])(00_data|02_analysis|03_results|docs/_internal)/'
# The two directories symlink_create_helper_lib mounts into 02_analysis/helpers/
# (symlinks.sh's literal `for libdir in figure-style interactive-style`), by
# their mounted shim spelling, plus the figure-style contract entry points a
# skill can only call if that shim is present.
_VCHECK_COUPLE_RE_TOOLKIT='(figure_style|interactive_style|project_theme\(|set_paper_style\(|save_figure\(|save_overview\()'

# _vcheck_skill_code_files <skilldir>
# Emit every CODE file under a skill, one per line. `.md` is excluded by
# construction (see the header: prose is not evidence). Build/venv residue is
# pruned — skills/mllmcelltype-consensus-annotation/.venv alone is 171MB.
_vcheck_skill_code_files() {
    local d="$1"
    find "$d" \
        \( -name '.venv' -o -name '__pycache__' -o -name 'node_modules' -o -name '.git' \) -prune -o \
        -type f \( -name '*.R' -o -name '*.r' -o -name '*.py' -o -name '*.sh' \
                   -o -name '*.bash' -o -name '*.js' -o -name '*.yaml' -o -name '*.yml' \) \
        -print 2>/dev/null | sort
}

# _vcheck_skill_evidence <skilldir> <ere>
# Echo `<file>:<line>` for the FIRST non-comment code line matching <ere>;
# return 1 when there is none. A line whose first non-blank characters are `#`
# (R/Python/shell/YAML) or `//` (JS) is skipped — see the header block.
# Known, accepted blind spot: a `#` inside a string literal, and Python
# docstrings, are not distinguished from real code (the same approximation
# _vcheck_stage_scan_defs and _lint_check_comment_intent already make).
_vcheck_skill_evidence() {
    local d="$1" ere="$2"
    local -a files=()
    local f
    while IFS= read -r f; do files+=("$f"); done < <(_vcheck_skill_code_files "$d")
    [[ ${#files[@]} -gt 0 ]] || return 1
    local hit
    hit=$(awk -v pat="$ere" '
        /^[ \t]*(#|\/\/)/ { next }
        $0 ~ pat { printf "%s:%d\n", FILENAME, FNR; exit }
    ' "${files[@]}" 2>/dev/null)
    [[ -n "$hit" ]] || return 1
    printf '%s' "$hit"
    return 0
}

# _vcheck_skill_dirrefs <skilldir> <prefix-ere>
# Echo `<name>|<file>:<line>` for every non-comment code line referencing
# `<prefix>/<name>/`, where <prefix-ere> is an alternation like `skills|\.\.`.
# The caller decides which names are meaningful (an existing sibling skill, an
# existing 01_modules entry) — this only harvests candidates. That split is why
# `skill-creator`'s `python utils/package_skill.py skills/public/my-skill` usage
# string is silent: `public` is not a skill directory, so the caller drops it.
_vcheck_skill_dirrefs() {
    local d="$1" prefix="$2"
    local -a files=()
    local f
    while IFS= read -r f; do files+=("$f"); done < <(_vcheck_skill_code_files "$d")
    [[ ${#files[@]} -gt 0 ]] || return 0
    awk -v pfx="$prefix" '
        /^[ \t]*(#|\/\/)/ { next }
        {
            line = $0
            re = "(" pfx ")/[A-Za-z0-9][A-Za-z0-9._-]*/"
            while (match(line, re)) {
                tok = substr(line, RSTART, RLENGTH)
                sub("^[^/]*/", "", tok)     # drop the prefix
                sub("/$", "", tok)          # drop the trailing slash
                if (!(tok in seen)) {
                    seen[tok] = 1
                    printf "%s|%s:%d\n", tok, FILENAME, FNR
                }
                line = substr(line, RSTART + RLENGTH)
            }
        }
    ' "${files[@]}" 2>/dev/null
    return 0
}

# _vcheck_compat_items <skill_md>
# Echo `<flavour><TAB><item>` for each declared item; nothing when the skill
# declares no `compatibility:`.
#
# This parses LENIENTLY and on purpose: grammar enforcement belongs to
# `sciagent validate` (`_validate_compat_parse`), which this file must not
# depend on — validate.sh depends on lint.sh, not the reverse, and inverting
# that would create the one source-time cycle bin/sciagent's load graph is
# built to make impossible. A malformed declaration is validate's finding; here
# it simply yields whatever splits cleanly, so the drift guard never
# double-reports a grammar defect in different words.
_vcheck_compat_items() {
    local skill_md="$1" fm val
    fm=$(_fm_extract "$skill_md")
    printf '%s\n' "$fm" | grep -q '^compatibility:' || return 0
    val=$(_fm_scalar compatibility <<< "$fm")
    [[ -n "$val" ]] || return 0
    printf '%s\n' "$val" | awk '
        {
            nc = split($0, cl, ";")
            for (i = 1; i <= nc; i++) {
                c = cl[i]
                sub(/^[ \t]+/, "", c); sub(/[ \t]+$/, "", c)
                p = index(c, ":")
                if (p == 0) continue
                fl = substr(c, 1, p - 1)
                it = substr(c, p + 1)
                ni = split(it, items, ",")
                for (j = 1; j <= ni; j++) {
                    x = items[j]
                    sub(/^[ \t]+/, "", x); sub(/[ \t]+$/, "", x)
                    if (x != "") printf "%s\t%s\n", fl, x
                }
            }
        }'
    return 0
}

# _vcheck_skill_mentions <skilldir> <needle>...
# True (0) when any <needle> appears as a literal substring anywhere under the
# skill — ANY file type, comments and prose included (see the header: rule 2
# acquits broadly by design). The declaration line itself is excluded, since a
# clause corroborating only itself is precisely the dead text being looked for.
# `grep -I` skips binaries (.pyc, .sam fixtures).
_vcheck_skill_mentions() {
    local d="$1"; shift
    local -a pats=()
    local n
    for n in "$@"; do pats+=(-e "$n"); done
    grep -rIF "${pats[@]}" "$d" \
        --exclude-dir='.venv' --exclude-dir='__pycache__' --exclude-dir='node_modules' \
        2>/dev/null | grep -qv '/SKILL\.md:compatibility:'
}

# _vcheck_compat_has_flavour <want> [<flavour>...]
# True (0) when <want> is among the declared flavours.
_vcheck_compat_has_flavour() {
    local want="$1"; shift
    local f
    for f in "$@"; do
        [[ "$f" == "$want" ]] && return 0
    done
    return 1
}

# _lint_check_skill_coupling <projdir-IGNORED> <strict-IGNORED> <quiet>
_lint_check_skill_coupling() {
    local quiet="$3"
    local tk_root
    tk_root="$(_vcheck_toolkit_root)"
    [[ -d "$tk_root/skills" ]] || return 0

    local d skill skill_md rel
    for d in "$tk_root"/skills/*/; do
        [[ -d "$d" ]] || continue
        skill="$(basename "$d")"
        # `_archive`, `_attic`, `_TEMPLATE` are not shipped skills.
        case "$skill" in _*) continue ;; esac
        skill_md="$d/SKILL.md"
        [[ -f "$skill_md" ]] || continue
        rel="skills/$skill"

        # --- read the declaration -----------------------------------------
        local -a decl_flavours=() decl_items=()
        local fl it
        while IFS=$'\t' read -r fl it; do
            [[ -n "$fl" ]] || continue
            decl_flavours+=("$fl")
            decl_items+=("$fl$(printf '\t')$it")
        done < <(_vcheck_compat_items "$skill_md")

        # --- rule 1: undeclared coupling ----------------------------------
        local hit
        if ! _vcheck_compat_has_flavour sciagent-scaffold ${decl_flavours[@]+"${decl_flavours[@]}"}; then
            if hit=$(_vcheck_skill_evidence "$d" "$_VCHECK_COUPLE_RE_SCAFFOLD"); then
                _vcheck_emit 0 "$quiet" skill-coupling \
                    "$skill: undeclared sciagent-scaffold coupling — ${hit#$tk_root/} references the analysis-repo layout, but SKILL.md declares no 'sciagent-scaffold:' clause"
            fi
        fi
        if ! _vcheck_compat_has_flavour sciagent-toolkit ${decl_flavours[@]+"${decl_flavours[@]}"}; then
            if hit=$(_vcheck_skill_evidence "$d" "$_VCHECK_COUPLE_RE_TOOLKIT"); then
                _vcheck_emit 0 "$quiet" skill-coupling \
                    "$skill: undeclared sciagent-toolkit coupling — ${hit#$tk_root/} uses a helper the toolkit mounts (figure-style/interactive-style), but SKILL.md declares no 'sciagent-toolkit:' clause"
            fi
        fi
        local ref name where
        if ! _vcheck_compat_has_flavour sibling-skill ${decl_flavours[@]+"${decl_flavours[@]}"}; then
            while IFS= read -r ref; do
                [[ -n "$ref" ]] || continue
                name="${ref%%|*}"; where="${ref#*|}"
                [[ "$name" == "$skill" ]] && continue
                [[ -f "$tk_root/skills/$name/SKILL.md" ]] || continue
                _vcheck_emit 0 "$quiet" skill-coupling \
                    "$skill: undeclared sibling-skill coupling — ${where#$tk_root/} reaches into skills/$name/, but SKILL.md declares no 'sibling-skill:' clause"
            done < <(_vcheck_skill_dirrefs "$d" 'skills|\.\.')
        fi
        if ! _vcheck_compat_has_flavour external-module ${decl_flavours[@]+"${decl_flavours[@]}"}; then
            while IFS= read -r ref; do
                [[ -n "$ref" ]] || continue
                name="${ref%%|*}"; where="${ref#*|}"
                _vcheck_emit 0 "$quiet" skill-coupling \
                    "$skill: undeclared external-module coupling — ${where#$tk_root/} reads 01_modules/$name, but SKILL.md declares no 'external-module:' clause"
            done < <(_vcheck_skill_dirrefs "$d" '01_modules')
        fi

        # --- rules 2 and 3: per declared item ------------------------------
        local pair item root
        for pair in ${decl_items[@]+"${decl_items[@]}"}; do
            fl="${pair%%$'\t'*}"
            item="${pair#*$'\t'}"

            # rule 3: scaffold item outside the scaffolded roots.
            if [[ "$fl" == "sciagent-scaffold" ]]; then
                root="${item%%/*}"
                if [[ "$root" == "01_modules" ]]; then
                    _vcheck_emit 0 "$quiet" skill-coupling \
                        "$skill: sciagent-scaffold '$item' names a sibling submodule — SciAgent does not ship it; declare it as 'external-module' (docs/packaged-skills.md §6)"
                else
                    local known=0 r
                    for r in "${_VCHECK_SCAFFOLD_ROOTS[@]}"; do
                        [[ "$root" == "$r" ]] && { known=1; break; }
                    done
                    if [[ "$known" -eq 0 ]]; then
                        _vcheck_emit 0 "$quiet" skill-coupling \
                            "$skill: sciagent-scaffold '$item' is rooted at '$root/', which 'sciagent new project --type analysis' does not create — if the path is project-local, 'external-module' is truer"
                    fi
                fi
            fi

            # rule 2: the declared item is mentioned nowhere in the skill.
            # Three spellings are accepted, each for a corpus reason:
            #   * the item with any trailing '/' dropped, so a directory item
            #     matches a mention that omits it;
            #   * for `sciagent-toolkit`, the mounted underscore spelling —
            #     lib/figure-style arrives as 02_analysis/helpers/figure_style.R;
            #   * for a file item, the extension-less stem. `figure-style`
            #     documents its two shims as `02_analysis/helpers/
            #     figure_style.{R,py}` in brace shorthand, which no literal
            #     search for the `.py` member can ever match. Without this the
            #     check's very first run reported a stale declaration that is
            #     plainly documented one line below it.
            local needle alt stem
            needle="${item%/}"
            alt="$needle"
            [[ "$fl" == "sciagent-toolkit" ]] && alt="${needle//-/_}"
            stem="$needle"
            [[ "$(basename "$needle")" == *.* ]] && stem="${needle%.*}"
            if ! _vcheck_skill_mentions "$d" "$needle" "$alt" "$stem"; then
                _vcheck_emit 0 "$quiet" skill-coupling \
                    "$skill: stale declaration '$fl: $item' — '$needle' appears nowhere under $rel/ outside the declaration itself"
            fi
        done
    done

    # Never returns non-zero: this check must not be able to block anything,
    # including under --strict. See the header block.
    return 0
}
